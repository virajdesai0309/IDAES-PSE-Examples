# Import objects from pyomo package
from pyomo.environ import ConcreteModel, Block, Constraint, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.base.property_base import StateBlockData

# Import the feed and product unit models and the tray column model
from idaes.models.unit_models import Feed, Product
from idaes.models_extra.column_models import TrayColumn
from idaes.models_extra.column_models.condenser import CondenserType, TemperatureSpec

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for the flowsheet
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the benzene-toluene property package
from benzene_toluene import configuration

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import propagate_state to pass the initialized states between unit models
from idaes.core.util.initialization import propagate_state

# Import Arc to connect the unit models
from pyomo.network import Arc

import re
import zipfile

import pandas as pd

COMPONENTS = ["benzene", "toluene"]


def setup_distillation_model(
    flow_mass=1.0,  # kg/s
    mole_frac={"benzene": 0.5, "toluene": 0.5},
    temperature=300,  # K
    pressure=101325,  # Pa
    number_of_stages=13,  # DWSIM numbering: condenser = 1, reboiler = last
    feed_stage=8,  # DWSIM numbering
    reflux_ratio=1.5,
    bottoms_flow=5.8,  # mol/s
    condenser_pressure=101325,  # Pa
    column_pressure_drop=10,  # Pa, condenser to reboiler
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a rigorous benzene-toluene distillation column.

    DWSIM counts the condenser and reboiler as stages, IDAES counts only the
    trays between them: IDAES tray k is DWSIM stage k + 1.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)

    number_of_trays = number_of_stages - 2
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.column = TrayColumn(
        number_of_trays=number_of_trays,
        feed_tray_location=feed_stage - 1,
        condenser_type=CondenserType.totalCondenser,
        condenser_temperature_spec=TemperatureSpec.atBubblePoint,
        has_pressure_change=True,
        property_package=m.fs.properties,
    )
    m.fs.distillate = Product(property_package=m.fs.properties)
    m.fs.bottoms = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.column.feed)
    m.fs.s02 = Arc(source=m.fs.column.condenser.distillate, destination=m.fs.distillate.inlet)
    m.fs.s03 = Arc(source=m.fs.column.reboiler.bottoms, destination=m.fs.bottoms.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    mw = sum(x * value(m.fs.properties.get_component(j).mw) for j, x in mole_frac.items())  # kg/mol
    m.fs.feed.flow_mol.fix(flow_mass / mw)  # mol/s
    m.fs.feed.temperature.fix(temperature)  # K
    m.fs.feed.pressure.fix(pressure)  # Pa
    for j, x in mole_frac.items():
        m.fs.feed.mole_frac_comp[0, j].fix(x)

    c = m.fs.column
    c.condenser.reflux_ratio.fix(reflux_ratio)
    c.condenser.condenser_pressure.fix(condenser_pressure)

    # Linear pressure profile, as DWSIM: equal pressure rise from stage to stage.
    # Tray pressure = pressure of the liquid entering from above - deltaP
    dp_stage = column_pressure_drop / (number_of_stages - 1)  # Pa
    for tray in c.tray_index:
        tray_block = c.feed_tray if tray == c.config.feed_tray_location else (
            c.rectification_section[tray] if tray < c.config.feed_tray_location else c.stripping_section[tray])
        tray_block.deltaP.fix(-dp_stage)
    c.reboiler.deltaP.fix(dp_stage)

    # Start from a boilup ratio guess; the DWSIM bottoms flow spec replaces it at the end
    c.reboiler.boilup_ratio.fix(3.0)
    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    # The tray, condenser and reboiler initialization routines use the bubble and
    # dew temperatures of their state blocks. Build them before initializing, so
    # that the property package initializes them too.
    states = [b for b in m.fs.component_data_objects(Block, descend_into=True) if isinstance(b, StateBlockData)]
    for state in states:
        state.temperature_bubble
        state.temperature_dew

    # 1. Initialize with the feed near its bubble point (365 K), which the
    #    tray-by-tray initialization handles well
    m.fs.feed.temperature.fix(365)
    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    c.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    propagate_state(arc=m.fs.s03)
    m.fs.distillate.initialize(outlvl=outlvl)
    m.fs.bottoms.initialize(outlvl=outlvl)

    solver = get_solver(options={"max_iter": 500})
    solver.solve(m, tee=False)

    # 2. Tighten the smoothing of the vapor-liquid complementarity. The default
    #    (1 mol/s) lets a vanishing phase linger, e.g. vapor in the total
    #    condenser outlet at its bubble point.
    for eps in (1e-1, 1e-2, 1e-3, 1e-4):
        for state in states:
            for name in ("eps_t_Vap_Liq", "eps_z_Vap_Liq"):
                if hasattr(state, name):
                    getattr(state, name).value = eps
        solver.solve(m, tee=False)

    # 3. Step the feed temperature down to the actual (subcooled) feed
    for T in (350, 335, 320, temperature):
        m.fs.feed.temperature.fix(T)
        solver.solve(m, tee=False)

    # 4. Switch to the DWSIM specification: bottoms molar flow instead of boilup ratio
    c.reboiler.boilup_ratio.unfix()
    m.fs.bottoms_flow_spec = Constraint(expr=c.reboiler.bottoms.flow_mol[0] == bottoms_flow)  # mol/s
    result = solver.solve(m, tee=False)

    return m, result


def report_distillation_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    distillate_prop = model.fs.distillate.report()
    bottoms_prop = model.fs.bottoms.report()

    return feed_prop, distillate_prop, bottoms_prop


def stage_blocks(model):
    """Stage blocks in DWSIM numbering: condenser, trays, reboiler."""
    c = model.fs.column
    stages = [c.condenser]
    for tray in c.tray_index:
        if tray == c.config.feed_tray_location:
            stages.append(c.feed_tray)
        elif tray < c.config.feed_tray_location:
            stages.append(c.rectification_section[tray])
        else:
            stages.append(c.stripping_section[tray])
    stages.append(c.reboiler)
    return stages


def profiles(model):
    """Temperature, pressure, flows and benzene mole fractions on every stage (DWSIM numbering)."""
    rows = []
    stages = stage_blocks(model)
    for n, stage in enumerate(stages, start=1):
        if n == 1:  # total condenser: liquid leaves at the bubble point; liquid flow is the reflux, as DWSIM
            state = stage.control_volume.properties_out[0]
            liquid = stage.reflux.flow_mol[0]
            vapor_out = 0
            x = stage.reflux.mole_frac_comp[0, "benzene"]
            y = stages[1].vap_out.mole_frac_comp[0, "benzene"]  # vapor entering the condenser
        elif n == len(stages):  # reboiler
            state = stage.control_volume.properties_out[0]
            liquid = stage.bottoms.flow_mol[0]
            vapor_out = stage.vapor_reboil.flow_mol[0]
            x = stage.bottoms.mole_frac_comp[0, "benzene"]
            y = stage.vapor_reboil.mole_frac_comp[0, "benzene"]
        else:
            state = stage.properties_out[0]
            liquid = stage.liq_out.flow_mol[0]
            vapor_out = stage.vap_out.flow_mol[0]
            x = stage.liq_out.mole_frac_comp[0, "benzene"]
            y = stage.vap_out.mole_frac_comp[0, "benzene"]
        rows.append({
            "Stage": n,
            "Temperature (K)": value(state.temperature),
            "Pressure (Pa)": value(state.pressure),
            "Liquid flow (mol/s)": value(liquid),
            "Vapor flow (mol/s)": value(vapor_out),
            "x benzene": value(x),
            "y benzene": value(y),
        })
    return pd.DataFrame(rows).set_index("Stage")


def column_results(model):
    """Key column results (mol/s, K, kW)."""
    c = model.fs.column
    dist, bott = c.condenser.distillate, c.reboiler.bottoms
    return {
        "Distillate flow (mol/s)": value(dist.flow_mol[0]),
        "Bottoms flow (mol/s)": value(bott.flow_mol[0]),
        "Distillate benzene mole fraction": value(dist.mole_frac_comp[0, "benzene"]),
        "Bottoms benzene mole fraction": value(bott.mole_frac_comp[0, "benzene"]),
        "Distillate temperature (K)": value(dist.temperature[0]),
        "Bottoms temperature (K)": value(bott.temperature[0]),
        "Condenser duty (kW)": -value(c.condenser.heat_duty[0]) / 1000,
        "Reboiler duty (kW)": value(c.reboiler.heat_duty[0]) / 1000,
        "Boilup ratio": value(c.reboiler.boilup_ratio),
    }


def read_dwsim_results(path):
    """Read the column results and stage profiles from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Streams are found
    by their tag; the converged stage profiles are stored with the column.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    column = xml[xml.find("UnitOperations.DistillationColumn<"):]

    # Stream tags (names shown on the flowsheet) are stored with the graphic objects
    tags = dict(re.findall(r"<Name>(MAT-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        streams[tags["MAT-" + stream[:36]]] = stream

    def mixture(stream):
        return stream[stream.find("<ComponentName>Mixture<"):]

    def benzene(stream):
        values = re.findall(r"<ComponentName>Benzene</ComponentName>.*?<MoleFraction>([^<]+)<", stream, re.S)
        return next(float(v) for v in values if v != "NaN")

    def vector(tag):
        return [float(v) for v in re.search(r"<{0}>\{{([^}}]*)\}}".format(tag), column).group(1).split(";")]

    def compositions(tag):
        block = re.search(r"<{0}>(.*?)</{0}>".format(tag), column, re.S).group(1)
        return [float(d.split(";")[0]) for d in re.findall(r"<data>\{([^}]*)\}</data>", block)]

    dist, bott = streams["distillate"], streams["residue"]
    results = {
        "Distillate flow (mol/s)": first("molarflow", mixture(dist)),
        "Bottoms flow (mol/s)": first("molarflow", mixture(bott)),
        "Distillate benzene mole fraction": benzene(dist),
        "Bottoms benzene mole fraction": benzene(bott),
        "Distillate temperature (K)": first("temperature", mixture(dist)),
        "Bottoms temperature (K)": first("temperature", mixture(bott)),
        "Condenser duty (kW)": first("CondenserDuty", column),
        "Reboiler duty (kW)": -first("ReboilerDuty", column),
    }
    V = vector("Vf")
    results["Boilup ratio"] = V[-1] / results["Bottoms flow (mol/s)"]

    stages = pd.DataFrame({
        "Temperature (K)": vector("Tf"),
        "Liquid flow (mol/s)": vector("Lf"),
        "Vapor flow (mol/s)": V,
        "x benzene": compositions("xf"),
        "y benzene": compositions("yf"),
    }, index=pd.RangeIndex(1, len(V) + 1, name="Stage"))
    return results, stages


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": column_results(model), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"] * 100
    return df


def compare_profiles(model, dwsim_stages):
    """Stage temperatures and benzene mole fractions, IDAES vs DWSIM."""
    p = profiles(model)
    return pd.DataFrame({
        "T IDAES (K)": p["Temperature (K)"],
        "T DWSIM (K)": dwsim_stages["Temperature (K)"],
        "T difference (K)": p["Temperature (K)"] - dwsim_stages["Temperature (K)"],
        "x IDAES": p["x benzene"],
        "x DWSIM": dwsim_stages["x benzene"],
    })


if __name__ == "__main__":
    model, result = setup_distillation_model()
    print(result.solver.termination_condition)
    report_distillation_properties(model)
    dwsim, dwsim_stages = read_dwsim_results("DWSim Distillation Column Model.dwxmz")
    print(compare_with_dwsim(model, dwsim))
    print(compare_profiles(model, dwsim_stages))
