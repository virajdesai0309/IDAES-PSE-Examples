# Import objects from pyomo package
from pyomo.environ import ConcreteModel, Constraint, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product, translator, mixer and equilibrium reactor unit models
from idaes.models.unit_models import (Feed, Product, Translator, Mixer,
                                      MomentumMixingType, EquilibriumReactor)

# Import the IAPWS-95 steam table property package for the steam feed
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property and reaction packages
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock

# Import the syngas property package and the water-gas shift reaction package
from syngas import configuration as gas_configuration
from wgs_reaction import configuration as reaction_configuration

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

COMPONENTS = ["carbon_monoxide", "water", "carbon_dioxide", "hydrogen"]
MW_WATER = 0.01801528  # kg/mol

DEFAULT_STEAM = {"flow_mass": 5.0, "temperature": 500, "pressure": 101325}  # kg/s, K, Pa
DEFAULT_HYDROGEN = {"flow_mass": 0.5, "temperature": 300, "pressure": 101325, "mole_frac": {"hydrogen": 1.0}}
DEFAULT_FLUE_GAS = {"flow_mass": 1.0, "temperature": 300, "pressure": 101325,
                    "mole_frac": {"carbon_monoxide": 0.5, "carbon_dioxide": 0.5}}


def setup_equilibrium_reactor_model(
    steam=DEFAULT_STEAM,
    hydrogen=DEFAULT_HYDROGEN,
    flue_gas=DEFAULT_FLUE_GAS,
    reactor_temperature=700,  # K
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a water-gas shift equilibrium reactor at a fixed temperature."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.steam_properties = iapws95.Iapws95ParameterBlock()
    m.fs.properties = GenericParameterBlock(**gas_configuration)
    m.fs.reactions = GenericReactionParameterBlock(
        property_package=m.fs.properties, **reaction_configuration
    )

    m.fs.steam = Feed(property_package=m.fs.steam_properties)
    # The translator hands the steam over from the steam tables to the syngas package
    m.fs.translator = Translator(
        inlet_property_package=m.fs.steam_properties,
        outlet_property_package=m.fs.properties,
        outlet_state_defined=True,
        has_phase_equilibrium=False,
    )
    m.fs.hydrogen = Feed(property_package=m.fs.properties)
    m.fs.flue_gas = Feed(property_package=m.fs.properties)
    m.fs.mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=["steam_inlet", "hydrogen_inlet", "flue_gas_inlet"],
        momentum_mixing_type=MomentumMixingType.minimize,
    )
    m.fs.reactor = EquilibriumReactor(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_rate_reactions=False,
        has_equilibrium_reactions=True,
        # Heats of formation are in the enthalpies, so no separate heat of reaction term
        has_heat_of_reaction=False,
        # Heat duty needed to hold the outlet temperature, as DWSIM "Outlet Temperature" mode
        has_heat_transfer=True,
        has_pressure_change=False,
    )
    m.fs.syngas = Product(property_package=m.fs.properties)

    """Translator constraints: same flow, temperature and pressure; pure water"""
    t_in = m.fs.translator.properties_in[0]
    t_out = m.fs.translator.properties_out[0]
    m.fs.translator.eq_flow = Constraint(expr=t_out.flow_mol == t_in.flow_mol)
    m.fs.translator.eq_temperature = Constraint(expr=t_out.temperature == t_in.temperature)
    m.fs.translator.eq_pressure = Constraint(expr=t_out.pressure == t_in.pressure)
    for j in COMPONENTS:
        # Small non-zero values avoid log(0) in the Peng-Robinson expressions
        t_out.mole_frac_comp[j].fix(1.0 if j == "water" else 1e-10)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.steam.outlet, destination=m.fs.translator.inlet)
    m.fs.s02 = Arc(source=m.fs.translator.outlet, destination=m.fs.mixer.steam_inlet)
    m.fs.s03 = Arc(source=m.fs.hydrogen.outlet, destination=m.fs.mixer.hydrogen_inlet)
    m.fs.s04 = Arc(source=m.fs.flue_gas.outlet, destination=m.fs.mixer.flue_gas_inlet)
    m.fs.s05 = Arc(source=m.fs.mixer.outlet, destination=m.fs.reactor.inlet)
    m.fs.s06 = Arc(source=m.fs.reactor.outlet, destination=m.fs.syngas.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.steam.flow_mol.fix(steam["flow_mass"] / MW_WATER)  # mol/s
    m.fs.steam.pressure.fix(steam["pressure"])  # Pa
    m.fs.steam.enth_mol.fix(
        value(htpx(T=steam["temperature"] * pyunits.K, P=steam["pressure"] * pyunits.Pa))
    )  # J/mol

    for feed, spec in ((m.fs.hydrogen, hydrogen), (m.fs.flue_gas, flue_gas)):
        mw = sum(
            spec["mole_frac"].get(j, 0) * value(m.fs.properties.get_component(j).mw)
            for j in COMPONENTS
        )  # kg/mol
        feed.flow_mol.fix(spec["flow_mass"] / mw)  # mol/s
        feed.temperature.fix(spec["temperature"])  # K
        feed.pressure.fix(spec["pressure"])  # Pa
        for j in COMPONENTS:
            feed.mole_frac_comp[0, j].fix(max(spec["mole_frac"].get(j, 0), 1e-10))

    m.fs.reactor.outlet.temperature.fix(reactor_temperature)

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.steam.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    m.fs.translator.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    m.fs.hydrogen.initialize(outlvl=outlvl)
    m.fs.flue_gas.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s03)
    propagate_state(arc=m.fs.s04)
    m.fs.mixer.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s05)
    m.fs.reactor.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s06)
    m.fs.syngas.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_equilibrium_reactor_properties(model):
    """Report properties of the solved model."""
    steam_prop = model.fs.steam.report()
    hydrogen_prop = model.fs.hydrogen.report()
    flue_gas_prop = model.fs.flue_gas.report()
    mixer_prop = model.fs.mixer.report()
    reactor_prop = model.fs.reactor.report()
    syngas_prop = model.fs.syngas.report()

    return steam_prop, hydrogen_prop, flue_gas_prop, mixer_prop, reactor_prop, syngas_prop


def reactor_results(model):
    """Key results (K, kW, mol/s)."""
    cv = model.fs.reactor.control_volume
    inlet, outlet = cv.properties_in[0], cv.properties_out[0]
    results = {
        "Mixed feed temperature (K)": value(inlet.temperature),
        "Reactor temperature (K)": value(outlet.temperature),
        "Reactor heat duty (kW)": value(cv.heat[0]) / 1000,
        "Syngas molar flow (mol/s)": value(outlet.flow_mol),
    }
    for j in COMPONENTS:
        results["Syngas " + j + " (mol/s)"] = value(outlet.flow_mol_comp[j])
    results["CO conversion"] = 1 - value(outlet.flow_mol_comp["carbon_monoxide"] / inlet.flow_mol_comp["carbon_monoxide"])
    return results


def read_dwsim_results(path):
    """Read the mixer and reactor results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Streams are found
    by their tag. Each compound is listed once per phase; the first numeric
    molar flow is the overall stream value.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    names = {
        "carbon_monoxide": "Carbon monoxide",
        "water": "Water",
        "carbon_dioxide": "Carbon dioxide",
        "hydrogen": "Hydrogen",
    }

    # Stream tags (names shown on the flowsheet) are stored with the graphic objects
    tags = dict(re.findall(r"<Name>(MAT-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        streams[tags["MAT-" + stream[:36]]] = stream

    def mixture(stream):
        return stream[stream.find("<ComponentName>Mixture<"):]

    def component_flow(stream, name):
        flows = re.findall(
            r"<ComponentName>{0}</ComponentName>.*?<MolarFlow>([^<]+)<".format(name), stream, re.S
        )
        return next(float(f) for f in flows if f != "NaN")

    reactor = xml[xml.find("Reactors.Reactor_Equilibrium<"):]
    feed, syngas = streams["feed"], streams["syngas"]
    results = {
        "Mixed feed temperature (K)": first("temperature", mixture(feed)),
        "Reactor temperature (K)": first("temperature", mixture(syngas)),
        "Reactor heat duty (kW)": first("DeltaQ", reactor),
        "Syngas molar flow (mol/s)": first("molarflow", mixture(syngas)),
    }
    for j, name in names.items():
        results["Syngas " + j + " (mol/s)"] = component_flow(syngas, name)
    results["CO conversion"] = 1 - component_flow(syngas, names["carbon_monoxide"]) / component_flow(
        feed, names["carbon_monoxide"])
    return results


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": reactor_results(model), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"].where(df["DWSIM"].abs() > 1e-6) * 100
    return df


if __name__ == "__main__":
    model, result = setup_equilibrium_reactor_model()
    print(result.solver.termination_condition)
    report_equilibrium_reactor_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Equilibrium Reactor Model.dwxmz")))
