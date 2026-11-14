# Import objects from pyomo package
from pyomo.environ import ConcreteModel, Constraint, Var, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the unit models: feeds and products for the water streams, and the refrigeration cycle units
from idaes.models.unit_models import Feed, Product, HeatExchanger
from idaes.models.unit_models.pressure_changer import Compressor
from idaes.models.unit_models.valve import Valve
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_smooth_callback

# Import the property packages: R134a (Helmholtz reference equation of state) and water (IAPWS-95)
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
)
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx

# Import idaes logger to set output levels
import idaes.logger as idaeslog

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

MW_WATER = 0.01801528  # kg/mol

# Design point of the DWSIM model
DESIGN = {
    "refrigerant_flow": 0.1375,  # kg/s R134a
    "compressor_efficiency": 0.75,
    "superheat": 5.0,  # K at the evaporator outlet (expansion valve setting)
    "subcooling": 3.0,  # K at the condenser outlet
    "condenser_U": 1000.0,  # W/m2/K, clean
    "condenser_area": 1.97,  # m2
    "evaporator_U": 800.0,  # W/m2/K, clean
    "evaporator_area": 5.249,  # m2
    "chilled_water_flow": 1.0,  # kg/s
    "chilled_water_return": 12.0,  # C
    "cooling_water_flow": 1.22,  # kg/s
    "cooling_water_inlet": 30.0,  # C
    "water_pressure": 101325,  # Pa
    # Initial guesses (the DWSIM operating pressures)
    "evaporator_pressure": 3.15e5,  # Pa
    "condenser_pressure": 10.17e5,  # Pa
}


def build_chiller(design=DESIGN):
    """Build the closed loop chiller flowsheet (not yet initialized or solved).

    Refrigerant loop: compressor -> condenser -> expansion valve -> evaporator -> compressor
    Chilled water flows through the evaporator, cooling water through the condenser.

    The condensing and evaporating pressures are not fixed: they float to give
    the specified condenser subcooling and evaporator superheat, as in a real
    chiller with an expansion valve. Fouling enters through the overall heat
    transfer coefficient, U = 1 / (1/U_clean + R_f).
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.r134a = HelmholtzParameterBlock(
        pure_component="r134a",
        phase_presentation=PhaseType.MIX,
        state_vars=StateVars.PH,
        amount_basis=AmountBasis.MOLE,
    )
    m.fs.water = iapws95.Iapws95ParameterBlock()

    # Refrigeration cycle
    m.fs.compressor = Compressor(property_package=m.fs.r134a)
    m.fs.condenser = HeatExchanger(
        delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
        hot_side={"property_package": m.fs.r134a},
        cold_side={"property_package": m.fs.water},
    )
    m.fs.expansion_valve = Valve(property_package=m.fs.r134a)
    m.fs.evaporator = HeatExchanger(
        delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
        hot_side={"property_package": m.fs.water},
        cold_side={"property_package": m.fs.r134a},
    )

    # Water streams
    m.fs.chilled_water_in = Feed(property_package=m.fs.water)
    m.fs.chilled_water_out = Product(property_package=m.fs.water)
    m.fs.cooling_water_in = Feed(property_package=m.fs.water)
    m.fs.cooling_water_out = Product(property_package=m.fs.water)

    """Connecting unit operations"""
    m.fs.s02 = Arc(source=m.fs.compressor.outlet, destination=m.fs.condenser.hot_side_inlet)
    m.fs.s05 = Arc(source=m.fs.condenser.hot_side_outlet, destination=m.fs.expansion_valve.inlet)
    m.fs.s07 = Arc(source=m.fs.expansion_valve.outlet, destination=m.fs.evaporator.cold_side_inlet)
    m.fs.cw_in = Arc(source=m.fs.cooling_water_in.outlet, destination=m.fs.condenser.cold_side_inlet)
    m.fs.cw_out = Arc(source=m.fs.condenser.cold_side_outlet, destination=m.fs.cooling_water_out.inlet)
    m.fs.chw_in = Arc(source=m.fs.chilled_water_in.outlet, destination=m.fs.evaporator.hot_side_inlet)
    m.fs.chw_out = Arc(source=m.fs.evaporator.hot_side_outlet, destination=m.fs.chilled_water_out.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Closing the refrigerant loop: the evaporator outlet is the compressor inlet.
    # Only enthalpy and pressure are linked; the circulating flow is set at the
    # compressor inlet (linking the flow too would repeat the loop mass balance).
    suction = m.fs.compressor.control_volume.properties_in[0]
    evap_out = m.fs.evaporator.cold_side.properties_out[0]
    m.fs.loop_enthalpy = Constraint(expr=suction.enth_mol == evap_out.enth_mol)
    m.fs.loop_pressure = Constraint(expr=suction.pressure == evap_out.pressure)

    # Fouling: U = 1 / (1/U_clean + R_f)
    for hx, key in ((m.fs.condenser, "condenser"), (m.fs.evaporator, "evaporator")):
        hx.U_clean = Var(initialize=design[key + "_U"], units=pyunits.W / pyunits.m**2 / pyunits.K)
        hx.fouling_resistance = Var(initialize=0, units=pyunits.m**2 * pyunits.K / pyunits.W)
        hx.fouling_eq = Constraint(
            expr=hx.overall_heat_transfer_coefficient[0] * (1 + hx.U_clean * hx.fouling_resistance) == hx.U_clean)

    # Control specs: condenser subcooling and evaporator superheat
    cond_out = m.fs.condenser.hot_side.properties_out[0]
    m.fs.subcooling = Var(initialize=design["subcooling"], units=pyunits.K)
    m.fs.superheat = Var(initialize=design["superheat"], units=pyunits.K)
    m.fs.subcooling_spec = Constraint(expr=cond_out.temperature == cond_out.temperature_sat - m.fs.subcooling)
    m.fs.superheat_spec = Constraint(expr=evap_out.temperature == evap_out.temperature_sat + m.fs.superheat)

    # Optional control mode: hold the chilled water supply temperature by adjusting
    # the refrigerant flow (compressor capacity). Inactive until set_operating_mode().
    m.fs.chilled_water_supply = Var(initialize=7.0 + 273.15, units=pyunits.K)
    m.fs.chilled_water_supply_spec = Constraint(
        expr=m.fs.evaporator.hot_side.properties_out[0].temperature == m.fs.chilled_water_supply)
    m.fs.chilled_water_supply_spec.deactivate()

    return m


def set_operating_mode(m, mode, chilled_water_supply=7.0, refrigerant_flow=None):
    """Switch between the two ways the chiller is operated.

    * "fixed_flow": the refrigerant flow is fixed (as in the DWSIM model); the
      chilled water supply temperature is a result.
    * "fixed_supply": the chilled water supply temperature is held (C) by
      adjusting the refrigerant flow, as a capacity-controlled chiller does.
    """
    flow = m.fs.compressor.control_volume.properties_in[0].flow_mol
    if mode == "fixed_supply":
        m.fs.chilled_water_supply.fix(chilled_water_supply + 273.15)
        m.fs.chilled_water_supply_spec.activate()
        flow.unfix()
    elif mode == "fixed_flow":
        m.fs.chilled_water_supply_spec.deactivate()
        m.fs.chilled_water_supply.unfix()
        if refrigerant_flow is not None:
            flow.fix(refrigerant_flow / value(m.fs.r134a.mw))
        else:
            flow.fix()
    else:
        raise ValueError("mode must be 'fixed_flow' or 'fixed_supply'")


def fix_design(m, design=DESIGN):
    """Fix the design inputs of the chiller (refrigerant flow mode)."""
    mw_r = value(m.fs.r134a.mw)
    m.fs.compressor.control_volume.properties_in[0].flow_mol.fix(design["refrigerant_flow"] / mw_r)
    m.fs.compressor.efficiency_isentropic.fix(design["compressor_efficiency"])
    m.fs.subcooling.fix(design["subcooling"])
    m.fs.superheat.fix(design["superheat"])
    for hx, key in ((m.fs.condenser, "condenser"), (m.fs.evaporator, "evaporator")):
        hx.U_clean.fix(design[key + "_U"])
        hx.fouling_resistance.fix(0)
        hx.area.fix(design[key + "_area"])
    m.fs.expansion_valve.valve_opening.fix(1)
    m.fs.expansion_valve.Cv.unfix()  # calculated from the flow and pressure drop

    p_w = design["water_pressure"]
    for feed, flow, T in ((m.fs.chilled_water_in, "chilled_water_flow", "chilled_water_return"),
                          (m.fs.cooling_water_in, "cooling_water_flow", "cooling_water_inlet")):
        feed.flow_mol.fix(design[flow] / MW_WATER)
        feed.pressure.fix(p_w)
        feed.enth_mol.fix(value(htpx(T=(design[T] + 273.15) * pyunits.K, P=p_w * pyunits.Pa)))


def initialize_chiller(m, design=DESIGN, outlvl=idaeslog.WARNING):
    """Initialize the cycle unit by unit with the loop open, then close it and solve."""
    solver = get_solver(options={"max_iter": 300})
    suction = m.fs.compressor.control_volume.properties_in[0]
    p_evap, p_cond = design["evaporator_pressure"], design["condenser_pressure"]

    # Open loop: fix the suction state and both operating pressures
    for c in (m.fs.loop_enthalpy, m.fs.loop_pressure, m.fs.subcooling_spec, m.fs.superheat_spec):
        c.deactivate()
    suction.pressure.fix(p_evap)
    suction.enth_mol.fix(value(m.fs.r134a.htpx(p=p_evap * pyunits.Pa, x=1)) + 400)  # slightly superheated vapor
    m.fs.compressor.outlet.pressure.fix(p_cond)
    m.fs.expansion_valve.outlet.pressure.fix(p_evap)

    m.fs.compressor.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    m.fs.cooling_water_in.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.cw_in)
    m.fs.condenser.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s05)
    m.fs.expansion_valve.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s07)
    m.fs.chilled_water_in.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.chw_in)
    m.fs.evaporator.initialize(outlvl=outlvl)
    for arc, product in ((m.fs.cw_out, m.fs.cooling_water_out), (m.fs.chw_out, m.fs.chilled_water_out)):
        propagate_state(arc=arc)
        product.initialize(outlvl=outlvl)
    solver.solve(m)

    # Close the loop and let the pressures float to the subcooling and superheat specs
    suction.pressure.unfix()
    suction.enth_mol.unfix()
    m.fs.compressor.outlet.pressure.unfix()
    m.fs.expansion_valve.outlet.pressure.unfix()
    for c in (m.fs.loop_enthalpy, m.fs.loop_pressure, m.fs.subcooling_spec, m.fs.superheat_spec):
        c.activate()
    return solver.solve(m)


def setup_chiller_model(design=DESIGN, outlvl=idaeslog.WARNING):
    """Build, initialize and solve the chiller at its design point."""
    m = build_chiller(design)
    fix_design(m, design)
    print("The DOF is {0}".format(degrees_of_freedom(m)))
    result = initialize_chiller(m, design, outlvl)
    return m, result


def chiller_results(m):
    """Key chiller results (C, bar, kW)."""
    r134a_mw = value(m.fs.r134a.mw)
    suction = m.fs.compressor.control_volume.properties_in[0]
    discharge = m.fs.compressor.control_volume.properties_out[0]
    cond_out = m.fs.condenser.hot_side.properties_out[0]
    valve_out = m.fs.expansion_valve.control_volume.properties_out[0]
    chw_out = m.fs.evaporator.hot_side.properties_out[0]
    cw_out = m.fs.condenser.cold_side.properties_out[0]
    Q_evap = value(m.fs.evaporator.heat_duty[0]) / 1000
    W = value(m.fs.compressor.work_mechanical[0]) / 1000
    return {
        "Evaporator pressure (bar)": value(suction.pressure) / 1e5,
        "Evaporating temperature (C)": value(suction.temperature_sat) - 273.15,
        "Condenser pressure (bar)": value(discharge.pressure) / 1e5,
        "Condensing temperature (C)": value(discharge.temperature_sat) - 273.15,
        "Suction temperature (C)": value(suction.temperature) - 273.15,
        "Discharge temperature (C)": value(discharge.temperature) - 273.15,
        "Condenser outlet temperature (C)": value(cond_out.temperature) - 273.15,
        "Valve outlet vapor fraction": value(valve_out.vapor_frac),
        "Chilled water supply (C)": value(chw_out.temperature) - 273.15,
        "Cooling water outlet (C)": value(cw_out.temperature) - 273.15,
        "Refrigerant flow (kg/s)": value(suction.flow_mol) * r134a_mw,
        "Cooling duty (kW)": Q_evap,
        "Condenser duty (kW)": value(m.fs.condenser.heat_duty[0]) / 1000,
        "Compressor power (kW)": W,
        "COP": Q_evap / W,
    }


def cycle_points(m):
    """Refrigerant state points of the cycle as (enthalpy kJ/kg, pressure bar), in loop order.

    Order: compressor suction, discharge, condenser outlet, valve outlet, back to suction.
    """
    mw = value(m.fs.r134a.mw)
    states = [
        m.fs.compressor.control_volume.properties_in[0],
        m.fs.compressor.control_volume.properties_out[0],
        m.fs.condenser.hot_side.properties_out[0],
        m.fs.expansion_valve.control_volume.properties_out[0],
    ]
    points = [(value(s.enth_mol) / mw / 1000, value(s.pressure) / 1e5) for s in states]
    return points + [points[0]]


def plot_ph(ax, cycles):
    """Draw the R134a saturation dome and one or more cycles on a pressure-enthalpy diagram.

    cycles: dict of label -> list of (h kJ/kg, P bar) from cycle_points().
    """
    import numpy as np
    from CoolProp.CoolProp import PropsSI

    T = np.linspace(250, PropsSI("Tcrit", "R134a") - 0.5, 200)
    for q in (0, 1):
        ax.plot([PropsSI("H", "T", t, "Q", q, "R134a") / 1000 for t in T],
                [PropsSI("P", "T", t, "Q", q, "R134a") / 1e5 for t in T], color="grey", linewidth=1)
    for label, points in cycles.items():
        h, p = zip(*points)
        ax.plot(h, p, "o-", label=label)
    ax.set_yscale("log")
    ax.set_xlabel("Enthalpy (kJ/kg)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_title("R134a cycle on the P-h diagram")
    ax.grid(alpha=0.3, which="both")
    ax.legend()


def read_dwsim_results(path):
    """Read the chiller results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Streams are found
    by their tag.
    """
    from CoolProp.CoolProp import PropsSI

    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([-+0-9.eE]+)<".format(tag), text).group(1))

    tags = dict(re.findall(r"<Name>(MAT-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        tag = tags.get("MAT-" + stream[:36])
        if tag:
            streams[tag] = stream[stream.find("<ComponentName>Mixture<"):]

    def T(tag):
        return first("temperature", streams[tag]) - 273.15

    def P(tag):
        return first("pressure", streams[tag])

    compressor = xml[xml.find("UnitOperations.Compressor<"):]
    W = first("DeltaQ", compressor)
    Q_evap = first("massflow", streams["Chilled Water In"]) * (
        first("enthalpy", streams["Chilled Water In"]) - first("enthalpy", streams["Chilled Water Out"]))
    Q_cond = first("massflow", streams["Cooling Water In"]) * (
        first("enthalpy", streams["Cooling Water Out"]) - first("enthalpy", streams["Cooling Water In"]))
    return {
        "Evaporator pressure (bar)": P("1") / 1e5,
        "Evaporating temperature (C)": PropsSI("T", "P", P("1"), "Q", 1, "R134a") - 273.15,
        "Condenser pressure (bar)": P("2") / 1e5,
        "Condensing temperature (C)": PropsSI("T", "P", P("2"), "Q", 0, "R134a") - 273.15,
        "Suction temperature (C)": T("1"),
        "Discharge temperature (C)": T("2"),
        "Condenser outlet temperature (C)": T("5"),
        "Valve outlet vapor fraction": PropsSI("Q", "P", P("7"), "H", first("enthalpy", streams["7"]) * 1000, "R134a"),
        "Chilled water supply (C)": T("Chilled Water Out"),
        "Cooling water outlet (C)": T("Cooling Water Out"),
        "Refrigerant flow (kg/s)": first("massflow", streams["1"]),
        "Cooling duty (kW)": Q_evap,
        "Condenser duty (kW)": Q_cond,
        "Compressor power (kW)": W,
        "COP": Q_evap / W,
    }


def compare_with_dwsim(m, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": chiller_results(m), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_chiller_model()
    print(result.solver.termination_condition)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Refrigeraion Chilled Water Skid.dwxmz")))
