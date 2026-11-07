# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product and valve unit models
from idaes.models.unit_models import Feed, Product
from idaes.models.unit_models.valve import Valve, ValveFunctionType

# Import the general Helmholtz equation of state framework. Nitrogen is added to
# it by nitrogen_helmholtz.py, which writes the files in helmholtz_data/
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
    set_parameter_path,
)

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import propagate_state to pass the initialized states between unit models
from idaes.core.util.initialization import propagate_state

# Import Arc to connect the unit models
from pyomo.network import Arc

import os
import re
import zipfile

import pandas as pd

HELMHOLTZ_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "helmholtz_data")


def setup_valve_model(
    flow_mass=1.0,  # kg/s (3600 kg/h)
    pressure_in=5e5,  # Pa
    vapor_frac_in=0.0,  # 0 = saturated liquid
    pressure_drop=0.5e5,  # Pa
    valve_opening=0.5,  # fraction open, as DWSIM 50 %
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a valve throttling saturated liquid nitrogen.

    The inlet temperature (-179.155 C at 5 bar) is the saturation temperature,
    so temperature and pressure alone do not fix the state. The inlet is
    therefore specified by pressure and vapor fraction.
    """
    # Point IDAES at the nitrogen parameter files, which also registers nitrogen
    set_parameter_path(HELMHOLTZ_DATA)

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = HelmholtzParameterBlock(
        pure_component="nitrogen",
        phase_presentation=PhaseType.MIX,
        state_vars=StateVars.PH,
        amount_basis=AmountBasis.MOLE,
    )
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.valve = Valve(
        property_package=m.fs.properties,
        valve_function_callback=ValveFunctionType.linear,
    )
    m.fs.product = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.valve.inlet)
    m.fs.s02 = Arc(source=m.fs.valve.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    mw = value(m.fs.properties.mw)  # kg/mol
    m.fs.feed.flow_mol.fix(flow_mass / mw)  # mol/s
    m.fs.feed.pressure.fix(pressure_in)
    m.fs.feed.enth_mol.fix(
        value(m.fs.properties.htpx(p=pressure_in * pyunits.Pa, x=vapor_frac_in))
    )  # J/mol

    # Pressure drop calculation mode: fix deltaP and the opening, and let the
    # valve flow coefficient Cv be calculated from the flow equation
    m.fs.valve.deltaP.fix(-pressure_drop)
    m.fs.valve.valve_opening.fix(valve_opening)
    m.fs.valve.Cv.unfix()

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    m.fs.valve.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    m.fs.product.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_valve_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    valve_prop = model.fs.valve.report()
    product_prop = model.fs.product.report()

    return feed_prop, valve_prop, product_prop


def valve_results(model):
    """Key valve results (K, Pa, kJ/kg)."""
    v = model.fs.valve.control_volume
    inlet, outlet = v.properties_in[0], v.properties_out[0]
    mw = value(model.fs.properties.mw)
    return {
        "Inlet temperature (K)": value(inlet.temperature),
        "Outlet temperature (K)": value(outlet.temperature),
        "Temperature change (K)": value(outlet.temperature - inlet.temperature),
        "Outlet pressure (Pa)": value(outlet.pressure),
        "Outlet vapor fraction": value(outlet.vapor_frac),
        "Outlet enthalpy (kJ/kg)": value(outlet.enth_mol) / mw / 1000,
    }


def read_dwsim_results(path):
    """Read the valve results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Each stream lists
    its phases (Mixture, Vapor, ...) followed by their properties.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    valve = xml[xml.find("UnitOperations.Valve<"):]
    pressure_out = first("OutletPressure", valve)  # Pa

    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        mixture = stream[stream.find("<ComponentName>Mixture<"):]
        vapor = stream[stream.find("<ComponentName>Vapor<"):]
        name = "outlet" if abs(first("pressure", mixture) - pressure_out) < 1 else "inlet"
        streams[name] = (mixture, vapor)

    (inlet, _), (outlet, outlet_vapor) = streams["inlet"], streams["outlet"]
    return {
        "Inlet temperature (K)": first("temperature", inlet),
        "Outlet temperature (K)": first("temperature", outlet),
        "Temperature change (K)": first("DeltaT", valve),
        "Outlet pressure (Pa)": first("pressure", outlet),
        # Pure component, so the vapor mass fraction equals the molar vapor fraction
        "Outlet vapor fraction": first("massflow", outlet_vapor) / first("massflow", outlet),
        "Outlet enthalpy (kJ/kg)": first("enthalpy", outlet),
    }


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the relative error in %."""
    df = pd.DataFrame({"IDAES": valve_results(model), "DWSIM": dwsim})
    df["Error (%)"] = (df["IDAES"] - df["DWSIM"]) / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_valve_model()
    print(result.solver.termination_condition)
    report_valve_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Valve Model.dwxmz")))
