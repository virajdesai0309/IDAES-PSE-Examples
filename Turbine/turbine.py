# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product and turbine unit models
from idaes.models.unit_models import Feed, Product
from idaes.models.unit_models.pressure_changer import Turbine

# Import the IAPWS-95 steam table property package and the htpx helper
# htpx returns the molar enthalpy of water from (T, P) or (P, vapor fraction x)
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import propagate_state to pass the initialized feed state to the turbine inlet
from idaes.core.util.initialization import propagate_state

# Import Arc to connect the unit models
from pyomo.network import Arc

import re
import zipfile

import pandas as pd

MW_WATER = 0.01801528  # kg/mol


def setup_turbine_model(
    flow_mass=10,  # kg/s (36000 kg/h)
    pressure_in=40e5,  # Pa
    vapor_frac_in=1.0,  # 1 = saturated steam, 0 = saturated liquid
    pressure_out=15e5,  # Pa
    efficiency=0.75,
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a steam turbine model.

    The inlet temperature (250.36 C at 40 bar) is the saturation temperature,
    so temperature and pressure alone do not fix the state. The inlet is
    therefore specified by pressure and vapor fraction through htpx.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.turbine = Turbine(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.turbine.inlet)
    m.fs.s02 = Arc(source=m.fs.turbine.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.flow_mol.fix(flow_mass / MW_WATER)  # mol/s
    m.fs.feed.pressure.fix(pressure_in)
    m.fs.feed.enth_mol.fix(
        value(htpx(P=pressure_in * pyunits.Pa, x=vapor_frac_in))
    )  # J/mol

    m.fs.turbine.efficiency_isentropic.fix(efficiency)
    m.fs.turbine.outlet.pressure.fix(pressure_out)

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    m.fs.turbine.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    m.fs.product.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_turbine_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    turbine_prop = model.fs.turbine.report()
    product_prop = model.fs.product.report()

    return feed_prop, turbine_prop, product_prop


def turbine_results(model):
    """Key turbine results in DWSIM units (C, kJ/kg, kW)."""
    t = model.fs.turbine
    out = t.control_volume.properties_out[0]
    return {
        "Outlet temperature (C)": value(out.temperature) - 273.15,
        "Outlet vapor fraction": value(out.vapor_frac),
        "Outlet enthalpy (kJ/kg)": value(out.enth_mol) / MW_WATER / 1000,
        "Power (kW)": -value(t.work_mechanical[0]) / 1000,
    }


def read_dwsim_results(path):
    """Read the turbine results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Each stream lists
    its phases (Mixture, OverallLiquid, Vapor, ...) followed by their properties.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    expander = xml[xml.find("UnitOperations.Expander"):]
    power = first("DeltaQ", expander)  # kW
    pressure_out = first("POut", expander)  # Pa

    # The outlet stream is the material stream at the turbine outlet pressure
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        mixture = stream[stream.find("<ComponentName>Mixture<"):]
        if abs(first("pressure", mixture) - pressure_out) > 1:
            continue
        vapor = stream[stream.find("<ComponentName>Vapor<"):]
        return {
            "Outlet temperature (C)": first("temperature", mixture) - 273.15,
            "Outlet vapor fraction": first("massflow", vapor) / first("massflow", mixture),
            "Outlet enthalpy (kJ/kg)": first("enthalpy", mixture),
            "Power (kW)": power,
        }
    raise ValueError("Turbine outlet stream not found in {0}".format(path))


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the relative error in %."""
    df = pd.DataFrame({"IDAES": turbine_results(model), "DWSIM": dwsim})
    df["Error (%)"] = (df["IDAES"] - df["DWSIM"]) / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_turbine_model()
    print(result.solver.termination_condition)
    report_turbine_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Turbine Model.dwxmz")))
