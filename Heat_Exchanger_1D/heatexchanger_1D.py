# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product and 1D heat exchanger unit models
from idaes.models.unit_models import Feed, Product
from idaes.models.unit_models.heat_exchanger_1D import (HeatExchanger1D,
                                                        HeatExchangerFlowPattern)

# Import the IAPWS-95 steam table property package for the cooling water
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for ethanol
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the ethanol property package to create a configuration file for the GenericParameterBlock
from ethanol import configuration

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

MW_ETHANOL = 0.046068  # kg/mol
MW_WATER = 0.01801528  # kg/mol

DEFAULT_HOT = {"flow_mass": 1.0, "temperature": 75.0, "pressure": 101325}  # ethanol: kg/s, C, Pa
DEFAULT_COLD = {"flow_mass": 1.0, "temperature": 25.0, "pressure": 101325}  # water: kg/s, C, Pa


def setup_heatexchanger_1d_model(
    hot=DEFAULT_HOT,
    cold=DEFAULT_COLD,
    area=1.0,  # m2, total heat transfer area
    length=1.0,  # m, only sets the length axis of the 1D model
    u=1000.0,  # W/m2/K, overall heat transfer coefficient
    finite_elements=10,
    collocation_points=3,
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a countercurrent 1D heat exchanger.

    Hot ethanol flows in the tube and cold water in the shell. The local heat
    flow is U * A / L * (T_hot - T_cold), so the result depends on U * A only.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.ethanol_properties = GenericParameterBlock(**configuration)
    m.fs.water_properties = iapws95.Iapws95ParameterBlock()

    m.fs.ethanol_inlet = Feed(property_package=m.fs.ethanol_properties)
    m.fs.ethanol_outlet = Product(property_package=m.fs.ethanol_properties)
    m.fs.water_inlet = Feed(property_package=m.fs.water_properties)
    m.fs.water_outlet = Product(property_package=m.fs.water_properties)

    # Orthogonal collocation keeps the discretized energy balances accurate on a coarse mesh
    discretization = {
        "transformation_method": "dae.collocation",
        "transformation_scheme": "LAGRANGE-RADAU",
    }
    m.fs.heat_exchanger = HeatExchanger1D(
        hot_side_name="tube",
        cold_side_name="shell",
        tube={"property_package": m.fs.ethanol_properties, **discretization},
        shell={"property_package": m.fs.water_properties, **discretization},
        finite_elements=finite_elements,
        collocation_points=collocation_points,
        flow_type=HeatExchangerFlowPattern.countercurrent,
    )

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.ethanol_inlet.outlet, destination=m.fs.heat_exchanger.hot_side_inlet)
    m.fs.s02 = Arc(source=m.fs.water_inlet.outlet, destination=m.fs.heat_exchanger.cold_side_inlet)
    m.fs.s03 = Arc(source=m.fs.heat_exchanger.hot_side_outlet, destination=m.fs.ethanol_outlet.inlet)
    m.fs.s04 = Arc(source=m.fs.heat_exchanger.cold_side_outlet, destination=m.fs.water_outlet.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.ethanol_inlet.flow_mol.fix(hot["flow_mass"] / MW_ETHANOL)  # mol/s
    m.fs.ethanol_inlet.mole_frac_comp[0, "ethanol"].fix(1.0)
    m.fs.ethanol_inlet.temperature.fix(hot["temperature"] + 273.15)  # K
    m.fs.ethanol_inlet.pressure.fix(hot["pressure"])  # Pa

    m.fs.water_inlet.flow_mol.fix(cold["flow_mass"] / MW_WATER)  # mol/s
    m.fs.water_inlet.pressure.fix(cold["pressure"])  # Pa
    m.fs.water_inlet.enth_mol.fix(
        value(htpx(T=(cold["temperature"] + 273.15) * pyunits.K, P=cold["pressure"] * pyunits.Pa))
    )  # J/mol

    m.fs.heat_exchanger.area.fix(area)
    m.fs.heat_exchanger.length.fix(length)
    m.fs.heat_exchanger.heat_transfer_coefficient.fix(u)

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.ethanol_inlet.initialize(outlvl=outlvl)
    m.fs.water_inlet.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    propagate_state(arc=m.fs.s02)
    m.fs.heat_exchanger.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s03)
    propagate_state(arc=m.fs.s04)
    m.fs.ethanol_outlet.initialize(outlvl=outlvl)
    m.fs.water_outlet.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_heatexchanger_properties(model):
    """Report properties of the solved model."""
    ethanol_inlet_prop = model.fs.ethanol_inlet.report()
    water_inlet_prop = model.fs.water_inlet.report()
    heat_exchanger_prop = model.fs.heat_exchanger.report()
    ethanol_outlet_prop = model.fs.ethanol_outlet.report()
    water_outlet_prop = model.fs.water_outlet.report()

    return ethanol_inlet_prop, water_inlet_prop, heat_exchanger_prop, ethanol_outlet_prop, water_outlet_prop


def heatexchanger_results(model):
    """Key heat exchanger results (K, kW)."""
    hx = model.fs.heat_exchanger
    hot_in = hx.hot_side_inlet
    hot_out = hx.hot_side_outlet
    duty = value(hx.hot_side.properties[0, 0].flow_mol) * (
        value(hx.hot_side.properties[0, 0].enth_mol)
        - value(hx.hot_side.properties[0, 1].enth_mol)
    )  # W
    return {
        "Hot (ethanol) outlet temperature (K)": value(hx.hot_side.properties[0, 1].temperature),
        "Cold (water) outlet temperature (K)": value(hx.cold_side.properties[0, 0].temperature),
        "Duty (kW)": duty / 1000,
    }


def temperature_profiles(model):
    """Hot and cold side temperatures along the exchanger length (C)."""
    hx = model.fs.heat_exchanger
    length = value(hx.length)
    return pd.DataFrame(
        [
            {
                "Position (m)": x * length,
                "Hot (ethanol) temperature (C)": value(hx.hot_side.properties[0, x].temperature) - 273.15,
                "Cold (water) temperature (C)": value(hx.cold_side.properties[0, x].temperature) - 273.15,
            }
            for x in hx.hot_side.length_domain
        ]
    )


def read_dwsim_results(path):
    """Read the heat exchanger results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. The values are read
    from the heat exchanger unit itself, not from the outlet streams.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    hx = xml[xml.find("UnitOperations.HeatExchanger<"):]
    return {
        "Hot (ethanol) outlet temperature (K)": first("HotSideOutletTemperature", hx),
        "Cold (water) outlet temperature (K)": first("ColdSideOutletTemperature", hx),
        "Duty (kW)": first("Q", hx),
    }


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the relative error in %."""
    df = pd.DataFrame({"IDAES": heatexchanger_results(model), "DWSIM": dwsim})
    df["Error (%)"] = (df["IDAES"] - df["DWSIM"]) / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_heatexchanger_1d_model()
    print(result.solver.termination_condition)
    report_heatexchanger_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Heat Exchanger Model.dwxmz")))
