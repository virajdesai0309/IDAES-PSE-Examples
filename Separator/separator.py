# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, separator and product unit models
from idaes.models.unit_models import Feed, Separator, Product

# Import the option to set the basis for splitting
from idaes.models.unit_models.separator import SplittingType

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for the flowsheet
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the gas property package to create a configuration file for the GenericParameterBlock
from gas import configuration

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

OUTLETS = ["outlet_1", "outlet_2"]


def setup_separator_model(
    flow_mass=1.0,  # kg/s (3600 kg/h)
    pressure=101325,  # Pa
    temperature=25.0,  # C
    split_fraction=0.75,  # fraction of the inlet sent to outlet_1
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a splitter model for a pure methane stream."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.separator = Separator(
        property_package=m.fs.properties,
        # Split the total flow, as DWSIM "Split Ratios" mode
        split_basis=SplittingType.totalFlow,
        outlet_list=OUTLETS,
        ideal_separation=False,
        has_phase_equilibrium=False,
    )
    m.fs.product_1 = Product(property_package=m.fs.properties)
    m.fs.product_2 = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.separator.inlet)
    m.fs.s02 = Arc(source=m.fs.separator.outlet_1, destination=m.fs.product_1.inlet)
    m.fs.s03 = Arc(source=m.fs.separator.outlet_2, destination=m.fs.product_2.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    mw = value(m.fs.properties.methane.mw)  # kg/mol
    m.fs.feed.flow_mol.fix(flow_mass / mw)  # mol/s
    m.fs.feed.mole_frac_comp[0, "methane"].fix(1.0)
    m.fs.feed.pressure.fix(pressure)
    m.fs.feed.temperature.fix(temperature + 273.15)  # K

    m.fs.separator.split_fraction[0, "outlet_1"].fix(split_fraction)

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    m.fs.separator.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    propagate_state(arc=m.fs.s03)
    m.fs.product_1.initialize(outlvl=outlvl)
    m.fs.product_2.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_separator_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    separator_prop = model.fs.separator.report()
    product_1_prop = model.fs.product_1.report()
    product_2_prop = model.fs.product_2.report()

    return feed_prop, separator_prop, product_1_prop, product_2_prop


def separator_results(model):
    """Key outlet results for both outlets (K, Pa, kg/h, mol/s)."""
    results = {}
    for o in OUTLETS:
        state = getattr(model.fs.separator, o + "_state")[0]
        results[o + " temperature (K)"] = value(state.temperature)
        results[o + " pressure (Pa)"] = value(state.pressure)
        results[o + " mass flow (kg/h)"] = value(state.flow_mol * state.mw) * 3600
        results[o + " molar flow (mol/s)"] = value(state.flow_mol)
    return results


def read_dwsim_results(path):
    """Read the splitter outlet streams from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. The inlet is the
    material stream with the largest mass flow; the outlets follow in order of
    decreasing mass flow (outlet_1 is the larger split).
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    mixtures = [
        stream[stream.find("<ComponentName>Mixture<"):]
        for stream in xml.split("<ComponentName>MAT-")[1:]
    ]
    mixtures.sort(key=lambda s: first("massflow", s), reverse=True)

    results = {}
    for o, mixture in zip(OUTLETS, mixtures[1:]):
        results[o + " temperature (K)"] = first("temperature", mixture)
        results[o + " pressure (Pa)"] = first("pressure", mixture)
        results[o + " mass flow (kg/h)"] = first("massflow", mixture) * 3600
        results[o + " molar flow (mol/s)"] = first("molarflow", mixture)
    return results


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the relative error in %."""
    df = pd.DataFrame({"IDAES": separator_results(model), "DWSIM": dwsim})
    df["Error (%)"] = (df["IDAES"] - df["DWSIM"]) / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_separator_model()
    print(result.solver.termination_condition)
    report_separator_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Splitter Model.dwxmz")))
