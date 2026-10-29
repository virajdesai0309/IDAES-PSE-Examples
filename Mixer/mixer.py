# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product and mixer unit models
from idaes.models.unit_models import Feed, Product, Mixer, MomentumMixingType

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for the flowsheet
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the light hydrocarbon property package to create a configuration file for the GenericParameterBlock
from light_hydrocarbons import configuration

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import propagate_state to pass the initialized feed states to the mixer inlets
from idaes.core.util.initialization import propagate_state

# Import Arc to connect the unit models
from pyomo.network import Arc

import re
import zipfile

import pandas as pd

COMPONENTS = ["methane", "ethane", "propane"]

# Each feed is a pure component stream
DEFAULT_STREAMS = {
    "methane": {"flow_mass": 1.0, "temperature": 80.0, "pressure": 101325},  # kg/s, C, Pa
    "ethane": {"flow_mass": 1.0, "temperature": 50.0, "pressure": 101325},
    "propane": {"flow_mass": 1.0, "temperature": 25.0, "pressure": 101325},
}


def setup_mixer_model(streams=DEFAULT_STREAMS, outlvl=idaeslog.WARNING):
    """Creates and solves a mixer model with one pure feed per component."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=[c + "_inlet" for c in COMPONENTS],
        # Outlet pressure = minimum inlet pressure, as DWSIM "Minimum" pressure calculation
        momentum_mixing_type=MomentumMixingType.minimize,
    )
    m.fs.product = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    for c in COMPONENTS:
        feed = Feed(property_package=m.fs.properties)
        m.fs.add_component(c + "_feed", feed)
        m.fs.add_component(
            c + "_arc", Arc(source=feed.outlet, destination=getattr(m.fs.mixer, c + "_inlet"))
        )
    m.fs.outlet_arc = Arc(source=m.fs.mixer.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    for c in COMPONENTS:
        feed = getattr(m.fs, c + "_feed")
        mw = value(m.fs.properties.get_component(c).mw)  # kg/mol
        feed.flow_mol.fix(streams[c]["flow_mass"] / mw)  # mol/s
        feed.temperature.fix(streams[c]["temperature"] + 273.15)  # K
        feed.pressure.fix(streams[c]["pressure"])  # Pa
        for j in COMPONENTS:
            # A small non-zero value avoids log(0) in the Peng-Robinson fugacity terms
            feed.mole_frac_comp[0, j].fix(1.0 if j == c else 1e-8)

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    for c in COMPONENTS:
        getattr(m.fs, c + "_feed").initialize(outlvl=outlvl)
        propagate_state(arc=getattr(m.fs, c + "_arc"))
    m.fs.mixer.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.outlet_arc)
    m.fs.product.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_mixer_properties(model):
    """Report properties of the solved model."""
    feed_props = [getattr(model.fs, c + "_feed").report() for c in COMPONENTS]
    mixer_prop = model.fs.mixer.report()
    product_prop = model.fs.product.report()

    return feed_props, mixer_prop, product_prop


def mixer_results(model):
    """Key mixer outlet results (K, Pa, kg/h, mol/s)."""
    out = model.fs.mixer.mixed_state[0]
    results = {
        "Outlet temperature (K)": value(out.temperature),
        "Outlet pressure (Pa)": value(out.pressure),
        "Outlet mass flow (kg/h)": value(out.flow_mol * out.mw) * 3600,
        "Outlet molar flow (mol/s)": value(out.flow_mol),
    }
    for c in COMPONENTS:
        results["Mole fraction " + c] = value(out.mole_frac_comp[c])
    return results


def read_dwsim_results(path):
    """Read the mixer outlet stream from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. The outlet is the
    material stream with the largest mass flow.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    streams = []
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        mixture = stream[stream.find("<ComponentName>Mixture<"):]
        streams.append((first("massflow", mixture), stream, mixture))
    _, outlet, mixture = max(streams, key=lambda s: s[0])

    results = {
        "Outlet temperature (K)": first("temperature", mixture),
        "Outlet pressure (Pa)": first("pressure", mixture),
        "Outlet mass flow (kg/h)": first("massflow", mixture) * 3600,
        "Outlet molar flow (mol/s)": first("molarflow", mixture),
    }
    # Overall composition: the first numeric MoleFraction listed for each compound
    for c in COMPONENTS:
        fractions = re.findall(
            r"<ComponentName>{0}<.*?<MoleFraction>([^<]+)<".format(c.capitalize()),
            outlet,
            re.S,
        )
        results["Mole fraction " + c] = next(float(f) for f in fractions if f != "NaN")
    return results


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the relative error in %."""
    df = pd.DataFrame({"IDAES": mixer_results(model), "DWSIM": dwsim})
    df["Error (%)"] = (df["IDAES"] - df["DWSIM"]) / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_mixer_model()
    print(result.solver.termination_condition)
    report_mixer_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Mixer Model.dwxmz")))
