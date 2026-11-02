# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product and CSTR unit models
from idaes.models.unit_models import Feed, Product, CSTR

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property and reaction packages
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock

# Import the methanol carbonylation property and reaction packages
from methanol_carbonylation import configuration as property_configuration
from carbonylation_reaction import configuration as reaction_configuration

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

COMPONENTS = ["methanol", "water", "acetic_acid", "carbon_monoxide"]

# Feed of the DWSIM model, component molar flows in mol/s (1 kg/s in total)
DEFAULT_FEED = {
    "flow_mol_comp": {"methanol": 20.4600949, "water": 18.4691703, "acetic_acid": 0, "carbon_monoxide": 0.41746463},
    "temperature": 328,  # K
    "pressure": 101325,  # Pa
}


def setup_cstr_model(feed=DEFAULT_FEED, volume=100, outlvl=idaeslog.WARNING):
    """Creates and solves an adiabatic CSTR for methanol carbonylation."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**property_configuration)
    m.fs.reactions = GenericReactionParameterBlock(
        property_package=m.fs.properties, **reaction_configuration
    )

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.cstr = CSTR(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_equilibrium_reactions=False,
        # Heats of formation are in the enthalpies, so no separate heat of reaction term
        has_heat_of_reaction=False,
        # Adiabatic reactor, as the DWSIM model
        has_heat_transfer=False,
        has_pressure_change=False,
        has_phase_equilibrium=True,
    )
    m.fs.product = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.cstr.inlet)
    m.fs.s02 = Arc(source=m.fs.cstr.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    flow_mol = sum(feed["flow_mol_comp"].values())
    m.fs.feed.flow_mol.fix(flow_mol)  # mol/s
    m.fs.feed.temperature.fix(feed["temperature"])  # K
    m.fs.feed.pressure.fix(feed["pressure"])  # Pa
    for j in COMPONENTS:
        # Small non-zero values avoid log(0) in the phase equilibrium expressions
        m.fs.feed.mole_frac_comp[0, j].fix(max(feed["flow_mol_comp"][j], 1e-8) / flow_mol)

    m.fs.cstr.volume.fix(volume)  # m3

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    m.fs.cstr.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    m.fs.product.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_cstr_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    cstr_prop = model.fs.cstr.report()
    product_prop = model.fs.product.report()

    return feed_prop, cstr_prop, product_prop


def cstr_results(model):
    """Key CSTR results (K, mol/s)."""
    cv = model.fs.cstr.control_volume
    inlet, outlet = cv.properties_in[0], cv.properties_out[0]
    results = {
        "Feed vapor flow (mol/s)": value(inlet.flow_mol_phase["Vap"]),
        "Product temperature (K)": value(outlet.temperature),
        "Temperature rise (K)": value(outlet.temperature - inlet.temperature),
        "Product molar flow (mol/s)": value(outlet.flow_mol),
    }
    for j in COMPONENTS:
        results["Product " + j + " (mol/s)"] = value(outlet.flow_mol_comp[j])
    results["CO conversion"] = 1 - value(outlet.flow_mol_comp["carbon_monoxide"] / inlet.flow_mol_comp["carbon_monoxide"])
    return results


def read_dwsim_results(path):
    """Read the feed and product results from a solved DWSIM .dwxmz file.

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
        "methanol": "Methanol",
        "water": "Water",
        "acetic_acid": "Acetic acid",
        "carbon_monoxide": "Carbon monoxide",
    }

    # Stream tags (names shown on the flowsheet) are stored with the graphic objects
    tags = dict(re.findall(r"<Name>(MAT-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        streams[tags["MAT-" + stream[:36]]] = stream

    def phase(stream, name):
        return stream[stream.find("<ComponentName>{0}<".format(name)):]

    def component_flow(stream, name):
        flows = re.findall(
            r"<ComponentName>{0}</ComponentName>.*?<MolarFlow>([^<]+)<".format(name), stream, re.S
        )
        return next(float(f) for f in flows if f != "NaN")

    feed, product = streams["feed"], streams["product"]
    feed_vapor = re.search(r"<molarflow>([^<]+)<", phase(feed, "Vapor"))
    results = {
        "Feed vapor flow (mol/s)": float(feed_vapor.group(1)) if feed_vapor else 0.0,
        "Product temperature (K)": first("temperature", phase(product, "Mixture")),
        "Temperature rise (K)": first("temperature", phase(product, "Mixture"))
        - first("temperature", phase(feed, "Mixture")),
        "Product molar flow (mol/s)": first("molarflow", phase(product, "Mixture")),
    }
    for j, name in names.items():
        results["Product " + j + " (mol/s)"] = component_flow(product, name)
    results["CO conversion"] = 1 - component_flow(product, names["carbon_monoxide"]) / component_flow(
        feed, names["carbon_monoxide"])
    return results


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": cstr_results(model), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"].where(df["DWSIM"].abs() > 1e-6) * 100
    return df


if __name__ == "__main__":
    model, result = setup_cstr_model()
    print(result.solver.termination_condition)
    report_cstr_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim CSTR Model.dwxmz")))
