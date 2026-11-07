# Import objects from pyomo package
from pyomo.environ import ConcreteModel, Constraint, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product and plug flow reactor unit models
from idaes.models.unit_models import Feed, Product, PFR

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property and reaction packages
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock

# Import the reformer gas property package and the steam reforming reaction package
from reforming_gas import configuration as gas_configuration, viscosity_vapor
from reforming_reaction import configuration as reaction_configuration

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

COMPONENTS = ["methane", "water", "hydrogen", "carbon_monoxide", "carbon_dioxide"]

DEFAULT_FEED = {
    "flow_mass": 0.0833,  # kg/s
    "temperature": 973.15,  # K
    "pressure": 101325,  # Pa
    "mole_frac": {"methane": 0.4975, "water": 0.4975, "hydrogen": 0.005},
}

DEFAULT_BED = {
    "length": 1.0,  # m
    "area": 1.0,  # m2 (1 m3 reactor volume)
    "void_fraction": 0.4,
    "particle_diameter": 0.002,  # m
    "sphericity": 1.0,
}


def setup_pfr_model(feed=DEFAULT_FEED, bed=DEFAULT_BED, finite_elements=20, collocation_points=3, outlvl=idaeslog.WARNING):
    """Creates and solves an isothermal packed bed PFR for steam methane reforming."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**gas_configuration)
    m.fs.reactions = GenericReactionParameterBlock(
        property_package=m.fs.properties, **reaction_configuration
    )

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.pfr = PFR(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_equilibrium_reactions=False,
        # Heats of formation are in the enthalpies, so no separate heat of reaction term
        has_heat_of_reaction=False,
        # Heat duty needed to keep the reactor isothermal, as DWSIM "Isothermic"
        has_heat_transfer=True,
        # Packed bed pressure drop (Ergun equation)
        has_pressure_change=True,
        # Orthogonal collocation resolves the steep rate near the inlet on a coarse mesh
        transformation_method="dae.collocation",
        transformation_scheme="LAGRANGE-RADAU",
        finite_elements=finite_elements,
        collocation_points=collocation_points,
    )
    m.fs.product = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pfr.inlet)
    m.fs.s02 = Arc(source=m.fs.pfr.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    cv = m.fs.pfr.control_volume
    x_domain = cv.length_domain
    x_in = x_domain.first()

    @m.fs.pfr.Constraint(m.fs.time, x_domain)
    def isothermal(b, t, x):
        if x == x_in:
            return Constraint.Skip
        return cv.properties[t, x].temperature == cv.properties[t, x_in].temperature

    # Ergun equation for the pressure drop per unit length of the packed bed
    eps = bed["void_fraction"]
    dp = bed["particle_diameter"] * bed["sphericity"]  # effective particle diameter, m

    @m.fs.pfr.Constraint(m.fs.time, x_domain)
    def ergun(b, t, x):
        state = cv.properties[t, x]
        u = state.flow_vol / m.fs.pfr.area  # superficial velocity
        mu = viscosity_vapor(state)
        return cv.deltaP[t, x] == -(
            150 * mu * u * (1 - eps) ** 2 / (eps**3 * (dp * pyunits.m) ** 2)
            + 1.75 * state.dens_mass_phase["Vap"] * u**2 * (1 - eps) / (eps**3 * dp * pyunits.m)
        )

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    mw = sum(
        feed["mole_frac"].get(j, 0) * value(m.fs.properties.get_component(j).mw) for j in COMPONENTS
    )  # kg/mol
    m.fs.feed.flow_mol.fix(feed["flow_mass"] / mw)  # mol/s
    m.fs.feed.temperature.fix(feed["temperature"])  # K
    m.fs.feed.pressure.fix(feed["pressure"])  # Pa
    for j in COMPONENTS:
        # Small non-zero values avoid log(0) in the Peng-Robinson expressions
        m.fs.feed.mole_frac_comp[0, j].fix(max(feed["mole_frac"].get(j, 0), 1e-10))

    m.fs.pfr.length.fix(bed["length"])  # m
    m.fs.pfr.area.fix(bed["area"])  # m2

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    m.fs.pfr.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    m.fs.product.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_pfr_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    pfr_prop = model.fs.pfr.report()
    product_prop = model.fs.product.report()

    return feed_prop, pfr_prop, product_prop


def heat_duty(model):
    """Total heat duty of the reactor [W] from the overall energy balance, H_out - H_in."""
    cv = model.fs.pfr.control_volume
    inlet = cv.properties[0, cv.length_domain.first()]
    outlet = cv.properties[0, cv.length_domain.last()]
    return value(outlet.flow_mol * outlet.enth_mol - inlet.flow_mol * inlet.enth_mol)


def pfr_results(model):
    """Key PFR results (Pa, kW, mol/s)."""
    cv = model.fs.pfr.control_volume
    x_in, x_out = cv.length_domain.first(), cv.length_domain.last()
    inlet, outlet = cv.properties[0, x_in], cv.properties[0, x_out]
    results = {
        "Outlet temperature (K)": value(outlet.temperature),
        "Outlet pressure (Pa)": value(outlet.pressure),
        "Pressure drop (Pa)": value(inlet.pressure - outlet.pressure),
        "Heat duty (kW)": heat_duty(model) / 1000,
        "Outlet molar flow (mol/s)": value(outlet.flow_mol),
    }
    for j in COMPONENTS:
        results["Outlet mole fraction " + j] = value(outlet.mole_frac_comp[j])
    results["CH4 conversion"] = 1 - value(outlet.flow_mol_comp["methane"] / inlet.flow_mol_comp["methane"])
    return results


def profiles(model):
    """Mole fractions and pressure along the reactor length."""
    cv = model.fs.pfr.control_volume
    L = value(model.fs.pfr.length)
    rows = []
    for x in cv.length_domain:
        state = cv.properties[0, x]
        row = {"Position (m)": x * L, "Pressure (Pa)": value(state.pressure)}
        for j in COMPONENTS:
            row[j] = value(state.mole_frac_comp[j])
        rows.append(row)
    return pd.DataFrame(rows)


def read_dwsim_results(path):
    """Read the PFR results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Streams are found
    by their tag. Each compound is listed once per phase; the first numeric
    mole fraction is the overall stream value.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    names = {
        "methane": "Methane",
        "water": "Water",
        "hydrogen": "Hydrogen",
        "carbon_monoxide": "Carbon monoxide",
        "carbon_dioxide": "Carbon dioxide",
    }

    # Stream tags (names shown on the flowsheet) are stored with the graphic objects
    tags = dict(re.findall(r"<Name>(MAT-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        streams[tags["MAT-" + stream[:36]]] = stream

    def mixture(stream):
        return stream[stream.find("<ComponentName>Mixture<"):]

    def component(stream, name, prop):
        values = re.findall(
            r"<ComponentName>{0}</ComponentName>.*?<{1}>([^<]+)<".format(name, prop), stream, re.S
        )
        return next(float(v) for v in values if v != "NaN")

    reactor = xml[xml.find("Reactors.Reactor_PFR<"):]
    feed, product = streams["feed"], streams["product"]
    results = {
        "Outlet temperature (K)": first("temperature", mixture(product)),
        "Outlet pressure (Pa)": first("pressure", mixture(product)),
        "Pressure drop (Pa)": first("DeltaP", reactor),
        "Heat duty (kW)": first("DeltaQ", reactor),
        "Outlet molar flow (mol/s)": first("molarflow", mixture(product)),
    }
    for j, name in names.items():
        results["Outlet mole fraction " + j] = component(product, name, "MoleFraction")
    results["CH4 conversion"] = 1 - component(product, "Methane", "MolarFlow") / component(feed, "Methane", "MolarFlow")
    return results


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": pfr_results(model), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"].where(df["DWSIM"].abs() > 1e-9) * 100
    return df


if __name__ == "__main__":
    model, result = setup_pfr_model()
    print(result.solver.termination_condition)
    report_pfr_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Plug Flow Reactor.dwxmz")))
