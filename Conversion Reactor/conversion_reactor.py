# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product, mixer and stoichiometric (conversion) reactor unit models
from idaes.models.unit_models import (Feed, Product, Mixer, MomentumMixingType,
                                      StoichiometricReactor)

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property and reaction packages
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock

# Import the combustion gas property package and the combustion reaction package
from combustion_gas import configuration as gas_configuration
from combustion_reaction import configuration as reaction_configuration

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

COMPONENTS = ["methane", "oxygen", "nitrogen", "carbon_dioxide", "water"]
STOICHIOMETRY = {"methane": -1, "oxygen": -2, "carbon_dioxide": 1, "water": 2}

DEFAULT_FUEL = {"flow_mass": 1.0, "temperature": 300, "pressure": 101325, "mole_frac": {"methane": 1.0}}
DEFAULT_AIR = {"flow_mass": 5.0, "temperature": 300, "pressure": 101325, "mole_frac": {"nitrogen": 0.79, "oxygen": 0.21}}


def limiting_extent(feed_flows, conversion, base="methane"):
    """Extent of reaction for a base component conversion, limited by the reactant that runs out first.

    Args:
        feed_flows: dict of component molar flows entering the reactor (mol/s)
        conversion: specified fractional conversion of the base component
    """
    extents = {base: conversion * feed_flows[base]}
    for j, nu in STOICHIOMETRY.items():
        if nu < 0 and j != base:
            extents[j] = feed_flows.get(j, 0) / -nu
    limiting = min(extents, key=extents.get)
    return extents[limiting], limiting


def setup_conversion_reactor_model(
    fuel=DEFAULT_FUEL,
    air=DEFAULT_AIR,
    conversion=0.95,  # specified conversion of methane, the base component
    outlvl=idaeslog.WARNING,
):
    """Creates and solves an adiabatic methane combustor (conversion reactor)."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**gas_configuration)
    m.fs.reactions = GenericReactionParameterBlock(
        property_package=m.fs.properties, **reaction_configuration
    )

    m.fs.fuel = Feed(property_package=m.fs.properties)
    m.fs.air = Feed(property_package=m.fs.properties)
    m.fs.mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=["fuel_inlet", "air_inlet"],
        momentum_mixing_type=MomentumMixingType.minimize,
    )
    m.fs.combustor = StoichiometricReactor(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        # Heats of formation are in the enthalpies, so no separate heat of reaction term
        has_heat_of_reaction=False,
        # Adiabatic combustor, as DWSIM "Adiabatic" operation mode
        has_heat_transfer=False,
        has_pressure_change=False,
    )
    m.fs.flue_gas = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.fuel.outlet, destination=m.fs.mixer.fuel_inlet)
    m.fs.s02 = Arc(source=m.fs.air.outlet, destination=m.fs.mixer.air_inlet)
    m.fs.s03 = Arc(source=m.fs.mixer.outlet, destination=m.fs.combustor.inlet)
    m.fs.s04 = Arc(source=m.fs.combustor.outlet, destination=m.fs.flue_gas.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    feed_flows = {}
    for feed, spec in ((m.fs.fuel, fuel), (m.fs.air, air)):
        mw = sum(
            spec["mole_frac"].get(j, 0) * value(m.fs.properties.get_component(j).mw)
            for j in COMPONENTS
        )  # kg/mol
        flow_mol = spec["flow_mass"] / mw  # mol/s
        feed.flow_mol.fix(flow_mol)
        feed.temperature.fix(spec["temperature"])  # K
        feed.pressure.fix(spec["pressure"])  # Pa
        for j in COMPONENTS:
            # Small non-zero values avoid log(0) in the Peng-Robinson expressions
            x = max(spec["mole_frac"].get(j, 0), 1e-10)
            feed.mole_frac_comp[0, j].fix(x)
            feed_flows[j] = feed_flows.get(j, 0) + flow_mol * x

    # The conversion reactor needs the extent of reaction. The specified methane
    # conversion is capped by the reactant that runs out first (oxygen here).
    extent, limiting = limiting_extent(feed_flows, conversion)
    m.fs.combustor.rate_reaction_extent[0, "combustion"].fix(extent)
    print("Limiting reactant: {0}, extent of reaction: {1:.4f} mol/s, methane conversion: {2:.2%}".format(
        limiting, extent, extent / feed_flows["methane"]))

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.fuel.initialize(outlvl=outlvl)
    m.fs.air.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    propagate_state(arc=m.fs.s02)
    m.fs.mixer.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s03)
    m.fs.combustor.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s04)
    m.fs.flue_gas.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_conversion_reactor_properties(model):
    """Report properties of the solved model."""
    fuel_prop = model.fs.fuel.report()
    air_prop = model.fs.air.report()
    mixer_prop = model.fs.mixer.report()
    combustor_prop = model.fs.combustor.report()
    flue_gas_prop = model.fs.flue_gas.report()

    return fuel_prop, air_prop, mixer_prop, combustor_prop, flue_gas_prop


def reactor_results(model):
    """Key combustor results (K, kg/s, mol/s)."""
    cv = model.fs.combustor.control_volume
    inlet, outlet = cv.properties_in[0], cv.properties_out[0]
    results = {
        "Flue gas temperature (K)": value(outlet.temperature),
        "Temperature rise (K)": value(outlet.temperature - inlet.temperature),
        "Flue gas mass flow (kg/s)": value(outlet.flow_mol * outlet.mw),
        "Flue gas molar flow (mol/s)": value(outlet.flow_mol),
    }
    for j in COMPONENTS:
        results["Flue gas " + j + " (mol/s)"] = value(outlet.flow_mol_comp[j])
    return results


def read_dwsim_results(path):
    """Read the combustor results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. The flue gas is the
    hottest material stream with a non-zero mass flow. Each compound is listed
    once per phase; the first numeric molar flow is the overall stream value.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    reactor = xml[xml.find("Reactors.Reactor_Conversion<"):]
    # DWSIM compound names (ChemSep and CoolProp databases name CO2 differently)
    names = {
        "methane": "Methane",
        "oxygen": "Oxygen",
        "nitrogen": "Nitrogen",
        "carbon_dioxide": "Carbon dioxide|CarbonDioxide",
        "water": "Water",
    }

    streams = []
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        mixture = stream[stream.find("<ComponentName>Mixture<"):]
        if first("massflow", mixture) > 0:
            streams.append((first("temperature", mixture), stream, mixture))
    _, flue, mixture = max(streams, key=lambda s: s[0])

    results = {
        "Flue gas temperature (K)": first("temperature", mixture),
        "Temperature rise (K)": first("DeltaT", reactor),
        "Flue gas mass flow (kg/s)": first("massflow", mixture),
        "Flue gas molar flow (mol/s)": first("molarflow", mixture),
    }
    for j, name in names.items():
        flows = re.findall(
            r"<ComponentName>(?:{0})</ComponentName>.*?<MolarFlow>([^<]+)<".format(name), flue, re.S
        )
        results["Flue gas " + j + " (mol/s)"] = next(float(f) for f in flows if f != "NaN")
    return results


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": reactor_results(model), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"].where(df["DWSIM"].abs() > 1e-6) * 100
    return df


if __name__ == "__main__":
    model, result = setup_conversion_reactor_model()
    print(result.solver.termination_condition)
    report_conversion_reactor_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Conversion Reactor Model.dwxmz")))
