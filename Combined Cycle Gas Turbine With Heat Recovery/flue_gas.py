#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# Gas phase property package for the combined cycle gas turbine: air, natural
# gas (methane, ethane, propane) and flue gas, using the IDAES Generic Property
# Package Framework.
#################################################################################
"""
Air, fuel gas and flue gas with the Peng-Robinson equation of state (vapor only).

This extends ``Conversion Reactor/combustion_gas.py`` with ethane, propane and
argon:

* Enthalpies include the ideal gas heats of formation at 25 C, so the heat of
  combustion is accounted for by the enthalpy balance of the combustor itself.
* The Poling ideal gas heat capacity correlations are only valid up to 1000 K,
  well below the turbine inlet temperature, so the ideal gas heat capacities are
  fitted to CoolProp's reference values over 298-2000 K.
* The binary interaction parameters (kij) are the Peng-Robinson values DWSIM uses
  for these compounds (read from the DWSIM flowsheet), so that the IDAES model can
  be validated against the DWSIM model.

The flue gas leaves the stack at about 105 C and 1 bar with ~8 mol % water, well
above its dew point (~40 C), so a vapor-only package is sufficient.
"""

__author__ = "Viraj Desai"

import chemicals
import numpy as np
import CoolProp.CoolProp as CP
from idaes.core import Component
from idaes.core import PhaseType as PT
from idaes.core import VaporPhase
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.pure import RPP5
from idaes.models.properties.modular_properties.state_definitions import FTPx
from pyomo.environ import units as pyunits

# Component name: (CAS number, CoolProp name, elemental composition)
COMPONENTS = {
    "methane": ("74-82-8", "Methane", {"C": 1, "H": 4}),
    "ethane": ("74-84-0", "Ethane", {"C": 2, "H": 6}),
    "propane": ("74-98-6", "n-Propane", {"C": 3, "H": 8}),
    "oxygen": ("7782-44-7", "Oxygen", {"O": 2}),
    "nitrogen": ("7727-37-9", "Nitrogen", {"N": 2}),
    "argon": ("7440-37-1", "Argon", {"Ar": 1}),
    "carbon_dioxide": ("124-38-9", "CarbonDioxide", {"C": 1, "O": 2}),
    "water": ("7732-18-5", "Water", {"H": 2, "O": 1}),
}

# Peng-Robinson binary interaction parameters (kij) from the DWSIM flowsheet
KIJ = {
    ("nitrogen", "argon"): -0.0026,
    ("nitrogen", "oxygen"): -0.0159,
    ("nitrogen", "methane"): 0.0289,
    ("nitrogen", "ethane"): 0.0344,
    ("nitrogen", "carbon_dioxide"): -0.0222,
    ("nitrogen", "propane"): 0.0878,
    ("argon", "oxygen"): 0.0107,
    ("argon", "methane"): 0.023,
    ("methane", "ethane"): -0.0033,
    ("methane", "carbon_dioxide"): 0.0793,
    ("methane", "propane"): 0.0119,
    ("carbon_dioxide", "ethane"): 0.1326,
    ("carbon_dioxide", "propane"): 0.1315,
    ("carbon_dioxide", "water"): -0.12155,
    ("ethane", "propane"): 0.0011,
    ("water", "methane"): 0.5,
    ("water", "ethane"): 0.5,
    ("water", "propane"): 0.5,
}


def cp_ig_coefficients(coolprop_name, T_min=275.0, T_max=2000.0):
    """Fit Cp/R = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4 to CoolProp's ideal gas Cp."""
    T = np.linspace(T_min, T_max, 200)
    cp = [CP.PropsSI("Cp0molar", "T", t, "P", 1e3, coolprop_name) for t in T]
    cp_over_R = np.array(cp) / CP.PropsSI("gas_constant", coolprop_name)
    return np.polynomial.polynomial.Polynomial.fit(T, cp_over_R, 4).convert().coef


def gas_component(name):
    """Vapor-only component with fitted ideal gas Cp and critical properties."""
    cas, coolprop_name, elements = COMPONENTS[name]
    a = cp_ig_coefficients(coolprop_name)
    return {
        "type": Component,
        "elemental_composition": elements,
        "enth_mol_ig_comp": RPP5,
        "entr_mol_ig_comp": RPP5,
        "valid_phase_types": PT.vaporPhase,
        "parameter_data": {
            "mw": (chemicals.identifiers.MW(ID=cas) / 1000, pyunits.kg / pyunits.mol),
            "pressure_crit": (chemicals.critical.Pc(CASRN=cas), pyunits.Pa),
            "temperature_crit": (chemicals.critical.Tc(CASRN=cas), pyunits.K),
            "omega": (chemicals.acentric.omega(CASRN=cas), None),
            "cp_mol_ig_comp_coeff": {
                "a0": (a[0], None),
                "a1": (a[1], pyunits.K**-1),
                "a2": (a[2], pyunits.K**-2),
                "a3": (a[3], pyunits.K**-3),
                "a4": (a[4], pyunits.K**-4),
            },
            "entr_mol_form_vap_comp_ref": (chemicals.reaction.S0g(CASRN=cas), pyunits.J / pyunits.mol / pyunits.K),
            "enth_mol_form_vap_comp_ref": (chemicals.reaction.Hfg(CASRN=cas), pyunits.J / pyunits.mol),
        },
    }


def kij(i, j):
    return KIJ.get((i, j), KIJ.get((j, i), 0.0))


configuration = {
    "components": {name: gas_component(name) for name in COMPONENTS},
    "phases": {
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        }
    },
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 1000, 10000, pyunits.mol / pyunits.s),
        "temperature": (250, 300, 2500, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },
    # Heats of formation are referenced to 25 C
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    "parameter_data": {
        "PR_kappa": {(i, j): kij(i, j) for i in COMPONENTS for j in COMPONENTS},
    },
}
