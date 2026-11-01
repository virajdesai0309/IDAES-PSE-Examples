#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a syngas mixture (carbon monoxide,
# water, carbon dioxide and hydrogen) for water-gas shift reactor modelling.
# The gas is modelled using the IDAES Generic Property Package Framework.
#################################################################################
"""
Gas phase property package for the water-gas shift reaction using the
Peng-Robinson equation of state.

Enthalpies and entropies include the ideal gas heats and entropies of
formation at 25 C, so the heat of reaction and the equilibrium constant both
follow from the component properties.

The ideal gas heat capacities are fitted to CoolProp's reference values over
298-2000 K, since the Poling correlations are only valid up to 1000 K.

Binary interaction parameters (kij) are the Peng-Robinson defaults used by
DWSIM, so that the IDAES model can be validated against the DWSIM model.
"""

__author__ = "Viraj Desai"

import chemicals
import numpy as np
import CoolProp.CoolProp as CP
# Import IDAES cores
from idaes.core import Component
from idaes.core import PhaseType as PT
from idaes.core import VaporPhase
from idaes.models.properties.modular_properties.eos.ceos import (Cubic,
                                                                 CubicType)
from idaes.models.properties.modular_properties.pure import RPP5
from idaes.models.properties.modular_properties.state_definitions import FTPx
# Import Pyomo units
from pyomo.environ import units as pyunits


'''
# Data Sources:
From python chemicals database → https://chemicals.readthedocs.io/index.html
Ideal gas heat capacities fitted to CoolProp → http://www.coolprop.org
'''

# Component name: (CAS number, CoolProp name, elemental composition)
COMPONENTS = {
    "carbon_monoxide": ("630-08-0", "CarbonMonoxide", {"C": 1, "O": 1}),
    "water": ("7732-18-5", "Water", {"H": 2, "O": 1}),
    "carbon_dioxide": ("124-38-9", "CarbonDioxide", {"C": 1, "O": 2}),
    "hydrogen": ("1333-74-0", "Hydrogen", {"H": 2}),
}

# Peng-Robinson binary interaction parameters (kij) from DWSIM
KIJ = {
    ("carbon_dioxide", "water"): -0.12155,
    ("hydrogen", "carbon_dioxide"): -0.1622,
    ("hydrogen", "carbon_monoxide"): 0.0919,
}


def cp_ig_coefficients(coolprop_name, T_min=298.15, T_max=2000):
    """Fit Cp/R = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4 to CoolProp ideal gas Cp."""
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
    # Specifying components
    "components": {name: gas_component(name) for name in COMPONENTS},
    # Specifying phases
    "phases": {
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        }
    },
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (250, 400, 1500, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },
    # Heats and entropies of formation are referenced to 25 C and 1 atm
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    "parameter_data": {
        "PR_kappa": {(i, j): kij(i, j) for i in COMPONENTS for j in COMPONENTS},
    },
}
