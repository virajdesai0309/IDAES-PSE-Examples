#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a light hydrocarbon gas mixture
# (methane, ethane and propane) for mixer modelling.
# The gas is modelled using the IDAES Generic Property Package Framework.
#################################################################################
"""
Gas phase property package using the Peng-Robinson equation of state.

Binary interaction parameters (kij) are the Peng-Robinson defaults used by
DWSIM, so that the IDAES model can be validated against the DWSIM model.
"""

__author__ = "Viraj Desai"

import chemicals
from chemicals.heat_capacity import Cp_data_Poling
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
'''

CAS = {
    "methane": "74-82-8",
    "ethane": "74-84-0",
    "propane": "74-98-6",
}

ELEMENTS = {
    "methane": {"C": 1, "H": 4},
    "ethane": {"C": 2, "H": 6},
    "propane": {"C": 3, "H": 8},
}

# Peng-Robinson binary interaction parameters (kij) from DWSIM
KIJ = {
    ("methane", "ethane"): -0.0033,
    ("methane", "propane"): 0.0119,
    ("ethane", "propane"): 0.0011,
}


def gas_component(name):
    """Vapor-only component with ideal gas Cp (Poling) and critical properties."""
    cas = CAS[name]
    cp = Cp_data_Poling.loc[cas]
    return {
        "type": Component,
        "elemental_composition": ELEMENTS[name],
        "enth_mol_ig_comp": RPP5,
        "entr_mol_ig_comp": RPP5,
        "valid_phase_types": PT.vaporPhase,
        "parameter_data": {
            "mw": (chemicals.identifiers.MW(ID=cas) / 1000, pyunits.kg / pyunits.mol),
            "pressure_crit": (chemicals.critical.Pc(CASRN=cas), pyunits.Pa),
            "temperature_crit": (chemicals.critical.Tc(CASRN=cas), pyunits.K),
            "omega": (chemicals.acentric.omega(CASRN=cas), None),
            "cp_mol_ig_comp_coeff": {
                "a0": (cp.iloc[3], None),
                "a1": (cp.iloc[4], pyunits.K**-1),
                "a2": (cp.iloc[5], pyunits.K**-2),
                "a3": (cp.iloc[6], pyunits.K**-3),
                "a4": (cp.iloc[7], pyunits.K**-4),
            },
            "entr_mol_form_vap_comp_ref": (chemicals.reaction.S0g(CASRN=cas), pyunits.J / pyunits.mol / pyunits.K),
            "enth_mol_form_vap_comp_ref": (chemicals.reaction.Hfg(CASRN=cas), pyunits.J / pyunits.mol),
        },
    }


def kij(i, j):
    return KIJ.get((i, j), KIJ.get((j, i), 0.0))


configuration = {
    # Specifying components
    "components": {name: gas_component(name) for name in CAS},
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
        "temperature": (250, 300, 500, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
    "parameter_data": {
        "PR_kappa": {(i, j): kij(i, j) for i in CAS for j in CAS},
    },
}
