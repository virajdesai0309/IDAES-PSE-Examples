#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a reaction package for methane
# combustion using the IDAES Generic Reaction Package Framework.
#################################################################################
"""
Methane combustion reaction: CH4 + 2 O2 -> CO2 + 2 H2O (vapor phase).

A stoichiometric (conversion) reactor only needs the stoichiometry: the extent
of reaction is specified directly, so no rate expression is required. The heat
of reaction is not given here because the property package enthalpies already
include the heats of formation.
"""

__author__ = "Viraj Desai"

# Import Pyomo units
from pyomo.environ import units as pyunits


configuration = {
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "rate_reactions": {
        "combustion": {
            "stoichiometry": {
                ("Vap", "methane"): -1,
                ("Vap", "oxygen"): -2,
                ("Vap", "carbon_dioxide"): 1,
                ("Vap", "water"): 2,
            },
        },
    },
}
