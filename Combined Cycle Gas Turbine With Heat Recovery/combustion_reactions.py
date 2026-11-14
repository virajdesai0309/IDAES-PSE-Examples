#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# Reaction package for the complete combustion of natural gas in the gas
# turbine combustor, using the IDAES Generic Reaction Package Framework.
#################################################################################
"""
Complete combustion of the fuel gas (vapor phase), as in the DWSIM conversion reactor:

    CH4  + 2   O2 ->   CO2 + 2 H2O
    C2H6 + 3.5 O2 -> 2 CO2 + 3 H2O
    C3H8 + 5   O2 -> 3 CO2 + 4 H2O

A stoichiometric (conversion) reactor only needs the stoichiometry: the extents of
reaction are specified directly (100 % conversion of each hydrocarbon), so no rate
expressions are required. The heat of reaction is not given here because the
property package enthalpies already include the heats of formation.
"""

__author__ = "Viraj Desai"

from pyomo.environ import units as pyunits

# Reaction name: (fuel component, stoichiometry)
REACTIONS = {
    "methane_combustion": ("methane", {"methane": -1, "oxygen": -2, "carbon_dioxide": 1, "water": 2}),
    "ethane_combustion": ("ethane", {"ethane": -1, "oxygen": -3.5, "carbon_dioxide": 2, "water": 3}),
    "propane_combustion": ("propane", {"propane": -1, "oxygen": -5, "carbon_dioxide": 3, "water": 4}),
}

configuration = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "rate_reactions": {
        name: {"stoichiometry": {("Vap", c): nu for c, nu in stoich.items()}}
        for name, (_, stoich) in REACTIONS.items()
    },
}
