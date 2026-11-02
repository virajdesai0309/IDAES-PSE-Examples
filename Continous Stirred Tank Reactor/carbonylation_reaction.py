#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a kinetic reaction package for
# methanol carbonylation using the IDAES Generic Reaction Package Framework.
#################################################################################
"""
Methanol carbonylation: CH3OH + CO -> CH3COOH

Rate law (as in the DWSIM model), per unit reactor volume:

    r = k * C_MeOH^1 * C_CO^0.5        [mol/m3/s]
    k = 3.5e6 * exp(-83680 / (R * T))  [m^1.5/mol^0.5/s]

DWSIM evaluates the concentrations on the "Mixture" phase: component molar
flow divided by the total volumetric flow of the stream. The custom rate
form below does the same, so the rate does not depend on how the stream is
split between vapor and liquid.
"""

__author__ = "Viraj Desai"

# Import Pyomo units
from pyomo.environ import units as pyunits

from idaes.models.properties.modular_properties.base.utility import ConcentrationForm
from idaes.models.properties.modular_properties.reactions.rate_constant import arrhenius
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn


class mixture_power_law_rate:
    """Power law rate with mixture concentrations C_j = F_j / Q_total.

    Uses the reaction orders from parameter_data["reaction_order"].
    """

    @staticmethod
    def build_parameters(rblock, config):
        pass

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        state = b.state_ref
        rate = b.k_rxn[r_idx]
        for (p, j), order in rblock.reaction_order.items():
            if order.value != 0:
                rate = rate * (state.flow_mol_comp[j] / state.flow_vol) ** order
        return rate


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
        "carbonylation": {
            "stoichiometry": {
                ("Liq", "methanol"): -1,
                ("Vap", "carbon_monoxide"): -1,
                ("Liq", "acetic_acid"): 1,
            },
            "rate_constant": arrhenius,
            "rate_form": mixture_power_law_rate,
            # Molar concentrations (mol/m3), so the Arrhenius constant has units m^1.5/mol^0.5/s
            "concentration_form": ConcentrationForm.molarity,
            # Only used to report the heat of reaction; the energy balance uses
            # the heats of formation in the property package enthalpies
            "heat_of_reaction": constant_dh_rxn,
            "parameter_data": {
                "arrhenius_const": (3.5e6, pyunits.m**1.5 * pyunits.mol**-0.5 * pyunits.s**-1),
                "energy_activation": (83680, pyunits.J / pyunits.mol),
                "dh_rxn_ref": (-121330, pyunits.J / pyunits.mol),
                "reaction_order": {("Liq", "methanol"): 1, ("Vap", "carbon_monoxide"): 0.5},
            },
        },
    },
}
