#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a reaction package for the water-gas
# shift equilibrium using the IDAES Generic Reaction Package Framework.
#################################################################################
"""
Water-gas shift equilibrium: CO + H2O <-> CO2 + H2 (vapor phase).

The equilibrium constant is calculated from the standard Gibbs energy of
reaction at the reactor temperature,

    ln K(T) = -sum(nu_j * (h_j(T) - T * s_j(T))) / (R * T),

using the ideal gas enthalpies and entropies of the property package. This
accounts for the temperature dependence of the heat of reaction, like DWSIM's
Gibbs energy based equilibrium constant. IDAES's built-in gibbs_energy and
van_t_hoff forms assume a constant heat of reaction instead.
"""

__author__ = "Viraj Desai"

# Import Pyomo units and functions
from pyomo.environ import exp, units as pyunits
from idaes.core.util.constants import Constants

from idaes.models.properties.modular_properties.base.utility import ConcentrationForm
from idaes.models.properties.modular_properties.pure import RPP5
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    power_law_equil,
)


class gibbs_energy_ig:
    """Equilibrium constant from the ideal gas Gibbs energy of reaction at T."""

    @staticmethod
    def build_parameters(rblock, config):
        # Keep the stoichiometry to sum the component Gibbs energies
        rblock.gibbs_stoichiometry = dict(config.stoichiometry)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        state = b.state_ref
        dg_rxn = sum(
            nu
            * (
                RPP5.enth_mol_ig_comp.return_expression(state, state.params.get_component(j), T)
                - T * RPP5.entr_mol_ig_comp.return_expression(state, state.params.get_component(j), T)
            )
            for (p, j), nu in rblock.gibbs_stoichiometry.items()
        )
        return exp(-dg_rxn / (Constants.gas_constant * T))

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1


configuration = {
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "water_gas_shift": {
            "stoichiometry": {
                ("Vap", "carbon_monoxide"): -1,
                ("Vap", "water"): -1,
                ("Vap", "carbon_dioxide"): 1,
                ("Vap", "hydrogen"): 1,
            },
            "equilibrium_constant": gibbs_energy_ig,
            "equilibrium_form": power_law_equil,
            # Equal moles on both sides, so K is the same on a mole fraction or partial pressure basis
            "concentration_form": ConcentrationForm.moleFraction,
        },
    },
}
