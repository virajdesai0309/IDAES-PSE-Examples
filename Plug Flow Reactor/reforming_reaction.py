#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a heterogeneous catalytic reaction
# package for steam methane reforming using the IDAES Generic Reaction Package
# Framework.
#################################################################################
"""
Steam methane reforming: CH4 + H2O <-> 3 H2 + CO (vapor phase, catalytic)

Rate per kg of catalyst, in the form of Xu & Froment (1989) as entered in the
DWSIM model, with partial pressures p_j in Pa:

    r = k / p_H2^2.5 * (p_CH4 * p_H2O - p_H2^3 * p_CO / K) / DEN^2   [mol/kg/s]
    k   = 4.22e15 * exp(-240100 / (R*T))
    K   = exp(30.42 - 27106 / T)
    DEN = 1 + K_CH4*p_CH4 + K_CO*p_CO + K_H2*p_H2 + K_H2O*p_H2O / p_H2

IDAES reactors need a rate per unit reactor volume, so the catalytic rate is
multiplied by the catalyst loading (kg of catalyst per m3 of reactor).
"""

__author__ = "Viraj Desai"

# Import Pyomo units and functions
from pyomo.environ import exp, units as pyunits


class xu_froment_rate:
    """Volumetric rate [mol/m3/s] = catalyst loading * Xu-Froment rate per kg of catalyst."""

    @staticmethod
    def build_parameters(rblock, config):
        rblock.catalyst_loading = config.parameter_data["catalyst_loading"]  # kg/m3

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        state = b.state_ref
        t = pyunits.convert(T, to_units=pyunits.K) / pyunits.K
        P = pyunits.convert(state.pressure, to_units=pyunits.Pa) / pyunits.Pa

        def p(j):  # partial pressure [Pa]
            return P * state.mole_frac_comp[j]

        R = 8.314
        k = 4.22e15 * exp(-240100 / (R * t))
        K = exp(30.42 - 27106 / t)
        den = (
            1
            + 6.65e-4 * exp(38280 / (R * t)) * p("methane")
            + 8.23e-5 * exp(70650 / (R * t)) * p("carbon_monoxide")
            + 6.12e-9 * exp(82900 / (R * t)) * p("hydrogen")
            + 1.77e5 * exp(-88680 / (R * t)) * p("water") / p("hydrogen")
        )
        rate = k / p("hydrogen") ** 2.5 * (
            p("methane") * p("water") - p("hydrogen") ** 3 * p("carbon_monoxide") / K
        ) / den**2  # mol/kg/s
        return rblock.catalyst_loading * rate * pyunits.mol / pyunits.m**3 / pyunits.s


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
        "steam_reforming": {
            "stoichiometry": {
                ("Vap", "methane"): -1,
                ("Vap", "water"): -1,
                ("Vap", "hydrogen"): 3,
                ("Vap", "carbon_monoxide"): 1,
            },
            "rate_form": xu_froment_rate,
            "parameter_data": {
                "catalyst_loading": 0.366,  # kg catalyst / m3 reactor
            },
        },
    },
}
