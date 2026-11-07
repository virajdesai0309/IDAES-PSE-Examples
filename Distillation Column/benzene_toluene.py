#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a benzene-toluene property package for
# distillation modelling.
# The mixture is modelled using the IDAES Generic Property Package Framework.
#################################################################################
"""
Benzene-toluene property package using the Peng-Robinson equation of state.

The pure component constants (critical properties, acentric factors, ideal
gas heats of formation and the DIPPR 16 ideal gas heat capacities) are the
ChemSep values used by the DWSIM model.

Three configurations are built from the same component data:

* configuration: vapor and liquid in equilibrium, for streams
* liquid_configuration: liquid phase only
* vapor_configuration: vapor phase only

The single phase packages evaluate a phase at a given temperature, pressure
and composition, for example the K-values K_j = phi_j(liquid) / phi_j(vapor)
of a liquid feed, or the enthalpy of a saturated vapor.
"""

__author__ = "Viraj Desai"

import copy

import chemicals
import numpy as np
# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.models.properties.modular_properties.eos.ceos import (Cubic,
                                                                 CubicType)
from idaes.models.properties.modular_properties.phase_equil.smooth_VLE_2 import CubicComplementarityVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.pure import RPP5
from idaes.models.properties.modular_properties.state_definitions import FTPx
# Import Pyomo units and functions
from pyomo.environ import exp, log, units as pyunits


'''
# Data Sources:
ChemSep database constants as stored in the DWSIM model file
(DWSim Shortcut Column Model.dwxmz, the same data as the distillation column model). Standard entropies from the python
chemicals database → https://chemicals.readthedocs.io/index.html
'''

# ChemSep data from the DWSIM model
# ig_cp: DIPPR 16, Cp = A + exp(B/T + C + D*T + E*T^2) [J/kmol/K]
# psat: DIPPR 101, ln(P) = A + B/T + C*ln(T) + D*T^E [Pa]
CHEMSEP = {
    "benzene": {
        "cas": "71-43-2", "mw": 78.11184, "Tc": 562.05, "Pc": 4895000, "omega": 0.209,
        "elements": {"C": 6, "H": 6},
        "Hf_ig": 1061.0427305258715 * 78.11184,  # J/mol
        "ig_cp": (34010.24, -588.0978, 12.81777, -0.000197306, 5.142899e-08),
        "psat": (88.368, -6712.9, -10.022, 7.694e-06, 2),
    },
    "toluene": {
        "cas": "108-88-3", "mw": 92.13843, "Tc": 591.75, "Pc": 4108000, "omega": 0.264,
        "elements": {"C": 7, "H": 8},
        "Hf_ig": 544.5067818064623 * 92.13843,
        "ig_cp": (47225, -565.85, 12.856, 5.535e-06, -1.998e-08),
        "psat": (32.89891, -5013.81, -1.348918, -1.869928e-06, 2),
    },
}


class dippr101:
    """Saturation pressure from DIPPR equation 101: ln(P) = A + B/T + C*ln(T) + D*T^E.

    Used only for the initial guesses of the bubble and dew point calculations;
    the phase equilibrium itself comes from the Peng-Robinson fugacities.
    """

    @staticmethod
    def build_parameters(cobj):
        A, B, C, D, E = CHEMSEP[cobj.local_name]["psat"]
        cobj.psat_coeff = {"A": A, "B": B, "C": C, "D": D, "E": E}

    @staticmethod
    def return_expression(b, cobj, T, dT=False):
        c = cobj.psat_coeff
        t = pyunits.convert(T, to_units=pyunits.K) / pyunits.K
        p = exp(c["A"] + c["B"] / t + c["C"] * log(t) + c["D"] * t ** c["E"]) * pyunits.Pa
        if dT:
            return p * (-c["B"] / t**2 + c["C"] / t + c["D"] * c["E"] * t ** (c["E"] - 1)) / pyunits.K
        units = b.params.get_metadata().derived_units
        return pyunits.convert(p, to_units=units.PRESSURE)

EOS = {"type": CubicType.PR}


def cp_ig_coefficients(c, T_min=250, T_max=600):
    """Fit Cp/R = a0 + ... + a4*T^4 to the ChemSep DIPPR 16 ideal gas heat capacity."""
    A, B, C, D, E = c
    T = np.linspace(T_min, T_max, 200)
    cp = (A + np.exp(B / T + C + D * T + E * T**2)) / 1000  # J/mol/K
    return np.polynomial.polynomial.Polynomial.fit(T, cp / 8.314462618, 4).convert().coef


def component(name):
    d = CHEMSEP[name]
    a = cp_ig_coefficients(d["ig_cp"])
    return {
        "type": Component,
        "elemental_composition": d["elements"],
        "enth_mol_ig_comp": RPP5,
        "entr_mol_ig_comp": RPP5,
        "pressure_sat_comp": dippr101,
        "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
        "parameter_data": {
            "mw": (d["mw"] / 1000, pyunits.kg / pyunits.mol),
            "pressure_crit": (d["Pc"], pyunits.Pa),
            "temperature_crit": (d["Tc"], pyunits.K),
            "omega": (d["omega"], None),
            "cp_mol_ig_comp_coeff": {
                "a0": (a[0], None),
                "a1": (a[1], pyunits.K**-1),
                "a2": (a[2], pyunits.K**-2),
                "a3": (a[3], pyunits.K**-3),
                "a4": (a[4], pyunits.K**-4),
            },
            "enth_mol_form_vap_comp_ref": (d["Hf_ig"], pyunits.J / pyunits.mol),
            "entr_mol_form_vap_comp_ref": (chemicals.reaction.S0g(CASRN=d["cas"]), pyunits.J / pyunits.mol / pyunits.K),
        },
    }


configuration = {
    # Specifying components
    "components": {name: component(name) for name in CHEMSEP},
    # Specifying phases
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Cubic, "equation_of_state_options": EOS},
        "Vap": {"type": VaporPhase, "equation_of_state": Cubic, "equation_of_state_options": EOS},
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
        "flow_mol": (0, 10, 1000, pyunits.mol / pyunits.s),
        "temperature": (250, 360, 500, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    # Heats of formation are referenced to 25 C
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): CubicComplementarityVLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {
        "PR_kappa": {(i, j): 0.0 for i in CHEMSEP for j in CHEMSEP},
    },
}


def single_phase_configuration(phase):
    """Copy of the configuration with only one phase ("Liq" or "Vap") and no phase equilibrium."""
    config = copy.deepcopy(configuration)
    config["phases"] = {phase: configuration["phases"][phase]}
    for key in ("phases_in_equilibrium", "phase_equilibrium_state", "bubble_dew_method"):
        config.pop(key)
    for comp in config["components"].values():
        comp.pop("phase_equilibrium_form")
    return config


liquid_configuration = single_phase_configuration("Liq")
vapor_configuration = single_phase_configuration("Vap")
