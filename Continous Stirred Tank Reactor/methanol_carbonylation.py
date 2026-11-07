#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a vapor-liquid property package for
# methanol carbonylation (methanol, water, acetic acid and carbon monoxide).
# The mixture is modelled using the IDAES Generic Property Package Framework.
#################################################################################
"""
Vapor-liquid property package for methanol carbonylation.

The DWSIM model uses the UNIQUAC activity coefficient model, which is not
available in the IDAES generic framework. This package approximates it with
an ideal liquid (Raoult's law) and takes the pure component data from the
same ChemSep database DWSIM uses:

* vapor pressure: ChemSep DIPPR equation 101 coefficients
* ideal gas heat capacity: polynomial fitted to the ChemSep DIPPR equation 16
* liquid density: ChemSep DIPPR equation 105 coefficients (Perry's form)
* liquid enthalpy: ideal gas enthalpy minus the heat of vaporization from the
  Pitzer corresponding states correlation, mirroring DWSIM's Lee-Kesler
  (corresponding states) liquid enthalpy

Carbon monoxide is far above its critical temperature and is treated as a
non-condensable gas (vapor phase only).
"""

__author__ = "Viraj Desai"

import numpy as np
# Import Pyomo units and functions
from pyomo.environ import exp, log, units as pyunits
from idaes.core.util.constants import Constants

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core import PhaseType as PT

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure import Perrys, RPP5


'''
# Data Sources:
ChemSep database coefficients as stored in the DWSIM model file
(DWSim CSTR Model.dwxmz). Water liquid density: Perry's Chemical Engineers'
Handbook 8th Ed., via the python chemicals database.
'''

# ChemSep data from the DWSIM model
# ig_cp: DIPPR 16, Cp = A + exp(B/T + C + D*T + E*T^2) [J/kmol/K]
# psat: DIPPR 101, ln(P) = A + B/T + C*ln(T) + D*T^E [Pa]
# rho_liq: DIPPR 105, rho = A / B^(1 + (1 - T/C)^D) [kmol/m3]
CHEMSEP = {
    "methanol": {
        "mw": 32.04186, "Tc": 512.64, "Pc": 8097000, "omega": 0.565,
        "elements": {"C": 1, "H": 4, "O": 1},
        "Hf_ig": -6271.1715237504941 * 32.04186,  # J/mol
        "Sf_ig": -4.0425903710759021 * 32.04186,  # J/mol/K
        "ig_cp": (36313.16, -680.4577, 11.10203, 0.000756766, -2.902645e-07),
        "psat": (73.40342, -6548.076, -7.409987, 5.72492e-06, 2),
        "rho_liq": (1.7918, 0.23929, 512.64, 0.21078),
    },
    "water": {
        "mw": 18.01528, "Tc": 647.14, "Pc": 22064000, "omega": 0.344,
        "elements": {"H": 2, "O": 1},
        "Hf_ig": -13422.716716032168 * 18.01528,
        "Sf_ig": -2.4619941145635971 * 18.01528,
        "ig_cp": (33200, -878.9001, 8.436956, 0.00207627, -6.467085e-07),
        "psat": (74.55502, -7295.586, -7.442448, 4.2881e-06, 2),
        "rho_liq": (5.459, 0.30542, 647.13, 0.081),  # Perry's
    },
    "acetic_acid": {
        "mw": 60.05196, "Tc": 594.45, "Pc": 5790000, "omega": 0.46,
        "elements": {"C": 2, "H": 4, "O": 2},
        "Hf_ig": -7207.091991668548 * 60.05196,
        "Sf_ig": -3.2561661291633914 * 60.05196,
        "ig_cp": (40110, -588.24, 12.017, 0.00016249, -8.6918e-08),
        "psat": (87.50607, -7603.906, -9.655308, 7.168835e-06, 2),
        "rho_liq": (1.0627, 0.22174, 594.76, 0.22566),
    },
    "carbon_monoxide": {
        "mw": 28.0101, "Tc": 132.85, "Pc": 3494000, "omega": 0.045,
        "elements": {"C": 1, "O": 1},
        "Hf_ig": -3946.0765938000932 * 28.0101,
        "Sf_ig": 3.1875615399872093 * 28.0101,
        "ig_cp": (29100, -1979.753, 10.58274, -7.90406e-05, -1.99685e-07),
    },
}


def cp_ig_coefficients(c, T_min=250, T_max=600):
    """Fit Cp/R = a0 + ... + a4*T^4 to the ChemSep DIPPR 16 ideal gas heat capacity."""
    A, B, C, D, E = c
    T = np.linspace(T_min, T_max, 200)
    cp = (A + np.exp(B / T + C + D * T + E * T**2)) / 1000  # J/mol/K
    return np.polynomial.polynomial.Polynomial.fit(T, cp / 8.314462618, 4).convert().coef


class dippr101:
    """Saturation pressure from DIPPR equation 101: ln(P) = A + B/T + C*ln(T) + D*T^E."""

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


class liquid_enthalpy_pitzer:
    """Liquid enthalpy = ideal gas enthalpy - heat of vaporization (Pitzer correlation).

    dHvap/(R*Tc) = 7.08*(1 - Tr)^0.354 + 10.95*omega*(1 - Tr)^0.456
    """

    @staticmethod
    def build_parameters(cobj):
        pass

    @staticmethod
    def return_expression(b, cobj, T):
        tr = 1 - T / cobj.temperature_crit
        dh_vap = Constants.gas_constant * cobj.temperature_crit * (
            7.08 * tr**0.354 + 10.95 * cobj.omega * tr**0.456
        )
        units = b.params.get_metadata().derived_units
        return pyunits.convert(
            RPP5.enth_mol_ig_comp.return_expression(b, cobj, T) - dh_vap,
            to_units=units.ENERGY_MOLE,
        )


def component(name):
    d = CHEMSEP[name]
    a = cp_ig_coefficients(d["ig_cp"])
    comp = {
        "type": Component,
        "elemental_composition": d["elements"],
        "enth_mol_ig_comp": RPP5,
        "entr_mol_ig_comp": RPP5,
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
            "entr_mol_form_vap_comp_ref": (d["Sf_ig"], pyunits.J / pyunits.mol / pyunits.K),
        },
    }
    if "psat" in d:
        A, B, C, D = d["rho_liq"]
        comp.update({
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": liquid_enthalpy_pitzer,
            "pressure_sat_comp": dippr101,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
        })
        comp["parameter_data"]["dens_mol_liq_comp_coeff"] = {
            "eqn_type": 1,
            "1": (A, pyunits.kmol * pyunits.m**-3),
            "2": (B, None),
            "3": (C, pyunits.K),
            "4": (D, None),
        }
    else:
        # Supercritical gas: present in the vapor phase only
        comp["valid_phase_types"] = PT.vaporPhase
    return comp


configuration = {
    # Specifying components
    "components": {name: component(name) for name in CHEMSEP},
    # Specifying phases
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
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
        "flow_mol": (0, 40, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 330, 450, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    # Heats of formation are referenced to 25 C
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}
