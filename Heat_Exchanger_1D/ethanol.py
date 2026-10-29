#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a liquid ethanol stream for heat
# exchanger modelling. The stream comprises of ethanol only.
# The stream is modelled using the IDAES Generic Property Package Framework.
#################################################################################
"""
Liquid-only ethanol property package.

The ethanol stays below its normal boiling point (78.3 C) in the heat
exchanger, so no vapor phase or phase equilibrium is needed. For a single
component, DWSIM's NRTL package reduces to an ideal liquid, which is what
this package models.
"""

__author__ = "Viraj Desai"


# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, Component
from idaes.core import PhaseType as PT

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import Perrys


'''
# Data Sources:
[1]  Perry's Chemical Engineers' Handbook 7th Ed.
[2]  https://webbook.nist.gov/cgi/cbook.cgi?Name=ethanol&Units=SI
[3]  https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Units=SI&Mask=2#Thermo-Condensed
'''


configuration = {
    # Specifying components
    "components": {
        "ethanol": {
            "type": Component,
            "elemental_composition": {"C": 2, "H": 6, "O": 1},
            "valid_phase_types": PT.liquidPhase,
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "parameter_data": {
                "mw": (46.068e-3, pyunits.kg / pyunits.mol),  # [2]
                "pressure_crit": (63.4e5, pyunits.Pa),  # [2]
                "temperature_crit": (514.7, pyunits.K),  # [2]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (1.648, pyunits.kmol * pyunits.m**-3),  # [1]
                    "2": (0.27627, None),
                    "3": (513.92, pyunits.K),
                    "4": (0.2331, None),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (1.0264e5, pyunits.J / pyunits.kmol / pyunits.K),  # [1]
                    "2": (-1.3963e2, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (-3.0341e-2, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (2.0386e-03, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (-276.2e3, pyunits.J / pyunits.mol),  # [3]
            },
        },
    },
    # Specifying phases
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
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
        "flow_mol": (0, 20, 100, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 400, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
}
