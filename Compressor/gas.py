#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a gas stream for compressor modelling 
# The gas comprises of methane only.
# The gas stream is modelled using the IDAES Generic Property Package Framework.
# The property package includes thermodynamic properties
# and phase equilibrium calculations for the gas mixture.
#################################################################################
"""
Gas phase equilibrium package using Peng-Robinson equation of state.

Example property package using the Generic Property Package Framework
This example shows how to set up a property package to do gas phase equilibrium 
in the generic framework using ideal gas assumptions along with methods
libraries.
"""

__author__ = "Viraj Desai"

# Import Python libraries
import logging

import chemicals
import idaes.logger as idaeslog
from chemicals.critical import Pc, Tc
from chemicals.heat_capacity import Cp_data_Poling
from chemicals.identifiers import MW
from chemicals.vapor_pressure import Psat_data_AntoinePoling
from chemicals.volume import rho_data_Perry_8E_105_l
# Import IDAES cores
from idaes.core import Component, LiquidPhase
from idaes.core import PhaseType as PT
from idaes.core import VaporPhase
from idaes.models.properties.modular_properties.eos.ceos import (Cubic,
                                                                 CubicType)
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew, LogBubbleDew)
from idaes.models.properties.modular_properties.phase_equil.forms import (
    fugacity, log_fugacity)
from idaes.models.properties.modular_properties.pure import (NIST, RPP4, RPP5,
                                                             Perrys)
from idaes.models.properties.modular_properties.state_definitions import FTPx
# Import Pyomo units
from pyomo.environ import units as pyunits

# Set up logger
_log = idaeslog.getLogger(__name__)


'''
# Data Sources:
From python chemicals database → https://chemicals.readthedocs.io/index.html 
'''


configuration = {
    # Specifying components
    "components": {
        "methane": {
            "type": Component,
            "elemental_composition": {"C": 1, "H": 4},
            "enth_mol_ig_comp": RPP5,
            "entr_mol_ig_comp": RPP5,
            "valid_phase_types": PT.vaporPhase,
            "parameter_data": {
                "mw": (chemicals.identifiers.MW(ID="74-82-8") / 1000, pyunits.kg / pyunits.mol), 
                "pressure_crit": (chemicals.critical.Pc(CASRN="74-82-8"), pyunits.Pa),  
                "temperature_crit": (chemicals.critical.Tc(CASRN="74-82-8"), pyunits.K),  
                "omega": (chemicals.acentric.omega(CASRN="74-82-8"), None),
                "cp_mol_ig_comp_coeff": {
                    "a0": (Cp_data_Poling.loc["74-82-8"][3], None),  
                    "a1": (Cp_data_Poling.loc["74-82-8"][4], pyunits.K**-1),
                    "a2": (Cp_data_Poling.loc["74-82-8"][5], pyunits.K**-2),
                    "a3": (Cp_data_Poling.loc["74-82-8"][6], pyunits.K**-3),
                    "a4": (Cp_data_Poling.loc["74-82-8"][7], pyunits.K**-4),
                },
                "entr_mol_form_vap_comp_ref": (chemicals.reaction.S0g(CASRN="74-82-8"), pyunits.J / pyunits.mol / pyunits.K),
                "enth_mol_form_vap_comp_ref": (chemicals.reaction.Hfg(CASRN="74-82-8"), pyunits.J / pyunits.mol),     
            },
        },
    },
    # Specifying phases
    "phases": {
        "Vap": {
            "type": VaporPhase, 
            "equation_of_state": Cubic, 
            "equation_of_state_options":{"type":CubicType.PR}
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
        "temperature": (273.15, 300, 1500, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
    "parameter_data":{
        "PR_kappa": {
            ("methane", "methane"): 0.000,
    }
},
}