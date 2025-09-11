#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a gas stream for compressor modelling 
# The gas comprises of methane, ethane, propane, n-butane
# The gas stream is modelled using the IDAES Generic Property Package Framework.
# The property package includes thermodynamic properties
# and phase equilibrium calculations for the gas mixture.
#################################################################################
"""
Gas phase equilibrium package using ideal gas law.

Example property package using the Generic Property Package Framework
This example shows how to set up a property package to do gas phase equilibrium 
in the generic framework using ideal gas assumptions along with methods
libraries.
"""

__author__ = "Viraj Desai"

import chemicals
import idaes.logger as idaeslog
from chemicals.critical import Pc, Tc
from chemicals.heat_capacity import Cp_data_Poling
from chemicals.identifiers import MW
from chemicals.vapor_pressure import Psat_data_AntoinePoling
from chemicals.volume import rho_data_Perry_8E_105_l
# Import IDAES cores
from idaes.core import Component, LiquidPhase, VaporPhase
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import \
    IdealBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import \
    fugacity
from idaes.models.properties.modular_properties.pure import NIST, RPP5, Perrys
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
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP5,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (chemicals.identifiers.MW(ID="74-82-8") / 1000, pyunits.kg / pyunits.mol), 
                "pressure_crit": (chemicals.critical.Pc(CASRN="74-82-8"), pyunits.Pa),  
                "temperature_crit": (chemicals.critical.Tc(CASRN="74-82-8"), pyunits.K),  
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (rho_data_Perry_8E_105_l.loc["74-82-8"][1] / 1000, pyunits.kmol * pyunits.m**-3),  
                    "2": (rho_data_Perry_8E_105_l.loc["74-82-8"][2], None),
                    "3": (rho_data_Perry_8E_105_l.loc["74-82-8"][3], pyunits.K),
                    "4": (rho_data_Perry_8E_105_l.loc["74-82-8"][2], None),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (Cp_data_Poling.loc["74-82-8"][3], pyunits.J / pyunits.mol / pyunits.K),  
                    "B": (Cp_data_Poling.loc["74-82-8"][4], pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (Cp_data_Poling.loc["74-82-8"][5], pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (Cp_data_Poling.loc["74-82-8"][6], pyunits.J / pyunits.mol / pyunits.K**4),
                    "E": (Cp_data_Poling.loc["74-82-8"][7], pyunits.J / pyunits.mol / pyunits.K**5),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (1.058e5, pyunits.J / pyunits.kmol / pyunits.K),  
                    "2": (-3.6223e2, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (9.3790e-01, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (chemicals.reaction.Hfl(CASRN="74-82-8"), pyunits.J / pyunits.mol),  
                "enth_mol_form_vap_comp_ref": (chemicals.reaction.Hfg(CASRN="74-82-8"), pyunits.J / pyunits.mol),  
                "pressure_sat_comp_coeff": {
                    "A": (Psat_data_AntoinePoling.loc["74-82-8"][1] - 5, None),  
                    "B": (Psat_data_AntoinePoling.loc["74-82-8"][2], pyunits.K),
                    "C": (Psat_data_AntoinePoling.loc["74-82-8"][3] + 273.15.K),
                },
            },
        },
    },
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
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 450, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}