#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# This file demonstrates the creation of a methanol property package
# using the IDAES Generic Property Package Framework.
#
# The methanol property package includes thermodynamic properties
# and phase equilibrium calculations for pure methanol and ethanol.
#################################################################################
"""
Methanol-Ethanol phase equilibrium package using ideal liquid and vapor.

Example property package using the Generic Property Package Framework.
This example shows how to set up a property package to do methanol-ethanol
phase equilibrium in the generic framework using ideal liquid and vapor
assumptions along with methods drawn from the pre-built IDAES property
libraries.
"""
# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
import idaes.logger as idaeslog

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure import Perrys
from idaes.models.properties.modular_properties.pure import RPP4
from idaes.models.properties.modular_properties.pure import NIST


# Set up logger
_log = idaeslog.getLogger(__name__)


'''
# Data Sources:
[1]  The Properties of Gases and Liquids (1987) 4th edition, Chemical Engineering Series - Robert C. Reid
[2]  Perry's Chemical Engineers' Handbook 7th Ed.
[3]  https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=4#Thermo-Phase
[4]  https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=2#Thermo-Condensed
[5]  https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=1#Thermo-Gas
[6]  https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=4#Thermo-Phase
[7]  https://webbook.nist.gov/cgi/cbook.cgi?Name=ethanol&Units=SI
[8]  https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Units=SI&Mask=4#Thermo-Phase
[9]  https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Units=SI&Mask=1#Thermo-Gas
[10] https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Units=SI&Mask=2#Thermo-Condensed

'''


configuration = {
    # Specifying components
    "components": {
        "methanol": {
            "type": Component,
            "elemental_composition": {"C": 1, "H": 4, "O": 1},
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (32.0419e-3, pyunits.kg / pyunits.mol),  # [4]
                "pressure_crit": (81.1e5, pyunits.Pa),  # [4]
                "temperature_crit": (513.1, pyunits.K),  # [4]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (2.288, pyunits.kmol * pyunits.m**-3),  # [2]
                    "2": (0.2685, None),
                    "3": (512.64, pyunits.K),
                    "4": (0.2453, None),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (2.11e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (7.092e-2, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (2.587e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (-2.852e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (1.058e5, pyunits.J / pyunits.kmol / pyunits.K),  # [2]
                    "2": (-3.6223e2, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (9.3790e-01, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (-250.6e3, pyunits.J / pyunits.mol),  # [4]
                "enth_mol_form_vap_comp_ref": (-205.0e3, pyunits.J / pyunits.mol),  # [5]
                "pressure_sat_comp_coeff": {
                    "A": (5.15853, None),  # [6] Temp range 353.5 to 512.63 K
                    "B": (1569.613, pyunits.K),
                    "C": (-34.846, pyunits.K),
                },
            },
        },
        "ethanol": {
            "type": Component,
            "elemental_composition": {"C": 2, "H": 6, "O": 1},
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (46.068e-3, pyunits.kg / pyunits.mol),  # [7]
                "pressure_crit": (63.4e5, pyunits.Pa),  # [7]
                "temperature_crit": (514.7, pyunits.K),  # [7]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (1.648, pyunits.kmol * pyunits.m**-3),  # [2] pg. 2-98 To be updated
                    "2": (0.27627, None),
                    "3": (513.92, pyunits.K),
                    "4": (0.2331, None),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (9.014e0, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (2.141e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (-8.390e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (1.373e-9, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (1.0264e5, pyunits.J / pyunits.kmol / pyunits.K),  # [2] To be updated
                    "2": (-1.3963e2, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (-3.0341e-2, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (2.0386e-03, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (-276.2e3, pyunits.J / pyunits.mol),  # [10]
                "enth_mol_form_vap_comp_ref": (-234.2e3, pyunits.J / pyunits.mol),  # [9]
                "pressure_sat_comp_coeff": {
                    "A": (4.92531, None),  # [8] Temp range 364.8 to 513.91 K
                    "B": (1432.526, pyunits.K),
                    "C": (-61.819, pyunits.K),
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