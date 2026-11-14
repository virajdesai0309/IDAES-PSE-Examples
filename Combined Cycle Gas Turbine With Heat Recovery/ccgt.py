#################################################################################
# IDAES PSE Examples Repository
# https://github.com/virajdesai0309/IDAES-PSE-Examples
#
# Combined cycle gas turbine (CCGT) with a single pressure heat recovery steam
# generator (HRSG), validated against the DWSIM model of the same plant.
#################################################################################
"""
Combined cycle gas turbine with a single pressure HRSG (~27 MW).

Gas turbine (Brayton cycle, Peng-Robinson flue gas package ``flue_gas.py``):
    Air -> compressor -> combustor (stoichiometric reactor, with fuel) -> gas turbine
Heat recovery steam generator (counter-current, gas path superheater -> evaporator -> economizer):
    flue gas: GT exhaust -> superheater -> evaporator -> economizer -> stack
    water (IAPWS-95): feedwater -> feed pump -> economizer -> evaporator -> drum -> superheater
Steam cycle:
    superheater -> steam turbine -> condenser (saturated liquid)

As in the DWSIM model the water side is open: feedwater enters at the condensate
state, and the condensate and drum blowdown leave the flowsheet. The flowsheet
has no recycle in IDAES because the whole model is solved simultaneously.

Operating modes (see ``set_mode``):

* ``"design"``: the DWSIM specifications. Fuel and feedwater flows are fixed; the
  economizer approach, evaporator outlet quality and steam temperature are specs,
  and the heat exchanger areas are results.
* ``"targets"``: as design, but the fuel flow is found for a turbine inlet
  temperature and the feedwater flow for an HRSG pinch.
* ``"offdesign"``: the HRSG areas are fixed at their design values and the steam
  flow, economizer outlet and steam temperatures float. The compressor draws a fixed
  volumetric flow, so the air mass flow follows the ambient density, and the gas
  turbine's swallowing capacity (choked first stage nozzle, m*sqrt(T)/P = constant)
  sets the pressure ratio. The load is set by the turbine inlet temperature (or the
  fuel flow). Component efficiencies and the condenser pressure are held constant.
"""

__author__ = "Viraj Desai"

import os
import re
import zipfile

import pandas as pd
from CoolProp.CoolProp import PropsSI
from pyomo.environ import (ConcreteModel, Constraint, Var, Expression, TransformationFactory, value, sqrt,
                           units as pyunits)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.models.properties import iapws95
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import GenericReactionParameterBlock
from idaes.models.unit_models import (Feed, Product, Mixer, MomentumMixingType, StoichiometricReactor,
                                      HeatExchanger, Heater)
from idaes.models.unit_models.pressure_changer import Compressor, Turbine, Pump
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_smooth_callback
from idaes.models_extra.power_generation.unit_models.helm.phase_separator import HelmPhaseSeparator
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

from flue_gas import configuration as gas_configuration, COMPONENTS
from combustion_reactions import configuration as reaction_configuration, REACTIONS

HERE = os.path.dirname(os.path.abspath(__file__))
DWSIM_FILE = os.path.join(HERE, "DWSim CCGT Combined Cycle.dwxmz")
MW_WATER = 0.01801528  # kg/mol

# Design basis: identical to build_dwsim_ccgt.py (temperatures in C, pressures in Pa, flows in kg/s)
DESIGN = {
    "ambient_T": 15.0,
    "ambient_P": 101325.0,
    "air_flow": 40.0,
    "air": {"oxygen": 0.21, "nitrogen": 0.78, "argon": 0.01},
    "pressure_ratio": 17.0,
    "compressor_eff": 0.86,
    "fuel": {"methane": 0.90, "ethane": 0.06, "propane": 0.04},
    "fuel_T": 25.0,
    "fuel_flow": 0.9777,
    "combustor_dP": 0.5e5,
    "turbine_eff": 0.89,
    "turbine_outlet_P": 1.06e5,
    "steam_flow": 8.074,  # feedwater
    "condenser_P": 0.08e5,
    "pump_outlet_P": 62e5,
    "pump_eff": 0.80,
    "evaporator_P_drop": 1e5,  # water side, economizer + evaporator
    "superheater_P_drop": 1e5,
    "steam_T": 500.0,
    "approach": 10.0,  # economizer outlet below saturation, K
    "evaporator_quality": 0.98,
    "steam_turbine_eff": 0.85,
    "gas_dP_per_section": 1000.0,
    # Overall heat transfer coefficients of the finned tube HRSG sections (gas side controlled), W/m2/K.
    # They only scale the areas; the off-design behaviour depends on UA.
    "U": {"superheater": 60.0, "evaporator": 55.0, "economizer": 50.0},
    # Design targets used in "targets" mode
    "TIT": 1350.0,
    "pinch": 10.0,
}

# Lower heating values at 25 C, kJ/mol (water as vapor in the products)
LHV_MOLAR = {"methane": 802.3, "ethane": 1428.6, "propane": 2043.1}
HX_NAMES = ("superheater", "evaporator", "economizer")


def saturation_T(P):
    """Saturation temperature of water (K) at pressure P (Pa), as used in the DWSIM build script."""
    return PropsSI("T", "P", P, "Q", 0, "Water")


def mixture_mw(m, fractions):
    return sum(x * value(m.fs.gas.get_component(j).mw) for j, x in fractions.items())


def fuel_lhv(m, fractions):
    """LHV of the fuel gas in kJ/kg."""
    return sum(x * LHV_MOLAR[j] for j, x in fractions.items()) / mixture_mw(m, fractions)


def build_ccgt(design=DESIGN):
    """Build the CCGT flowsheet (not yet specified, initialized or solved)."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.gas = GenericParameterBlock(**gas_configuration)
    m.fs.reactions = GenericReactionParameterBlock(property_package=m.fs.gas, **reaction_configuration)
    m.fs.water = iapws95.Iapws95ParameterBlock()

    # ---------------- Gas turbine ----------------
    m.fs.air = Feed(property_package=m.fs.gas)
    m.fs.fuel = Feed(property_package=m.fs.gas)
    m.fs.compressor = Compressor(property_package=m.fs.gas)
    m.fs.mixer = Mixer(property_package=m.fs.gas, inlet_list=["air_inlet", "fuel_inlet"],
                       momentum_mixing_type=MomentumMixingType.minimize)
    m.fs.combustor = StoichiometricReactor(
        property_package=m.fs.gas,
        reaction_package=m.fs.reactions,
        has_heat_of_reaction=False,  # heats of formation are in the enthalpies
        has_heat_transfer=False,  # adiabatic
        has_pressure_change=True,
    )
    m.fs.gas_turbine = Turbine(property_package=m.fs.gas)

    # ---------------- HRSG ----------------
    def hrsg_section():
        return HeatExchanger(
            delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
            hot_side={"property_package": m.fs.gas, "has_pressure_change": True},
            cold_side={"property_package": m.fs.water, "has_pressure_change": True},
        )

    m.fs.superheater = hrsg_section()
    m.fs.evaporator = hrsg_section()
    m.fs.economizer = hrsg_section()
    m.fs.stack = Product(property_package=m.fs.gas)

    # ---------------- Steam cycle ----------------
    m.fs.feedwater = Feed(property_package=m.fs.water)
    m.fs.feed_pump = Pump(property_package=m.fs.water)
    m.fs.drum = HelmPhaseSeparator(property_package=m.fs.water)
    m.fs.blowdown = Product(property_package=m.fs.water)
    m.fs.steam_turbine = Turbine(property_package=m.fs.water)
    m.fs.condenser = Heater(property_package=m.fs.water)
    m.fs.condensate = Product(property_package=m.fs.water)

    # Streams carry the DWSIM stream names in the comments
    arcs = {
        "s_air": (m.fs.air.outlet, m.fs.compressor.inlet),  # Air
        "s_compressed_air": (m.fs.compressor.outlet, m.fs.mixer.air_inlet),  # Compressed Air
        "s_fuel": (m.fs.fuel.outlet, m.fs.mixer.fuel_inlet),  # Fuel
        "s_combustor_feed": (m.fs.mixer.outlet, m.fs.combustor.inlet),  # Combustor Feed
        "s_hot_gas": (m.fs.combustor.outlet, m.fs.gas_turbine.inlet),  # Hot Gas
        "s_gt_exhaust": (m.fs.gas_turbine.outlet, m.fs.superheater.hot_side_inlet),  # GT Exhaust
        "s_gas_after_sh": (m.fs.superheater.hot_side_outlet, m.fs.evaporator.hot_side_inlet),  # Gas After Superheater
        "s_gas_after_ev": (m.fs.evaporator.hot_side_outlet, m.fs.economizer.hot_side_inlet),  # Gas After Evaporator
        "s_stack": (m.fs.economizer.hot_side_outlet, m.fs.stack.inlet),  # Stack Gas
        "s_feedwater": (m.fs.feedwater.outlet, m.fs.feed_pump.inlet),  # Feedwater
        "s_pumped_water": (m.fs.feed_pump.outlet, m.fs.economizer.cold_side_inlet),  # Pumped Water
        "s_eco_water_out": (m.fs.economizer.cold_side_outlet, m.fs.evaporator.cold_side_inlet),  # Economizer Water Out
        "s_evap_outlet": (m.fs.evaporator.cold_side_outlet, m.fs.drum.inlet),  # Evaporator Outlet
        "s_drum_steam": (m.fs.drum.vap_outlet, m.fs.superheater.cold_side_inlet),  # Drum Steam
        "s_blowdown": (m.fs.drum.liq_outlet, m.fs.blowdown.inlet),  # Blowdown
        "s_superheated_steam": (m.fs.superheater.cold_side_outlet, m.fs.steam_turbine.inlet),  # Superheated Steam
        "s_turbine_exhaust": (m.fs.steam_turbine.outlet, m.fs.condenser.inlet),  # Turbine Exhaust Steam
        "s_condensate": (m.fs.condenser.outlet, m.fs.condensate.inlet),  # Condensate
    }
    for name, (src, dst) in arcs.items():
        setattr(m.fs, name, Arc(source=src, destination=dst))
    TransformationFactory("network.expand_arcs").apply_to(m)

    # ---------------- Specification variables and constraints ----------------
    t = m.fs.time.first()
    eco_out = m.fs.economizer.cold_side.properties_out[t]
    evap_out = m.fs.evaporator.cold_side.properties_out[t]
    sh_out = m.fs.superheater.cold_side.properties_out[t]
    drum = m.fs.drum.mixed_state[t]

    m.fs.approach = Var(initialize=design["approach"], units=pyunits.K, doc="Economizer approach to saturation")
    m.fs.evaporator_quality = Var(initialize=design["evaporator_quality"], doc="Evaporator outlet vapor fraction")
    m.fs.steam_T = Var(initialize=design["steam_T"] + 273.15, units=pyunits.K, doc="Superheated steam temperature")
    m.fs.TIT = Var(initialize=design["TIT"] + 273.15, units=pyunits.K, doc="Turbine inlet temperature")
    m.fs.pinch = Var(initialize=design["pinch"], units=pyunits.K, doc="HRSG pinch")

    m.fs.approach_spec = Constraint(expr=eco_out.temperature == eco_out.temperature_sat - m.fs.approach)
    m.fs.quality_spec = Constraint(expr=evap_out.vapor_frac == m.fs.evaporator_quality)
    m.fs.steam_T_spec = Constraint(expr=sh_out.temperature == m.fs.steam_T)
    m.fs.TIT_spec = Constraint(expr=m.fs.combustor.control_volume.properties_out[t].temperature == m.fs.TIT)
    m.fs.pinch_spec = Constraint(
        expr=m.fs.evaporator.hot_side.properties_out[t].temperature == drum.temperature_sat + m.fs.pinch)
    m.fs.TIT_spec.deactivate()
    m.fs.pinch_spec.deactivate()

    # Complete combustion: the extent of each reaction equals the inlet flow of its fuel
    # (less a negligible remainder, which keeps the Peng-Robinson log terms finite)
    feed = m.fs.combustor.control_volume.properties_in[t]
    m.fs.complete_combustion = Constraint(
        list(REACTIONS),
        rule=lambda fs, r: fs.combustor.rate_reaction_extent[t, r]
        == (1 - 1e-7) * feed.flow_mol * feed.mole_frac_comp[REACTIONS[r][0]])

    # Compressor inlet volumetric flow (fixed in off-design: the air mass flow follows density)
    air = m.fs.air.properties[t]
    m.fs.air_volumetric_flow = Var(initialize=33.0, units=pyunits.m**3 / pyunits.s)
    m.fs.air_volume_spec = Constraint(expr=air.flow_mol == m.fs.air_volumetric_flow * air.dens_mol)
    m.fs.air_volume_spec.deactivate()

    # Gas turbine swallowing capacity (off-design): flow * sqrt(T) / P is constant at the turbine inlet,
    # and the fuel is supplied at the combustor pressure
    gt_in = m.fs.gas_turbine.control_volume.properties_in[t]
    m.fs.turbine_flow_coefficient = Var(initialize=0.034, units=pyunits.mol / pyunits.s * pyunits.K**0.5 / pyunits.Pa)
    m.fs.turbine_swallowing = Constraint(
        expr=gt_in.flow_mol * sqrt(gt_in.temperature) == m.fs.turbine_flow_coefficient * gt_in.pressure)
    m.fs.fuel_pressure_link = Constraint(
        expr=m.fs.fuel.properties[t].pressure == m.fs.compressor.control_volume.properties_out[t].pressure)
    m.fs.turbine_swallowing.deactivate()
    m.fs.fuel_pressure_link.deactivate()

    # Key results
    m.fs.net_power = Expression(expr=-(m.fs.compressor.work_mechanical[t] + m.fs.gas_turbine.work_mechanical[t]
                                       + m.fs.steam_turbine.work_mechanical[t] + m.fs.feed_pump.work_mechanical[t]))
    m.fs.fuel_heat_input = Expression(
        expr=sum(m.fs.fuel.properties[t].flow_mol * m.fs.fuel.mole_frac_comp[t, j] * LHV_MOLAR[j] * 1000
                 for j in LHV_MOLAR))  # W
    m.fs.efficiency = Expression(expr=m.fs.net_power / m.fs.fuel_heat_input)
    return m


def scale_ccgt(m):
    """Scaling factors for the gas states, heat duties, areas and powers (W, m2, mol/s, K, Pa)."""
    for name, sf in (("flow_mol", 1e-3), ("flow_mol_phase", 1e-3), ("temperature", 1e-2), ("pressure", 1e-5)):
        m.fs.gas.set_default_scaling(name, sf)
    for n in HX_NAMES:
        hx = getattr(m.fs, n)
        iscale.set_scaling_factor(hx.area, 1e-3)
        iscale.set_scaling_factor(hx.heat_duty, 1e-6)
        iscale.set_scaling_factor(hx.overall_heat_transfer_coefficient, 1e-2)
    for u in (m.fs.compressor, m.fs.gas_turbine, m.fs.steam_turbine, m.fs.feed_pump):
        iscale.set_scaling_factor(u.control_volume.work, 1e-6)
    iscale.calculate_scaling_factors(m)


def get_ccgt_solver(max_iter=300):
    """IPOPT using the model's scaling factors."""
    return get_solver(options={"max_iter": max_iter, "nlp_scaling_method": "user-scaling"})


def fix_design(m, design=DESIGN):
    """Fix the DWSIM design inputs ("design" mode)."""
    t = m.fs.time.first()
    small = 1e-10  # zero mole fractions are replaced by a small value (Peng-Robinson log terms)

    def gas_feed(feed, flow_mass, T, P, fractions):
        feed.flow_mol.fix(flow_mass / mixture_mw(m, fractions))
        feed.temperature.fix(T + 273.15)
        feed.pressure.fix(P)
        for j in COMPONENTS:
            feed.mole_frac_comp[t, j].fix(max(fractions.get(j, 0.0), small))

    p_comb = design["ambient_P"] * design["pressure_ratio"]
    gas_feed(m.fs.air, design["air_flow"], design["ambient_T"], design["ambient_P"], design["air"])
    gas_feed(m.fs.fuel, design["fuel_flow"], design["fuel_T"], p_comb, design["fuel"])
    m.fs.compressor.outlet.pressure.fix(p_comb)
    m.fs.compressor.efficiency_isentropic.fix(design["compressor_eff"])
    m.fs.combustor.control_volume.deltaP.fix(-design["combustor_dP"])
    m.fs.gas_turbine.outlet.pressure.fix(design["turbine_outlet_P"])
    m.fs.gas_turbine.efficiency_isentropic.fix(design["turbine_eff"])

    # Feedwater: condensate just subcooled (as DWSIM)
    p_cond = design["condenser_P"]
    m.fs.feedwater.flow_mol.fix(design["steam_flow"] / MW_WATER)
    m.fs.feedwater.pressure.fix(p_cond)
    m.fs.feedwater.enth_mol.fix(value(m.fs.water.htpx(T=(saturation_T(p_cond) - 0.5) * pyunits.K, p=p_cond * pyunits.Pa)))
    m.fs.feed_pump.outlet.pressure.fix(design["pump_outlet_P"])
    m.fs.feed_pump.efficiency_pump.fix(design["pump_eff"])

    dp_water = {"economizer": design["evaporator_P_drop"] / 2, "evaporator": design["evaporator_P_drop"] / 2,
                "superheater": design["superheater_P_drop"]}
    for name in HX_NAMES:
        hx = getattr(m.fs, name)
        hx.hot_side.deltaP.fix(-design["gas_dP_per_section"])
        hx.cold_side.deltaP.fix(-dp_water[name])
        hx.overall_heat_transfer_coefficient.fix(design["U"][name])
    m.fs.approach.fix(design["approach"])
    m.fs.evaporator_quality.fix(design["evaporator_quality"])
    m.fs.steam_T.fix(design["steam_T"] + 273.15)
    m.fs.TIT.fix(design["TIT"] + 273.15)
    m.fs.pinch.fix(design["pinch"])

    m.fs.steam_turbine.outlet.pressure.fix(p_cond)
    m.fs.steam_turbine.efficiency_isentropic.fix(design["steam_turbine_eff"])
    m.fs.condenser.outlet.enth_mol.fix(value(m.fs.water.htpx(p=p_cond * pyunits.Pa, x=0)))  # saturated liquid


def set_mode(m, mode):
    """Switch the specifications: "design", "targets" or "offdesign" (see the module docstring)."""
    t = m.fs.time.first()
    fuel_flow = m.fs.fuel.flow_mol[t]
    feedwater_flow = m.fs.feedwater.flow_mol[t]
    air_flow = m.fs.air.flow_mol[t]
    areas = [getattr(m.fs, n).area for n in HX_NAMES]

    p_comb = m.fs.compressor.outlet.pressure[t]
    fuel_pressure = m.fs.fuel.pressure[t]
    if mode in ("design", "targets"):
        m.fs.air_volume_spec.deactivate()
        air_flow.fix()
        m.fs.turbine_swallowing.deactivate()
        m.fs.fuel_pressure_link.deactivate()
        p_comb.fix()
        fuel_pressure.fix()
        for a in areas:
            a.unfix()
        m.fs.approach_spec.activate()
        m.fs.steam_T_spec.activate()
        if mode == "design":
            m.fs.TIT_spec.deactivate()
            m.fs.pinch_spec.deactivate()
            fuel_flow.fix()
            feedwater_flow.fix()
        else:
            m.fs.TIT_spec.activate()
            m.fs.pinch_spec.activate()
            fuel_flow.unfix()
            feedwater_flow.unfix()
    elif mode == "offdesign":
        # Areas at their current (design) values; the steam side floats
        for a in areas:
            a.fix()
        m.fs.approach_spec.deactivate()
        m.fs.steam_T_spec.deactivate()
        m.fs.pinch_spec.deactivate()
        feedwater_flow.unfix()
        # Fixed compressor inlet volumetric flow at the current value
        air = m.fs.air.properties[t]
        m.fs.air_volumetric_flow.fix(value(air.flow_mol / air.dens_mol))
        air_flow.unfix()
        m.fs.air_volume_spec.activate()
        # Turbine swallowing capacity at its current value sets the pressure ratio
        gt_in = m.fs.gas_turbine.control_volume.properties_in[t]
        m.fs.turbine_flow_coefficient.fix(value(gt_in.flow_mol * sqrt(gt_in.temperature) / gt_in.pressure))
        m.fs.turbine_swallowing.activate()
        p_comb.unfix()
        m.fs.fuel_pressure_link.activate()
        fuel_pressure.unfix()
        # Load set by the turbine inlet temperature
        m.fs.TIT.fix(value(m.fs.combustor.control_volume.properties_out[t].temperature))
        m.fs.TIT_spec.activate()
        fuel_flow.unfix()
    else:
        raise ValueError("mode must be 'design', 'targets' or 'offdesign'")
    m.fs.mode = mode


def init_unit(unit, outlvl=idaeslog.WARNING, T_out=None):
    """Initialize a unit by solving it alone with its inlet fixed.

    With an outlet temperature guess (used for the compressor and turbines) the unit is solved directly:
    the IDAES isentropic initialization starts from the inlet temperature, which is far off for a pressure
    ratio of 17. Otherwise the IDAES routine is tried first.
    """
    t = unit.flowsheet().time.first()
    if T_out is not None:
        for name in ("properties_isentropic",):
            if hasattr(unit, name):
                getattr(unit, name)[t].temperature.value = T_out
        unit.control_volume.properties_out[t].temperature.value = T_out
    if T_out is None:
        try:
            unit.initialize(outlvl=outlvl)
            return
        except Exception:  # noqa: BLE001  fall back to a direct solve
            pass
    fixed = [v for v in unit.inlet.vars.values() for v in v.values() if not v.fixed]
    for v in fixed:
        v.fix()
    result = get_solver(options={"max_iter": 300}).solve(unit)
    for v in fixed:
        v.unfix()
    if str(result.solver.termination_condition) != "optimal":
        raise RuntimeError("{0} failed to initialize".format(unit.name))


def initialize_ccgt(m, design=DESIGN, outlvl=idaeslog.WARNING):
    """Initialize unit by unit along the gas path, then the water path, and solve the design point.

    The HRSG sections are first initialized at the DWSIM duties (a guess), with the
    steam side inlets estimated where the flowsheet is counter-current.
    """
    t = m.fs.time.first()
    solver = get_ccgt_solver()

    m.fs.air.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s_air)
    init_unit(m.fs.compressor, outlvl, T_out=690.0)
    propagate_state(arc=m.fs.s_compressed_air)
    m.fs.fuel.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s_fuel)
    m.fs.mixer.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s_combustor_feed)
    feed = m.fs.combustor.control_volume.properties_in[t]
    for r, (fuel, _) in REACTIONS.items():
        m.fs.combustor.rate_reaction_extent[t, r].fix((1 - 1e-7) * value(feed.flow_mol * feed.mole_frac_comp[fuel]))
    m.fs.complete_combustion.deactivate()
    m.fs.combustor.initialize(outlvl=outlvl)
    m.fs.combustor.rate_reaction_extent.unfix()
    m.fs.complete_combustion.activate()
    propagate_state(arc=m.fs.s_hot_gas)
    init_unit(m.fs.gas_turbine, outlvl, T_out=940.0)
    propagate_state(arc=m.fs.s_gt_exhaust)

    m.fs.feedwater.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s_feedwater)
    init_unit(m.fs.feed_pump, outlvl)
    propagate_state(arc=m.fs.s_pumped_water)

    # Steam side inlet estimates for the counter-current sections
    flow = value(m.fs.feedwater.flow_mol[t])
    p_eco = design["pump_outlet_P"] - design["evaporator_P_drop"] / 2
    p_drum = design["pump_outlet_P"] - design["evaporator_P_drop"]
    sh_in = m.fs.superheater.cold_side_inlet
    sh_in.flow_mol[t].value = flow * design["evaporator_quality"]
    sh_in.pressure[t].value = p_drum
    sh_in.enth_mol[t].value = value(m.fs.water.htpx(p=p_drum * pyunits.Pa, x=1))
    ev_in = m.fs.evaporator.cold_side_inlet
    ev_in.flow_mol[t].value = flow
    ev_in.pressure[t].value = p_eco
    ev_in.enth_mol[t].value = value(m.fs.water.htpx(T=(saturation_T(p_eco) - design["approach"]) * pyunits.K,
                                                    p=p_eco * pyunits.Pa))

    duty_guess = {"superheater": 5.06e6, "evaporator": 12.77e6, "economizer": 8.00e6}
    for name, arc in (("superheater", m.fs.s_gas_after_sh), ("evaporator", m.fs.s_gas_after_ev),
                      ("economizer", m.fs.s_stack)):
        hx = getattr(m.fs, name)
        hx.area.fix(1000.0)
        hx.heat_duty.fix(duty_guess[name] * flow * MW_WATER / DESIGN["steam_flow"])
        hx.area.unfix()
        hx.initialize(outlvl=outlvl)
        hx.heat_duty.unfix()
        propagate_state(arc=arc)
    m.fs.stack.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s_eco_water_out)
    propagate_state(arc=m.fs.s_evap_outlet)
    m.fs.drum.initialize(outlvl=outlvl)
    for arc in (m.fs.s_drum_steam, m.fs.s_blowdown):
        propagate_state(arc=arc)
    m.fs.blowdown.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s_superheated_steam)
    init_unit(m.fs.steam_turbine, outlvl, T_out=315.0)
    propagate_state(arc=m.fs.s_turbine_exhaust)
    m.fs.condenser.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s_condensate)
    m.fs.condensate.initialize(outlvl=outlvl)
    return solver.solve(m)


def setup_ccgt_model(design=DESIGN, outlvl=idaeslog.WARNING):
    """Build, initialize and solve the CCGT at its design point (DWSIM specifications)."""
    m = build_ccgt(design)
    fix_design(m, design)
    set_mode(m, "design")
    scale_ccgt(m)
    print("The DOF is {0}".format(degrees_of_freedom(m)))
    result = initialize_ccgt(m, design, outlvl)
    return m, result


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
GAS_STREAMS = {  # DWSIM tag: (unit, port)
    "Air": ("air", "outlet"),
    "Compressed Air": ("compressor", "outlet"),
    "Fuel": ("fuel", "outlet"),
    "Hot Gas": ("combustor", "outlet"),
    "GT Exhaust": ("gas_turbine", "outlet"),
    "Gas After Superheater": ("superheater", "hot_side_outlet"),
    "Gas After Evaporator": ("evaporator", "hot_side_outlet"),
    "Stack Gas": ("economizer", "hot_side_outlet"),
}
WATER_STREAMS = {
    "Feedwater": ("feedwater", "outlet"),
    "Pumped Water": ("feed_pump", "outlet"),
    "Economizer Water Out": ("economizer", "cold_side_outlet"),
    "Evaporator Outlet": ("evaporator", "cold_side_outlet"),
    "Drum Steam": ("drum", "vap_outlet"),
    "Blowdown": ("drum", "liq_outlet"),
    "Superheated Steam": ("superheater", "cold_side_outlet"),
    "Turbine Exhaust Steam": ("steam_turbine", "outlet"),
    "Condensate": ("condenser", "outlet"),
}


def _state(m, unit, port):
    """State block behind a unit port."""
    u = getattr(m.fs, unit)
    t = m.fs.time.first()
    if unit in ("air", "fuel", "feedwater"):
        return u.properties[t]
    if unit == "drum":
        return (u.vap_state if port == "vap_outlet" else u.liq_state)[t]
    if port.startswith("hot_side"):
        return u.hot_side.properties_out[t]
    if port.startswith("cold_side"):
        return u.cold_side.properties_out[t]
    return u.control_volume.properties_out[t]


def stream_table(m):
    """Temperature (C), pressure (bar), mass flow (kg/s) and vapor fraction of the DWSIM streams."""
    rows = {}
    for tag, (unit, port) in {**GAS_STREAMS, **WATER_STREAMS}.items():
        s = _state(m, unit, port)
        if tag in GAS_STREAMS:
            mw = sum(value(s.mole_frac_comp[j] * m.fs.gas.get_component(j).mw) for j in COMPONENTS)
            flow, vf = value(s.flow_mol) * mw, 1.0
        else:
            flow, vf = value(s.flow_mol) * MW_WATER, value(s.vapor_frac)
        rows[tag] = {"Temperature (C)": value(s.temperature) - 273.15, "Pressure (bar)": value(s.pressure) / 1e5,
                     "Mass flow (kg/s)": flow, "Vapor fraction": vf}
    return pd.DataFrame(rows).T


def raw_results(m):
    """Stream conditions and powers (kW) in the same structure as read_dwsim_results()."""
    t = m.fs.time.first()
    streams = stream_table(m)
    power = {
        "Compressor Power": value(m.fs.compressor.work_mechanical[t]) / 1e3,
        "GT Power": -value(m.fs.gas_turbine.work_mechanical[t]) / 1e3,
        "ST Power": -value(m.fs.steam_turbine.work_mechanical[t]) / 1e3,
        "Pump Power": value(m.fs.feed_pump.work_mechanical[t]) / 1e3,
        "Condenser Duty": -value(m.fs.condenser.heat_duty[t]) / 1e3,
    }
    duties = {name.capitalize(): value(getattr(m.fs, name).heat_duty[t]) / 1e3 for name in HX_NAMES}
    return {"streams": streams, "power": power, "duties": duties,
            "fuel_lhv": fuel_lhv(m, DESIGN["fuel"]), "fuel_carbon": fuel_carbon_ratio(DESIGN["fuel"])}


def fuel_carbon_ratio(fuel):
    """kg CO2 per kg of fuel for complete combustion."""
    carbon = {"methane": 1, "ethane": 2, "propane": 3}
    mw = {"methane": 16.043, "ethane": 30.069, "propane": 44.096}
    return 44.0095 * sum(x * carbon[j] for j, x in fuel.items()) / sum(x * mw[j] for j, x in fuel.items())


def kpis(raw):
    """Plant performance from raw stream/power results (same definitions for IDAES and DWSIM)."""
    s, p = raw["streams"], raw["power"]
    T = lambda tag: s.loc[tag, "Temperature (C)"]  # noqa: E731
    q_fuel = s.loc["Fuel", "Mass flow (kg/s)"] * raw["fuel_lhv"]
    gt_net = p["GT Power"] - p["Compressor Power"]
    st_net = p["ST Power"] - p["Pump Power"]
    net = gt_net + st_net
    drum_T = T("Drum Steam")
    co2 = s.loc["Fuel", "Mass flow (kg/s)"] * raw["fuel_carbon"]
    return {
        "Fuel LHV input (kW)": q_fuel,
        "Air flow (kg/s)": s.loc["Air", "Mass flow (kg/s)"],
        "Fuel flow (kg/s)": s.loc["Fuel", "Mass flow (kg/s)"],
        "Turbine inlet temperature (C)": T("Hot Gas"),
        "GT exhaust temperature (C)": T("GT Exhaust"),
        "Gas turbine gross (kW)": p["GT Power"],
        "Air compressor (kW)": p["Compressor Power"],
        "Gas turbine net (kW)": gt_net,
        "Steam turbine (kW)": p["ST Power"],
        "Feed pump (kW)": p["Pump Power"],
        "Steam cycle net (kW)": st_net,
        "Plant net power (kW)": net,
        "GT simple cycle efficiency LHV (%)": 100 * gt_net / q_fuel,
        "Combined cycle efficiency LHV (%)": 100 * net / q_fuel,
        "Heat rate LHV (kJ/kWh)": 3600 * q_fuel / net,
        "Steam flow (kg/s)": s.loc["Superheated Steam", "Mass flow (kg/s)"],
        "Steam temperature (C)": T("Superheated Steam"),
        "Drum saturation temperature (C)": drum_T,
        "HRSG pinch (K)": T("Gas After Evaporator") - drum_T,
        "Economizer approach (K)": drum_T - T("Economizer Water Out"),
        "Stack temperature (C)": T("Stack Gas"),
        "Steam turbine exhaust quality": s.loc["Turbine Exhaust Steam", "Vapor fraction"],
        "Condenser duty (kW)": p["Condenser Duty"],
        "CO2 emissions (kg/MWh)": co2 * 3600 / (net / 1000),
    }


def ccgt_results(m):
    """Plant performance of the solved IDAES model."""
    return kpis(raw_results(m))


def hrsg_areas(m):
    t = m.fs.time.first()
    return pd.DataFrame({name.capitalize(): {
        "Duty (kW)": value(getattr(m.fs, name).heat_duty[t]) / 1e3,
        "U (W/m2/K)": value(getattr(m.fs, name).overall_heat_transfer_coefficient[t]),
        "Area (m2)": value(getattr(m.fs, name).area),
        "UA (kW/K)": value(getattr(m.fs, name).area * getattr(m.fs, name).overall_heat_transfer_coefficient[t]) / 1e3,
        "LMTD (K)": value(getattr(m.fs, name).delta_temperature[t]),
    } for name in HX_NAMES}).T


def energy_balance(m):
    """Gas side, water side and overall energy balance closure (kW).

    The IDAES gas enthalpies include heats of formation, so the combustion heat
    release appears through the stream enthalpies themselves.
    """
    t = m.fs.time.first()
    W = lambda u: value(getattr(m.fs, u).work_mechanical[t]) / 1e3  # noqa: E731
    H = lambda sb: value(sb.flow_mol * sb.enth_mol) / 1e3  # noqa: E731
    q_hrsg = sum(value(getattr(m.fs, n).heat_duty[t]) for n in HX_NAMES) / 1e3
    gas_in = H(m.fs.air.properties[t]) + H(m.fs.fuel.properties[t]) + W("compressor")
    gas_out = H(m.fs.economizer.hot_side.properties_out[t]) - W("gas_turbine") + q_hrsg
    water_in = H(m.fs.feedwater.properties[t]) + W("feed_pump") + q_hrsg
    water_out = (H(m.fs.condenser.control_volume.properties_out[t]) + H(m.fs.drum.liq_state[t])
                 - W("steam_turbine") - value(m.fs.condenser.heat_duty[t]) / 1e3)
    return {"Gas side imbalance (kW)": gas_in - gas_out, "Water side imbalance (kW)": water_in - water_out}


def hrsg_tq(m):
    """Temperature vs cumulative heat duty (from the stack end) for the HRSG gas and water.

    The gas temperatures are at the section boundaries; the water profile has the
    economizer, the subcooled entry of the evaporator up to saturation, the
    boiling plateau, and the superheater.
    """
    t = m.fs.time.first()
    Tc = lambda sb: value(sb.temperature) - 273.15  # noqa: E731
    q = {n: value(getattr(m.fs, n).heat_duty[t]) / 1e3 for n in HX_NAMES}
    eco, ev, sh = m.fs.economizer, m.fs.evaporator, m.fs.superheater
    Q1 = q["economizer"]
    Q2 = Q1 + q["evaporator"]
    Q3 = Q2 + q["superheater"]
    gas = pd.DataFrame({"Q (kW)": [0, Q1, Q2, Q3],
                        "T (C)": [Tc(eco.hot_side.properties_out[t]), Tc(ev.hot_side.properties_out[t]),
                                  Tc(sh.hot_side.properties_out[t]), Tc(sh.hot_side.properties_in[t])]})
    # Heat to bring the evaporator inlet water to saturation (flow x enthalpy difference)
    ev_in, ev_out = ev.cold_side.properties_in[t], ev.cold_side.properties_out[t]
    h_sat = value(m.fs.water.htpx(p=value(ev_out.pressure) * pyunits.Pa, x=0))
    q_sub = max(value(ev_in.flow_mol) * (h_sat - value(ev_in.enth_mol)) / 1e3, 0.0)
    water = pd.DataFrame({"Q (kW)": [0, Q1, Q1 + q_sub, Q2, Q3],
                          "T (C)": [Tc(eco.cold_side.properties_in[t]), Tc(eco.cold_side.properties_out[t]),
                                    Tc(ev_out), Tc(ev_out), Tc(sh.cold_side.properties_out[t])]})
    return gas, water


def plot_tq(ax, m, label="", style="-"):
    """Plot the HRSG T-Q diagram on a matplotlib axis, marking the pinch."""
    gas, water = hrsg_tq(m)
    ax.plot(gas["Q (kW)"] / 1e3, gas["T (C)"], "o" + style, color="tab:red", label="Flue gas " + label)
    ax.plot(water["Q (kW)"] / 1e3, water["T (C)"], "s" + style, color="tab:blue", label="Water / steam " + label)
    q_pinch, T_gas, T_water = gas["Q (kW)"][1] / 1e3, gas["T (C)"][1], water["T (C)"][2]
    ax.annotate("pinch {0:.1f} K".format(T_gas - T_water), xy=(q_pinch, (T_gas + T_water) / 2),
                xytext=(q_pinch - 6, (T_gas + T_water) / 2 + 80), arrowprops={"arrowstyle": "->"})
    ax.set_xlabel("Heat transferred from the stack end (MW)")
    ax.set_ylabel("Temperature (C)")
    ax.grid(alpha=0.3)
    ax.legend()


# ---------------------------------------------------------------------------
# DWSIM comparison
# ---------------------------------------------------------------------------
def read_dwsim_results(path=DWSIM_FILE):
    """Read the solved CCGT from the DWSIM .dwxmz file (a zip archive holding the flowsheet as XML)."""
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([-+0-9.eE]+)<".format(tag), text).group(1))

    tags = dict(re.findall(r"<Name>((?:MAT|EN)-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    blocks = {}
    for kind in ("MAT", "EN"):
        for block in xml.split("<ComponentName>{0}-".format(kind))[1:]:
            tag = tags.get(kind + "-" + block[:36])
            if tag:
                blocks[tag] = block

    rows = {}
    for tag in list(GAS_STREAMS) + list(WATER_STREAMS):
        mix = blocks[tag][blocks[tag].find("<ComponentName>Mixture<"):]
        vap = blocks[tag][blocks[tag].find("<ComponentName>Vapor<"):]
        rows[tag] = {"Temperature (C)": first("temperature", mix) - 273.15, "Pressure (bar)": first("pressure", mix) / 1e5,
                     "Mass flow (kg/s)": first("massflow", mix),
                     "Vapor fraction": 1.0 if tag in GAS_STREAMS else first("molarfraction", vap)}
    power = {tag: first("EnergyFlow", blocks[tag]) for tag in
             ("Compressor Power", "GT Power", "ST Power", "Pump Power", "Condenser Duty")}

    def h(tag):
        mix = blocks[tag][blocks[tag].find("<ComponentName>Mixture<"):]
        return first("enthalpy", mix)
    m_ = lambda tag: rows[tag]["Mass flow (kg/s)"]  # noqa: E731
    duties = {
        "Superheater": m_("Superheated Steam") * (h("Superheated Steam") - h("Drum Steam")),
        "Evaporator": m_("Evaporator Outlet") * (h("Evaporator Outlet") - h("Economizer Water Out")),
        "Economizer": m_("Economizer Water Out") * (h("Economizer Water Out") - h("Pumped Water")),
    }
    fuel = DESIGN["fuel"]
    mw = {"methane": 16.043, "ethane": 30.069, "propane": 44.096}
    lhv = sum(x * LHV_MOLAR[j] for j, x in fuel.items()) / sum(x * mw[j] for j, x in fuel.items()) * 1000
    return {"streams": pd.DataFrame(rows).T, "power": power, "duties": duties,
            "fuel_lhv": lhv, "fuel_carbon": fuel_carbon_ratio(fuel)}


def compare_with_dwsim(m, dwsim):
    """Side by side comparison of the IDAES and DWSIM results.

    Returns (streams, performance) DataFrames; differences are IDAES - DWSIM.
    """
    idaes_raw = raw_results(m)
    rows = []
    for tag in dwsim["streams"].index:
        for col in ("Temperature (C)", "Pressure (bar)", "Mass flow (kg/s)"):
            d, i = dwsim["streams"].loc[tag, col], idaes_raw["streams"].loc[tag, col]
            rows.append({"Stream": tag, "Property": col, "DWSIM": d, "IDAES": i, "Difference": i - d})
    streams = pd.DataFrame(rows)

    k_d, k_i = kpis(dwsim), kpis(idaes_raw)
    perf = pd.DataFrame({"DWSIM": k_d, "IDAES": k_i})
    for name in ("Superheater", "Evaporator", "Economizer"):
        perf.loc[name + " duty (kW)"] = [dwsim["duties"][name], idaes_raw["duties"][name]]
    perf["Difference"] = perf["IDAES"] - perf["DWSIM"]
    perf["Difference (%)"] = 100 * perf["Difference"] / perf["DWSIM"].abs()
    return streams, perf


if __name__ == "__main__":
    m, res = setup_ccgt_model()
    print(res.solver.termination_condition)
    pd.set_option("display.width", 200)
    print(stream_table(m).round(3))
    print(hrsg_areas(m).round(2))
    streams, perf = compare_with_dwsim(m, read_dwsim_results())
    print(perf.round(3))
    print(energy_balance(m))
