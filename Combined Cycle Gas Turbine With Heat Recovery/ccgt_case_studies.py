"""Case studies for the combined cycle gas turbine (see ccgt.py for the model).

1. Ambient temperature: off-design derating of power and efficiency from 0 to 45 C
2. Part load: efficiency vs load with turbine inlet temperature control or inlet guide vanes
3. Economics: levelized cost of electricity, fuel and carbon price sensitivity, simple vs combined cycle
4. Design optimization: drum pressure and steam temperature, and the HRSG pinch trade-off
"""

import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
from pyomo.environ import Objective, Constraint, value, maximize, units as pyunits
from ccgt import DESIGN, HX_NAMES, MW_WATER, ccgt_results, set_mode, get_ccgt_solver, saturation_T

SOLVER = get_ccgt_solver(max_iter=100)  # a failed step is halved rather than iterated for long
MMBTU = 1.055056e6  # kJ per MMBtu


def try_solve(m):
    """Solve the model; keep the previous solution if the solve fails. Returns True on success."""
    result = SOLVER.solve(m, load_solutions=False)
    if str(result.solver.termination_condition) == "optimal":
        m.solutions.load_from(result)
        return True
    return False


def solve(m):
    if not try_solve(m):
        raise RuntimeError("CCGT model did not converge")


def step_to(m, apply, start, end, depth=0, max_depth=8):
    """Move a model input from start to end, halving the step whenever a solve fails.

    apply(m, x) sets the input to x. Large changes are reached as a series of
    smaller, warm started steps (continuation).
    """
    apply(m, end)
    if try_solve(m):
        return
    if depth >= max_depth or start == end:
        raise RuntimeError("CCGT model did not converge stepping to {0}".format(end))
    mid = (start + end) / 2
    step_to(m, apply, start, mid, depth + 1, max_depth)
    step_to(m, apply, mid, end, depth + 1, max_depth)


def sweep(m, name, values, apply, start, extra=None):
    """Solve the model for each value (warm started from the previous one) and tabulate the results.

    start is the current value of the input; extra(m) may add columns to each row.
    """
    rows = []
    current = start
    for v in values:
        step_to(m, apply, current, v)
        current = v
        row = {name: v}
        row.update(ccgt_results(m))
        if extra is not None:
            row.update(extra(m))
        rows.append(row)
    return pd.DataFrame(rows).set_index(name)


def two_way_sweep(m, name, values, apply, start, extra=None):
    """Sweep the values below and above start (each walked away from start and back), sorted by value."""
    parts = []
    for side in ([v for v in values if v <= start], [v for v in values if v > start]):
        if not side:
            continue
        parts.append(sweep(m, name, side, apply, start, extra))
        # Walk back to the start through the same points (warm started, no large jumps)
        back = list(reversed(side[:-1])) + [start]
        current = side[-1]
        for v in back:
            step_to(m, apply, current, v)
            current = v
    return pd.concat(parts).sort_index()


def hrsg_extra(m):
    """HRSG areas and the gas turbine pressure ratio for the sweep tables."""
    t = m.fs.time.first()
    row = {"{0} area (m2)".format(n.capitalize()): value(getattr(m.fs, n).area) for n in HX_NAMES}
    row["HRSG area (m2)"] = sum(value(getattr(m.fs, n).area) for n in HX_NAMES)
    row["Pressure ratio"] = value(m.fs.compressor.outlet.pressure[t] / m.fs.air.pressure[t])
    return row


# Condenser saturation temperature above ambient that reproduces the design condenser pressure (0.08 bar) at 15 C
CONDENSER_DT = saturation_T(DESIGN["condenser_P"]) - 273.15 - DESIGN["ambient_T"]


def set_condenser_pressure(m, P):
    """Condenser pressure in Pa: steam turbine back pressure, condensate and (open loop) feedwater state."""
    t = m.fs.time.first()
    m.fs.steam_turbine.outlet.pressure.fix(P)
    m.fs.condenser.outlet.enth_mol.fix(value(m.fs.water.htpx(p=P * pyunits.Pa, x=0)))
    m.fs.feedwater.pressure.fix(P)
    m.fs.feedwater.enth_mol.fix(value(m.fs.water.htpx(T=(saturation_T(P) - 0.5) * pyunits.K, p=P * pyunits.Pa)))


def set_ambient(m, T, condenser_follows_ambient=False):
    """Ambient (compressor inlet) temperature in C.

    With condenser_follows_ambient the condenser saturation temperature is held CONDENSER_DT above
    ambient (a cooling system with a constant temperature rise), otherwise the condenser pressure is fixed.
    """
    m.fs.air.temperature.fix(T + 273.15)
    if condenser_follows_ambient:
        set_condenser_pressure(m, PropsSI("P", "T", T + 273.15 + CONDENSER_DT, "Q", 0, "Water"))


# ---------------------------------------------------------------------------
# 1. Ambient temperature
# ---------------------------------------------------------------------------
def ambient_sweep(m, temperatures=(15, 10, 5, 0, 20, 25, 30, 35, 40, 45), condenser_follows_ambient=False):
    """Off-design performance vs ambient temperature at full load (design TIT).

    The model must be at the design point in "offdesign" mode. The sweep walks
    away from 15 C in both directions; the result is sorted by temperature.
    """
    t = m.fs.time.first()
    apply = lambda m, T: set_ambient(m, T, condenser_follows_ambient)  # noqa: E731

    def extra(m):
        row = hrsg_extra(m)
        row["Condenser pressure (bar)"] = value(m.fs.steam_turbine.outlet.pressure[t]) / 1e5
        return row
    return two_way_sweep(m, "Ambient temperature (C)", temperatures, apply, DESIGN["ambient_T"], extra)


def seasonal_output(ambient_df, monthly_T, hours_per_month=730.0, availability=1.0):
    """Monthly energy (MWh) and average efficiency from the ambient curve and monthly mean temperatures (C)."""
    T = ambient_df.index.values
    rows = []
    for month, Tm in monthly_T.items():
        power = np.interp(Tm, T, ambient_df["Plant net power (kW)"])
        eff = np.interp(Tm, T, ambient_df["Combined cycle efficiency LHV (%)"])
        rows.append({"Month": month, "Ambient (C)": Tm, "Net power (MW)": power / 1e3,
                     "Efficiency (%)": eff, "Energy (MWh)": power / 1e3 * hours_per_month * availability})
    return pd.DataFrame(rows).set_index("Month")


# ---------------------------------------------------------------------------
# 2. Part load
# ---------------------------------------------------------------------------
def _fuel_flow(m):
    return m.fs.fuel.flow_mol[m.fs.time.first()]


def set_fuel_control(m):
    """Load set by the fuel flow (the turbine inlet temperature becomes a result)."""
    m.fs.TIT_spec.deactivate()
    _fuel_flow(m).fix()


def part_load_tit(m, fractions=(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)):
    """Part load by reducing the fuel flow at a fixed compressor flow (no inlet guide vanes).

    The turbine inlet and exhaust temperatures fall with load. fractions are of the design fuel flow.
    """
    set_fuel_control(m)
    f0 = value(_fuel_flow(m))
    apply = lambda m, x: _fuel_flow(m).fix(x * f0)  # noqa: E731
    return _add_load(two_way_sweep(m, "Fuel flow fraction", fractions, apply, 1.0, hrsg_extra).sort_index(ascending=False))


def part_load_igv(m, fractions=(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4), min_air=0.70):
    """Part load with inlet guide vanes: the air flow is reduced to hold the GT exhaust temperature
    (maximizing the steam cycle output), down to min_air of the design volumetric flow. Below that the
    guide vanes are at their limit and the turbine inlet temperature falls.
    """
    t = m.fs.time.first()
    set_fuel_control(m)
    f0 = value(_fuel_flow(m))
    v0 = value(m.fs.air_volumetric_flow)
    exhaust_T = m.fs.gas_turbine.control_volume.properties_out[t].temperature
    T_exh = value(exhaust_T)

    def igv_mode():
        m.fs.air_volumetric_flow.unfix()
        exhaust_T.fix(T_exh)

    def limit_mode():
        exhaust_T.unfix()
        m.fs.air_volumetric_flow.fix(min_air * v0)

    apply = lambda m, x: _fuel_flow(m).fix(x * f0)  # noqa: E731
    igv_mode()
    rows, current, at_limit = [], 1.0, False
    for x in fractions:
        step_to(m, apply, current, x)
        current = x
        if not at_limit and value(m.fs.air_volumetric_flow) < min_air * v0:
            # Guide vanes at their limit: walk the air flow back to the limit at this fuel flow
            start = value(m.fs.air_volumetric_flow) / v0
            limit_mode()
            step_to(m, lambda m, a: m.fs.air_volumetric_flow.fix(a * v0), start, min_air)
            at_limit = True
        row = {"Fuel flow fraction": x}
        row.update(ccgt_results(m))
        row.update(hrsg_extra(m))
        row["Air flow fraction"] = value(m.fs.air_volumetric_flow) / v0
        rows.append(row)

    # Back to the design point along the same path: fuel flow up at the guide vane limit, then
    # guide vanes open at constant exhaust temperature
    limit_reached_at = next((r["Fuel flow fraction"] for r in rows if r["Air flow fraction"] <= min_air + 1e-6), None)
    for x in reversed([r["Fuel flow fraction"] for r in rows]):
        if at_limit and limit_reached_at is not None and x > limit_reached_at:
            # Leave the limit when the exhaust temperature is back at its setpoint
            exhaust_T.fix(value(exhaust_T))
            m.fs.air_volumetric_flow.unfix()
            step_to(m, lambda m, T: exhaust_T.fix(T), value(exhaust_T), T_exh)
            at_limit = False
        step_to(m, apply, current, x)
        current = x
    exhaust_T.unfix()
    m.fs.air_volumetric_flow.fix(v0)
    solve(m)
    return _add_load(pd.DataFrame(rows).set_index("Fuel flow fraction"))


def _add_load(df):
    df["Load (%)"] = 100 * df["Plant net power (kW)"] / df["Plant net power (kW)"].iloc[0]
    return df


def restore_tit_control(m, TIT=DESIGN["TIT"]):
    """Back to load control by the turbine inlet temperature."""
    t = m.fs.time.first()
    m.fs.TIT.fix(value(m.fs.combustor.control_volume.properties_out[t].temperature))
    m.fs.TIT_spec.activate()
    _fuel_flow(m).unfix()
    step_to(m, lambda m, T: m.fs.TIT.fix(T + 273.15), value(m.fs.TIT) - 273.15, TIT)


# ---------------------------------------------------------------------------
# 3. Economics
# ---------------------------------------------------------------------------
ECONOMICS = {
    # Combined cycle (CC) and simple cycle (SC) gas turbine plant, small industrial size
    "capex_cc": 1100.0,  # $/kW net
    "capex_sc": 750.0,  # $/kW net
    "fixed_om_cc": 18.0,  # $/kW-yr
    "fixed_om_sc": 12.0,  # $/kW-yr
    "var_om_cc": 3.5,  # $/MWh
    "var_om_sc": 5.0,  # $/MWh
    "discount_rate": 0.08,
    "lifetime": 25,  # years
    "fuel_price": 6.0,  # $/MMBtu (LHV)
    "capacity_factor": 0.60,
    "carbon_price": 0.0,  # $/t CO2
}


def crf(rate, years):
    """Capital recovery factor."""
    return rate * (1 + rate) ** years / ((1 + rate) ** years - 1)


def plant_cases(kpi):
    """Net power, heat rate and CO2 intensity of the combined cycle and the gas turbine alone (simple cycle)."""
    gt_net = kpi["Gas turbine net (kW)"]
    return {
        "Combined cycle": {"Net power (kW)": kpi["Plant net power (kW)"], "Heat rate (kJ/kWh)": kpi["Heat rate LHV (kJ/kWh)"],
                           "CO2 (kg/MWh)": kpi["CO2 emissions (kg/MWh)"]},
        "Simple cycle": {"Net power (kW)": gt_net, "Heat rate (kJ/kWh)": 3600 * kpi["Fuel LHV input (kW)"] / gt_net,
                         "CO2 (kg/MWh)": kpi["CO2 emissions (kg/MWh)"] * kpi["Plant net power (kW)"] / gt_net},
    }


def lcoe(heat_rate, co2, capex, fixed_om, var_om, econ=ECONOMICS, **override):
    """Levelized cost of electricity breakdown in $/MWh.

    heat_rate in kJ/kWh (LHV), co2 in kg/MWh, capex in $/kW, fixed_om in $/kW-yr, var_om in $/MWh.
    Any economics key (fuel_price, capacity_factor, carbon_price, ...) can be overridden.
    """
    e = dict(econ, **override)
    mwh_per_kw = 8.760 * e["capacity_factor"]  # MWh per kW-yr
    parts = {
        "Capital": capex * crf(e["discount_rate"], e["lifetime"]) / mwh_per_kw,
        "Fixed O&M": fixed_om / mwh_per_kw,
        "Fuel": heat_rate * 1000 / MMBTU * e["fuel_price"],
        "Variable O&M": var_om,
        "Carbon": co2 / 1000 * e["carbon_price"],
    }
    parts["LCOE"] = sum(parts.values())
    return parts


def lcoe_table(kpi, econ=ECONOMICS, **override):
    """LCOE breakdown of the combined cycle and the simple cycle gas turbine."""
    cases = plant_cases(kpi)
    out = {}
    for name, key in (("Combined cycle", "cc"), ("Simple cycle", "sc")):
        c = cases[name]
        out[name] = lcoe(c["Heat rate (kJ/kWh)"], c["CO2 (kg/MWh)"], econ["capex_" + key], econ["fixed_om_" + key],
                         econ["var_om_" + key], econ, **override)
    return pd.DataFrame(out)


def lcoe_heatmap(kpi, fuel_prices=(3, 4, 6, 8, 10, 12), capacity_factors=(0.3, 0.5, 0.7, 0.9), econ=ECONOMICS):
    """Combined cycle LCOE ($/MWh) over fuel price and capacity factor."""
    c = plant_cases(kpi)["Combined cycle"]
    return pd.DataFrame(
        {cf: [lcoe(c["Heat rate (kJ/kWh)"], c["CO2 (kg/MWh)"], econ["capex_cc"], econ["fixed_om_cc"], econ["var_om_cc"],
                   econ, fuel_price=fp, capacity_factor=cf)["LCOE"] for fp in fuel_prices] for cf in capacity_factors},
        index=pd.Index(fuel_prices, name="Fuel price ($/MMBtu)")).rename_axis(columns="Capacity factor")


def tornado(kpi, econ=ECONOMICS, swing=0.25):
    """Change of the combined cycle LCOE for +/- swing (fraction) on each input."""
    c = plant_cases(kpi)["Combined cycle"]
    base_inputs = {"Heat rate": c["Heat rate (kJ/kWh)"], "CAPEX": econ["capex_cc"]}

    def f(hr=None, capex=None, **ov):
        return lcoe(hr or c["Heat rate (kJ/kWh)"], c["CO2 (kg/MWh)"], capex or econ["capex_cc"], econ["fixed_om_cc"],
                    econ["var_om_cc"], econ, **ov)["LCOE"]
    base = f()
    rows = {}
    for name, fn in (
            ("Fuel price", lambda k: f(fuel_price=econ["fuel_price"] * k)),
            ("Capacity factor", lambda k: f(capacity_factor=min(econ["capacity_factor"] * k, 1.0))),
            ("CAPEX", lambda k: f(capex=base_inputs["CAPEX"] * k)),
            ("Discount rate", lambda k: f(discount_rate=econ["discount_rate"] * k)),
            ("Heat rate", lambda k: f(hr=base_inputs["Heat rate"] * k)),
            ("Carbon price (+25 $/t)", lambda k: f(carbon_price=econ["carbon_price"] + (25 if k > 1 else 0)))):
        rows[name] = {"Low": fn(1 - swing) - base, "High": fn(1 + swing) - base}
    df = pd.DataFrame(rows).T
    return df.reindex(df.abs().max(axis=1).sort_values().index), base


def carbon_price_study(kpi, carbon_prices=(0, 25, 50, 75, 100, 150), econ=ECONOMICS):
    """LCOE of the combined and simple cycle vs carbon price."""
    rows = {cp: lcoe_table(kpi, econ, carbon_price=cp).loc["LCOE"] for cp in carbon_prices}
    return pd.DataFrame(rows).T.rename_axis("Carbon price ($/t CO2)")


def dispatch_study(kpi, power_prices=np.arange(20, 121, 10), econ=ECONOMICS, **override):
    """Spark spread and short run margin vs the electricity price ($/MWh).

    The plant runs when the power price exceeds its marginal cost (fuel + variable O&M + carbon).
    """
    rows = []
    t = lcoe_table(kpi, econ, **override)
    for name in t.columns:
        marginal = t.loc[["Fuel", "Variable O&M", "Carbon"], name].sum()
        for p in power_prices:
            rows.append({"Plant": name, "Power price ($/MWh)": p, "Spark spread ($/MWh)": p - t.loc["Fuel", name],
                         "Margin ($/MWh)": p - marginal, "Dispatched": p > marginal})
    return pd.DataFrame(rows)


def cc_vs_sc_savings(kpi, econ=ECONOMICS, hours=None):
    """Annual fuel, cost and CO2 savings of the combined cycle against generating the same energy in simple cycle."""
    hours = 8760 * econ["capacity_factor"] if hours is None else hours
    cases = plant_cases(kpi)
    energy = cases["Combined cycle"]["Net power (kW)"] / 1e3 * hours  # MWh
    fuel = {k: c["Heat rate (kJ/kWh)"] * energy * 1000 / MMBTU for k, c in cases.items()}  # MMBtu
    co2 = {k: c["CO2 (kg/MWh)"] * energy / 1000 for k, c in cases.items()}  # t
    return {
        "Annual energy (MWh)": energy,
        "Fuel saved (MMBtu/yr)": fuel["Simple cycle"] - fuel["Combined cycle"],
        "Fuel cost saved ($/yr)": (fuel["Simple cycle"] - fuel["Combined cycle"]) * econ["fuel_price"],
        "CO2 avoided (t/yr)": co2["Simple cycle"] - co2["Combined cycle"],
        "Extra CAPEX of the steam cycle ($)": econ["capex_cc"] * cases["Combined cycle"]["Net power (kW)"]
        - econ["capex_sc"] * cases["Simple cycle"]["Net power (kW)"],
    }


# ---------------------------------------------------------------------------
# 4. Design optimization
# ---------------------------------------------------------------------------
def set_drum_pressure(m, p_bar):
    """Feed pump outlet pressure in bar (the drum is 1 bar lower, the superheater outlet 2 bar lower)."""
    m.fs.feed_pump.outlet.pressure.fix(p_bar * 1e5)


def drum_pressure_sweep(m, pressures=(62, 50, 40, 30, 20, 70, 80, 90, 100, 120)):
    """Steam cycle design vs pressure, in "targets" mode (fixed TIT, pinch, approach and steam temperature).

    The steam flow follows the pinch; the HRSG areas are results.
    """
    set_mode(m, "targets")
    solve(m)
    return two_way_sweep(m, "Pump outlet pressure (bar)", pressures, set_drum_pressure, DESIGN["pump_outlet_P"] / 1e5,
                         hrsg_extra)


def optimize_steam_conditions(m, max_steam_T=540.0, min_exhaust_quality=0.88, min_sh_approach=25.0,
                              pressure_bounds=(20e5, 140e5)):
    """Maximize the net power over the pump outlet pressure and the superheated steam temperature.

    Constraints: steam turbine exhaust quality (blade erosion), steam temperature below a material
    limit and at least min_sh_approach below the GT exhaust. TIT, pinch and approach stay fixed.
    Returns the optimum results; the model is left at the optimum.
    """
    t = m.fs.time.first()
    set_mode(m, "targets")
    solve(m)
    p = m.fs.feed_pump.outlet.pressure[t]
    p.unfix()
    p.setlb(pressure_bounds[0])
    p.setub(pressure_bounds[1])
    m.fs.steam_T.unfix()
    m.fs.steam_T.setub(max_steam_T + 273.15)
    st_out = m.fs.steam_turbine.control_volume.properties_out[t]
    gt_out = m.fs.gas_turbine.control_volume.properties_out[t]
    m.fs.opt_quality = Constraint(expr=st_out.vapor_frac >= min_exhaust_quality)
    m.fs.opt_sh_approach = Constraint(expr=m.fs.steam_T <= gt_out.temperature - min_sh_approach)
    m.fs.opt_objective = Objective(expr=m.fs.net_power / 1e6, sense=maximize)
    ok = try_solve(m)
    result = ccgt_results(m)
    result["Pump outlet pressure (bar)"] = value(p) / 1e5
    result["Converged"] = ok
    # Clean up: leave the optimum's inputs fixed
    m.fs.del_component(m.fs.opt_objective)
    m.fs.del_component(m.fs.opt_quality)
    m.fs.del_component(m.fs.opt_sh_approach)
    p.fix()
    m.fs.steam_T.fix()
    m.fs.steam_T.setub(None)
    return result


def restore_steam_design(m, design=DESIGN):
    """Return the steam conditions to the design values (targets mode)."""
    t = m.fs.time.first()
    p0 = value(m.fs.feed_pump.outlet.pressure[t]) / 1e5
    T0 = value(m.fs.steam_T) - 273.15
    step_to(m, lambda m, T: m.fs.steam_T.fix(T + 273.15), T0, design["steam_T"])
    step_to(m, set_drum_pressure, p0, design["pump_outlet_P"] / 1e5)


def pinch_study(m, pinches=(10, 8, 6, 5, 4, 3, 12, 15, 20, 25, 30, 40)):
    """HRSG performance and surface vs the pinch, in "targets" mode (fixed TIT, approach and steam conditions).

    A smaller pinch makes more steam (power) but needs more heat transfer surface; see pinch_economics().
    """
    set_mode(m, "targets")
    solve(m)
    apply = lambda m, x: m.fs.pinch.fix(x)  # noqa: E731
    return two_way_sweep(m, "Pinch (K)", pinches, apply, DESIGN["pinch"], hrsg_extra)


def pinch_economics(pinch_df, area_cost=600.0, power_price=60.0, hours=8760 * 0.6, econ=ECONOMICS):
    """Annual net benefit of each pinch relative to the design pinch ($/yr).

    Extra energy x power price - annualized extra HRSG cost. area_cost is the installed HRSG cost per m2 of the
    heat transfer area used in the model (U of 50-60 W/m2/K); 600 $/m2 puts the design HRSG at about $4M.
    """
    ref = pinch_df.loc[DESIGN["pinch"]]
    annuity = crf(econ["discount_rate"], econ["lifetime"])
    out = pd.DataFrame(index=pinch_df.index)
    out["Extra revenue ($/yr)"] = (pinch_df["Plant net power (kW)"] - ref["Plant net power (kW)"]) / 1e3 * hours * power_price
    out["Extra HRSG cost ($/yr)"] = (pinch_df["HRSG area (m2)"] - ref["HRSG area (m2)"]) * area_cost * annuity
    out["Net benefit ($/yr)"] = out["Extra revenue ($/yr)"] - out["Extra HRSG cost ($/yr)"]
    return out


def optimal_pinch(pinch_df, area_costs=(200, 400, 600, 800, 1000), **kwargs):
    """Pinch with the highest net benefit for a range of HRSG area costs (interpolated on a fine grid)."""
    fine = np.linspace(pinch_df.index.min(), pinch_df.index.max(), 400)
    rows = {}
    for c in area_costs:
        benefit = pinch_economics(pinch_df, area_cost=c, **kwargs)["Net benefit ($/yr)"]
        # Smooth interpolation of the benefit curve
        coeffs = np.polyfit(pinch_df.index.values, benefit.values, 4)
        b = np.polyval(coeffs, fine)
        rows[c] = {"Optimal pinch (K)": fine[b.argmax()], "Benefit vs 10 K design ($/yr)": b.max()}
    return pd.DataFrame(rows).T.rename_axis("HRSG area cost ($/m2)")
