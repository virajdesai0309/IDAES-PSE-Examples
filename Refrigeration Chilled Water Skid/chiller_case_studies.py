"""Case studies for the chiller skid digital twin (see chiller.py for the model).

1. Fouling sweep: effect of condenser fouling on pressures, power, capacity and COP
2. Energy, cost and CO2 penalty of fouling, and the optimal cleaning interval
3. Fouling soft sensor: estimating the fouling resistance from plant measurements
4. Operating optimization and part load performance
"""

import math

import numpy as np
import pandas as pd
from pyomo.environ import Objective, value, units as pyunits, minimize
from idaes.core.solvers import get_solver
from idaes.models.properties.iapws95 import htpx

from chiller import MW_WATER, DESIGN, chiller_results, set_operating_mode

SOLVER = get_solver(options={"max_iter": 300})


def try_solve(m):
    """Solve the model; keep the previous solution if the solve fails. Returns True on success."""
    result = SOLVER.solve(m, load_solutions=False)
    if str(result.solver.termination_condition) == "optimal":
        m.solutions.load_from(result)
        return True
    return False


def solve(m):
    if not try_solve(m):
        raise RuntimeError("Chiller model did not converge")


def step_to(m, apply, start, end, depth=0, max_depth=8):
    """Move a model input from start to end, halving the step whenever a solve fails.

    apply(m, x) sets the input to x. Large changes are reached as a series of
    smaller, warm started steps (continuation).
    """
    apply(m, end)
    if try_solve(m):
        return
    if depth >= max_depth or start == end:
        raise RuntimeError("Chiller model did not converge stepping to {0}".format(end))
    mid = (start + end) / 2
    step_to(m, apply, start, mid, depth + 1, max_depth)
    step_to(m, apply, mid, end, depth + 1, max_depth)


def set_fouling(m, condenser=None, evaporator=None):
    """Set the fouling resistances (m2.K/W)."""
    if condenser is not None:
        m.fs.condenser.fouling_resistance.fix(condenser)
    if evaporator is not None:
        m.fs.evaporator.fouling_resistance.fix(evaporator)


def set_water_inlet(feed, flow=None, temperature=None, pressure=DESIGN["water_pressure"]):
    """Change a water feed: flow in kg/s, temperature in C."""
    if flow is not None:
        feed.flow_mol.fix(flow / MW_WATER)
    if temperature is not None:
        feed.enth_mol.fix(value(htpx(T=(temperature + 273.15) * pyunits.K, P=pressure * pyunits.Pa)))


def restore_design(m):
    """Return to the clean design point (fixed refrigerant flow, design water inlets) and solve."""
    set_fouling(m, condenser=0, evaporator=0)
    set_water_inlet(m.fs.cooling_water_in, flow=DESIGN["cooling_water_flow"], temperature=DESIGN["cooling_water_inlet"])
    set_water_inlet(m.fs.chilled_water_in, flow=DESIGN["chilled_water_flow"], temperature=DESIGN["chilled_water_return"])
    set_operating_mode(m, "fixed_flow", refrigerant_flow=DESIGN["refrigerant_flow"])
    solve(m)


def sweep(m, name, values, apply, start=None):
    """Solve the model for each value and tabulate the results.

    Each point is warm started from the previous one; start is the current
    value of the input (defaults to the first value).
    """
    rows = []
    current = values[0] if start is None else start
    for v in values:
        step_to(m, apply, current, v)
        current = v
        row = {name: v}
        row.update(chiller_results(m))
        rows.append(row)
    return pd.DataFrame(rows).set_index(name)


# ---------------------------------------------------------------------------
# 1. Fouling sweep
# ---------------------------------------------------------------------------
def fouling_sweep(m, fouling_values, side="condenser", mode="fixed_flow", chilled_water_supply=7.0):
    """Effect of fouling on one exchanger.

    mode "fixed_flow": constant refrigerant flow, the cooling capacity drops.
    mode "fixed_supply": the chilled water supply is held, so the cooling load
    is met and the compressor works harder (the real fouling penalty).
    """
    def apply(m, rf):
        set_fouling(m, **{side: rf})

    set_operating_mode(m, mode, chilled_water_supply=chilled_water_supply)
    df = sweep(m, "Fouling resistance (m2K/W)", fouling_values, apply, start=0.0)
    # Return to the clean design point, so the next study starts from a solved state
    step_to(m, apply, fouling_values[-1], 0.0)
    restore_design(m)
    df["U (W/m2K)"] = [1 / (1 / DESIGN[side + "_U"] + rf) for rf in df.index]
    return df


# ---------------------------------------------------------------------------
# 2. Energy, cost and CO2 penalty, and the cleaning interval
# ---------------------------------------------------------------------------
def annual_penalty(sweep_df, hours=6000, tariff=8.0, co2_factor=0.7):
    """Annual extra energy, cost and CO2 caused by fouling, relative to the clean case.

    hours: operating hours per year; tariff: cost per kWh; co2_factor: kg CO2 per kWh.
    """
    W = sweep_df["Compressor power (kW)"]
    extra = W - W.iloc[0]
    return pd.DataFrame({
        "Compressor power (kW)": W,
        "Extra power (kW)": extra,
        "Extra energy (kWh/yr)": extra * hours,
        "Extra cost (per yr)": extra * hours * tariff,
        "Extra CO2 (t/yr)": extra * hours * co2_factor / 1000,
        "COP": sweep_df["COP"],
    })


def power_surrogate(sweep_df, degree=2):
    """Polynomial fit of compressor power (kW) vs fouling resistance from an IDAES sweep."""
    coeffs = np.polyfit(sweep_df.index.values, sweep_df["Compressor power (kW)"].values, degree)
    return np.poly1d(coeffs)


def fouling_growth(t_days, rf_max=5e-4, tau_days=60.0):
    """Asymptotic (Kern-Seaton) fouling growth: R_f(t) = R_f,max (1 - exp(-t / tau))."""
    return rf_max * (1 - np.exp(-np.asarray(t_days, dtype=float) / tau_days))


def cleaning_interval_study(power, intervals_days, rf_max=5e-4, tau_days=60.0, hours_per_day=24 * 6000 / 8760,
                            tariff=8.0, cleaning_cost=5000.0, operating_days=365):
    """Annual cost of fouling (extra energy) plus cleaning, for each cleaning interval.

    Fouling restarts from zero after each cleaning. power is the surrogate W(R_f) in kW.
    """
    W_clean = power(0.0)
    rows = []
    for interval in intervals_days:
        t = np.linspace(0, interval, 400)
        extra_kw = np.trapezoid(power(fouling_growth(t, rf_max, tau_days)) - W_clean, t) / interval  # cycle average
        energy_cost = extra_kw * hours_per_day * operating_days * tariff
        cleanings = operating_days / interval
        rows.append({
            "Cleaning interval (days)": interval,
            "Average extra power (kW)": extra_kw,
            "Energy penalty (per yr)": energy_cost,
            "Cleaning cost (per yr)": cleanings * cleaning_cost,
            "Total cost (per yr)": energy_cost + cleanings * cleaning_cost,
        })
    return pd.DataFrame(rows).set_index("Cleaning interval (days)")


# ---------------------------------------------------------------------------
# 3. Fouling soft sensor
# ---------------------------------------------------------------------------
MEASUREMENTS = {
    # name: (model expression, measurement noise standard deviation)
    "Condenser pressure (bar)": (lambda m: m.fs.compressor.control_volume.properties_out[0].pressure / 1e5, 0.02),
    "Discharge temperature (C)": (lambda m: m.fs.compressor.control_volume.properties_out[0].temperature - 273.15, 0.2),
    "Cooling water outlet (C)": (lambda m: m.fs.condenser.cold_side.properties_out[0].temperature - 273.15, 0.05),
    "Compressor power (kW)": (lambda m: m.fs.compressor.work_mechanical[0] / 1000, 0.03),
}


def synthetic_measurements(m, true_fouling, seed=0):
    """Plant data from the twin itself at known condenser fouling, plus measurement noise.

    true_fouling should be in increasing order (fouling grows between cleanings).
    """
    rng = np.random.default_rng(seed)
    rows = []
    current = 0.0

    def apply(m, rf):
        set_fouling(m, condenser=rf)

    for rf in true_fouling:
        step_to(m, apply, current, rf)
        current = rf
        row = {"True fouling (m2K/W)": rf, "True COP": chiller_results(m)["COP"]}
        for name, (expr, sigma) in MEASUREMENTS.items():
            row[name] = value(expr(m)) + rng.normal(0, sigma)
        rows.append(row)
    step_to(m, apply, current, 0.0)
    return pd.DataFrame(rows)


def estimate_fouling(m, measured):
    """Least-squares estimate of the condenser fouling resistance from one set of measurements.

    The fouling resistance is freed and chosen to minimize the noise-weighted
    squared differences between the model and the measurements.
    """
    rf = m.fs.condenser.fouling_resistance
    rf.unfix()
    rf.setlb(0)
    rf.setub(5e-3)
    m.fs.soft_sensor = Objective(
        expr=sum(((expr(m) - measured[name]) / sigma) ** 2 for name, (expr, sigma) in MEASUREMENTS.items()),
        sense=minimize,
    )
    try:
        solve(m)
        estimate = value(rf)
        cop = chiller_results(m)["COP"]
    finally:
        m.fs.del_component(m.fs.soft_sensor)
        rf.fix()
        rf.setlb(None)
        rf.setub(None)
    return estimate, cop


def soft_sensor_study(m, measurements, cop_clean, alert_cop_loss=0.05):
    """Estimate fouling for each measurement set and flag when the COP loss exceeds the threshold."""
    rows = []
    for _, meas in measurements.iterrows():
        rf_est, cop_est = estimate_fouling(m, meas)
        loss = 1 - cop_est / cop_clean
        rows.append({
            "True fouling (m2K/W)": meas["True fouling (m2K/W)"],
            "Estimated fouling (m2K/W)": rf_est,
            "True COP": meas["True COP"],
            "Estimated COP": cop_est,
            "COP loss (%)": loss * 100,
            "Cleaning alert": loss > alert_cop_loss,
        })
    step_to(m, lambda m, rf: set_fouling(m, condenser=rf), rows[-1]["Estimated fouling (m2K/W)"], 0.0)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. Operating optimization and part load
# ---------------------------------------------------------------------------
def setpoint_reset_study(m, supplies=(7, 8, 9, 10, 6), delta_t=5.0):
    """Chilled water supply setpoint at constant load (constant flow, return = supply + delta_t)."""
    def apply(m, Ts):
        set_water_inlet(m.fs.chilled_water_in, temperature=Ts + delta_t)
        m.fs.chilled_water_supply.fix(Ts + 273.15)

    set_operating_mode(m, "fixed_supply", chilled_water_supply=7.0)
    df = sweep(m, "Chilled water supply setpoint (C)", list(supplies), apply, start=7.0).sort_index()
    step_to(m, apply, supplies[-1], 7.0)
    restore_design(m)
    return df


PUMP_POWER_REF = 0.26  # kW at the design cooling water flow (15 m head, 70 % efficiency; an assumption)


def pump_power(flow, flow_ref=DESIGN["cooling_water_flow"], power_ref=PUMP_POWER_REF):
    """Cooling water pump power (kW) by the affinity laws, P ~ Q^3."""
    return power_ref * (flow / flow_ref) ** 3


def cooling_water_flow_study(m, flows, chilled_water_supply=7.0):
    """Compressor + cooling water pump power vs cooling water flow, at constant load (7 C supply)."""
    def apply(m, f):
        set_water_inlet(m.fs.cooling_water_in, flow=f)

    set_operating_mode(m, "fixed_supply", chilled_water_supply=chilled_water_supply)
    df = sweep(m, "Cooling water flow (kg/s)", list(flows), apply, start=DESIGN["cooling_water_flow"])
    step_to(m, apply, flows[-1], DESIGN["cooling_water_flow"])
    restore_design(m)
    df["Pump power (kW)"] = [pump_power(f) for f in df.index]
    df["Total power (kW)"] = df["Compressor power (kW)"] + df["Pump power (kW)"]
    return df


def optimize_cooling_water_flow(m, chilled_water_supply=7.0, bounds=(0.6, 3.0)):
    """Cooling water flow that minimizes compressor + pump power, by IDAES/Pyomo optimization."""
    set_operating_mode(m, "fixed_supply", chilled_water_supply=chilled_water_supply)
    solve(m)
    flow = m.fs.cooling_water_in.flow_mol[0]
    flow.unfix()
    flow.setlb(bounds[0] / MW_WATER)
    flow.setub(bounds[1] / MW_WATER)
    flow_kg = flow * MW_WATER / (pyunits.kg / pyunits.mol)
    m.fs.total_power = Objective(
        expr=m.fs.compressor.work_mechanical[0] / 1000
        + PUMP_POWER_REF * (flow_kg / DESIGN["cooling_water_flow"]) ** 3,
        sense=minimize,
    )
    try:
        solve(m)
        best_flow = value(flow_kg)
        best = {"Optimal cooling water flow (kg/s)": best_flow}
        best.update(chiller_results(m))
        best["Pump power (kW)"] = pump_power(best_flow)
        best["Total power (kW)"] = best["Compressor power (kW)"] + best["Pump power (kW)"]
    finally:
        m.fs.del_component(m.fs.total_power)
        flow.setlb(None)
        flow.setub(None)
        flow.fix()
    step_to(m, lambda m, f: set_water_inlet(m.fs.cooling_water_in, flow=f), best_flow, DESIGN["cooling_water_flow"])
    restore_design(m)
    return best


def part_load_map(m, loads=(1.0, 0.75, 0.5, 0.25), cooling_water_temps=(30, 26, 34), fouling=0.0,
                  chilled_water_supply=7.0, design_delta_t=5.0):
    """COP vs cooling load and cooling water temperature at a given condenser fouling.

    Load is set by the chilled water return temperature at constant flow and
    7 C supply. The compressor isentropic efficiency is kept constant (a
    simplification: real compressors lose efficiency at part load).
    """
    def set_load(m, load):
        set_water_inlet(m.fs.chilled_water_in, temperature=chilled_water_supply + load * design_delta_t)

    def set_tcw(m, T):
        set_water_inlet(m.fs.cooling_water_in, temperature=T)

    def set_rf(m, rf):
        set_fouling(m, condenser=rf)

    set_operating_mode(m, "fixed_supply", chilled_water_supply=chilled_water_supply)
    step_to(m, set_rf, 0.0, fouling)
    rows = []
    current_tcw, current_load = DESIGN["cooling_water_inlet"], 1.0
    for Tcw in cooling_water_temps:
        step_to(m, set_load, current_load, 1.0)
        step_to(m, set_tcw, current_tcw, Tcw)
        current_tcw, current_load = Tcw, 1.0
        for load in loads:
            step_to(m, set_load, current_load, load)
            current_load = load
            r = chiller_results(m)
            rows.append({"Cooling water inlet (C)": Tcw, "Load (%)": load * 100, "Condenser fouling (m2K/W)": fouling,
                         "Cooling duty (kW)": r["Cooling duty (kW)"],
                         "Compressor power (kW)": r["Compressor power (kW)"], "COP": r["COP"]})
    # Back to the design point
    step_to(m, set_load, current_load, 1.0)
    step_to(m, set_tcw, current_tcw, DESIGN["cooling_water_inlet"])
    step_to(m, set_rf, fouling, 0.0)
    restore_design(m)
    return pd.DataFrame(rows)
