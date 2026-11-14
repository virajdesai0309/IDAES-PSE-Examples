"""Build the Combined Cycle Gas Turbine (CCGT) DWSIM model through the DWSIM Automation API.

Gas turbine (Brayton cycle):
    Air -> Air Compressor -> Combustor (conversion reactor, with Fuel) -> Gas Turbine -> HRSG gas path
Heat recovery steam generator (single pressure, counter-current; gas path SH -> EVAP -> ECO):
    Feedwater -> Feed Pump -> Economizer -> Evaporator -> Steam Drum -> Superheater -> Steam Turbine -> Condenser
The water side is open (feedwater source and condensate sink at the same state), so no recycle is needed
around the steam cycle. Two recycles on the gas path (between superheater/evaporator and evaporator/economizer)
break the counter-current coupling of the HRSG sections.

Run with the project venv (needs pythonnet and DWSIM installed):
    python build_dwsim_ccgt.py          # build, solve, report and save the .dwxmz
    python build_dwsim_ccgt.py --tune   # also re-tune fuel flow (TIT 1350 C) and feedwater flow (10 K pinch)
"""

import os
import sys

DWSIM_PATH = os.environ.get("DWSIM_PATH", r"C:\Users\VirajDesai\AppData\Local\DWSIM") + "\\"
sys.path.append(DWSIM_PATH)

import clr  # noqa: E402  (pythonnet)

for dll in ("CapeOpen.dll", "DWSIM.Automation.dll", "DWSIM.Interfaces.dll", "DWSIM.GlobalSettings.dll",
            "DWSIM.SharedClasses.dll", "DWSIM.Thermodynamics.dll", "DWSIM.UnitOperations.dll",
            "DWSIM.Inspector.dll", "System.Buffers.dll"):
    clr.AddReference(DWSIM_PATH + dll)

from System.IO import Directory  # noqa: E402
from System.Collections.Generic import Dictionary  # noqa: E402
from System import String, Double  # noqa: E402
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType  # noqa: E402
from DWSIM.UnitOperations import UnitOperations, Reactors  # noqa: E402
from DWSIM.Automation import Automation3  # noqa: E402
from DWSIM.GlobalSettings import Settings  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "DWSim CCGT Combined Cycle.dwxmz")

# ---------------------------------------------------------------------------
# Design basis
# ---------------------------------------------------------------------------
DESIGN = {
    # Gas turbine
    "ambient_T": 288.15,  # K (ISO 15 C)
    "ambient_P": 101325.0,  # Pa
    "air_flow": 40.0,  # kg/s (sized for ~27 MW net combined cycle)
    "air": {"Oxygen": 0.21, "Nitrogen": 0.78, "Argon": 0.01},  # mol
    "pressure_ratio": 17.0,
    "compressor_eff": 86.0,  # % isentropic
    "fuel": {"Methane": 0.90, "Ethane": 0.06, "Propane": 0.04},  # mol
    "fuel_T": 298.15,  # K
    "fuel_flow": 0.9777,  # kg/s (tuned for a 1350 C turbine inlet temperature, see tune())
    "combustor_dP": 0.5e5,  # Pa
    "turbine_eff": 89.0,  # % isentropic
    "turbine_outlet_P": 1.06e5,  # Pa (HRSG back pressure)
    # HRSG and steam cycle
    "steam_flow": 8.074,  # kg/s feedwater (tuned for a 10 K HRSG pinch, see tune())
    "condenser_P": 0.08e5,  # Pa
    "pump_outlet_P": 62e5,  # Pa
    "pump_eff": 80.0,  # %
    "evaporator_P_drop": 1e5,  # Pa, water side economizer + evaporator
    "superheater_P_drop": 1e5,  # Pa
    "steam_T": 773.15,  # K (500 C)
    "approach": 10.0,  # K subcooling at the economizer outlet
    "evaporator_quality": 0.98,
    "steam_turbine_eff": 85.0,  # %
    "gas_dP_per_section": 1000.0,  # Pa, HRSG gas side per section
}


def build(design=DESIGN):
    Directory.SetCurrentDirectory(DWSIM_PATH)
    interf = Automation3()
    sim = interf.CreateFlowsheet()

    for c in ("Methane", "Ethane", "Propane", "Oxygen", "Nitrogen", "Argon", "Carbon dioxide", "Water"):
        sim.AddCompound(c)
    pr = sim.CreateAndAddPropertyPackage("Peng-Robinson (PR)")
    steam = sim.CreateAndAddPropertyPackage("Steam Tables (IAPWS-IF97)")

    # Combustion reactions, 100 % conversion of each fuel component (excess air)
    rset = sim.CreateReactionSet("Combustion", "Complete combustion of the fuel gas")
    sim.AddReactionSet(rset)
    reactions = {
        "Methane": {"Methane": -1, "Oxygen": -2, "Carbon dioxide": 1, "Water": 2},
        "Ethane": {"Ethane": -2, "Oxygen": -7, "Carbon dioxide": 4, "Water": 6},
        "Propane": {"Propane": -1, "Oxygen": -5, "Carbon dioxide": 3, "Water": 4},
    }
    for rank, (base, stoich) in enumerate(reactions.items()):
        coeffs = Dictionary[String, Double]()
        for k, v in stoich.items():
            coeffs[k] = float(v)
        rx = sim.CreateConversionReaction(base + " combustion", base + " combustion", coeffs, base, "Vapor", "100")
        sim.AddReaction(rx)
        sim.AddReactionToSet(rx.ID, rset.ID, True, rank)

    objs = {}

    def add(kind, tag, x, y, package=None):
        o = sim.AddObject(kind, x, y, tag)
        if package is not None:
            set_property_package(o.GetAsObject(), package)
        objs[tag] = o
        return o.GetAsObject()

    def connect(a, b, ia=-1, ib=-1):
        sim.ConnectObjects(objs[a].GraphicObject, objs[b].GraphicObject, ia, ib)

    # ---------------- Gas turbine ----------------
    air = add(ObjectType.MaterialStream, "Air", 0, 100, pr)
    add(ObjectType.Compressor, "Air Compressor", 100, 100, pr)
    add(ObjectType.MaterialStream, "Compressed Air", 200, 100, pr)
    add(ObjectType.EnergyStream, "Compressor Power", 100, 180)
    fuel = add(ObjectType.MaterialStream, "Fuel", 200, 20, pr)
    add(ObjectType.Mixer, "Combustor Inlet Mixer", 280, 100, pr)
    add(ObjectType.MaterialStream, "Combustor Feed", 340, 100, pr)
    add(ObjectType.RCT_Conversion, "Combustor", 400, 100, pr)
    add(ObjectType.MaterialStream, "Hot Gas", 480, 100, pr)
    add(ObjectType.MaterialStream, "Combustor Liquid", 480, 160, pr)
    add(ObjectType.EnergyStream, "Combustor Duty", 400, 180)
    add(ObjectType.Expander, "Gas Turbine", 560, 100, pr)
    add(ObjectType.MaterialStream, "GT Exhaust", 640, 100, pr)
    add(ObjectType.EnergyStream, "GT Power", 560, 180)

    connect("Air", "Air Compressor")
    connect("Air Compressor", "Compressed Air")
    connect("Compressor Power", "Air Compressor")
    connect("Compressed Air", "Combustor Inlet Mixer", -1, 0)
    connect("Fuel", "Combustor Inlet Mixer", -1, 1)
    connect("Combustor Inlet Mixer", "Combustor Feed")
    connect("Combustor Feed", "Combustor")
    connect("Combustor", "Hot Gas", 0, -1)
    connect("Combustor", "Combustor Liquid", 1, -1)
    connect("Combustor Duty", "Combustor")
    connect("Hot Gas", "Gas Turbine")
    connect("Gas Turbine", "GT Exhaust")
    connect("Gas Turbine", "GT Power")

    # Stream and unit specifications
    air.SetTemperature(design["ambient_T"])
    air.SetPressure(design["ambient_P"])
    air.SetMassFlow(design["air_flow"])
    set_composition(air, sim, design["air"])

    p_comb = design["ambient_P"] * design["pressure_ratio"]
    fuel.SetTemperature(design["fuel_T"])
    fuel.SetPressure(p_comb)
    fuel.SetMassFlow(design["fuel_flow"])
    set_composition(fuel, sim, design["fuel"])

    comp = objs["Air Compressor"].GetAsObject()
    comp.CalcMode = UnitOperations.Compressor.CalculationMode.OutletPressure
    comp.POut = p_comb
    comp.AdiabaticEfficiency = design["compressor_eff"]

    objs["Combustor Inlet Mixer"].GetAsObject().PressureCalculation = UnitOperations.Mixer.PressureBehavior.Minimum

    rx = objs["Combustor"].GetAsObject()
    rx.ReactionSetID = rset.ID
    rx.ReactorOperationMode = Reactors.OperationMode.Adiabatic
    rx.DeltaP = design["combustor_dP"]

    gt = objs["Gas Turbine"].GetAsObject()
    gt.CalcMode = UnitOperations.Expander.CalculationMode.OutletPressure
    gt.POut = design["turbine_outlet_P"]
    gt.AdiabaticEfficiency = design["turbine_eff"]

    # ---------------- HRSG and steam cycle ----------------
    # Gas path: GT Exhaust -> Superheater -> Evaporator -> Economizer -> Stack
    # Water path: Feedwater -> Feed Pump -> Economizer -> Evaporator -> Steam Drum -> Superheater -> Steam Turbine -> Condenser
    water = {"Water": 1.0}
    feedwater = add(ObjectType.MaterialStream, "Feedwater", 640, 420, steam)
    add(ObjectType.Pump, "Feed Pump", 720, 420, steam)
    add(ObjectType.EnergyStream, "Pump Power", 720, 500)
    add(ObjectType.MaterialStream, "Pumped Water", 800, 420, steam)
    add(ObjectType.HeatExchanger, "Economizer", 880, 300, steam)
    add(ObjectType.MaterialStream, "Economizer Water Out", 960, 420, steam)
    add(ObjectType.HeatExchanger, "Evaporator", 800, 200, steam)
    add(ObjectType.MaterialStream, "Evaporator Outlet", 880, 200, steam)
    add(ObjectType.Vessel, "Steam Drum", 960, 200, steam)
    add(ObjectType.MaterialStream, "Drum Steam", 1040, 160, steam)
    add(ObjectType.MaterialStream, "Blowdown", 1040, 240, steam)
    add(ObjectType.HeatExchanger, "Superheater", 720, 100, steam)
    add(ObjectType.MaterialStream, "Superheated Steam", 1200, 100, steam)
    add(ObjectType.Expander, "Steam Turbine", 1280, 100, steam)
    add(ObjectType.EnergyStream, "ST Power", 1280, 180)
    add(ObjectType.MaterialStream, "Turbine Exhaust Steam", 1360, 100, steam)
    add(ObjectType.Cooler, "Condenser", 1440, 100, steam)
    add(ObjectType.EnergyStream, "Condenser Duty", 1440, 180)
    add(ObjectType.MaterialStream, "Condensate", 1520, 100, steam)
    # Flue gas through the HRSG
    # The recycles (tear streams) sit on the gas path: the water side then solves straight through and the
    # counter-current coupling of the three sections is converged on the gas temperatures.
    add(ObjectType.MaterialStream, "Gas After Superheater", 760, 160, pr)
    add(ObjectType.OT_Recycle, "Recycle Gas 1", 780, 180)
    add(ObjectType.MaterialStream, "Evaporator Gas In", 790, 190, pr)
    add(ObjectType.MaterialStream, "Gas After Evaporator", 840, 260, pr)
    add(ObjectType.OT_Recycle, "Recycle Gas 2", 860, 280)
    add(ObjectType.MaterialStream, "Economizer Gas In", 870, 290, pr)
    add(ObjectType.MaterialStream, "Stack Gas", 920, 360, pr)

    # Heat exchanger ports: inlet 0 = hot, 1 = cold; outlet 0 = hot, 1 = cold
    connect("GT Exhaust", "Superheater", -1, 0)
    connect("Drum Steam", "Superheater", -1, 1)
    connect("Superheater", "Gas After Superheater", 0, -1)
    connect("Superheater", "Superheated Steam", 1, -1)
    connect("Gas After Superheater", "Recycle Gas 1")
    connect("Recycle Gas 1", "Evaporator Gas In")
    connect("Evaporator Gas In", "Evaporator", -1, 0)
    connect("Economizer Water Out", "Evaporator", -1, 1)
    connect("Evaporator", "Gas After Evaporator", 0, -1)
    connect("Evaporator", "Evaporator Outlet", 1, -1)
    connect("Gas After Evaporator", "Recycle Gas 2")
    connect("Recycle Gas 2", "Economizer Gas In")
    connect("Economizer Gas In", "Economizer", -1, 0)
    connect("Pumped Water", "Economizer", -1, 1)
    connect("Economizer", "Stack Gas", 0, -1)
    connect("Economizer", "Economizer Water Out", 1, -1)
    connect("Evaporator Outlet", "Steam Drum", -1, 0)
    connect("Steam Drum", "Drum Steam", 0, -1)
    connect("Steam Drum", "Blowdown", 1, -1)
    connect("Feedwater", "Feed Pump")
    connect("Pump Power", "Feed Pump")
    connect("Feed Pump", "Pumped Water")
    connect("Superheated Steam", "Steam Turbine")
    connect("Steam Turbine", "Turbine Exhaust Steam")
    connect("Steam Turbine", "ST Power")
    connect("Turbine Exhaust Steam", "Condenser")
    connect("Condenser", "Condensate")
    connect("Condenser", "Condenser Duty")

    p_pump = design["pump_outlet_P"]
    p_eco = p_pump - design["evaporator_P_drop"] / 2
    p_drum = p_pump - design["evaporator_P_drop"]
    m_steam = design["steam_flow"]

    feedwater.SetPressure(design["condenser_P"])
    feedwater.SetTemperature(saturation_T(design["condenser_P"]) - 0.5)  # condensate, just subcooled
    feedwater.SetMassFlow(m_steam)
    set_composition(feedwater, sim, water)

    pump = objs["Feed Pump"].GetAsObject()
    pump.CalcMode = UnitOperations.Pump.CalculationMode.OutletPressure
    pump.Pout = p_pump
    pump.Efficiency = design["pump_eff"]

    def hx(tag, mode, cold_dp, **specs):
        h = objs[tag].GetAsObject()
        h.FlowDir = UnitOperations.FlowDirection.CounterCurrent
        h.CalculationMode = mode
        h.HotSidePressureDrop = design["gas_dP_per_section"]
        h.ColdSidePressureDrop = cold_dp
        for k, v in specs.items():
            setattr(h, k, v)

    mode = UnitOperations.HeatExchangerCalcMode
    hx("Economizer", mode.CalcTempHotOut, design["evaporator_P_drop"] / 2,
       ColdSideOutletTemperature=saturation_T(p_eco) - design["approach"])
    hx("Evaporator", mode.OutletVaporFraction2, design["evaporator_P_drop"] / 2,
       OutletVaporFraction2=design["evaporator_quality"])
    hx("Superheater", mode.CalcTempHotOut, design["superheater_P_drop"],
       ColdSideOutletTemperature=design["steam_T"])

    drum = objs["Steam Drum"].GetAsObject()
    drum.PressureCalculation = UnitOperations.Vessel.PressureBehavior.Minimum

    st = objs["Steam Turbine"].GetAsObject()
    st.CalcMode = UnitOperations.Expander.CalculationMode.OutletPressure
    st.POut = design["condenser_P"]
    st.AdiabaticEfficiency = design["steam_turbine_eff"]

    cond = objs["Condenser"].GetAsObject()
    cond.CalcMode = UnitOperations.Cooler.CalculationMode.OutletVaporFraction
    cond.OutletVaporFraction = 0.0
    cond.DeltaP = 0.0

    # Initial estimates for the gas-side tear streams (approximate flue gas composition)
    flue_guess = {"Nitrogen": 0.75, "Oxygen": 0.13, "Argon": 0.01, "Carbon dioxide": 0.035, "Water": 0.075}
    m_gas = design["air_flow"] + design["fuel_flow"]
    for tag, T in (("Evaporator Gas In", 850.0), ("Economizer Gas In", 600.0)):
        g = objs[tag].GetAsObject()
        g.SetTemperature(T)
        g.SetPressure(design["turbine_outlet_P"])
        g.SetMassFlow(m_gas)
        set_composition(g, sim, flue_guess)

    return interf, sim, objs


def set_property_package(obj, package):
    """Assign a property package to a stream or unit operation.

    Through pythonnet, ``obj.PropertyPackage = pp`` hits the interface property, which DWSIM ignores, so every
    object would silently use the flowsheet's first package. Setting the concrete class property works.
    """
    from System.Reflection import BindingFlags
    flags = BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic
    prop = next(p for p in obj.GetType().GetProperties(flags)
                if p.Name == "PropertyPackage" and p.DeclaringType.Name != "BaseClass")
    prop.SetValue(obj, package)


def saturation_T(P):
    """Saturation temperature of water (K) at pressure P (Pa)."""
    from CoolProp.CoolProp import PropsSI
    return PropsSI("T", "P", P, "Q", 0, "Water")


def set_composition(stream, sim, fractions):
    """Set the overall molar composition of a stream from a dict of compound mole fractions."""
    names = list(sim.SelectedCompounds.Keys)
    from System import Array
    x = Array[Double]([float(fractions.get(n, 0.0)) for n in names])
    stream.SetOverallComposition(x)


def solve(interf, sim):
    Settings.SolverMode = 0
    errors = interf.CalculateFlowsheet4(sim)
    return [str(e) for e in errors] if errors else []


def stream(objs, tag):
    return objs[tag].GetAsObject()


# Lower heating values at 25 C, kJ/mol, and molar masses, g/mol
LHV_MOLAR = {"Methane": 802.3, "Ethane": 1428.6, "Propane": 2043.1}
MOLAR_MASS = {"Methane": 16.043, "Ethane": 30.069, "Propane": 44.096}


def fuel_lhv(fuel):
    """LHV of the fuel gas in kJ/kg from its molar composition."""
    return sum(x * LHV_MOLAR[c] for c, x in fuel.items()) / sum(x * MOLAR_MASS[c] for c, x in fuel.items()) * 1000


def performance(objs, design=DESIGN):
    """Key plant performance figures from the solved flowsheet (powers in kW)."""
    e = lambda tag: stream(objs, tag).EnergyFlow  # noqa: E731
    T = lambda tag: stream(objs, tag).GetTemperature() - 273.15  # noqa: E731
    fuel = stream(objs, "Fuel")
    q_fuel = fuel.GetMassFlow() * fuel_lhv(design["fuel"])
    gt_net = e("GT Power") - e("Compressor Power")
    st_net = e("ST Power") - e("Pump Power")
    net = gt_net + st_net
    drum_T = T("Drum Steam")
    co2 = stream(objs, "Stack Gas")
    co2_kg_s = co2.GetMassFlow() * co2.GetOverallMassComposition()[list(co2.FlowSheet.SelectedCompounds.Keys).index("Carbon dioxide")]         if hasattr(co2, "FlowSheet") else float("nan")
    return {
        "Fuel LHV input (kW)": q_fuel,
        "Turbine inlet temperature (C)": T("Hot Gas"),
        "GT exhaust temperature (C)": T("GT Exhaust"),
        "Gas turbine gross (kW)": e("GT Power"),
        "Air compressor (kW)": e("Compressor Power"),
        "Gas turbine net (kW)": gt_net,
        "Steam turbine (kW)": e("ST Power"),
        "Feed pump (kW)": e("Pump Power"),
        "Steam cycle net (kW)": st_net,
        "Plant net power (kW)": net,
        "GT simple cycle efficiency LHV (%)": 100 * gt_net / q_fuel,
        "Combined cycle efficiency LHV (%)": 100 * net / q_fuel,
        "Heat rate LHV (kJ/kWh)": 3600 * q_fuel / net,
        "Steam flow (kg/s)": stream(objs, "Superheated Steam").GetMassFlow(),
        "Drum saturation temperature (C)": drum_T,
        "HRSG pinch (K)": T("Gas After Evaporator") - drum_T,
        "Economizer approach (K)": drum_T - T("Economizer Water Out"),
        "Stack temperature (C)": T("Stack Gas"),
        "Condenser duty (kW)": e("Condenser Duty"),
        "CO2 emissions (kg/MWh)": co2_kg_s * 3600 / (net / 1000),
    }


def tune(interf, sim, objs, design=DESIGN, TIT=1350.0, pinch=10.0, tol=0.05, max_iter=12):
    """Secant search: fuel flow for the turbine inlet temperature, then steam flow for the HRSG pinch.

    Returns the tuned (fuel_flow, steam_flow) in kg/s; the flowsheet is left solved at that point.
    """
    def secant(set_x, get_y, x0, x1, target):
        set_x(x0); y0 = get_y()
        for _ in range(max_iter):
            set_x(x1); y1 = get_y()
            if abs(y1 - target) < tol:
                return x1
            x0, y0, x1 = x1, y1, x1 + (target - y1) * (x1 - x0) / (y1 - y0)
        raise RuntimeError("tuning did not converge")

    T = lambda tag: stream(objs, tag).GetTemperature() - 273.15  # noqa: E731
    fuel, feedwater = stream(objs, "Fuel"), stream(objs, "Feedwater")

    def run(errs_ok=False):
        errs = solve(interf, sim)
        if errs:
            raise RuntimeError(errs[0].split(" --->")[0])

    def set_fuel(m): fuel.SetMassFlow(m); run()  # noqa: E704
    def set_steam(m): feedwater.SetMassFlow(m); run()  # noqa: E704
    m_fuel = secant(set_fuel, lambda: T("Hot Gas"), design["fuel_flow"], design["fuel_flow"] * 1.02, TIT)
    m_steam = secant(set_steam, lambda: T("Gas After Evaporator") - T("Drum Steam"),
                     design["steam_flow"], design["steam_flow"] * 1.05, pinch)
    return m_fuel, m_steam


def energy_balance(objs, design=DESIGN):
    """Energy balance closure of the gas side and the water side (kW).

    DWSIM stream enthalpies exclude heats of formation, so the combustion heat release enters the gas side
    as fuel flow x LHV (the reaction heat at 25 C with water as vapour).
    """
    H = lambda tag: stream(objs, tag).GetMassFlow() * stream(objs, tag).GetMassEnthalpy()  # noqa: E731
    e = lambda tag: stream(objs, tag).EnergyFlow  # noqa: E731
    q_hrsg = sum(stream(objs, t).Q for t in ("Superheater", "Evaporator", "Economizer"))
    q_comb = stream(objs, "Fuel").GetMassFlow() * fuel_lhv(design["fuel"])
    gas_in = H("Air") + H("Fuel") + q_comb + e("Compressor Power")
    gas_out = H("Stack Gas") + e("GT Power") + q_hrsg
    water_in = H("Feedwater") + e("Pump Power") + q_hrsg
    water_out = H("Condensate") + H("Blowdown") + e("ST Power") + e("Condenser Duty")
    return {"Gas side in (kW)": gas_in, "Gas side out (kW)": gas_out,
            "Gas side closure (%)": 100 * (gas_in - gas_out) / gas_in,
            "Water side in (kW)": water_in, "Water side out (kW)": water_out,
            "Water side closure (%)": 100 * (water_in - water_out) / water_in}


def report(objs, design=DESIGN):
    print("\nStreams")
    for tag in ("Air", "Compressed Air", "Fuel", "Hot Gas", "GT Exhaust", "Gas After Superheater", "Gas After Evaporator",
                "Stack Gas", "Feedwater", "Pumped Water", "Economizer Water Out", "Evaporator Outlet", "Drum Steam",
                "Blowdown", "Superheated Steam", "Turbine Exhaust Steam", "Condensate"):
        s = stream(objs, tag)
        print("  %-22s T=%8.2f C  P=%8.3f bar  m=%8.3f kg/s" % (tag, s.GetTemperature() - 273.15, s.GetPressure() / 1e5, s.GetMassFlow()))
    print("\nHRSG duties")
    for tag in ("Superheater", "Evaporator", "Economizer"):
        print("  %-22s %10.1f kW" % (tag, stream(objs, tag).Q))
    print("\nPerformance")
    for k, v in performance(objs, design).items():
        print("  %-38s %10.2f" % (k, v))
    print("\nEnergy balance")
    for k, v in energy_balance(objs, design).items():
        print("  %-38s %10.2f" % (k, v))


if __name__ == "__main__":
    interf, sim, objs = build()
    errs = solve(interf, sim)
    print("errors:", errs or "none")
    if "--tune" in sys.argv:
        print("tuned fuel and feedwater flow (kg/s):", tune(interf, sim, objs))
    report(objs)
    interf.SaveFlowsheet(sim, OUTPUT, True)
    print("\nsaved", OUTPUT)
