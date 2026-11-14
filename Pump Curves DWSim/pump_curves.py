# Import objects from pyomo package
from pyomo.environ import ConcreteModel, Constraint, Expression, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed, product and pump unit models
from idaes.models.unit_models import Feed, Product, Pump

# Import the IAPWS-95 steam table property package
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx

# Import a smooth maximum function to build the piecewise linear curves
from idaes.core.util.math import smooth_max

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import propagate_state to pass the initialized states between unit models
from idaes.core.util.initialization import propagate_state

# Import Arc to connect the unit models
from pyomo.network import Arc

import re
import zipfile

import pandas as pd

MW_WATER = 0.01801528  # kg/mol
GRAVITY = 9.81  # m/s2, as used by DWSIM for the head

# Pump curves of the DWSIM model (impeller 200 mm, 1450 rpm), flow in m3/h
FLOW_POINTS = [0, 3, 6, 9, 12, 13.5, 15, 16.5, 18, 19.5]  # m3/h
CURVES = {
    "head": [44.526, 44.154, 43.783, 42.67, 40.815, 38.96, 37.105, 35.249, 32.652, 29.684],  # m
    "efficiency": [0, 30, 48.8, 63.8, 71.3, 73.5, 75, 73.5, 71.3, 67.5],  # %
    "npsh_required": [0.506, 0.506, 0.563, 0.675, 0.844, 0.957, 1.126, 1.351, 1.688, 2.026],  # m
}


def piecewise_linear(x, xs, ys, eps=1e-6):
    """Linear interpolation of the points (xs, ys) at x, written as a smooth expression.

    y(x) = y0 + s0 (x - x0) + sum_k (s_k - s_(k-1)) max(0, x - x_k)

    Each max is smoothed (smooth_max with a tiny eps), so the expression is
    valid for any x during the solve, and is linear beyond the end points.
    """
    slopes = [(ys[k + 1] - ys[k]) / (xs[k + 1] - xs[k]) for k in range(len(xs) - 1)]
    y = ys[0] + slopes[0] * (x - xs[0])
    for k in range(1, len(slopes)):
        y += (slopes[k] - slopes[k - 1]) * smooth_max(0, x - xs[k], eps)
    return y


def setup_pump_curve_model(
    flow_mass=4.1547394372007966,  # kg/s (15 m3/h)
    temperature=298.15199942972197,  # K
    pressure=288159.28901000967,  # Pa
    flow_points=FLOW_POINTS,
    curves=CURVES,
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a water pump whose head and efficiency follow its pump curves."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.pump = Pump(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.s02 = Arc(source=m.fs.pump.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    """Pump curves: head, efficiency and NPSH required as functions of the volumetric flow"""
    inlet = m.fs.pump.control_volume.properties_in[0]
    flow_m3h = pyunits.convert(inlet.flow_vol, to_units=pyunits.m**3 / pyunits.hr) / (pyunits.m**3 / pyunits.hr)
    m.fs.pump.head = Expression(
        expr=piecewise_linear(flow_m3h, flow_points, curves["head"]) * pyunits.m, doc="Pump head from the curve")
    m.fs.pump.npsh_required = Expression(
        expr=piecewise_linear(flow_m3h, flow_points, curves["npsh_required"]) * pyunits.m, doc="NPSH required from the curve")
    m.fs.pump.npsh_available = Expression(
        expr=(inlet.pressure - inlet.pressure_sat) / (inlet.dens_mass * GRAVITY * pyunits.m / pyunits.s**2),
        doc="NPSH available at the pump inlet")

    # The curves set the pressure rise and the efficiency
    m.fs.pump.head_curve = Constraint(
        expr=m.fs.pump.deltaP[0] == inlet.dens_mass * GRAVITY * pyunits.m / pyunits.s**2 * m.fs.pump.head)
    m.fs.pump.efficiency_curve = Constraint(
        expr=m.fs.pump.efficiency_pump[0] == piecewise_linear(flow_m3h, flow_points, curves["efficiency"]) / 100)

    m.fs.feed.flow_mol.fix(flow_mass / MW_WATER)  # mol/s
    m.fs.feed.pressure.fix(pressure)  # Pa
    m.fs.feed.enth_mol.fix(
        value(htpx(T=temperature * pyunits.K, P=pressure * pyunits.Pa)))  # J/mol

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    # Initialize the pump at a pressure rise and efficiency guess, then release them to the curves
    m.fs.pump.deltaP.fix(3.5e5)
    m.fs.pump.efficiency_pump.fix(0.75)
    m.fs.pump.head_curve.deactivate()
    m.fs.pump.efficiency_curve.deactivate()
    m.fs.pump.initialize(outlvl=outlvl)
    m.fs.pump.deltaP.unfix()
    m.fs.pump.efficiency_pump.unfix()
    m.fs.pump.head_curve.activate()
    m.fs.pump.efficiency_curve.activate()
    propagate_state(arc=m.fs.s02)
    m.fs.product.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_pump_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    pump_prop = model.fs.pump.report()
    product_prop = model.fs.product.report()

    return feed_prop, pump_prop, product_prop


def pump_results(model):
    """Key pump results (m3/h, m, %, bar, kW, K)."""
    p = model.fs.pump
    inlet, outlet = p.control_volume.properties_in[0], p.control_volume.properties_out[0]
    return {
        "Volumetric flow (m3/h)": value(pyunits.convert(inlet.flow_vol, to_units=pyunits.m**3 / pyunits.hr)),
        "Head (m)": value(p.head),
        "Efficiency (%)": value(p.efficiency_pump[0]) * 100,
        "Pressure increase (bar)": value(p.deltaP[0]) / 1e5,
        "Outlet pressure (bar)": value(outlet.pressure) / 1e5,
        "Power (kW)": value(p.work_mechanical[0]) / 1000,
        "Hydraulic power (kW)": value(p.work_fluid[0]) / 1000,
        "Temperature rise (K)": value(outlet.temperature - inlet.temperature),
        "NPSH available (m)": value(p.npsh_available),
        "NPSH required (m)": value(p.npsh_required),
    }


def read_dwsim_results(path):
    """Read the pump results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Streams are found
    by their tag.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        # Numeric values only: <CurveHead> is also the name of the curve data block
        return float(re.search(r"<{0}>([-+0-9.eE]+)<".format(tag), text).group(1))

    pump = xml[xml.find("UnitOperations.Pump<"):]

    # Stream tags (names shown on the flowsheet) are stored with the graphic objects
    tags = dict(re.findall(r"<Name>(MAT-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        tag = tags.get("MAT-" + stream[:36])
        if tag:
            streams[tag] = stream[stream.find("<ComponentName>Mixture<"):]

    inlet, outlet = streams["Pump Inlet"], streams["Pump Outlet"]
    return {
        "Volumetric flow (m3/h)": first("CurveFlow", pump) * 3600,
        "Head (m)": first("CurveHead", pump),
        "Efficiency (%)": first("CurveEff", pump),
        "Pressure increase (bar)": first("DeltaP", pump) / 1e5,
        "Outlet pressure (bar)": first("pressure", outlet) / 1e5,
        "Power (kW)": first("CurvePower", pump),
        # DWSIM adds only the hydraulic work to the fluid enthalpy
        "Hydraulic power (kW)": first("massflow", outlet) * (first("enthalpy", outlet) - first("enthalpy", inlet)),
        "Temperature rise (K)": first("temperature", outlet) - first("temperature", inlet),
        "NPSH available (m)": first("NPSH", pump),
        "NPSH required (m)": first("CurveNPSHr", pump),
    }


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": pump_results(model), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_pump_curve_model()
    print(result.solver.termination_condition)
    report_pump_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Pump Curves Model.dwxmz")))
