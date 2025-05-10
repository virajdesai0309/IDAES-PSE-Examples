from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Feed, Product, HeatExchanger
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_smooth_callback
import idaes.logger as idaeslog
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from methanol_ethanol import configuration
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.network import Arc
from idaes.core.util.tables import create_stream_table_dataframe

def setup_heatexchanger_model(
    feed_stream,
    water_stream,
    hx_area=1.0,
    hx_u=100.0
):
    """Creates and solves a heat exchanger model with fixed U and area."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    
    # Add property packages
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.steam_properties = iapws95.Iapws95ParameterBlock()
    
    # Build units
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.steam = Feed(property_package=m.fs.steam_properties)
    m.fs.condensate = Product(property_package=m.fs.steam_properties)
    m.fs.heat_exchanger = HeatExchanger(
        delta_temperature_callback=delta_temperature_lmtd_smooth_callback,
        hot_side_name="shell",
        cold_side_name="tube",
        shell={"property_package": m.fs.steam_properties},
        tube={"property_package": m.fs.properties},
    )
    
    # Connect units with arcs
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.heat_exchanger.cold_side_inlet)
    m.fs.s02 = Arc(source=m.fs.steam.outlet, destination=m.fs.heat_exchanger.hot_side_inlet)
    m.fs.s03 = Arc(source=m.fs.heat_exchanger.cold_side_outlet, destination=m.fs.product.inlet)
    m.fs.s04 = Arc(source=m.fs.heat_exchanger.hot_side_outlet, destination=m.fs.condensate.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    
    # Check initial DOF
    print(f"Initial Degrees of Freedom: {degrees_of_freedom(m)}")
    
    # Fix feed conditions
    m.fs.feed.flow_mol[0].fix(feed_stream["flow_mol"])
    m.fs.feed.mole_frac_comp[0, "methanol"].fix(feed_stream["mole_frac_methanol"])
    m.fs.feed.mole_frac_comp[0, "ethanol"].fix(feed_stream["mole_frac_ethanol"])
    m.fs.feed.pressure[0].fix(feed_stream["pressure"])
    m.fs.feed.temperature[0].fix(feed_stream["temperature"])
    
    # Fix steam conditions
    m.fs.steam.outlet.flow_mol[0].fix(water_stream["flow_mol"])
    m.fs.steam.outlet.pressure[0].fix(water_stream["pressure"])
    steam_enthalpy = htpx(
        T=water_stream["temperature"] * pyunits.K,
        P=water_stream["pressure"] * pyunits.Pa
    )
    m.fs.steam.outlet.enth_mol[0].fix(steam_enthalpy)
    
    # Fix HX parameters
    m.fs.heat_exchanger.area.fix(hx_area)
    m.fs.heat_exchanger.overall_heat_transfer_coefficient[0].fix(hx_u)
    
    # Check final DOF
    print(f"Final Degrees of Freedom: {degrees_of_freedom(m)}")
    
    # Solve model
    solver = get_solver()
    results = solver.solve(m, tee=False)
    
    return m, results

def report_heatexchanger_properties(model):
    """Generates a report with key metrics and stream tables."""
    report = {}
    
    # Performance metrics
    report["heat_duty"] = value(model.fs.heat_exchanger.heat_duty[0])
    report["delta_T"] = value(model.fs.heat_exchanger.delta_temperature[0])
    
    # Stream table
    report["stream_table"] = create_stream_table_dataframe({
        "Cold Inlet": model.fs.feed.outlet,
        "Cold Outlet": model.fs.product.inlet,
        "Hot Inlet": model.fs.steam.outlet,
        "Hot Outlet": model.fs.condensate.inlet
    })
    
    return report

# Example Usage
if __name__ == "__main__":
    feed_stream = {
        "flow_mol": 1.0,
        "mole_frac_methanol": 0.5,
        "mole_frac_ethanol": 0.5,
        "pressure": 101325,
        "temperature": 300
    }
    water_stream = {
        "flow_mol": 1.0,
        "pressure": 101325,
        "temperature": 380
    }
    
    # Create and solve model
    m, results = setup_heatexchanger_model(feed_stream, water_stream, hx_area=1.0, hx_u=100.0)
    
    # Generate report
    report = report_heatexchanger_properties(m)
    print("\nHeat Duty:", report["heat_duty"], "W")
    print("ΔT Driving Force:", report["delta_T"], "K")
    print("\nStream Table:")
    print(report["stream_table"])