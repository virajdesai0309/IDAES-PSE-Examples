# Import necessary Pyomo components
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver utility from IDAES
from idaes.core.solvers import get_solver

# Import the FlowsheetBlock for creating process flowsheets
from idaes.core import FlowsheetBlock

# Import unit models for feed, product, and heat exchanger
from idaes.models.unit_models import Feed, Product, HeatExchanger

# Import property packages for steam and generic properties
from idaes.models.properties import iapws95
from idaes.models.properties.iapws95 import htpx

# Import the LMTD temperature difference callback for the heat exchanger
from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_smooth_callback

# Import the IDAES logger for logging purposes
import idaes.logger as idaeslog

# Import the generic property package base class
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the configuration for methanol-ethanol property package
from methanol_ethanol import configuration

# Import utility for checking degrees of freedom in the model
from idaes.core.util.model_statistics import degrees_of_freedom

# Import the Arc component for connecting unit models
from pyomo.network import Arc

# Import utility for creating stream tables
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