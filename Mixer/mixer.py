# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed unit model
from idaes.models.unit_models import Feed, Product, Mixer, MomentumMixingType

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for the flowsheet
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the methanol_ethanol property package to create a configuration file for the GenericParameterBlock
from methanol_ethanol import configuration

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Arc to connect the unit models
from pyomo.network import Arc

def setup_mixer_model(
    methanol_stream={
        "flow_mol": 100,  # mol/s
        "mole_frac_methanol": 1.0,
        "mole_frac_ethanol": 1e-15,
        "pressure": 101325,  # Pa
        "temperature": 298,  # K
    },
    ethanol_stream={
        "flow_mol": 100,  # mol/s
        "mole_frac_methanol": 1e-15,
        "mole_frac_ethanol": 1.0,
        "pressure": 101325,  # Pa
        "temperature": 298,  # K
    }
):
    """Creates and solve a pump model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.methanol_feed = Feed(property_package=m.fs.properties)
    m.fs.ethanol_feed = Feed(property_package=m.fs.properties)
    m.fs.mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list = ["Methanol_Inlet", "Ethanol_Inlet"],
        momentum_mixing_type=MomentumMixingType.minimize)
    m.fs.product = Product(property_package=m.fs.properties)
    
    """Connectiing unit operations"""
    m.fs.s01 = Arc(source=m.fs.methanol_feed.outlet, destination=m.fs.mixer.Methanol_Inlet)
    m.fs.s02 = Arc(source=m.fs.ethanol_feed.outlet, destination=m.fs.mixer.Ethanol_Inlet)
    m.fs.s03 = Arc(source=m.fs.mixer.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print(degrees_of_freedom(m))

    m.fs.methanol_feed.flow_mol.fix(methanol_stream["flow_mol"])
    m.fs.methanol_feed.mole_frac_comp[0, "methanol"].fix(methanol_stream["mole_frac_methanol"])
    m.fs.methanol_feed.mole_frac_comp[0, "ethanol"].fix(methanol_stream["mole_frac_ethanol"])
    m.fs.methanol_feed.pressure.fix(methanol_stream["pressure"])
    m.fs.methanol_feed.temperature.fix(methanol_stream["temperature"])

    m.fs.ethanol_feed.flow_mol.fix(ethanol_stream["flow_mol"])
    m.fs.ethanol_feed.mole_frac_comp[0, "methanol"].fix(ethanol_stream["mole_frac_methanol"])
    m.fs.ethanol_feed.mole_frac_comp[0, "ethanol"].fix(ethanol_stream["mole_frac_ethanol"])
    m.fs.ethanol_feed.pressure.fix(ethanol_stream["pressure"])
    m.fs.ethanol_feed.temperature.fix(ethanol_stream["temperature"])

    print(degrees_of_freedom(m))
    
    solver = get_solver()
    result = solver.solve(m, tee=True)  # Set tee=True to see solver output
    
    return m, result

def report_mixer_properties(model):
    """Report properties of the solved model."""
    methanol_feed_prop = model.fs.methanol_feed.report()
    ethanol_feed_prop = model.fs.ethanol_feed.report() 
    mixer_prop = model.fs.mixer.report()
    product_prop = model.fs.product.report()
    
    return methanol_feed_prop, ethanol_feed_prop, mixer_prop, product_prop