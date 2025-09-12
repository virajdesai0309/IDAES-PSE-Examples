# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed unit model
from idaes.models.unit_models import Feed, Compressor, Product

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for the flowsheet
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the gas property package to create a configuration file for the GenericParameterBlock
from gas import configuration

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Arc to connect the unit models
from pyomo.network import Arc

def setup_compressor_model(
    flow_mol=100,  # mol/s
    mole_frac_methane=1.0,
    pressure=101325,  # Pa
    temperature=298,  # K
):
    """Creates and solve a compressor model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.compressor = Compressor(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    """Connectiing unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.compressor.inlet)
    m.fs.s02 = Arc(source=m.fs.compressor.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print(degrees_of_freedom(m))

    m.fs.feed.flow_mol.fix(flow_mol)
    m.fs.feed.mole_frac_comp[0, "methane"].fix(mole_frac_methane)
    m.fs.feed.pressure.fix(pressure)
    m.fs.feed.temperature.fix(temperature)

    m.fs.compressor.deltaP.fix(250000) # Pa
    m.fs.compressor.efficiency_isentropic.fix(0.85)

    print(degrees_of_freedom(m))
    
    solver = get_solver()
    result = solver.solve(m, tee=True)  # Set tee=True to see solver output
    
    return m, result

def report_compressor_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    compressor_prop = model.fs.compressor.report()
    product_prop = model.fs.product.report()

    return feed_prop, compressor_prop, product_prop