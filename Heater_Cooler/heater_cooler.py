# Import objects from pyomo package
from pyomo.environ import ConcreteModel, SolverFactory, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed unit model
from idaes.models.unit_models import Feed, Product, Heater

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

def setup_heater_cooler_model(
    flow_mol=100,  # mol/s
    mole_frac_methanol=0.6,
    mole_frac_ethanol=0.4,
    pressure=101325,  # Pa
    temperature=298.15,  # K
    heater_duty=200000,  # W
    cooler_outlet_temperature=300,  # K
):
    """Creates and solve a heater cooler model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.HT_1 = Heater(property_package=m.fs.properties)
    m.fs.CL_1 = Heater(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    
    """Connectiing unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.HT_1.inlet)
    m.fs.s02 = Arc(source=m.fs.HT_1.outlet, destination=m.fs.CL_1.inlet)
    m.fs.s03 = Arc(source=m.fs.CL_1.outlet, destination=m.fs.product.inlet)
    
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Print the degrees of freedom before fixing any variables
    print("Degrees of freedom before fixing variables:", degrees_of_freedom(m))

    m.fs.feed.flow_mol.fix(flow_mol)
    m.fs.feed.mole_frac_comp[0, "methanol"].fix(mole_frac_methanol)
    m.fs.feed.mole_frac_comp[0, "ethanol"].fix(mole_frac_ethanol)
    m.fs.feed.pressure.fix(pressure)
    m.fs.feed.temperature.fix(temperature)

    m.fs.HT_1.heat_duty.fix(heater_duty) # W
    m.fs.CL_1.outlet.temperature.fix(cooler_outlet_temperature) # K

    # Print the degrees of freedom after fixing some variables
    print("Degrees of freedom after fixing variables:", degrees_of_freedom(m))

    # Initialize the model
    m.fs.feed.initialize()
    m.fs.HT_1.initialize()
    m.fs.CL_1.initialize()
    m.fs.product.initialize()

    # Solve the model
    solver = SolverFactory("ipopt")
    result = solver.solve(m, tee=True)  # Set tee=True to see solver output
    
    return m, result

def report_heater_cooler_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    heater_prop = model.fs.HT_1.report()
    heater_prop2 = model.fs.CL_1.report()
    product_prop = model.fs.product.report()

    return feed_prop, heater_prop, heater_prop2, product_prop