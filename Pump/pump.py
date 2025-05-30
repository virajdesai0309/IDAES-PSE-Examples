# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, SolverFactory, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed unit model
from idaes.models.unit_models import Feed, Product
from idaes.models.unit_models.pressure_changer import Pump

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

# Import the IAPWS Property Package for water properties
from idaes.models.properties import iapws95

# Import the PhaseType from the Helmholtz equation of state
from idaes.models.properties.helmholtz.helmholtz import PhaseType

def setup_pump_methanol_ethanol_model(
    flow_mol=100,  # mol/s
    mole_frac_methanol=0.6,
    mole_frac_ethanol=0.4,
    pressure=101325,  # Pa
    temperature=298,  # K
    pressure_increase=101325,  # Pa (equivalent to 1 atm)
    efficiency_pump_one=0.75,  # Efficiency of the pump
    outlet_pressure=202650,  # Pa (2 atm)
    efficiency_pump_two=0.75,  # Efficiency of the pump
):
    """Creates and solve a pump model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.pump_one = Pump(property_package=m.fs.properties)
    m.fs.pump_two = Pump(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump_one.inlet)
    m.fs.s02 = Arc(source=m.fs.pump_one.outlet, destination=m.fs.pump_two.inlet)
    m.fs.s03 = Arc(source=m.fs.pump_two.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Print the degrees of freedom before fixing any variables
    print("Degrees of freedom before fixing variables:", degrees_of_freedom(m))

    m.fs.feed.flow_mol.fix(flow_mol)
    m.fs.feed.mole_frac_comp[0, "methanol"].fix(mole_frac_methanol)
    m.fs.feed.mole_frac_comp[0, "ethanol"].fix(mole_frac_ethanol)
    m.fs.feed.pressure.fix(pressure)
    m.fs.feed.temperature.fix(temperature)

    m.fs.pump_one.deltaP.fix(pressure_increase) # Pa
    m.fs.pump_one.efficiency_pump.fix(efficiency_pump_one)
    m.fs.pump_two.outlet.pressure.fix(outlet_pressure)  # Pa
    m.fs.pump_two.efficiency_pump.fix(efficiency_pump_two)

    # Print the degrees of freedom after fixing some variables
    print("Degrees of freedom after fixing variables:", degrees_of_freedom(m))
    
    # Initialize the model
    m.fs.feed.initialize()
    m.fs.pump_one.initialize()
    m.fs.pump_two.initialize()
    m.fs.product.initialize()

    # Solve the model   
    solver = SolverFactory("ipopt")
    result = solver.solve(m, tee=True)  # Set tee=True to see solver output
    
    return m, result

def report_pump_methanol_ethanol_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    pump_one_prop = model.fs.pump_one.report()
    pump_two_prop = model.fs.pump_two.report()
    product_prop = model.fs.product.report()

    return feed_prop, pump_one_prop, pump_two_prop, product_prop
