# Import pyomo package
from pyomo.environ import ConcreteModel, Constraint, value, units

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the IAPWS property package to create a properties block for steam in the flowsheet
from idaes.models.properties import iapws95

from idaes.models.properties.iapws95 import htpx

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
from idaes.core.util.model_statistics import degrees_of_freedom

# Import the default IPOPT solver
from idaes.core.solvers import get_solver

# Import a compressor unit
from idaes.models.unit_models.pressure_changer import Turbine, ThermodynamicAssumption

# Import a feed and product stream
from idaes.models.unit_models import Feed, Product, Heater, Flash
from idaes.models.unit_models.pressure_changer import Pump

# Creating the model block
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# Creating property packages for individual unit operations
m.fs.properties = iapws95.Iapws95ParameterBlock()

# Assigning property packages
m.fs.one = Feed(property_package=m.fs.properties)
m.fs.two = Product(property_package=m.fs.properties)
m.fs.Pump_1 = Pump(
    dynamic = False,
    property_package=m.fs.properties
)
m.fs.Heater_1 = Heater(
    dynamic = False,
    property_package = m.fs.properties
)
# Import necessary Pyomo components
from pyomo.environ import TransformationFactory

# Import the Arc component for connecting unit models
from pyomo.network import Arc

# Connecting unit operations
m.fs.s01 = Arc(source=m.fs.one.outlet, destination = m.fs.Pump_1.inlet)
m.fs.s02 = Arc(source=m.fs.Pump_1.outlet, destination = m.fs.Heater_1.inlet)
m.fs.s03 = Arc(source=m.fs.Heater_1.outlet, destination = m.fs.two.inlet)
TransformationFactory("network.expand_arcs").apply_to(m)

# Call the degrees_of_freedom function, get initial DOF
print("The DOF after connecting the unit operations is", degrees_of_freedom(m))

# Fix the stream inlet conditions
m.fs.one.flow_mol[0].fix(100) # mol/s

# Use htpx method to obtain the molar enthalpy of inlet stream at the given temperature and pressure conditions 
m.fs.one.enth_mol[0].fix(value(htpx(T=298.15*units.K, P=101325*units.Pa))) # T in K, P in Pa
m.fs.one.pressure[0].fix(101325)

# Fix pump conditions
m.fs.Pump_1.deltaP.fix(101325*5)
m.fs.Pump_1.efficiency_pump.fix(0.85)

# Fix Heater conditions
m.fs.Heater_1.heat_duty.fix(4000000)

# Call the degrees_of_freedom function, get initial DOF
print("The DOF after connecting the unit operations and defining the variables is", degrees_of_freedom(m))

# Initialize the flowsheet
m.fs.one.initialize(outlvl=idaeslog.INFO)
m.fs.Pump_1.initialize(outlvl=idaeslog.INFO)
m.fs.Heater_1.initialize(outlvl=idaeslog.INFO)
m.fs.two.initialize(outlvl=idaeslog.INFO)

# Solve the simulation using ipopt
# NOTE: If the degrees of freedom = 0, we have a square problem
from pyomo.environ import SolverFactory
opt = SolverFactory('ipopt')
solve_status = opt.solve(m, tee=True)

print(m.fs.Pump_1.report())
print(m.fs.Heater_1.report())