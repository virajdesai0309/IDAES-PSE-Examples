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
from idaes.models.unit_models.pressure_changer import Pump ,Turbine, ThermodynamicAssumption

# Import a feed and product stream
from idaes.models.unit_models import Feed, Product

# Import Flash unit model
from idaes.models.unit_models import Flash

# Creating an instance of the ConcreteModel class
m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

# Creating property packages for individual unit operations
m.fs.properties = iapws95.Iapws95ParameterBlock()

# Assigning property packages to unit operations
m.fs.one = Feed(property_package=m.fs.properties)
m.fs.Pump_1 = Pump(
    dynamic = False,
    property_package=m.fs.properties,
    compressor=True,
    thermodynamic_assumption=ThermodynamicAssumption.isentropic)
m.fs.Flash_1 = Flash(
    dynamic = False,
    property_package=m.fs.properties,
    has_heat_transfer=True,
    has_pressure_change=True
)
m.fs.two = Product(property_package=m.fs.properties)
m.fs.three = Product(property_package=m.fs.properties)

# Call the degrees_of_freedom function, get initial DOF
DOF_initial = degrees_of_freedom(m)
print("The initial DOF is {0}".format(DOF_initial))

# Import necessary Pyomo components
from pyomo.environ import TransformationFactory

# Import the Arc component for connecting unit models
from pyomo.network import Arc

m.fs.s01 = Arc(source=m.fs.one.outlet, destination = m.fs.Pump_1.inlet)
m.fs.s02 = Arc(source=m.fs.Pump_1.outlet, destination = m.fs.Flash_1.inlet)
m.fs.s03 = Arc(source=m.fs.Flash_1.vap_outlet, destination = m.fs.two.inlet)
m.fs.s04 = Arc(source=m.fs.Flash_1.liq_outlet, destination = m.fs.three.inlet)

TransformationFactory("network.expand_arcs").apply_to(m)

# Call the degrees_of_freedom function, get initial DOF
degrees_of_freedom(m)

# Fix the stream inlet conditions
m.fs.Pump_1.inlet.flow_mol[0].fix(100) # mol/s

# Use htpx method to obtain the molar enthalpy of inlet stream at the given temperature and pressure conditions 
m.fs.Pump_1.inlet.enth_mol[0].fix(value(htpx(T=308.15*units.K, P=101325*units.Pa))) # T in K, P in Pa
m.fs.Pump_1.inlet.pressure[0].fix(101325)

# Fix pump conditions
m.fs.Pump_1.deltaP.fix(101325*5)
m.fs.Pump_1.efficiency_isentropic.fix(0.85)

# Fix the flash conditions
m.fs.Flash_1.heat_duty[0].fix(0)  # No heat transfer
m.fs.Flash_1.deltaP[0].fix(0)  # Pressure drop across the flash unit

# Call the degrees_of_freedom function, get final DOF
DOF_final = degrees_of_freedom(m)
print("The final DOF is {0}".format(DOF_final))

"""
Diagnosing the model
"""
from idaes.core.util import DiagnosticsToolbox
dt = DiagnosticsToolbox(m)
print(dt.report_structural_issues())

"""
Commenting the initalization and solving steps
# Initialize the flowsheet, and set the output at INFO level
m.fs.Pump_1.initialize(outlvl=idaeslog.INFO)
m.fs.one.initialize(outlvl=idaeslog.INFO)
m.fs.two.initialize(outlvl=idaeslog.INFO)
m.fs.three.initialize(outlvl=idaeslog.INFO)

# Solve the simulation using ipopt
# Note: If the degrees of freedom = 0, we have a square problem
from pyomo.environ import SolverFactory
opt = SolverFactory('ipopt')
solve_status = opt.solve(m, tee=True)
"""