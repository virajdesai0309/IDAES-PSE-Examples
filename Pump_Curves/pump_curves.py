from pyomo.environ import ConcreteModel, SolverFactory, units, value
from idaes.core import FlowsheetBlock
from idaes.models.unit_models.pressure_changer import (
    PressureChanger,
    IsentropicPressureChangerInitializer,
    ThermodynamicAssumption
)
from idaes.models.properties import iapws95

# Define a callback function for performance curves
def pump_perf_callback(blk):
    # `blk` is the performance_curve block
    parent = blk.parent_block()  # Main pressure changer block

    # Add efficiency constraint (use `efficiency_pump` for pumps)
    @blk.Constraint(parent.flowsheet().time)
    def pump_efficiency_eqn(b, t):
        return parent.efficiency_isentropic[t] == 0.75  # Fixed efficiency

    # Add head constraint (use `head` instead of `head_isentropic`)
    @blk.Constraint(parent.flowsheet().time)
    def pump_head_eqn(b, t):
        return b.head_isentropic[t] == 100000 * units.J/units.kg  # Example head

# Create model and flowsheet
solver = SolverFactory('ipopt')
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = iapws95.Iapws95ParameterBlock()

# Configure a custom pump using PressureChanger
m.fs.unit = PressureChanger(
    property_package=m.fs.properties,
    compressor=True,  # Pump is a "compressor" (pressure increases)
    thermodynamic_assumption=ThermodynamicAssumption.isentropic,  # Enable curves
    support_isentropic_performance_curves=True,
    isentropic_performance_curves={"build_callback": pump_perf_callback},
)

# Set inlet conditions
m.fs.unit.inlet.flow_mol[0].fix(1000)  # mol/s
Pin = 1e5  # Inlet pressure (Pa)
Pout = 2e5  # Outlet pressure (Pa)
hin = iapws95.htpx(T=300*units.K, P=Pin*units.Pa)
m.fs.unit.inlet.enth_mol[0].fix(hin)
m.fs.unit.inlet.pressure[0].fix(Pin)

# Fix pressure ratio (outlet/inlet)
# m.fs.unit.ratioP[0].fix(2.0)

# Initialize and solve
initializer = IsentropicPressureChangerInitializer()
initializer.initialize(m.fs.unit)
results = solver.solve(m, tee=True)

# Check results
print("Outlet Pressure:", value(m.fs.unit.outlet.pressure[0]))
print("Pressure Ratio:", value(m.fs.unit.ratioP[0]))

