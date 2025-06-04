from pyomo.environ import ConcreteModel, SolverFactory, units, value, Var, Param, Constraint
from idaes.core import FlowsheetBlock
from idaes.models.unit_models.pressure_changer import (
    PressureChanger,
    IsentropicPressureChangerInitializer,
    ThermodynamicAssumption
)
from idaes.models.properties import iapws95
import numpy as np

# Pump curve data (converted to SI units where needed)
flow_data_m3h = np.array([0, 0.725, 1.451, 2.176, 2.902, 3.265, 3.627, 3.99, 4.353, 4.716])  # m³/h
head_data_m = np.array([271.162, 268.903, 266.643, 259.864, 248.565, 237.267, 225.969, 214.67, 198.852, 180.775])  # m
eff_data = np.array([0, 0.3, 0.488, 0.638, 0.713, 0.735, 0.75, 0.735, 0.713, 0.675])  # fraction

# Convert flow from m³/h to mol/s for water
water_density = 996.557  # kg/m³
water_mw = 0.01801528  # kg/mol
flow_data_mol_s = flow_data_m3h * (water_density / water_mw) / 3600

# Convert head from m to J/kg
g = 9.81  # m/s²
head_data_J_kg = head_data_m * g

def pump_perf_callback(blk):
    parent = blk.parent_block()
    
    # Store curve data
    blk.flow_points = Param(range(len(flow_data_mol_s)), 
                      initialize=dict(enumerate(flow_data_mol_s)), 
                      units=units.mol/units.s)
    blk.head_points = Param(range(len(head_data_J_kg)), 
                   initialize=dict(enumerate(head_data_J_kg)), 
                   units=units.J/units.kg)
    blk.eff_points = Param(range(len(eff_data)), 
                 initialize=dict(enumerate(eff_data)), 
                 units=units.dimensionless)
    
    # Variables
    blk.flow_mol = Var(initialize=44.56, bounds = (0, 100), units=units.mol/units.s)
    blk.head_interp = Var(initialize=2500,  bounds = (0, 3000), units=units.J/units.kg)
    blk.eff_interp = Var(initialize=0.7, bounds = (0, 0.75), units=units.dimensionless)  # Added bounds
    
    @blk.Constraint()
    def flow_link(b):
        return b.flow_mol == parent.control_volume.properties_in[0].flow_mol
    
    @blk.Constraint()
    def head_interp_eqn(b):
        # Find segment
        for i in range(len(flow_data_mol_s)-1):
            if value(b.flow_mol) >= value(b.flow_points[i]) and value(b.flow_mol) <= value(b.flow_points[i+1]):
                slope = (b.head_points[i+1] - b.head_points[i]) / (b.flow_points[i+1] - b.flow_points[i])
                return b.head_interp == b.head_points[i] + slope * (b.flow_mol - b.flow_points[i])
        # Clip to bounds if outside
        if value(b.flow_mol) < value(b.flow_points[0]):
            return b.head_interp == b.head_points[0]
        else:
            return b.head_interp == b.head_points[-1]
    
    @blk.Constraint()
    def eff_interp_eqn(b):
        # Find segment
        for i in range(len(flow_data_mol_s)-1):
            if value(b.flow_mol) >= value(b.flow_points[i]) and value(b.flow_mol) <= value(b.flow_points[i+1]):
                slope = (b.eff_points[i+1] - b.eff_points[i]) / (b.flow_points[i+1] - b.flow_points[i])
                return b.eff_interp == b.eff_points[i] + slope * (b.flow_mol - b.flow_points[i])
        # Clip to bounds if outside
        if value(b.flow_mol) < value(b.flow_points[0]):
            return b.eff_interp == 0
        else:
            return b.eff_interp == b.eff_points[-1]
    
    @blk.Constraint(parent.flowsheet().time)
    def pump_efficiency_eqn(b, t):
        return parent.efficiency_isentropic[t] == b.eff_interp
    
    @blk.Constraint(parent.flowsheet().time)
    def pump_head_eqn(b, t):
        return b.head_isentropic[t] == b.head_interp
    
# Create model and flowsheet
solver = SolverFactory('ipopt')
solver.options = {
    'tol': 1e-6,
    'bound_relax_factor': 0.1,
    'max_iter': 1000,
    'print_level': 5
}
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = iapws95.Iapws95ParameterBlock()

# Configure pump with performance curve
m.fs.unit = PressureChanger(
    property_package=m.fs.properties,
    compressor=True,
    thermodynamic_assumption=ThermodynamicAssumption.isentropic,
    support_isentropic_performance_curves=True,
    isentropic_performance_curves={"build_callback": pump_perf_callback},
)

# Set inlet conditions (using water properties)
m.fs.unit.inlet.flow_mol[0].fix(44.5612)  # mol/s (~ 2.9 m³/h)
Pin = 1e5  # Pa
Pout = 3e5  # Pa (target - will be adjusted by pump curve)
hin = iapws95.htpx(T=300*units.K, P=Pin*units.Pa)
m.fs.unit.inlet.enth_mol[0].fix(hin)
m.fs.unit.inlet.pressure[0].fix(Pin)

# Initialize and solve
initializer = IsentropicPressureChangerInitializer()
initializer.config.solver_options = {
    'tol': 1e-6,
    'bound_relax_factor': 0.1,
    'max_iter': 1000
}
try:
    initializer.initialize(m.fs.unit)
    results = solver.solve(m, tee=True)
    
    # Results
    print("\nResults:")
    print(f"Flow rate: {value(m.fs.unit.inlet.flow_mol[0]):.1f} mol/s")
    print(f"Inlet pressure: {value(m.fs.unit.inlet.pressure[0]/1e5):.1f} bar")
    print(f"Outlet pressure: {value(m.fs.unit.outlet.pressure[0]/1e5):.1f} bar")
    print(f"Head: {value(m.fs.unit.performance_curve.head_interp)/9.81:.1f} m")
    print(f"Efficiency: {value(m.fs.unit.efficiency_isentropic[0])*100:.1f}%")
except Exception as e:
    print(f"Initialization failed: {str(e)}")
    # Try solving without initializer
    results = solver.solve(m, tee=True)
    if results.solver.termination_condition == 'optimal':
        print("\nManaged to solve without initializer:")
        print(f"Outlet pressure: {value(m.fs.unit.outlet.pressure[0]/1e5):.1f} bar")


import matplotlib.pyplot as plt

# Create figure and primary axis
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot Head curve (primary Y-axis)
color = 'tab:blue'
ax1.set_xlabel('Flow Rate (m³/h)')
ax1.set_ylabel('Head (m)', color=color)
ax1.plot(flow_data_m3h, head_data_m, 'o-', color=color, label='Head')
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid(True)

# Create secondary Y-axis for Efficiency
ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('Efficiency (%)', color=color)
ax2.plot(flow_data_m3h, eff_data*100, 's-', color=color, label='Efficiency')
ax2.tick_params(axis='y', labelcolor=color)

# Add operating point if solved successfully
if 'results' in locals() and results.solver.termination_condition == 'optimal':
    op_flow = value(m.fs.unit.inlet.flow_mol[0]) * water_mw / water_density * 3600  # Convert back to m³/h
    op_head = value(m.fs.unit.performance_curve.head_interp)/9.81
    op_eff = value(m.fs.unit.efficiency_isentropic[0])*100
    
    # Plot operating points
    ax1.plot(op_flow, op_head, 'k*', markersize=12, label='Operating Head')
    ax2.plot(op_flow, op_eff, 'k*', markersize=12, label='Operating Eff.')

# Add title and legends
plt.title('Pump Performance Curves')
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.tight_layout()
#plt.show()
plt.savefig('/home/viraj/Documents/Github/IDAES PSE Examples/Pump_Curves/pump_curves.png', dpi=300, bbox_inches='tight')
print("Pump curves saved as 'pump_curves.png'")