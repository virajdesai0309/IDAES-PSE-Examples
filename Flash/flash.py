# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed unit model
from idaes.models.unit_models import Feed, Flash, Product

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

def setup_flash_model(
    flow_mol=100,  # mol/s
    mole_frac_methanol=0.6,
    mole_frac_ethanol=0.4,
    pressure=101325,  # Pa
    temperature=298,  # K
):
    """Creates and solve a flash model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.flash = Flash(property_package=m.fs.properties)
    m.fs.vapor = Product(property_package=m.fs.properties)
    m.fs.liquid = Product(property_package=m.fs.properties)

    """Connectiing unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.flash.inlet)
    m.fs.s02 = Arc(source=m.fs.flash.vap_outlet, destination=m.fs.vapor.inlet)
    m.fs.s03 = Arc(source=m.fs.flash.liq_outlet, destination=m.fs.liquid.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print(degrees_of_freedom(m))

    m.fs.feed.flow_mol.fix(flow_mol)
    m.fs.feed.mole_frac_comp[0, "methanol"].fix(mole_frac_methanol)
    m.fs.feed.mole_frac_comp[0, "ethanol"].fix(mole_frac_ethanol)
    m.fs.feed.pressure.fix(pressure)
    m.fs.feed.temperature.fix(temperature)

    m.fs.flash.heat_duty.fix(2500000) # W
    m.fs.flash.deltaP.fix(0)

    print(degrees_of_freedom(m))
    
    solver = get_solver()
    result = solver.solve(m, tee=True)  # Set tee=True to see solver output
    
    return m, result

def report_flash_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    flash_prop = model.fs.flash.report()
    vapor_prop = model.fs.vapor.report()
    liquid_prop = model.fs.liquid.report()

    return feed_prop, flash_prop, vapor_prop, liquid_prop