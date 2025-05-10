# Import objects from pyomo package
from pyomo.environ import ConcreteModel, value, TransformationFactory, units as pyunits

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the feed unit model
from idaes.models.unit_models import Feed, heat_exchanger, Product

# Import IAPWS property package to create a properties block for steam in the flowsheet
from idaes.models.properties import iapws95

from idaes.models.properties.iapws95 import htpx

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


feed_stream={
    "flow_mol": 100,  # mol/s
    "mole_frac_methanol": 0.5,
    "mole_frac_ethanol": 0.5,
    "pressure": 101325,  # Pa
    "temperature": 298,  # K
},

water_stream={
    "flow_mol": 100,  # mol/s
    "pressure": 101325,  # Pa
    "temperature": 450,  # K
}

m = ConcreteModel()

m.fs = FlowsheetBlock(dynamic=False)

m.fs.properties = GenericParameterBlock(**configuration)

m.fs.properties_shell = iapws95.Iapws95ParameterBlock()

m.fs.properties_tube = GenericParameterBlock(**configuration)

m.fs.feed = Feed(property_package=m.fs.properties)
m.fs.product = Product(property_package=m.fs.properties)

m.fs.steam = Feed(property_package=m.fs.iapws95.Iapws95ParameterBlock())
m.fs.condensate = Product(property_package=m.fs.iapws95.Iapws95ParameterBlock())
"""
def setup_separator_model(
    feed_stream={
        "flow_mol": 100,  # mol/s
        "mole_frac_methanol": 0.5,
        "mole_frac_ethanol": 0.5,
        "pressure": 101325,  # Pa
        "temperature": 298,  # K
    },
    water_stream={
        "flow_mol": 100,  # mol/s
        "pressure": 101325,  # Pa
        "temperature": 450,  # K
        "enthalphy": htpx(T=water_stream["temperature"]*pyunits.K, P=water_stream["pressure"]*pyunits.Pa)
        }
    ):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.separator = Separator(
        property_package=m.fs.properties,
        split_basis=SplittingType.totalFlow,
        outlet_list = ["To_Vessel_1", "To_Vessel_2",],
        ideal_separation = False,
        has_phase_equilibrium = False,
        )
    m.fs.to_vessel_1 = Product(property_package=m.fs.properties)
    m.fs.to_vessel_2 = Product(property_package=m.fs.properties)
    
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.separator.inlet)
    m.fs.s02 = Arc(source=m.fs.separator.To_Vessel_1, destination=m.fs.to_vessel_1.inlet)
    m.fs.s03 = Arc(source=m.fs.separator.To_Vessel_2, destination=m.fs.to_vessel_2.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print(degrees_of_freedom(m))

    m.fs.feed.flow_mol.fix(flow_mol)
    m.fs.feed.mole_frac_comp[0, "methanol"].fix(mole_frac_methanol)
    m.fs.feed.mole_frac_comp[0, "ethanol"].fix(mole_frac_ethanol)
    m.fs.feed.pressure.fix(pressure)
    m.fs.feed.temperature.fix(temperature)

    m.fs.separator.split_fraction[0, "To_Vessel_1"].fix(0.2)
    
    print(degrees_of_freedom(m))
    
    solver = get_solver()
    result = solver.solve(m, tee=True)  # Set tee=True to see solver output
    
    return m, result

def report_separator_properties(model):
    feed_prop = model.fs.feed.report()
    seprator_prop = model.fs.separator.report()
    to_vessel_1 = model.fs.to_vessel_1.report()
    to_vessel_2 = model.fs.to_vessel_1.report()
    
    return feed_prop, seprator_prop, to_vessel_1, to_vessel_2
"""