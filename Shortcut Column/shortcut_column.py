# Import objects from pyomo package
from pyomo.environ import (ConcreteModel, Constraint, Var, value, log,
                           TransformationFactory, units as pyunits)
from pyomo.common.config import ConfigValue, ListOf

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock and the base classes for a custom unit model
from idaes.core import FlowsheetBlock, declare_process_block_class, UnitModelBlockData, useDefault
from idaes.core.util.config import is_physical_parameter_block

# Import the feed and product unit models
from idaes.models.unit_models import Feed, Product

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for the flowsheet
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the benzene-toluene property package (two phase, liquid only and vapor only)
from benzene_toluene import configuration, liquid_configuration, vapor_configuration

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

# Import propagate_state to pass the initialized states between unit models
from idaes.core.util.initialization import propagate_state

# Import Arc to connect the unit models
from pyomo.network import Arc

import re
import zipfile

import pandas as pd

VLE = ("Vap", "Liq")


@declare_process_block_class("ShortcutColumn")
class ShortcutColumnData(UnitModelBlockData):
    """Fenske-Underwood-Gilliland shortcut distillation column.

    * Fenske: minimum number of stages at total reflux
    * Underwood: minimum reflux ratio
    * Gilliland (Eduljee correlation): number of stages at the actual reflux ratio
    * Kirkbride: feed stage location

    The relative volatilities are evaluated from the Peng-Robinson K-values at
    the feed conditions, as in DWSIM's shortcut column. The distillate and
    bottoms leave as saturated liquids (total condenser).
    """

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault, domain=is_physical_parameter_block,
        doc="Two phase property package for the feed and products"))
    CONFIG.declare("liquid_property_package", ConfigValue(
        domain=is_physical_parameter_block, doc="Liquid only property package"))
    CONFIG.declare("vapor_property_package", ConfigValue(
        domain=is_physical_parameter_block, doc="Vapor only property package"))
    CONFIG.declare("light_key", ConfigValue(domain=str, doc="Light key component"))
    CONFIG.declare("heavy_key", ConfigValue(domain=str, doc="Heavy key component"))
    CONFIG.declare("light_non_keys", ConfigValue(
        default=[], domain=ListOf(str), doc="Components lighter than the light key (all to distillate)"))
    CONFIG.declare("heavy_non_keys", ConfigValue(
        default=[], domain=ListOf(str), doc="Components heavier than the heavy key (all to bottoms)"))

    def build(self):
        super().build()
        t = self.flowsheet().time.first()
        pp = self.config.property_package
        comps = pp.component_list
        lk, hk = self.config.light_key, self.config.heavy_key

        def state(name, package, defined):
            sb = package.build_state_block([t], defined_state=defined)
            self.add_component(name, sb)
            return sb[t]

        feed = state("properties_in", pp, True)
        dist = state("properties_distillate", pp, True)
        bott = state("properties_bottoms", pp, True)
        # Single phase states: feed liquid and vapor (K-values), feed bubble and
        # dew points (feed quality q) and the vapor leaving the top stage
        feed_liq = state("feed_liquid", self.config.liquid_property_package, True)
        feed_vap = state("feed_vapor", self.config.vapor_property_package, True)
        bubble_liq = state("feed_bubble_liquid", self.config.liquid_property_package, True)
        dew_vap = state("feed_dew_vapor", self.config.vapor_property_package, True)
        top_vap = state("top_vapor", self.config.vapor_property_package, True)

        self.add_port("inlet", self.properties_in)
        self.add_port("distillate", self.properties_distillate)
        self.add_port("bottoms", self.properties_bottoms)

        # Specifications
        self.reflux_ratio = Var(initialize=1.5, doc="Reflux ratio L/D")
        self.mole_frac_hk_distillate = Var(initialize=0.01, doc="Heavy key mole fraction in the distillate")
        self.mole_frac_lk_bottoms = Var(initialize=0.01, doc="Light key mole fraction in the bottoms")
        self.stage_height = Var(initialize=0.5, units=pyunits.m, doc="Height of one stage")

        # Results
        flow_units = pyunits.mol / pyunits.s
        self.alpha = Var(comps, initialize=2.0, doc="Relative volatility to the heavy key at the feed")
        self.q = Var(initialize=1.0, doc="Feed quality (liquid fraction added to the stripping section)")
        self.theta = Var(initialize=1.5, doc="Underwood root between the key volatilities")
        self.min_reflux_ratio = Var(initialize=1.0, doc="Underwood minimum reflux ratio")
        self.min_stages = Var(initialize=8.0, doc="Fenske minimum number of stages")
        self.stages = Var(initialize=12.0, doc="Number of theoretical stages (Gilliland, Eduljee)")
        self.rectifying_stages = Var(initialize=6.0, doc="Stages above the feed (Kirkbride)")
        self.stripping_stages = Var(initialize=6.0, doc="Stages below the feed (Kirkbride)")
        self.L = Var(initialize=10, units=flow_units, doc="Liquid flow in the rectifying section")
        self.V = Var(initialize=10, units=flow_units, doc="Vapor flow in the rectifying section")
        self.L_stripping = Var(initialize=10, units=flow_units, doc="Liquid flow in the stripping section")
        self.V_stripping = Var(initialize=10, units=flow_units, doc="Vapor flow in the stripping section")
        self.condenser_duty = Var(initialize=-4e5, units=pyunits.W, doc="Condenser heat duty (negative: removed)")
        self.reboiler_duty = Var(initialize=5e5, units=pyunits.W, doc="Reboiler heat duty")
        self.height = Var(initialize=7, units=pyunits.m, doc="Column height, (N + 2) stages")

        F, z = feed.flow_mol, feed.mole_frac_comp
        D, xD = dist.flow_mol, dist.mole_frac_comp
        B, xB = bott.flow_mol, bott.mole_frac_comp

        # Material balances and product specifications
        self.component_balance = Constraint(comps, rule=lambda b, j: F * z[j] == D * xD[j] + B * xB[j])
        self.distillate_sum = Constraint(expr=sum(xD[j] for j in comps) == 1)
        self.bottoms_sum = Constraint(expr=sum(xB[j] for j in comps) == 1)
        self.hk_spec = Constraint(expr=xD[hk] == self.mole_frac_hk_distillate)
        self.lk_spec = Constraint(expr=xB[lk] == self.mole_frac_lk_bottoms)
        self.light_non_key_split = Constraint(self.config.light_non_keys, rule=lambda b, j: xB[j] == 1e-10)
        self.heavy_non_key_split = Constraint(self.config.heavy_non_keys, rule=lambda b, j: xD[j] == 1e-10)

        # Products leave as saturated liquids
        self.distillate_temperature = Constraint(expr=dist.temperature == dist.temperature_bubble[VLE])
        self.bottoms_temperature = Constraint(expr=bott.temperature == bott.temperature_bubble[VLE])

        # Helper states follow the feed (or distillate) conditions
        def link(helper, source, T):
            yield helper.temperature == T
            yield helper.pressure == source.pressure
            yield helper.flow_mol == 1 * flow_units
            for j in comps:
                yield helper.mole_frac_comp[j] == source.mole_frac_comp[j]

        links = []
        links += list(link(feed_liq, feed, feed.temperature))
        links += list(link(feed_vap, feed, feed.temperature))
        links += list(link(bubble_liq, feed, feed.temperature_bubble[VLE]))
        links += list(link(dew_vap, feed, feed.temperature_dew[VLE]))
        links += list(link(top_vap, dist, dist.temperature_dew[VLE]))
        self.helper_links = Constraint(range(len(links)), rule=lambda b, i: links[i])

        # Relative volatility from the K-values at the feed, K = phi_liquid / phi_vapor
        def k_value(j):
            return feed_liq.fug_coeff_phase_comp["Liq", j] / feed_vap.fug_coeff_phase_comp["Vap", j]

        self.alpha_eq = Constraint(comps, rule=lambda b, j: b.alpha[j] * k_value(hk) == k_value(j))

        # Feed quality from the feed, saturated liquid and saturated vapor enthalpies
        self.q_eq = Constraint(expr=self.q * (dew_vap.enth_mol - bubble_liq.enth_mol)
                               == dew_vap.enth_mol - feed_liq.enth_mol)

        # Fenske
        self.fenske = Constraint(expr=self.min_stages * log(self.alpha[lk])
                                 == log(xD[lk] / xD[hk] * xB[hk] / xB[lk]))

        # Underwood
        self.underwood_theta = Constraint(expr=sum(
            self.alpha[j] * z[j] / (self.alpha[j] - self.theta) for j in comps) == 1 - self.q)
        self.underwood_reflux = Constraint(expr=self.min_reflux_ratio + 1 == sum(
            self.alpha[j] * xD[j] / (self.alpha[j] - self.theta) for j in comps))

        # Gilliland, Eduljee correlation: Y = 0.75 (1 - X^0.5668)
        X = (self.reflux_ratio - self.min_reflux_ratio) / (self.reflux_ratio + 1)
        Y = (self.stages - self.min_stages) / (self.stages + 1)
        self.gilliland = Constraint(expr=Y == 0.75 * (1 - X**0.5668))

        # Kirkbride feed stage location
        self.kirkbride = Constraint(expr=(self.rectifying_stages / self.stripping_stages) ** (1 / 0.206)
                                    == z[hk] / z[lk] * (xB[lk] / xD[hk]) ** 2 * B / D)
        self.stage_sum = Constraint(expr=self.rectifying_stages + self.stripping_stages == self.stages)

        # Internal flows
        self.L_eq = Constraint(expr=self.L == self.reflux_ratio * D)
        self.V_eq = Constraint(expr=self.V == self.L + D)
        self.L_stripping_eq = Constraint(expr=self.L_stripping == self.L + self.q * F)
        self.V_stripping_eq = Constraint(expr=self.V_stripping == self.V - (1 - self.q) * F)

        # Energy balances: total condenser, then overall balance for the reboiler.
        # The products are saturated liquids, so their liquid phase enthalpies are
        # used: exactly at the bubble point the phase split of the two phase
        # states is not unique.
        h_dist = dist.enth_mol_phase["Liq"]
        h_bott = bott.enth_mol_phase["Liq"]
        self.condenser_eq = Constraint(expr=self.condenser_duty == -self.V * (top_vap.enth_mol - h_dist))
        self.reboiler_eq = Constraint(expr=self.reboiler_duty
                                      == D * h_dist + B * h_bott - F * feed.enth_mol - self.condenser_duty)

        self.height_eq = Constraint(expr=self.height == (self.stages + 2) * self.stage_height)

    def initialize_build(self, outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        """Initialize the states from Raoult's law estimates, then solve the unit."""
        t = self.flowsheet().time.first()
        lk, hk = self.config.light_key, self.config.heavy_key
        feed = self.properties_in[t]
        self.properties_in.initialize(outlvl=outlvl)

        # Product estimates from the key specifications (binary key split)
        F = value(feed.flow_mol)
        z = {j: value(feed.mole_frac_comp[j]) for j in feed.mole_frac_comp}
        xhk, xlk = value(self.mole_frac_hk_distillate), value(self.mole_frac_lk_bottoms)
        D = F * (z[lk] - xlk) / (1 - xhk - xlk)
        for sb, flow, x in ((self.properties_distillate, D, {lk: 1 - xhk, hk: xhk}),
                            (self.properties_bottoms, F - D, {lk: xlk, hk: 1 - xlk})):
            s = sb[t]
            s.flow_mol.value = flow
            for j in z:
                s.mole_frac_comp[j].value = x.get(j, 1e-10)
            s.temperature.fix(value(feed.temperature_bubble[VLE]) if hasattr(feed, "temperature_bubble") else 360)
            sb.initialize(outlvl=outlvl)
            s.temperature.unfix()

        for sb, src, T in ((self.feed_liquid, feed, feed.temperature),
                           (self.feed_vapor, feed, feed.temperature),
                           (self.feed_bubble_liquid, feed, feed.temperature_bubble[VLE]),
                           (self.feed_dew_vapor, feed, feed.temperature_dew[VLE]),
                           (self.top_vapor, self.properties_distillate[t], self.properties_distillate[t].temperature_dew[VLE])):
            s = sb[t]
            s.temperature.value = value(T)
            s.pressure.value = value(src.pressure)
            s.flow_mol.value = 1
            for j in z:
                s.mole_frac_comp[j].value = value(src.mole_frac_comp[j])
            sb.initialize(outlvl=outlvl)

        # Seed the shortcut variables from Raoult's law relative volatility at the feed
        K = {j: value(self.feed_liquid[t].fug_coeff_phase_comp["Liq", j] / self.feed_vapor[t].fug_coeff_phase_comp["Vap", j])
             for j in z}
        for j in z:
            self.alpha[j].value = K[j] / K[hk]
        a = K[lk] / K[hk]
        self.theta.value = (1 + a) / 2
        self.min_stages.value = 9.19 / log(a)
        self.stages.value = 2 * value(self.min_stages)
        self.rectifying_stages.value = self.stripping_stages.value = value(self.stages) / 2

        solver_obj = get_solver(solver, optarg)
        result = solver_obj.solve(self)
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(result)))
        return result


def setup_shortcut_model(
    flow_mass=1.0,  # kg/s
    mole_frac={"benzene": 0.5, "toluene": 0.5},
    temperature=300,  # K
    pressure=101325,  # Pa
    reflux_ratio=1.5,
    mole_frac_lk_bottoms=0.01,
    mole_frac_hk_distillate=0.01,
    condenser_pressure=101325,  # Pa
    reboiler_pressure=101325,  # Pa
    stage_height=0.5,  # m
    outlvl=idaeslog.WARNING,
):
    """Creates and solves a benzene-toluene shortcut column."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.liquid_properties = GenericParameterBlock(**liquid_configuration)
    m.fs.vapor_properties = GenericParameterBlock(**vapor_configuration)

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.column = ShortcutColumn(
        property_package=m.fs.properties,
        liquid_property_package=m.fs.liquid_properties,
        vapor_property_package=m.fs.vapor_properties,
        light_key="benzene",
        heavy_key="toluene",
    )
    m.fs.distillate = Product(property_package=m.fs.properties)
    m.fs.bottoms = Product(property_package=m.fs.properties)

    """Connecting unit operations"""
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.column.inlet)
    m.fs.s02 = Arc(source=m.fs.column.distillate, destination=m.fs.distillate.inlet)
    m.fs.s03 = Arc(source=m.fs.column.bottoms, destination=m.fs.bottoms.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    print("The initial DOF is {0}".format(degrees_of_freedom(m)))

    mw = sum(x * value(m.fs.properties.get_component(j).mw) for j, x in mole_frac.items())  # kg/mol
    m.fs.feed.flow_mol.fix(flow_mass / mw)  # mol/s
    m.fs.feed.temperature.fix(temperature)  # K
    m.fs.feed.pressure.fix(pressure)  # Pa
    for j, x in mole_frac.items():
        m.fs.feed.mole_frac_comp[0, j].fix(x)

    c = m.fs.column
    c.reflux_ratio.fix(reflux_ratio)
    c.mole_frac_lk_bottoms.fix(mole_frac_lk_bottoms)
    c.mole_frac_hk_distillate.fix(mole_frac_hk_distillate)
    c.stage_height.fix(stage_height)
    c.properties_distillate[0].pressure.fix(condenser_pressure)
    c.properties_bottoms[0].pressure.fix(reboiler_pressure)

    print("The final DOF is {0}".format(degrees_of_freedom(m)))

    m.fs.feed.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s01)
    m.fs.column.initialize(outlvl=outlvl)
    propagate_state(arc=m.fs.s02)
    propagate_state(arc=m.fs.s03)
    m.fs.distillate.initialize(outlvl=outlvl)
    m.fs.bottoms.initialize(outlvl=outlvl)

    solver = get_solver()
    result = solver.solve(m, tee=False)

    return m, result


def report_shortcut_properties(model):
    """Report properties of the solved model."""
    feed_prop = model.fs.feed.report()
    distillate_prop = model.fs.distillate.report()
    bottoms_prop = model.fs.bottoms.report()

    return feed_prop, distillate_prop, bottoms_prop


def shortcut_results(model):
    """Key shortcut column results (stages, mol/s, K, kW, m)."""
    c = model.fs.column
    dist, bott = c.properties_distillate[0], c.properties_bottoms[0]
    return {
        "Relative volatility (feed)": value(c.alpha["benzene"]),
        "Feed quality q": value(c.q),
        "Minimum stages (Fenske)": value(c.min_stages),
        "Minimum reflux ratio (Underwood)": value(c.min_reflux_ratio),
        "Stages (Gilliland)": value(c.stages),
        "Distillate flow (mol/s)": value(dist.flow_mol),
        "Bottoms flow (mol/s)": value(bott.flow_mol),
        "Distillate temperature (K)": value(dist.temperature),
        "Bottoms temperature (K)": value(bott.temperature),
        "L rectifying (mol/s)": value(c.L),
        "V rectifying (mol/s)": value(c.V),
        "L stripping (mol/s)": value(c.L_stripping),
        "V stripping (mol/s)": value(c.V_stripping),
        "Condenser duty (kW)": -value(c.condenser_duty) / 1000,
        "Reboiler duty (kW)": value(c.reboiler_duty) / 1000,
        "Column height (m)": value(c.height),
    }


def read_dwsim_results(path):
    """Read the shortcut column results from a solved DWSIM .dwxmz file.

    A .dwxmz is a zip archive holding the flowsheet as XML. Streams are found
    by their tag.
    """
    with zipfile.ZipFile(path) as z:
        xml_name = next(n for n in z.namelist() if n.endswith(".xml"))
        xml = z.read(xml_name).decode("utf-8", "ignore")

    def first(tag, text):
        return float(re.search(r"<{0}>([^<]+)<".format(tag), text).group(1))

    column = xml[xml.find("UnitOperations.ShortcutColumn<"):]

    # Stream tags (names shown on the flowsheet) are stored with the graphic objects
    tags = dict(re.findall(r"<Name>(MAT-[^<]+)</Name>\s*<Tag>([^<]*)</Tag>", xml))
    streams = {}
    for stream in xml.split("<ComponentName>MAT-")[1:]:
        streams[tags["MAT-" + stream[:36]]] = stream[stream.find("<ComponentName>Mixture<"):]

    nmin, rmin, n = first("m_Nmin", column), first("m_Rmin", column), first("m_N", column)
    return {
        # Relative volatility implied by DWSIM's Fenske result
        "Relative volatility (feed)": (99 * 99) ** (1 / nmin),
        "Feed quality q": (first("L_", column) - first("L", column)) / first("molarflow", streams["feed"]),
        "Minimum stages (Fenske)": nmin,
        "Minimum reflux ratio (Underwood)": rmin,
        "Stages (Gilliland)": n,
        "Distillate flow (mol/s)": first("molarflow", streams["distillate"]),
        "Bottoms flow (mol/s)": first("molarflow", streams["residue"]),
        "Distillate temperature (K)": first("temperature", streams["distillate"]),
        "Bottoms temperature (K)": first("temperature", streams["residue"]),
        "L rectifying (mol/s)": first("L", column),
        "V rectifying (mol/s)": first("V", column),
        "L stripping (mol/s)": first("L_", column),
        "V stripping (mol/s)": first("V_", column),
        "Condenser duty (kW)": first("m_Qc", column),
        "Reboiler duty (kW)": first("m_Qb", column),
        "Column height (m)": (n + 2) * first("StageHeight", column),
    }


def compare_with_dwsim(model, dwsim):
    """Table of IDAES vs DWSIM results with the absolute and relative difference."""
    df = pd.DataFrame({"IDAES": shortcut_results(model), "DWSIM": dwsim})
    df["Difference"] = df["IDAES"] - df["DWSIM"]
    df["Error (%)"] = df["Difference"] / df["DWSIM"] * 100
    return df


if __name__ == "__main__":
    model, result = setup_shortcut_model()
    print(result.solver.termination_condition)
    report_shortcut_properties(model)
    print(compare_with_dwsim(model, read_dwsim_results("DWSim Shortcut Column Model.dwxmz")))
