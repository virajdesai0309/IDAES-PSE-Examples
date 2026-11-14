# IDAES PSE Examples Repository

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Pyomo](https://img.shields.io/badge/Pyomo-6.x-brightgreen)](http://www.pyomo.org/)
[![Stars](https://img.shields.io/github/stars/virajdesai0309/IDAES-PSE-Examples?style=social)](https://github.com/virajdesai0309/IDAES-PSE-Examples/stargazers)
[![Forks](https://img.shields.io/github/forks/virajdesai0309/IDAES-PSE-Examples?style=social)](https://github.com/virajdesai0309/IDAES-PSE-Examples/network/members)
[![Repo Size](https://img.shields.io/github/repo-size/virajdesai0309/IDAES-PSE-Examples)](https://github.com/virajdesai0309/IDAES-PSE-Examples)

This repository contains example notebooks and scripts using the **IDAES PSE (Process Systems Engineering)** framework for process modeling, thermodynamics, unit operations, and optimization using Pyomo.

---

## Setting Up Your Environment

### 1. Clone This Repository
```bash
git clone https://github.com/virajdesai0309/idaes-pse-examples.git
cd idaes-pse-examples
```

### 2. Create a Virtual Environment

Tested with Python 3.12.

```
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

### 3. Install Dependencies
```
pip install --upgrade pip
pip install -r requirements.txt
idaes get-extensions
```

## Tutorials and Examples

> [!TIP]
> Click on the link to access the tutorial and example notebooks.

| **Tutorial**           | **Description**                                                                 |
|-------------------------|---------------------------------------------------------------------------------|
| [Compound Creation](Compound_Creation) | Learn how to create and configure compounds for process modeling. |
| [Flash](Flash) | Learn how to create and configure flash unit operation.  |
| [Heater](Heater_Cooler) | Learn how to create and configure heater unit operation and process modelling. |
| [Pump](Pump) | Learn how to create and configure pump unit operation for process modeling. |
| [Pump Curves](Pump_Curves) | Learn how to create and configure pump unit operation using pump curves for process modeling. |
| [Pump Curves (DWSIM validated)](Pump%20Curves%20DWSim) | Learn how to model a pump from its head, efficiency and NPSH curves with smooth piecewise linear interpolation, validated against DWSIM. |
| [HeatExchanger 0D](Heat_Exchanger_0D) | Learn how to create and configure heat exchanger 0D unit operation and process modelling. |
| [HeatExchanger 1D](Heat_Exchanger_1D) | Learn how to create and configure heat exchanger 1D unit operation with temperature profiles and process modelling. |
| [Compressor](Compressor) | Learn how to create and configure compressor unit operation and process modelling. |
| [Turbine](Turbine) | Learn how to create and configure steam turbine unit operation using steam tables and process modelling. |
| [Mixer](Mixer) | Learn how to create and configure mixer unit operation for mixing gas streams using Peng-Robinson and process modelling. |
| [Separator(Splitter)](Separator) | Learn how to create and configure separator aka splitter unit operation and process modelling. |
| [Valve](Valve) | Learn how to create and configure valve unit operation for flashing liquid nitrogen, including adding a new component to the IDAES Helmholtz equation of state. |
| [Conversion Reactor](Conversion%20Reactor) | Learn how to create and configure a conversion (stoichiometric) reactor for methane combustion, including reaction packages and limiting reactants. |
| [Equilibrium Reactor](Equilibrium%20Reactor) | Learn how to create and configure an equilibrium reactor for the water-gas shift reaction, including a temperature-dependent equilibrium constant and a Translator between property packages. |
| [CSTR](Continous%20Stirred%20Tank%20Reactor) | Learn how to create and configure a continuous stirred tank reactor with reaction kinetics for methanol carbonylation (approximate comparison: UNIQUAC in DWSIM vs ideal liquid in IDAES). |
| [PFR](Plug%20Flow%20Reactor) | Learn how to create and configure a packed bed plug flow reactor for steam methane reforming, including Xu-Froment catalytic kinetics, Ergun pressure drop and isothermal operation. |
| [Shortcut Column](Shortcut%20Column) | Learn how to build a custom IDAES unit model: a Fenske-Underwood-Gilliland shortcut distillation column for benzene-toluene. |
| [Distillation Column](Distillation%20Column) | Learn how to create and configure a rigorous tray distillation column for benzene-toluene with Peng-Robinson, validated stage by stage, and check a shortcut column design. |
| [Chiller Skid](Refrigeration%20Chilled%20Water%20Skid/chiller_flowsheet_visualizer.ipynb) | Learn how to model a closed-loop R134a vapor compression chiller with floating condensing and evaporating pressures, validated against DWSIM. |
| [Chiller Digital Twin Case Studies](Refrigeration%20Chilled%20Water%20Skid) | Condenser fouling vs COP and compressor power, energy cost and optimal cleaning interval, a fouling soft sensor, and operating optimization (setpoint reset, cooling water flow, part load). |
| [Combined Cycle Gas Turbine](Combined%20Cycle%20Gas%20Turbine%20With%20Heat%20Recovery/ccgt_flowsheet_visualizer.ipynb) | Learn how to model a ~27 MW combined cycle power plant (gas turbine, single pressure HRSG and steam cycle) with natural gas combustion, Peng-Robinson flue gas and IAPWS-95 steam, validated against a DWSIM model built through the DWSIM Automation API. |
| [CCGT Case Studies](Combined%20Cycle%20Gas%20Turbine%20With%20Heat%20Recovery) | Ambient temperature derating, part load with TIT control vs inlet guide vanes, economics (LCOE, fuel and carbon price, simple vs combined cycle, dispatch), and steam cycle design optimization (drum pressure, steam temperature, HRSG pinch). |

> [!IMPORTANT]
> The following tutorials are in progress

| **Tutorial**           | **Description**                                                                 |
|-------------------------|---------------------------------------------------------------------------------|
| [Statejunction]() | Learn how to create and configure statejunction unit operation and process modelling. |

## Documentation and links

- [IDAES PSE Documentation](https://idaes-pse.readthedocs.io/en/stable/)
- [IDAES PSE GitHub Repository](https://github.com/IDAES/idaes-pse)
