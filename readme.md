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

```
python3 -m venv idaes-env
source idaes-env/bin/activate  # On Windows: idaes-env\Scripts\activate
```

### 3. Install Dependencies
```
pip install --upgrade pip
pip install idaes-pse
idaes get-extensions
```

## Tutorials and Examples

> [!TIP]
> Click on the link to access the tutorial and example notebooks.

| **Tutorial**           | **Description**                                                                 |
|-------------------------|---------------------------------------------------------------------------------|
| [Compound Creation](Compound_Creation) | Learn how to create and configure compounds for process modeling.              |
| [Flash](Flash)             | Learn how to create and configure flash unit opeation.  |
| [Pump](Pump) | Learn how to create and configure pump unit opeation for process modeling.              |
| [Heater](Heater_Cooler) | Learn how to create and configure heater unit opeation and process modelling. |
| [Mixer](Mixer) | Learn how to create and configure mixer unit opeation and process modelling. |
| [Separator(Splitter)](Separator) | Learn how to create and configure separator aka splitter unit opeation and process modelling. |

> [!IMPORTANT]
> The following tutorials are in progress

| **Tutorial**           | **Description**                                                                 |
|-------------------------|---------------------------------------------------------------------------------|
| [Compressor]() | Learn how to create and configure compressor unit opeation and process modelling. |
| [Turbine]() | Learn how to create and configure turbine unit opeation and process modelling. |
| [Valve]() | Learn how to create and configure valve unit opeation and process modelling. |
| [Translator]() | Learn how to create and configure translator unit opeation and process modelling. |
| [Statejunction]() | Learn how to create and configure statejunction unit opeation and process modelling. |
| [HeatExchanger 0D]() | Learn how to create and configure heat exchanger 0D unit opeation and process modelling. |
| [HeatExchanger 1D]() | Learn how to create and configure heat exchanger 1D unit opeation and process modelling. |

## Documentation and links

- [IDAES PSE Documentation](https://idaes-pse.readthedocs.io/en/stable/)
- [IDAES PSE GitHub Repository](https://github.com/IDAES/idaes-pse)
