This repository contains all code and relevant data used in the manuscript "Nonlinear vegetation controls on urban water availability emerge from precipitation–evapotranspiration coupling".

## Overview

This repository provides:

1. The R code used to fit Conditional Vine Copula models linking precipitation (P), evapotranspiration (ET), urban vegetation (LAIur), and peri-urban vegetation (LAIperi).
2. The input datasets used for scenario simulations.
3. The source data underlying the main figures presented in the manuscript.

Urban water availability (WA) is defined as:

WA = P − ET

The analysis includes 1,029 cities worldwide and covers the period 1981–2023 at daily temporal resolution.

The Vine Copula framework was used to characterize the dependence structure among hydroclimatic variables and to simulate urban water availability under four vegetation scenarios:

* **Both Dynamic scenario**: observed urban and peri-urban LAI trajectories.
* **Both Fixed scenario**: urban and peri-urban LAI fixed at baseline conditions.
* **Urban Change scenario**: urban LAI fixed at +1 standard deviation relative to the baseline period, while peri-urban LAI remains at baseline.
* **Peri-urban Change scenario**: peri-urban LAI fixed at +1 standard deviation relative to the baseline period, while urban LAI remains at baseline.



---

## Input Data

The following files are required to run the simulation model:

- P19812023.xlsx       : Daily precipitation time series (1981–2023)      

- ET19812023.xlsx      : Daily evapotranspiration time series (1981–2023) 

- LAI19812023ur.xlsx   : Daily urban LAI time series (1981–2023)          

- LAI19812023peri.xlsx : Daily peri-urban LAI time series (1981–2023)     


Each column represents one city and each row represents one day.

---

## Code

`Urban_WA_Scenario_Simulation.R`

Main workflow for:

* Harmonic reconstruction of daily LAI dynamics;
* Empirical cumulative distribution function (CDF) transformation;
* Conditional Vine Copula fitting;
* Scenario-based simulation of precipitation and evapotranspiration;
* Calculation of urban water availability (WA).

### Required R Packages

VineCopula
CDVineCopulaConditional
openxlsx
readxl


---

## Source Data for Figures

The following files contain the numerical source data used to generate figures in the main manuscript:

- SourceData_Fig.2.xlsx : Figure 2 
- SourceData_Fig.3.xlsx : Figure 3 
- SourceData_Fig.4.xlsx : Figure 4 
- SourceData_Fig.5.xlsx : Figure 5 
- SourceData_Fig.6.xlsx : Figure 6 

Figure 1 presents the conceptual framework and therefore does not contain source data.

---
## Computing Environment

The analyses were conducted under the following environment:

- Operating System: Windows 11
- R version: 4.4.3


## Reproducibility

### Study Period

1981–2023

### Baseline Period

1981–1990

### Copula Model

Conditional Vine Copula

### Model Selection Criterion

Bayesian Information Criterion (BIC)

### Variable Order

1. Precipitation (P)
2. Evapotranspiration (ET)
3. Urban LAI (LAIur)
4. Peri-urban LAI (LAIperi)

---
## Running the Code

All required input datasets are included in this repository.

- P19812023.xlsx
- ET19812023.xlsx
- LAI19812023ur.xlsx
- LAI19812023peri.xlsx

Run:

source("Urban_WA_Scenario_Simulation.R")


The script will generate simulated daily precipitation (P), evapotranspiration (ET), and water availability (WA) for all cities under four vegetation scenarios:

- Both_Dynamic
- Urban_Change
- Periurban_Change
- Both_Fixed
---
## License

This repository is distributed under the MIT License.
