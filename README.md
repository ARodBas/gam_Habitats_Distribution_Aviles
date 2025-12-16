  
  # Species distribution and density modelling using GAMs

This repository contains the R scripts and data used to reproduce the analyses
and figures presented in the associated manuscript.

The objective of the study is to model species presence–absence and density
patterns using Generalized Additive Models (GAMs) and environmental predictors.

---

## Repository structure

├── scripts/
│ ├── 01_models_zone1.R
│ ├── 02_models_zone2.R
│
├── data/
│ ├── pres_abs/ Presence–absence data
│ ├── densities/ Species density data
│ └── predictors/ Environmental predictors
│
├── outputs/
│ ├── pres_abs/ Binary presence–absence maps
│ └── delta/ Delta maps
│
└── README.md

----


---

## How to reproduce the analysis

1. Clone or download this repository.
2. Open the R project file (`.Rproj`) in RStudio.
3. Run the scripts in the `scripts/` folder in numerical order:
   - `01_data_preparation.R`
   - `02_model_fitting.R`
   - `03_model_evaluation.R`

All outputs will be written to the `outputs/` directory.

---

## Data

All input data required to reproduce the analyses are stored in the `data/`
directory.

Subfolders contain:
- presence–absence observations,
- species density data,
- environmental predictor variables.

Data are provided in R-native formats (`.rds`) when possible to ensure
reproducibility.

---

## Requirements

The analysis was developed and tested using:

- R (>= 4.2)
- Packages:
  - mgcv
  - terra
  - sf
  - dismo
  - dplyr
  - DHARMa
  - mgcViz

---


