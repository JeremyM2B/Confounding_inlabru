# Getting started

This project includes various scripts and datasets used for simulations and applications related to Cadmium (Cd) concentrations in France. Below is an overview of the available directories and files:

- **Simulation.R**: Contains the code to reproduce the simulation results presented in the manuscript.

- **Application.R**: Includes the code for reproducing the application results from the manuscript.

- **Predition_maps.R**: Contains the code for generating the maps results featured in the manuscript.

- **function.R**: Contains functions used in the Simulation.R script.

- **data folder contains**:

  - regbiofr: Contains the .shp shapefiles of France.
  
  - df_FRANCE.Rdata: Contains data on Cadmium (Cd) concentrations in mosses throughout France, along with their coordinates.
  
  - df_france_cd_air_2026: Contains data on Cd concentrations in the air, modelled using the EMEP physical model.

  - simu_data.Rdata: Contains simulated datasets, allowing for the avoidance of code compilation in the Simulation.R script.


# Usage
To reproduce the results, run the scripts in the following order:
  1. **Simulation.R**
  2. **Application.R**
  3. **Prediction_maps.R**

Before running the scripts, make sure to load the necessary data from the data/folder.

# Requirements

Please ensure that the necessary packages are installed to run the R scripts.
