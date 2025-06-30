# Getting started

This project includes various scripts and datasets used for simulations and applications related to Cadmium (Cd) concentrations in France. Below is an overview of the available directories and files:

- **Simulation.R**: Contains the code to reproduce the simulation results presented in the manuscript.

- **Application.R**: Includes the code for reproducing the application results from the manuscript.

- **Predition_maps.R**: Contains the code for generating the maps results featured in the manuscript.

- **function.R**: Contains functions used in the Simulation.R script.

- **data** folder contains:

  - **regbiofr.zip**: This file can be unzipped to extract the .shp shapefiles containing geographical data of France.
  
  - **df_FRANCE.Rdata**: Contains data on Cadmium (Cd) concentrations in mosses throughout France, along with their coordinates.
  
  - **df_france_cd_air_2026**: Contains data on Cd concentrations in the air, modelled using the EMEP physical model.

  - **simu_data.Rdata**: Contains simulated datasets, allowing for the avoidance of code compilation in the Simulation.R script.


# Usage
To reproduce the results, run the scripts in the following order:
  1. **Simulation.R**
  2. **Application.R**
  3. **Prediction_maps.R**

Before running the scripts, make sure to load the necessary data from the data/folder.

# Requirements

Please ensure that the necessary packages are installed to run the R scripts.

# DOI
DOI: 10.5281/zenodo.15748177 

# Data availability statement
Data on Cd concentrations in mosses are publicly available on request from Sébastien Leblond (sebastien.leblond@mnhn.fr) or Caroline Meyer (caroline.meyer@mnhn.fr). Data on Cd concentrations in the air are publicly available on the website https://www.emep.int/ and on request from msc-east@ijs.si (Ilyin and Travnikov, 2005). All models used for this work are publicly available. The localisation of the data can be accessed via the following dataset on the Global Biodiversity Information Facility (GBIF):  https://doi.org/10.15468/c2bqic. Information regarding the Cd measurements in mosses is available at: https://bramm.mnhn.fr/cadmium-cd/ (Inventaire National du Patrimoine Naturel, 2017). The repositories were archived on Zenodo at the time of this publication (https://doi.org/10.5281/zenodo.15748177)

# References

Ilyin, I. and O. Travnikov (2005). Regional Model MSCE-HM of Heavy Metal Transboundary Air Pollution in Europe Technical report 8, 2005

Inventaire National du Patrimoine Naturel (2017). Données issues du dispositif de Biosurveillance des Retombées Atmosphériques Métalliques par les Mousses (BRAMM). UMS PatriNat (OFB-CNRS-MNHN), Paris. Occurrence dataset. https://doi.org/10.15468/c2bqic accessed via GBIF.org on 2025-06-30.
