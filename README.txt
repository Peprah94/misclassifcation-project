R scripts
a. Data Simulation
The code should be run in this order for the results presented in the paper:
- simulation.R is used to simulate the data that is needed for the model.
- nimbleSimulation.R is a function used to fit the model and estimate the parameters in NIMBLE
- constant.R, covariate.R, fixed_intercov.R, fixed_covariate.R, main.R, intercept.R fit the study scenario models in the main article
- plotSimulations.R is used to summarize the results from the model fitting in tables and graphs

b. Application to Gull dataset from GBIF
The code should be run in this order for the results presented in the paper:
- gullDataFormat.R is used to format the data that is needed for the model. It needs occurrence data from gullGBIFdata folder.
- nimble.R is used to fit the model with NIMBLE.
- plotForData.R is used to summarise the results from nimble.R and present results in Tables and graphs.
