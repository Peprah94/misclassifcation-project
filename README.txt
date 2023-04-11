
README.txt



The following explains what the various scripts, data and supplementary information is used for or were derived. The Rscripts and resulta are divided into tow main parts: The simulation study and analysis of the gull dataset.


A. Simulation study:

i. simulationData.R - Rscript used to simulate the data

ii. nimbleSimulatedData.R - Rscript used to run the MCMC using NIMBLE.

iii. constant.R, fixedIntercov.R, intercept.R, main.R, onlyCov.R, variable.R - Rscripts used to run a specific case study in the main paper. Depends on the simulated data from simulationData.R and the function defined in nimbleSimulatedData.R.

iv. plotSimulations - Rscript to plot and summarise the results from the simulations.

v. allDataResults[1:16].csv - Saved results from MCMC results for the simulation study

vi. estimateBiasPars.csv - Bias in all the parameters used in the simulation study.


B. Case study : Gull dataset

i. gull_data_formatting.R : Format the data downloaded from GBIF. The results of the formatting is saved as gull_df.RData.

ii. fullMLFormatting.R : Format the data and ML output for use in MCMC.

iii. nimble.R : Rscript used to run the MCMC using NIMBLE.

iv. constantData.R, fixedIntercovData.R, interceptData.R, mainData.R, onlyCovData.R, variableData.R, ML_estimatesData.R - Rscripts used to run a specific case study in the main paper. Depends on the simulated data from simulationData.R and the function defined in nimbleSimulatedData.R.

v. plot_for_data.R : Rscript to plot and summarise the results from the application of model to gull dataset.

vi. predPerformanceGullDataset.csv : Table of predictive performance for the various case study scenarios




