# casecrossoverME
Code to run simulation study for "Corrections for Measurement Error Due to Delayed Onset of Illness for Case-Crossover Designs"

## Code Files
#### full_sim_singleboot0.R
Runs simulations for given exposure effects. Defines functions for fitting marginal likelihood, regression calibration, and conditional score methods. along with function to generate events based on true exposure time series. Loops over RR simulated datasets, each with BB bootstrap resamples. Will run for beta=0 (uncomment for other levels). Lower RR and BB for speed (fewer datasets and resamples).

#### MeasErr_tables.R
Takes results of a simulation (from full_sim_singleboot0.R) and produces summary tables and boxplots. Will run for beta=0 (uncomment for other exposure levels).

#### real_data_analysis_singleboot.R
Code for data analysis. NOTE: This will not run, since it requires data that are not publicly available.

## Data Files
#### MLmatrix.txt
Hourly exposure time series at all possible event times and matched control times. Number of matches varies depending on month/day.

#### brentpmhourly.csv
Hourly exposure time series (pre-processed). Used to create MLmatrix.txt.

#### delays.csv
Distribution of true delay (lag) times in hours.
