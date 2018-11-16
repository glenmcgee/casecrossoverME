# casecrossoverME
Code to run simulation study for "Corrections for Measurement Error Due to Delayed Onset of Illness for Case-Crossover Designs"

## Code Files
#### full_sim_singleboot0.R
Generates events based on true exposure time series. Defines functions for fitting marginal likelihood, regression calibration, and conditional score methods. along with function to generate events based on true exposure time series. Loops over RR simulated datasets, each with BB bootstrap resamples.

#### combineResults.R
Used only if running a job array via sbatch.

## Data Files
#### MLmatrix.txt
Hourly exposure time series at all possible event times and matched control times. Number of matches varies depending on month/day.

#### brentpmhourly.csv
Hourly exposure time series (pre-processed). Used to create MLmatrix.txt.

#### delays.csv
Distribution of true delay (lag) times in hours.
