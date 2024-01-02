# Simulation results

The code for reproducing the simulation results are stored in the subdirectories `model_i, i=1,2,3`.

Each of those subdirectories contains the following:

- `run_simulation.R` contains the code for running the simulation, see instructions below.
- the `simulation_results` subdirectory contains the raw output from the simulations
- `generate_summary.Rmd` contains to code for producing the simulations summaries presented in the paper, such as confidence intervals and error plots. Knitting the file will produce a `generate_summary_files` directory containing the plots presented in the paper as PDF. It will also generate a `generate_summary.pdf` and `generate_summary.tex` file, containing the table with the values for the confidence intervals.

## Re-running the simulations

The simulation data in the subdirectories `model_i/simulation_results` can be reproduced by running the `run_simulation.R` files. The file will look for the environment variable `SGE_TASK_ID` and uses its value for setting the seed (line 75). Running the file once will do 20 simulation runs each for each sample size. To produce the full 1000 simulation runs, the file is meant to be run 50 times, each time with a different seed.

# Simulation results for ER-Model with Covariates

The code for reproducing the simulation results for the ERC-model are stored in the subdirectories `erc/model_i, i=1,2,3`.

Each of those subdirectories contains the following:

- `ER-C.R` contains the code for running the simulation, see instructions below. As with the other simulation results, the script will lokk for the `SGE_TASK_ID` environment variable for setting the seed (line 34). Running the script once will run 20 simulations for each sample size and to reproduce the full simulation results the script is meant to be run 50 times, each time with a different seed.
- `ER-C_generate_summary.Rmd` contains to code for producing the simulations summaries presented in the paper. Knitting the file will produce a `ER-C_generate_summary.pdf` and `ER-C_generate_summary.tex` file, containing the table with the values for the confidence intervals.

# Lawyer network

The results from analyzing Lazega's lawyer data can be reproduced by running `lawyer_friendship_network.R`. Before running the code, please download the necessary data from the [original data owner's website](https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm) and place it in the `LazegaLawyers` subdirectory.

# World trade netwrok

The results from the analysis of the world trade network can be reproduced by running `trade_network.R`. Before running the code, please download the necessary data from the [original data owner's website](http://personal.lse.ac.uk/tenreyro/LGW.html) and place it in the `trade_network` subdirectory.