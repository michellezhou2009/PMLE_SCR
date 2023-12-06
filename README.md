This repository provides the R code for the simulation studies and data analysis in the manuscript "Two-Stage Pseudo Maximum Likelihood Estimation of Semiparametric Copula-based Regression Models for Semi-Competing Risks Data" by Sakie J. Arachchige, Xinyuan Chen and Qian M. Zhou. 

The "R" folder includes all the R files:

- "_RunSimulation.R_" includes the R code for running simulation studies and produces individual ".RData" files, saved in "sim_results" folder, for each replication under each simulation setting. The R code requires R functions defined in three R files: "_helpers.R_", "_fitSPT.R_", and "_SemiCompCopSTPFuns.R_", and R libraries: _survival_, _dplyr_, _tidyr_, _purrr_, _VineCopula_, and _trust_. In addition, it requires three R libraries: _foreach_, _doSNOW_, _parallel_ for parallel computing. 

- A .RData file, named "_SimRes.RData_" in the "RData" folder, collects all the results from the replications where the maximization of the log-likelihood function and pseudo log-likelihood function converged. 

- "_ResToSummary.R_" includes the R code for calculating the summary statistics using the results in "_SimRes.RData_", and the summary statistics are stored in the "_Summary.RData_" in the "RData" folder.

- "_CreateTables.R_" includes the R code for creating Tables 1 - 5 of the manuscript using "_Summary.RData_".

- "_data_analysis.R_" includes the R code for the analysis of Bone Marrow Transplant data and creating Table 6 of the manuscript. The R code requires three R files: "helpers.R", "_fitSPT.R_", and "_SemiCompCopSTPFuns.R_", and R libraries: _survival_, _dplyr_, _tidyr_, _purrr_, _VineCopula_, and _trust_.

All the generated tables are stored in the .csv files in the "Tables" folder.
