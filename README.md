# Semiparametric-RL-for-Mental-Disorder

R code for the analysis from the paper "A Semiparametric Inverse Reinforcement Learning Approach to Characterize Decision Making for Mental Disorders". by Xingche Guo, Donglin Zeng, and Yuanjia Wang. 

## Overview

### algorithm
Folder **"algorithm"** includes .R funcs for simulating RL data and estimating the parameters using our proposed optimization algorithms. Before running the R code, make sure the R packages ("splines", "splines2", "statmod") are installed.

* The .R file **"algorithm/simulate_RL_data_funcs.R"** contains functions to simulate RL data.

* The .R file **"algorithm/gauss_hermite_2d.R"** is the function to produce bivariate Gauss-Hermite quadrature knots and weights.

* The .R files **"algorithm/joint_likelihood_funcs.R"** contains functions to compute joint log-likelihood.

* The .R file **"algorithm/SRL_algo.R"** is the R function that performs our proposed optimization algorithm to estimate RL parameters.

* The .R files **"algorithm/joint_likelihood_funcs_linear.R"** contains functions to compute joint log-likelihood of a Linear RL model.

* The .R file **"algorithm/LRL_algo.R"** is the R function that performs our proposed optimization algorithm to estimate Linear RL parameters.





