#############################################################
## Run this .R file in local PC to simulate RL data 
## and estimate the RL data using the proposed algorithm.
#############################################################
source("algorithm/simulate_RL_data_funcs.R")
source("algorithm/joint_likelihood_funcs.R")
source("algorithm/gauss_hermite_2d.R")
source("algorithm/SRL_algo.R")

##########################################################
#################### Simulate Data ####################
##########################################################
## hyperparameters
# trial length
N <- 100
# probability of getting rewards
#p1 <- 0.3
#p2 <- 0.75
# number of subjects
n <- 100
# group-level learning rate & sensitivity
mu <- c(-2.5, 1)
Sigma <- matrix(c(0.25, -0.1, -0.1, 0.16), nrow = 2)
# uncertainty
w <- 0.8
# alpha (initial Q value)
alpha <- 2

# suppose the true function f(.) is semiparameteric: f(x) = sgn(x) * log(2|x|+1)
ff <- function(x){
  return( sign(x) * log(2*abs(x)+1) )
}


## initialize knots for the I-splines
x <- seq(-6, 6, 0.01)
knots <- c(-2.4, -1.2, -0.6, -0, 0.6, 1.2, 2.4)


## simulate parameters (subject-level random effects)
Theta1 <- t( mu + t(chol(Sigma)) %*% matrix(rnorm(2*n), nrow=2) )
Beta1 <- Theta1[,1]
Rho1 <- Theta1[,2] 

## simulate S, A, R (i.e. state, action, rewards)
A_mat <- matrix(0, N, n)
S_mat <- matrix(0, N, n)
R_mat <- matrix(0, N, n)
for (i in 1:n){
  S_mat[, i] <- rbinom(n=N, size=1, p=0.5) + 1
  experiment <- RL_updates(S_mat[,i], Beta1[i], Rho1[i], alpha, w, ff)
  A_mat[, i] <- experiment$actions
  R_mat[, i] <- experiment$rewards
}




##########################################################
#################### Estimation ####################
##########################################################

## initial values
l0 <- c(0.3, 0.3, 0)
L0 <- matrix(0, 2, 2)
L0[1,1] <- l0[1]
L0[2,2] <- l0[2]
L0[2,1] <- l0[3]
Sigma_hat0 <- L0 %*% t(L0)
start_val <- list(Sigma0 = Sigma_hat0,
                  alpha0 = c(1.8),
                  w0 = 0.75, mu_nu0 = -3)

## optimization
## change cores (CPU cores) to fit your local PC/cluster
## expect less than an hour to complete if you use multiple CPU cores

#---------------------------------------------
## S_mat: N x n state matrix
## A_mat: N x n action matrix
## R_mat: N x n reward matrix
## cores: number of CPU cores
## knots: I-spline knots
## degree: I-spline degree
## x: I-spline evaluate points
## start_f: starting spline coefficients for the f, (if NULL start at f(x)=x)
## start_val: initial value of parameters (except f)
## Miter: Maximum number of iterations
## nn: 2D Gauss-Hermite quadrature grids (i.e. c(3,3) means there are 9 knots totally)
#---------------------------------------------
output <- SRL_algorithm(S_mat, A_mat, R_mat, cores=10, knots=knots, degree=1, 
                        x=x, start_f = NULL, start_val=start_val, 
                        Miter=150, nn=c(3,3))


saveRDS(output, file ="output_one_sample_rep.rds")
