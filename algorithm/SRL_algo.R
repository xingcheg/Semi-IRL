library(splines)
library(splines2)

###############################################################
####       the proposed semiparametric RL algorithm        ####
###############################################################

SRL_algorithm <- function(S_mat, A_mat, R_mat, cores=10, knots, degree=2, x, 
                          equal_start=FALSE, start_val, start_f=NULL, 
                          Miter=50, nn=c(3,3)){
  n <- ncol(S_mat)
  N <- nrow(S_mat)
  ######### Initial Estimate #########
  I_mat <- iSpline(x, knots = knots, degree = degree, intercept = TRUE)
  I0 <- predict(I_mat, 0)
  I0_mat <- sweep(I_mat, 2, I0, "-")
  if (is.null(start_f)){
    c0 <- coef(lm(x~0+I0_mat))
  } else {
    c0 <- start_f
  }
  
  ## starting values
  Sigma0 <- start_val$Sigma0
  L0 <- t(chol(Sigma0))
  l0 <- c(L0[1,1], L0[2,2], L0[2,1])
  alpha0 <- start_val$alpha0
  w0 <- start_val$w0
  mu0 <- c(start_val$mu_nu0, 1)
  
  
  ###### method 2 ######

  st1 <- Sys.time()
  
  ## estimate (f, alpha, w) given (mu, Sigma)
  opt1 <- optim(par = c(mu0, l0, alpha0, w0, c0), fn = joint_neg_log_lik2_parallel, 
                lower = c(-4.5, 0.999, 0.05, 0.05, -1, 
                          rep(0, length(alpha0)), 0.5, rep(0, ncol(I0_mat))), 
                upper = c(0, 1.001, 3, 3, 1, 
                          rep(4, length(alpha0)), 1, rep(8, ncol(I0_mat))), 
                method = "L-BFGS-B", control = list(maxit = Miter, trace=1), nn = nn, 
                I_mat = I_mat, S_mat = S_mat, A_mat = A_mat, 
                R_mat = R_mat, cores = cores, equal_start=equal_start)
  
  st2 <- Sys.time()
  st12 <- difftime(st2, st1, units = "mins")
  cat(".. optimization using:", st12, "minutes, done! \n")
  
  Mu_hat <- opt1$par[1:2]
  l1 <- opt1$par[3:5]
  L1 <- matrix(0, 2, 2)
  L1[1,1] <- l1[1]
  L1[2,2] <- l1[2]
  L1[2,1] <- l1[3]
  Sigma_hat <- L1 %*% t(L1)
  
  if (equal_start==TRUE){
    alpha_hat <- opt1$par[6]
    w_hat <- opt1$par[7] 
    c_hat <- opt1$par[-c(1:7)]
  } else {
    alpha_hat <- opt1$par[6:7]
    w_hat <- opt1$par[8] 
    c_hat <- opt1$par[-c(1:8)]
  }

  
  df <- 5+length(alpha_hat)+ncol(I0_mat)
  output <- list(Mu_hat = Mu_hat, Sigma_hat = Sigma_hat, C_mat = c_hat, 
                 alpha_vec = alpha_hat, w_vec = w_hat, knots=knots, degree=degree,
                 AIC=2*opt1$value+2*df, BIC=2*opt1$value+log(n)*df  )
  
  return(output)
}
