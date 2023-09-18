library(parallel)


###########################################
#####      likelihood functions      ######
###########################################

## compute Q
compute_Q <- function(states, actions, rewards, beta, rho, alpha){
  N <- length(states)
  Q <- matrix(0, 2, 2)
  Q_chain <- matrix(0, N, 4)
  if (length(alpha)==1){
    Q[1,1] <- alpha
    Q[2,2] <- alpha
    Q_chain[1,1] <- alpha
    Q_chain[1,4] <- alpha
  } else if (length(alpha)==2){
    Q[1,1] <- alpha[1]
    Q[2,2] <- alpha[2]
    Q_chain[1,1] <- alpha[1]
    Q_chain[1,4] <- alpha[2]
  } else {
    cat("Error, length of alpha should be 1 or 2.")
  }
  
  for (i in 1:N){
    s <- states[i]
    a <- actions[i]
    R <- rewards[i]
    Q[a,s] <- Q[a,s] + beta * (rho * R - Q[a,s])
    if (i < N){
      Q_chain[i+1,] <- c(Q)
    }
  }
  return(Q_chain)
}


## compute g
compute_g1 <- function(states, Q11, Q22, w, c, I_mat){
  N <- length(states)
  g0 <- rep(0, N)
  for (t in 1:N){
    s <- states[t]
    if ( s==2 ){
      g0[t] <- w * Q22[t] - (1-w) * Q11[t]
    }
    if ( s==1 ){
      g0[t] <- (1-w) * Q22[t] - w * Q11[t]
    }
  }
  
  I <- predict(I_mat, g0)
  i0 <- predict(I_mat, 0)
  I0 <- sweep(I, 2, i0, "-")
  g <- as.numeric( I0 %*% c )
  g <- pmax(g, -6)
  g <- pmin(g, 6)
  return(g)
}




## negative log likelihood function for one subject: theta = (beta1, rho1)
neg_log_lik1 <- function(theta, alpha, w, states, actions, rewards, c, I_mat){
  beta1 <- theta[1]
  rho1 <- theta[2]
  beta <- exp(beta1) / (1 + exp(beta1))
  rho <- exp(rho1)
  Q_chain <- compute_Q(states, actions, rewards, beta, rho, alpha)
  Q11 <- Q_chain[,1]
  Q22 <- Q_chain[,4]
  g <- compute_g1(states, Q11, Q22, w, c, I_mat)
  out <- - sum( ( actions - 1 ) * g - log( 1 + exp(g) ) )
  return(out)
}






## marginal negative log likelihood function for one subject (using bivariate Gauss-Hermit quadrature)
neg_log_lik2 <- function(mu, L, alpha, w, nn, states, actions, rewards, c, I_mat){
  pts <- gauss_hermite_2d(nn, mu, L, prune = NULL)
  m <- nrow(pts$nodes)
  if (m > 1){
    f_vec <- -apply(pts$nodes, 1, neg_log_lik1, alpha = alpha, w = w,
                    states=states, actions=actions, rewards=rewards, 
                    c=c, I_mat=I_mat)
  }
  if (m == 1){
    f_vec <- -neg_log_lik1(pts$nodes, alpha = alpha, w = w, 
                           states=states, actions=actions, rewards=rewards, 
                           c=c, I_mat=I_mat)
  }
  f0 <- max(f_vec)
  h_vec <- exp(f_vec - f0)
  int_h <- sum( pts$weights * h_vec )
  out <- -log(int_h) - f0
  
  return(out)
}






## joint negative log likelihood function for all subjects
joint_neg_log_lik2_parallel <- function(pp, nn, I_mat, S_mat, A_mat, R_mat, 
                                        equal_start=FALSE, cores = 15){
  mu <- pp[1:2]
  l <- pp[3:5]
  L <- matrix(0, 2, 2)
  L[1,1] <- l[1]
  L[2,2] <- l[2]
  L[2,1] <- l[3]
  if (equal_start){
    alpha <- pp[6]
    w <- pp[7]
    c <- pp[-c(1:7)]
  } 
  if (!equal_start){
    alpha <- pp[6:7]
    w <- pp[8]
    c <- pp[-c(1:8)]
  } 
  n <- ncol(S_mat)
  out <- mclapply(1:n, FUN = function(i){
    neg_log_lik2(mu, L, alpha, w, nn, S_mat[,i], A_mat[,i], R_mat[,i], c, I_mat)
  }, mc.cores = cores)
  out <- sum(unlist(out))
  return(out)
}





