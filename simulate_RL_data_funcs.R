###########################################
######        simulate RL data       ######
###########################################


## mechanism of getting rewards
reward_mechanism <- function(a, s, case=1){
  R <- 0
  if (case==1){
    if (a==s){
      if (a==1){
        R <- rbinom(1, 1, 0.3)
      }
      if (a==2){
        R <- rbinom(1, 1, 0.75)
      }
    }
  }
  if (case==2){
    if (a==s){
      if (a==1){
        R <- rbeta(1, 1, 3)
      }
      if (a==2){
        R <- rbeta(1, 3, 1)
      }
    }
  }
  return(R)
}



## mechanism of updating Q: (action x state matrix)
update_Q <- function(a, s, R, Q, beta, rho){
  Q1 <- Q
  Q1[a,s] <- Q[a,s] + beta * (rho * R - Q[a,s])
  return(Q1)
}


## mechanism of taking actions
g_func <- function(Q, s, w, ff){
  g0 <- w * (Q[2,s] - Q[1,s]) + (1-w) *(Q[2,-s] - Q[1,-s])
  g <- ff(g0)
  return(g)
}

action_mechanism <- function(Q, s, w, ff){
  g <- g_func(Q, s, w, ff)
  prob <- 1 / (1 + exp(-g))
  a <- rbinom(1, 1, prob) + 1
  return(list(a = a, p = prob))
}




## function to simulate RL data
RL_updates <- function(states, beta1, rho1, alpha, w, ff, case=1){
  
  beta <- exp(beta1) / (1 + exp(beta1))
  rho <- exp(rho1)
  
  N <- length(states)
  actions <- rep(0, N)
  rewards <- rep(0, N)
  probs <- rep(0, N)
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
    
    ## take action
    a_out <- action_mechanism(Q, s, w, ff)
    a <- a_out$a
    p <- a_out$p
    actions[i] <- a
    probs[i] <- ifelse(s==2, p, 1-p)
    
    ## compute reward
    R <- reward_mechanism(a, s, case=case)
    rewards[i] <- R
    
    ## update Q
    Q <- update_Q(a, s, R, Q, beta, rho)
    if (i < N){
      Q_chain[i+1,] <- c(Q)
    }
  }
  return(list(actions=actions, rewards=rewards, Q_chain=Q_chain, probs=probs))
}





