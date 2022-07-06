HMM_Metropolis_Hastings <- function(x, iterations, Gamma0, tau0, alpha_gamma, alpha_tau, beta_tau){
  acc <- 0 # initialise acceptance count
  m <- dim(Gamma0)[1] # number of states in HMM
  
  Gammas <- array(NA, dim = c(m, m, iterations + 1)) # empty 3D array to store gamma simulations
  Gammas[,,1] <- Gamma0 # setting first array as initial gamma matrix
  
  taus <- array(NA, dim = c(1, m, iterations+1)) # empty 3D array to store tau simulations
  taus[,,1] <- tau0 # setting first array as initial tau vector
  
  for (i in 1:iterations) {
    Gamma_star <- matrix(0, nrow = m, ncol = m)
    tau_star <- matrix(0, nrow = 1, ncol = m)
    for (j in 1:m){       # simulate each row of gamma independently from a gamma distribution
      for (k in 1:m){     # with beta = 1. For each row, dividing each entry over the sum makes it equivalent to
        # drawing from a Dirichlet Distribution
        Gamma_star[j, k] <- rgamma(n = 1, shape =  (alpha_gamma + Gammas[j,k,i]), rate = 1)
      }
      
      Gamma_star[j,] <-  Gamma_star[j,]/ sum(Gamma_star[j,])
      
      tau_star[1, j] <-rgamma(1, shape = (beta_tau + taus[1,j,i]), rate = beta_tau) #draw Tau from Gamma Distribution
    }
    
    
    
    
    p <- HMM_Log_Likelihood(x, Gamma_star, cumsum(tau_star)) - # likelihood of proposal
      HMM_Log_Likelihood(x, Gammas[,,i], cumsum(taus[,,i])) # likelihood of current
    
    for (j in 1:m){
      
      for (k in 1:m){
        p <- p +
          (dgamma(Gammas[j,k,i], shape = alpha_gamma + Gamma_star[j, k], rate = 1, log = TRUE)) +
          (dgamma(Gamma_star[j, k], shape = alpha_gamma, rate = 1, log = TRUE)) -
          (dgamma(Gamma_star[j, k], shape =  alpha_gamma + Gammas[j,k,i], rate = 1, log = TRUE)) -
          (dgamma(Gammas[j,k,i], shape = alpha_gamma, rate = 1, log = TRUE))
        
      }
      
      p <- p + 
        dgamma(tau_star[1, j], shape = alpha_tau , rate = beta_tau, log = TRUE) + # prior
        dgamma(taus[1,j,i], shape = beta_tau + tau_star[1, j], rate = beta_tau, log = TRUE)  - # proposal
        dgamma(taus[1,j,i], shape = alpha_tau, rate = beta_tau, log = TRUE) - # prior
        dgamma(tau_star[1, j], shape = (beta_tau + taus[1,j,i]), rate = beta_tau, log = TRUE) # proposal
    }
    
    
    if (log(runif(1)) < min(0, p)) { # if we accept the proposal
      Gammas[,,i+1] <- Gamma_star
      taus[,,i+1] <- tau_star
      acc <- acc + 1
    } 
    else { # if we reject the proposal
      Gammas[,,i+1] <- Gammas[,,i]
      taus[,,i+1] <- taus[,,i]
    }
    print(i)
  }
  print(acc/iterations)
  
  thetas <- array(0, dim = c(iterations + 1, m + (m**2))) # converting it into this format 
  # in order to be able to return it and calculate the effective sample size
  
  for (i in 1:(iterations+1)){
    thetas[i,] <- as.vector(c(taus[,,i], as.vector(t(Gammas[,,i]))))
  }
  return(thetas)
}



