library(mcmcse)
library(MASS)
library(DirichletReg)

HMM_Metropolis_Hastings <- function(x, iterations, gamma0, tau0, alpha_gamma, alpha_tau, beta_tau){
  acc <- 0 # initialise acceptance count
  m <- dim(gamma0)[1] # number of states in HMM
  
  gammas <- array(NA, dim = c(m, m, iterations + 1)) # empty 3D array to store gamma simulations
  gammas[,,1] <- gamma0 # setting first array as initial gamma matrix
  
  taus <- array(NA, dim = c(1, m, iterations+1)) # empty 3D array to store tau simulations
  taus[,,1] <- tau0 # setting first array as initial tau vector
  
  for (i in 1:iterations) {
    gamma_star <- matrix(0, nrow = m, ncol = m)
    tau_star <- matrix(0, nrow = 1, ncol = m)
    for (j in 1:m){       # simulate each row of gamma independently from a gamma distribution
      for (k in 1:m){     # with beta = 1. For each row, dividing each entry over the sum makes it equivalent to
                          # drawing from a Dirichlet Distribution
        gamma_star[j, k] <- rgamma(n = 1, shape = alpha_gamma * gammas[j,k,i], rate = 1)
      }

      gamma_star[j,] <-  gamma_star[j,]/ sum( gamma_star[j,])

      # gamma_star[j, ] <- rdirichlet(n = 1, alpha = alpha_gamma * gammas[j,,i])
      tau_star[1, j] <-rgamma(1, shape = beta_tau * taus[1,j,i], rate = beta_tau) #draw Tau from Gamma Distribution
    }
    
    

    
    p <- log(HMM_Likelihood(x, gamma_star, cumsum(tau_star))) + # likelihood of proposal
        - log(HMM_Likelihood(x, gammas[,,i], cumsum(taus[,,i]))) # likelihood of current
    print(p)

    for (j in 1:m){
      for (k in 1:m){
        p <- p
        + (dgamma(gammas[j,k,i], shape = alpha_gamma * gamma_star[j, k], rate = 1, log = TRUE))
        + (dgamma(gamma_star[j, k], shape = alpha_gamma, rate = 1, log = TRUE))
        - (dgamma(gamma_star[j, k], shape = alpha_gamma * gammas[j,k,i], rate = 1, log = TRUE))
        - (dgamma(gammas[j,k,i], shape = alpha_gamma, rate = 1, log = TRUE))
        # print(p)
      }

      p <- p + 
              # ddirichlet(as.matrix(gammas[j,,i]), alpha = alpha_gamma * as.matrix(gamma_star[j, ]), sum.up = TRUE, log = T)
              # - ddirichlet(as.matrix(gamma_star[j, ]), alpha = alpha_gamma * as.matrix(gammas[j,,i]),  sum.up = TRUE, log = T) +
              dgamma(tau_star[1, j], shape = alpha_tau , rate = beta_tau, log = TRUE) + # prior
              (dgamma(taus[1,j,i], shape = beta_tau *  tau_star[1, j], rate = beta_tau, log = TRUE))  - # posterior
              (dgamma(taus[1,j,i], shape = alpha_tau, rate = beta_tau, log = TRUE)) - # prior
              (dgamma(tau_star[1, j], shape = beta_tau * taus[1,j,i], rate = beta_tau, log = TRUE)) # posterior
      # print(p)
    }
    
    
    
    if (log(runif(1)) < min(0, p)) {
      gammas[,,i+1] <- gamma_star
      taus[,,i+1] <- tau_star
      acc <- acc + 1
    } 
    else {
      gammas[,,i+1] <- gammas[,,i]
      taus[,,i+1] <- taus[,,i]
    }
    print(p)
    print(i)
  }
  print(acc/iterations)

  thetas <- array(0, dim = c(iterations + 1, m + (m**2))) # converting it into this format 
  # in order to be able to return it and calculate the effective sample size
  
  for (i in 1:(iterations+1)){
    thetas[i,] <- as.vector(c(taus[,,i], as.vector(t(gammas[,,i]))))
  }
  return(thetas)
}

