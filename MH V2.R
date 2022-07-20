HMM_MH <- function(x, sigma, iterations, Gamma0, tau0, alpha_tau, beta_tau, alpha_gamma, beta_gamma){
  m <- dim(Gamma0)[1] # number of states in HMM
  acc_gam <- rep(0, m) # acceptance rate of each row of Gamma and Taus
  acc_tau <- rep(0, m)
  Thetas <- array(NA, dim = c(m+1, m, iterations + 1)) # empty 3D array to store Tau and Gamma simulations
  Thetas[1,,1] <- tau0 # setting first row as the Taus
  Thetas[2:(m+1),,1] <- Gamma0 # remaining rows as the initial gamma matrix
  
  for (i in 1:iterations) {
    if (i %% 1000 == 0){print(i)}
    
    # simulate each row of gamma independently from a gamma distribution
    
    proposal <- Thetas[,,i]
    current <- Thetas[,,i]
    
    for (j in 2:(m+1)) { # for the Gamma Matrix
      for (k in 1:m){
        # proposal[j, k] <- rgamma(n = 1, shape = alpha_gamma + current[j,k]/5, rate = 1)
        proposal[j, k] <- rgamma(n = 1, shape = alpha_gamma + current[j,k], rate = 1)
      }
      
      proposal[j, ] <- proposal[j, ] / sum(proposal[j, ])
      
      p <- HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(proposal[2:(m+1),]), mu = as.vector(cumsum(proposal[1,])) ) - # likelihood of proposal
        HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(current[2:(m+1),]), mu = as.vector(cumsum(current[1,])) )   # likelihood of current
      
      # a0_1 <- (m*alpha_gamma) + sum(proposal[j, ])
      # a0_2 <- (m*alpha_gamma) + sum(current[j,])
      for (k in 1:m) {
        
        p <- p + 
          dbeta(proposal[j, k], shape1 = alpha_gamma, shape2 =  beta_gamma, log = TRUE) + # prior
          dbeta(current[j,k], shape1 = alpha_gamma + proposal[j, k], shape2 = beta_gamma, log = TRUE) - # proposal
          dbeta(current[j,k], shape1 = alpha_gamma, shape2 = beta_gamma , log = TRUE) - # prior
          dbeta(proposal[j, k], shape1 = alpha_gamma + current[j, k], shape2 = beta_gamma, log = TRUE) # proposal
        
        # (dbeta(Thetas[j,k,i], shape1 = alpha_gamma + proposal[j, k], shape2 = a0_1 - alpha_gamma - proposal[j, k], log = TRUE)) +  # proposal
        # (dbeta(proposal[j, k], shape1 = alpha_gamma, shape2 = (m-1)*alpha_gamma, log = TRUE)) - # prior
        # (dbeta(proposal[j, k], shape1 = alpha_gamma + current[j,k], shape2 = a0_2 - alpha_gamma - current[j,k], log = TRUE)) -  # proposal
        # (dbeta(Thetas[j,k,i], shape1 = alpha_gamma, shape2 = (m-1)*alpha_gamma , log = TRUE)) # prior
      }
      
      if (min(p, 0) < log(runif(1))) { # if we fail to accept the proposal
        proposal[j,] <- current[j,] # use previous
      }
      else {
        acc_gam[j-1] <- acc_gam[j-1] + 1
        current[j, ] <- proposal[j,]
      }
      
      
    } # end of gamma
    
    p <- 0 # setting p to 0 for the Tau simulation [TAU]
    
    for (k in 1:m) {     # Simulating the Taus using gamma distribution
      proposal[1, k] <- rnorm(1, mean = current[1, k], sd = 2*beta_tau[k])
      
      p <- p +
        dnorm(proposal[1,k], mean = alpha_tau[k], sd = 2*beta_tau[k], log = TRUE) + # prior
        dnorm(current[1,k], mean =  proposal[1, k], sd = 2*beta_tau[k], log = TRUE) - # proposal
        dnorm(current[1,k], mean = alpha_tau[k], sd = 2*beta_tau[k], log = TRUE) - # prior
        dnorm(proposal[1, k], mean = current[1,k], sd = 2*beta_tau[k], log = TRUE) # proposal   
      + HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(proposal[2:(m+1),]), mu = as.vector(cumsum(proposal[1,])))  - # likelihood of proposal
        HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(current[2:(m+1),]), mu = as.vector(cumsum(current[1,] )))   # likelihood of current
      
      if (min(p, 0) < log(runif(1))) { # if we fail to accept the proposal
        proposal[1, k] <- current[1,k] # use previous
      }
      else {
        acc_tau[k] <- acc_tau[k] + 1
        current[1, k] <- proposal[1, k]
        
      }
    } # end of tau
    
    Thetas[,,i+1] <- proposal
    
  }
  print("Gamma acceptance:")
  print( acc_gam / iterations)
  print("Tau acceptance:")
  print( acc_tau / iterations)
  return(Thetas)
}