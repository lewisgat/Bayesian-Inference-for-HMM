HMM_MH <- function(x, iterations, Gamma0, tau0, alpha_gamma, alpha_tau, beta_tau, beta_gamma){
  m <- dim(Gamma0)[1] # number of states in HMM
  acc <- rep(0, m+1) # acceptance rate of each row of Gamma and Taus
  
  Thetas <- array(NA, dim = c(m+1, m, iterations + 1)) # empty 3D array to store Tau and Gamma simulations
  Thetas[1,,1] <- tau0 # setting first row as the Taus
  Thetas[2:(m+1),,1] <- Gamma0 # remaining rows as the initial gamma matrix
  
  for (i in 1:iterations) {
    print(i)

      # simulate each row of gamma independently from a gamma distribution
      
    proposal <- Thetas[,,i]
      
      for (j in 2:(m+1)) { # for the Gamma Matrix
        for (k in 1:m){
          proposal[j, k] <- rgamma(n = 1, shape =  alpha_gamma + Thetas[j,k,i], rate = 1)
        }
        
        proposal[j, ] <- proposal[j, ] / sum(proposal[j, ])
        
        p <- HMM_Log_Likelihood(x, as.matrix(proposal[2:(m+1),]), as.vector(cumsum(proposal[1,])) ) - # likelihood of proposal
             HMM_Log_Likelihood(x, as.matrix(Thetas[2:(m+1),,i]), as.vector(cumsum(Thetas[1,,i])) )   # likelihood of current
       
        # a0_1 <- (m*alpha_gamma) + sum(proposal[j, ])
        # a0_2 <- (m*alpha_gamma) + sum(Thetas[j,,i]) 
        for (k in 1:m) {
          
          p <- p + 
          # dbeta(proposal[j, k], shape1 = alpha_gamma, shape2 = (m-1)*alpha_gamma, log = TRUE) + # prior
          # dbeta(Thetas[j,k,i], shape1 = alpha_gamma + proposal[j, k], shape2 = a0_1 - alpha_gamma - (proposal[j, k]), log = TRUE) + # proposal
          # dbeta(Thetas[j,k,i], shape1 = alpha_gamma, shape2 = (m-1)*alpha_gamma , log = TRUE) - # prior
          # dbeta(proposal[j, k], shape1 = alpha_gamma + Thetas[j,k,i], shape2 = a0_2 - alpha_gamma - (Thetas[j,k,i]), log = TRUE) # current
            
          dbeta(proposal[j, k], shape1 = alpha_gamma, shape2 =  beta_gamma, log = TRUE) + # prior
          dbeta(Thetas[j,k,i], shape1 = alpha_gamma + proposal[j, k], shape2 = beta_gamma, log = TRUE) + # proposal
          dbeta(Thetas[j,k,i], shape1 = alpha_gamma, shape2 = beta_gamma , log = TRUE) - # prior
          dbeta(proposal[j, k], shape1 = alpha_gamma + Thetas[j,k,i], shape2 = beta_gamma, log = TRUE) # current
        
        }
        
        if (min(p, 0) < log(runif(1))) { # if we fail to accept the proposal
          proposal[j,] <- Thetas[j,,i] # use previous
        }
        else {
          acc[j] <- acc[j] + 1
        }
        
        p <- 0 # setting p to 0 for the Tau simulation
        
        for (k in 1:m) {     # Simulating the Taus using gamma distribution
          
          proposal[1, k] <- rgamma(1, shape = (beta_tau[k] + Thetas[1, k, i]), rate = beta_tau[k]) # draw Tau from Gamma Distribution
          
          
          p <- p +
            dgamma(proposal[1, k], shape = alpha_tau[k] , rate = beta_tau[k], log = TRUE) + # prior
            dgamma(Thetas[1,k,i], shape = beta_tau[k] + proposal[1, k], rate = beta_tau[k], log = TRUE)  - # proposal
            dgamma(Thetas[1,k,i], shape = alpha_tau[k], rate = beta_tau[k], log = TRUE) - # prior
            dgamma(proposal[1, k], shape = (beta_tau[k] + Thetas[1,k,i]), rate = beta_tau[k], log = TRUE) # proposal
          
        }
        
        p <- p + HMM_Log_Likelihood(x, as.matrix(proposal[2:(m+1),]), as.vector(cumsum(proposal[1,])) ) - # likelihood of proposal
          HMM_Log_Likelihood(x, as.matrix(Thetas[2:(m+1),,i]), as.vector(cumsum(Thetas[1,,i])) )   # likelihood of current
        
        if (min(p, 0) < log(runif(1))) { # if we fail to accept the proposal
          proposal[1,] <- Thetas[1,,i] # use previous
        }
        else {
          acc[1] <- acc[1] + 1
        }
        
        
        Thetas[,,i+1] <- proposal
      }
  }
  
  print( acc / iterations)
  return(Thetas)
}






