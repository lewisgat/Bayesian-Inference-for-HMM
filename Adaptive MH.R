library(mvtnorm)
library(gtools)

Adaptive_MH <- function(x, sigma, iterations, Gamma0, mu0, Sigma0, epsilon, t0, alpha_gamma){
  m <- dim(Gamma0)[1] # number of states in HMM
  acc_gam <- rep(0, m) # acceptance rate of each row of Gamma and Taus
  # acc_tau <- rep(0, m)
  acc_mu <- 0
  Thetas <- array(NA, dim = c(m+1, m, iterations + 1)) # empty 3D array to store Tau and Gamma simulations
  Thetas[1,,1] <- mu0 # setting first row as the Taus
  Thetas[2:(m+1),,1] <- Gamma0 # remaining rows as the initial gamma matrix

  s <- (2.38**2)/m # for adaptive part

  for (i in 1:iterations) {
    # if (i %% 10000 == 0){print(i)}

    # simulate each row of gamma independently from a gamma distribution

    proposal <- Thetas[,,i]
    current <- Thetas[,,i]

    for (j in 2:(m+1)) { # for the Gamma Matrix

      proposal[j, ] <- rdirichlet(1, alpha = alpha_gamma+current[j, ])
      
      p <- HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(proposal[2:(m+1),]), mu = as.vector((proposal[1,])) ) - HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(current[2:(m+1),]), mu = as.vector((current[1,])) ) + log(ddirichlet(current[j, ], alpha =alpha_gamma+ proposal[j, ])) + log(ddirichlet(proposal[j, ], alpha = alpha_gamma )) - log(ddirichlet(proposal[j, ], alpha = alpha_gamma+current[j, ])) -  log(ddirichlet(current[j, ], alpha = alpha_gamma))

      if (min(p, 0) < log(runif(1))) { # if we fail to accept the proposal
        proposal[j,] <- current[j,] # use previous
      }
      else {
        acc_gam[(j-1)] <- acc_gam[(j-1)] + 1
        current[j, ] <- proposal[j,]
      }
    } # end of Gamma

    if (i < t0) {C <- Sigma0}

    # else { C <- (s * cov(t(Thetas[1,1:m,1:i])) ) + s*epsilon*diag(m) }
    else if (i >= t0) { C <- (s * cov(t(Thetas[1,1:m,1:i])) ) + s*epsilon*diag(m) }

      else {
        C <- (C * (i-1)/(i)) + (s/i) *( (i* colMeans(t(Thetas[1,1:m,1:(i-1)])) %*% t(colMeans(t(Thetas[1,1:m,1:(i-1)])))) - ((i+1)* colMeans(t(Thetas[1,1:m,1:i])) %*% t(colMeans(t(Thetas[1,1:m,1:i]))))
                   + (Thetas[1,1:m,i] %*% t(Thetas[1,1:m,i])) + (epsilon*diag(m)) )
    }
    proposal[1,] <- rmvnorm(1, mean = current[1,], sigma = C)


    p <- dmvnorm(current[1,], mean = proposal[1,], sigma = C ,log = TRUE) + dmvnorm(proposal[1,], mean = mu0, sigma = Sigma0, log = TRUE) - dmvnorm(proposal[1,], mean = current[1,], sigma = C, log = TRUE) - dmvnorm(current[1,], mean = mu0, sigma = Sigma0, log = TRUE) + HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(proposal[2:(m+1),]), mu = as.vector((proposal[1,])))- HMM_Log_Likelihood(x, sigma, Gamma = as.matrix(current[2:(m+1),]), mu = as.vector((current[1,] )))

    if (min(p, 0) < log(runif(1))) { # if we fail to accept the proposal
      proposal[1, ] <- current[1, ] # use previous
    }
    else {
      acc_mu[1] <- acc_mu[1] + 1
      current[1, ] <- proposal[1, ]
    } # end of tau

    Thetas[,,i+1] <- proposal

  }
  print("Gamma acceptance:")
  print( acc_gam / iterations)
  print("Mu acceptance:")
  print( acc_mu / iterations)
  return(Thetas)
}

