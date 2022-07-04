# gamma - transition matrix
# mu - vector of (increasing) means of each state
# m - number of states
# delta - the initial probabilities of each state
# N - length of time series we are generating (uses T in the book but replace with N to 
# avoid conflict with boolean operator)

# MU's MUST BE IN INCREASING ORDER

HMM_Log_Likelihood <- function(x, Gamma, mu){ # returns the log likelihood for the HMM given mu and Gamma
  
  N <- length(x) # length of observed time series
  m <- dim(Gamma)[1] # number of states 
  U <- matrix(1, nrow = m, ncol = m)  # m x m matrix of ones
  delta <- matrix(1, nrow = 1, ncol = m) %*% solve(diag(m) - Gamma + U) # obtaining initial probabilities
  
  alpha <- delta * dnorm(x[1], mean = mu, sd = mu)
  lscale <- log(sum(alpha))
  alpha <- alpha/sum(alpha)
  
  for (i in 2:N){
    alpha <- alpha %*% Gamma * dnorm(x[i], mean = mu, sd = mu)
    lscale <- lscale + log(sum(alpha))
    alpha <- alpha/sum(alpha)
  }
  
  return(lscale)

}

