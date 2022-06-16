# gamma - transition matrix
# mu - vector of (increasing) means of each state
# m - number of states
# delta - the initial probabilities of each state
# N - length of time series we are generating (uses T in the book but replace with N to 
# avoid conflict with boolean operator)

# MU's MUST BE IN INCREASING ORDER

Emission_Matrix <- function(x, mu, m){ # function creates the emission probability matrix for a given observation
  P <- matrix(0, nrow = m, ncol = m) # m x m matrix of 0's
  for (i in 1:m){
    P[i,i] <- pnorm(x, mean = mu[i], sd = mu[i]) # fills in the diagonals using state dependent distribution
  }
  return(P)
}

HMM_Likelihood <- function(x, gamma, mu){ # function returns the likelihood for the HMM given mu and gamma
  
  N <- length(x) # length of observed time series
  m <- dim(gamma)[1] # number of states 
  U <- matrix(1, nrow = m, ncol = m)  # m x m matrix of ones
  delta <- matrix(1, nrow = 1, ncol = m) %*% solve(diag(m) - gamma + U) # obtaining initial probabilities
  
  # delta <- eigen(gamma)$vector[,1]
  L <- delta %*% Emission_Matrix(x[1], mu, m) # first iteration for observation 1
  
  for (s in 2:N){
    L <- L %*% gamma %*% Emission_Matrix(x[s], mu, m) # iterating for remaining observations
  }
  
  L <- L %*% matrix(1, nrow = m, ncol = 1) #summing up all the columns in L to obtain final Likelihood value
  
  return(L)
}


##### 
# Testing the likelihood
gamma = matrix(c(1/3, 1/3, 1/3, 2/3, 0, 1/3, 1/2, 1/2, 0), nrow = 3, byrow = TRUE)
mu <- c(10, 50, 100)
df <- DataGenerator(gamma, mu, 100)
x <- as.vector(df[,2])


HMM_Likelihood(x, gamma, mu)

