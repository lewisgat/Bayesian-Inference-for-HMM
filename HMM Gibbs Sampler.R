# gamma - transition matrix
# mu - vector of (increasing) means of each state
# m - number of states
# delta - the initial probabilities of each state
# N - length of time series we are generating (uses T in the book but replace with N to 
# avoid conflict with boolean operator)

# MU's MUST BE IN INCREASING ORDER

Sample_Path_Generator <- function(x, Gamma, mu){
  N <- length(x) # length of observed time series
  U <- matrix(1, nrow = m, ncol = m)  # m x m matrix of ones
  C <- rep(0, m)
  delta <- matrix(1, nrow = 1, ncol = m) %*% solve(diag(m) - Gamma + U) # obtaining initial probabilities
  Alphas <- matrix(NA, nrow = N, ncol = m)
  
  alpha <- delta * dnorm(x[1], mean = mu, sd = mu)
  Alphas[1,] <- alpha
  for (i in 2:N){
    alpha <- alpha %*% Gamma * dnorm(x[i], mean = mu, sd = mu)
    Alphas[i,] <- alpha
  }
  
  C[N] <- sample(c(1:m), size=1, prob=as.vector(alpha/sum(alpha)))
  
  for (i in (N-1):1){
    alpha <- Alphas[i,]
    alpha <- alpha * Gamma[,C[i+1]]
    C[i] <- sample(c(1:m), size=1, prob=as.vector(alpha/sum(alpha)))
  }
  return(C)
}

Transition_Counts <- function(C, m){
  N <- length(C)
  transitions <- matrix(0, nrow = m, ncol = m)
  
  for (i in 1:(N-1)){
    transitions[C[i],C[i+1]] <- transitions[C[i],C[i+1]] + 1
  }
  return(transitions)
}

Gamma_Simulator <- function(C, m){
  
  transitions <- Transition_Counts(C, m)
  
  Gamma <- matrix(0, nrow = m, ncol = m) # empty Gamma matrix
  
  for (i in 1:m) { # for the Gamma Matrix
    for (j in 1:m){
      Gamma[i, j] <- rgamma(n = 1, shape =  1 + transitions[i,j], rate = 1) # simulate each entry 
    }
    Gamma[i, ] <- Gamma[i, ] / sum(Gamma[i, ]) # so each row sums to 1, simulates a Dirichlet Distribution
  }
  return(Gamma)
}

Mu_Simulator <- function(C, x, m, sigma){
  
  mu <- rep(0, m)
  # S <- rep(0, m)
  # n <- rep(0, m)
  
  xi <- (min(x) + max(x))/2
  R <- max(x) - min(x) 
  kappa <- 1/(R**2)
  
  for (i in 1:m){
    index <- which(C == i)
    S <- sum(x[index]) # sum of number of obvs under each state
    n <- length(index) # number of states 
    
    mu[i] <- rnorm(1, mean = (S + kappa*xi*(sigma**2))/(n+kappa*(sigma**2)), # simulating Mu as Ryden did in his paper
                   sd = sqrt( (sigma**2)/(n + kappa*(sigma**2)) ))
  }
  
  return(mu)
  
}

Gibbs_Sampler <- function(Gamma0, mu0, x, iterations, sigma){
  m <- length(mu0) # number of states 
  Thetas <- array(NA, dim = c(m+1, m, iterations + 1))
  Thetas[1,,1] <- mu0 # setting first row as the Taus
  Thetas[2:(m+1),,1] <- Gamma0 # remaining rows as the initial gamma matrix
  
  for (i in 1:iterations){
    mu_prior <- as.vector(Thetas[1,,i])
    Gamma_prior <- as.matrix(Thetas[2:(m+1),,i])
    C <- Sample_Path_Generator(x, Gamma = Gamma_prior, mu = mu_prior)
    Thetas[2:(m+1),,i+1] <- Gamma_Simulator(C, m)
    Thetas[1,,i+1] <- Mu_Simulator(C, x, m, sigma)
  }
  return(Thetas)
}
