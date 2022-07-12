# gamma - transition matrix
# mu - vector of (increasing) means of each state
# m - number of states
# delta - the initial probabilities of each state
# N - length of time series we are generating (uses T in the book but replace with N to 
# avoid conflict with boolean operator)

# MU's MUST BE IN INCREASING ORDER


DataGenerator2 <- function(gamma, mu, N, sigma)
{
  
  m <- dim(gamma)[1] # number of states 
  
  U <- matrix(1, nrow = m, ncol = m)  # m x m matrix of ones
  
  delta <- matrix(1, nrow = 1, ncol = m) %*% solve(diag(m) - gamma + U) # solving this equation to find delta
  # delta <- eigen(gamma)$vector[,1]/sum(eigen(gamma)$vector[,1])        
  
  # for ease of notation in the for loop, we will denote delta and 
  # all future transition probabilities transition_probs
  transition_probs <- delta
  
  df = data.frame(matrix(vector(), N, 2, dimnames=list(c(), c("State","Data")))) # creating an empty N x 2 
  # dataframe to store states and observations in 
  
  for (i in 1:N) {
    state <- sample(c(1:m), size=1, prob=as.vector(transition_probs)) # sample a state using transition 
    # probability, given previous state
    
    mu_state <- mu[state] # obtaining the mean associated with the new state we have sampled
    
    x <- rnorm(1, mean = mu_state, sd = sigma) # generating an observation using this state-dependent 
    # distribution
    df[i, ] <- c(state, x) # storing the state and the observation
    
    transition_probs <- gamma[state,] # obtaining the transition probabilities associated with this state
    # repeats the loop for N iterations to generate N observations 
  }
  return(df)
}

#####
# Testing the data generation matrix

gamma = matrix(c(1/3, 1/3, 1/3, 2/3, 0, 1/3, 1/2, 1/2, 0), nrow = 3, byrow = TRUE)
mu <- c(10, 50, 100)

df <- DataGenerator(gamma, mu, 100)
head(df)

hist(df[,2], breaks = 100)


