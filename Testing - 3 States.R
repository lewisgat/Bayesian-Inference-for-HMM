setwd("C:/Users/lewis/OneDrive/Desktop/Dissertation")
library(mcmcse)

#####
# Dispersed means and Gamma matrix, 200 Observations

mu <- c(5, 20, 50)
Gamma <- matrix(c(0.9, 0.1, 0,
                  0.1, 0.8, 0.1,
                  0.1, 0.2, 0.7), nrow = 3, ncol = 3, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 200)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/3, nrow = 3, ncol = 3)
tau0 <- c(10,20,30)
# acc rate  0.33956 0.27397 0.28696 0.25108
Thetas1 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(10,20,30), beta_tau = c(1,1,1), alpha_gamma = 0.3, beta_gamma = 1)

mean(Thetas1[4,3,10000:100001])

mu1_1 <- Thetas1[1,1,]
mu2_1 <- Thetas1[1,1,] + Thetas1[1,2,]
mu3_1 <- Thetas1[1,1,] + Thetas1[1,2,] + Thetas1[1,3,]
mean(mu1_1[10000:100001])
mean(mu2_1[10000:100001])
mean(mu3_1[10000:100001])

# means and first row of Gamma accurately estimated. The remaining two rows of Gamma are not estimated so well (bias). 

#####
# Dispersed means and Gamma matrix, 100 Observations

mu <- c(5, 20, 50)
Gamma <- matrix(c(0.9, 0.1, 0,
                  0.1, 0.8, 0.1,
                  0.1, 0.2, 0.7), nrow = 3, ncol = 3, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 100)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/3, nrow = 3, ncol = 3)
tau0 <- c(10,20,30)

# acc rate 0.46077 0.34829 0.38302 0.43843
Thetas2 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(10,20,30), beta_tau = c(1,1,1), alpha_gamma = 0.3, beta_gamma = 1)

mean(Thetas2[2,1,10000:100001])

mu1_2 <- Thetas2[1,1,]
mu2_2 <- Thetas2[1,1,] + Thetas2[1,2,]
mu3_2 <- Thetas2[1,1,] + Thetas2[1,2,] + Thetas2[1,3,]
mean(mu1_2[10000:100001])
mean(mu2_2[10000:100001])
mean(mu3_2[10000:100001])

# very similar to the one with 200 observations, first row of Gamma estimated accurately. Means are estimated to be larger than they
# actually are

##### 
# Not dispersed means, dispersed Gamma matrix 100 observations

mu <- c(2, 5, 7)
Gamma <- matrix(c(0.9, 0.1, 0,
                  0.05, 0.9, 0.05,
                  0.1, 0.1, 0.8), nrow = 3, ncol = 3, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 100)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/3, nrow = 3, ncol = 3)
tau0 <- c(3,3,3)
#acc rate  0.27188 0.50565 0.59215 0.63226
Thetas3 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(3,3,3), beta_tau = c(1,1,1), alpha_gamma = 0.5, beta_gamma = 1)

mean(Thetas3[2,1,10000:100001])

mu1_3 <- Thetas3[1,1,]
mu2_3 <- Thetas3[1,1,] + Thetas3[1,2,]
mu3_3 <- Thetas3[1,1,] + Thetas3[1,2,] + Thetas3[1,3,]
mean(mu1_3[10000:100001])
mean(mu2_3[10000:100001])
mean(mu3_3[10000:100001])

# overestimation of the means. Underestimation of the high probabilities in the matrix, overestimation of the small probabilities
plot(Thetas3[3,1,10000:100001])
# mixes well in all examples


#####
# Not dispersed means, dispersed Gamma matrix 200 observations
mu <- c(2, 5, 7)
Gamma <- matrix(c(0.9, 0.1, 0,
                  0.05, 0.9, 0.05,
                  0.1, 0.1, 0.8), nrow = 3, ncol = 3, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 200)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/3, nrow = 3, ncol = 3)
tau0 <- c(2,2,2)
# acc 0.35747 0.38980 0.53370 0.42436
Thetas4 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(2,2,2), beta_tau = 2*c(1,1,1), alpha_gamma = 0.3, beta_gamma = 1)

mu1_4 <- Thetas4[1,1,]
mu2_4 <- Thetas4[1,1,] + Thetas4[1,2,]
mu3_4 <- Thetas4[1,1,] + Thetas4[1,2,] + Thetas4[1,3,]
mean(mu1_4[10000:100001])
mean(mu2_4[10000:100001])
mean(mu3_4[10000:100001])

mean(Thetas4[1,2,10000:100001])

# small bias in the means
# Gamma matrix converges to a uniform matrix 
# seems that non-dispersed means leads to high bias Gamma

par(mfrow = c(4,3))

for (i in 1:4){
  for (j in 1:3){
    plot(Thetas4[i,j,])
  }
}

for (i in 1:4){
  for (j in 1:3){
    acf(Thetas4[i,j,])
  }
}


#####
# Dispersed Means, Gamma which encourages large amount of model switching. 100 observations

mu <- c(2, 20, 40)
Gamma <- matrix(c(0.1, 0.8, 0.1,
                  0.1, 0.1, 0.8,
                  0.8, 0.1, 0.1), nrow = 3, ncol = 3, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 100)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/3, nrow = 3, ncol = 3)
tau0 <- c(5, 15, 25)

# acc rate  0.40161 0.40348 0.33697 0.42012
Thetas5 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(5, 15, 25), beta_tau = c(1,1,1), alpha_gamma = 0.5, beta_gamma = 1)

mu1_5 <- Thetas5[1,1,]
mu2_5 <- Thetas5[1,1,] + Thetas5[1,2,]
mu3_5 <- Thetas5[1,1,] + Thetas5[1,2,] + Thetas5[1,3,]
mean(mu1_5[10000:100001])
mean(mu2_5[10000:100001])
mean(mu3_5[10000:100001])

mean(Thetas5[4,1,10000:100001])

#0.8 estimated to be about 0.6-0.7. 0.1's estimated to be around 0.1-0.2. Mu's estimated more accurately but off by a small amount

#####
# Dispersed Means, Gamma which encourages large amount of model switching. 200 observations

mu <- c(2, 20, 40)
Gamma <- matrix(c(0.1, 0.8, 0.1,
                  0.1, 0.1, 0.8,
                  0.8, 0.1, 0.1), nrow = 3, ncol = 3, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 200)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/3, nrow = 3, ncol = 3)
tau0 <- c(5, 15, 25)

# acc rate  0.21677 0.27920 0.30968 0.29082
Thetas6 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(5, 15, 25), beta_tau = c(1,1,1), alpha_gamma = 0.5, beta_gamma = 1)

mu1_6 <- Thetas6[1,1,]
mu2_6 <- Thetas6[1,1,] + Thetas6[1,2,]
mu3_6 <- Thetas6[1,1,] + Thetas6[1,2,] + Thetas6[1,3,]
mean(mu1_6[10000:100001])
mean(mu2_6[10000:100001])
mean(mu3_6[10000:100001])

mean(Thetas6[2,2,10000:100001])

# far more accurate once the number of observations have been doubled. Closely estimate all of the parameters in the Gamma and Tau
# largest improvement is in the Gamma matrix


# Flat_Thetas6<- array(NA, dim = c(100001, 12 ))
# for (i in 1:100001){
#   Flat_Thetas6[i,] <- as.vector(t(Thetas6[,,i]))
# }

par(mfrow = c(4,3))

for (i in 1:4){
  for (j in 1:3){
    plot(Thetas6[i,j,])
  }
}

for (i in 1:4){
  for (j in 1:3){
    acf(Thetas6[i,j,])
  }
}





