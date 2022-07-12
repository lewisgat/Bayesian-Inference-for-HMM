library(mcmcse)

#####
# Dispersed Means and Gamma 100 observations
mu <- c(5, 30)
Gamma <- matrix(c(0.9, 0.1,
                  0.2, 0.8), nrow = 2, ncol = 2, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 100)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/2, nrow = 2, ncol = 2)
tau0 <- c(10,15)
# acc rate  0.24526 0.18522 0.38127
Thetas1_2 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(10,15), beta_tau = c(1,1), alpha_gamma = 1, beta_gamma = 1)

mu11 <- Thetas1_2[1,1,]
mu21 <- Thetas1_2[1,1,] + Thetas1_2[1,2,]
mean(mu11[10000:100001])
mean(mu21[10000:100001])
mean(Thetas1_2[2,1,10000:100001])
mean(Thetas1_2[3,2,10000:100001])

# well estimated the larger probability, overestimated the smaller probability

#####
# Dispersed means, dispersed Gamma 200 observations

mu <- c(5, 30)
Gamma <- matrix(c(0.9, 0.1,
                  0.2, 0.8), nrow = 2, ncol = 2, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 200)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/2, nrow = 2, ncol = 2)
tau0 <- c(10,15)
# acc rate  0.16140 0.16392 0.21463
Thetas2_2 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(10, 15), beta_tau = c(1,1), alpha_gamma = 1, beta_gamma = 1)

mean(Thetas2_2[3,2,10000:100001])
mu12<- Thetas2_2[1,1,]
mu22 <- Thetas2_2[1,1,] + Thetas2_2[1,2,]
mean(mu12[10000:100001])
mean(mu22[10000:100001])
mean(Thetas2_2[2,1,10000:100001])
mean(Thetas2_2[3,2,10000:100001])

#####
# Dispersed means, dispersed Gamma 150 observations
mu <- c(5, 30)
Gamma <- matrix(c(0.9, 0.1,
                  0.2, 0.8), nrow = 2, ncol = 2, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 150)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/2, nrow = 2, ncol = 2)
tau0 <- c(10,15)
# acc rate  0.21287 0.17538 0.13813
Thetas2_3 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(10, 15), beta_tau = c(1,1), alpha_gamma = 1, beta_gamma = 1)

mean(Thetas2_3[3,2,10000:100001])
mu13<- Thetas2_3[1,1,]
mu23 <- Thetas2_3[1,1,] + Thetas2_3[1,2,]
mean(mu13[10000:100001])
mean(mu23[10000:100001])
mean(Thetas2_3[2,1,10000:100001])
mean(Thetas2_3[3,2,10000:100001])

Flat_Thetas2_3 <- array(NA, dim = c(100001, 6 ))
for (i in 1:100001){
  Flat_Thetas2_3[i,] <- as.vector(t(Thetas2_3[,,i]))
}
multiESS(Flat_Thetas2_3)
#####
# dispersed means, mixing encouraging Gamma matrix, 100 obvs
mu <- c(2, 20)
Gamma <- matrix(c(0.15, 0.85,
                  0.85, 0.15), nrow = 2, ncol = 2, byrow = TRUE)

df <- DataGenerator(Gamma, mu, 100)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/2, nrow = 2, ncol = 2)
tau0 <- c(5,10)

Thetas2_4 <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(5,10), beta_tau = c(1,1), alpha_gamma = 1, beta_gamma = 1)

mu14<- Thetas2_4[1,1,]
mu24 <- Thetas2_4[1,1,] + Thetas2_4[1,2,]
mean(mu14[10000:100001])
mean(mu24[10000:100001])
mean(Thetas2_4[2,2,10000:100001])
mean(Thetas2_4[3,1,10000:100001])

par(mfrow = c(3,2))

for (i in 1:3){
  for (j in 1:2){
    plot(Thetas2_4[i,j,])
  }
}

par(mfrow = c(3,2))

for (i in 1:3){
  for (j in 1:2){
    acf(Thetas2_4[i,j,])
  }
}


###### 
# TRYING WITH CONSTANT VARIANCE

mu <- c(5, 30)
Gamma <- matrix(c(0.9, 0.1,
                  0.2, 0.8), nrow = 2, ncol = 2, byrow = TRUE)

df <- DataGenerator2(Gamma, mu, 100, sigma = 2)
x <- as.vector(df[,2])

Gamma0 <- matrix(1/2, nrow = 2, ncol = 2)
tau0 <- c(10,15)
# acc rate  0.16363 0.19184 0.20305
Thetas1_2_fixed_sigma <- HMM_MH(x, 100000, Gamma0, tau0, alpha_tau = c(10,15), beta_tau = c(1,1), alpha_gamma = 1, beta_gamma = 1)

mu1_s <- Thetas1_2[1,1,]
mu2_s <- Thetas1_2[1,1,] + Thetas1_2[1,2,]
mean(mu1_s[10000:100001])
mean(mu2_s[10000:100001])
mean(Thetas1_2_fixed_sigma[2,1,10000:100001])
mean(Thetas1_2_fixed_sigma[3,2,10000:100001])
