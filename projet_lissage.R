# PROJET LISSAGE   

rm(list = ls()) #We use this to clear all variables in memory
library(stats) #load packages (stats normally loaded by default)
library(datasets) #We use some Inbuilt Data later

set.seed(22.1)
#Size of our data
n <- 100
# Sequence of n numbers between 0 and 1 . These points are fixed
X <- seq(from = 1/n, to = 1, by = 1/n)

# this is the m(x) function 
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05))

# theta squared function 
theta_squared <- function(x) 2 + sin(2* pi *x)

#Epsilon 
epsilon <- rnorm(n, mean = 0, sd = 1)
#function Y 
set.seed(22.2)
Y <- m(X) + 0.5 * sqrt(theta_squared(X)) * epsilon

# Plot of the true function
plot(X, Y, pch = "+")
rug(X)
x <- seq(0,1, length.out = 300)
lines(x, m(x), col=2)
lines(X,m(X), col=4)

################################################
#   COMPARAISON OF DIFFERENT ESTIMATORS FOR FUNCTION M(X)
################################################

# bandwith h for the estimation
h = 0.02
# Nadaraya-Watson  estimator
nWEstim <- function(x, X, Y, h, K) sum(Y * K((x - X)/h))/sum(K((x - X)/h))

# Kernels CDF
Knorm <- function(u) dnorm(u) #Gaussian kernel
Kunif <- function(u) (abs(u) <= 1) * 0.5 #Uniform kernel
Kepan <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1) #Epanechnikov kernel
Ktria <- function(u) (1 - abs(u)) * (abs(u) <= 1) #Triangle kernel


NWnormEst <- sapply(x, function(x) nwRegEst(x, X, Y, h, Knorm))

plot(X, Y, pch="+", col = "black")
lines(x, NWnormEst, type = "l", col = 4)
lines(x, m(x), type = "l", col = 3)
legend("topleft", legend = c("NW", "True"),
       col = c( 4, 3),
       lty = c( 1, 1), cex = 0.5)



################################################
#   MONTE CARLO SIMULATION 
################################################

#Here I consider a fixed n
#and I compute MSE for a sequence of points x's
Knorm <- function(u) dnorm(u) #Gaussian kernel
## NW regression estimator
nwRegEst <- function(x, X, Y, h, K) sum(Y * K((x - X)/h))/sum(K((x - X)/h))
# monte carlo estimation
K <- 1000 # number of replications
n <- 400 # number of observations
h=0.05
x <- seq(from = 0, to = 1, length = 20) #a sequence of points
m <- function(x) (sin(2 * pi * x^3))^3 #true regression model
NWnormEst <- matrix(NA, nrow = K, ncol = length(x))
for (k in 1:K) {
  set.seed(k)
  X <- runif(n, 0, 1)
  set.seed(k + 1)
  Y <- m(X) + 0.5 * rnorm(n, 0, 1)
  NWnormEst[k,] <- sapply(x, function(x) nwRegEst(x, X, Y, h, Knorm))
}
NW<- colMeans(NWnormEst)
bias = NW - m(x)
SE = matrix(NA, nrow = K, ncol = length(x))
for (i in 1:length(x)) {
  SE[,i] = (NWnormEst[,i]-m(x[i]))^2
}
MSE = colMeans(SE)
MSSE = mean(MSE) # MSSE as an estimation of MISE in simulation













