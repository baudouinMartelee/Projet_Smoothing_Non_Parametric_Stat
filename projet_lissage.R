# PROJET SMOOTHING TECHNIQUE

rm(list = ls()) #We use this to clear all variables in memory
library(stats) #load packages (stats normally loaded by default)
library(datasets) #We use some Inbuilt Data later

# We create our data
set.seed(22.1)
n <- 100
X <- seq(from = 1/n, to = 1, by = 1/n)
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05))
theta_squared <- function(x) 2 + sin(2* pi *x)
epsilon <- rnorm(n, mean = 0, sd = 1)
set.seed(22.2)
Y <- m(X) + 0.5 * sqrt(theta_squared(X)) * epsilon
Y
# Plot of the true function
x <- seq(0,1, length.out = 100)
plot(X, Y, pch = "+")
rug(X)
lines(X, m(X), col=3)

# Plot the true  variance theta^2
plot(X,theta_squared(X), pch = "+" ,type = "l", ylab = "Variance")


################################################
#   CHOOSING A GOOD BANDWITH
################################################

#Rule of Thumb estimator
h2 <- bw.nrd(X)

#cross validation estimator
h3 <- bw.ucv(X)

# SHeater & Jones bandwith selector 
h_dpi <- bw.SJ(X,method = "dpi")
h_ste <- bw.SJ(X)

plot(range(X), c(0,1), type="n")
lines(density(X),col = 1)
lines(density(X, n=100, from=min(X), to=max(X), bw="nrd", kernel = "gaussian"), col=2)
lines(density(X, n=100, from=min(X), to=max(X), bw="ucv", kernel = "gaussian"), col=4)
lines(density(X, n=100, from=min(X), to=max(X), bw="SJ", kernel = "gaussian"), col=5)
legend("center", legend=c("default", "nrd", "ucv", "SJ"),col = c(1,2, 4, 5), lty=c(1,1,1), rug(X),cex = 0.5)
rug(X)
################################################
#   COMPARAISON OF DIFFERENT ESTIMATORS FOR FUNCTION M(X)
################################################

# A KERNEL ESTIMATOR
############################

h = h_ste
# Nadaraya-Watson  estimator
nwRegEstim <- function(x, X, Y, h, K) sum(Y * K((x - X)/h))/sum(K((x - X)/h))

# Kernels CDF
Knorm <- function(u) dnorm(u) #Gaussian kernel
Kunif <- function(u) (abs(u) <= 1) * 0.5 #Uniform kernel
Kepan <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1) #Epanechnikov kernel
Ktria <- function(u) (1 - abs(u)) * (abs(u) <= 1) #Triangle kernel


NWnormEst <- sapply(x, function(x) nwRegEstim(x, X, Y, h, Knorm))

plot(X, Y, pch="+", col = "black")
lines(X, NWnormEst, type = "l", col = 4)
lines(x, m(x), type = "l", col = 3)
legend("topleft", legend = c("NW", "True"),
       col = c( 4, 3),
       lty = c( 1, 1), cex = 0.5)


# A NON-KERNEL ESTIMATOR (SPLINE SMOOTHING)
#####################################


################################################
#  ESTIMATORS FOR VARIANCE FUNCTION
################################################

# FAN et YAO
epsilonEstim = Y - NWnormEst

varRegEstim <- function(x, X, Y, h, K) sum(epsilonEstim^2 * K((x - X)/h))/sum(K((x - X)/h))

varNormRegEstim <- sapply(x,function(x)  varRegEstim(x, X, Y, h, Knorm))

VarianceX <- theta_squared(X)

plot(X, VarianceX, type = "l", col = "black")
plot(X, varNormRegEstim, type = "l", col = "black")
lines(X, varNormRegEstim , type = "l", col = 4)
lines(x, varNormRegEstim, type = "l", col = 3)
legend("topleft", legend = c("TRUE VARIANCE", "VARIANCE ESTIMATOR"),
       col = c( 4, 3),
       lty = c( 1, 1), cex = 0.5)

################################################
# plot the estimate function
################################################




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













