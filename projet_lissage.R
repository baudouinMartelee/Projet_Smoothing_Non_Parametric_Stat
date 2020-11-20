# PROJET SMOOTHING TECHNIQUE


install.packages("lokern")
library("lokern")
rm(list = ls())


# We create our data
set.seed(22.1)
n <- 100
X <- seq(from = 1/n, to = 1, by = 1/n)
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05))
sigma_squared <- function(x) 2 + sin(2* pi *x)
epsilon <- rnorm(n, mean = 0, sd = 1)
set.seed(22.2)
Y <- m(X) + 0.5 * sqrt(sigma_squared(X)) * epsilon
Y
# Plot of observations and m(x) function.
x <- seq(0,1, length.out = 100)
plot(X, Y, pch =20, lwd=0.2)
rug(X)
lines(X, m(X), col=3)
legend("top", legend = "m(x)", col = "green",lty=c(1,1),cex=1)

# Plot the true  variance theta^2
plot(X, sigma_squared(X), type = "l", ylab = "Variance sigma2")


################################################
#   CHOOSING A GOOD BANDWITH
################################################

#global bandwith
h_gl <- glkerns(X, Y, hetero = TRUE, sig = sigma_squared(X))
h_gl$bandwidth
#Local bandwith
h_lo <- lokerns(X, Y, hetero = TRUE, sig = sigma_squared(X))
h_lo$bandwidth

#Rule of Thumb estimator
h2 <- bw.nrd(X)

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


NWnormEst <- sapply(x, function(x) nwRegEstim(x, X, Y, h, Knorm))

plot(X, Y, pch="+", col = "black")
lines(X, NWnormEst, type = "l", col = 4)
lines(x, m(x), type = "l", col = 3)
legend("topright", legend = c("NW", "True"),
       col = c( 4, 3),
       lty = c( 1, 1), cex = 0.5)


# A NON-KERNEL ESTIMATOR (SPLINE SMOOTHING)
#####################################
spl0 <- smooth.spline(X, Y) # default (Cross-Validation)
spl1 <- smooth.spline(X, Y, spar=0) # smoothing parameter = 0 (equi. lambda = 0) (interp
# lambda = r * 256^(3*spar - 1)
# r = tr(X' W X) / tr(Sigma)
spl2 <- smooth.spline(X, Y, spar=0.5)# smoothing parameter = 0.5
spl3 <- smooth.spline(X, Y, spar=2) # smoothing parameter = 2 (linear function) #(equi. lambda = \infty)
spl4 <- smooth.spline(X, Y, spar=0.32)

## Visualization
plot(X, Y, main="Spline Smoothing")
lines(x, m(x), col=1, lty=2)
lines(predict(spl0,x), col=2)
lines(predict(spl1,x), col=3)
lines(predict(spl2,x), col=4)
lines(predict(spl3,x), col=5)
lines(predict(spl4,x), col=6)
legend("topleft", legend=c("True","spline default (Cross-Validation)", "spline(spar=0)","spline(spar=0.5)","spline(spar=2)","spline(spar=1)"),
       col = c(1,2,3,4,5,6),lty=c(2,1,1,1,1,1), cex=0.5)

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
legend("topright", legend = c("TRUE VARIANCE", "VARIANCE ESTIMATOR"),
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
K <- 1000# number of replications
n <- 100 # number of observations
h= h_ste
x <- seq(from = 0, to = 1, length = 300) #a sequence of points
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05)) #true regression model
theta_squared <- function(x) 2 + sin(2* pi *x) # true Variance
epsilon <- rnorm(n, mean = 0, sd = 1)

NWnormEst <- matrix(NA, nrow = K, ncol = length(x))

for (k in 1:K) {
  set.seed(k)
  X <- seq(from = 1/n, to = 1, by = 1/n)
  set.seed(k + 1)
  Y <- m(X) + 0.5 * sqrt(theta_squared(X)) * epsilon
  NWnormEst[k,] <- sapply(x, function(x) nwRegEst(x, X, Y, h, Knorm))
}
NW<- colMeans(NWnormEst)
bias = NW - m(x)
SE = matrix(NA, nrow = K, ncol = length(x))
for (i in 1:length(x)) {
  SE[,i] = (NWnormEst[,i]-m(x[i]))^2
}

variance1<-colMeans((NWnormEst-NW)^2)
plot(x,variance1 ,ylab="Variance",type="l",col="#C72C48")

MSE = colMeans(SE)
MSSE = mean(MSE) # MSSE as an estimation of MISE in simulation


plot(x, bias, type="l", col = 2)
abline(h=0)









