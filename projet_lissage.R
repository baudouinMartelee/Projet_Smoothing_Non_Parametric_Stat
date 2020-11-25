# PROJET SMOOTHING TECHNIQUE


#install.packages("lokern")
library("lokern")
rm(list = ls())

# We create our data
set.seed(22.01)
n <- 100

X <- seq(from = 1/n, to = 1, by = 1/n)
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05))
variance <- function(X) { 
  var <- 2 + sin(2*pi*X) 
  return(var)
}
variance
sigma_squared <- function(x) 2 + sin(2*pi*x)
epsilon <- rnorm(n, mean = 0, sd = 1)
set.seed(22.02)
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

# In the course, we have seen that global and local bandwith + rule of thumb are good fits

#global bandwith
h_gl <- glkerns(X, Y, hetero = TRUE, sig = sigma_squared(X), is.rand = F)
h_gl$bandwidth
#Local bandwith
h_lo <- lokerns(X, Y, hetero = TRUE, sig = sigma_squared(X), n.out = 100)
h_lo$bandwidth
#Rule of Thumb estimator
h_RT <- bw.nrd(X)

#Pour  prendre la bandwidth

################################################
#   COMPARAISON OF DIFFERENT ESTIMATORS FOR FUNCTION M(X)
################################################

# A KERNEL ESTIMATOR
############################

# Nadaraya-Watson  estimator
nwRegEstim <- function(x, X, Y, h, K) sum(Y * K((x - X)/h))/sum(K((x - X)/h))

# Kernels CDF
Knorm <- function(u) dnorm(u) #Gaussian kernel


NWnormEst_h_gl <- sapply(x, function(x) nwRegEstim(x, X, Y, h_gl$bandwidth, Knorm))
NWnormEst_h_lo <- sapply(x, function(x) nwRegEstim(x, X, Y, h_lo$bandwidth, Knorm))
NWnormEst_h_RT <- sapply(x, function(x) nwRegEstim(x, X, Y, h_RT, Knorm))

plot(X, Y, pch=20, col = "black", main = "Comparaison m(x) avec différents bandwidths pour Nadarya-Watson")
lines(x, m(x), type = "l", col = 3)
lines(X, NWnormEst_h_gl, type = "l", col = 4)
lines(X, NWnormEst_h_lo, type = "l", col = 2)
lines(X, NWnormEst_h_RT, type = "l", col = 5)
legend("bottomright", legend = c("True m(x)", "NW, h global = 0.05","NW, h local","NW, h Thumb of rule = 0.12"),
       col = c( 3, 4, 2, 5),
       lty = c( 1, 1))

#Gass and Muller estimator
cdfKnorm <- function(u) pnorm(u) #Gaussian kernel cdf


gmRegEstim <- function(x, X, Y, h, cdfK) { 
  IntK <- NULL
  n <- length(X)
  S <- -Inf
  for (i in 1:(n - 1)) {
    S <- rbind(S, 0.5 * (X[i] + X[i + 1])) 
    }
    S <- rbind(S, Inf)
    l <- (x - S)/h 
    for (i in 1:n) {
      IntK[i] <- cdfK(l[i]) - cdfK(l[i + 1]) 
      }
      sum(IntK * Y) 
}

GMnormEst_h_gl <- sapply(x, function(x) gmRegEstim(x, X, Y, h_gl$bandwidth, cdfKnorm ))
GMnormEst_h_lo <- sapply(x, function(x) gmRegEstim(x, X, Y, h_lo$bandwidth, cdfKnorm ))
GMnormEst_h_RT <- sapply(x, function(x) gmRegEstim(x, X, Y, h_RT, cdfKnorm ))


plot(X, Y, pch=20, col = "black", main = "Comparaison m(x) avec différents bandwidths pour Gasser & Muller")
lines(x, m(x), type = "l", col = 3)
lines(X, GMnormEst_h_gl, type = "l", col = 4)
lines(X, GMnormEst_h_lo, type = "l", col = 2)
lines(X, GMnormEst_h_RT, type = "l", col = 5)
legend("bottomright", legend = c("True m(x)", "GM, h global = 0.05","GM, h local","GM, h Thumb of rule = 0.12"),
       col = c( 3, 4, 2),
       lty = c( 1, 1), cex = 0.8)

# A NON-KERNEL ESTIMATOR (SPLINE SMOOTHING)
#####################################


spl0 <- smooth.spline(X, Y) # default (Cross-Validation)
spl2 <- smooth.spline(X, Y, spar=h_gl$bandwidth)
spl3 <- smooth.spline(X, Y, spar=h_RT)# smoothing parameter = h_RT

## Visualization
plot(X, Y, pch=20, col = "black", main="Spline Smoothing")
lines(x, m(x), col=1, lty=2)
lines(predict(spl0,x), col=2)
lines(predict(spl2,x), col=3)
lines(predict(spl3,x), col=4)
legend("bottomright", legend=c("True m(x)","spline default (Cross-Validation)","spline(spar=h_gl)","spline(spar=h_RT)"),
       col = c(1,2,3,4),lty=c(2,1,1,1,1), cex=1)

################################################
#  ESTIMATORS FOR VARIANCE FUNCTION
################################################

# FAN YAO  use 2 differents smoothing parameters
# 1for estimate m(x) and 1 for estimate sigma2.

estimateVariance <- function(x, X, Y, h1, h2, RegEstim, K){
  epsilonEstim = 2*(Y - sapply(x, function(x) RegEstim(x, X, Y, h1, K)))
  fanAndYaoEstim <- sum(epsilonEstim^2 * K((x - X)/h2))/sum(K((x - X)/h2))
  return(fanAndYaoEstim)
}

varianceNW1 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, nwRegEstim, Knorm))
varianceNW2 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim, Knorm))
varianceNW2 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim, Knorm))
varianceNW3 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_gl$bandwidth, nwRegEstim, Knorm))
varianceNW4 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_lo$bandwidth, nwRegEstim, Knorm))
varianceNW5 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_RT, nwRegEstim, Knorm))
varianceNW6 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_gl$bandwidth, nwRegEstim, Knorm))
varianceNW7 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_lo$bandwidth, nwRegEstim, Knorm))

varianceGM1 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim, cdfKnorm))
varianceGM2 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_lo$bandwidth, gmRegEstim, cdfKnorm))
varianceGM3 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_gl$bandwidth, gmRegEstim, cdfKnorm))
varianceGM32 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_lo$bandwidth, gmRegEstim, cdfKnorm))
varianceGM4 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_RT, gmRegEstim, cdfKnorm))
varianceGM5 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_gl$bandwidth, gmRegEstim, cdfKnorm))
varianceGM6 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_lo$bandwidth, gmRegEstim, cdfKnorm))

VarianceX <- sigma_squared(X)

# PLOT VARIANCE
plot(X, VarianceX, type = "l",ylim = c(0,15), col = "black", main= paste("Comparaison de Variance en utilisant la méthode Nadarya-Watson"))
lines(x, varianceNW1, type = "l", col = 2)
lines(x, varianceNW2, type = "l", col = 3)
lines(x, varianceNW3, type = "l", col = 4)
lines(x, varianceNW32, type = "l", col = 5)


legend("topright", legend = c("True Variance", "variance NW, h1 = gl,h2 = lo","variance NW, h1 = gl,h2 = RT","variance NW, h1 = lo,h2 = gl","variance NW, h1 = lo,h2 = RT","variance NW, h1 = RT,h2 = gl","variance NW, h1 = RT,h2 = lo"),
       col = c( 1,2,3,4,5),
       lty = c( 1, 1), cex = 0.5)

plot(X, VarianceX, type = "l",ylim = c(0,30), col = "black", main= paste("Comparaison de Variance en utilisant la méthode Gasser and Muller"))
lines(x, varianceGM1, type = "l", col = 2)
lines(x, varianceGM2, type = "l", col = 3)
lines(x, varianceGM3, type = "l", col = 4)
lines(x, varianceGM32, type = "l", col = 5)

legend("topright", legend = c("True Variance", "variance NW, h1 = gl,h2 = lo","variance NW, h1 = gl,h2 = RT","variance NW, h1 = lo,h2 = gl","variance NW, h1 = lo,h2 = RT","variance NW, h1 = RT,h2 = gl","variance NW, h1 = RT,h2 = lo"),
       col = c( 1,2,3,4,5),
       lty = c( 1, 1), cex = 0.5)

# variance for spline smoothing

estimateVarianceSMCV <- function(x, X, Y,h2, K){
  epsilonEstim = 2*(Y - smooth.spline(X, Y)$y)
  cat("epsilon = ",epsilonEstim)
  fanAndYaoEstim <- sum(epsilonEstim^2 * K((x - X)/h2))/sum(K((x - X)/h2))
  return(fanAndYaoEstim)
}

estimateVarianceSM <- function(x, X, Y, h1, h2, K){
  
  epsilonEstim = 2*(Y - smooth.spline(X, Y, spar = h1)$y)
  cat("epsilon = ",epsilonEstim)
  fanAndYaoEstim <- sum(epsilonEstim^2 * K((x - X)/h2))/sum(K((x - X)/h2))
  return(fanAndYaoEstim)
}

varianceSM_CV_h_gl <-sapply(x, function(x) estimateVarianceSMCV(x, X, Y, h_gl$bandwidth, Knorm))
varianceSM_CV_h_lo <-sapply(x, function(x) estimateVarianceSMCV(x, X, Y, h_lo$bandwidth, Knorm))
varianceSM_CV_h_RT <-sapply(x, function(x) estimateVarianceSMCV(x, X, Y, h_RT, Knorm))

varianceSM_hgl_hgl <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, Knorm))

varianceSM_hgl_hlo <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_lo$bandwidth, Knorm))
varianceSM_hgl_hrt <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_RT, Knorm))
varianceSM_hgl_hgl <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, Knorm))
varianceSM_hgl_hgl <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, Knorm))
varianceSM_hgl_hgl <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, Knorm))

epsilonEstim = 2*(Y -  sapply(x, function(x) smooth.spline(X, Y))$y)
varianceSM1 <- sum((epsilonEstim^2) * spl0 )/sum(spl0)

VarianceX <- sigma_squared(X)

plot(X, VarianceX, type = "l",ylim = c(-5,20), col = "black", main="Variance using Spline Smoothing")
lines(x, varianceSM_CV_h_gl, col=2)
lines(x, varianceSM_CV_h_lo, col=3)
lines(x, varianceSM_CV_h_RT, col=4)
legend("topright", legend=c("spline(spar=h_gl)","spline(spar=h_RT)"),
       col = c(3,5),lty=c(1,1), cex=0.8)

# in this notation the first h is for estimate m(x) and the second sigma2
# for example varNormRegEstim_h_gl_h_lo  we use h_gl for estimate m(x) and h_lo to estimate h_lo
varNormRegEstim_h_gl_h_RT <- sapply(x,function(x)  varRegEstim(x, X, Y, h_RT, Knorm, epsilonEstim_h_gl))
varNormRegEstim_h_RT_h_gl <- sapply(x,function(x)  varRegEstim(x, X, Y, h_gl$bandwidth, Knorm, epsilonEstim_h_RT))
VarianceX <- sigma_squared(X)


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









