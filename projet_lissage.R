# PROJET SMOOTHING TECHNIQUE


#install.packages("lokern")
library("lokern")
rm(list = ls())

# We create our data
set.seed(22.01)
n <- 100

X <- seq(from = 1/n, to = 1, by = 1/n)
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05))
variance <- function(x)  2 + sin(2*pi*x) 
variance
sigma_squared <- function(x) 2 + sin(2*pi*x)
epsilon <- rnorm(n, mean = 0, sd = 1)
#set.seed(22.02)
Y <- m(X) + 0.5 * sqrt(sigma_squared(X)) * epsilon
Y
# Plot of observations and m(x) function.
x <- seq(0,1, length.out = 100)

plot(X, Y, pch =20, lwd=0.2)
rug(X)
lines(X, m(X), col=1)
legend("bottomright", legend = "m(x)", col = "black",lty=c(1,1),cex=1)

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
h_RT
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

plot(X, Y, pch=20, col = "black", main = "Estimateurs de m(x) avec Nadarya-Watson")
lines(x, m(x), type = "l", col = 1)
lines(x, NWnormEst_h_gl, type = "l", col = 2)
lines(x, NWnormEst_h_lo, type = "l", col = 3)
lines(x, NWnormEst_h_RT, type = "l", col = 4)
legend("bottomright", legend = c("True m(x)", "NW, h global = 0.05","NW, h local","NW, h Thumb of rule = 0.12"),
       col = c( 1, 2, 3, 4),
       lty = c( 1, 1), cex = 0.8)

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


plot(X, Y, pch=20, col = "black", main = "Estimateurs de m(x) avec Gasser & Muller")
lines(x, m(x), type = "l", col = 1)
lines(x, GMnormEst_h_gl, type = "l", col = 2)
lines(x, GMnormEst_h_lo, type = "l", col = 3)
lines(x, GMnormEst_h_RT, type = "l", col = 4)
legend("bottomright", legend = c("True m(x)", "GM, h global = 0.05","GM, h local","GM, h Thumb of rule = 0.12"),
       col = c(1, 2, 3, 4),
       lty = c( 1, 1), cex = 0.8)

# A NON-KERNEL ESTIMATOR (SPLINE SMOOTHING)
#####################################


spl0 <- smooth.spline(X, Y) # default (Cross-Validation)
spl2 <- smooth.spline(X, Y, spar=h_gl$bandwidth)
spl3 <- smooth.spline(X, Y, spar=h_RT)# smoothing parameter = h_RT

## Visualization
plot(X, Y, pch=20, col = "black", main="Estimateurs de m(x) avec Spline Smoothing")
lines(x, m(x), col=1)
lines(predict(spl0,x), col=2)
lines(predict(spl2,x), col=3)
lines(predict(spl3,x), col=4)
legend("bottomright", legend=c("True m(x)","spline default (Cross-Validation)","spline(spar=h_gl)","spline(spar=h_RT)"),
       col = c(1,2,3,4),lty=c(1,1,1,1,1), cex=0.8)

################################################
#  ESTIMATORS FOR VARIANCE FUNCTION
################################################

# FAN YAO  use 2 differents smoothing parameters
# 1for estimate m(x) and 1 for estimate sigma2.

estimateVariance <- function(x, X, Y, h1, h2, RegEstim, K){
  epsilonEstim = Y - sapply(x, function(x) RegEstim(x, X, Y, h1, K))
  fanAndYaoEstim <- sum((epsilonEstim^2) * K((x - X)/h2))/sum(K((x - X)/h2))
  return(fanAndYaoEstim)
}

varianceNW11 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, nwRegEstim, Knorm))
varianceNW12 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim, Knorm))
varianceNW12 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_RT, nwRegEstim, Knorm))
varianceNW21 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_gl$bandwidth, nwRegEstim, Knorm))
varianceNW22 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_lo$bandwidth, nwRegEstim, Knorm))
varianceNW23 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_RT, nwRegEstim, Knorm))
varianceNW31 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_gl$bandwidth, nwRegEstim, Knorm))
varianceNW32 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_lo$bandwidth, nwRegEstim, Knorm))
varianceNW33 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_RT, nwRegEstim, Knorm))


varianceGM11 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim, cdfKnorm))
varianceGM12 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_lo$bandwidth, gmRegEstim, cdfKnorm))
varianceGM12 <-sapply(x, function(x) estimateVariance(x, X, Y, h_gl$bandwidth, h_RT, gmRegEstim, cdfKnorm))
varianceGM21 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_gl$bandwidth, gmRegEstim, cdfKnorm))
varianceGM22 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_lo$bandwidth, gmRegEstim, cdfKnorm))
varianceGM23 <-sapply(x, function(x) estimateVariance(x, X, Y, h_lo$bandwidth, h_RT, gmRegEstim, cdfKnorm))
varianceGM31 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_gl$bandwidth, gmRegEstim, cdfKnorm))
varianceGM32 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_lo$bandwidth, gmRegEstim, cdfKnorm))
varianceGM33 <-sapply(x, function(x) estimateVariance(x, X, Y, h_RT, h_RT, gmRegEstim, cdfKnorm))


VarianceX <- sigma_squared(X)

# PLOT VARIANCE
plot(X, VarianceX, type = "l",ylim = c(0,5),ylab = "Variance de X", col = "black", main= paste("Comparaison de variances en utilisant la méthode Nadarya-Watson"))
lines(x, varianceNW11, type = "l", col = 2)
lines(x, varianceNW12, type = "l", col = 3)
lines(x, varianceNW13, type = "l", col = 4)
lines(x, varianceNW21, type = "l", col = 5)
lines(x, varianceNW22, type = "l", col = 6)
lines(x, varianceNW23, type = "l", col = 7)
lines(x, varianceNW31, type = "l", col = 8)
lines(x, varianceNW32, type = "l", col = 9)
lines(x, varianceNW33, type = "l", col = 10)
legend("topright", legend = c("True Variance","h1 =h_gl, h2 =h_gl","h1 =h_gl, h2 =h_lo","h1 =h_gl, h2 =h_RT","h1 =h_lo, h2 =h_gl","h1 =h_lo, h2 =h_lo","h1 =h_lo, h2 =h_RT","h1 =h_RT, h2 =h_gl","h1 =h_RT, h2 =h_lo","h1 =h_RT, h2 =h_RT"),
       col = c(1,2,3,4,5,6,7,8,9,10),
       lty = c(1, 1), cex = 0.8)


plot(X, VarianceX, type = "l", ylim=c(0,7), ylab = "Variance de X", col = "black", main= paste("Comparaison de variances en utilisant la méthode Gasser and Muller"))
lines(x, varianceGM11, type = "l", col = 2)
lines(x, varianceGM12, type = "l", col = 3)
lines(x, varianceGM13, type = "l", col = 4)
lines(x, varianceGM21, type = "l", col = 5)
lines(x, varianceGM22, type = "l", col = 6)
lines(x, varianceGM23, type = "l", col = 7)
lines(x, varianceGM31, type = "l", col = 8)
lines(x, varianceGM32, type = "l", col = 9)
lines(x, varianceGM33, type = "l", col = 10)
legend("topright", legend = c("True Variance","h1 =h_gl, h2 =h_gl","h1 =h_gl, h2 =h_lo","h1 =h_gl, h2 =h_RT","h1 =h_lo, h2 =h_gl","h1 =h_lo, h2 =h_lo","h1 =h_lo, h2 =h_RT","h1 =h_RT, h2 =h_gl","h1 =h_RT, h2 =h_lo","h1 =h_RT, h2 =h_RT"),
       col = c( 1,2,3,4,5,6,7,8,9,10),
       lty = c( 1, 1), cex = 0.8)

# variance for spline smoothing

estimateVarianceSMCV <- function(x, X, Y, h2, K){
  epsilonEstim = Y - smooth.spline(X, Y)$y
  fanAndYaoEstim <- sum(epsilonEstim^2 * K((x - X)/h2))/sum(K((x - X)/h2))
  return(fanAndYaoEstim)
}

estimateVarianceSM <- function(x, X, Y, h1, h2, K){
  
  epsilonEstim = Y - smooth.spline(X, Y, spar = h1)$y
  fanAndYaoEstim <- sum(epsilonEstim^2 * K((x - X)/h2))/sum(K((x - X)/h2))
  return(fanAndYaoEstim)
}

varianceSM_CV_h_gl <-sapply(x, function(x) estimateVarianceSMCV(x, X, Y, h_gl$bandwidth, Knorm))
varianceSM_CV_h_lo <-sapply(x, function(x) estimateVarianceSMCV(x, X, Y, h_lo$bandwidth, Knorm))
varianceSM_CV_h_RT <-sapply(x, function(x) estimateVarianceSMCV(x, X, Y, h_RT, Knorm))
varianceSM_hgl_hgl <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_gl$bandwidth, Knorm))
varianceSM_hgl_hlo <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_lo$bandwidth, Knorm))
varianceSM_hgl_hrt <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_gl$bandwidth, h_RT, Knorm))
varianceSM_hRT_hRT <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_RT, h_RT, Knorm))
varianceSM_hRT_hgl <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_RT, h_gl$bandwidth, Knorm))
varianceSM_hRT_hlo <-sapply(x, function(x) estimateVarianceSM(x, X, Y, h_RT, h_lo$bandwidth, Knorm))

VarianceX <- sigma_squared(X)

plot(X, VarianceX, type = "l",ylim = c(0,3.5) ,ylab = "Variance de X" , col = "black", main="Comparaison de variance en utilisant la méthode Spline Smoothing")
lines(x, varianceSM_CV_h_gl, col=2)
lines(x, varianceSM_CV_h_lo, col=3)
lines(x, varianceSM_CV_h_RT, col=4)
lines(x, varianceSM_hgl_hgl, col=5)
lines(x, varianceSM_hgl_hlo, col=6)
lines(x, varianceSM_hgl_hrt, col=7)
lines(x, varianceSM_hRT_hRT, col=8)
lines(x, varianceSM_hRT_hgl, col=9)
lines(x, varianceSM_hRT_hlo, col=10)
legend("topright", legend=c("true Variance","h1 = CV, h2 = h_gl","h1 = CV, h2 = h_lo","h1 = CV, h2 = h_RT","h1 = h_gl, h2 = h_gl","h1 = h_gl, h2 = h_lo","h1 = h_gl, h2 = h_RT","h1 = h_RT, h2 = h_RT","h1 = h_RT, h2 = h_gl","h1 = h_RT, h2 = h_lo"),
       col = c(1,2,3,4,5,6,7,8,9,10),lty=c(1,1,1,1,1,1,1,1,1,1), cex=0.8)

################################################
#   MONTE CARLO SIMULATION 
################################################
# we have to do the simulation for k = 25,50,100,500,1000
# monte carlo estimation
# number of replications

n <- 100 # number of observations
x <- seq(from = 0, to = 1, length = 100) #a sequence of points
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05)) #true regression model
epsilon <- rnorm(n, mean = 0, sd = 1)

GM_simulation <- function(x, sigma_squared, m, K, iter, h1, h2, regEstim){
  
  GMvarianceEstim <- matrix(NA, nrow = iter, ncol = length(x))
  GMvariance <- matrix(NA, nrow = iter, ncol = length(x))
  
  for (k in 1:iter) {
    set.seed(k)
    epsilon <- rnorm(n, mean = 0, sd = 1)
    set.seed(k+1)
    Y <- m(X) + 0.5 * sqrt(sigma_squared(X)) * epsilon
    for (i in 1:length(x)) {
      GMvarianceEstim[k,i] <- estimateVariance (x[i], X, Y, h1, h2, regEstim, K)
      GMvariance[k,i] <- sigma_squared(x[i])
    }
  }
  SE = matrix(NA, nrow = iter, ncol = length(x))
  for (i in 1:length(x)) {
    SE[,i] = (GMvarianceEstim[,i]-GMvariance[,i])^2
  }  

  meanGMvarianceEstim <- colMeans(GMvarianceEstim)
  meanGMvariance <- colMeans(GMvariance)
  biais = meanGMvarianceEstim - meanGMvariance
  variance<- colMeans(GMvarianceEstim^2)-meanGMvarianceEstim^2
  
  MSE = colMeans(SE)
  MSSE = mean(MSE) # MSSE as an estimation of MISE in simulation
  my_list <- list("biais" = biais,"IMSE" =  MSSE,"variance" = variance)
  return(my_list)
}


#GM biais 
GMresult_25 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 25, h_RT, h_RT, gmRegEstim)
GMresult_50 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 50, h_RT, h_RT, gmRegEstim)
GMresult_100 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 100, h_RT, h_RT, gmRegEstim)
GMresult_500 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 500, h_RT, h_RT, gmRegEstim)
GMresult_1000 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 1000, h_RT, h_RT, gmRegEstim)

cat("IMSE with Gasser-Müller : ",GMresult_1000$IMSE)


plot_stat <- function(stat, method,stat25, stat50, stat100, stat500, stat1000){
  plot(X, stat25, type = "l", col = "black", main=paste("", stat," de l'estimateur ", method," de la variance ", sep=" "))
  lines(x, stat50, col=2)
  lines(x, stat100, col=3)
  lines(x, stat500, col=4)
  lines(x, stat1000, col=5)
  legend("topright", legend=c("25 iterations","50 iterations","100 iterations","500 iterations","1000 iterations"),
         col = c(1,2,3,4,5),lty=c(1,1,1,1,1), cex=0.8)
}

plot_stat("Biais","GM",GMresult_25$biais, GMresult_50$biais, GMresult_100$biais, GMresult_500$biais, GMresult_1000$biais)
plot_stat("Variance","GM",GMresult_25$variance, GMresult_50$variance, GMresult_100$variance, GMresult_500$variance, GMresult_1000$variance)

#NW biais 
NWresult_25 <- GM_simulation(x, sigma_squared, m, Knorm, 25, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim)
NWresult_50 <- GM_simulation(x, sigma_squared, m, Knorm, 50, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim)
NWresult_100 <- GM_simulation(x, sigma_squared, m, Knorm, 100, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim)
NWresult_500 <- GM_simulation(x, sigma_squared, m, Knorm, 500, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim)
NWresult_1000 <- GM_simulation(x, sigma_squared, m, Knorm, 1000, h_gl$bandwidth, h_lo$bandwidth, nwRegEstim)

NWresult_25 <- GM_simulation(x, sigma_squared, m, Knorm, 25, h_RT, h_RT, nwRegEstim)
NWresult_50 <- GM_simulation(x, sigma_squared, m, Knorm, 50, h_RT, h_RT, nwRegEstim)
NWresult_100 <- GM_simulation(x, sigma_squared, m, Knorm, 100, h_RT, h_RT, nwRegEstim)
NWresult_500 <- GM_simulation(x, sigma_squared, m, Knorm, 500, h_RT, h_RT, nwRegEstim)
NWresult_1000 <- GM_simulation(x, sigma_squared, m, Knorm, 1000, h_RT, h_RT, nwRegEstim)


plot_stat("Biais","NW",NWresult_25$biais, NWresult_50$biais, NWresult_100$biais, NWresult_500$biais, NWresult_1000$biais)
plot_stat("Variance","NW",NWresult_25$variance, NWresult_50$variance, NWresult_100$variance, NWresult_500$variance, NWresult_1000$variance)

cat("IMSE with GNadarya-Watson : ",NWresult_1000$IMSE)


# Spline smoothing simulation avec CV
n <- 100 # number of observations
x <- seq(from = 0, to = 1, length = 100) #a sequence of points
m <- function(x) 3 * (x)^(1/4) * cos(1.2/(x + 0.05)) #true regression model
epsilon <- rnorm(n, mean = 0, sd = 1)
SM_simulation <- function(x, sigma_squared, m, K, iter, h2){
  
  SMvarianceEstim <- matrix(NA, nrow = iter, ncol = length(x))
  SMvariance <- matrix(NA, nrow = iter, ncol = length(x))
  for (k in 1:iter) {
    set.seed(k)
    epsilon <- rnorm(n, mean = 0, sd = 1)
    set.seed(k+1)
    Y <- m(X) + 0.5 * sqrt(sigma_squared(X)) * epsilon
    for (i in 1:length(x)) {
      SMvarianceEstim[k,i] <- estimateVarianceSMCV (x[i], X, Y, h2, K)
      SMvariance[k,i] <- sigma_squared(x[i])
    }
  }
  
  SE = matrix(NA, nrow = iter, ncol = length(x))
  for (i in 1:length(x)) {
    SE[,i] = (SMvarianceEstim[,i]-SMvariance[,i])^2
  }  
  
  meanSMvarianceEstim <- colMeans(SMvarianceEstim)
  meanSMvariance <- colMeans(SMvariance)
  biais = meanSMvarianceEstim - meanSMvariance
  variance<- colMeans(SMvarianceEstim^2)-meanSMvarianceEstim^2
  MSE = colMeans(SE)
  MSSE = mean(MSE) # MSSE as an estimation of MISE in simulation
  my_list <- list("biais" = biais,"IMSE" =  MSSE,"variance" = variance)
  return(my_list)
}

SMresult_25 <- SM_simulation(x, sigma_squared, m, Knorm, 25, h_gl$bandwidth)
SMresult_50 <- SM_simulation(x, sigma_squared, m, Knorm, 50, h_gl$bandwidth)
SMresult_100 <- SM_simulation(x, sigma_squared, m, Knorm, 100, h_gl$bandwidth)
SMresult_500 <- SM_simulation(x, sigma_squared, m, Knorm, 500, h_gl$bandwidth)
SMresult_1000 <- SM_simulation(x, sigma_squared, m, Knorm, 1000, h_gl$bandwidth)

plot_stat("Biais","SM",SMresult_25$biais, SMresult_50$biais, SMresult_100$biais, SMresult_500$biais, SMresult_1000$biais)
plot_stat("Variance","SM",SMresult_25$variance, SMresult_50$variance, SMresult_100$variance, SMresult_500$variance, SMresult_1000$variance)


cat("IMSE with Spline smoothing : ",SMresult_1000$IMSE)
