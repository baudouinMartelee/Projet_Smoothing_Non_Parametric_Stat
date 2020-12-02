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
  fanAndYaoEstim <- sum((epsilonEstim^2) * K((x - X)/h2))/sum(K((x - X)/h2))
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
lines(x, varianceNW4, type = "l", col = 5)
lines(x, varianceNW5, type = "l", col = 2)
lines(x, varianceNW6, type = "l", col = 3)
lines(x, varianceNW7, type = "l", col = 4)


legend("topright", legend = c("True Variance", "variance NW, h1 = gl,h2 = lo","variance NW, h1 = gl,h2 = RT","variance NW, h1 = lo,h2 = gl","variance NW, h1 = lo,h2 = RT","variance NW, h1 = RT,h2 = gl","variance NW, h1 = RT,h2 = lo"),
       col = c( 1,2,3,4,5),
       lty = c( 1, 1), cex = 0.5)

plot(X, VarianceX, type = "l",ylim = c(0,30), col = "black", main= paste("Comparaison de Variance en utilisant la méthode Gasser and Muller"))
lines(x, varianceGM1, type = "l", col = 2)
lines(x, varianceGM2, type = "l", col = 3)
lines(x, varianceGM3, type = "l", col = 4)
lines(x, varianceGM32, type = "l", col = 5)
lines(x, varianceGM4, type = "l", col = 5)
lines(x, varianceGM5, type = "l", col = 5)
lines(x, varianceGM6, type = "l", col = 5)

legend("topright", legend = c("True Variance", "variance NW, h1 = gl,h2 = lo","variance NW, h1 = gl,h2 = RT","variance NW, h1 = lo,h2 = gl","variance NW, h1 = lo,h2 = RT","variance NW, h1 = RT,h2 = gl","variance NW, h1 = RT,h2 = lo"),
       col = c( 1,2,3,4,5),
       lty = c( 1, 1), cex = 0.5)

# variance for spline smoothing

estimateVarianceSMCV <- function(x, X, Y, h2, K){
  epsilonEstim = 2*(Y - smooth.spline(X, Y)$y)
  fanAndYaoEstim <- sum(epsilonEstim^2 * K((x - X)/h2))/sum(K((x - X)/h2))
  return(fanAndYaoEstim)
}

estimateVarianceSM <- function(x, X, Y, h1, h2, K){
  
  epsilonEstim = 2*(Y - smooth.spline(X, Y, spar = h1)$y)
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

plot(X, VarianceX, type = "l",ylim = c(0,5), col = "black", main="Variance using Spline Smoothing")
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
  #cat(class(meanGMvarianceEstim))
  #print.table(meanGMvarianceEstim)
  #print.table(meanGMvariance)
  biais = meanGMvarianceEstim - meanGMvariance
  variance<- colMeans(GMvarianceEstim^2)-meanGMvarianceEstim^2
  
  MSE = colMeans(SE)
  MSSE = mean(MSE) # MSSE as an estimation of MISE in simulation
  #plot(x, biais, type="l", col = 2)
  #abline(h=0)
  my_list <- list("biais" = biais,"IMSE" =  MSSE,"variance" = variance)
  return(my_list)
}


#GM biais 
GMresult_25 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 25, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim)
GMresult_50 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 50, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim)
GMresult_100 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 100, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim)
GMresult_500 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 500, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim)
GMresult_1000 <- GM_simulation(x, sigma_squared, m, cdfKnorm, 1000, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim)


GMresult_25$biais
GMresult_25$variance
cat("IMSE with Gasser-Müller : ",GMresult_1000$IMSE)

GMresult_1000$biais

#MSE <- GM_simulation(x, sigma_squared, m, cdfKnorm, 25, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim)[[2]]

#variance <- GM_simulation(x, sigma_squared, m, cdfKnorm, 25, h_gl$bandwidth, h_gl$bandwidth, gmRegEstim)[[3]]


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
NWresult_25 <- GM_simulation(x, sigma_squared, m, Knorm, 25, h_gl$bandwidth, h_gl$bandwidth, nwRegEstim)
NWresult_50 <- GM_simulation(x, sigma_squared, m, Knorm, 50, h_gl$bandwidth, h_gl$bandwidth, nwRegEstim)
NWresult_100 <- GM_simulation(x, sigma_squared, m, Knorm, 100, h_gl$bandwidth, h_gl$bandwidth, nwRegEstim)
NWresult_500 <- GM_simulation(x, sigma_squared, m, Knorm, 500, h_gl$bandwidth, h_gl$bandwidth, nwRegEstim)
NWresult_1000 <- GM_simulation(x, sigma_squared, m, Knorm, 1000, h_gl$bandwidth, h_gl$bandwidth, nwRegEstim)


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
  #cat(class(meanGMvarianceEstim))
  print.table(meanSMvarianceEstim)
  print.table(meanSMvariance)
  biais = meanSMvarianceEstim - meanSMvariance
  cat(biais)
  #variance1<-colMeans((NWnormEst-NW)^2)
  #plot(x,variance1 ,ylab="Variance",type="l",col="#C72C48")
  variance<- colMeans(SMvarianceEstim^2)-meanSMvarianceEstim^2
  MSE = colMeans(SE)
  MSSE = mean(MSE) # MSSE as an estimation of MISE in simulation
  
  
  #plot(x, biais, type="l", col = 2)
  #abline(h=0)
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
