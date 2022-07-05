library(MASS)
# parameters changing
n=200  # {200, 1000}
VIO.str =0.5 # {0.5, 2}
sim.round = 1 # seq(1,20)
nsim = 25 # the number of simulations

####################
case = "homo"
IV.str = 0.5
pi.value = IV.str*VIO.str
beta = 1
L = 100
alpha = c(rep(0, 5), rep(pi.value, 2), rep(0, L-7))
gamma = c(rep(IV.str, 7), rep(0, L-7))

px = 150
p=L+px # p stands for the number of total exogeneous variables
phi<-rep(0,px)
psi<-rep(0,px)
phi[1:10] <- 1/10*seq(1, 10)+0.5
psi[1:10] <- 1/10*seq(1, 10)+1

rho=0.5
A1gen <- function(rho, p){
  A1 = matrix(0, nrow=p, ncol=p)
  for(i in 1:p) for(j in 1:p) A1[i, j] = rho^(abs(i-j))
  return(A1)
}
Cov<-(A1gen(rho,p))

########### Part-2: Simulations ##########
CI.sear.m1 = CI.sear.m2 = matrix(NA, nrow=nsim, ncol=2)
CI.samp.m1 = CI.sear.m2 = matrix(NA, nrow=nsim, ncol=2)

for(i.sim in 1:nsim){
  set.seed(i.sim+(sim.round-1)*nsim)
  print(i.sim)
  W = mvrnorm(n, rep(0, p), Cov)
  Z = W[, 1:L]
  X = W[, (L+1):p]
  if(case=="homo"){
    # epsilonSigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)
    epsilonSigma = matrix(c(1.5, 0.75, 0.75, 1.5), 2, 2)
    epsilon = mvrnorm(n, rep(0, 2), epsilonSigma)
    epsilon1 = epsilon[,1]
    epsilon2 = epsilon[,2]
  }
  D = 0.5 + Z %*% gamma+ X%*% psi + epsilon1
  Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon2
  cat("DeLasso method... \n")
  # out1 = SearchingSampling(Y, D, Z, X, intercept=TRUE, method="DeLasso", Sampling=FALSE)
  out2 = SearchingSampling(Y, D, Z, X, intercept=TRUE, method="DeLasso", Sampling=TRUE)
  # CI.sear.m1[i.sim, ] = out1$ci
  CI.samp.m1[i.sim, ] = out2$ci
  cat("Fast.DeLasso method... \n")
  # out3 = SearchingSampling(Y, D, Z, X, intercept=TRUE, method="Fast.DeLasso", Sampling=FALSE)
  out4 = SearchingSampling(Y, D, Z, X, intercept=TRUE, method="Fast.DeLasso", Sampling=TRUE)
  # CI.sear.m2[i.sim, ] = out3$ci
  CI.samp.m2[i.sim, ] = out4$ci
}

rm(Y, D, Z, X)
filename <- paste("Highd-n",n,"-pz",L,"-Violation",VIO.str,"-SimRound",sim.round,".RData",sep="")
save.image(filename)
