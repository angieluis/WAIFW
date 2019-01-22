
########################################################################
# Functions to model and estimate parameters with 
# transmission rates varying by sex
########################################################################

library(zoo)
library(deSolve)

WAIFW.model <- function(t, x, f.spline, m.spline, K.spline, params){
  # parameter estimates from Luis et al. 2015
  mu <- 0.085           # disease-induced mortality rate
  b <- 0.3154089        # max birth rate
  d <- 3.655774e-05     # min death rate
  a <- 0.614247         # how density-dependence is partitioned across births and deaths

  # transmission rates
  beta.ff=params["beta.ff"] # from female to female
  beta.mf=params["beta.mf"] # from female to male
  beta.fm=params["beta.fm"] # from male to female
  beta.mm=params["beta.mm"] # from male to male
  
  phi=params["phi"]         # infected immigration
  K=predict(K.spline,t)$y
  Nf=predict(f.spline,t)$y
  Nm=predict(m.spline,t)$y
  If=x["If"]
  Im=x["Im"]
  Sf=Nf-If
  Sm=Nm-Im
  N=Nf+Nm
  
  
  dIf=Sf*(beta.ff*If+beta.fm*Im)-If*(mu+d+(1-a)*(b-d)*N/K)+phi
  
  dIm=Sm*(beta.mf*If+beta.mm*Im)-Im*(mu+d+(1-a)*(b-d)*N/K)+phi
  
  list(c(dIf,dIm))
}



null.model <- function(t, x, f.spline, m.spline, K.spline, params){
  # parameter estimates from Luis et al. 2015
  mu <- 0.085           # disease-induced mortality rate
  b <- 0.3154089        # max birth rate
  d <- 3.655774e-05     # min death rate
  a <- 0.614247         # how density-dependence is partitioned across births and deaths
  
  # tranmsision rates the same
  beta.ff=beta.mf=beta.fm=beta.mm=params["beta"]
  
  phi=params["phi"]         # infected immigration
  K=predict(K.spline,t)$y
  Nf=predict(f.spline,t)$y
  Nm=predict(m.spline,t)$y
  If=x["If"]
  Im=x["Im"]
  Sf=Nf-If
  Sm=Nm-Im
  N=Nf+Nm
  
  
  dIf=Sf*(beta.ff*If+beta.fm*Im)-If*(mu+d+(1-a)*(b-d)*N/K)+phi
  
  dIm=Sm*(beta.mf*If+beta.mm*Im)-Im*(mu+d+(1-a)*(b-d)*N/K)+phi
  
  list(c(dIf,dIm))
  
}


trajmatchSexfun <- function(model, #WAIFW.model or null.model
                         data,  # data frame like Cascade11
                         init){    # initial values for parameters
  MNA.f <- data$MNAf
  MNA.m <- data$MNAm
  MNI.f <- data$MNIf
  MNI.m <- data$MNIm
  MNA <- data$MNA
  
  N.int <- na.approx(MNA,1:length(MNA))
  f.int <- na.approx(MNA.f,1:length(MNA.f))
  m.int <- na.approx(MNA.m,1:length(MNA.m))
  f.spline <- smooth.spline(f.int, df=length(f.int))
  m.spline <- smooth.spline(m.int, df=length(m.int))
  K.spline <- smooth.spline(N.int[c(4:length(N.int),length(N.int),length(N.int), length(N.int))],df=length(N.int)) 
  

  ind <- min(which(MNI.f+MNI.m>0))
  
  # trajectory matching
  trajmatch <- function(p, ...){
    params=p
    
    out <- matrix(NA,nrow=length(MNA.f[ind:length(MNA.f)]),ncol=3)
    xstart <- c( If=MNI.f[ind],Im=MNI.m[ind])
    out <- lsoda(xstart, ind:length(MNA.f), model, params, f.spline=f.spline, m.spline=m.spline, K.spline=K.spline)  
    Lf <- sum(dpois(floor(MNI.f[ind:length(MNA.f)]),out[,2],log=TRUE),na.rm=TRUE)
    Lm <- sum(dpois(floor(MNI.m[ind:length(MNA.f)]),out[,3],log=TRUE),na.rm=TRUE)
    return(-(Lf+Lm))
  }
  fit1 <- optim(init, trajmatch, f.spline=f.spline, m.spline=m.spline, K.spline=K.spline, method="L-BFGS-B",lower=rep(0,length(init)),upper=rep(1,length(init)),control=list(trace=2),hessian=TRUE)
  
  return(fit1)
}




##################################################################
# Model and plot fits for Montana sites
##################################################################
MT.site.names <- c("Cascade11","Cascade12","CascadeMITH","PetersonUpper","PetersonLower","Smith")
SW.site.names <- c("HesperusA","HesperusB")

MT.WAIFW.trajmatch <- list()
MT.null.trajmatch <- list()

for(i in 1:length(MT.site.names)){
  data=get(MT.site.names[i])
  
  MT.WAIFW.trajmatch[[i]] <- trajmatchSexfun(model=WAIFW.model, data=data, init=c(beta.ff=0.001,beta.mf=0.005,beta.fm=0.001,beta.mm=0.005,phi=0.02))

  MT.null.trajmatch[[i]] <- trajmatchSexfun(model=null.model, data=data, init=c(beta=0.001,phi=0.02))

  cat(i,"\n")
}



