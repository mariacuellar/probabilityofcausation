# # Estimating the probability of causation for the RCT described (and with data) here 
# # https://www.povertyactionlab.org/evaluation/cleaning-springs-kenya

 # rm(list=ls())


# -----------------------
# Packages
# -----------------------

# run this first
expit = function(x){exp(x)/(1+exp(x))}



# -----------------------
# 1) Simulate data ----
# -----------------------

fun.dat.sim <- function(n=100){

  # simulated data to train model
  n = n
  index = 1:n
  beta = 0.5

  # potential outcome when a=1
  y1 = rbinom(n = n, size = 1, prob = beta)

  # x's for correctly specified model
  x1 = rnorm(n, 0, 0.1)
  x2 = rnorm(n, 0, 0.1)
  x3 = rnorm(n, 0, 0.1)
  x4 = rnorm(n, 0, 0.1)
  x = as.matrix(cbind(x1, x2, x3, x4))

  # x's for misspecified model, from Kang and Schafer (2007)
  x1star = exp(x1/2)
  x2star = x2/(1 + exp(x1)) + 10
  x3star = (x1*x3/25 + 0.6)^3
  x4star = (x2 + x4 + 20)^2
  xstar = as.matrix(cbind(x1star, x2star, x3star, x4star))

  # Kang and Schafer vector
  ks = c(-1, 0.5, -0.25, -0.1)

  # propensity score
  pi.n = expit(x %*% ks)
  pi = as.data.frame(expit(x %*% ks))
  colnames(pi) = "pi"

  # treatment
  a = rbinom(n, 1, pi.n)

  # nuisance functions
  mu0 = beta/(1+exp(x %*% ks))
  mu1 = rep(beta, n)
  gamma = as.data.frame(expit(x %*% ks)) # same as pi
  colnames(gamma) = "gamma"

  df = cbind(index, x, xstar, pi, gamma, a, y1)
  head(df)

  # generate y0 conditional on combinations of values of a and y1
  dfy11 = df[which(df$y1==1),]
  dfy11$y0 = rbinom(n = nrow(dfy11) , size = 1, prob = (1-pi.n))
  dfy10 = df[which(df$y1==0),]
  dfy10$y0 = 0

  # add y0 to dataframe
  df.wy0 = as.data.frame(rbind(dfy11, dfy10))

  # apply consistency to get y
  df.wy0$y = ifelse(df.wy0$a==1, df.wy0$y1, df.wy0$y0)
  y = df.wy0$y

  # ordering data so it's as it was at the beginning
  dff = df.wy0[order(df.wy0$index),]
  head(dff)

  return(dff)
}







# -----------------------
# 2) Function to estimate PC ----
# -----------------------

# function to estimate the probability of causation
pcausation = function(y, a, x, xtest, nsplits=2, start.list=c(rep(0,ncol(x))), tracetf=FALSE, printres=TRUE, showprogress=TRUE,
                      sl.lib = c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger","SL.rpart")){
  
  # packages
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("rpart")
  
  # setting things up
  n <- nrow(x)
  avals <- names(table(a))
  n.avals <- length(avals)
  s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
  
  # progress bar
  if(showprogress==TRUE){ pb <- txtProgressBar(min=0, max=2*nsplits*n.avals, style=3) }
  
  muhat <- as.data.frame(matrix(NA,nrow=n,ncol=n.avals))
  colnames(muhat) <- paste("a",avals,sep="")
  pihat <- muhat
  muhat.xtest <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=n.avals))
  pihat.xtest <- muhat.xtest
  
  # estimating nuisance parameters: propensity score and outcome regressions
  pbcount <- 0
  for (i in 1:n.avals){
    if (i==1){ Sys.sleep(0.1); if(showprogress==TRUE){setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }}
    
    # starting cross validation
    for (vfold in 1:nsplits){
      
      # split into train and test sets
      train <- s!=vfold; test <- s==vfold
      if (nsplits==1){ train <- test }
      
      # estimate propensity score (pi)
      if (i != n.avals){
        
        pifit.xtest <- SuperLearner(as.numeric(a==avals[i])[train],as.data.frame(x[train,]),
                                    newX=xtest, SL.library=sl.lib, family=binomial)
        
        pihat.xtest <- as.numeric(pifit.xtest$SL.predict) # evaluated at pre-selected set of x's 
        
        Sys.sleep(0.1)
        if(showprogress==TRUE){ setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }
      }
      
    }
    
    # estimate outcome regression function (mu), one for each value of a
    mufit.xtest <- SuperLearner(y[a==avals[i] & train],
                                as.data.frame(x[a==avals[i] & train,]),
                                newX=xtest, SL.library=sl.lib)
    
    muhat.xtest[,i] <- mufit.xtest$SL.predict # evaluated at pre-selected set of x's
    
    Sys.sleep(0.1)
    if(showprogress==TRUE){ setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }
    
  }
  if (i == n.avals){ pihat[,i] <- 1 - apply(pihat,1,sum, na.rm=T) }
  
  # estimates of E{Y(0)} and E{Y(1)}
  names(muhat.xtest) <- c("a0", "a1")
  amat <- matrix(rep(a,n.avals),nrow=n,byrow=F)
  alevel <- matrix(rep(as.numeric(avals), rep(n,n.avals)),nrow=n,byrow=F)
  ymat <- matrix(rep(y,n.avals),nrow=n,byrow=F)
  
  # estimates of probability of causation - plugin estimation
  est.pi <- abs(1 - muhat.xtest[1]/muhat.xtest[2]); names(est.pi)<-"est.pi" # se and CI's are not of interest
  
  # estimates of probability of causation - influence function estimation
  include.x = names(as.data.frame(x))
  include.b = gsub("x", "b", include.x)
  include.xb = diag(sapply( include.b, FUN=function(x, y) paste(x, include.x, sep="*") ))
  formula.x = paste("y ~ ", paste(include.x, collapse="+"), sep="")
  formula.pi = paste("a ~ ", paste(include.x, collapse="+"), sep="")
  formula.nls.x = paste("y ~ expit(", paste(include.xb, collapse="+"), ")", sep = "")
  start.values = setNames(as.list(start.list), c(include.b))
  ystar = (1/muhat.xtest[2])*((muhat.xtest[1]/muhat.xtest[2])*(1/pihat.xtest)*a*(y-muhat.xtest[2]) -
                                (1/(1-pihat.xtest))*(1-a)*(y-muhat.xtest[1])) + (1 - muhat.xtest[1]/muhat.xtest[2])
  ifvals <- ystar; names(ifvals) = "ystar"
  
  # fit nls model - gives NA if the NLS doesn't converge
  mod = tryCatch({
    nls( as.formula(formula.nls.x), start=start.values, data=as.data.frame(cbind(x,y,a)), nls.control(maxiter = 500), algorithm = "port",
        lower=c(rep(-5, length(include.x))), upper=c(rep(5, length(include.x))), trace = tracetf )
  }, warning = function(warning_condition) {
    mod <- NA
  }, error = function(error_condition) {
    mod <- NA
  })
  
  # estimate covariance matrix for betahat
  preds <- try( predict(mod) , silent = TRUE )
  xmat <- try( as.matrix(x) , silent = TRUE ); try( xtestmat <- as.matrix(xtest) , silent = TRUE )
  wts <- try( mod$weights , silent = TRUE ); try( if (is.null(wts)){ wts <- 1 } , silent = TRUE )
  bread <- try( solve( (t(xmat *(preds*(1-preds)*wts)) %*% xmat)/n ) , silent = TRUE )
  meat <- try( (t(xmat *( ((y-preds) * wts)^2)) %*% xmat)/n , silent = TRUE )
  vcov <- try( bread %*% meat %*% bread / n , silent = TRUE )
  coefs <- try( coef(mod) , silent = TRUE )
  res.betas <- tryCatch({ data.frame(Estimate=coefs, Robust.SE=sqrt(diag(vcov)), 
                          z.val=coefs/sqrt(diag(vcov)),	p.val= round(2*(1-pnorm(abs(coefs/sqrt(diag(vcov))))),3) , 
                          ci.ll= coefs-1.96*sqrt(diag(vcov)) , ci.ul=coefs+1.96*sqrt(diag(vcov)) )
                          }, warning = function(warning_condition){ res.betas })
  
  # get predicted value / CI for probability at specific x0
  x0 <- xtestmat
  est.if <- try( expit(x0 %*% coefs) , silent = TRUE )
  se.if <- try( as.data.frame(apply(  x0, 1, function(M){(t(M) %*% vcov %*% M)}  )) , silent = TRUE ); names(se.if)="se.if"
  ci.ll.if = try( as.data.frame(apply(  x0, 1, function(M){expit(M %*% coefs - 1.96*(t(M) %*% vcov %*% M))}  )) , silent = TRUE ); names(ci.ll.if)="ci.ll.if"
  ci.ul.if = try( as.data.frame(apply(  x0, 1, function(M){expit(M %*% coefs + 1.96*(t(M) %*% vcov %*% M))}  )) , silent = TRUE ); names(ci.ul.if)="ci.ul.if"
  
  # get SEs/CIs/pvals for betas
  coefs <- try( coef(mod), silent=TRUE)
  res.coefs <- try( data.frame(Estimate=coefs, Robust.SE=sqrt(diag(vcov)), 
                        z.val=coefs/sqrt(diag(vcov)),	p.val= round(2*(1-pnorm(abs(coefs/sqrt(diag(vcov))))),3) , 
                        ci.ll= coefs-1.96*sqrt(diag(vcov)) , ci.ul=coefs+1.96*sqrt(diag(vcov)) ), silent=TRUE)
  
  # table of values for influence function estimator and plugin estimator
  res.if <- tryCatch({
     data.frame( est.pc = est.if, se = se.if, ci.ll = ci.ll.if, ci.ul = ci.ul.if )
  }, warning = function(warning_condition) {
    empty = rep(NA,nrow(dat.eval))
    data.frame( est.pc = empty, se = empty, ci.ll = empty, ci.ul = empty )
  }, error = function(error_condition) {
    empty = rep(NA,nrow(dat.eval))
    data.frame( est.pc = empty, se = empty, ci.ll = empty, ci.ul = empty )
  })
  
  res.pi <- data.frame( est.pc = est.pi )
  
  Sys.sleep(0.1)
  if(showprogress==TRUE){ setTxtProgressBar(pb,pbcount) ; close(pb) }
  
  nuis <- as.data.frame(cbind(pihat.xtest,muhat.xtest))
  colnames(nuis) <- c("pihat", "mu0hat", "mu1hat")
  
  if(printres==TRUE){ print(res.if) }
  return(invisible(list(res.if=res.if, res.pi=res.pi, nuis=nuis, ifvals=ifvals, res.coefs=res.coefs)))
  
}

# simulated data - this would be the researcher's data
dat <- fun.dat.sim(n=500)

# evaluation data with multiple observations - this is for simulations, to evaluate at same X matrix at every iteration
dat.eval <- fun.dat.sim(n=5)

# testing function - using X's
out <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1,x2=dat$x2,x3=dat$x3,x4=dat$x4)), 
                  xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)), 
                  nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE, showprogress = FALSE,
                  sl.lib = c("SL.ranger"))
out

# testing function - using Xstar's
out1 <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1star,x2=dat$x2star,x3=dat$x3star,x4=dat$x4star)), 
                  xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)), 
                  nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE,
                  sl.lib = c("SL.ranger"))
out1





# -----------------------
# 3) Simulation to get RMSE, bias, and coverage ----
# -----------------------
# NOTE THAT THE FOR LOOP HAS TO BE RUN SEPARATELY FOR PARAMETRIC AND NONPARAMETRIC ESTIMATIONS

# number of simulations ( enumerated by 1,...,j,...,J )
n.sims = 20

# three sample sizes
sample.sizes.for.sim = c(200,1000,5000)

# use cross-validation?
nsplits <- 1 #nope

# how to estimate the nuisance functions (parametric vs. nonparametric)?
sl.lib <- c("SL.glm") # nonparametric
# sl.lib <- c("SL.glm") # parametric

# data for evaluation - will be the same for all estimators
dat.eval <- fun.dat.sim(5); xtest <- dat.eval

# where to store values
sim.df <- as.data.frame(matrix(NA,nrow=n.sims,ncol=24))
names(sim.df) <- c(
                    "plugin.correct.bias.500", "plugin.correct.bias.5000", "plugin.correct.bias.10000",
                    "plugin.correct.RMSE.500", "plugin.correct.RMSE.5000", "plugin.correct.RMSE.10000",
                    "ifb.correct.bias.500", "ifb.correct.bias.5000", "ifb.correct.bias.10000",
                    "ifb.correct.RMSE.500", "ifb.correct.RMSE.5000", "ifb.correct.RMSE.10000",
                    "plugin.misspecified.bias.500", "plugin.misspecified.bias.5000", "plugin.misspecified.bias.10000",
                    "plugin.misspecified.RMSE.500", "plugin.misspecified.RMSE.5000", "plugin.misspecified.RMSE.10000",
                    "ifb.misspecified.bias.500", "ifb.misspecified.bias.5000", "ifb.misspecified.bias.10000",
                    "ifb.misspecified.RMSE.500", "ifb.misspecified.RMSE.5000", "ifb.misspecified.RMSE.10000"
                   )

sim.df.cov <- as.data.frame(matrix(NA, nrow=n.sims, ncol=6))
names(sim.df.cov) <- c("ifb.correct.cov.500", "ifb.correct.cov.5000", "ifb.correct.cov.10000", 
                       "ifb.misspecified.cov.500", "ifb.misspecified.cov.5000", "ifb.misspecified.cov.10000" )

# simulation for loop
pb <- txtProgressBar(min=0, max=6*n.sims-1, style=3)
pbcount <- 0
Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount)

# j=i=1
for(j in 1:n.sims){
  
  # first sample size
  dat <- fun.dat.sim(sample.sizes.for.sim[1])

  # 1/6 correct.500 (plugin.correct.bias.500, plugin.correct.RMSE.500, ifb.correct.bias.500, ifb.correct.RMSE.500)
  pc.df <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=4))
  names(pc.df) <- c("PC.plugin", "PC.ifb", "PC.ci.ll", "PC.ci.ul")
  
  for(i in 1:nrow(dat.eval)){
    out <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1,x2=dat$x2,x3=dat$x3,x4=dat$x4)), 
                      xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)), 
                      nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE, showprogress = FALSE,
                      sl.lib = c("SL.ranger"))
    
    pc.df$PC.ifb <- out$res.if[,1]
    pc.df$PC.plugin <- out$res.pi[,1]
    pc.df$PC.ci.ll <- out$res.if[,3]
    pc.df$PC.ci.ul <- out$res.if[,4]
    
    }
  
  # mean bias
  sim.df$plugin.correct.bias.500[j] <- mean((pc.df$PC.plugin-dat.eval$gamma )^2, na.rm=TRUE)
  sim.df$ifb.correct.bias.500[j] <- mean((pc.df$PC.ifb-dat.eval$gamma)^2, na.rm=TRUE)
  # single RMSE
  sim.df$plugin.correct.RMSE.500[j] <- mean(abs(pc.df$PC.plugin-dat.eval$gamma), na.rm=TRUE)
  sim.df$ifb.correct.RMSE.500[j] <- mean(abs(pc.df$PC.ifb-dat.eval$gamma), na.rm=TRUE)
  # coverage
  sim.df.cov$ifb.correct.cov.500[j] <- length(which( dat.eval$gamma >= pc.df$PC.ci.ll & dat.eval$gamma <= pc.df$PC.ci.ul) ) / nrow(dat.eval)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
  


  # 2/6 misspecified.500 (plugin.misspecified.bias.500, plugin.misspecified.RMSE.500, ifb.misspecified.bias.500, ifb.misspecified.RMSE.500)
  pc.df <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=2))
  names(pc.df) <- c("PC.plugin", "PC.ifb")
  for(i in 1:nrow(dat.eval)){
    out <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1star,x2=dat$x2star,x3=dat$x3star,x4=dat$x4star)), 
                      xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)),
                      nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE, showprogress = FALSE,
                      sl.lib = c("SL.ranger"))
    
    pc.df$PC.ifb <- out$res.if[,1]
    pc.df$PC.plugin <- out$res.pi[,1]
    pc.df$PC.ci.ll <- out$res.if[,3]
    pc.df$PC.ci.ul <- out$res.if[,4]
    
    }
  # mean bias
  sim.df$plugin.misspecified.bias.500[j] <- mean((pc.df$PC.plugin-dat.eval$gamma)^2, na.rm=TRUE)
  sim.df$ifb.misspecified.bias.500[j] <- mean((pc.df$PC.ifb-dat.eval$gamma)^2, na.rm=TRUE)
  # mean RMSE
  sim.df$plugin.misspecified.RMSE.500[j] <- mean(abs(pc.df$PC.plugin-dat.eval$gamma), na.rm=TRUE)
  sim.df$ifb.misspecified.RMSE.500[j] <- mean(abs(pc.df$PC.ifb-dat.eval$gamma), na.rm=TRUE)
  # coverage
  sim.df.cov$ifb.misspecified.cov.500[j] <- length(which( dat.eval$gamma >= pc.df$PC.ci.ll & dat.eval$gamma <= pc.df$PC.ci.ul) ) / nrow(dat.eval)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
  


  # second sample size
  dat <- fun.dat.sim(sample.sizes.for.sim[2])

  # 3/6 correct.5000 (plugin.correct.bias.5000, plugin.correct.RMSE.5000, ifb.correct.bias.5000, ifb.correct.RMSE.5000)
  pc.df <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=2))
  names(pc.df) <- c("PC.plugin", "PC.ifb")
  for(i in 1:nrow(dat.eval)){
    out <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1,x2=dat$x2,x3=dat$x3,x4=dat$x4)), 
                      xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)), 
                      nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE, showprogress = FALSE,
                      sl.lib = c("SL.ranger"))

    pc.df$PC.ifb <- out$res.if[,1]
    pc.df$PC.plugin <- out$res.pi[,1]
    pc.df$PC.ci.ll <- out$res.if[,3]
    pc.df$PC.ci.ul <- out$res.if[,4]

    }
  # mean bias
  sim.df$plugin.correct.bias.5000[j] <- mean((pc.df$PC.plugin-dat.eval$gamma)^2, na.rm=TRUE)
  sim.df$ifb.correct.bias.5000[j] <- mean((pc.df$PC.ifb-dat.eval$gamma)^2, na.rm=TRUE)
  # mean RMSE
  sim.df$plugin.correct.RMSE.5000[j] <- mean(abs(pc.df$PC.plugin-dat.eval$gamma), na.rm=TRUE)
  sim.df$ifb.correct.RMSE.5000[j] <- mean(abs(pc.df$PC.ifb-dat.eval$gamma), na.rm=TRUE)
  # coverage
  sim.df.cov$ifb.correct.cov.5000[j] <- length(which( dat.eval$gamma >= pc.df$PC.ci.ll & dat.eval$gamma <= pc.df$PC.ci.ul) ) / nrow(dat.eval)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
  

  # 4/6 misspecified.5000 (plugin.misspecified.bias.5000, plugin.misspecified.RMSE.5000, ifb.misspecified.bias.5000, ifb.misspecified.RMSE.5000)
  pc.df <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=2))
  names(pc.df) <- c("PC.plugin", "PC.ifb")
  for(i in 1:nrow(dat.eval)){
    out <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1star,x2=dat$x2star,x3=dat$x3star,x4=dat$x4star)), 
                      xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)),
                      nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE, showprogress = FALSE,
                      sl.lib = c("SL.ranger"))

    pc.df$PC.ifb <- out$res.if[,1]
    pc.df$PC.plugin <- out$res.pi[,1]
    pc.df$PC.ci.ll <- out$res.if[,3]
    pc.df$PC.ci.ul <- out$res.if[,4]
    
    }
  # mean bias
  sim.df$plugin.misspecified.bias.5000[j] <- mean((pc.df$PC.plugin-dat.eval$gamma)^2, na.rm=TRUE)
  sim.df$ifb.misspecified.bias.5000[j] <- mean((pc.df$PC.ifb-dat.eval$gamma)^2, na.rm=TRUE)
  # mean RMSE
  sim.df$plugin.misspecified.RMSE.5000[j] <- mean(abs(pc.df$PC.plugin-dat.eval$gamma), na.rm=TRUE)
  sim.df$ifb.misspecified.RMSE.5000[j] <- mean(abs(pc.df$PC.ifb-dat.eval$gamma), na.rm=TRUE)
  # coverage
  sim.df.cov$ifb.misspecified.cov.5000[j] <- length(which( dat.eval$gamma >= pc.df$PC.ci.ll & dat.eval$gamma <= pc.df$PC.ci.ul) ) / nrow(dat.eval)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
  
  


  # third sample size
  dat <- fun.dat.sim(sample.sizes.for.sim[3])
  
  # 5/6 correct.10000 (plugin.correct.bias.10000, plugin.correct.RMSE.10000, ifb.correct.bias.10000, ifb.correct.RMSE.10000)
  pc.df <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=2))
  names(pc.df) <- c("PC.plugin", "PC.ifb")
  for(i in 1:nrow(dat.eval)){
    out <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1,x2=dat$x2,x3=dat$x3,x4=dat$x4)), 
                      xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)), 
                      nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE, showprogress = FALSE,
                      sl.lib = c("SL.ranger"))
    
    pc.df$PC.ifb <- out$res.if[,1]
    pc.df$PC.plugin <- out$res.pi[,1]
    pc.df$PC.ci.ll <- out$res.if[,3]
    pc.df$PC.ci.ul <- out$res.if[,4]

    }
  # mean bias
  sim.df$plugin.correct.bias.10000[j] <- mean((pc.df$PC.plugin-dat.eval$gamma)^2, na.rm=TRUE)
  sim.df$ifb.correct.bias.10000[j] <- mean((pc.df$PC.ifb-dat.eval$gamma)^2, na.rm=TRUE)
  # mean RMSE
  sim.df$plugin.correct.RMSE.10000[j] <- mean(abs(pc.df$PC.plugin-dat.eval$gamma), na.rm=TRUE)
  sim.df$ifb.correct.RMSE.10000[j] <- mean(abs(pc.df$PC.ifb-dat.eval$gamma), na.rm=TRUE)
  # coverage
  sim.df.cov$ifb.correct.cov.10000[j] <- length(which( dat.eval$gamma >= pc.df$PC.ci.ll & dat.eval$gamma <= pc.df$PC.ci.ul) ) / nrow(dat.eval)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
  

  # 6/6 misspecified.10000 (plugin.misspecified.bias.10000, plugin.misspecified.RMSE.10000, ifb.misspecified.bias.10000, ifb.misspecified.RMSE.10000)
  pc.df <- as.data.frame(matrix(NA,nrow=nrow(xtest),ncol=2))
  names(pc.df) <- c("PC.plugin", "PC.ifb")
  for(i in 1:nrow(dat.eval)){
    out <- pcausation(y=dat$y, a=dat$a, x=as.data.frame(cbind(x1=dat$x1star,x2=dat$x2star,x3=dat$x3star,x4=dat$x4star)), 
                      xtest=as.data.frame(cbind(x1=dat.eval$x1,x2=dat.eval$x2,x3=dat.eval$x3,x4=dat.eval$x4)),
                      nsplits=1, start.list=c(rep(0,4)), tracetf=FALSE, printres=FALSE, showprogress = FALSE,
                      sl.lib = c("SL.ranger"))
    
    pc.df$PC.ifb <- out$res.if[,1]
    pc.df$PC.plugin <- out$res.pi[,1]
    pc.df$PC.ci.ll <- out$res.if[,3]
    pc.df$PC.ci.ul <- out$res.if[,4]

    }
  
  # mean bias
  sim.df$plugin.misspecified.bias.10000[j] <- mean((pc.df$PC.plugin-dat.eval$gamma)^2, na.rm=TRUE)
  sim.df$ifb.misspecified.bias.10000[j] <- mean((pc.df$PC.ifb-dat.eval$gamma)^2, na.rm=TRUE)
  # mean RMSE
  sim.df$plugin.misspecified.RMSE.10000[j] <- mean(abs(pc.df$PC.plugin-dat.eval$gamma), na.rm=TRUE)
  sim.df$ifb.misspecified.RMSE.10000[j] <- mean(abs(pc.df$PC.ifb-dat.eval$gamma), na.rm=TRUE)
  # coverage
  sim.df.cov$ifb.misspecified.cov.10000[j] <- length(which( dat.eval$gamma >= pc.df$PC.ci.ll & dat.eval$gamma <= pc.df$PC.ci.ul) ) / nrow(dat.eval)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

}
Sys.sleep(0.1)
close(pb)

head(sim.df)

# NOTE THAT THE FOR LOOP HAS TO BE RUN SEPARATELY FOR PARAMETRIC AND NONPARAMETRIC ESTIMATIONS

# save(sim.df, file="sim.df.np-2018-08-25.rda")
# save(sim.df, file="sim.df.p-2018-08-25.rda")


# take mean across all j simulations - RMSE and bias
sim.df.m <- as.data.frame(t(colMeans(sim.df, dims = 1, na.rm=TRUE)))
colnames(sim.df.m) = names(sim.df.m)
sim.df.m = round(sim.df.m, 2)
sim.df.m

# sim.df.m.np <- format(round(sim.df.m, digits=2), nsmall = 2) 
# sim.df.m.p <- format(round(sim.df.m, digits=2), nsmall = 2) 


# take mean across all j simulations - coverage
sim.df.cov.m <- as.data.frame(t(colMeans(sim.df.cov, dims = 1, na.rm=TRUE)))
colnames(sim.df.cov.m) = names(sim.df.cov.m)
sim.df.cov.m = round(sim.df.cov.m, 2)
sim.df.cov.m

# sim.df.cov.m.np <- format(round(sim.df.cov.m, digits=2), nsmall = 2) 
# sim.df.cov.m.p <- format(round(sim.df.cov.m, digits=2), nsmall = 2) 




# -----------------------
# 4) Preparing for publication ----
# -----------------------

# make table for publication - RMSE and bias
publish.df.biasRMSE <- as.data.frame(matrix(NA,nrow=6,ncol=6))
names(publish.df.biasRMSE) <- c("Sample size", "Method", "Cor P", "Mis P", "Cor NP", "Mis NP")
publish.df.biasRMSE$`Sample size` <- c(500, "", 5000, "", 10000, "")
publish.df.biasRMSE$`Method` <- c("Plug-in", "Proposed", "Plug-in", "Proposed", "Plug-in", "Proposed")
publish.df.biasRMSE$`Cor P` <- c( paste( sim.df.m.p$plugin.correct.bias.500, " (", sim.df.m.p$plugin.correct.RMSE.500,")", sep=""), 
                          paste( sim.df.m.p$ifb.correct.bias.500, " (", sim.df.m.p$ifb.correct.RMSE.500,")", sep=""),
                          paste( sim.df.m.p$plugin.correct.bias.5000, " (", sim.df.m.p$plugin.correct.RMSE.5000,")", sep=""),
                          paste( sim.df.m.p$ifb.correct.bias.5000, " (", sim.df.m.p$ifb.correct.RMSE.5000,")", sep=""),
                          paste( sim.df.m.p$plugin.correct.bias.10000, " (", sim.df.m.p$plugin.correct.RMSE.10000,")", sep=""),
                          paste( sim.df.m.p$ifb.correct.bias.10000, " (", sim.df.m.p$ifb.correct.RMSE.10000,")", sep="")
                        )
publish.df.biasRMSE$`Mis P` <- c( paste( sim.df.m.p$plugin.misspecified.bias.500, " (", sim.df.m.p$plugin.misspecified.RMSE.500,")", sep=""),
                          paste( sim.df.m.p$ifb.misspecified.bias.500, " (", sim.df.m.p$ifb.misspecified.RMSE.500,")", sep=""),
                          paste( sim.df.m.p$plugin.misspecified.bias.5000, " (", sim.df.m.p$plugin.misspecified.RMSE.5000,")", sep=""),
                          paste( sim.df.m.p$ifb.misspecified.bias.5000, " (", sim.df.m.p$ifb.misspecified.RMSE.5000,")", sep=""),
                          paste( sim.df.m.p$plugin.misspecified.bias.10000, " (", sim.df.m.p$plugin.misspecified.RMSE.10000,")", sep=""),
                          paste( sim.df.m.p$ifb.misspecified.bias.10000, " (", sim.df.m.p$ifb.misspecified.RMSE.10000,")", sep="")
                        )
publish.df.biasRMSE$`Cor NP` <- c( paste( sim.df.m.np$plugin.correct.bias.500, " (", sim.df.m.np$plugin.correct.RMSE.500,")", sep=""), 
                          paste( sim.df.m.np$ifb.correct.bias.500, " (", sim.df.m.np$ifb.correct.RMSE.500,")", sep=""),
                          paste( sim.df.m.np$plugin.correct.bias.5000, " (", sim.df.m.np$plugin.correct.RMSE.5000,")", sep=""),
                          paste( sim.df.m.np$ifb.correct.bias.5000, " (", sim.df.m.np$ifb.correct.RMSE.5000,")", sep=""),
                          paste( sim.df.m.np$plugin.correct.bias.10000, " (", sim.df.m.np$plugin.correct.RMSE.10000,")", sep=""),
                          paste( sim.df.m.np$ifb.correct.bias.10000, " (", sim.df.m.np$ifb.correct.RMSE.10000,")", sep="")
)
publish.df.biasRMSE$`Mis NP` <- c( paste( sim.df.m.np$plugin.misspecified.bias.500, " (", sim.df.m.np$plugin.misspecified.RMSE.500,")", sep=""),
                          paste( sim.df.m.np$ifb.misspecified.bias.500, " (", sim.df.m.np$ifb.misspecified.RMSE.500,")", sep=""),
                          paste( sim.df.m.np$plugin.misspecified.bias.5000, " (", sim.df.m.np$plugin.misspecified.RMSE.5000,")", sep=""),
                          paste( sim.df.m.np$ifb.misspecified.bias.5000, " (", sim.df.m.np$ifb.misspecified.RMSE.5000,")", sep=""),
                          paste( sim.df.m.np$plugin.misspecified.bias.10000, " (", sim.df.m.np$plugin.misspecified.RMSE.10000,")", sep=""),
                          paste( sim.df.m.np$ifb.misspecified.bias.10000, " (", sim.df.m.np$ifb.misspecified.RMSE.10000,")", sep="")
)
publish.df.biasRMSE


library(xtable)
print(xtable(publish.df.biasRMSE, caption="RMSE and bias from simulation with 500 repetitions.", label="tab:simresults"), include.rownames=FALSE)





# make plots for publication
dfsumRMSE <- as.data.frame ( c( rep("Parametric", 24), rep("Nonparametric", 24) )) ; names(dfsumRMSE) <- "PorNP"
dfsumRMSE$Algorithm_type <- (c( rep("Plug-in", 12), rep("Proposed", 12) ))
dfsumRMSE$X_Specification <- c( rep("Cor", 6), rep("Mis", 6))
dfsumRMSE$Stat <- c( rep("Bias", 3), rep("RMSE", 3))
dfsumRMSE$Sample_Sizes <- c( "100", "1000", "10000" )
dfsumRMSE <- dfsumRMSE %>% mutate_if(is.character,as.factor)
dfsumRMSE$PorNP <- factor(dfsumRMSE$PorNP, levels = c("Parametric", "Nonparametric"))
dfsumRMSE$Sample_Sizes <- factor(dfsumRMSE$Sample_Sizes, levels = c("100", "1000", "10000"))
dfsumRMSE$Value <- as.numeric(format(round(c ( as.numeric(t(sim.df.m.p)), as.numeric(t(sim.df.m.np)) ), digits=2), nsmall = 2) )

library(dplyr)
library(ggplot2)
library(ggthemes)

setwd("/Users/mariacuellar/Desktop/")
pdf(file="errors-RMSE-2018-08-26.pdf", width=8, height=3.5)
newdata_6 = dfsumRMSE[ which(dfsumRMSE$Stat=="RMSE"), ]
lp6 = ggplot(data=newdata_6, aes(x=Sample_Sizes, y=Value, colour=Algorithm_type, group = interaction(Algorithm_type, X_Specification)) ) + 
  geom_blank() + geom_line(aes(linetype=X_Specification), size=1) + facet_grid(Stat ~ PorNP, scales = "free_y", switch = "y") +
  ylim(0, max(dfsumRMSE$Value)) + xlab("Sample size") + labs(colour = "Algorithm type:", linetype = "X specification:") +
  ylab(NULL) + geom_point(size = 2)
lp6 = lp6 + theme_minimal() + theme(legend.position="bottom")
lp6
dev.off()






# make table for publication - coverage
publish.df.cov <- as.data.frame(matrix(NA,nrow=3,ncol=5))
names(publish.df.cov) <- c("Sample size", "Cor P", "Mis P", "Cor NP", "Mis NP")
publish.df.cov$`Sample size` <- c("500", "5000", "10000")
publish.df.cov$`Cor P` <- c( sim.df.cov.m.p$ifb.correct.cov.500, 
                             sim.df.cov.m.p$ifb.correct.cov.5000, 
                             sim.df.cov.m.p$ifb.correct.cov.1000
                            )
publish.df.cov$`Mis P` <- c( sim.df.cov.m.p$ifb.misspecified.cov.500, 
                             sim.df.cov.m.p$ifb.misspecified.cov.5000, 
                             sim.df.cov.m.p$ifb.misspecified.cov.1000
                            )
publish.df.cov$`Cor NP` <- c( sim.df.cov.m.np$ifb.correct.cov.500, 
                              sim.df.cov.m.np$ifb.correct.cov.5000, 
                              sim.df.cov.m.np$ifb.correct.cov.1000
                            )
publish.df.cov$`Mis NP` <- c( sim.df.cov.m.np$ifb.misspecified.cov.500, 
                              sim.df.cov.m.np$ifb.misspecified.cov.5000, 
                              sim.df.cov.m.np$ifb.misspecified.cov.1000
                            )
publish.df.cov


library(xtable)
print(xtable(publish.df.cov, caption="Coverage from simulation with 500 repetitions using nonparametric estimation for nuisance parameters.", 
             label="tab:simresultscoverage"), include.rownames=FALSE)





# TO DO:
# WHY ARE SOME OF THE PLUGINS PC'S STILL NEGATIVE? RIGHT NOW I'M USING ABS() IN THE FUNCTION.
# COVERAGE GETS SMALLER WITH SAMPLE SIZE. WHY IS THIS HAPPENING?
