rm(list = ls())
graphics.off()

# install.packages('DEoptim')
# install.packages('truncnorm')

## Set the seed (for reproducibility)
set.seed(1234)

#==========================================================================================================================================================
# Objective Function                                                                                                                                      #
# Objective function for use with DEoptim                                                                                                                 #
# Carl Fredrick G. Aquino / cga5030 / 994085920                                                                                                           #
#==========================================================================================================================================================
findShape <- function(pars,quantiles,myshape){

  peg <- c(.05,.17,.50,.83,.95)
  
  library(truncnorm)
  
  # can do IF test later to specify this by bringing in another variable from driver
  # ie, if (test==1/2/3) do norm/lnorm/truncnorm
  
    if (myshape==1){test <- qnorm(peg,mean = pars[1],sd = pars[2])}
    if (myshape==2){test <- qlnorm(peg,meanlog = pars[1],sdlog = pars[2])}
    if (myshape==3){test <- qtruncnorm(peg,a=-.25,b=Inf,mean = pars[1],sd = pars[2])}
  
  # RMSE = sqrt(  sum(observation - known)^2 / N    )
  N=length(quantiles)
  
  RMSE = sqrt(sum((quantiles-test)^2)/N)
  
  return(RMSE)
}
#==========================================================================================================================================================

if (TRUE){
#==========================================================================================================================================================
# Driver                                                                                                                                   #
# Driver for shape sensitivity - IN METERS!
# Carl Fredrick G. Aquino / cga5030 / 994085920                                                                                                           #
#==========================================================================================================================================================

# Configure !!! - ENTER TRUE or FALSE FOR DESIRED ANALYSIS
shape.norm = TRUE # true if normal distribution
shape.lnorm = TRUE # true if lognormal distribution
shape.truncnorm = TRUE # true if truncated normal distribution

# specify quantiles
quantiles_2100H <- c(-.11,.2,.21,.57,1.32)
  
library(DEoptim)
library(truncnorm)

if (shape.norm){
  
  # for normal distribution
  myshape=1
  
  ## Set up DEoptim
  npars=2
  lower.bounds = rep(NA,npars)
  upper.bounds = rep(NA,npars)
  NP_scale = 11  
  n_iter = 200 # Tony ran 100 # Vivek tends to do 2000 - 5000. # Start with 10 and then see what happens
  
  lower.bounds <- c(-100,0)
  upper.bounds <- c(100,100)
  
  DEoptim.output <- DEoptim(findShape,lower.bounds,upper.bounds,control=list(NP=NP_scale*npars, itermax=n_iter, parallelType=1,
                                                                             parVar = list("quantiles_2100H")),
                            quantiles=quantiles_2100H,
                            myshape=myshape
  )
  
  grid = seq(-.50,3.00,.01)
  
  plot(grid,pnorm(grid,mean=DEoptim.output$optim$bestmem[1],sd=DEoptim.output$optim$bestmem[2]),type="l",xlab="x",ylab="CDF")
  points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
  
  OUTPUT_norm_pars <- DEoptim.output$optim$bestmem
  OUTPUT_norm_rmse <- DEoptim.output$optim$bestval
  
}

if (shape.lnorm){
  
  # for normal distribution
  myshape=2
  
  ## Set up DEoptim
  npars=2
  lower.bounds = rep(NA,npars)
  upper.bounds = rep(NA,npars)
  NP_scale = 11  
  n_iter = 300 # Tony ran 100 # Vivek tends to do 2000 - 5000. # Start with 10 and then see what happens
  
  lower.bounds <- c(-1000,0)
  upper.bounds <- c(1000,1000)
  
  DEoptim.output <- DEoptim(findShape,lower.bounds,upper.bounds,control=list(NP=NP_scale*npars, itermax=n_iter, parallelType=1,
                                                                             parVar = list("quantiles_2100H")),
                            quantiles=quantiles_2100H,
                            myshape=myshape
  )
  
  grid = seq(-.50,3.00,.01)
  
  plot(grid,plnorm(grid,meanlog=DEoptim.output$optim$bestmem[1],sdlog=DEoptim.output$optim$bestmem[2]),type="l",xlab="x",ylab="CDF")
  points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
  
  OUTPUT_lnorm_pars <- DEoptim.output$optim$bestmem
  OUTPUT_lnorm_rmse <- DEoptim.output$optim$bestval
  
}


if (shape.truncnorm){
  
  # for truncated normal distribution
  myshape=3
  
  ## Set up DEoptim
  npars=2
  lower.bounds = rep(NA,npars)
  upper.bounds = rep(NA,npars)
  NP_scale = 11  
  n_iter = 700 # Tony ran 100 # Vivek tends to do 2000 - 5000. # Start with 10 and then see what happens
  
  lower.bounds <- c(-4000,-1000)
  upper.bounds <- c(1000,1000)
  
  DEoptim.output <- DEoptim(findShape,lower.bounds,upper.bounds,control=list(NP=NP_scale*npars, itermax=n_iter, parallelType=1,
                                                                             parVar = list("quantiles_2100H")),
                            quantiles=quantiles_2100H,
                            myshape=myshape
                            )
  
  grid = seq(-.50,3.00,.01)
  
  plot(grid,ptruncnorm(grid,a=-.25,b=Inf,mean=DEoptim.output$optim$bestmem[1],sd=DEoptim.output$optim$bestmem[2]),type="l",xlab="x",ylab="CDF")
  points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
  
  OUTPUT_trunc_pars <- DEoptim.output$optim$bestmem
  OUTPUT_trunc_rmse <- DEoptim.output$optim$bestval
}

#==========================================================================================================================================================
}
  
  
# PROCESS FINDINGS
if (TRUE){

  grid = seq(-.50,3.00,.01)
  
  singleplot=FALSE
  if (singleplot){
    
    normal = c(74,50.767063)
    RMSE_norm = 16.389841
    plot(grid,pnorm(grid,92.5,52),type="l",xlab="x",ylab="f(x)")
    points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
    
    log_normal = c(3.925, 0.774)
    RMSE_lnorm = 4.494364
    plot(grid,plnorm(grid,meanlog = log_normal[1],sdlog = log_normal[2]),type="l")
    points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
    
    trunc_normal = c(-127.7569, 125.5338)
    RMSE_truncnorm = 3.929016
    plot(grid,ptruncnorm(grid,a=0,b=Inf,mean=trunc_normal[1],sd=trunc_normal[2]), type="l",
         xlab="SLR [cm]",ylab="CDF")
    points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
    legend("bottomright",
           legend=c("expert quantiles","truncated normal"),
           col=c("black","black"),
           lty=c(0,1),
           pch=c(19,NA)    
           )
    
    trunc_normal_v2 = c(43.526701, 71.975001)
    RMSE_truncnorm_v2 = 11.602644
    plot(grid,ptruncnorm(grid,a=-25,b=Inf,mean=trunc_normal_v2[1],sd=trunc_normal_v2[2]), type="l",
         xlab="SLR [cm]",ylab="CDF")
    points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
    legend("bottomright",
           legend=c("expert quantiles","truncated normal_v2"),
           col=c("black","black"),
           lty=c(0,1),
           pch=c(19,NA)    
          )
  }
  
  multiplot=TRUE
  if (multiplot){
    normal = OUTPUT_norm_pars
    log_normal = OUTPUT_lnorm_pars
    trunc_normal = OUTPUT_trunc_pars
    
    plot(grid,pnorm(grid,mean=normal[1],sd=normal[2]),type="l",xlab="SLR [m]",ylab="CDF",lwd=2,lty=1)
    lines(grid,plnorm(grid,meanlog = log_normal[1],sdlog = log_normal[2]),col="red",lwd=2)
    lines(grid,ptruncnorm(grid,a=-.25,b=Inf,mean=trunc_normal[1],sd=trunc_normal[2]),col="blue",lwd=2, lty=2)
    points(quantiles_2100H,c(.05,.17,.50,.83,.95),pch=19)
    
    legend("bottomright",
           legend=c("expert quantiles","normal","log-normal","truncated normal"),
           col=c("black","black","red","blue"),
           lty=c(0,1,1,2),
           pch=c(19,NA,NA,NA),
           lwd=2
           )
    }
}
