##==============================================================================
## AR1 model for errors (centered at 0 -- X(t)=rho*X(t-1)+eps(t), eps(t)~N(0,sigma1)
##		For clarity -- sigma.proc is for AR(1) process;
##					-- sigma1 (sampled) is the whitened innovation sigma
## Estimate the log likelihood of the AR1 process

if(TRUE){
  ## APPROX AR1?
  logl.ar1 = function(r,sigma1,rho1,eps1=0) # default obs error is 0
  {
    n = length(r) # r is the residuals
    if(length(eps1)==1) eps1 = rep(eps1,n)
    
    logl=0
    if(n>1) {
      w = r[2:n] - rho1*r[1:(n-1)] # this process whitens the residuals
      logl = logl + sum(dnorm(w,sd=sqrt((sigma1)^2+(eps1[c(-1)])^2),log=TRUE)) # add in the sum of
      # density of the whitened residuals with a standard deviation of the
      # variance and the obs. errors
    }
    return(logl)
  }
  ##
} else {
  ## EXACT AR1?
  library(mvtnorm)
  logl.ar1 <-
    function(r,sigma1,rho1,eps1) # default obs error is 0. sigma1 is standard error.
    {
      library(mvtnorm)
      
      n = length(r)
      sigma.proc = sigma1/sqrt(1-rho1^2) # stationary process variance sigma.proc^2
      if(all(eps1==0)){
        logl = dnorm(r[1],sd=sigma.proc,log=TRUE)
        if(n>1) {
          w = r[2:n] - rho1*r[1:(n-1)] # whitened residuals
          # logl = logl + sum(dnorm(w,sd=sigma1,log=TRUE))
          # This is what we had before to make the computation faster.
          # This approximation should not change the result, but it is worth trying
          logl = logl + sum(dnorm(w,sd=sqrt(sigma1^2+eps1[-1]^2),log=TRUE))
        }
      }else{
        H <- abs(outer(1:n, 1:n, "-"))
        v = sigma.proc^2*rho1^H
        if(length(eps1)>1) {v = v+diag(eps1^2)
        } else {v = v+diag(rep(eps1^2,n))}
        # Need R package "mvtnorm"
        logl = dmvnorm(r,sigma=v,log=TRUE)
      }
      return(logl)
    }
  ##
}

##==============================================================================
## (log of the) prior probability for BRICK
log.pri.BRICK = function(parameters.in , parnames.in, bound.lower.in, bound.upper.in,
                   shape.in, scale.in )
{
  
  # Pluck off the model and statistical parameters (only ones without uniform prior)
  ind.S      = match("S",parnames.in)
  ind.invtau = match("invtau.te",parnames.in)
  ind.rho.gsic = match("rho.gsic",parnames.in)
  lpri.S      = 0
  lpri.invtau = 0
  
  in.range = all(parameters.in > bound.lower.in) & all(parameters.in < bound.upper.in)
  
  if(in.range){
    lpri.uni = 0									# Sum of all uniform priors (log(1)=0)
    #    lpri.S = log(dcauchy(parameters.in[ind.S],location=3,scale=2) / 	# S has truncated Cauchy(3,2) prior
    #					(pcauchy(bound.upper[ind.S],location=3,scale=2)-pcauchy(bound.lower[ind.S],location=3,scale=2)))
    if(!is.na(ind.invtau)) {lpri.invtau = dgamma( parameters.in[ind.invtau], shape=shape.in, scale=scale.in, log=TRUE)}
    lpri = lpri.uni + lpri.S + lpri.invtau
    #    lpri = lpri.uni
  } else {
    lpri = -Inf
  }
  
  return(lpri)
}

##==============================================================================
## Log-prior for DAIS
log.pri.DAIS = function( parameters.in,
                    parnames.in,
                    alpha.var.in,
                    beta.var.in,
                    shape.lambda=NULL,
                    rate.lambda=NULL,
                    shape.Tcrit=NULL,
                    rate.Tcrit=NULL,
                    l.aisfastdy=l.aisfastdy,
                    bound.lower.in,
                    bound.upper.in
)
{
  
  if(l.aisfastdy) {
    ind.lambda = match("lambda",parnames.in)
    lambda =parameters.in[match("lambda",parnames.in)]
    Tcrit =parameters.in[match("Tcrit",parnames.in)]
    ind.var = match("var.dais",parnames.in)
    var.dais =parameters.in[ind.var]
    # var.dais has inverse gamma prior, so there is a lower bound at 0 but no
    # upper bound
    if(fd.priors=='uniform') {
      in.range = all(parameters.in[1:(ind.var-1)] < bound.upper.in[1:(ind.var-1)]) &
        all(parameters.in > bound.lower.in)
    }
    if(fd.priors=='gamma') {
      in.range = all(parameters.in[1:(ind.lambda-1)] < bound.upper.in[1:(ind.lambda-1)]) &
        all(parameters.in[1:(ind.lambda-1)] > bound.lower.in[1:(ind.lambda-1)]) &
        all(parameters.in[ind.var] > bound.lower.in[ind.var]) &
        all(parameters.in[ind.lambda] > 0) &
        all(parameters.in[ind.lambda+1] < 0)
    }
    var_pri = 0
    lambda_pri = 0
    Tcrit_pri = 0
    if(in.range) {
      var_pri = (-alpha.var.in - 1)*log(var.dais) + (-beta.var.in/var.dais)
      if(fd.priors=='gamma') {
        lambda_pri = dgamma(x=lambda, shape=shape.lambda, rate=rate.lambda, log=TRUE)
        Tcrit_pri = dgamma(x=-Tcrit, shape=shape.Tcrit, rate=rate.Tcrit, log=TRUE)
      }
      lpri=0 + var_pri + lambda_pri + Tcrit_pri
    } else {
      lpri = -Inf
    }
  } else {
    ind.var = match("var.dais",parnames.in)
    var.dais =parameters.in[ind.var]
    # var.dais has inverse gamma prior, so there is a lower bound at 0 but no
    # upper bound
    in.range = all(parameters.in[1:(ind.var-1)] < bound.upper.in[1:(ind.var-1)]) &
      all(parameters.in > bound.lower.in)
    var_pri = 0
    if(in.range) {
      var_pri = (-alpha.var.in - 1)*log(var.dais) + (-beta.var.in/var.dais)
      lpri=0 + var_pri
    } else {
      lpri = -Inf
    }
  }
  lpri
}

##==============================================================================
## (log of the) density of expert assessment for log-normal shape 

log_exp_lognorm <- function(model_out){
  # expert distribution - imputed from Bamber et al. 2019, Table 2 "2100 H"
  meanlog = 0.147
  sdlog = 0.444
  
  # get BRICK outputs for 2100
  slr_2000 <- model_out$slr.out[151]
  slr_2100 <- model_out$slr.out[251]
  test_slr <- slr_2100 - slr_2000
  
  # log-likelihood given expert assessment
  dlnorm(test_slr, meanlog = meanlog, sdlog=sdlog, log=TRUE)
}

##==============================================================================
## Log-likelihood
log.lik.DAISpaleo = function( parameters.in,
                              parnames.in,
                              obs.in,
                              obs.err.in,
                              obs.step.in,
                              ind.norm.in,
                              SL.in,
                              dSL.in,
                              trends.ais.in,
                              trends.err.in,
                              ind.trends.in,
                              slope.Ta2Tg.in=1,
                              intercept.Ta2Tg.in=0,
                              l.aisfastdy=TRUE,
                              Tg.in=NULL,
                              Ta.in=NULL,
                              Toc.in=NULL
)
{
  
  anto.a=parameters.in[match("anto.a",parnames.in)]
  anto.b=parameters.in[match("anto.b",parnames.in)]
  gamma =parameters.in[match("gamma" ,parnames.in)]
  alpha =parameters.in[match("alpha.dais" ,parnames.in)]
  mu =parameters.in[match("mu" ,parnames.in)]
  nu =parameters.in[match("nu" ,parnames.in)]
  P0 =parameters.in[match("P0" ,parnames.in)]
  kappa =parameters.in[match("kappa.dais" ,parnames.in)]
  f0 =parameters.in[match("f0" ,parnames.in)]
  h0 =parameters.in[match("h0" ,parnames.in)]
  c =parameters.in[match("c" ,parnames.in)]
  b0 =parameters.in[match("b0" ,parnames.in)]
  slope =parameters.in[match("slope" ,parnames.in)]
  var.dais =parameters.in[match("var.dais" ,parnames.in)]
  if(l.aisfastdy) {
    Tcrit =parameters.in[match("Tcrit" ,parnames.in)]
    lambda =parameters.in[match("lambda" ,parnames.in)]
  } else {
    Tcrit = NULL
    lambda = NULL
  }
  
  if(TRUE){
    dais.out = daisanto_fastdynF(
      anto.a=anto.a, anto.b=anto.b,
      slope.Ta2Tg=slope.Ta2Tg.in, intercept.Ta2Tg=intercept.Ta2Tg.in,
      gamma=gamma  , alpha=alpha  ,
      mu=mu        , nu=nu        ,
      P0=P0        , kappa=kappa  ,
      f0=f0        , h0=h0        ,
      c=c          , b0=b0        ,
      slope=slope  ,
      Tcrit=Tcrit  , lambda=lambda,
      l.aisfastdy=l.aisfastdy,
      Tg=Tg.in     ,
      SL=SL.in     , dSL=dSL.in)
  }
  # if(FALSE){
  #   dais.out = dais_fastdynF(
  #     gamma=gamma  , alpha=alpha  ,
  #     mu=mu        , nu=nu        ,
  #     P0=P0        , kappa=kappa  ,
  #     f0=f0        , h0=h0        ,
  #     c=c          , b0=b0        ,
  #     slope=slope  , SL=SL.in     ,
  #     dSL=dSL.in   , Ta=Ta.in     ,
  #     Toc=Toc.in   , includes_dSLais=1,
  #     Tcrit=Tcrit  , lambda=lambda,
  #     l.aisfastdy=l.aisfastdy
  #   )
  # }
  
  dais.out.norm = dais.out$Vais - mean(dais.out$Vais[ind.norm.in])
  
  llik0  = 0
  
  # Calculate the residuals
  
  # First check if the modeled AIS-only SLR exceeds total
  # Use first 20 years and only check residuals after that for exceeding total SLR
  # because there will be noise around the trend in the first 20 years; don't want
  # to throw out the run because of a little noise.
  
  resid.sl = ( SL.new-mean(SL.new[1:20]) ) - ( dais.out$Vais[midx.sl]-mean(dais.out$Vais[midx.sl[1:20]]) )
  
  # Only admit runs which fit the Pliocene window, to 1-sigma
  resid.lig = abs(obs.in[1]+0-dais.out.norm[obs.step.in[1]])
  
  if(!is.na(resid.lig)){
    if(all(resid.sl[20:length(resid.sl)]>0) & (resid.lig < 2.0*obs.err.in[1])){
      
      # Calculate the modeled trend for the IPCC periods
      # 1993-2010
      x=seq(ind.trends[1,1],ind.trends[1,2]); y=dais.out.norm[ind.trends[1,1]:ind.trends[1,2]]
      barx=mean(x); bary=mean(y)
      trend.1993 = sum( (x-barx)*(y-bary) )/sum( (x-barx)^2 )
      # 1992-2001
      x=seq(ind.trends[2,1],ind.trends[2,2]); y=dais.out.norm[ind.trends[2,1]:ind.trends[2,2]]
      barx=mean(x); bary=mean(y)
      trend.1992 = sum( (x-barx)*(y-bary) )/sum( (x-barx)^2 )
      # 2002-2011
      x=seq(ind.trends[3,1],ind.trends[3,2]); y=dais.out.norm[ind.trends[3,1]:ind.trends[3,2]]
      barx=mean(x); bary=mean(y)
      trend.2002 = sum( (x-barx)*(y-bary) )/sum( (x-barx)^2 )
      
      resid = c( obs.in[1] - dais.out.norm[obs.step.in[1]] ,
                 obs.in[2] - dais.out.norm[obs.step.in[2]] ,
                 obs.in[3] - dais.out.norm[obs.step.in[3]] ,
                 obs.in[4] - (dais.out.norm[obs.step.in[4]]-dais.out.norm[i1992]) ,
                 trends.ais.in[1] - trend.1993 ,
                 trends.ais.in[2] - trend.1992 ,
                 trends.ais.in[3] - trend.2002
      )
      
      # Calculate the likelihood. The observations are not correlated. They are independent
      llik0 = sum (dnorm(resid, mean=rep(0,length(resid)), sd = c(sqrt(var.dais + obs.err.in[1:4]^2), trends.err.in), log=TRUE))
      #    llik0 = sum (dnorm(resid, mean=rep(0,length(resid)), sd = c(sqrt(var.dais + obs.err.in[1:4]^2)), log=TRUE))
      
      llik.DAISpaleo = llik0 # assume residuals are independent
    } else {
      llik.DAISpaleo = -Inf
    }
  } else {
    llik.DAISpaleo=-Inf
  }
  llik.DAISpaleo
}

##==============================================================================
## rest of the statistical model
##==============================================================================
log.lik.BRICK = function( 
                          parameters.in,
                          parnames.in,
                          forcing.in,
                          bound.lower.in,
                          bound.upper.in,
                          l.project=FALSE,
                          rho.simple.in=NULL,
                          sigma.simple.in=NULL,
                          shape.in,
                          scale.in,
                          slope.Ta2Tg.in=1,
                          intercept.Ta2Tg.in=0,
                          mod.time,
                          ind.norm.data,
                          ind.norm.sl,
                          midx,
                          oidx,
                          obs,
                          obs.err,
                          trends.te,
                          luse.brick,
                          i0,
                          l.aisfastdy=TRUE,
                          
                          brick.out
  
){
  
  ## Calculate contribution from DOECLIM temperature
  llik.temp = 0
  if(!is.null(oidx$temp) & luse.brick[,"luse.doeclim"]) {
    
    ## Grab the DOECLIM parameters
    parnames.doeclim   =c("S" ,"kappa.doeclim","alpha.doeclim","T0"  ,"H0" ,"sigma.T","sigma.H","rho.T","rho.H")	# parameters names
    p0.doeclim         =c(3.1 , 3.5   , 1.1           , -0.06, -33 , 0.1     , 2       , 0.55  , 0.9   )	# initial parameter guesses
    S            =p0.doeclim[match("S"            ,parnames.doeclim)]
    kappa.doeclim=p0.doeclim[match("kappa.doeclim",parnames.doeclim)]
    alpha.doeclim=p0.doeclim[match("alpha.doeclim",parnames.doeclim)]
    T0           =p0.doeclim[match("T0"           ,parnames.doeclim)]
    H0           =p0.doeclim[match("H0"           ,parnames.doeclim)]
    sigma.T      =p0.doeclim[match("sigma.T"      ,parnames.doeclim)]
    sigma.H      =p0.doeclim[match("sigma.H"      ,parnames.doeclim)]
    rho.T        =p0.doeclim[match("rho.T"        ,parnames.doeclim)]
    rho.H        =p0.doeclim[match("rho.H"        ,parnames.doeclim)]
    
    # Calculate the DOECLIM temperature residuals; apply AR1 error model
    resid.temp= obs$temp[oidx$temp] - (brick.out$doeclim.out$temp[midx$temp]+T0)
    llik.temp = logl.ar1(resid.temp, sigma.T, rho.T, obs.err$temp[oidx$temp]) # AR(1)
    
  }
  
  ## Calculate contribution from DOECLIM ocean heat
  llik.ocheat = 0
  if(!is.null(oidx$ocheat) & luse.brick[,"luse.doeclim"]) {
    
    ## Grab the DOECLIM parameters
    parnames.doeclim   =c("S" ,"kappa.doeclim","alpha.doeclim","T0"  ,"H0" ,"sigma.T","sigma.H","rho.T","rho.H")	# parameters names
    p0.doeclim         =c(3.1 , 3.5   , 1.1           , -0.06, -33 , 0.1     , 2       , 0.55  , 0.9   )	# initial parameter guesses
    S            =p0.doeclim[match("S"            ,parnames.doeclim)]
    kappa.doeclim=p0.doeclim[match("kappa.doeclim",parnames.doeclim)]
    alpha.doeclim=p0.doeclim[match("alpha.doeclim",parnames.doeclim)]
    T0           =p0.doeclim[match("T0"           ,parnames.doeclim)]
    H0           =p0.doeclim[match("H0"           ,parnames.doeclim)]
    sigma.T      =p0.doeclim[match("sigma.T"      ,parnames.doeclim)]
    sigma.H      =p0.doeclim[match("sigma.H"      ,parnames.doeclim)]
    rho.T        =p0.doeclim[match("rho.T"        ,parnames.doeclim)]
    rho.H        =p0.doeclim[match("rho.H"        ,parnames.doeclim)]
    
    # Calculate the DOECLIM ocean heat residuals; apply AR1 error model
    resid.ocheat= obs$ocheat[oidx$ocheat] - (brick.out$doeclim.out$ocheat[midx$ocheat]+H0)
    llik.ocheat = logl.ar1(resid.ocheat, sigma.H, rho.H, obs.err$ocheat[oidx$ocheat]) # AR(1)
    
  }
  
  ## Calculate contribution from GSIC SLR
  llik.gsic = 0
  if(!is.null(oidx$gsic) & luse.brick[,"luse.gsic"]) {
    
    # Grab the GSIC statistical parameters
    sigma.gsic=parameters.in[match("sigma.gsic",parnames.in)]
    rho.gsic  =parameters.in[match("rho.gsic"  ,parnames.in)]
    
    # Calculate the GSIC residuals; apply AR1 error model
    resid.gsic= obs$gsic[oidx$gsic] - brick.out$gsic.out[midx$gsic]
    llik.gsic = logl.ar1(resid.gsic, sigma.gsic, rho.gsic, obs.err$gsic[oidx$gsic]) # AR(1)
  }
  
  ## Calculate contribution from thermosteric expansion
  llik.te = 0
  if(luse.brick[,"luse.te"] | luse.brick[,"luse.tee"]) {
    
    # Calculate the SLR residuals - only proceed if all TE SLR < total SLR
    # (all after the first 20 years, that is, because they are the 0 point)
    
    resid.sl.te = (obs$sl[oidx$sl]-mean(obs$sl[oidx$sl[1:20]])) -
      (brick.out$te.out[midx$sl]-mean(brick.out$te.out[midx$sl[1:20]]))
    
    if(all(resid.sl.te[20:length(resid.sl.te)]>0)){
      
      # Note 1: the trends from IPCC are in mm/year, and model output is m
      # Note 2: these calculate the least squares regression slope coefficients. It
      # is more than twice as fast to calcualte by hand like this than to use R's
      # "lm(...)" function.
      # Note 3: Need 1000*trends.mod because they're in meters, but trends.te is mm
      
      trends.mod = rep(0, nrow(trends.te))
      for (i in 1:nrow(trends.te)) {
        x = seq(trends.te[i,6],trends.te[i,7]);              barx = mean(x);
        y = brick.out$te.out[trends.te[i,6]:trends.te[i,7]]; bary = mean(y);
        trends.mod[i] = sum( (x-rep(barx,length(x)))*(y-rep(bary,length(y))))/sum( (x-rep(barx,length(x)))^2 )
      }
      resid.trends = 1000*trends.mod - trends.te[,1]
      err.trends   = 0.5*(trends.te[,3]-trends.te[,2])
      llik.te = sum (dnorm(resid.trends, mean=rep(0,length(resid.trends)), sd = sqrt(err.trends^2), log=TRUE))
      
    } else {
      llik.te = -Inf
    }
  }
  
  ## Calculate contribution from SIMPLE (Greenland Ice Sheet)
  llik.simple = 0
  if(!is.null(oidx$gis) & luse.brick[,"luse.simple"]) {
    
    # Grab the SIMPLE statistical parameters
    sigma.simple=parameters.in[match("sigma.simple",parnames.in)]
    rho.simple  =parameters.in[match("rho.simple"  ,parnames.in)]
    
    # Overwrite the SIMPLE statistical parameters if values were fed into MCMC
    if(!is.null(rho.simple.in  )) rho.simple  =rho.simple.in
    if(!is.null(sigma.simple.in)) sigma.simple=sigma.simple.in
    
    # Calibrate SIMPLE based on GIS data alone?
    resid.simple = obs$gis[oidx$gis] - brick.out$simple.out$sle.gis[midx$gis] #Calculate the residuals
    
    if(!all(is.finite(resid.simple))) {
      llik.simple = -Inf
    } else {
      llik.simple  = logl.ar1(r=resid.simple, sigma1=sigma.simple,
                              rho1=rho.simple, eps1=obs.err$gis) # AR(1) #Set up the likelihood function
    }
  }
  
  ## Calculate contribution from Antarctic ice sheet
  llik.dais = 0
  if(luse.brick[,"luse.dais"]) {
    
    # Calculate the SLR residuals - only proceed if all AIS SLR < total SLR
    # (all after the first 20 years, that is, because they are the 0 point)
    
    if(is.na(brick.out$dais.out)){
        
      llik.dais = -Inf
      
    } else {
    
      resid.sl.ais = (obs$sl[oidx$sl]-mean(obs$sl[oidx$sl[1:20]])) -
        #(brick.out[[5]][1][midx$sl]-mean(brick.out[[5]][1][midx$sl[1:20]]))
        (brick.out$dais.out$Vais[midx$sl]-mean(brick.out$dais.out$Vais[midx$sl[1:20]]))
      
      if(all(resid.sl.ais[20:length(resid.sl.ais)]>0)){
        llik.dais = 0
      } else {
        llik.dais = -Inf
      }
    
    }
    
    
  }
  
  # Calculate contribution from total sea level rise
  llik.sl = 0
  if(FALSE) {
    if(!is.null(oidx$sl)) {
      
      resid.sl = brick.out$slr.out[midx$sl] - obs$sl[oidx$sl]
      llik.sl  = sum (dnorm(resid.sl, mean=rep(0,length(resid.sl)), sd = sqrt(obs.err$sl[oidx$sl]^2), log=TRUE))
      
    }
  }

  # Assume residual time series are independent
  llik.BRICK = llik.temp + llik.ocheat + llik.gsic + llik.te + llik.simple + llik.dais + llik.sl
  
  return(llik.BRICK)

}

##==============================================================================
## (log of the) posterior distribution:  posterior ~ likelihood * prior
log.post = function(  parameters.in,
                      parnames.in,
                      forcing.in,
                      bound.lower.in,
                      bound.upper.in,
                      l.project=FALSE,
                      rho.simple.in=NULL,
                      sigma.simple.in=NULL,
                      shape.in,
                      scale.in,
                      slope.Ta2Tg.in=1,
                      intercept.Ta2Tg.in=0,
                      mod.time,
                      ind.norm.data,
                      ind.norm.sl,
                      midx,
                      oidx,
                      obs,
                      obs.err,
                      trends.te,
                      luse.brick,
                      i0,
                      l.aisfastdy=TRUE,
                      obs.in, obs.err.in, obs.step.in,
                      trends.ais.in , trends.err.in , ind.trends.in ,
                      ind.norm.in,
                      alpha.var,
                      beta.var,
                      shape.lambda=NULL,
                      rate.lambda=NULL,
                      shape.Tcrit=NULL,
                      rate.Tcrit=NULL,
                      Tg.in=NULL,
                      Toc.in=NULL,
                      Ta.in=NULL,
                      SL.in,
                      dSL.in,
                      experts,
                      alldata
){
  
  # saveRDS(parameters.in, file="parameters_in.RDS")
  
  lpri.BRICK = log.pri.BRICK( parameters.in=parameters.in,
                              parnames.in=parnames.in,
                              bound.lower.in=bound.lower.in,
                              bound.upper.in=bound.upper.in,
                              shape.in=shape.in,
                              scale.in=scale.in)
  
  lpri.DAIS = log.pri.DAIS(parameters.in                , parnames.in=parnames.in      ,
                           alpha.var.in=alpha.var       , beta.var.in=beta.var         ,
                           bound.lower.in=bound.lower.in, bound.upper.in=bound.upper.in,
                           shape.lambda=shape.lambda    , rate.lambda=rate.lambda      ,
                           shape.Tcrit=shape.Tcrit      , rate.Tcrit=rate.Tcrit        ,
                           l.aisfastdy=l.aisfastdy                                      )
  
  lpri = lpri.BRICK + lpri.DAIS
  
  if(is.finite(lpri)){
    lpost <- lpri
      
      # if expert assessment inversion is included, add log-expert assessment density
      
      brick.out = brick_model(  parameters.in=parameters.in,
                                parnames.in=parnames.in,
                                forcing.in=forcing.in,
                                l.project=l.project,
                                slope.Ta2Tg.in=slope.Ta2Tg.in,
                                intercept.Ta2Tg.in=intercept.Ta2Tg.in,
                                mod.time=mod.time,
                                ind.norm.data=ind.norm.data,
                                ind.norm.sl=ind.norm.sl,
                                luse.brick=luse.brick,
                                i0=i0,
                                l.aisfastdy=l.aisfastdy
      )
        
      if (experts){
        lpost <- lpost + log_exp_lognorm(brick.out)
      }
      
      if (alldata){
      
        llik.BRICK = log.lik.BRICK( 
          
                                parameters.in=parameters.in,
                                parnames.in=parnames.in,
                                forcing.in=forcing.in,
                                l.project=l.project,
                                rho.simple.in=rho.simple.in,
                                sigma.simple.in=sigma.simple.in,
                                slope.Ta2Tg.in=slope.Ta2Tg.in,
                                intercept.Ta2Tg.in=intercept.Ta2Tg.in,
                                mod.time=mod.time,
                                ind.norm.data=ind.norm.data,
                                ind.norm.sl=ind.norm.sl,
                                midx=midx,
                                oidx=oidx,
                                obs=obs,
                                obs.err=obs.err,
                                trends.te=trends.te,
                                luse.brick=luse.brick,
                                i0=i0,
                                l.aisfastdy=l.aisfastdy,
          
                                brick.out=brick.out
                                 )  
        
        llik.DAISpaleo =         log.lik.DAISpaleo( 
          
                                  parameters.in=parameters.in  , parnames.in=parnames.in ,
                                  obs.in=obs.in                , obs.err.in=obs.err.in   ,
                                  obs.step.in=obs.step.in      , ind.trends.in=ind.trends.in ,
                                  trends.ais.in=trends.ais.in  , trends.err.in = trends.err.in ,
                                  slope.Ta2Tg.in=slope.Ta2Tg.in, intercept.Ta2Tg.in=intercept.Ta2Tg.in,
                                  Tg.in=Tg.in                  , l.aisfastdy=l.aisfastdy,
                                  #Ta.in=Ta.in                  , Toc.in=Toc.in     ,
                                  SL.in=SL.in                  , dSL.in=dSL.in     ,
                                  ind.norm.in=ind.norm.in )
        
      lpost<-lpost+llik.BRICK+llik.DAISpaleo
          
      }
      
      
  

} else {
    lpost = -Inf
  }
  
  # print(paste0("lpost: ", lpost))
  return(lpost)
}

##==============================================================================
## (log of the) posterior distribution:  posterior ~ likelihood * prior
neg.log.post = function(  parameters.in,
                      parnames.in,
                      forcing.in,
                      bound.lower.in,
                      bound.upper.in,
                      l.project=FALSE,
                      rho.simple.in=NULL,
                      sigma.simple.in=NULL,
                      shape.in,
                      scale.in,
                      slope.Ta2Tg.in=1,
                      intercept.Ta2Tg.in=0,
                      mod.time,
                      ind.norm.data,
                      ind.norm.sl,
                      midx,
                      oidx,
                      obs,
                      obs.err,
                      trends.te,
                      luse.brick,
                      i0,
                      l.aisfastdy=TRUE,
                      obs.in, obs.err.in, obs.step.in,
                      trends.ais.in , trends.err.in , ind.trends.in ,
                      ind.norm.in,
                      alpha.var,
                      beta.var,
                      shape.lambda=NULL,
                      rate.lambda=NULL,
                      shape.Tcrit=NULL,
                      rate.Tcrit=NULL,
                      Tg.in=NULL,
                      Toc.in=NULL,
                      Ta.in=NULL,
                      SL.in,
                      dSL.in,
                      experts,
                      alldata
){
  
  # evaluate log posterior
  lp = log.post (  parameters.in,
                        parnames.in,
                        forcing.in,
                        bound.lower.in,
                        bound.upper.in,
                        l.project,
                        rho.simple.in,
                        sigma.simple.in,
                        shape.in,
                        scale.in,
                        slope.Ta2Tg.in,
                        intercept.Ta2Tg.in,
                        mod.time,
                        ind.norm.data,
                        ind.norm.sl,
                        midx,
                        oidx,
                        obs,
                        obs.err,
                        trends.te,
                        luse.brick,
                        i0,
                        l.aisfastdy,
                        obs.in, obs.err.in, obs.step.in,
                        trends.ais.in , trends.err.in , ind.trends.in ,
                        ind.norm.in,
                        alpha.var,
                        beta.var,
                        shape.lambda,
                        rate.lambda,
                        shape.Tcrit,
                        rate.Tcrit,
                        Tg.in,
                        Toc.in,
                        Ta.in,
                        SL.in,
                        dSL.in,
                        experts,
                        alldata
  )
  
  if (is.na(lp)){lp=-Inf}
  
  # return negative log posterior
  return(-1*lp)
}


##==============================================================================
## End
##==============================================================================