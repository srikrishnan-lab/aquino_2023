rm(list=ls())                        # Clear all previous variables

## Set up MCMC stuff here so that it can be automated for HPC
nnode_mcmc000 <- 1
niter_mcmc000 <- 4e7 # 2e7 for full inversion, 4e7 for expert-only inversion
gamma_mcmc000 <- 0.51    # rate of adaptation (between 0.5 and 1, lower is faster adaptation)
accept_mcmc000 <- 0.15  # changed to 0.15 to help the sampler explore the space more carefully
                        # from Tony's 0.234: "Optimal as # parameters->infinity (Gelman et al, 1996; Roberts et al, 1997)" (Wong et al., 2017)
experts=FALSE # do you wish to invert expert assessment?
alldata=FALSE # do you wish to invert paleo and instrumental data?

if(alldata){
  if(experts){
  configure <- '_Complete_0921_2e7_' # configure for your run
  }
  else{
    configure <- '_Standard_0714_2e7_' # configure for your run
  }
} else{
  if (experts){
    configure <- '_Expert_0115_4e7_' # configure for your run
  }
  else{
    configure <- '_Priors_0426_4e7_'
  }
}
print(configure)
## Show plots? (probably want FALSE on HPC, non-interactive)
l.doplots <- FALSE

## Set up a filename for saving RData images along the way
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.saveprogress <- paste('CARL_calib_MCMC',configure,today,'.RData',sep='')

## Use the AIS fast dynamics emulator (Wong et al., 2017)
l.aisfastdy <- TRUE
fd.priors <- 'gamma'			# which prior distirbutions to use?

## Set the seed (for reproducibility)
set.seed(1234)
#set.seed(as.double(Sys.time())) # should yield same distributions... (a good test!)

## Needed libraries
library(adaptMCMC)
library(DEoptim)
# library(ncdf4)
library(compiler)
library(invgamma)
library(lhs)
enableJIT(3)

## Do you want to use RCP8.5 to make projections? (l.project=TRUE)
## Or do you want to use historical data to make hindcasts? (l.project=FALSE)
## Note -- l.project=FALSE => Kriegler (2005) data, same as Urban and Keller (2010)
l.project = TRUE
#begyear = 1765  # SNEASY start date
begyear = 1850  # DOECLIM start date
endyear = 2100
tstep   = 1
mod.time= seq(from=begyear, to=endyear, by=tstep)
begyear.norm <- 1986
endyear.norm <- 2005
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

## Source the models
source('../fortran/R/sneasyF.R')        # the SNEASY model (includes DOECLIM, and carbon cycle)
source('../fortran/R/doeclimF.R')       # the DOECLIM model
source('../fortran/R/GSIC_magiccF.R')   # the GSIC model
source('../fortran/R/brick_te_F.R')     # TE (thermosteric expansion) model
source('../fortran/R/brick_tee_F.R')    # TEE (explicit thermosteric expansion) model
source('../fortran/R/simpleF.R')        # GIS (Greenland Ice Sheet) model
source('../fortran/R/daisanto_fastdynF.R') # DAIS (Antarctic Ice Sheet) model
source('../R/brick_lws.R')              # LWS (land water storage)

## Source some useful functions for manipulating data
source('../R/forcing_total_carl.R')          # function to add up the total forcing
# source('../R/forcing_total.R')
source('../R/compute_indices.R')        # function to determine the model and
# data indices for comparisons

##==============================================================================
## BRICK:

## Read the model calibration data sets

## TODO
## TODO -- revise SNEASY_readData.R to match other components
## TODO

#source('../calibration/SNEASY_readData.R')    # read SNEASY calibration data
source('../calibration/DOECLIM_readData.R')   # read DOECLIM calibration data
source('../calibration/GSIC_readData.R')      # read GSIC calibration data
source('../calibration/TE_readData.R')        # read TE data
source('../calibration/SIMPLE_readData.R')    # GIS data, and trends in mass balance
#source('../calibration/DAIS_readData.R')     # DAIS forcing data (if at all uncoupled)


## TODO
## TODO -- add SNEASY calibration data to midx, oidx, obs, obs.err,
## TODO -- ind.norm.data, and i0
## TODO

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.sl)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "sl" )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
  c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
  c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
  c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Set the indices of the initial condition for each sub-model
i0 = vector("list",nrow(ind.norm.data)); names(i0)=as.character(ind.norm.data[,1])

## GSIC initial conditions are actually relative to 1990 (Wigley and Raper 2005)
## Re-set these. The simulation is relative to 1990, but results and comparison
## to data is relative to 1960.
i0$gsic = which(mod.time==1990)

## GIS initial conditions are relative to 1961-1990
i0$gis = which(mod.time==1961)

## DAIS: 

## Read the data forcing for hindcasts and projections. Yields:
##  Ta (Antarctic temperature reduced to sea-level)
##  Toc (High latitude subsurface ocean temperature)
##  SL, obs.sl (Reconstructed sea level, Church and White (2011) modern-era sea level)
##  dSL (time rate of change of sea level)
##  date (240 Kyr before present to 2100 AD at 1-year intervals of forcings)

source('../calibration/DAIS_readData.R')

## For prior distributions:
alpha.var <- 2     # alpha parameter for inverse gamma for var (E[x]=beta/(alpha+1))
beta.var <- 1      # beta parameter for inverse gamma for var (uncertainty parameter)
# note that the upper bound on var.dais is not actually imposed; just for plotting
shape.lambda <- 8.1              # gives 5% quantile at lambda=0.005 and
rate.lambda <- 100*shape.lambda  # gives mean at 0.01 m/yr, DeConto and Pollard (2016)
rate.Tcrit <- 1.37               # gives 5% quantile at Tcrit = -10 deg C
shape.Tcrit <- 15*rate.Tcrit     # gives mean at -15 deg C (negative requires multiplication of Tcrit by -1)

source('../fortran/R/daisanto_fastdynF.R')

##==============================================================================
## Which model(s) will you use?
## If you want to plug in your own model, insert the relevant "luse.XXX" line
## below, as well as into the "luse.brick = ..." command.
## Exactly one of te or tee must be TRUE.
luse.sneasy   = FALSE    # Simple Nonlinear EArth SYstem model (DOECLIM+CCM)
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermosteric expansion contribution to SLR
luse.tee      = FALSE   # explicit thermosteric expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.lws      = FALSE    # land water storage
luse.brick = cbind(luse.sneasy, luse.doeclim, luse.gsic, luse.te, luse.tee,
                   luse.simple, luse.dais, luse.lws)

## If you are using DAIS, include the fast dynamics emulator?
if (luse.dais) {
  l.aisfastdy <- TRUE
} else {
  l.aisfastdy <- FALSE  # force FALSE if not using DAIS
  slope.Ta2Tg <- NULL
  intercept.Ta2Tg <- NULL
}

if(luse.te & luse.tee) {
  luse.tee <- FALSE
  print('Only use 1 thermosteric expansion model; switching off explicit model.')
}

##==============================================================================
## Define parameters and their prior ranges
## -> Note: 'parnames' is defined here, which establishes how the parameters
##    are passed around into DEoptim, MCMC, likelihood functions, and the models

source('../calibration/CARL_BRICK_parameterSetup_v2.R')
# source('../calibration/CARL_BRICK_parameterSetup_v2_RCP85.R')


##==============================================================================
## Get the forcing data

if(luse.sneasy) {
  
  ## SNEASY
  setup.sneasy() # call this to initialize SNEASY
  
  ## TODO
  ## TODO -- (YG, TW) likely will need modified, to make sure forcing is as we want it
  ## TODO
  
  forcing <- vector('list',3)
  names(forcing) <- c('co2','aero','other')
  
  # load emissions data time series.
  rcp8.5.emis = read.csv("../data/sneasy_tmp/RCP85_EMISSIONS.csv")
  ibeg = which(rcp8.5.emis[,1]==begyear)
  iend = which(rcp8.5.emis[,1]==endyear)
  if(length(iend)*length(ibeg)==0) {print('ERROR - emissions data does not span mod.time')}
  emis = rcp8.5.emis[ibeg:iend,2] + rcp8.5.emis[ibeg:iend,3] # fossil + land use
  emisdata = data.frame(cbind(mod.time, emis))
  colnames(emisdata) = c("year","co2")
  forcing$co2 = emisdata$co2
  
  forcing.dat = read.table("../data/sneasy_tmp/forcing_rcp85.txt", header=TRUE)
  ibeg = which(forcing.dat[,1]==begyear)
  iend = which(forcing.dat[,1]==endyear)
  if(length(iend)*length(ibeg)==0) {print('ERROR - other radiative forcing data does not span mod.time')}
  aero.tmp  <- forcing.dat$aerosol.direct + forcing.dat$aerosol.indirect # aerosol radiative forcing
  other.tmp <- forcing.dat$ghg.nonco2 + forcing.dat$solar + forcing.dat$volcanic + forcing.dat$other
  forcing$aero  <- aero.tmp[ibeg:iend]
  forcing$other <- other.tmp[ibeg:iend]
  
} else {
  
  if(l.project) {
    forcing = read.csv( '../data/forcing_high6000.csv', header=TRUE )
    # forcing = read.csv( '../data/forcing_rcp85.csv', header=TRUE )
  } else {
    forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )
  }
  
}

##==============================================================================
## Define the coupled model
## -> need it defined before DEoptim, so we can calculate the objective function

source('../R/CARL_BRICK_coupledModel.R')
# source('../R/BRICK_coupledModel.R')


##==============================================================================

## Get a standard hindcast and future projections (a "control" model) for DAIS
if(l.aisfastdy) {
  Tcrit0 <- p0[match('Tcrit',parnames)]
  lambda0 <- p0[match('lambda',parnames)]
} else {
  Tcrit0 <- NULL
  lambda0 <- NULL
}

AIS_melt <- daisanto_fastdynF(
  anto.a = p0[match("anto.a",parnames)],
  anto.b = p0[match("anto.b",parnames)],
  gamma  = p0[match("gamma",parnames)],
  alpha  = p0[match("alpha.dais",parnames)],
  mu     = p0[match("mu",parnames)],
  nu     = p0[match("nu",parnames)],
  P0     = p0[match("P0",parnames)],
  kappa  = p0[match("kappa.dais",parnames)],
  f0     = p0[match("f0",parnames)],
  h0     = p0[match("h0",parnames)],
  c      = p0[match("c",parnames)],
  b0     = p0[match("b0",parnames)],
  slope  = p0[match("slope",parnames)],
  Tcrit  = Tcrit0,
  lambda = lambda0,
  l.aisfastdy = l.aisfastdy,
  Tg     = Tg.recon,
  slope.Ta2Tg = slope.Ta2Tg,
  intercept.Ta2Tg = intercept.Ta2Tg,
  SL     = SL,
  dSL    = dSL)

##==============================================================================
## DAIS: Set up (pre-)calibration windows around the data
## The point: We run the model at many parameter values and see which ones
##            send the simulation through a window (width determined by the
##            observational errors) around the observational data.
## These windows are presented in Shaffer (2014) and Shepherd et al. (2012)
## 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
## We want the cumulative sea-level equivalent in meters for the year 2002
## Note the conversion: 360Gt = 1mm SLE
## A fifth window is added to match IPCC AR5 Ch13 (page 1151) AIS SLR trend:
## 0.27 +/- 0.11 mm/year (convert to m/year here)

## Precal windows 5-?:
## Last "precalibration window" is 1993-2010 mean trend, from the IPCC AR5 Ch13
## (Page 1151), for AIS SLR contribution: 0.27 +- 0.11 mm/year
## Note that model output is in meters SLE and these trends are mm, so a
## conversion is necessary.

trends.ais <- c(0.27 , 0.08 , 0.40 )/1000   # m/year (after the /1000)
trends.err <- c(0.11 , 0.185, 0.205)/1000   # m/year (after the /1000)
trends.2up <- trends.ais+2*trends.err
trends.2dn <- trends.ais-2*trends.err
ind.trends <- mat.or.vec( length(trends.ais), 2)
ind.trends[1,] <- c(which(date==-7) , which(date==10)) # 1993-2010
ind.trends[2,] <- c(which(date==-8) , which(date== 1)) # 1992-2001
ind.trends[3,] <- c(which(date== 2) , which(date==11)) # 2002-2011

## Precal window 4:
## Adding observational constraint
estimate.SLE.rate <- abs(-71/360)/1000
time.years <- 2002-1992      # using the midpoint of the 19-year interval
mid.cum.SLE_2002 <- estimate.SLE.rate*time.years
i1992 <- which(date==-8)

estimate.SLE.rate.error <- abs(-53/360)/1000     #1-sigma error
estimate.SLE.error <- sqrt(time.years)*estimate.SLE.rate.error #1-sigma error
# (*sqrt(10) because 10 years of potentially accumulated error:
#  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
#                = 10*year X error^2)
SE2_2002 <- estimate.SLE.error*2 #2-sigma error

positive_2SE <- mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
negative_2SE <- mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

## Precal windows 1-3:
## from Shaffer (2014). modified by Kelsey
upper.wind <- c(7.4, -6.9, -1.25, positive_2SE) # Windows 2-3 from Kelsey, Window 1 from DeConto and Pollard 2016
lower.wind <- c(3.6, -15.8, -4.0, negative_2SE)
#upper.wind <- c(6.0, -6.9, -1.25, positive_2SE) # Windows 1-3 from Kelsey
#lower.wind <- c(1.8, -15.8, -4.0, negative_2SE)
#upper.wind <- c(5.5, -8 , -2, positive_2SE) # Windows 1-3 fFrom Shaffer 2014, p 1809
#lower.wind <- c(2.5, -17, -4, negative_2SE)

windows <- matrix(c(lower.wind, upper.wind), nrow = length(upper.wind), ncol=2)
obs.targets <- (windows[,2]+windows[,1])*.5  # middle of window = obs to compare model to
obs.err.dais <- (windows[,2]-windows[,1])*.5      # half-width of window = uncertainty
obs.err.dais <- 0.5*obs.err.dais                       # assume all windows are 2*stdErr (last two actually are)

## Create a vector with each observation year
## 120kyr, 20Kyr, 6kyr (before present), 2002, and 1993 (first year of the IPCC trend)
obs.years <- c(120000, 220000, 234000, 240002)

##==============================================================================


##==============================================================================
## Use differential optimization (DEoptim) algorithm to find suitable initial
## parameters for the MCMC chains

##TODO
##TODO -- TW, YG -- need to incorporate SNEASY CO2 (and MOC data?) into optimization
##TODO -- Also note that some of the SNEASY parameter bounds are infinite, which
##TODO -- causes problems in optimization.
##TODO
if(FALSE){ # pull starting parameters from a previous run; fix this later
  source('../calibration/CARL_BRICK_DEoptim.R')
  # source('../calibration/BRICK_DEoptim.R')
  p0.deoptim=p0                          # initialize optimized initial parameters
  niter.deoptim= 100                      # number of iterations for DE optimization
  NP.deoptim=11*length(index.model)      # population size for DEoptim (do at least 10*[N parameters])
  F.deoptim=0.8                          # as suggested by Storn et al (2006)
  CR.deoptim=0.9                        # as suggested by Storn et al (2006)
  outDEoptim <- DEoptim(minimize_residuals_brick, bound.lower[index.model], bound.upper[index.model],
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames.in=parnames[index.model], forcing.in=forcing        , l.project=l.project      ,
                        slope.Ta2Tg.in=slope.Ta2Tg       , intercept.Ta2Tg.in=intercept.Ta2Tg,
                        ind.norm.data=ind.norm.data      , ind.norm.sl=ind.norm      , mod.time=mod.time        ,
                        tstep=tstep                      , oidx = oidx.all           , midx = midx.all          ,
                        obs=obs.all                      , obs.err = obs.err.all     , trends.te = trends.te    ,
                        luse.brick = luse.brick           , i0 = i0                   , l.aisfastdy = l.aisfastdy )
  p0.deoptim[index.model] = outDEoptim$optim$bestmem
  
  p0.start <- p0.deoptim
}
# DAIS values from Ben Vega Westhoff # NVM, fix this later
if(FALSE){
p0.start[14] = 0.268933326
p0.start[15] = 0.013722995
p0.start[16] = 3.431541681
p0.start[17] = 0.297776878
p0.start[18] = 10.51496029
p0.start[19] = 0.008511003
p0.start[20] = 0.674262404
p0.start[21] = 0.057917442
p0.start[22] = 1.728871822
p0.start[23] = 1839.953735
p0.start[24] = 97.35321808
p0.start[25] = 801.059021
p0.start[26] = 0.000619136
p0.start[27] = 0.009839232
p0.start[28] = -15.8764534
p0.start[33] = 0.201209769
}

# for now pull starting values from a past calibration
amcmc.out<-readRDS("CARL_calib_amcmc_out_Standard_0714_extend_to_4e7_18Dec2021.rds",)
p0.start <- amcmc.out$samples[4e7,]
rm(amcmc.out)
## Run the model and examine output at these parameter values
brick.out = brick_model(
                        # parameters.in=p0.deoptim,
                        parameters.in=p0.start,
                        parnames.in=parnames,
                        forcing.in=forcing,
                        l.project=l.project,
                        #slope.Ta2Tg.in=slope.Ta2Tg,
                        #intercept.Ta2Tg.in=intercept.Ta2Tg,
                        mod.time=mod.time,
                        tstep=tstep,
                        ind.norm.data = ind.norm.data,
                        ind.norm.sl = ind.norm,
                        luse.brick = luse.brick,
                        i0 = i0,
                        l.aisfastdy = l.aisfastdy)

# p0.doplots = p0.deoptim
p0.doplots = p0.start

##TODO -- TW -- modify plotting for only the components that are used (incl SNEASY)
if(l.doplots) {
  if(TRUE){
  par(mfrow=c(3,2))
  # plot 1 -- DOECLIM, temperature match
  plot(obs.temp.time[oidx.temp], obs.temp[oidx.temp], pch=20, ylab='surface temperature anomaly [deg C]', xlab='year')
  lines(brick.out$doeclim.out$time[midx.temp], brick.out$doeclim.out$temp[midx.temp]+p0.start[4], col='red', lwd=2)
  # plot 2 -- DOECLIM, ocean heat match
  plot(obs.ocheat.time[oidx.ocheat], obs.ocheat[oidx.ocheat], pch=20, ylab='ocean heat uptake [10^22 J]', xlab='year')
  lines(brick.out$doeclim.out$time[midx.ocheat], brick.out$doeclim.out$ocheat[midx.ocheat]+p0.start[5], col='red', lwd=2)
  # plot 3 -- GSIC match
  plot(obs.gsic.time, obs.gsic, pch=20, ylab = 'sea-level equivalence [m]', xlab='year')
  lines(obs.gsic.time, brick.out$gsic.out[midx.gsic], col="blue", lwd=2)
  # plot 4 -- GIS match
  plot(obs.gis.time[oidx.gis], obs.gis[oidx.gis], pch=20, ylab = 'GIS mass balance change (SLE [m])', xlab='year')
  lines(obs.gis.time[oidx.gis], brick.out$simple.out$sle.gis[midx.gis], col="blue", lwd=2)
  # plot 5 -- SLR match
  plot(obs.sl.time[oidx.sl], obs.sl[oidx.sl], pch=20, ylab = 'sea level [m]', xlab='year')
  lines(brick.out$doeclim.out$time[midx.sl], brick.out$slr.out[midx.sl], col="purple", lwd=2)
  }
}


##==============================================================================
## Establish a gamma prior for drawing 1/tau.
## gamma distribution is the conjugate prior for the uncertain parameter beta
## in exponentially distributed random variable V(t):
##    dV/dt = -beta*V => V(t) = exp(-beta*t)
## For the TE model, this is not quite accurate, but similar enough (and
## previous testing with uniform priors on 1/tau confirm) that the gamma is an
## excellent approximation of the distribution of 1/tau (or tau)

if(luse.te) {
  # Original code did optimization to determine these parameters.
  # This, and future versions, will hand-pick.
  shape.invtau <- 1.81
  scale.invtau <- 0.00275
}

##==============================================================================
## Now set up the coupled model calibration
##  -- will need to whip up a coupled model likelihood function file
##   -- will calibrate with log.post(...) (in the above LL file) calling 'brick_model'
## Notes:
##  -- With DOECLIM+GSIC+TE+SIMPLE, all R models, requires ~1900s for 1e6 iterations.
##==============================================================================
## Set up and run the MCMC calibration

##TODO
##TODO -- TW, YG -- need to incorporate SNEASY CO2 (and MOC data?) into calibration
##TODO

## Source the statistical models
source('../calibration/CARL_assimLikelihood_priors.R')

## MCMC calibration
niter.mcmc <- niter_mcmc000       # number of iterations for MCMC
nnode.mcmc <- nnode_mcmc000        # number of nodes for parallel MCMC
accept.mcmc <- accept_mcmc000
gamma.mcmc <- gamma_mcmc000
stopadapt.mcmc <- round(niter.mcmc*0.1) # stop adapting after ?? iterations? (niter*1 => don't stop)

if(FALSE){ # manually step through likelihood for debugging
  if(TRUE){
    parameters.in = p0.start
    parnames.in=parnames
    forcing.in=forcing
    bound.lower.in = bound.lower
    bound.upper.in = bound.upper
    l.project=TRUE
    shape.in=shape.invtau       
    scale.in=scale.invtau
    slope.Ta2Tg.in=slope.Ta2Tg
    intercept.Ta2Tg.in=intercept.Ta2Tg    
    mod.time=mod.time
    ind.norm.data=ind.norm.data
    ind.norm.sl=ind.norm
    oidx = oidx.all                
    midx = midx.all
    obs=obs.all 
    obs.err = obs.err.all
    trends.te = trends.te
    luse.brick=luse.brick 
    i0=i0
    l.aisfastdy=l.aisfastdy
    
    obs.in=obs.targets 
    obs.err.in=obs.err.dais
    obs.step.in=obs.years
    trends.ais.in=trends.ais 
    trends.err.in=trends.err 
    ind.trends.in=ind.trends
    ind.norm.in = ind.relative
    alpha.var=alpha.var
    beta.var=beta.var
    SL.in=SL
    dSL.in=dSL
    Tg.in=Tg.recon
    
    rho.simple.in=NULL
    sigma.simple.in=NULL
  }
}

##==============================================================================
## Actually run the calibration

if(nnode.mcmc == 1) {
  t.beg <- proc.time()
  amcmc.out <- MCMC(log.post, niter.mcmc, p0.start, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
                    gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
                    
                    # BRICK:
                    parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                    slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg,
                    ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                    oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                    obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                    bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                    luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                    
                    # DAIS:
                    obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                    trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends  ,
                    ind.norm=ind.relative          , alpha.var=alpha.var        , beta.var=beta.var         ,
                    Tg.in=Tg.recon                 ,
                    shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit   ,
                    rate.Tcrit=rate.Tcrit          , 
                    #Ta.in=Ta  , Toc.in=Toc ,
                    SL.in=SL  , dSL.in=dSL                 , 
                    
                    # Configure
                    experts=experts                , alldata=alldata
                    
                    )
  
  t.end <- proc.time()
  chain1 <- amcmc.out$samples
  
  
} else if(nnode.mcmc > 1) {
  t.beg <- proc.time()
  amcmc.out <- MCMC.parallel(log.post, niter.mcmc, p0.start, n.chain=nnode.mcmc, n.cpu=nnode.mcmc,
                             dyn.libs=c('../fortran/doeclim.so','../fortran/brick_te.so','../fortran/brick_tee.so','../fortran/gsic_magicc.so','../fortran/simple.so','../fortran/dais_fastdyn.so'),
                             scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc, gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
                             
                             # BRICK:
                             parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                             slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
                             ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                             oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                             obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                             bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                             luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                             
                             # DAIS:
                             obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                             trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
                             ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
                             #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
                             Tg.in=Tg.recon                 ,
                             shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
                             rate.Tcrit=rate.Tcrit          , 
                             #l.aisfastdy=l.aisfastdy, # already accounted for above 
                             #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
                             SL.in=SL                       , dSL.in=dSL                 , 
                             
                             # Configure
                             experts=experts                , alldata=alldata
                             
                             )
  t.end <- proc.time()
  
  chain1 <- amcmc.out[[1]][1]$samples
}

if(luse.sneasy) {cleanup.sneasy()}  # deallocates memory after SNEASY is done

## Extend and run more MCMC samples?
if(FALSE){
  t.beg <- proc.time()
  amcmc.extend <- MCMC.add.samples(amcmc.out, niter.mcmc,
                                   
                                   # BRICK:
                                   parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                                   slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
                                   ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                                   oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                                   obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                                   bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                                   luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                                   
                                   # DAIS:
                                   obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                                   trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
                                   ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
                                   #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
                                   Tg.in=Tg.recon                 ,
                                   shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
                                   rate.Tcrit=rate.Tcrit          , 
                                   #l.aisfastdy=l.aisfastdy, # already accounted for above 
                                   #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
                                   SL.in=SL                       , dSL.in=dSL                 , 
                                   
                                   # Configure
                                   experts=experts                , alldata=alldata
                                   
                                   )
  t.end <- proc.time()
  chain1 <- amcmc.extend$samples
}

filename.RDS <- paste('CARL_calib_amcmc_out',configure,today,'.rds',sep='')
saveRDS(amcmc.out,file <- filename.RDS)
# saveRDS(time.elapsed,file <- "CALIB_amcmc_out_0415_1e7.rds")

##############
# Make some trace plots
parameters.posterior <- chain1
sequence1 = c(26,23,24,22,20,14,15,17,19,21,16,18,25,27,28)
sequence2 = c(5,6,7,8)
sequence3 = c(9,10,11,12,13)
sequence4 = c(1,2,3,4)
sequence5 = c(29,30,31,32,33)

# Trace plots for Antarctic parameters
if(TRUE){
  print("beginning trace plot 1")
  filename.trace1<- paste('../output_calibration/TracePlot',configure,'Set01','.png',sep='')
  png(filename=filename.trace1, width=1920, height=1080, units ="px")
  par(mfrow=c(3,5))
  for (pp in sequence1) {
    plot(parameters.posterior[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="", cex.lab=2.0)
    abline(v=stopadapt.mcmc,col="red",lwd=3)
  }
  dev.off()
} # End of if true

# Trace Plots for the rest of the parameters.
if(TRUE){
  print("beginning trace plot 2")
  filename.trace2<- paste('../output_calibration/TracePlot',configure,'Set02','.png',sep='')
  png(filename=filename.trace2, width=1920, height=1080, units ="px")
  par(mfrow=c(3,6))
  for (pp in sequence2){
    
    plot(parameters.posterior[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="", cex.lab=2.0)
    abline(v=stopadapt.mcmc,col="red",lwd=3)
  }
  
  for (pp in sequence3){
    plot(parameters.posterior[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="", cex.lab=2.0)
    abline(v=stopadapt.mcmc,col="red",lwd=3)
  }
  
  for (pp in sequence4){
    plot(parameters.posterior[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="", cex.lab=2.0)
    abline(v=stopadapt.mcmc,col="red",lwd=3)
  }
  
  for (pp in sequence5){
    plot(parameters.posterior[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="", cex.lab=2.0)
    abline(v=stopadapt.mcmc,col="red",lwd=3)
  }
  dev.off()
} # End of if true

# Individual Trace Plots
if(FALSE){
  print("beginning individual trace plots")
  for (pp in 1:length(parnames)){
    print(pp)
    filename.individual = paste('../output_calibration/individual_plots/TraceIndividual',configure,pp,'.png',sep='')
    png(filename=filename.individual, width=1920, height=1080, units ="px")
    
    plot(parameters.posterior[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
    abline(v=stopadapt.mcmc,col="red",lwd=3)
    
    dev.off()
  } 
}

## Save workspace image
t.end <- proc.time()
time.elapsed <- t.end - t.beg
save.image(file=filename.saveprogress)
print("line 705")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")
    
