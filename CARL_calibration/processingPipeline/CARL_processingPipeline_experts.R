rm(list=ls())

t.beg <- proc.time()

##==============================================================================
## IMPORTANT SETTINGS FOR YOU TO MODIFY
##              |
##              |
##              V
##==============================================================================

## Define the files you want to process/read/create
amcmc.out <- readRDS(file="CARL_calib_amcmc_out_Expert_0426_4e7_26Apr2022.rds",refhook=NULL) # load your calibration file here

# Set up a filename for saving RData images along the way
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
configure <- '_Expert_0426_4e7_' 
filename.saveprogress <- paste('CARL_pipeline_hindcasts_projections_',configure,today,'.RData',sep='')

print(configure)

appen <- ''        ## Append file name? In case you process multiple files in one day

l.aisfastdy <- TRUE        # including AIS fast dynamics in the DAIS version used? (must be consistent with how DAIS_calib_driver.R was run)
chain.length <- nrow(amcmc.out$samples)      # total proposed ensemble before burn-in
burn <- chain.length*0.1 # burn in amount
n.ensemble <- 1e6 #1e6 # total proposed ensemble before rejection sampling
                  # 1e6 for 1 core with 200GB RAM 
n.ensemble.final <- 40e3 # after rejection sampling, subsample to this total number

plotdir <- '../plots/'

# How many DAIS paleo simulations in your sample? (If you do a huge BRICK ensemble,
# it may be infeasible to do that number of 240,000 time step long simulations for
# DAIS paleo runs.)
n.samplepaleo <- chain.length

## Add heteroscedastic (observational) error into the hindcasts?
l.ar1.hetero <- TRUE

## Fingerprints of sea-level rise sources on local sea-level rise, if you want
## obtain local sea-level projections
lat.fp <- 29.95      # latitude of location to fingerprint local sea level rise (>0 is North, <0 is South)
lon.fp <- -90.07      # longitude of location ... (>0 is East, <0 is West)

## Mean and standard deviation for sampling land water storage (LWS)
## contributions to GMSL
## -- Using IPCC AR5 (Church et al. 2013) --
#lws.mean <- 0.38           # mm/y
#lws.sd   <- (0.49-.26)/4   # mm/y (take the IPCC 5-95% range as +/-2sigma)
## -- Using Dieng et al 2015 (doi:10.1088/1748-9326/10/12/124010)--
lws.mean <- 0.30           # mm/y
lws.sd   <- 0.18           # mm/y

##==============================================================================
##==============================================================================
##
##  NOTHING BELOW HERE SHOULD NEED MODIFIED
##
##==============================================================================
##==============================================================================

##==============================================================================
##==============================================================================
# simulate stationary AR(1) process (approximate - faster, better convergence, and
# results not sensitive to use of this as opposed to exact AR1)
ar1.sim <- function(N,rho1,sigma) {
  x <- rep(NA,N)
  if(length(sigma) > 1) {
    x[1] = rnorm(n=1,sd=sigma[1]/sqrt(1-rho1^2))
    for (i in 2:N) {
      x[i] <- rho1*x[i-1] + rnorm(1,sd=sigma[i])
    }
  } else {
    x[1] = rnorm(n=1,sd=sigma/sqrt(1-rho1^2))
    for (i in 2:N) {
      x[i] <- rho1*x[i-1] + rnorm(1,sd=sigma)
    }
  }
  return(x)
}
##==============================================================================
##==============================================================================

##==============================================================================
##==============================================================================
## Make posterior parameter draws to run an ensemble of simulations and make
## projections of SLR (include DAIS, from calibrated parameters)

# burn-in
chain.burned <- amcmc.out$samples[burn:chain.length,] 

# subsample
parameters <- matrix(0, n.ensemble, ncol(chain.burned))
sampleme <- seq(from=1, to=nrow(chain.burned),by=1)
sample.index <- sample(sampleme, n.ensemble, replace=TRUE)

for (i in 1:length(sample.index)){
  k = sample.index[i]
  parameters[i,] <- chain.burned[k,]
}

## Set up the model for hindcasts
l.project <- FALSE
begyear <- 1850
endyear <- 2009
mod.time <- begyear:endyear
begyear.norm <- 1986
endyear.norm <- 2005
ind.norm <- which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time <- length(mod.time)

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

## Initialize matrix to store model ensemble output
brick.out <- vector("list", n.ensemble)

## Initialize flag for possibly bad runs (DAIS+BRICK parameters could go wrong,
## because the other model components were calibrated without DAIS, and vice
## versa)
badruns <- rep(0, n.ensemble)

## Run the sample, and enjoy a nice progress bar
print(paste('Starting ',n.ensemble,' model hindcasts...',sep=''))
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {
  brick.out[[i]] <- brick_model(parameters.in     = as.numeric(parameters[i,]),
                                parnames.in       = parnames,
                                forcing.in        = forcing,
                                l.project         = l.project,
                                slope.Ta2Tg.in    = slope.Ta2Tg,
                                intercept.Ta2Tg.in= intercept.Ta2Tg,
                                mod.time          = mod.time,
                                ind.norm.data     = ind.norm.data,
                                ind.norm.sl       = ind.norm,
                                luse.brick        = luse.brick,
                                i0                = i0,
                                l.aisfastdy       = l.aisfastdy)
  setTxtProgressBar(pb, i)
}
close(pb)
print(paste(' ... done running model hindcasts'))

## Before rejection sampling to global mean sea level data, need to add the
## modeled statistical noise back in. Only using sea-level rise data, so only
## need to modify GSIC, GIS. Using the statistical parameters for AR1, AR1 and
## Gaussian noise, respectively. Do not do for AIS, because var.dais was fit to
## paleo data-model mismatch, not representative of the current era.

## Gather the fields for each simulation (easy referencing for plotting and
## analysis)
slr.out    <- mat.or.vec(n.ensemble,length(mod.time))
temp.out   <- mat.or.vec(n.ensemble,length(mod.time))
ocheat.out <- mat.or.vec(n.ensemble,length(mod.time))
gsic.out   <- mat.or.vec(n.ensemble,length(mod.time))
gis.out    <- mat.or.vec(n.ensemble,length(mod.time))
ais.out    <- mat.or.vec(n.ensemble,length(mod.time))
ais.disint.out <- mat.or.vec(n.ensemble,length(mod.time))
te.out     <- mat.or.vec(n.ensemble,length(mod.time))

## Will also normalize the output to "ind.norm" (1961-1990? 1986-2005? (Mengel, IPCC))
slr.out.norm    <- slr.out
temp.out.norm   <- temp.out
ocheat.out.norm <- ocheat.out
gsic.out.norm   <- gsic.out
gis.out.norm    <- gis.out
ais.out.norm    <- ais.out
te.out.norm     <- te.out

## And add statistical noise
slr.norm.stat    <- slr.out
temp.norm.stat   <- temp.out
ocheat.norm.stat <- ocheat.out
gsic.norm.stat   <- gsic.out
gis.norm.stat    <- gis.out
ais.norm.stat    <- ais.out
te.norm.stat     <- te.out

# Need to grab fixed DOECLIM parameters since we bypassed the forcings
## DOECLIM (Urban and Keller, 2010, values
parnames.doeclim   =NULL; p0.doeclim       =NULL; bound.lower.doeclim=NULL;
bound.upper.doeclim=NULL; step.mcmc.doeclim=NULL; index.model.doeclim=NULL;
if (luse.doeclim) {
  parnames.doeclim   =c("S" ,"kappa.doeclim","alpha.doeclim","T0"  ,"H0" ,"sigma.T","sigma.H","rho.T","rho.H")	# parameters names
  p0.doeclim		     =c(3.1 , 3.5   , 1.1           , -0.06, -33 , 0.1     , 2       , 0.55  , 0.9   )	# initial parameter guesses
  bound.lower.doeclim=c(0.1 , 0.1   , 0             , -0.3 , -50 , 0.05    , 0.1     , 0     , 0     )	# prior range lower bounds
  bound.upper.doeclim=c(10  , 4     , 2             ,  0.3 ,   0 , 5       , 10      , 0.999 , 0.999 )	# prior range upper bounds
  step.mcmc.doeclim	 =c(0.16, 0.17  ,0.025          ,0.003 , 0.9 , 5e-4    , 0.025   , 0.007 , 0.006 )	# step size for parameters in MCMC (proposals)
  index.model.doeclim=c(1,2,3,4,5)		# which are model parameters? (index within parnames.doeclim)
}

## Go through each simulation and collect, normalize and add modeled error to
## the fields

print(paste('Starting to add up total sea level rise from model hindcasts...',sep=''))

pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {
  
  T0 <- p0.doeclim[match("T0"           ,parnames.doeclim)]
  H0 <- p0.doeclim[match("H0"           ,parnames.doeclim)]
  
  # set the results (note that we already have slr.out)
  temp.out[i,]   <- brick.out[[i]]$doeclim.out$temp + T0
  ocheat.out[i,] <- brick.out[[i]]$doeclim.out$ocheat + H0
  
  # Normalize the output to "ind.norm.data"
  temp.out.norm[i,]   <- temp.out[i,]  -mean(temp.out[i,1:20])
  ocheat.out.norm[i,] <- ocheat.out[i,]#-mean(ocheat.out[i,ind.norm.data[which(ind.norm.data[,1]=='ocheat'),2]:ind.norm.data[which(ind.norm.data[,1]=='ocheat'),3]])

  # # Add the statistcal model for AR1 (or otherwise) noise
  sigma.T    <- p0.doeclim[match("sigma.T"           ,parnames.doeclim)]
  rho.T      <- p0.doeclim[match("rho.T"           ,parnames.doeclim)]
  sigma.H    <- p0.doeclim[match("sigma.H"           ,parnames.doeclim)]
  rho.H      <- p0.doeclim[match("rho.H"           ,parnames.doeclim)]
  err.temp   <- rep(sigma.T,n.time); if(l.ar1.hetero) {err.temp[midx.temp]=sqrt(sigma.T^2 + obs.temp.err[oidx.temp]^2)}
  err.ocheat <- rep(sigma.H,n.time); if(l.ar1.hetero) {err.ocheat[midx.ocheat]=sqrt(sigma.H^2+obs.ocheat.err[oidx.ocheat]^2)}
  temp.norm.stat[i,]   <- temp.out.norm[i,]   + ar1.sim(n.time, rho.T, err.temp)
  ocheat.norm.stat[i,] <- ocheat.out.norm[i,] + ar1.sim(n.time, rho.H, err.ocheat)

  gsic.out[i,]   <- brick.out[[i]]$gsic.out
  gis.out[i,]    <- brick.out[[i]]$simple.out$sle.gis
  ais.out[i,]    <- brick.out[[i]]$dais.out$Vais
  ais.disint.out[i,] <- brick.out[[i]]$dais.out$Vdisint
  te.out[i,]     <- brick.out[[i]]$te.out
  
  slr.out[i,]    <- brick.out[[i]]$slr.out
  
  # Normalize the output to "ind.norm.data"
  gsic.out.norm[i,] <- gsic.out[i,]  -mean(gsic.out[i,ind.norm.data[which(ind.norm.data[,1]=='gsic'),2]:ind.norm.data[which(ind.norm.data[,1]=='gsic'),3]])
  gis.out.norm[i,]  <- gis.out[i,]   -mean(gis.out[i,ind.norm.data[which(ind.norm.data[,1]=='gis'),2]:ind.norm.data[which(ind.norm.data[,1]=='gis'),3]])
  ais.out.norm[i,]  <- ais.out[i,]   -mean(ais.out[i,ind.norm.data[which(ind.norm.data[,1]=='ais'),2]:ind.norm.data[which(ind.norm.data[,1]=='ais'),3]])
  te.out.norm[i,]   <- te.out[i,]    -mean(te.out[i,ind.norm.data[which(ind.norm.data[,1]=='te'),2]:ind.norm.data[which(ind.norm.data[,1]=='te'),3]])

  # Add the statistical model for AR1 (or otherwise) noise
  sigma.gsic   <- parameters[i,match("sigma.gsic"  ,parnames)]
  rho.gsic     <- parameters[i,match("rho.gsic"    ,parnames)]
  sigma.simple <- parameters[i,match("sigma.simple",parnames)]
  rho.simple   <- parameters[i,match("rho.simple"  ,parnames)]
  var.dais     <- parameters[i,match("var.dais"    ,parnames)]

  err.gsic <- rep(sigma.gsic,n.time); if(l.ar1.hetero) {err.gsic[midx.gsic] <- sqrt(sigma.gsic^2+obs.gsic.err[oidx.gsic]^2)}
  err.gis  <- rep(sigma.simple,n.time); if(l.ar1.hetero) {err.gis[midx.gis] <- sqrt(sigma.simple^2+obs.gis.err^2)}

  gsic.norm.stat[i,] <- gsic.out.norm[i,] + ar1.sim(n.time, rho.gsic, err.gsic)
  gis.norm.stat[i,]  <- gis.out.norm[i,]  + ar1.sim(n.time, rho.simple, err.gis)
  ais.norm.stat[i,]  <- ais.out.norm[i,] #+ rnorm(  n.time, mean=0,sd=sqrt(var.dais))
  te.norm.stat[i,]   <- te.out.norm[i,]

  slr.norm.stat[i,] <- gsic.out.norm[i,] +
    gis.out.norm[i,]  +
    ais.out.norm[i,]  +
    te.out.norm[i,]

  slr.norm.stat[i,] <- slr.norm.stat[i,] - mean(slr.norm.stat[i,ind.norm.data[which(ind.norm.data[,1]=='sl'),2]:ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])


  
  setTxtProgressBar(pb, i)
}
close(pb)
print(paste(' ... done adding up model hindcast sea level rise and contributions'))

# Rejection sampling, with target distribution as the likelihood of the sea-
# level rise data, proposing uniformly across the ensemble members.
survive <- rep(0, n.ensemble)

## Make sure SLR data are also normalized
ibeg <- which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),2]])
iend <- which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])
obs.sl <- obs.sl - mean(obs.sl[ibeg:iend])

# calibrate to the Church and White data with land water subtracted out
# and uncertainties added in quadrature
# assumed budget: TE+AIS+GIS+GSIC+LWS = GMSL
# 1901-1990: –0.11 [–0.16 to –0.06] (5-95% range)
lw.time.1900 <- 1900:1989
i1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
lw.1900 <- (-0.11/1000)*(lw.time.1900 - 1900)
lw.err.1900 <- (0.25*(-0.06--0.16)/1000)*sqrt(lw.time.1900 - lw.time.1900[1])
# 1971-2010: 0.12 [0.03 to 0.22]
lw.time.1970 <- 1970:2009
i1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
lw.1970 <- (0.12/1000)*(lw.time.1970 - lw.time.1970[1])
lw.err.1970 <- (0.25*(0.22-0.03)/1000)*sqrt(lw.time.1970 - lw.time.1970[1])
# 1993-2010: 0.38 [0.26 to 0.49]
lw.time.1992 <- 1992:2009
i1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])
lw.1992 <- (0.38/1000)*(lw.time.1992 - lw.time.1992[1])
lw.err.1992 <- (0.25*(0.49-0.26)/1000)*sqrt(lw.time.1992 - lw.time.1992[1])

# normalize, subtract and add error in quadrature
obs.sl.lw.1900 <- obs.sl[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])] - obs.sl[which(obs.sl.time==lw.time.1900[1])]
obs.sl.lw.1970 <- obs.sl[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])] - obs.sl[which(obs.sl.time==lw.time.1970[1])]
obs.sl.lw.1992 <- obs.sl[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])] - obs.sl[which(obs.sl.time==lw.time.1992[1])]

obs.sl.lw.1900 <- obs.sl.lw.1900 - lw.1900
obs.sl.lw.1970 <- obs.sl.lw.1970 - lw.1970
obs.sl.lw.1992 <- obs.sl.lw.1992 - lw.1992

obs.sl.lw.err.1900 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])]^2 + lw.err.1900^2)
obs.sl.lw.err.1970 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])]^2 + lw.err.1970^2)
obs.sl.lw.err.1992 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])]^2 + lw.err.1992^2)

# calculate likelihood as the product of the three independent likelihoods
resid.1900 <- obs.sl.lw.1900 - obs.sl.lw.1900
llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
resid.1970 <- obs.sl.lw.1970 - obs.sl.lw.1970
llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
resid.1992 <- obs.sl.lw.1992 - obs.sl.lw.1992
llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
lik.max <- (llik.1900 + llik.1970 + llik.1992)/n.ensemble

imod.1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
imod.1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
imod.1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])

uni.rnd <- log(runif(n.ensemble))
for (i in 1:n.ensemble) {
  resid.1900 <- obs.sl.lw.1900 - (slr.norm.stat[i,imod.1900]-slr.norm.stat[i,imod.1900[1]])
  resid.1970 <- obs.sl.lw.1970 - (slr.norm.stat[i,imod.1970]-slr.norm.stat[i,imod.1970[1]])
  resid.1992 <- obs.sl.lw.1992 - (slr.norm.stat[i,imod.1992]-slr.norm.stat[i,imod.1992[1]])
  llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
  llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
  llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
  lik.mem <- llik.1900 + llik.1970 + llik.1992
  if( uni.rnd[i] <= lik.mem-lik.max) {survive[i]=1}
}
ind.survive <- which( as.logical(survive))
print(paste('Calibration to sea level data by rejection sampling leaves ',length(ind.survive),' full calibrated ensemble members',sep=''))

slr.out.good <- slr.norm.stat[ind.survive,]
parameters.good <- parameters[ind.survive,]
colnames(parameters.good) <- parnames

# subsample from post-rejection sample down to n.ensemble.final (40,000 for 1 core, 200 GB RAM)
if (nrow(parameters.good)>n.ensemble.final){
  
  # subsample
  parameters.subsample <- matrix(0, n.ensemble.final, ncol(parameters.good))
  sampleme <- seq(from=1, to=nrow(parameters.good),by=1)
  sample.index <- sample(sampleme, n.ensemble.final, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters.subsample[i,] <- parameters.good[k,]
  }
  
  parameters.good <- parameters.subsample
  
}

saveRDS(parameters.good, file="experts_0426_parameters_good.rds")

## Save progress so far via workspace image
# save.image(file=filename.saveprogress)

##==============================================================================
##==============================================================================

## DAIS paleo runs with the post-calibrated parameters
## Already have dSL, SL, Tg.recon.

## How many members do you want in your ensemble?
n.ensemble <- nrow(parameters.good)
n.parameters <- ncol(parameters.good)
n.sample <- min(n.samplepaleo, n.ensemble)
ind.sample <- sample( 1:n.ensemble, size=n.sample, replace=FALSE)
parameters.sample <- parameters.good[ind.sample,]

n.paleo <- length(SL)
dais.paleo <- mat.or.vec(n.sample, n.paleo)
date <- seq(-239999,16,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
norm.period <- c(1961,1990)
ibeg <- which(date==(norm.period[1]-2000))
iend <- which(date==(norm.period[2]-2000))
ind.norm.paleo <- ibeg:iend
t.paleo <- date

print(paste('Starting ',n.sample,' DAIS paleo hindcast simulations, using the post-calibrated parameters ...',sep=''))
pb <- txtProgressBar(min=0,max=n.sample,initial=0,style=3);
for (i in 1:n.sample) {
  
  anto.a <- parameters.sample[i,match("anto.a",parnames)]
  anto.b <- parameters.sample[i,match("anto.b",parnames)]
  gamma  <- parameters.sample[i,match("gamma" ,parnames)]
  alpha  <- parameters.sample[i,match("alpha.dais" ,parnames)]
  mu     <- parameters.sample[i,match("mu" ,parnames)]
  nu     <- parameters.sample[i,match("nu" ,parnames)]
  P0     <- parameters.sample[i,match("P0" ,parnames)]
  kappa  <- parameters.sample[i,match("kappa.dais" ,parnames)]
  f0     <- parameters.sample[i,match("f0" ,parnames)]
  h0     <- parameters.sample[i,match("h0" ,parnames)]
  c      <- parameters.sample[i,match("c" ,parnames)]
  b0     <- parameters.sample[i,match("b0" ,parnames)]
  slope  <- parameters.sample[i,match("slope" ,parnames)]
  var.dais <- parameters.sample[i,match("var.dais" ,parnames)]
  if(l.aisfastdy) {
    Tcrit  <- parameters.sample[i,match("Tcrit" ,parnames)]
    lambda <- parameters.sample[i,match("lambda" ,parnames)]
  } else {
    Tcrit  <- NULL
    lambda <- NULL
  }
  
  dais.tmp <- daisanto_fastdynF(
    anto.a=anto.a, anto.b=anto.b,
    gamma=gamma  , alpha=alpha  ,
    mu=mu        , nu=nu        ,
    P0=P0        , kappa=kappa  ,
    f0=f0        , h0=h0        ,
    c=c          , b0=b0        ,
    slope=slope  , l.aisfastdy=l.aisfastdy,
    Tcrit=Tcrit  , lambda=lambda,
    slope.Ta2Tg=slope.Ta2Tg, intercept.Ta2Tg=intercept.Ta2Tg,
    Tg=Tg.recon  , SL=SL ,
    dSL=dSL      , includes_dSLais=1
  )
  
  # Subtract off the 1961-1990 normalization period
  dais.norm <- dais.tmp$Vais - mean(dais.tmp$Vais[ind.norm.paleo])
  
  # Add the modeled error back in
  dais.paleo[i,] <- dais.norm + rnorm(n.paleo, mean=0,sd=sqrt(var.dais))
  
  setTxtProgressBar(pb, i)
}
close(pb)

## Save progress so far via workspace image
# save.image(file=filename.saveprogress)

## Get 5-95% CI for hindcasts
## Initialize arrays for the output
dais.paleo.05  <- rep(NA,length(date)); dais.paleo.50  <- rep(NA,length(date)); dais.paleo.95 <- rep(NA,length(date))
dais.paleo.max <- rep(NA,length(date)); dais.paleo.min <- rep(NA,length(date));

## Source a useful script, to allow for assigning multiple outputs at once
source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
pb <- txtProgressBar(min=0,max=length(date),initial=0,style=3);
for (t in 1:length(date)){
  c(dais.paleo.05[t] , dais.paleo.50[t] , dais.paleo.95[t] , dais.paleo.max[t], dais.paleo.min[t]) := quantile(dais.paleo[,t],c(0.05,.50,.95,1,0), na.rm=TRUE)
  setTxtProgressBar(pb, t)
}
close(pb)
print(paste(' ... done with the paleo simulations!',sep=''))

## Make smoothed version of the AIS paleo results, so the plots are not massive
n.avg <- 100  # number of years in averaging period
n.time.avg <- ceiling(length(date)/n.avg)
dais.paleo.05.avg <- rep(NA,n.time.avg)
dais.paleo.50.avg <- rep(NA,n.time.avg)
dais.paleo.95.avg <- rep(NA,n.time.avg)
dais.paleo.max.avg <- rep(NA,n.time.avg)
dais.paleo.min.avg <- rep(NA,n.time.avg)
date.avg <- seq(date[1],date[length(date)],by=n.avg)

print(paste('Smoothing the paleo simulations with ',n.avg,'-year averages ...',sep=''))
pb <- txtProgressBar(min=0,max=n.time.avg,initial=0,style=3);
for (t in 1:(n.time.avg-1)){
  dais.paleo.05.avg[t] <- mean(dais.paleo.05[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.50.avg[t] <- mean(dais.paleo.50[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.95.avg[t] <- mean(dais.paleo.95[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.max.avg[t] <- mean(dais.paleo.max[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.min.avg[t] <- mean(dais.paleo.min[((t-1)*n.avg+1) : (t*n.avg)])
  setTxtProgressBar(pb, t)
}
dais.paleo.05.avg[n.time.avg] <- mean(dais.paleo.05[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.50.avg[n.time.avg] <- mean(dais.paleo.50[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.95.avg[n.time.avg] <- mean(dais.paleo.95[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.max.avg[n.time.avg] <- mean(dais.paleo.max[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.min.avg[n.time.avg] <- mean(dais.paleo.min[((n.time.avg-1)*n.avg+1) : length(date)])
close(pb)
print(paste(' ... done smoothing the paleo simulations!',sep=''))

##==============================================================================
##==============================================================================

## Save the hindcasts, trimming down first for ind.survive (Church and White, 2011)
## and then for AIS Vmin
gsic.hind <- t(gsic.norm.stat[ind.survive,])
te.hind <- t(te.norm.stat[ind.survive,])
gis.hind <- t(gis.norm.stat[ind.survive,])
ais.hind <- t(ais.norm.stat[ind.survive,])
temp.hind <- t(temp.norm.stat[ind.survive,])
ocheat.hind <- t(ocheat.norm.stat[ind.survive,])
gsl.hind <- t(slr.out.good)
t.hind <- mod.time

saveRDS(t.hind, file="experts_0426_t_hind.rds")

## Save progress so far via workspace image
# save.image(file=filename.saveprogress)

#END HINDCASTS
##==============================================================================
##==============================================================================










##==============================================================================
##==============================================================================

## Make High Scenario projections to 2100. Need:
## (1) global total sea level, (2) global sea level without fast dynamics,
## (3) local sea level, (4) local sea level without fast dynamics
n.ensemble <- nrow(parameters.good)      # ... or take all of them
n.parameters <- ncol(parameters.good)    # number of distinct model parameters

## Display output to let the user know what is going on, and how many ensemble
## members will be processed into projections.
print(paste('Starting projections for ensemble of ',n.ensemble,' members...',sep=''))

## Set up the model projections
l.project <- TRUE
begyear <- 1850
endyear <- 2100; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
mod.time <- begyear:endyear
t.proj <- mod.time # save for output

## NOTE: IPCC generally are relative to 1986-2005.
## Therefore use that period to be commensurate with their results
begyear.norm <- 2000
endyear.norm <- 2000
ind.norm <- which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time <- length(mod.time)

## Want model projections under RCP2.6, 4.5 and 8.5. Initialize a list, each
## will hold the ensemble model output for a different forcing scenario.
names.scen = c('highT')
n.scen <- length(names.scen)
proj.out <- vector("list", n.scen)
forc.scen <- vector("list", n.scen)
badruns.scen <- vector("list", n.scen)
names(proj.out) <- names.scen

## Get the forcing data
forc.scen[[1]] <- read.csv( '../data/forcing_high6000.csv', header=TRUE )

## Loop over forcing scenarios
for (ff in 1:n.scen) {
  
  forcing = forc.scen[[ff]]
  
  ## Initialize matrix to store model ensemble output and flag for bad runs
  ## (could go bad because of different forcing from hindcasts; different
  ## temperatures lead to different stability requirements for SIMPLE)
  brick.out <- vector("list", n.ensemble)
  badruns   <- rep(0, n.ensemble)
  
  ## Run the sample, and enjoy a nice progress bar
  pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3);
  for (i in 1:n.ensemble) {
    
    brick.out[[i]] <- brick_model(parameters.in     = parameters.good[i,],
                                  parnames.in       = parnames,
                                  forcing.in        = forcing,
                                  l.project         = l.project,
                                  slope.Ta2Tg.in    = slope.Ta2Tg,
                                  intercept.Ta2Tg.in= intercept.Ta2Tg,
                                  mod.time          = mod.time,
                                  ind.norm.data     = ind.norm.data,
                                  ind.norm.sl       = ind.norm,
                                  luse.brick        = luse.brick,
                                  i0                = i0,
                                  l.aisfastdy       = l.aisfastdy)
    
    # check if the run turned out bad
    if( is.na(brick.out[[i]]$slr.out[length(mod.time)]) ) {badruns[i]=1}
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  badruns.scen[[ff]] <- badruns
  proj.out[[ff]] <- brick.out
}

## Save progress so far via workspace image
# save.image(file=filename.saveprogress)

## Filter out any bad runs
ind.good <- which(badruns.scen[[1]]==0)
n.ensemble <- length(ind.good)
parameters <- parameters.good[ind.good,]

## Trim down the hindcasts to match (some runs might go off the rails with RCP
## forcings but didn't with the hindcast forcings)
gsl.hind <- gsl.hind[,ind.good]
gsic.hind <- gsic.hind[,ind.good]
te.hind <- te.hind[,ind.good]
gis.hind <- gis.hind[,ind.good]
ais.hind <- ais.hind[,ind.good]
temp.hind <- temp.hind[,ind.good]
ocheat.hind <- ocheat.hind[,ind.good]

brick.rcp26 <- proj.out[[1]][ind.good]
# brick.rcp45 <- proj.out[[2]][ind.good]
# brick.rcp60 <- proj.out[[3]][ind.good]
# brick.rcp85 <- proj.out[[4]][ind.good]

## Gather the fields for each simulation (easy referencing for plotting and
## analysis)
proj.rcp26 <- vector("list", n.scen)
# proj.rcp45 <- vector("list", n.scen)
# proj.rcp60 <- vector("list", n.scen)
# proj.rcp85 <- vector("list", n.scen)

## Make each projections list a list of the SLR contributions and other fields
## Get the 90% CI at each time step, of all model output fields

names.output <- c('slr','gsic','gis','ais','disint','te','lws','temp','ocheat','slr.local')
n.output <- length(names.output)

proj.rcp26 <- vector("list", n.output)
# proj.rcp45 <- vector("list", n.output)
# proj.rcp60 <- vector("list", n.output)
# proj.rcp85 <- vector("list", n.output)

names(proj.rcp26) <- names.output
# names(proj.rcp45) <- names.output
# names(proj.rcp60) <- names.output
# names(proj.rcp85) <- names.output

for (j in 1:n.output) {
  proj.rcp26[[j]] <- mat.or.vec(n.ensemble, n.time)
  # proj.rcp45[[j]] <- mat.or.vec(n.ensemble, n.time)
  # proj.rcp60[[j]] <- mat.or.vec(n.ensemble, n.time)
  # proj.rcp85[[j]] <- mat.or.vec(n.ensemble, n.time)
}

## Grab the DOECLIM parameters
parnames.doeclim   =c("S" ,"kappa.doeclim","alpha.doeclim","T0"  ,"H0" ,"sigma.T","sigma.H","rho.T","rho.H")	# parameters names
p0.doeclim         =c(3.1 , 3.5   , 1.1           , -0.06, -33 , 0.1     , 2       , 0.55  , 0.9   )	# initial parameter guesses
S            =p0.doeclim[match("S"            ,parnames.doeclim)]

## Go through each simulation and collect, normalize and add modeled error to
## the fields

pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)

for (i in 1:n.ensemble) {
  
  T0 <- p0.doeclim[match("T0"            ,parnames.doeclim)]
  H0 <- p0.doeclim[match("H0"            ,parnames.doeclim)]
  
  # set the results (note that we already have slr.out)
  proj.rcp26$temp[i,] <- brick.rcp26[[i]]$doeclim.out$temp + T0
  proj.rcp26$ocheat[i,] <- brick.rcp26[[i]]$doeclim.out$ocheat + H0
  
  # proj.rcp45$temp[i,] <- brick.rcp45[[i]]$doeclim.out$temp + T0
  # proj.rcp45$ocheat[i,] <- brick.rcp45[[i]]$doeclim.out$ocheat + H0
  # 
  # proj.rcp60$temp[i,] <- brick.rcp60[[i]]$doeclim.out$temp + T0
  # proj.rcp60$ocheat[i,] <- brick.rcp60[[i]]$doeclim.out$ocheat + H0
  # 
  # proj.rcp85$temp[i,] <- brick.rcp85[[i]]$doeclim.out$temp + T0
  # proj.rcp85$ocheat[i,] <- brick.rcp85[[i]]$doeclim.out$ocheat + H0
  
  # Normalize the output to "ind.norm".
  # Normalize ocean heat uptake too, for sake of plotting (it will be plotted as
  # the amount of heat taken up by ocean since 1986-2005 period)
  proj.rcp26$temp[i,] <- proj.rcp26$temp[i,] - mean( proj.rcp26$temp[i,ind.norm])
  proj.rcp26$ocheat[i,] <- proj.rcp26$ocheat[i,] - mean( proj.rcp26$ocheat[i,ind.norm])
  
  # proj.rcp45$temp[i,] <- proj.rcp45$temp[i,] - mean( proj.rcp45$temp[i,ind.norm])
  # proj.rcp45$ocheat[i,] <- proj.rcp45$ocheat[i,] - mean( proj.rcp45$ocheat[i,ind.norm])
  # 
  # proj.rcp60$temp[i,] <- proj.rcp60$temp[i,] - mean( proj.rcp60$temp[i,ind.norm])
  # proj.rcp60$ocheat[i,] <- proj.rcp60$ocheat[i,] - mean( proj.rcp60$ocheat[i,ind.norm])
  # 
  # proj.rcp85$temp[i,] <- proj.rcp85$temp[i,] - mean( proj.rcp85$temp[i,ind.norm])
  # proj.rcp85$ocheat[i,] <- proj.rcp85$ocheat[i,] - mean( proj.rcp85$ocheat[i,ind.norm])
  
  # Add the statistcal model for AR1 (or otherwise) noise?
  # Edit: not for DAIS, because that noise matches paleo data, which has
  # different uncertainties than modern data.
  sigma.T  <- p0.doeclim[match("sigma.T"            ,parnames.doeclim)]
  rho.T    <- p0.doeclim[match("rho.T"            ,parnames.doeclim)]
  sigma.H  <- p0.doeclim[match("sigma.H"            ,parnames.doeclim)]
  rho.H    <- p0.doeclim[match("rho.H"            ,parnames.doeclim)]
  
  proj.rcp26$temp[i,] <- proj.rcp26$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
  proj.rcp26$ocheat[i,] <- proj.rcp26$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)
  
  # proj.rcp45$temp[i,] <- proj.rcp45$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
  # proj.rcp45$ocheat[i,] <- proj.rcp45$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)
  # 
  # proj.rcp60$temp[i,] <- proj.rcp60$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
  # proj.rcp60$ocheat[i,] <- proj.rcp60$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)
  # 
  # proj.rcp85$temp[i,] <- proj.rcp85$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
  # proj.rcp85$ocheat[i,] <- proj.rcp85$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)
  
  # set the results (note that we already have slr.out)
  proj.rcp26$gsic[i,] <- brick.rcp26[[i]]$gsic.out
  proj.rcp26$gis[i,] <- brick.rcp26[[i]]$simple.out$sle.gis
  proj.rcp26$ais[i,] <- brick.rcp26[[i]]$dais.out$Vais
  proj.rcp26$disint[i,] <- brick.rcp26[[i]]$dais.out$Vdisint
  proj.rcp26$te[i,] <- brick.rcp26[[i]]$te.out
  
  # proj.rcp45$gsic[i,] <- brick.rcp45[[i]]$gsic.out
  # proj.rcp45$gis[i,] <- brick.rcp45[[i]]$simple.out$sle.gis
  # proj.rcp45$ais[i,] <- brick.rcp45[[i]]$dais.out$Vais
  # proj.rcp45$disint[i,] <- brick.rcp45[[i]]$dais.out$Vdisint
  # proj.rcp45$te[i,] <- brick.rcp45[[i]]$te.out
  # 
  # proj.rcp60$gsic[i,] <- brick.rcp60[[i]]$gsic.out
  # proj.rcp60$gis[i,] <- brick.rcp60[[i]]$simple.out$sle.gis
  # proj.rcp60$ais[i,] <- brick.rcp60[[i]]$dais.out$Vais
  # proj.rcp60$disint[i,] <- brick.rcp60[[i]]$dais.out$Vdisint
  # proj.rcp60$te[i,] <- brick.rcp60[[i]]$te.out
  # 
  # proj.rcp85$gsic[i,] <- brick.rcp85[[i]]$gsic.out
  # proj.rcp85$gis[i,] <- brick.rcp85[[i]]$simple.out$sle.gis
  # proj.rcp85$ais[i,] <- brick.rcp85[[i]]$dais.out$Vais
  # proj.rcp85$disint[i,] <- brick.rcp85[[i]]$dais.out$Vdisint
  # proj.rcp85$te[i,] <- brick.rcp85[[i]]$te.out
  
  # Normalize the output to "ind.norm" (1961-1990? 1986-2005 (Mengel)?).
  # Normalize ocean heat uptake too, for sake of plotting (it will be plotted as
  # the amount of heat taken up by ocean since 1986-2005 period)
  proj.rcp26$gsic[i,] <- proj.rcp26$gsic[i,] - mean( proj.rcp26$gsic[i,ind.norm])
  proj.rcp26$gis[i,] <- proj.rcp26$gis[i,] - mean( proj.rcp26$gis[i,ind.norm])
  proj.rcp26$ais[i,] <- proj.rcp26$ais[i,] - mean( proj.rcp26$ais[i,ind.norm])
  proj.rcp26$te[i,] <- proj.rcp26$te[i,] - mean( proj.rcp26$te[i,ind.norm])
  
  # proj.rcp45$gsic[i,] <- proj.rcp45$gsic[i,] - mean( proj.rcp45$gsic[i,ind.norm])
  # proj.rcp45$gis[i,] <- proj.rcp45$gis[i,] - mean( proj.rcp45$gis[i,ind.norm])
  # proj.rcp45$ais[i,] <- proj.rcp45$ais[i,] - mean( proj.rcp45$ais[i,ind.norm])
  # proj.rcp45$te[i,] <- proj.rcp45$te[i,] - mean( proj.rcp45$te[i,ind.norm])
  # 
  # proj.rcp60$gsic[i,] <- proj.rcp60$gsic[i,] - mean( proj.rcp60$gsic[i,ind.norm])
  # proj.rcp60$gis[i,] <- proj.rcp60$gis[i,] - mean( proj.rcp60$gis[i,ind.norm])
  # proj.rcp60$ais[i,] <- proj.rcp60$ais[i,] - mean( proj.rcp60$ais[i,ind.norm])
  # proj.rcp60$te[i,] <- proj.rcp60$te[i,] - mean( proj.rcp60$te[i,ind.norm])
  # 
  # proj.rcp85$gsic[i,] <-proj.rcp85$gsic[i,] - mean( proj.rcp85$gsic[i,ind.norm])
  # proj.rcp85$gis[i,] <- proj.rcp85$gis[i,] - mean( proj.rcp85$gis[i,ind.norm])
  # proj.rcp85$ais[i,] <- proj.rcp85$ais[i,] - mean( proj.rcp85$ais[i,ind.norm])
  # proj.rcp85$te[i,] <- proj.rcp85$te[i,] - mean( proj.rcp85$te[i,ind.norm])
  # 
  # Add the statistcal model for AR1 (or otherwise) noise?
  # Edit: not for DAIS, because that noise matches paleo data, which has
  # different uncertainties than modern data.
  sigma.gsic   <- parameters[i,match("sigma.gsic"  ,parnames)]
  rho.gsic     <- parameters[i,match("rho.gsic"    ,parnames)]
  sigma.simple <- parameters[i,match("sigma.simple",parnames)]
  rho.simple   <- parameters[i,match("rho.simple"  ,parnames)]
  var.dais     <- parameters[i,match("var.dais"    ,parnames)]
  
  proj.rcp26$gsic[i,] <- proj.rcp26$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
  proj.rcp26$gis[i,] <- proj.rcp26$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)
  
  # proj.rcp45$gsic[i,] <- proj.rcp45$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
  # proj.rcp45$gis[i,] <- proj.rcp45$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)
  # 
  # proj.rcp60$gsic[i,] <- proj.rcp60$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
  # proj.rcp60$gis[i,] <- proj.rcp60$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)
  # 
  # proj.rcp85$gsic[i,] <- proj.rcp85$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
  # proj.rcp85$gis[i,] <- proj.rcp85$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)
  
  # Add contributions to land water storage. /1000 to convert to meters.
  # This is done for each ensemble member and each RCP.
  # All other contributions are normalized to 1986-2005 (some have the
  # estimated noise added, so mean(1986-2005) not nec. equall to 0), so
  # normalize the lws.est contribution to this period.
  
  # RCP2.6
  proj.rcp26$lws[i,] <- cumsum(rnorm(n=n.time, mean=lws.mean, sd=lws.sd)) /1000
  proj.rcp26$lws[i,] <- proj.rcp26$lws[i,] - mean(proj.rcp26$lws[i,ind.norm])
  
  # # RCP4.5
  # proj.rcp45$lws[i,] <- cumsum(rnorm(n=n.time, mean=lws.mean, sd=lws.sd)) /1000
  # proj.rcp45$lws[i,] <- proj.rcp45$lws[i,] - mean(proj.rcp45$lws[i,ind.norm])
  # 
  # # RCP6.0
  # proj.rcp60$lws[i,] <- cumsum(rnorm(n=n.time, mean=lws.mean, sd=lws.sd)) /1000
  # proj.rcp60$lws[i,] <- proj.rcp60$lws[i,] - mean(proj.rcp60$lws[i,ind.norm])
  # 
  # # RCP8.5
  # proj.rcp85$lws[i,] <- cumsum(rnorm(n=n.time, mean=lws.mean, sd=lws.sd)) /1000
  # proj.rcp85$lws[i,] <- proj.rcp85$lws[i,] - mean(proj.rcp85$lws[i,ind.norm])
  
  # Add up to total sea-level rise
  proj.rcp26$slr[i,] <- proj.rcp26$gsic[i,] +
    proj.rcp26$gis[i,]  +
    proj.rcp26$ais[i,]  +
    proj.rcp26$te[i,]   +
    proj.rcp26$lws[i,]
  
  # proj.rcp45$slr[i,] <- proj.rcp45$gsic[i,] +
  #   proj.rcp45$gis[i,]  +
  #   proj.rcp45$ais[i,]  +
  #   proj.rcp45$te[i,]   +
  #   proj.rcp45$lws[i,]
  # 
  # proj.rcp60$slr[i,] <- proj.rcp60$gsic[i,] +
  #   proj.rcp60$gis[i,]  +
  #   proj.rcp60$ais[i,]  +
  #   proj.rcp60$te[i,]   +
  #   proj.rcp60$lws[i,]
  # 
  # proj.rcp85$slr[i,] <- proj.rcp85$gsic[i,] +
  #   proj.rcp85$gis[i,]  +
  #   proj.rcp85$ais[i,]  +
  #   proj.rcp85$te[i,]   +
  #   proj.rcp85$lws[i,]
  
  # And normalize sea-level rise
  proj.rcp26$slr[i,] <- proj.rcp26$slr[i,] - mean(proj.rcp26$slr[i,ind.norm])
  # proj.rcp45$slr[i,] <- proj.rcp45$slr[i,] - mean(proj.rcp45$slr[i,ind.norm])
  # proj.rcp60$slr[i,] <- proj.rcp60$slr[i,] - mean(proj.rcp60$slr[i,ind.norm])
  # proj.rcp85$slr[i,] <- proj.rcp85$slr[i,] - mean(proj.rcp85$slr[i,ind.norm])
  
  setTxtProgressBar(pb, i)
}
close(pb)

print(paste('... finished projections!',sep=''))

##==============================================================================
##==============================================================================

# SAVE THE FILES!

# for hindcasts
saveRDS(dais.paleo.05.avg, file="experts_0426_ais_paleo_05.rds")
saveRDS(dais.paleo.50.avg, file="experts_0426_ais_paleo_50.rds")
saveRDS(dais.paleo.95.avg, file="experts_0426_ais_paleo_95.rds")
saveRDS(date.avg, file="experts_0426_t_paleo.rds")

saveRDS(gsic.hind, file="experts_0426_gsic_hind.rds")
saveRDS(te.hind, file="experts_0426_te_hind.rds")
saveRDS(gis.hind, file="experts_0426_gis_hind.rds")
saveRDS(ais.hind, file="experts_0426_ais_hind.rds")
saveRDS(temp.hind, file="experts_0426_temp_hind.rds")
saveRDS(ocheat.hind, file="experts_0426_ocheat_hind.rds")
saveRDS(gsl.hind, file="experts_0426_gsl_hind.rds")

# for projections
saveRDS(proj.rcp26$slr,file="experts_0426_slr_rcp26.rds")
saveRDS(proj.rcp26$te, file="experts_0426_te_rcp26.rds") 
saveRDS(proj.rcp26$gis, file="experts_0426_gis_rcp26.rds")
saveRDS(proj.rcp26$gsic, file="experts_0426_gsic_rcp26.rds")
saveRDS(proj.rcp26$ais, file="experts_0426_ais_rcp26.rds")
saveRDS(proj.rcp26$temp, file="experts_0426_temp_rcp26.rds")
saveRDS(proj.rcp26$ocheat, file="experts_0426_ocheat_rcp26.rds")
saveRDS(t.proj, file="experts_0426_t_proj.rds")

## Save workspace image
t.end <- proc.time()
time.elapsed <- t.end - t.beg
# save.image(file=filename.saveprogress)
print("line 1044")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")
