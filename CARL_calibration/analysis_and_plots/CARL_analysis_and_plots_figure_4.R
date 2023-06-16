rm(list=ls())
graphics.off()

t.beg <- proc.time()

## Load pipeline files
t.hind <- readRDS(file="priors_0426_t_hind.rds")
t.proj <- readRDS(file="priors_0426_t_proj.rds")
t.paleo <- readRDS(file="priors_0426_t_paleo.rds")

## NOTE: THERE IS ONLY ONE SCENARIO: "high temperature" from SEJ2018.
## here, "rcp26" refers to "high temperature" from SEJ2018.

#####################################################################
## Priors
priors.ais.paleo.05 <- readRDS(file="priors_0426_ais_paleo_05.rds")
priors.ais.paleo.50 <- readRDS(file="priors_0426_ais_paleo_50.rds")
priors.ais.paleo.95 <- readRDS(file="priors_0426_ais_paleo_95.rds")

priors.gsic.hind <- readRDS(file="priors_0426_gsic_hind.rds")
priors.te.hind <- readRDS(file="priors_0426_te_hind.rds")
priors.gis.hind <- readRDS(file="priors_0426_gis_hind.rds")
priors.ais.hind <- readRDS(file="priors_0426_ais_hind.rds")
priors.temp.hind <- readRDS(file="priors_0426_temp_hind.rds")
priors.ocheat.hind <- readRDS(file="priors_0426_ocheat_hind.rds")
priors.gsl.hind <- readRDS(file="priors_0426_gsl_hind.rds")

# for projections
priors.slr.rcp26 <- readRDS(file="priors_0426_slr_rcp26.rds")
priors.te.rcp26 <- readRDS(file="priors_0426_te_rcp26.rds") 
priors.gis.rcp26 <- readRDS(file="priors_0426_gis_rcp26.rds")
priors.gsic.rcp26 <- readRDS(file="priors_0426_gsic_rcp26.rds")
priors.ais.rcp26 <- readRDS(file="priors_0426_ais_rcp26.rds")
priors.temp.rcp26 <- readRDS(file="priors_0426_temp_rcp26.rds")
priors.ocheat.rcp26 <- readRDS(file="priors_0426_ocheat_rcp26.rds")
#####################################################################

#####################################################################
## Experts
experts.ais.paleo.05 <- readRDS(file="experts_0426_ais_paleo_05.rds")
experts.ais.paleo.50 <- readRDS(file="experts_0426_ais_paleo_50.rds")
experts.ais.paleo.95 <- readRDS(file="experts_0426_ais_paleo_95.rds")

experts.gsic.hind <- readRDS(file="experts_0426_gsic_hind.rds")
experts.te.hind <- readRDS(file="experts_0426_te_hind.rds")
experts.gis.hind <- readRDS(file="experts_0426_gis_hind.rds")
experts.ais.hind <- readRDS(file="experts_0426_ais_hind.rds")
experts.temp.hind <- readRDS(file="experts_0426_temp_hind.rds")
experts.ocheat.hind <- readRDS(file="experts_0426_ocheat_hind.rds")
experts.gsl.hind <- readRDS(file="experts_0426_gsl_hind.rds")

# for projections
experts.slr.rcp26 <- readRDS(file="experts_0426_slr_rcp26.rds")
experts.te.rcp26 <- readRDS(file="experts_0426_te_rcp26.rds") 
experts.gis.rcp26 <- readRDS(file="experts_0426_gis_rcp26.rds")
experts.gsic.rcp26 <- readRDS(file="experts_0426_gsic_rcp26.rds")
experts.ais.rcp26 <- readRDS(file="experts_0426_ais_rcp26.rds")
experts.temp.rcp26 <- readRDS(file="experts_0426_temp_rcp26.rds")
experts.ocheat.rcp26 <- readRDS(file="experts_0426_ocheat_rcp26.rds")
#####################################################################

#####################################################################
## Standard
standard.ais.paleo.05 <- readRDS(file="standard_0426_ais_paleo_05.rds")
standard.ais.paleo.50 <- readRDS(file="standard_0426_ais_paleo_50.rds")
standard.ais.paleo.95 <- readRDS(file="standard_0426_ais_paleo_95.rds")

standard.gsic.hind <- readRDS(file="standard_0426_gsic_hind.rds")
standard.te.hind <- readRDS(file="standard_0426_te_hind.rds")
standard.gis.hind <- readRDS(file="standard_0426_gis_hind.rds")
standard.ais.hind <- readRDS(file="standard_0426_ais_hind.rds")
standard.temp.hind <- readRDS(file="standard_0426_temp_hind.rds")
standard.ocheat.hind <- readRDS(file="standard_0426_ocheat_hind.rds")
standard.gsl.hind <- readRDS(file="standard_0426_gsl_hind.rds")

# for projections
standard.slr.rcp26 <- readRDS(file="standard_0426_slr_rcp26.rds")
standard.te.rcp26 <- readRDS(file="standard_0426_te_rcp26.rds") 
standard.gis.rcp26 <- readRDS(file="standard_0426_gis_rcp26.rds")
standard.gsic.rcp26 <- readRDS(file="standard_0426_gsic_rcp26.rds")
standard.ais.rcp26 <- readRDS(file="standard_0426_ais_rcp26.rds")
standard.temp.rcp26 <- readRDS(file="standard_0426_temp_rcp26.rds")
standard.ocheat.rcp26 <- readRDS(file="standard_0426_ocheat_rcp26.rds")
#####################################################################

#####################################################################
## Complete
complete.ais.paleo.05 <- readRDS(file="complete_0426_ais_paleo_05.rds")
complete.ais.paleo.50 <- readRDS(file="complete_0426_ais_paleo_50.rds")
complete.ais.paleo.95 <- readRDS(file="complete_0426_ais_paleo_95.rds")

complete.gsic.hind <- readRDS(file="complete_0426_gsic_hind.rds")
complete.te.hind <- readRDS(file="complete_0426_te_hind.rds")
complete.gis.hind <- readRDS(file="complete_0426_gis_hind.rds")
complete.ais.hind <- readRDS(file="complete_0426_ais_hind.rds")
complete.temp.hind <- readRDS(file="complete_0426_temp_hind.rds")
complete.ocheat.hind <- readRDS(file="complete_0426_ocheat_hind.rds")
complete.gsl.hind <- readRDS(file="complete_0426_gsl_hind.rds")

# for projections
complete.slr.rcp26 <- readRDS(file="complete_0426_slr_rcp26.rds")
complete.te.rcp26 <- readRDS(file="complete_0426_te_rcp26.rds") 
complete.gis.rcp26 <- readRDS(file="complete_0426_gis_rcp26.rds")
complete.gsic.rcp26 <- readRDS(file="complete_0426_gsic_rcp26.rds")
complete.ais.rcp26 <- readRDS(file="complete_0426_ais_rcp26.rds")
complete.temp.rcp26 <- readRDS(file="complete_0426_temp_rcp26.rds")
complete.ocheat.rcp26 <- readRDS(file="complete_0426_ocheat_rcp26.rds")
#####################################################################



## Other useful scripts
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array
source('../Useful/MultipleOutput.R')    # defines the useful ":=" operator

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

## Set colors to use for control model, and observations/experimental model,
## in all plots. This is indexed within "mycol", from "colorblindPalette.R".
colmod <- 2
colobs <- 11

##==============================================================================
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
##==============================================================================

##==============================================================================
##==============================================================================
## setup hindcasts and projections

## Initialize arrays for the output
# priors
priors.slr.05 = rep(NA,length(t.hind));		  priors.slr.50 = rep(NA,length(t.hind));	  	priors.slr.95 = rep(NA,length(t.hind))
priors.gsic.05 = rep(NA,length(t.hind));		priors.gsic.50 = rep(NA,length(t.hind));		priors.gsic.95 = rep(NA,length(t.hind))
priors.gis.05 = rep(NA,length(t.hind));		  priors.gis.50 = rep(NA,length(t.hind));		  priors.gis.95 = rep(NA,length(t.hind))
priors.te.05 = rep(NA,length(t.hind));			priors.te.50 = rep(NA,length(t.hind));			priors.te.95 = rep(NA,length(t.hind))
priors.temp.05 = rep(NA,length(t.hind));		priors.temp.50 = rep(NA,length(t.hind));		priors.temp.95 = rep(NA,length(t.hind))
priors.ocheat.05 = rep(NA,length(t.hind));  priors.ocheat.50 = rep(NA,length(t.hind));	priors.ocheat.95 = rep(NA,length(t.hind))

# experts
experts.slr.05 = rep(NA,length(t.hind));		experts.slr.50 = rep(NA,length(t.hind));	  experts.slr.95 = rep(NA,length(t.hind))
experts.gsic.05 = rep(NA,length(t.hind));		experts.gsic.50 = rep(NA,length(t.hind));		experts.gsic.95 = rep(NA,length(t.hind))
experts.gis.05 = rep(NA,length(t.hind));		experts.gis.50 = rep(NA,length(t.hind));		experts.gis.95 = rep(NA,length(t.hind))
experts.te.05 = rep(NA,length(t.hind));			experts.te.50 = rep(NA,length(t.hind));			experts.te.95 = rep(NA,length(t.hind))
experts.temp.05 = rep(NA,length(t.hind));		experts.temp.50 = rep(NA,length(t.hind));		experts.temp.95 = rep(NA,length(t.hind))
experts.ocheat.05 = rep(NA,length(t.hind)); experts.ocheat.50 = rep(NA,length(t.hind));	experts.ocheat.95 = rep(NA,length(t.hind))

# standard
standard.slr.05 = rep(NA,length(t.hind));		  standard.slr.50 = rep(NA,length(t.hind));	  	standard.slr.95 = rep(NA,length(t.hind))
standard.gsic.05 = rep(NA,length(t.hind));		standard.gsic.50 = rep(NA,length(t.hind));		standard.gsic.95 = rep(NA,length(t.hind))
standard.gis.05 = rep(NA,length(t.hind));		  standard.gis.50 = rep(NA,length(t.hind));		  standard.gis.95 = rep(NA,length(t.hind))
standard.te.05 = rep(NA,length(t.hind));			standard.te.50 = rep(NA,length(t.hind));			standard.te.95 = rep(NA,length(t.hind))
standard.temp.05 = rep(NA,length(t.hind));		standard.temp.50 = rep(NA,length(t.hind));		standard.temp.95 = rep(NA,length(t.hind))
standard.ocheat.05 = rep(NA,length(t.hind));  standard.ocheat.50 = rep(NA,length(t.hind));	standard.ocheat.95 = rep(NA,length(t.hind))

# complete
complete.slr.05 = rep(NA,length(t.hind));		  complete.slr.50 = rep(NA,length(t.hind));	  	complete.slr.95 = rep(NA,length(t.hind))
complete.gsic.05 = rep(NA,length(t.hind));		complete.gsic.50 = rep(NA,length(t.hind));		complete.gsic.95 = rep(NA,length(t.hind))
complete.gis.05 = rep(NA,length(t.hind));		  complete.gis.50 = rep(NA,length(t.hind));		  complete.gis.95 = rep(NA,length(t.hind))
complete.te.05 = rep(NA,length(t.hind));			complete.te.50 = rep(NA,length(t.hind));			complete.te.95 = rep(NA,length(t.hind))
complete.temp.05 = rep(NA,length(t.hind));		complete.temp.50 = rep(NA,length(t.hind));		complete.temp.95 = rep(NA,length(t.hind))
complete.ocheat.05 = rep(NA,length(t.hind));  complete.ocheat.50 = rep(NA,length(t.hind));	complete.ocheat.95 = rep(NA,length(t.hind))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.hind)){
  # priors  
  c(priors.slr.05[t],    priors.slr.50[t],    priors.slr.95[t])				:= quantile(priors.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gsic.05[t],   priors.gsic.50[t],   priors.gsic.95[t])			:= quantile(priors.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gis.05[t],    priors.gis.50[t],    priors.gis.95[t])				:= quantile(priors.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.te.05[t],     priors.te.50[t],     priors.te.95[t])				:= quantile(priors.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.temp.05[t],   priors.temp.50[t],   priors.temp.95[t])			:= quantile(priors.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.ocheat.05[t], priors.ocheat.50[t], priors.ocheat.95[t])	  := quantile(priors.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # experts 
  c(experts.slr.05[t],    experts.slr.50[t],    experts.slr.95[t])		:= quantile(experts.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gsic.05[t],   experts.gsic.50[t],   experts.gsic.95[t])		:= quantile(experts.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gis.05[t],    experts.gis.50[t],    experts.gis.95[t])		:= quantile(experts.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.te.05[t],     experts.te.50[t],     experts.te.95[t])			:= quantile(experts.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.temp.05[t],   experts.temp.50[t],   experts.temp.95[t])		:= quantile(experts.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.ocheat.05[t], experts.ocheat.50[t], experts.ocheat.95[t])	:= quantile(experts.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # standard  
  c(standard.slr.05[t],    standard.slr.50[t], standard.slr.95[t])				:= quantile(standard.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gsic.05[t],   standard.gsic.50[t], standard.gsic.95[t])			:= quantile(standard.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gis.05[t],    standard.gis.50[t], standard.gis.95[t])				:= quantile(standard.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.te.05[t],     standard.te.50[t], standard.te.95[t])					:= quantile(standard.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.temp.05[t],   standard.temp.50[t], standard.temp.95[t])			:= quantile(standard.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.ocheat.05[t], standard.ocheat.50[t], standard.ocheat.95[t])	:= quantile(standard.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # complete  
  c(complete.slr.05[t],    complete.slr.50[t], complete.slr.95[t])				:= quantile(complete.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gsic.05[t],   complete.gsic.50[t], complete.gsic.95[t])			:= quantile(complete.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gis.05[t],    complete.gis.50[t], complete.gis.95[t])				:= quantile(complete.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.te.05[t],     complete.te.50[t], complete.te.95[t])					:= quantile(complete.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.temp.05[t],   complete.temp.50[t], complete.temp.95[t])			:= quantile(complete.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.ocheat.05[t], complete.ocheat.50[t], complete.ocheat.95[t])	:= quantile(complete.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
}

begyear = t.hind[1]
endyear = t.hind[length(t.hind)]
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time = length(mod.time)

source('../R/compute_indices.R')  # function to determine the model and
# data indices for comparisons

## Source the data for hindcast comparisons
source('../calibration/DOECLIM_readData.R')
source('../calibration/GSIC_readData.R')
source('../calibration/SIMPLE_readData.R')
source('../calibration/DAIS_readData.R')
source('../calibration/TE_readData.R')

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

## Also normalize the available data to this time period.
ibeg=which(obs.temp.time==mod.time[1])
iend=which(obs.temp.time==mod.time[20])
obs.temp.norm = obs.temp - mean(obs.temp[ibeg:iend])

ibeg=which(obs.ocheat.time==mod.time[ind.norm[1]])
iend=which(obs.ocheat.time==mod.time[ind.norm[length(ind.norm)]])
obs.ocheat.norm = obs.ocheat - mean(obs.ocheat[ibeg:iend])

ibeg=which(obs.gsic.time==mod.time[ind.norm[1]])
iend=which(obs.gsic.time==mod.time[ind.norm[length(ind.norm)]])
obs.gsic.norm = obs.gsic #- mean(obs.gsic[ibeg:iend]) # GSIC does not need normalized - already is normalized to 1960

ibeg=which(obs.gis.time==mod.time[ind.norm[1]])
iend=which(obs.gis.time==mod.time[ind.norm[length(ind.norm)]])
obs.gis.norm = obs.gis - mean(obs.gis[ibeg:iend])

ibeg=which(obs.sl.time==mod.time[ind.norm[1]])
iend=which(obs.sl.time==mod.time[ind.norm[length(ind.norm)]])
obs.sl.norm = obs.sl - mean(obs.sl[ibeg:iend])

## TE trends
# read in TE_readData.R -- "trends.te"

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

##
## 5-95% CI of hindcasts, with obs in there
##

#transpose:

# priors
priors.slr.rcp26 <- t(priors.slr.rcp26)
priors.te.rcp26 <- t(priors.te.rcp26)
priors.gis.rcp26 <- t(priors.gis.rcp26)
priors.gsic.rcp26 <- t(priors.gsic.rcp26)
priors.ais.rcp26 <- t(priors.ais.rcp26)
priors.temp.rcp26 <- t(priors.temp.rcp26)
priors.ocheat.rcp26 <- t(priors.ocheat.rcp26)

# experts
experts.slr.rcp26 <- t(experts.slr.rcp26)
experts.te.rcp26 <- t(experts.te.rcp26)
experts.gis.rcp26 <- t(experts.gis.rcp26)
experts.gsic.rcp26 <- t(experts.gsic.rcp26)
experts.ais.rcp26 <- t(experts.ais.rcp26)
experts.temp.rcp26 <- t(experts.temp.rcp26)
experts.ocheat.rcp26 <- t(experts.ocheat.rcp26)


# standard
standard.slr.rcp26 <- t(standard.slr.rcp26)
standard.te.rcp26 <- t(standard.te.rcp26)
standard.gis.rcp26 <- t(standard.gis.rcp26)
standard.gsic.rcp26 <- t(standard.gsic.rcp26)
standard.ais.rcp26 <- t(standard.ais.rcp26)
standard.temp.rcp26 <- t(standard.temp.rcp26)
standard.ocheat.rcp26 <- t(standard.ocheat.rcp26)

# complete
complete.slr.rcp26 <- t(complete.slr.rcp26)
complete.te.rcp26 <- t(complete.te.rcp26)
complete.gis.rcp26 <- t(complete.gis.rcp26)
complete.gsic.rcp26 <- t(complete.gsic.rcp26)
complete.ais.rcp26 <- t(complete.ais.rcp26)
complete.temp.rcp26 <- t(complete.temp.rcp26)
complete.ocheat.rcp26 <- t(complete.ocheat.rcp26)


## Initialize arrays for the output
# priors
priors.slr.rcp26.05 = rep(NA,length(t.proj)); priors.slr.rcp26.50 = rep(NA,length(t.proj)); priors.slr.rcp26.95 = rep(NA,length(t.proj))
priors.ais.rcp26.05 = rep(NA,length(t.proj)); priors.ais.rcp26.50 = rep(NA,length(t.proj)); priors.ais.rcp26.95 = rep(NA,length(t.proj))
priors.gis.rcp26.05 = rep(NA,length(t.proj)); priors.gis.rcp26.50 = rep(NA,length(t.proj)); priors.gis.rcp26.95 = rep(NA,length(t.proj))
priors.gsic.rcp26.05 = rep(NA,length(t.proj)); priors.gsic.rcp26.50 = rep(NA,length(t.proj)); priors.gsic.rcp26.95 = rep(NA,length(t.proj))
priors.te.rcp26.05 = rep(NA,length(t.proj)); priors.te.rcp26.50 = rep(NA,length(t.proj)); priors.te.rcp26.95 = rep(NA,length(t.proj))
priors.temp.rcp26.05 = rep(NA,length(t.proj)); priors.temp.rcp26.50 = rep(NA,length(t.proj)); priors.temp.rcp26.95 = rep(NA,length(t.proj))
priors.ocheat.rcp26.05 = rep(NA,length(t.proj)); priors.ocheat.rcp26.50 = rep(NA,length(t.proj)); priors.ocheat.rcp26.95 = rep(NA,length(t.proj))

# experts
experts.slr.rcp26.05 = rep(NA,length(t.proj)); experts.slr.rcp26.50 = rep(NA,length(t.proj)); experts.slr.rcp26.95 = rep(NA,length(t.proj))
experts.ais.rcp26.05 = rep(NA,length(t.proj)); experts.ais.rcp26.50 = rep(NA,length(t.proj)); experts.ais.rcp26.95 = rep(NA,length(t.proj))
experts.gis.rcp26.05 = rep(NA,length(t.proj)); experts.gis.rcp26.50 = rep(NA,length(t.proj)); experts.gis.rcp26.95 = rep(NA,length(t.proj))
experts.gsic.rcp26.05 = rep(NA,length(t.proj)); experts.gsic.rcp26.50 = rep(NA,length(t.proj)); experts.gsic.rcp26.95 = rep(NA,length(t.proj))
experts.te.rcp26.05 = rep(NA,length(t.proj)); experts.te.rcp26.50 = rep(NA,length(t.proj)); experts.te.rcp26.95 = rep(NA,length(t.proj))
experts.temp.rcp26.05 = rep(NA,length(t.proj)); experts.temp.rcp26.50 = rep(NA,length(t.proj)); experts.temp.rcp26.95 = rep(NA,length(t.proj))
experts.ocheat.rcp26.05 = rep(NA,length(t.proj)); experts.ocheat.rcp26.50 = rep(NA,length(t.proj)); experts.ocheat.rcp26.95 = rep(NA,length(t.proj))

# standard
standard.slr.rcp26.05 = rep(NA,length(t.proj)); standard.slr.rcp26.50 = rep(NA,length(t.proj)); standard.slr.rcp26.95 = rep(NA,length(t.proj))
standard.ais.rcp26.05 = rep(NA,length(t.proj)); standard.ais.rcp26.50 = rep(NA,length(t.proj)); standard.ais.rcp26.95 = rep(NA,length(t.proj))
standard.gis.rcp26.05 = rep(NA,length(t.proj)); standard.gis.rcp26.50 = rep(NA,length(t.proj)); standard.gis.rcp26.95 = rep(NA,length(t.proj))
standard.gsic.rcp26.05 = rep(NA,length(t.proj)); standard.gsic.rcp26.50 = rep(NA,length(t.proj)); standard.gsic.rcp26.95 = rep(NA,length(t.proj))
standard.te.rcp26.05 = rep(NA,length(t.proj)); standard.te.rcp26.50 = rep(NA,length(t.proj)); standard.te.rcp26.95 = rep(NA,length(t.proj))
standard.temp.rcp26.05 = rep(NA,length(t.proj)); standard.temp.rcp26.50 = rep(NA,length(t.proj)); standard.temp.rcp26.95 = rep(NA,length(t.proj))
standard.ocheat.rcp26.05 = rep(NA,length(t.proj)); standard.ocheat.rcp26.50 = rep(NA,length(t.proj)); standard.ocheat.rcp26.95 = rep(NA,length(t.proj))

# complete
complete.slr.rcp26.05 = rep(NA,length(t.proj)); complete.slr.rcp26.50 = rep(NA,length(t.proj)); complete.slr.rcp26.95 = rep(NA,length(t.proj))
complete.ais.rcp26.05 = rep(NA,length(t.proj)); complete.ais.rcp26.50 = rep(NA,length(t.proj)); complete.ais.rcp26.95 = rep(NA,length(t.proj))
complete.gis.rcp26.05 = rep(NA,length(t.proj)); complete.gis.rcp26.50 = rep(NA,length(t.proj)); complete.gis.rcp26.95 = rep(NA,length(t.proj))
complete.gsic.rcp26.05 = rep(NA,length(t.proj)); complete.gsic.rcp26.50 = rep(NA,length(t.proj)); complete.gsic.rcp26.95 = rep(NA,length(t.proj))
complete.te.rcp26.05 = rep(NA,length(t.proj)); complete.te.rcp26.50 = rep(NA,length(t.proj)); complete.te.rcp26.95 = rep(NA,length(t.proj))
complete.temp.rcp26.05 = rep(NA,length(t.proj)); complete.temp.rcp26.50 = rep(NA,length(t.proj)); complete.temp.rcp26.95 = rep(NA,length(t.proj))
complete.ocheat.rcp26.05 = rep(NA,length(t.proj)); complete.ocheat.rcp26.50 = rep(NA,length(t.proj)); complete.ocheat.rcp26.95 = rep(NA,length(t.proj))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.proj)){
  # priors
  c(priors.slr.rcp26.05[t], priors.slr.rcp26.50[t], priors.slr.rcp26.95[t]) := quantile(priors.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.ais.rcp26.05[t], priors.ais.rcp26.50[t], priors.ais.rcp26.95[t]) := quantile(priors.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gis.rcp26.05[t], priors.gis.rcp26.50[t], priors.gis.rcp26.95[t]) := quantile(priors.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gsic.rcp26.05[t], priors.gsic.rcp26.50[t], priors.gsic.rcp26.95[t]) := quantile(priors.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.te.rcp26.05[t], priors.te.rcp26.50[t], priors.te.rcp26.95[t]) := quantile(priors.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.temp.rcp26.05[t], priors.temp.rcp26.50[t], priors.temp.rcp26.95[t]) := quantile(priors.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.ocheat.rcp26.05[t], priors.ocheat.rcp26.50[t], priors.ocheat.rcp26.95[t]) := quantile(priors.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # experts
  c(experts.slr.rcp26.05[t], experts.slr.rcp26.50[t], experts.slr.rcp26.95[t]) := quantile(experts.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.ais.rcp26.05[t], experts.ais.rcp26.50[t], experts.ais.rcp26.95[t]) := quantile(experts.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gis.rcp26.05[t], experts.gis.rcp26.50[t], experts.gis.rcp26.95[t]) := quantile(experts.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gsic.rcp26.05[t], experts.gsic.rcp26.50[t], experts.gsic.rcp26.95[t]) := quantile(experts.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.te.rcp26.05[t], experts.te.rcp26.50[t], experts.te.rcp26.95[t]) := quantile(experts.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.temp.rcp26.05[t], experts.temp.rcp26.50[t], experts.temp.rcp26.95[t]) := quantile(experts.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.ocheat.rcp26.05[t], experts.ocheat.rcp26.50[t], experts.ocheat.rcp26.95[t]) := quantile(experts.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # standard
  c(standard.slr.rcp26.05[t], standard.slr.rcp26.50[t], standard.slr.rcp26.95[t]) := quantile(standard.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.ais.rcp26.05[t], standard.ais.rcp26.50[t], standard.ais.rcp26.95[t]) := quantile(standard.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gis.rcp26.05[t], standard.gis.rcp26.50[t], standard.gis.rcp26.95[t]) := quantile(standard.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gsic.rcp26.05[t], standard.gsic.rcp26.50[t], standard.gsic.rcp26.95[t]) := quantile(standard.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.te.rcp26.05[t], standard.te.rcp26.50[t], standard.te.rcp26.95[t]) := quantile(standard.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.temp.rcp26.05[t], standard.temp.rcp26.50[t], standard.temp.rcp26.95[t]) := quantile(standard.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.ocheat.rcp26.05[t], standard.ocheat.rcp26.50[t], standard.ocheat.rcp26.95[t]) := quantile(standard.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # complete
  c(complete.slr.rcp26.05[t], complete.slr.rcp26.50[t], complete.slr.rcp26.95[t]) := quantile(complete.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.ais.rcp26.05[t], complete.ais.rcp26.50[t], complete.ais.rcp26.95[t]) := quantile(complete.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gis.rcp26.05[t], complete.gis.rcp26.50[t], complete.gis.rcp26.95[t]) := quantile(complete.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gsic.rcp26.05[t], complete.gsic.rcp26.50[t], complete.gsic.rcp26.95[t]) := quantile(complete.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.te.rcp26.05[t], complete.te.rcp26.50[t], complete.te.rcp26.95[t]) := quantile(complete.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.temp.rcp26.05[t], complete.temp.rcp26.50[t], complete.temp.rcp26.95[t]) := quantile(complete.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.ocheat.rcp26.05[t], complete.ocheat.rcp26.50[t], complete.ocheat.rcp26.95[t]) := quantile(complete.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
}

iproj = which(t.proj==2000):which(t.proj==2100)
i2100 = which(t.proj==2100)

# compute SLR densities
# GIS
priors.gis.2100 <- priors.gis.rcp26[nrow(priors.gis.rcp26),]
experts.gis.2100 <- experts.gis.rcp26[nrow(experts.gis.rcp26),]
standard.gis.2100 <- standard.gis.rcp26[nrow(standard.gis.rcp26),]
complete.gis.2100 <- complete.gis.rcp26[nrow(complete.gis.rcp26),]

priors.density.gis.2100 <- density(priors.gis.2100)
experts.density.gis.2100 <- density(experts.gis.2100)
standard.density.gis.2100 <- density(standard.gis.2100)
complete.density.gis.2100 <- density(complete.gis.2100)

# AIS
priors.ais.2100 <- priors.ais.rcp26[nrow(priors.ais.rcp26),]
experts.ais.2100 <- experts.ais.rcp26[nrow(experts.ais.rcp26),]
standard.ais.2100 <- standard.ais.rcp26[nrow(standard.ais.rcp26),]
complete.ais.2100 <- complete.ais.rcp26[nrow(complete.ais.rcp26),]

priors.density.ais.2100 <- density(priors.ais.2100)
experts.density.ais.2100 <- density(experts.ais.2100)
standard.density.ais.2100 <- density(standard.ais.2100)
complete.density.ais.2100 <- density(complete.ais.2100)

# TOTAL SLR
priors.slr.2100 <- priors.slr.rcp26[nrow(priors.slr.rcp26),]
experts.slr.2100 <- experts.slr.rcp26[nrow(experts.slr.rcp26),]
standard.slr.2100 <- standard.slr.rcp26[nrow(standard.slr.rcp26),]
complete.slr.2100 <- complete.slr.rcp26[nrow(complete.slr.rcp26),]

priors.density.2100 <- density(priors.slr.2100)
experts.density.2100 <- density(experts.slr.2100)
standard.density.2100 <- density(standard.slr.2100)
complete.density.2100 <- density(complete.slr.2100)

# TEMPERATURE
priors.temp.2100 <- priors.temp.rcp26[nrow(priors.temp.rcp26),]
experts.temp.2100 <- experts.temp.rcp26[nrow(experts.temp.rcp26),]
standard.temp.2100 <- standard.temp.rcp26[nrow(standard.temp.rcp26),]
complete.temp.2100 <- complete.temp.rcp26[nrow(complete.temp.rcp26),]

priors.density.temp.2100 <- density(priors.temp.2100)
experts.density.temp.2100 <- density(experts.temp.2100)
standard.density.temp.2100 <- density(standard.temp.2100)
complete.density.temp.2100 <- density(complete.temp.2100)

# TE
priors.te.2100 <- priors.te.rcp26[nrow(priors.te.rcp26),]
experts.te.2100 <- experts.te.rcp26[nrow(experts.te.rcp26),]
standard.te.2100 <- standard.te.rcp26[nrow(standard.te.rcp26),]
complete.te.2100 <- complete.te.rcp26[nrow(complete.te.rcp26),]

# GSIC
priors.gsic.2100 <- priors.gsic.rcp26[nrow(priors.gsic.rcp26),]
experts.gsic.2100 <- experts.gsic.rcp26[nrow(experts.gsic.rcp26),]
standard.gsic.2100 <- standard.gsic.rcp26[nrow(standard.gsic.rcp26),]
complete.gsic.2100 <- complete.gsic.rcp26[nrow(complete.gsic.rcp26),]


##==============================================================================
##==============================================================================



##########################################################################################################
# Figure 4 - GIS parameters

# start post rejection-sampled parameter sets
chain_priors <- readRDS("priors_0426_parameters_good.rds")
chain_experts <- readRDS("experts_0426_parameters_good.rds")
chain_standard <- readRDS("standard_0426_parameters_good.rds")
chain_complete <- readRDS("complete_0426_parameters_good.rds")

# how many to subsample to?
n.ensemble <- 2e6 # subsample to one million

##==============================================================================
## subsample

# priors
if (nrow(chain_priors)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_priors))
  sampleme <- seq(from=1, to=nrow(chain_priors),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_priors[k,]
  }
  
  chain_priors <- parameters
}

# experts
if (nrow(chain_experts)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_experts))
  sampleme <- seq(from=1, to=nrow(chain_experts),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_experts[k,]
  }
  
  chain_experts <- parameters
}

# standard
if (nrow(chain_standard)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_standard))
  sampleme <- seq(from=1, to=nrow(chain_standard),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_standard[k,]
  }
  
  chain_standard <- parameters
}

# complete
if (nrow(chain_complete)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_complete))
  sampleme <- seq(from=1, to=nrow(chain_complete),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_complete[k,]
  }
  
  chain_complete <- parameters
}

if (exists("parameters")){
  rm(parameters) 
}

##==============================================================================
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
##==============================================================================


## get posterior densities
chain_priors_densities <- replicate(length(parnames), vector("list", 7), simplify = FALSE)
chain_experts_densities <- chain_priors_densities
chain_standard_densities <- chain_priors_densities
chain_complete_densities <- chain_priors_densities


for (pp in 1:length(parnames)){
  chain_priors_densities[[pp]] <- density(chain_priors[,pp])
}

for (pp in 1:length(parnames)){
  chain_experts_densities[[pp]] <- density(chain_experts[,pp])
}

for (pp in 1:length(parnames)){
  chain_standard_densities[[pp]] <- density(chain_standard[,pp])
}

for (pp in 1:length(parnames)){
  chain_complete_densities[[pp]] <- density(chain_complete[,pp])
}

## get nice colors
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array

## mycol colors, by row:
# 1 black
# 2 dark blue
# 3 aqua
# 4 pink
# 5 light orange
# 6 purple
# 7 blue
# 8 light purple
# 9 light blue
# 10 very light blue
# 11 dark red
# 12 brown
# 13 dark orange
# 14 neon green
# 15 neon yellow

mycol_experts <- 9 # originally 5
mycol_standard <- 1
mycol_complete <- 9
mycol_priors <- 1

## draw pdfs
sequence_GIS = c(9,10,11,12,13,31,32)
sequence_GIS_select = c(9,10,11,12,13)

cexlab = 1.2
mylwd = 2

if(TRUE){
  layout(cbind(c(1,3,5,7),c(2,4,6,7)))
  par(mai=c(.6,.8,.3,.3) #c(bottom, left, top, right)
  )
  for (pp in sequence_GIS_select){
    
    myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
    
    # experts
    plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
         ylim = c(0,myylim),
         type="l",
         xlab = "", #xlab=parnames[pp],
         ylab='Probability density', lwd=1, lty="dashed", 
         col=mycol.rgb[mycol_experts],
         yaxt='n', cex.lab=cexlab)
    

    
    if(pp==9){ # a.simple
      mtext(side=3, text=expression(bold('a')), line=.25, cex=.9, adj=0);
      mtext("equilibrium volume sensitivity to",side=1,line=2, cex=.8)
      mtext(expression(paste("temperature changes (mSLE","\u00B0C"^-1,")" )),side=1,line=3.4, cex=.8)
    }
    
    if(pp==10){ # b.simple
      mtext(side=3, text=expression(bold('b')), line=.25, cex=.9, adj=0);
      mtext("equilibrium volume for zero",side=1,line=2, cex=.8)
      mtext(expression(paste("temperature anomaly (mSLE)")),side=1,line=3.4, cex=.8)
    }
    
    if(pp==11){ # alpha.simple
      mtext(side=3, text=expression(bold('c')), line=.25, cex=.9, adj=0);
      mtext("sensitivity of e-folding timescale of GIS volume response",side=1,line=2, cex=.8)
      mtext(expression(paste("to changes in temperature (","\u00B0C"^-1,"yr"^-1,")")),side=1,line=3.4, cex=.8)
    }
    
    if(pp==12){ # beta.simple
      mtext(side=3, text=expression(bold('d')), line=.25, cex=.9, adj=0);
      mtext("equilibrium timescale of GIS volume response",side=1,line=2, cex=.8)
      mtext(expression(paste("to changes in temperature (yr"^-1,")")),side=1,line=3.4, cex=.8)
    }
    
    if(pp==13){ # V0
      mtext(side=3, text=expression(bold('e')), line=.25, cex=.9, adj=0);
      mtext("initial condition for GIS contributions",side=1,line=2, cex=.8)
      mtext(expression(paste("to GMSL (mSLE)")),side=1,line=3.4, cex=.8)
    }
    
    # priors
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
          col="black", lwd=1, lty=2)
    
    # standard
    lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,
          col=mycol.rgb[mycol_standard], 
          cex.lab=cexlab)
    
    # complete
    lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,
          col=mycol.rgb[mycol_complete],
          cex.lab=cexlab)
    
    if(pp==9){
      # 
      # text(x=-1.7,y=.39,labels=expression(bold("BRICK")), cex=1.4)
      # text(x=-0.8,y=.25,labels=expression(bold("BRICK+SEJ")), cex=1.4, col=mycol.rgb[mycol_complete])
      
    }
    
  }
  
  plot.new() # space
  
  ## LEGEND
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab=""
       # ,main="Global mean sea level at 2100 [m]"
  )
  
  legend(x = "top",inset = 0,

         legend = c("expert-data-model calibration" , 
                    "data-model calibration",
                    "expert-informed priors",
                    "standard priors"
         ), 
         
         lwd=c(2,2,1,1),
         
         lty=c("solid",
               "solid",
               "dashed",
               "dashed"
         ),
         bty='n', cex=cexlab,
         col=c(mycol.rgb[mycol_complete],
               mycol.rgb[mycol_standard],
               mycol.rgb[mycol_experts],
               mycol.rgb[mycol_priors]
               # "#c77d7dff"
         ), 
         horiz = FALSE,
         
  )

}

##########################################################################################################

t.end <- proc.time()
time.elapsed <- t.end - t.beg
print("line 800")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")
