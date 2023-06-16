#==========================================================================================================================================================
# Objective Function                                                                                                                                      #
# Objective function for use with DEoptim                                                                                                                 #
# Carl Fredrick G. Aquino / cga5030 / 994085920                                                                                                           #
#==========================================================================================================================================================
findRMSE <- function(
  pars,
  S ,    
  kappa ,
  alpha.doeclim, 
  mod.time,
  temp,
  hindcast.forcing,
  forcing,
  projections_start,
  endyear,
  duration,
  series.time) 
{
  #============================================================================= 
  # Run standard projections and calculate RMSE
  # Calculate Total Forcing
  
  source('fortran/R/doeclimF_carl.R')      # the DOECLIM model - resets the mod.time
  
  # initialize projection forcings
  calc.forcing = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 2))
  names(calc.forcing) = c('year','total_forcing')
  calc.forcing[,1] = forcing$year
  calc.forcing = subset(calc.forcing, year>=projections_start&year<=endyear) # trim to desired years
  
  # initialize pars
  # pars = c(-10,10,-11,11,-12,12,-13,13,-14,14,-15,15,-16,16,-17,17,-18,18,-19,19,-20,20) # DELETE THIS LATER
  
  t <- seq(from = 2000, to = 2100, by = 1)
  t_0 <- seq(from = 2000, to = 2100, by = 10)
  rows = length(t)
  nbases = length(t_0)
  npars = nbases*2
  
  basis.series = data.frame(matrix(NA, nrow = rows, ncol = nbases+2))
  names(basis.series) = c("basis1","basis2","basis3","basis4","basis5","basis6","basis7","basis8","basis9","basis10","basis11","total","time")
  basis.series$time=series.time
  
  deoptim.pars = data.frame(matrix(NA, nrow = 11, ncol = 3))
  names(deoptim.pars) = c('year','weight','phi')
  deoptim.pars$year=c(t_0)

  for (i in 1:nbases) {
    deoptim.pars[i,2] = pars[(2*i)-1]
    deoptim.pars[i,3] = pars[2*i]
  }
  
  y=0
  
  for (i in 1:nbases){
    for (j in 1:duration){
      y=deoptim.pars[i,2]*exp(-(abs(deoptim.pars[i,3])*abs(t[j]-t_0[i]))^2)
      basis.series[j,i]=y
    }
  }
  
  basis.series$total = rowSums(basis.series[,1:11])
  
  for (i in 1:duration){
    calc.forcing[i,2]=basis.series$total[i]
  }
  
  # rbind hindcast forcings to projection forcings
  hindcast.total_forcing = data.frame(matrix(NA, nrow = nrow(hindcast.forcing), ncol = 2))
  names(hindcast.total_forcing) = c('year','total_forcing')
  hindcast.total_forcing[,1] = hindcast.forcing$year
  hindcast.total_forcing$total_forcing=hindcast.forcing$non_aerosol+hindcast.forcing$aerosol

  
  # calculate total forcings
  forcing.total_df = rbind(hindcast.total_forcing,calc.forcing)
  forcing.total = forcing.total_df$total_forcing
  # deoptim.forcing[,3] = deoptim.forcing[,3]*alpha.doeclim
  # calc.forcing[,3] = calc.forcing[,3]*alpha.doeclim # for debugging
  
  # forcing.total = deoptim.forcing$non_aerosol + deoptim.forcing$aerosol
  # forcing.total = calc.forcing$non_aerosol + calc.forcing$aerosol # for debugging
  
  # Run DOECLIM
  projections = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=forcing.total, mod.time=mod.time)
  
  # Grab Temps
  doeclim.temp = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 2))
  names(doeclim.temp) = c('year','doeclim.temp')
  doeclim.temp[,1] = forcing[,1]
  doeclim.temp[,2] = projections$temp
  doeclim.temp = subset(doeclim.temp, year>=projections_start&year<=endyear) # trim to desired years
  
  # Calculate RMSE
  for (i in 1:duration) {
    temp[i,2] = (doeclim.temp[i,2]-temp[i,2])^2
  }
  
  RMSE = sqrt(sum(temp[,2])/duration)
  
  return(RMSE)
  
}
#==========================================================================================================================================================









#==========================================================================================================================================================
# temp_series Driver                                                                                                                                             #
# Create DOECLIM outputs consistent with Bamber et al. 2019                                                                                               #
# Carl Fredrick G. Aquino / cga5030 / 994085920                                                                                                           #
#==========================================================================================================================================================

#=============================================================================
# Intro and configuration

library(DEoptim)

## Source some useful functions for manipulating data
# source('../R/forcing_total.R')     # function to add up the total forcing

## Source the model
# source('../fortran/R/doeclimF.R')      # the DOECLIM model - resets the mod.time

#============================================================================= 
# Interpret global mean surface air temperatures SAT trajectories (Bamber et al. 2019 Figure S1)
initial_timestamp = Sys.time()
l.project = FALSE
begyear  = 1850
projections_start = 2000
endyear  = 2100
tstep    = 1
spinup.time = seq(from=begyear, to=projections_start-1, by=tstep)
series.time= seq(from=projections_start, to=endyear, by=tstep)
mod.time= seq(from=begyear, to=endyear, by=tstep)

duration = length(series.time)

forcing <- read.csv( 'data/forcing_rcp26.csv', header=TRUE )
# forcing85 <- read.csv( '/data/forcing_rcp85.csv', header=TRUE )


## Initialize dataframes for L and H temperature series

# L(+2deg C)
SAT_L <- data.frame(matrix(NA, nrow = nrow(forcing), ncol = 2))
names(SAT_L) = c('year','SAT')
SAT_L[,1] = forcing[,1]
SAT_L = subset(SAT_L, year>=projections_start&year<=endyear) # trim to desired years

# H(+5deg C)
SAT_H = SAT_L

# Populate SAT_L and SAT_H temperature series

# L(+2deg C)
# set intercepts
L_intercept_2000 = 0.65 # (+deg C)
L_intercept_2050 = 1.50 # (+deg C)
L_intercept_2100 = 2.00 # (+deg C)

# determine slopes of each of the four time segments:
L_slope1 = (L_intercept_2050 - L_intercept_2000) /50  # (+deg C/year)
L_slope2 = (L_intercept_2100 - L_intercept_2050) /50  # (+deg C/year)
L_slope3 = 0 # (+deg C/year)
L_slope4 = 0 # (+deg C/year)

# fill in dataframe through 2050
SAT_L[1,2] = L_intercept_2000

for (i in 2:51) {
  SAT_L[i,2] = SAT_L[i-1,2] + L_slope1
}

# fill in dataframe through 2100
if (duration > 51){
  for (i in 52:101) {
    SAT_L[i,2] = SAT_L[i-1,2] + L_slope2
  }
}

# fill in dataframe through 2200
if (duration>101){
  for (i in 102:duration) {
    SAT_L[i,2] = SAT_L[i-1,2] + L_slope3
  }
}

# fill in dataframe through 2300
if (duration>201){
  for (i in 202:duration) {
    SAT_L[i,2] = SAT_L[i-1,2] + L_slope3
  }
}

# H(+5deg C)
# set intercepts
H_intercept_2000 = L_intercept_2000 # (+deg C)
H_intercept_2050 = 2.00 # (+deg C)
H_intercept_2100 = 5.00 # (+deg C)

# determine slopes of each of the four time segments:
H_slope1 = (H_intercept_2050 - H_intercept_2000) / 50 # (+deg C/year)
H_slope2 = (H_intercept_2100 - H_intercept_2050) / 50 # (+deg C/year)
H_slope3 = 0 # (+deg C/year)
H_slope4 = 0 # (+deg C/year) 

# fill in dataframe through 2050
SAT_H[1,2] = H_intercept_2000

for (i in 2:51) {
  SAT_H[i,2] = SAT_H[i-1,2] + H_slope1
}

# fill in dataframe through 2100
if (duration>51){
  for (i in 52:101) {
    SAT_H[i,2] = SAT_H[i-1,2] + H_slope2
  }
}

# fill in dataframe through 2200
if (duration>101){
  for (i in 102:duration) {
    SAT_H[i,2] = SAT_H[i-1,2] + H_slope3
  }
}

# fill in dataframe through 2300
if (duration>201){
  for (i in 202:duration) {
    SAT_H[i,2] = SAT_H[i-1,2] + H_slope4
  }
}

#=============================================================================
# Create DOECLIM outputs consistent with Bamber et al. 2019
  
  # Set up DOECLIM #took this section from BRICK_parameterSetup.R
    parnames.in   =NULL; p0.doeclim       =NULL; bound.lower.doeclim=NULL;
    bound.upper.doeclim=NULL; step.mcmc.doeclim=NULL; index.model.doeclim=NULL;
    
    parnames.in   =c("S" ,"kappa.doeclim","alpha.doeclim","T0"  ,"H0" ,"sigma.T","sigma.H","rho.T","rho.H")	# parameters names
    parameters.in         =c(3.1 , 3.5   , 1.1           , -0.06, -33 , 0.1     , 2       , 0.55  , 0.9   )	# initial parameter guesses
    bound.lower.doeclim=c(0.1 , 0.1   , 0             , -0.3 , -50 , 0.05    , 0.1     , 0     , 0     )	# prior range lower bounds
    bound.upper.doeclim=c(10  , 4     , 2             ,  0.3 ,   0 , 5       , 10      , 0.999 , 0.999 )	# prior range upper bounds
    step.mcmc.doeclim  =c(0.16, 0.17  ,0.025          ,0.003 , 0.9 , 5e-4    , 0.025   , 0.007 , 0.006 )	# step size for parameters in MCMC (proposals)
    index.model.doeclim=c(1,2,3,4,5)		# which are model parameters? (index within parnames.doeclim)
    
    ## Grab the DOECLIM parameters
      S            =parameters.in[match("S"            ,parnames.in)]
      kappa.doeclim=parameters.in[match("kappa.doeclim",parnames.in)]
      alpha.doeclim=parameters.in[match("alpha.doeclim",parnames.in)]
      T0           =parameters.in[match("T0"           ,parnames.in)]
      H0           =parameters.in[match("H0"           ,parnames.in)]
    
#=============================================================================
# Use DEoptim to run projections out to 2300
      
  ## Source the model
  source('fortran/R/doeclimF_carl.R')      # the DOECLIM model - resets the mod.time
  
  ## Calculate forcings
    forcing = subset(forcing, year>=begyear&year<=endyear) # start at start year
    hindcast.forcing = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 3))
    names(hindcast.forcing) = c('year','non_aerosol','aerosol')
    hindcast.forcing[,1] = forcing$year
    hindcast.forcing[,2] = forcing$co2 + forcing$nonco2 + forcing$solar + forcing$volcanic + forcing$other
    hindcast.forcing[,3] = forcing$aerosol.direct + forcing$aerosol.indirect
    hindcast.forcing = subset(hindcast.forcing, year>=begyear&year<=projections_start-1) # trim to desired years
    
  ## Set up DEoptim
  npars=22
  lower.bounds = rep(NA,npars)
  upper.bounds = rep(NA,npars)
  NP_scale = 11  
  n_iter = 6000 # Tony ran 100 # Vivek tends to do 2000 - 5000. # Start with 10 and then see what happens
  
  for (i in 1:11) {
    lower.bounds[(i*2)-1]=-10
    lower.bounds[(i*2)]=-11
  }
  
  for (i in 1:11) {
    upper.bounds[(i*2)-1]=10
    upper.bounds[(i*2)]=11
  }

  ## Call DEoptim to infer some forcings
  DEoptim.output = DEoptim(findRMSE,lower.bounds,upper.bounds,control=list(NP=NP_scale*npars, itermax=n_iter, parallelType=1,
                                                                           parVar = list("S","kappa.doeclim","alpha.doeclim","mod.time","SAT_L",
                                                                                         "hindcast.forcing","forcing","projections_start","endyear",
                                                                                         "duration","series.time")),
                                                                  S=S, kappa=kappa.doeclim, alpha.doeclim=alpha.doeclim, mod.time=mod.time,
                                                                  temp=SAT_L, hindcast.forcing=hindcast.forcing, forcing=forcing,
                                                                  projections_start=projections_start, endyear=endyear, duration=duration, series.time=series.time)


  ## Save RDS
  saveRDS(DEoptim.output,file <- "DEoptim_output_RBF_low6000.rds")
  final_timestamp <- Sys.time()
  times <- c(initial_timestamp,final_timestamp)
  saveRDS(times,file <- "times_output_RBF_low6000.rds")
  
  