#==========================================================================================================================================================
# MAE                                                                                                                                                 #
# Calculates the mean absolute error                                                                                                               #
# Carl Fredrick G. Aquino / cga5030 / 994085920                                                                                                           #
#==========================================================================================================================================================
carl_MAE <- function(actual,
                      predicted) {
  N <- length(actual)
  RMSE <- sum(abs(predicted-actual))/N
}

#==========================================================================================================================================================
# RMSE                                                                                                                                                    #
# Calculates the root mean square error                                                                                                                   #
# Carl Fredrick G. Aquino / cga5030 / 994085920                                                                                                           #
#==========================================================================================================================================================
carl_RMSE <- function(actual,
                      predicted) {
  N <- length(actual)
  RMSE <- sqrt(sum((predicted-actual)^2)/N)
}


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
  duration) 
{
  #============================================================================= 
  # Run standard projections and calculate RMSE
  # Calculate Total Forcing
  
  # populate projection forcings
  calc.forcing = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 3))
  names(calc.forcing) = c('year','non_aerosol','aerosol')
  calc.forcing[,1] = forcing$year
  calc.forcing = subset(calc.forcing, year>=projections_start&year<=endyear) # trim to desired years
  
  for (i in 1:duration) {
    calc.forcing[i,2] = pars[(2*i)-1]
    calc.forcing[i,3] = pars[2*i]
  }
  
  # rbind hindcast forcings to projection forcings
  deoptim.forcing = rbind(hindcast.forcing,calc.forcing)
  
  # calculate total forcings
  deoptim.forcing[,3] = deoptim.forcing[,3]*alpha.doeclim
  # calc.forcing[,3] = calc.forcing[,3]*alpha.doeclim # for debugging
  
  forcing.total = deoptim.forcing$non_aerosol + deoptim.forcing$aerosol
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
  for (i in 1:length(duration)) {
    temp[i,2] = (doeclim.temp[i,2]-temp[i,2])^2
  }
  
  RMSE = sqrt(sum(temp[,2])/duration)
  
  return(RMSE)
  
}
#==========================================================================================================================================================



#==========================================================================================================================================================
# temp_series                                                                                                                                             #
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

l.project = TRUE
begyear  = 1850
projections_start = 2000
endyear  = 2100
tstep    = 1
spinup.time = seq(from=begyear, to=projections_start-1, by=tstep)
series.time= seq(from=projections_start, to=endyear, by=tstep)
mod.time= seq(from=begyear, to=endyear, by=tstep)

duration = length(series.time)

forcing <- read.csv( '../data/forcing_rcp26.csv', header=TRUE )
forcing85 <- read.csv( '../data/forcing_rcp85.csv', header=TRUE )
forcing26 <- forcing


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
  source('../fortran/R/doeclimF.R')      # the DOECLIM model - resets the mod.time
  source('../fortran/R/dais_fastdynF.R')  # DAIS model  
  
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
    NP_scale = 12  
    n_iter = 100 # Tony ran 100 # Vivek tends to do 2000 - 5000. # Start with 10 and then see what happens
    
    for (i in 1:11) {
      lower.bounds[(i*2)-1]=-10
      lower.bounds[(i*2)]=-11
    }
    
    for (i in 1:11) {
      upper.bounds[(i*2)-1]=10
      upper.bounds[(i*2)]=11
    }
  
  ## Call DEoptim to infer some forcings
    # # DEoptim.output <- readRDS(file.choose())
    # times <- readRDS(file.choose())
    DEoptim.output <- readRDS('DEoptim_low6000_0106.rds', refhook=NULL)
    times <- readRDS('DEoptim_low6000_0106-TIMES.rds', refhook=NULL)
    
    
  ## Extract / interpret forcings
  bestmem = DEoptim.output$optim$bestmem
  
  
  ## Populate forcings
  hindcast.total_forcing = data.frame(matrix(NA, nrow = nrow(hindcast.forcing), ncol = 2))
  names(hindcast.total_forcing) = c('year','total_forcing')
  hindcast.total_forcing[,1] = hindcast.forcing$year
  hindcast.total_forcing$total_forcing=hindcast.forcing$non_aerosol+hindcast.forcing$aerosol
  
  implied.forcing = data.frame(matrix(NA, nrow = nrow(SAT_H), ncol = 2))
  names(implied.forcing) = c('year','total_forcing')
  implied.forcing[,1] = SAT_H[,1]
  
  new.forcing = data.frame(matrix(NA, nrow = nrow(SAT_H), ncol = 2))
  names(new.forcing) = c('year','total_forcing')
  new.forcing[,1] = SAT_H[,1]
  
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
    deoptim.pars[i,2] = bestmem[(2*i)-1]
    deoptim.pars[i,3] = bestmem[2*i]
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
    implied.forcing[i,2]=basis.series$total[i]
  }
  
  ## Create new total forcings
  new.forcing_df = rbind(hindcast.total_forcing,implied.forcing)
  doeclim.forcing = new.forcing_df$total_forcing
  calibration.forcing_df = subset(new.forcing_df, year>=1850&year<=2009) # trim to desired years
  calibration.forcing = calibration.forcing_df$total_forcing
  
  ## save RDS files
  saveRDS(doeclim.forcing,file <- "doeclim.forcing_RBF_low6000_0106.rds")
  saveRDS(calibration.forcing,file <- "calibration_RBF_low6000_0106.rds")
  
  ## Run DOECLIM projections using inferred forcings
  # doeclim.out = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=doeclim.forcing, mod.time=mod.time)
  
  #=============================================================================
  # old plots
  oldplots = FALSE
  if (oldplots){
  # Plots and Analysis - forcings
  plot(doeclim.out$time,doeclim.out$temp, xlab="Year",ylab="Global Mean Surface Temperature Anomaly [deg C]", pch=16)
  lines(series.time,SAT_H[,2], lwd=2, col="blue")
  legend("topleft", c("DOECLIM Outputs (DEoptim - 6000 iterations, NP = 11x22 = 242)", "Prescribed Temperature Series (Bamber et al. 2019)"),
         col = c("black", "blue"), lwd = 2, bty = "n")
  
  plot (doeclim.out$time, doeclim.forcing, xlab="Year",ylab="Total Forcing[W/m^2]", main="Implied Forcings for DOECLIM (n_iter = 6000, NP = 11x22 = 242)")
  }
  
  # plot (doeclim.out$time, new.forcing$aerosol, xlab="Year",ylab="Aerosol Forcing [W/m^2]", main="Implied Aerosol Forcings for DOECLIM (n_iter = 2000, NP = 11x22 = 242)")
  # 
  # plot (doeclim.out$time, new.forcing$non_aerosol, xlab="Year",ylab="Non-Aerosol Forcing[W/m^2]", main="Implied Non-Aerosol Forcings for DOECLIM (n_iter = 2000, NP = 11x202)")
  # 
  # plot (new.forcing$aerosol, new.forcing$non_aerosol, xlab="Aerosol Forcing[W/m^2]",ylab="Non-Aerosol Forcing[W/m^2]", main="Aerosol vs. Non-Aerosol Forcings (n_iter = 2000, NP = 11x202)")
 
  #=============================================================================
  # New analysis / plots

    # RCP comparison 
      # Generate vectors
        forcing26 = subset(forcing26, year>=begyear&year<=endyear) # start at start year
        forcing26.total = forcing26$co2 + forcing26$nonco2 + alpha.doeclim*forcing26$aerosol.direct + alpha.doeclim*forcing26$aerosol.indirect + forcing26$solar + forcing26$volcanic + forcing26$other
        doeclim.out26 = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=forcing26.total, mod.time=mod.time)
        
        forcing85 = subset(forcing85, year>=begyear&year<=endyear) # start at start year
        forcing85.total = forcing85.total = forcing85$co2 + forcing85$nonco2 + alpha.doeclim*forcing85$aerosol.direct + alpha.doeclim*forcing85$aerosol.indirect + forcing85$solar + forcing85$volcanic + forcing85$other
        doeclim.out85 = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=forcing85.total, mod.time=mod.time)
        
        doeclim.forcing_high6000 <- readRDS('doeclim.forcing_RBF_high6000_0106.rds', refhook=NULL)
        doeclim.out_high6000_2 = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=doeclim.forcing_high6000, mod.time=mod.time)
        
        doeclim.forcing_low6000 <- readRDS('doeclim.forcing_RBF_low6000_0106.rds', refhook=NULL)
        doeclim.out_low6000 = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=doeclim.forcing_low6000, mod.time=mod.time)
      
        # make CSVs
        forcing_low6000_df = data.frame(matrix(NA, nrow = length(doeclim.forcing_low6000), ncol = 2))
        names(forcing_low6000_df) = c("year","total forcing")
        forcing_low6000_df$year = forcing[,1]
        forcing_low6000_df$`total forcing`=doeclim.forcing_low6000
        write.csv(forcing_low6000_df, file="../data/forcing_low6000.csv", row.names=FALSE)
        
        forcing_high6000_df = forcing_low6000_df
        forcing_high6000_df$`total forcing`=doeclim.forcing_high6000
        write.csv(forcing_high6000_df, file="../data/forcing_high6000.csv", row.names=FALSE)
        
      # Analysis - RCP
        # Calculate RMSE of RCP8.5
          rcp85temp = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 2))
          names(rcp85temp) = c('year','SAT')
          rcp85temp[,1] = forcing[,1]
          rcp85temp[,2] = doeclim.out85$temp
          rcp85temp = subset(rcp85temp, year>=projections_start&year<=endyear) # trim to desired years
          RMSE_rcp85 = carl_RMSE(actual = SAT_H$SAT,predicted = rcp85temp$SAT)
          
        # # Calculate MAE of RCP8.5
        #   MAE_rcp85 = carl_MAE(actual = SAT_H$SAT, predicted = rcp85temp$SAT)
          
        # Calculate RMSE of RCP2.6
          rcp26temp = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 2))
          names(rcp26temp) = c('year','SAT')
          rcp26temp[,1] = forcing[,1]
          rcp26temp[,2] = doeclim.out26$temp
          rcp26temp = subset(rcp26temp, year>=projections_start&year<=endyear) # trim to desired years
          RMSE_rcp26 = carl_RMSE(actual = SAT_L$SAT,predicted = rcp26temp$SAT)
          
        # # Calculate MAE of RCP2.6
        #   MAE_rcp26 = carl_MAE(actual = SAT_L$SAT,predicted = rcp26temp$SAT)
      
    # Inferred Results
        # Generate vectors
          high6000 = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 2))
          names(high6000) = c('year','SAT')
          high6000[,1] = forcing[,1]
          high6000[,2] = doeclim.out_high6000_2$temp
          high6000 = subset(high6000, year>=projections_start&year<=endyear) # trim to desired years
          
          low6000 = data.frame(matrix(NA, nrow = nrow(forcing), ncol = 2))
          names(low6000) = c('year','SAT')
          low6000[,1] = forcing[,1]
          low6000[,2] = doeclim.out_low6000$temp
          low6000 = subset(low6000, year>=projections_start&year<=endyear) # trim to desired years
        
        # Calculate RMSE 
          RMSE_high6000 = carl_RMSE(actual = SAT_H$SAT,predicted =high6000$SAT)
          RMSE_low6000 = carl_RMSE(actual = SAT_L$SAT,predicted =low6000$SAT)
        
        # # Calculate MAE
        #   MAE_high6000 = carl_MAE(actual = SAT_H$SAT, predicted = high6000$SAT)
        #   MAE_low6000 = carl_MAE(actual = SAT_L$SAT, predicted = low6000$SAT)
        
        # difference_high = abs(high6000$SAT - SAT_H$SAT)
        # logdifference_high = log10(abs(difference_high))
        
        
        newplot=TRUE
        if (newplot==TRUE) {
          # par(mfrow=c(2,1))
          
          # bitmap("plots/plot_RCP.tiff", height = 4, width = 5, units = 'in', type="tiff24nc", res=300)
          plot(x=doeclim.out26$time, y=doeclim.out26$temp, ylim=c(0,5.1), xlim=c(2000,2100), xlab="Year",
               ylab="Global Mean Surface Temperature Anomaly [deg C]", type="l", lwd=4, col="yellow")
          lines(x=doeclim.out85$time, y=doeclim.out85$temp, lwd=4, col="blue",xlim=c(2000,2100))
          lines(series.time,SAT_H[,2], lwd=3, col="black", lty="dotted")
          lines(series.time,SAT_L[,2], lwd=3, col="red", lty="dotted")
          legend(x="topleft", legend=c("Prescribed High Emission Scenario", 
                                       "Prescribed Low Emission Scenario",
                                       "RCP8.5 - DOECLIM Projection", 
                                       "RCP2.6 - DOECLIM Projection"), 
                 col=c("black","red","blue","yellow"), lwd=c(3,3,4,4), bty="n", lty=c(3,3,1,1))
          # dev.off()

          plot(x=low6000$year, y=low6000$SAT, ylim=c(0,5.1), xlim=c(2000,2100), xlab="Year", ylab="Global Mean Surface Temperature Anomaly [deg C]", type="l", lwd=4, col="yellow")
          lines(x=high6000$year, y=high6000$SAT, lwd=4, col="blue",xlim=c(2000,2100))
          lines(series.time,SAT_H[,2], lwd=3, col="black", lty="dotted")
          lines(series.time,SAT_L[,2], lwd=3, col="red", lty="dotted")
          legend(x="topleft", legend=c("Prescribed High Emission Scenario", 
                                       "Prescribed Low Emission Scenario",
                                       "Inferred DOECLIM Projection - High", 
                                       "Inferred DOECLIM Projection - Low"), 
                 col=c("black","red","blue","yellow"), lwd=c(3,3,4,4), bty="n", lty=c(3,3,1,1))
          # plot(x=high6000$year,y=difference_high,xlim=c(2000,2100),ylim=c(0,0.1))
          # plot(x=high6000$year,y=logdifference_high,xlim=c(2000,2100))
          
          plot(x=  doeclim.out_high6000_2$time, y=  doeclim.out_high6000_2$temp, ylim=c(0,5.1), xlim=c(2000,2100), xlab="Year",
               ylab="Global Mean Surface Temperature Anomaly [deg C]", type="l", lwd=4, col="blue")
          lines(series.time,SAT_H[,2], lwd=3, col="black", lty="dotted")
          legend(x="topleft", legend=c("Prescribed High Emission Scenario",
                                       "Inferred Projection - Spline Interpolation"), 
                 col=c("black","blue"), lwd=c(3,4), bty="n", lty=c(3,1))
          
          # Forcings Plots
          plot(x=doeclim.out_high6000_2$time, y=doeclim.forcing_high6000,xlim=c(2000,2100),
               type="l",lwd=3,col="blue",
               xlab="Year",ylab="Total Atmospheric Forcing [W/m^2]")
          lines(x=doeclim.out_high6000_2$time,y=doeclim.forcing_low6000,xlim=c(2000,2100),
                lwd=3,col="yellow")
          legend(x="topleft", legend=c("Inferred Forcings - High","Inferred Forcings - Low"),
                 col=c("blue","yellow"), lwd=c(3,3), bty="n")
          
          plot(x=series.time,SAT_H[,2], type="l", lwd=3, col="black",lty="dotted",xlab="Year",
               ylab="Global Mean Surface Temperature Anomaly [deg C]",ylim=c(0,5))
          lines(series.time,SAT_L[,2], lwd=3, col="red", lty="dotted")
          legend(x="topleft", legend=c("Prescribed High Emission Scenario", 
                                       "Prescribed Low Emission Scenario"),
                                        col=c("black","red"), lwd=c(3,3), bty="n")
          
          # Ocean Heat Uptake Plots
          plot(x=doeclim.out_high6000_2$time, y=doeclim.out_high6000_2$ocheat,xlim=c(2000,2100),
               type="l",lwd=3,col="blue",
               xlab="Year", ylab="Ocean Heat [10^22 J]")
          lines(x=doeclim.out_low6000$time, y=doeclim.out_low6000$ocheat,xlim=c(2000,2100),
                lwd=3, col="yellow")
          legend(x="topleft", legend=c("Inferred Ocean Heat - High","Inferred Ocean Heat - Low"),
                 col=c("blue","yellow"), lwd=c(3,3), bty="n")
          
        }
      
    # Results
 
  #==========================================================================================================================================================
  



