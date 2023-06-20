# plot_hindcasts.R: this file generates hindcasts from the calibration experiments  

# load packages
library(ggplot2)
library(ggnewscale)

cur_path = "R/plots"
main_path = setwd(cur_path)
rds_path = "../../output_calibration" # location of data files for plotting
plot_path = "../../plots" # output path for plots

# load scripts
source("../util/MultipleOutput.R") # define multiple assignment operator
source("../util/plot_utilities.R") # miscellaneous plot utilities

# load files
t_hind <- readRDS(file=file.path(rds_path, "priors_0426_t_hind.rds")) # hindcast years
t_paleo <- readRDS(file=file.path(rds_path, "priors_0426_t_paleo.rds")) # paleodata years


## years for the hindcast

begyear <- t_hind[1]
endyear <- t_hind[length(t_hind)]
mod_time <- begyear:endyear
begyear_norm <- 1961
endyear_norm <- 1990
ind_norm <- which(mod_time==begyear_norm):which(mod_time==endyear_norm) # normalizing BRICK output
n_time <- length(mod_time)

# Gather up all the data/model indices for comparisons

source('../util/compute_indices.R')  # function to determine the model and data indices for comparisons

## Source the data for hindcast comparisons
source("../calibration/DOECLIM_readData.R")
source("../calibration/GSIC_readData.R")
source("../calibration/SIMPLE_readData.R")
source("../calibration/DAIS_readData.R")
source("../calibration/TE_readData.R")

mod_idx <- list(temp = midx_temp,
                ocheat = midx_ocheat,
                gis = midx_gis,
                gsic = midx_gsic,
                sl = midx_sl)
obs_idx <- list(temp = oidx_temp,
                ocheat = oidx_ocheat,
                gis = oidx_gis,
                gsic = oidx_gsic,
                sl = oidx_sl)
obs_time <- list(temp = obs_temp_time,
                ocheat = obs_ocheat_time,
                gis = obs_gis_time,
                gsic = obs_gsic_time,
                sl = obs_sl_time)

## Gather up all the observations for comparisons
obs <- list(temp = obs_temp,
            ocheat = obs_ocheat,
            gis = obs_gis,
            gsic = obs_gsic,
            sl = obs_sl)
obs_err <- list(temp=obs_temp_err,
                ocheat=obs_ocheat_err,
                gis=obs_gis_err,
                gsic=obs_gsic_err,
                sl=obs_sl_err)

# normalize data to appropriate periods
ind_norm = data.frame(
  temp = c(which(mod_time == 1850), which(mod_time == 1870)),
  ocheat = c(which(mod_time == 1960), which(mod_time == 1990)),
  gsic = c(which(mod_time == 1961), which(mod_time == 1961)),
  gis = c(which(mod_time == 1961), which(mod_time == 1990)),
  sl = c(which(mod_time == 1961), which(mod_time == 1990))
)

for (j in seq_along(obs)) {
  ibeg <- which(obs_time[[j]] == mod_time[ind_norm[[names(obs)[j]]][1]])
  iend <- which(obs_time[[j]] == mod_time[ind_norm[[names(obs)[j]]][2]])
  obs[[j]] <- obs[[j]] - mean(obs[[j]][ibeg:iend])
} 

# load model hindcasts
## no-expert "standard" hindcasts

standard_ais_paleo_05 <- readRDS(file=file.path(rds_path, "standard_0426_ais_paleo_05.rds"))
standard_ais_paleo_50 <- readRDS(file=file.path(rds_path, "standard_0426_ais_paleo_50.rds"))
standard_ais_paleo_95 <- readRDS(file=file.path(rds_path, "standard_0426_ais_paleo_95.rds"))

standard_gsic_hind <- readRDS(file=file.path(rds_path, "standard_0426_gsic_hind.rds"))
standard_te_hind <- readRDS(file=file.path(rds_path, "standard_0426_te_hind.rds"))
standard_gis_hind <- readRDS(file=file.path(rds_path, "standard_0426_gis_hind.rds"))
standard_ais_hind <- readRDS(file=file.path(rds_path, "standard_0426_ais_hind.rds"))
standard_temp_hind <- readRDS(file=file.path(rds_path, "standard_0426_temp_hind.rds"))
standard_ocheat_hind <- readRDS(file=file.path(rds_path, "standard_0426_ocheat_hind.rds"))
standard_gsl_hind <- readRDS(file=file.path(rds_path, "standard_0426_gsl_hind.rds"))

complete_ais_paleo_05 <- readRDS(file=file.path(rds_path, "complete_0426_ais_paleo_05.rds"))
complete_ais_paleo_50 <- readRDS(file=file.path(rds_path, "complete_0426_ais_paleo_50.rds"))
complete_ais_paleo_95 <- readRDS(file=file.path(rds_path, "complete_0426_ais_paleo_95.rds"))

complete_gsic_hind <- readRDS(file=file.path(rds_path, "complete_0426_gsic_hind.rds"))
complete_te_hind <- readRDS(file=file.path(rds_path, "complete_0426_te_hind.rds"))
complete_gis_hind <- readRDS(file=file.path(rds_path, "complete_0426_gis_hind.rds"))
complete_ais_hind <- readRDS(file=file.path(rds_path, "complete_0426_ais_hind.rds"))
complete_temp_hind <- readRDS(file=file.path(rds_path, "complete_0426_temp_hind.rds"))
complete_ocheat_hind <- readRDS(file=file.path(rds_path, "complete_0426_ocheat_hind.rds"))
complete_gsl_hind <- readRDS(file=file.path(rds_path, "complete_0426_gsl_hind.rds"))

# get quantiles for BRICK output
standard_slr_05 = rep(NA,length(t_hind))
standard_slr_50 = rep(NA,length(t_hind))
standard_slr_95 = rep(NA,length(t_hind))
complete_slr_05 = rep(NA,length(t_hind))
complete_slr_50 = rep(NA,length(t_hind))
complete_slr_95 = rep(NA,length(t_hind))
standard_gis_05 = rep(NA,length(t_hind))
standard_gis_50 = rep(NA,length(t_hind))
standard_gis_95 = rep(NA,length(t_hind))
complete_gis_05 = rep(NA,length(t_hind))
complete_gis_50 = rep(NA,length(t_hind))
complete_gis_95 = rep(NA,length(t_hind))

for (t in seq_along(t_hind)) {
  c(standard_gis_05[t],    standard_gis_50[t], standard_gis_95[t]) := quantile(standard_gis_hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete_gis_05[t],    complete_gis_50[t], complete_gis_95[t]) := quantile(complete_gis_hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard_slr_05[t],    standard_slr_50[t], standard_slr_95[t])				:= quantile(standard_gsl_hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete_slr_05[t],    complete_slr_50[t], complete_slr_95[t])				:= quantile(complete_gsl_hind[t,], c(0.05,.50,.95), na.rm=TRUE)
}


# set up antartica precalibration windows

## Note that model output is in meters SLE and these trends are mm, so a
## conversion is necessary.

trends_ais <- c(0.27 , 0.08 , 0.40 )/1000   # m/year (after the /1000)
trends_err <- c(0.11 , 0.185, 0.205)/1000   # m/year (after the /1000)
trends_2up <- trends_ais+2*trends_err
trends_2dn <- trends_ais-2*trends_err
ind_trends <- mat.or.vec(length(trends_ais), 2)
ind_trends[1,] <- c(which(date==-7) , which(date==10)) # 1993-2010
ind_trends[2,] <- c(which(date==-8) , which(date== 1)) # 1992-2001
ind_trends[3,] <- c(which(date== 2) , which(date==11)) # 2002-2011

## Precal window 4:
## Adding observational constraint
estimate_SLE_rate <- abs(-71/360)/1000
time_years <- 2002-1992      # using the midpoint of the 19-year interval
mid_cum_SLE_2002 <- estimate_SLE_rate*time_years
i1992 <- which(date==-8)

estimate_SLE_rate_error <- abs(-53/360)/1000     #1-sigma error
estimate_SLE_error <- sqrt(time_years)*estimate_SLE_rate_error #1-sigma error
# (*sqrt(10) because 10 years of potentially accumulated error:
#  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
#                = 10*year X error^2)
SE2_2002 <- estimate_SLE_error*2 #2-sigma error

positive_2SE <- mid_cum_SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
negative_2SE <- mid_cum_SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

## Precal windows 1-3:
## from Shaffer (2014). modified by Kelsey
upper_wind <- c(7.4, -6.9, -1.25, positive_2SE) # Windows 2-3 from Kelsey, Window 1 from DeConto and Pollard 2016
lower_wind <- c(3.6, -15.8, -4.0, negative_2SE)
#upper.wind <- c(6.0, -6.9, -1.25, positive_2SE) # Windows 1-3 from Kelsey
#lower.wind <- c(1.8, -15.8, -4.0, negative_2SE)
#upper.wind <- c(5.5, -8 , -2, positive_2SE) # Windows 1-3 fFrom Shaffer 2014, p 1809
#lower.wind <- c(2.5, -17, -4, negative_2SE)

windows <- matrix(c(lower_wind, upper_wind), nrow = length(upper_wind), ncol=2)
obs_targets <- (windows[,2]+windows[,1])*.5  # middle of window = obs to compare model to
obs_err_dais <- (windows[,2]-windows[,1])*.5      # half-width of window = uncertainty
obs_err_dais <- 0.5*obs_err_dais                       # assume all windows are 2*stdErr (last two actually are)

## Create a vector with each observation year
## 120kyr, 20Kyr, 6kyr (before present), 2002, and 1993 (first year of the IPCC trend)
obs_years <- c(120000, 220000, 234000, 240002)
ipaleo=which(t_paleo==-149999):which(t_paleo==1)

standard_paleo_ais <- data.frame(
  year = t_paleo[ipaleo], 
  med = standard_ais_paleo_50[ipaleo], 
  lowerci = standard_ais_paleo_05[ipaleo], 
  upperci = standard_ais_paleo_95[ipaleo],
  case = "standard"
  )
complete_paleo_ais <- data.frame(
  year = t_paleo[ipaleo], 
  med = complete_ais_paleo_50[ipaleo], 
  lowerci = complete_ais_paleo_05[ipaleo], 
  upperci = complete_ais_paleo_95[ipaleo],
  case = "complete"
  )
paleo_ais <- rbind(standard_paleo_ais, complete_paleo_ais)
paleo_ais$case <- factor(paleo_ais$case, levels=c("standard", "complete"))
precal_windows <- data.frame(year_low = date[obs_years - 1000], year_high = pmin(date[obs_years + 1000], max(t_paleo[ipaleo]), na.rm = TRUE), lower = windows[, 1], upper = windows[, 2])

## Greenland Ice Sheet
standard_gis <- data.frame(
  year = mod_time, 
  med = standard_gis_50, 
  lowerci = standard_gis_05, 
  upperci = standard_gis_95,
  case = "standard"
  )
complete_gis <- data.frame(
  year = mod_time, 
  med = complete_gis_50, 
  lowerci = complete_gis_05, 
  upperci = complete_gis_95,
  case = "complete"
  )
gis <- rbind(standard_gis, complete_gis)
gis_obs <- data.frame(year = obs_gis_time, data = obs$gis)

## Global Mean Sea Level
standard_gmsl <- data.frame(
  year = mod_time, 
  med = standard_slr_50, 
  lowerci = standard_slr_05, 
  upperci = standard_slr_95,
  case = "standard"
  )
complete_gmsl <- data.frame(
  year = mod_time, 
  med = complete_slr_50, 
  lowerci = complete_slr_05, 
  upperci = complete_slr_95,
  case = "complete"
  )
gmsl <- rbind(standard_gmsl, complete_gmsl)
gmsl_obs <- data.frame(year = obs_sl_time, data = obs$sl)

# plot 
data_colors <- c(complete=mycol.rgb[13], standard=mycol.rgb[7])

## AIS
pais <- ggplot() + theme_classic(base_size = 16) + scale_color_manual("Calibration", values = data_colors, labels = c(complete="Data + Experts", standard = "Data Only")) + scale_fill_manual("Calibration", values = data_colors, labels = c(complete="Data + Experts", standard = "Data Only"))
pais <- pais +
        geom_line(data=paleo_ais, aes(x = year, y = med, color = case), lwd = 1) +
        geom_ribbon(data=paleo_ais, aes(x=year, ymin = lowerci, ymax = upperci, fill = case), alpha = 0.4) +
        new_scale_color() +
        new_scale_fill() +
        geom_rect(data=precal_windows, aes(xmin = year_low, xmax = year_high, ymin = lower, ymax = upper, fill = "Observations", color = "Observations")) +
        scale_fill_manual("", values = mycol.rgb[1]) +
        scale_color_manual("", values = mycol.rgb[1]) +
        geom_hline(yintercept = 0, lty = 2) +
        xlim(min(t_paleo[ipaleo]), max(t_paleo[ipaleo])) +
        xlab("Year") +
        ylab("[m SLE]") +
        ggtitle("(a) Antarctic Ice Sheet") +
        theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

leg <- extract_legend(pais) # extract legend from figure for common use
pais <- pais + theme(legend.position = "none") # remove extracted legend

## GIS
pgis <- ggplot() + theme_classic(base_size = 16) + scale_color_manual("Calibration", values = data_colors, labels = c(complete="Data + Experts", standard = "Data Only")) + scale_fill_manual("Calibration", values = data_colors, labels = c(complete="Data + Experts", standard = "Data Only"))
pgis <- pgis +
        geom_line(data=gis, aes(x = year, y = med, color = case), lwd = 1) +
        geom_ribbon(data=gis, aes(x=year, ymin = lowerci, ymax = upperci, fill = case), alpha = 0.4) +
        new_scale_color() +
        geom_point(data=gis_obs, aes(x=year,  y = data, color = "Observations"), size=1, alpha=0.8) +
        scale_color_manual("", values = mycol.rgb[1]) +
        geom_hline(yintercept = 0, lty = 2) +
        xlab("Year") +
        xlim(c(1950, 2009)) +
        ylab("[m SLE]") +
        ylim(c(-0.003, 0.012)) +
        ggtitle("(b) Greenland Ice Sheet") +
        theme(plot.title = element_text(face = "bold"), legend.position = "none")

## GMSL
pgmsl <- ggplot() + theme_classic(base_size = 16) + scale_color_manual("Calibration", values = data_colors, labels = c(complete="Data + Experts", standard = "Data Only")) + scale_fill_manual("Calibration", values = data_colors, labels = c(complete="Data + Experts", standard = "Data Only"))
pgmsl <- pgmsl +
        geom_line(data=gmsl, aes(x = year, y = med, color = case), lwd = 1) +
        geom_ribbon(data=gmsl, aes(x=year, ymin = lowerci, ymax = upperci, fill = case), alpha = 0.4) +
        new_scale_color() +
        geom_point(data=gmsl_obs, aes(x=year,  y = data, color = "Observations"), size=1, alpha=0.8) +
        scale_color_manual("", values = mycol.rgb[1]) +
        geom_hline(yintercept = 0, lty = 2) +
        xlab("Year") +
        xlim(c(1850, 2009)) +
        ylab("[m SLE]") +
        ylim(c(-0.25, 0.12)) +
        ggtitle("(c) Global Mean Sea Level") +
        theme(plot.title = element_text(face = "bold"), legend.position = "none")

## combine plots into one panel

png("../../plots/hindcast.png", width=180, height=180, units="mm", res=600)
p <- grid.arrange(arrangeGrob(pais, pgis, pgmsl, ncol=1), leg, ncol=1, heights=c(10, 1))
dev.off()
