## PRODUCE TRACE PLOTS FOR CALIBRATION 

rm(list=ls())                        # Clear all previous variables
graphics.off()                       # clear plots

# load("CARL_calib_MCMC_Complete_0614_2e7_14Jun2021.RData")

t.beg <- proc.time()

## Set up a filename for saving RData images along the way
configure <- 'COMPARE_v5_0313_'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.saveprogress <- paste(configure,today,'.RData',sep='')

# # Start with amcmc.out
# amcmc.out_priors   <- readRDS("CARL_calib_amcmc_out_Priors_0115_1e7_15Jan2022.rds")
# amcmc.out_experts  <- readRDS(file = "CARL_calib_amcmc_out_Expert_0115_4e7_15Jan2022.rds",refhook = NULL)
# amcmc.out_standard <- readRDS(file = "CARL_calib_amcmc_out_Standard_0115_extend_to_45e6_17Feb2022.rds",refhook=NULL)
# amcmc.out_complete <- readRDS(file = "CARL_calib_amcmc_out_Complete_0115_extend_to_45e6_17Feb2022.rds",refhook = NULL)
# 
# ## burn-in
# chain_priors <- amcmc.out_priors$samples[4e6:10e6,] # keep 4 million to 10 million
# chain_experts <- amcmc.out_experts$samples[5e6:4e7,] # keep 5 million to 40 million
# chain_standard <- amcmc.out_standard$samples[25e6:45e6,] # keep 25 million to 45 million
# chain_complete <- amcmc.out_complete$samples[25e6:45e6,] # keep 25 million to 45 million

# start post rejection-sampled parameter sets
chain_priors <- readRDS("priors_0303_parameters_good.rds")
chain_experts <- readRDS("experts_0115_parameters_good.rds")
chain_standard <- readRDS("standard_0115_parameters_good.rds")
chain_complete <- readRDS("complete_0115_parameters_good.rds")

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

mycol_experts <- 5
mycol_standard <- 1
mycol_complete <- 9
mycol_priors <- 1

## draw pdfs
sequence1 = c(14,15,16,17,18,19,20,21,22,23,24,25,26,27,28) # DAIS
sequence2 = c(5,6,7,8) # TE
sequence3 = c(9,10,11,12,13) # SIMPLE
sequence4 = c(1,2,3,4) # GSIC
sequence5 = c(29,30,31,32,33) # Statistical

sequence_antarctic = c(14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,33)
sequence_simple = c(9,10,11,12,13,31,32)

cexlab = 1.2
mylwd = 2

## ANTARCTICA
print("begin Antarctic PDFs")
filename.compare1<- paste('../output_calibration/',configure,'Set01','.png',sep='')
png(filename=filename.compare1, width=1920, height=1080, units ="px")
par(mfrow=c(3,5),
    mai=c(.6,.6,.6,.6) #c(bottom, left, top, right)
    )
for (pp in sequence1){
  
  myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
  
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       ylim = c(0,myylim),
       type="l",xlab=parnames[pp], ylab='Probability density', lwd=mylwd, lty="dashed", 
       col=mycol.rgb[mycol_experts],
       yaxt='n', cex.lab=cexlab)
  
  # priors
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=2, lty=2)
  
  # standard
  lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,
        col=mycol.rgb[mycol_standard], 
        cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,
        col=mycol.rgb[mycol_complete],
        cex.lab=cexlab)
  
}
dev.off()

## PRINT THE REST

print("begin rest of PDFs")
filename.compare2<- paste('../output_calibration/',configure,'Set03','.png',sep='')
png(filename=filename.compare2, width=1920, height=1080, units ="px")
par(mfrow=c(4,5),
    mai=c(.6,.6,.6,.6) #c(bottom, left, top, right)
    )
for (pp in sequence2){
  
  myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
  
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       ylim = c(0,myylim),
       type="l",xlab=parnames[pp], ylab='Probability density', lwd=mylwd, lty="dashed",
       col=mycol.rgb[mycol_experts], 
       yaxt='n', cex.lab=cexlab)
  
  # priors
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=2, lty=2)
  
  # standard
  lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_standard], cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_complete], cex.lab=cexlab)
  
}

plot.new() # space

for (pp in sequence3){
  myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
  
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       ylim = c(0,myylim),
       type="l",xlab=parnames[pp], ylab='Probability density', lwd=mylwd, lty="dashed", col=mycol.rgb[mycol_experts], yaxt='n', cex.lab=cexlab)
  
  # priors
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=2, lty=2)
  
  # standard
  lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_standard], cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_complete], cex.lab=cexlab)
  
  }



for (pp in sequence4){
  myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
  
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       ylim = c(0,myylim),
       type="l",xlab=parnames[pp], ylab='Probability density', lwd=mylwd, lty="dashed", col=mycol.rgb[mycol_experts], yaxt='n', cex.lab=cexlab)
  
  # priors
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=2, lty=2)
  
  # standard
  lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_standard], cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_complete], cex.lab=cexlab)
  
}

plot.new() # space

for (pp in sequence5){
  myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
  
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       ylim = c(0,myylim),
       type="l",xlab=parnames[pp], ylab='Probability density', lwd=mylwd, lty="dashed", col=mycol.rgb[mycol_experts], yaxt='n', cex.lab=cexlab)
  
  # priors
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=2, lty=2)
  
  # standard
  lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_standard], cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,col=mycol.rgb[mycol_complete], cex.lab=cexlab)
  
}
dev.off()


# save results in case you need to revisit later
save.image(file=filename.saveprogress)
t.end <- proc.time()
time.elapsed <- t.end - t.beg
print("line 128")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")

