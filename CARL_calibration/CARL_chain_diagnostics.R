## PRODUCE DIAGNOSTIC PDFs

rm(list=ls())                        # Clear all previous variables
graphics.off()                       # clear plots

## CONFIGURE
load("CARL_calib_MCMC_Complete_0614_2e7_14Jun2021.RData") # configure this
configure <- '_Complete_0614_2e7_' # configure this
## specify windows for marginal distributions
chain1 <- amcmc.out$samples[2e6:1e7,] # 5 million to 20 million   # configure this
chain2 <- amcmc.out$samples[1e7:2e7,] # 20 million to 40 million  # configure this

t.beg <- proc.time()

## get posterior densities
chain1_densities <- replicate(length(parnames), vector("list", 7), simplify = FALSE)
chain2_densities <- chain1_densities

for (pp in 1:length(parnames)){
  chain1_densities[[pp]] <- density(chain1[,pp])
}

for (pp in 1:length(parnames)){
  chain2_densities[[pp]] <- density(chain2[,pp])
}

## draw pdfs
sequence1 = c(26,23,24,22,20,14,15,17,19,21,16,18,25,27,28)
sequence2 = c(5,6,7,8)
sequence3 = c(9,10,11,12,13)
sequence4 = c(1,2,3,4)
sequence5 = c(29,30,31,32,33)

print("begin Antarctic PDFs")
filename.antarctic <- paste('../output_calibration/!!PDFs',configure,'Set01','.png',sep='')
png(filename=filename.antarctic, width=1280, height=720, units ="px")
par(mfrow=c(3,5))
for (pp in sequence1){
  plot(chain1_densities[[pp]]$x,chain1_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain2_densities[[pp]]$x,chain2_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}
dev.off()

print("begin rest of PDFs")
filename.therest <- paste('../output_calibration/!!PDFs',configure,'Set02','.png',sep='')
png(filename=filename.therest, width=1280, height=720, units ="px")
par(mfrow=c(3,6))
for (pp in sequence2){
  plot(chain1_densities[[pp]]$x,chain1_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain2_densities[[pp]]$x,chain2_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}

for (pp in sequence3){
  plot(chain1_densities[[pp]]$x,chain1_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain2_densities[[pp]]$x,chain2_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}

for (pp in sequence4){
  plot(chain1_densities[[pp]]$x,chain1_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain2_densities[[pp]]$x,chain2_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}

for (pp in sequence5){
  plot(chain1_densities[[pp]]$x,chain1_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain2_densities[[pp]]$x,chain2_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}
dev.off()


# save results in case you need to revisit later
t.end <- proc.time()
time.elapsed <- t.end - t.beg
print("line 78")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")

