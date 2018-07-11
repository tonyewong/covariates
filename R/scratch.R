#===============================================================================
# scratch.R
#
# Questions? Tony Wong (<anthony.e.wong@colorado.edu>)
#===============================================================================
# Copyright 2018 Tony Wong
#
# MESS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MESS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MESS.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

rm(list=ls())

setwd('~/codes/datalengths/R')

library(ncdf4)
library(extRemes)
library(Bolstad)

# get some colorblind-friendly colors to plot with
source('colorblindPalette.R')

# set useful directories -- assumes you are in the 'R' directory within the repo
# directory structure
plot.dir <- '../figures/'
output.dir <- '../output/'

# site names
site.names <- c('norfolk','delfzijl')

# calibrated parameter sets (samples; all should be same size)
# these are the results from 'calibration_dayPOT-experiments_driver.R'
filename.normalgamma <- vector('list', length(site.names)); names(filename.normalgamma) <- site.names
filename.normalgamma$norfolk <-  paste(output.dir,'calibratedParameters_datalengths_norfolk_normalgamma_decl3-pot99_08Jun2018.nc', sep='')
filename.normalgamma$delfzijl <- paste(output.dir,'calibratedParameters_datalengths_delfzijl_normalgamma_decl3-pot99_08Jun2018.nc',sep='')

filename.datacalib <- "../data/tidegauge_processed_norfolk-delfzijl_decl3-pot99-annual_07Jun2018.rds"
filename.priors <- '../data/surge_priors_normalgamma_ppgpd_decl3-pot99_28Mar2018.rds'

# file to save progress as you run
filename.saveprogress <- '../output/analysis_inprogress.RData'

# can get all experiment names from calibration data
data_calib <- readRDS(filename.datacalib)

# read results for posterior parameters, across both sites and all data lengths

parameters_posterior <- vector('list', length(data_calib)); names(parameters_posterior) <- names(data_calib)
for (site in names(parameters_posterior)) {
  parameters_posterior[[site]] <- vector('list', length(data_calib[[site]]))
  names(parameters_posterior[[site]]) <- names(data_calib[[site]])
}

for (site in site.names) {
  ncdata <- nc_open(filename.normalgamma[[site]])
  for (data.exp in names(parameters_posterior[[site]])) {
    parameters_posterior[[site]][[data.exp]] <- t(ncvar_get(ncdata, paste('parameters.',data.exp,sep='')))
  }
}


# can reformat return levels as a matrix!
# -> rows are different ensemble members
# -> columns are different data length experiments
# different matrix for each site, and for each parameter

# then return levels:


# calculate return levels

rl.years <- c(20, 100)
rl.names <- NULL
for (yy in rl.years) {rl.names <- c(rl.names, paste('y',yy,sep=''))}
rl <- vector('list', length(site.names)); names(rl) <- site.names
for (site in site.names) {
  rl[[site]] <- vector('list', length(rl.names))
  names(rl[[site]]) <- rl.names
}

for (site in site.names) {
  for (yy in 1:length(rl.years)) {
    # initialize return levels matrix with # ensemble members rows and # data length experiment columns
    rl[[site]][[yy]] <- mat.or.vec(nrow(parameters_posterior[[site]][[1]]), length(parameters_posterior[[site]]))
    colnames(rl[[site]][[yy]]) <- names(parameters_posterior[[site]])
    for (data.exp in names(parameters_posterior[[site]])) {
      rl[[site]][[yy]][, data.exp] <- rlevd(rl.years[yy], scale=parameters_posterior[[site]][[data.exp]][,2],
                                            shape=parameters_posterior[[site]][[data.exp]][,3],
                                            threshold=data_calib[[site]][[data.exp]]$threshold,
                                            type='GP', npy=365.25,
                                            rate=parameters_posterior[[site]][[data.exp]][,1])
    }
  }
}


# examples
ks.test(rl$norfolk$y20[,'y89'], rl$norfolk$y20[,'y86'])



library(fitdistrplus)
tmp <- fitdist(rl20$delfzijl$y137, "lnorm", start = c(meanlog=0, sdlog=1))

#===============================================================================
#===============================================================================

data_gpd <- readRDS('../data/tidegauge_processed_norfolk-delfzijl_decl3-pot99-annual_10Jul2018.rds')

years_norfolk <- rep(NA, length(data_gpd$norfolk))
threshold_norfolk <- rep(NA, length(data_gpd$norfolk))
for (yy in 1:length(data_gpd$norfolk)) {
  years_norfolk[yy] <- length(data_gpd$norfolk[[yy]]$year)
  threshold_norfolk[yy] <- data_gpd$norfolk[[yy]]$threshold
}
threshold0_norfolk <- threshold_norfolk[length(threshold_norfolk)]
relerr_norfolk <- (threshold_norfolk - threshold0_norfolk)/threshold0_norfolk

years_delfzijl <- rep(NA, length(data_gpd$delfzijl))
threshold_delfzijl <- rep(NA, length(data_gpd$delfzijl))
for (yy in 1:length(data_gpd$delfzijl)) {
  years_delfzijl[yy] <- length(data_gpd$delfzijl[[yy]]$year)
  threshold_delfzijl[yy] <- data_gpd$delfzijl[[yy]]$threshold
}
threshold0_delfzijl <- threshold_delfzijl[length(threshold_delfzijl)]
relerr_delfzijl <- (threshold_delfzijl - threshold0_delfzijl)/threshold0_delfzijl

# relative error to plot
re <- 0.01

par(mfrow=c(2,1))

plot(years_norfolk, relerr_norfolk, pch=4, ylim=c(-.05,.05),
     xlab='Years of data', ylab='RE', main='Norfolk')
lines(c(0,140),c(0,0), lty=1, col='red')
lines(c(0,140),c(re, re), lty=2, col='red')
lines(c(0,140),-c(re, re), lty=2, col='red')

plot(years_delfzijl, relerr_delfzijl, pch=4, ylim=c(-.05,.05),
     xlab='Years of data', ylab='RE', main='Delfzijl')
lines(c(0,140),c(0,0), lty=1, col='red')
lines(c(0,140),c(re, re), lty=2, col='red')
lines(c(0,140),-c(re, re), lty=2, col='red')



# plot in actual sea levels

#pdf(paste(plot.dir,'returnlevels_pdf_sf.pdf',sep=''),width=7,height=9.5,colormodel='rgb')
pdf(paste('../figures/p99_datalengths.pdf',sep=''),width=4,height=5,colormodel='rgb')

par(mfrow=c(2,1))
par(mai=c(.67,.8,.25,.1))

plot(years_norfolk, threshold_norfolk, pch=4, ylim=c(920,1020), xlab='', ylab='', main='', yaxt='n')
lines(c(0,140),c(threshold0_norfolk, threshold0_norfolk), lty=1, col='red')
mtext('Years of data', side=1, line=2, cex=1)
mtext('99th percentile [mm]', side=2, line=2.8, cex=1)
mtext('Norfolk', side=3, line=.1, cex=1)
#mtext(expression(bold(' c.')), side=3, line=.1, cex=0.75, adj=0)
axis(1, at=seq(10,140,10), labels=c('10','','30','','50','','70','','90','','110','','130',''), cex.axis=1)
axis(2, at=seq(940,1020,20), las=1, cex.axis=1, mgp=c(3,.7,0))

par(mai=c(.67,.8,.25,.1))

plot(years_delfzijl, threshold_delfzijl, pch=4, ylim=c(2380,2580), xlab='', ylab='', main='', yaxt='n')
lines(c(0,140),c(threshold0_delfzijl, threshold0_delfzijl), lty=1, col='red')
mtext('Years of data', side=1, line=2, cex=1)
mtext('99th percentile [mm]', side=2, line=2.8, cex=1)
mtext('Delfzijl', side=3, line=.1, cex=1)
#mtext(expression(bold(' c.')), side=3, line=.1, cex=0.75, adj=0)
axis(1, at=seq(10,140,10), labels=c('10','','30','','50','','70','','90','','110','','130',''), cex.axis=1)
axis(2, at=seq(2380,2580,40), las=1, cex.axis=1, mgp=c(3,.7,0))

dev.off()

#===============================================================================
# End
#===============================================================================
