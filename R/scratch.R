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
#filename.normalgamma$norfolk <-  paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_normalgamma_decl3-pot99_31Mar2018.nc', sep='')
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


# can reformat things as a matrix!
# -> rows are different ensemble members
# -> columns are different data length experiments
# different matrix for each site, and for each parameter

# then return levels:


# calculate return levels

# vectorized?
rl20 <- vector('list', length(site.names)); names(rl20) <- site.names
for (site in site.names) {
  rl20[[site]] <- vector('list', length(parameters_posterior[[site]])); names(rl20[[site]]) <- names(parameters_posterior[[site]])
  for (data.exp in names(parameters_posterior[[site]])) {
    n.ensemble <- nrow(parameters_posterior[[site]][[data.exp]])
    rl20[[site]][[data.exp]] <- rlevd(20, scale=parameters_posterior[[site]][[data.exp]][,2],
                                             shape=parameters_posterior[[site]][[data.exp]][,3],
                                             threshold=data_calib[[site]][[data.exp]]$threshold,
                                             type='GP', npy=365.25,
                                             rate=parameters_posterior[[site]][[data.exp]][,1])
  }
}


# example
ks.test(c(rl20$delfzijl$y137,rl20$delfzijl$y134,rl20$delfzijl$y131), c(rl20$delfzijl$y122) )

#===============================================================================
# End
#===============================================================================
