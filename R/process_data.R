#===============================================================================
# process_data.R
#
#
# Questions?  Tony Wong (<anthony.e.wong@colorado.edu>)
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

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/datalengths/R')
}

#===============================================================================

# install some preliminary packages, or load the libraries if you already have them
l.installpackages <- FALSE
if(l.installpackages) {
  install.packages('date')
  install.packages('zoo')
  install.packages('Hmisc')
  install.packages('ncdf4')
  install.packages('extRemes')
  install.packages('adaptMCMC')
  install.packages('lhs')
  install.packages('DEoptim')
  install.packages('foreach')
  install.packages('doParallel')
}
library(date)
library(zoo)
library(Hmisc)
library(ncdf4)
library(extRemes)
library(adaptMCMC)
library(DEoptim)
library(foreach)
library(doParallel)

#===============================================================================

dt.decluster <- 3
detrend.method <- 'linear'
pot.threshold <- 0.99

source('decluster_timeseries.R') # routine for declustering time series

source('processing_norfolk.R') # yields data_norfolk

source('processing_delfzijl.R') # yields data_delfzijl

data_gpd <- vector('list', 2); names(data_gpd) <- c('norfolk', 'delfzijl')
data_gpd$norfolk <- data_norfolk
data_gpd$delfzijl <- data_delfzijl

# write calibration data RDS file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_norfolk-delfzijl_decl',dt.decluster,'-pot',pot.threshold*100,'-',detrend.method,'_',today,'.rds', sep='')
saveRDS(data_gpd, file=filename.output)

#===============================================================================
# End
#===============================================================================
