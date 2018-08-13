#===============================================================================
# process_data.R
#
# Questions?  Tony Wong (<anthony.e.wong@colorado.edu>)
#===============================================================================

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/covariates/R')
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/covariates/R')
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
detrend.method <- 'annual' # MUST be annual
pot.threshold <- 0.99

source('decluster_timeseries.R') # routine for declustering time series

source('processing_norfolk.R') # yields data_norfolk

# write calibration data RDS file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_norfolk_decl',dt.decluster,'-pot',pot.threshold*100,'-',detrend.method,'_',today,'.rds', sep='')
saveRDS(data_norfolk, file=filename.output)

#===============================================================================
# End
#===============================================================================
