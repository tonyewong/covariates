#===============================================================================
# process_data.R
#
# Leads to:
#
# Starts with:  resulting data object from Wong et al 2018 (doi to follow...)
#               where detrending was an annual moving average, POT threshold of
#               99th percentile, and declustering time scale of 3 days
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

# read processed data from Wong et al 2018 (doi to follow...)
data_norfolk <- readRDS('../../EVT/data/tidegauge_processed_norfolk_decl3-pot99-annual_06Dec2017.rds')
data_delfzijl <- readRDS('../../EVT/data/tidegauge_processed_delfzijl_decl3-pot99-annual_20Dec2017.rds')

# determine the minimum block (# years) such that across both sites, each block
# will have at least one count in it
no_counts_norfolk <- which(data_norfolk$gpd89$counts==0)
no_counts_delfzijl <- which(data_delfzijl$gpd137$counts==0)

# start with blocks of size 2 years, increase from there if any
dt <- 2
l_gaps <- any(diff(no_counts_norfolk)==1) | any(diff(no_counts_delfzijl)==1)
while (l_gaps) {
  dt <- dt+1
  no_counts_norfolk <- diff(no_counts_norfolk)
  no_counts_delfzijl <- diff(no_counts_delfzijl)
  l_gaps <- any(diff(no_counts_norfolk)==0) | any(diff(no_counts_delfzijl)==0)
}

# now, dt is the shortest block of time (in years) for which every block will
# have at least one data point in it. shoudl find dt=3 years

# subsample smaller sets of the data
# starting from 11 years of data because that will allow an integer number of
# 3-year blocks for both Norfolk and Delfzijl

# norfolk lengths:
norfolk_total <- length(data_norfolk$gpd89$year)
norfolk_years <- c(seq(from=11, to=norfolk_total-dt, by=dt), norfolk_total)
len_norfolk <- NULL
for (yy in norfolk_years) {
  len_norfolk <- c(len_norfolk, paste('y', yy,sep=''))
}
names(norfolk_years) <- len_norfolk

# delfzijl lengths
delfzijl_total <- length(data_delfzijl$gpd137$year)
delfzijl_years <- c(seq(from=11, to=delfzijl_total-dt, by=dt), delfzijl_total)
len_delfzijl <- NULL
for (yy in delfzijl_years) {
  len_delfzijl <- c(len_delfzijl, paste('y', yy,sep=''))
}
names(delfzijl_years) <- len_delfzijl

data_gpd <- vector('list', 2); names(data_gpd) <- c('norfolk', 'delfzijl')

data_gpd$norfolk <- vector('list', length(len_norfolk)); names(data_gpd$norfolk) <- len_norfolk
for (len in len_norfolk) {
  data_gpd$norfolk[[len]] <- data_norfolk$gpd89  # to initialize everything with full data
  idx <- (length(data_norfolk$gpd89$year) - norfolk_years[[len]]+1):norfolk_total
  data_gpd$norfolk[[len]]$counts <- data_gpd$norfolk[[len]]$counts[idx]
  data_gpd$norfolk[[len]]$time_length <- data_gpd$norfolk[[len]]$time_length[idx]
  data_gpd$norfolk[[len]]$excesses <- data_gpd$norfolk[[len]]$excesses[idx]
  data_gpd$norfolk[[len]]$year <- data_gpd$norfolk[[len]]$year[idx]
  data_gpd$norfolk[[len]]$counts_all <- NULL
  data_gpd$norfolk[[len]]$time_length_all <- NULL
  data_gpd$norfolk[[len]]$excesses_all <- NULL
}

data_gpd$delfzijl <- vector('list', length(len_delfzijl)); names(data_gpd$delfzijl) <- len_delfzijl

for (len in len_delfzijl) {
  data_gpd$delfzijl[[len]] <- data_delfzijl$gpd137  # to initialize everything with full data
  idx <- (length(data_delfzijl$gpd137$year) - delfzijl_years[[len]]+1):delfzijl_total
  data_gpd$delfzijl[[len]]$counts <- data_gpd$delfzijl[[len]]$counts[idx]
  data_gpd$delfzijl[[len]]$time_length <- data_gpd$delfzijl[[len]]$time_length[idx]
  data_gpd$delfzijl[[len]]$excesses <- data_gpd$delfzijl[[len]]$excesses[idx]
  data_gpd$delfzijl[[len]]$year <- data_gpd$delfzijl[[len]]$year[idx]
  data_gpd$delfzijl[[len]]$counts_all <- NULL
  data_gpd$delfzijl[[len]]$time_length_all <- NULL
  data_gpd$delfzijl[[len]]$excesses_all <- NULL
}

# write calibration data RDS file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_norfolk-delfzijl_decl',data_norfolk$dt.decluster,'-pot',data_norfolk$gpd$p.threshold*100,'-',detrend.method,'_',today,'.rds', sep='')
saveRDS(data_gpd, file=filename.output)

#
#===============================================================================
# End
#===============================================================================
#
