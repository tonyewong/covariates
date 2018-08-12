#===============================================================================
# This file reads in the output files from the marginal likelihood estimator
# for each model and combines them into a combined .RDS file with bma weights
# and marginal likelihoods.
#
# Requires the RData files from the 'bridge_sample.R' routine as input
#
#===============================================================================
# Copyright 2017 Tony Wong
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

library(Hmisc)

appen <- '_threshold99'
path.ml <- paste('/home/scrim/axw322/codes/datalengths/output/bma',appen,sep='')
path.out <- '/home/scrim/axw322/codes/datalengths/output'
path.R <- '/home/scrim/axw322/codes/datalengths/R'

filename.likelihood <- paste('log_marginal_likelihood',appen,'.rds',sep='')
filename.weights <- paste('bma_weights',appen,'.rds',sep='')
filename.likelihood_all <- paste('log_marginal_likelihood_all',appen,'.rds',sep='')
filename.weights_all <- paste('bma_weights_all',appen,'.rds',sep='')

types.of.priors <- 'normalgamma'

site.names <- c('Norfolk')
n_sites <- length(site.names)
bma.weights <- vector('list', n_sites)
names(bma.weights) <- site.names

log.marg.lik <- vector('list', n_sites)
names(log.marg.lik) <- site.names


gpd.models <- c('gpd3','gpd4','gpd5','gpd6')
n_model <- length(gpd.models)

source(paste(path.R,'get_timeseries_covariates.R',sep='/'))

for (site in site.names) {

  bma.weights[[site]] <- vector('list', length(names_covariates))
  names(bma.weights[[site]]) <- names_covariates

  for (cc in names_covariates) {
    bma.weights[[site]][[cc]] <- rep(NA, n_model)
    names(bma.weights[[site]][[cc]]) <- gpd.models
  }
  log.marg.lik[[site]] <- vector('list', length(names_covariates))
  names(log.marg.lik[[site]]) <- names_covariates
  for (cc in names_covariates) {
      log.marg.lik[[site]][[cc]] <- rep(NA, n_model)
      names(log.marg.lik[[site]][[cc]]) <- gpd.models
  }
}

files <- list.files(path=path.ml, full.names=TRUE, recursive=FALSE)

for (file in files) {
  load(file)
  ##site <- capitalize(toString(station))
  station_name <- toString(station)
  site <- paste(toupper(substr(station_name, 1, 1)), substr(station_name, 2, nchar(station_name)), sep="")
  cc <- unlist(strsplit(file, split='_'))[4]
  # gpd.model is on the file, so just don't mess with it!
  log.marg.lik[[site]][[cc]][[gpd.model]] <- ml[length(ml)]
}

for (site in site.names) {
  for (cc in names_covariates) {
    ml <- log.marg.lik[[site]][[cc]]
    ml.scale <- ml - max(ml,na.rm=TRUE)
    for (model in gpd.models) {
      if (!is.na(log.marg.lik[[site]][[cc]][[model]])) {
        bma.weights[[site]][[cc]][[model]] <- exp(ml.scale[model])/sum(exp(ml.scale), na.rm=TRUE)
      }
    }
  }
}

saveRDS(log.marg.lik, paste(path.out,filename.likelihood, sep="/"))
saveRDS(bma.weights, paste(path.out,filename.weights, sep="/"))



#===============================================================================
#  Now, put ALL 16 of the candidate models into the mix. How are the weights
#  distributed among the various covariates?


for (site in site.names) {

  n_model_all <- length(names_covariates)*n_model
  bma.weights[[site]] <- vector('list', n_model_all)

  cnt <- 1
  for (cc in names_covariates) {
    for (model in gpd.models) {
      name_tmp <- paste(cc,model,sep="_")
      names(bma.weights[[site]])[cnt] <- name_tmp
      cnt <- cnt+1
    }
  }

  log.marg.lik[[site]] <- vector('list', n_model_all)
  names(log.marg.lik[[site]]) <- names(bma.weights[[site]])
}

files <- list.files(path=path.ml, full.names=TRUE, recursive=FALSE)

for (file in files) {
  load(file)
  ##site <- capitalize(toString(station))
  station_name <- toString(station)
  site <- paste(toupper(substr(station_name, 1, 1)), substr(station_name, 2, nchar(station_name)), sep="")
  cc <- unlist(strsplit(file, split='_'))[4]
  name_tmp <- paste(cc,gpd.model,sep="_")
  log.marg.lik[[site]][[name_tmp]] <- ml[length(ml)]
}

for (site in site.names) {
  ml <- unlist(log.marg.lik[[site]])
  ml.scale <- ml - max(ml,na.rm=TRUE)
  for (name_tmp in names(log.marg.lik[[site]])) {
    if (!is.na(log.marg.lik[[site]][[name_tmp]])) {
      bma.weights[[site]][[name_tmp]] <- exp(ml.scale[name_tmp])/sum(exp(ml.scale), na.rm=TRUE)
    }
  }
}

saveRDS(log.marg.lik, paste(path.out,filename.likelihood_all, sep="/"))
saveRDS(bma.weights, paste(path.out,filename.weights_all, sep="/"))


# a bit more exploratory analysis and draft plotting ===========================
if(FALSE) {

bw <- rev(sort(unlist(bma.weights)))
new_names <- NULL
for (i in 1:length(bw)) {
  cc_mod <- unlist(strsplit(names(bw)[i], split="[.]"))[3]
  cc_mod <- unlist(strsplit(cc_mod, split="[_]"))
  cc <- capitalize(cc_mod[1])
  if (cc=='Nao') {cc <- 'NAO'} else if (cc=='Sealevel') {cc <- 'Sea level'} else if (cc=='Temp') {cc <- 'Temperature'}
  mod <- cc_mod[2]
  if (mod=='gpd3') {mod <- 'ST'} else if (mod=='gpd4') {mod <- 'NS1'} else if (mod=='gpd5') {mod <- 'NS2'} else if (mod=='gpd6') {mod <- 'NS3'}
  new_names <- c(new_names, paste(cc,mod, sep=', '))
}
names(bw) <- new_names


# get totals for each covariate type:
bw_totals <- rep(0, length(names_covariates)); names(bw_totals) <- names_covariates
for (cc in names(bw)) {
  if (grepl('Time', cc)) {bw_totals['time'] <- bw_totals['time'] + bw[cc]
  } else if (grepl('Sea', cc)) {bw_totals['sea'] <- bw_totals['sea'] + bw[cc]
  } else if (grepl('Temp', cc)) {bw_totals['temp'] <- bw_totals['temp'] + bw[cc]
  } else if (grepl('NAO', cc)) {bw_totals['nao'] <- bw_totals['nao'] + bw[cc]
  }
}
#> bw_totals
#      time       temp   sealevel        nao
#0.09128698 0.29934667 0.43283316 0.17653319


par(las=1, mai=c(1,1.5,.2,.2))
barplot(bw, horiz=TRUE, names.arg=new_names, xlab='BMA weight')


#==========
# individual covariates
bma_weights <- readRDS('../output/bma_weights_threshold99.rds')

site <- 'Norfolk'
better_names <- c('ST','NS1','NS2','NS3')
covar_names <- c('Time', 'Temperature','Sea level', 'NAO'); names(covar_names) <- names_covariates

bw_cov <- vector('list', length(names_covariates)); names(bw_cov) <- names_covariates
for (cc in names_covariates) {
  bw_cov[[cc]] <- bma_weights[[site]][[cc]]
  names(bw_cov[[cc]]) <- better_names
}

par(mfrow=c(2,2), mai=c(.2,.5,.2,.2))
# Time
barplot(bw_cov$time, names.arg=better_names, ylab='BMA weight', space=1)
# Temperature
barplot(bw_cov$temp, names.arg=better_names, ylab='BMA weight', space=1)
# Sea level
barplot(bw_cov$sealevel, names.arg=better_names, ylab='BMA weight', space=1)
# NAO index
barplot(bw_cov$nao, names.arg=better_names, ylab='BMA weight', space=1)



} #=============================================================================


#===============================================================================
# End
#===============================================================================
