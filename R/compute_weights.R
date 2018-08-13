#===============================================================================
# compute_weights.R
#
# This file reads in the output files from the marginal likelihood estimator
# for each model and combines them into a combined .RDS file with bma weights
# and marginal likelihoods.
#
# Requires the RData files from the 'bridge_sample.R' routine as input
#
# Original code: Vivek Srikrishnan (Penn State) 2017
# Modified code: Tony Wong (CU Boulder) 2018
#===============================================================================

rm(list=ls())

library(Hmisc)

appen <- '_threshold99'

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/covariates/R')
  path.ml <- paste('/Users/tony/codes/covariates/output/bma',appen,sep='')
  path.out <- '/Users/tony/codes/covariates/output'
  path.R <- '/Users/tony/codes/covariates/R'
  nnode_mcmc_prod000 <- 1          # number of CPUs to use (PRODUCTION chains)
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/covariates/R')
  path.ml <- paste('/home/scrim/axw322/codes/covariates/output/bma',appen,sep='')
  path.out <- '/home/scrim/axw322/codes/covariates/output'
  path.R <- '/home/scrim/axw322/codes/covariates/R'
  nnode_mcmc_prod000 <- 10          # number of CPUs to use (PRODUCTION chains)
}


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


site <- "Norfolk"

# n_model_all corrected for multiple counting of stationary model (gpd3, or ST)
n_model_all <- length(names_covariates)*n_model - (length(names_covariates)-1)
bma.weights[[site]] <- vector('list', n_model_all)

cnt <- 1
for (cc in names_covariates) {
  for (model in gpd.models) {
    if (model=='gpd3') {
      if (cc=='time') {
        name_tmp <- paste(cc,model,sep="_")
        names(bma.weights[[site]])[cnt] <- name_tmp
        cnt <- cnt+1
      }
    } else {
      name_tmp <- paste(cc,model,sep="_")
      names(bma.weights[[site]])[cnt] <- name_tmp
      cnt <- cnt+1
    }
  }
}

log.marg.lik[[site]] <- vector('list', n_model_all)
names(log.marg.lik[[site]]) <- names(bma.weights[[site]])

files <- list.files(path=path.ml, full.names=TRUE, recursive=FALSE)

for (file in files) {
  load(file)
  ##site <- capitalize(toString(station))
  station_name <- toString(station)
  site <- paste(toupper(substr(station_name, 1, 1)), substr(station_name, 2, nchar(station_name)), sep="")
  cc <- unlist(strsplit(file, split='_'))[4]
  name_tmp <- paste(cc,gpd.model,sep="_")
  if (name_tmp %in% names(log.marg.lik[[site]])) {log.marg.lik[[site]][[name_tmp]] <- ml[length(ml)]}
}

ml <- unlist(log.marg.lik[[site]])
ml.scale <- ml - max(ml,na.rm=TRUE)
for (name_tmp in names(log.marg.lik[[site]])) {
  if (!is.na(log.marg.lik[[site]][[name_tmp]])) {
    bma.weights[[site]][[name_tmp]] <- exp(ml.scale[name_tmp])/sum(exp(ml.scale), na.rm=TRUE)
  }
}

saveRDS(log.marg.lik, paste(path.out,filename.likelihood_all, sep="/"))
saveRDS(bma.weights, paste(path.out,filename.weights_all, sep="/"))

#===============================================================================
# End
#===============================================================================
