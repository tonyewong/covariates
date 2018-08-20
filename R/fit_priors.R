#===============================================================================
# fit_priors.R
#
# Fit prior distributions for four candidate PP/GPD-based model structures,
# times four candidate covariates, using all UHSLC data (research quality
# versions) at least 90 years, plus Sewells Point (Norfolk).
# Makes for 29 stations total.
#
# 1. get tide gauge data objects for all UHSLC database sites with > 90 years
# 2. get ... for Norfolk, VA, USA (not in database)
# 3. calculate maximum likelihood pp/gpd parameters, for each of the candidate
#    model structures, for each of the sites
# 4. fit normal or gamma prior distributions to these parameter sets, for each
#    model parameter within each of the candidate model structures.
# 5. write this priors object to a file (rds) and save progress to revisit later
#    (rdata)
#
# Updated 11 Dec 2017 // revised processing // tony wong
# Updated 18 Aug 2018 // revised for covariates // tony wong
#
# Questions? Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================

rm(list=ls())

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/covariates/R')
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/covariates/R')
}

pot.threshold <- 0.99   # POT threshold (percentile)
dt.decluster <- 3       # declustering time-scale (days)

.NP.deoptim <- 100      # number of DE population members (at least 10*[# parameters])
.niter.deoptim <- 100   # number of DE iterations
output.dir <- '../data/'
l.doprocessing <- 'FALSE'  # true if you need to run the processing
                          # false -> read in some previous RDS processing results,
                          # with the filenames defined below

appen <- paste('ppgpd_decl',dt.decluster,'-pot',100*pot.threshold,sep='')

if(pot.threshold==0.99 & dt.decluster==3) {
  filename.many <- '../data/tidegauge_processed_manystations_decl3-pot99-annual_10Dec2017.rds'
  filename.norfolk <- '../data/tidegauge_processed_norfolk_decl3-pot99-annual_13Aug2018.rds'
} else {
  print('ERROR: unknown input for pot.threshold and/or dt.decluster')
}

#
#===============================================================================
# relevant libraries - do 'install.packages([library name])' if you do not have
# one yet
#===============================================================================
#

library(extRemes)
library(zoo)
library(adaptMCMC)
library(lhs)
library(DEoptim)
library(ncdf4)

#
#===============================================================================
# read and process data for forcing (auxiliary covariate for nonstationary
# parameters).
#===============================================================================
#

print('reading forcing data...')

# Load the covariate names and values
source('get_timeseries_covariates_priors.R')

print('...done.')

#
#===============================================================================
# read and process data for tide gauge stations (or read previous processing)
#===============================================================================
#

print('reading/processing data from tide gauge stations...')
print(' (if you do not have this already done, this might take a while)...')

if(l.doprocessing) {
  source('process_data.R')
} else {
  data_many <- readRDS(filename.many)
  data_norfolk <- readRDS(filename.norfolk)
}

# round them all up as one big data set
data_all <- data_many
data_all$norfolk <- data_norfolk$y89

print('...done.')

#
#===============================================================================
# set up PP-GPD model parameters
#===============================================================================
#

print('setting up PP-GPD model parameters for DE optimization...')

source('parameter_setup.R')

print('...done.')

#
#===============================================================================
# parameters for DE optim (for maximum likelihood/minimum negative likelihoood)
#===============================================================================
#

NP.deoptim <- .NP.deoptim
niter.deoptim <- .niter.deoptim
F.deoptim <- 0.8
CR.deoptim <- 0.9

#
#===============================================================================
# fit MLE PP-GPD model parameters for each candidate model at each tide gauge
#===============================================================================
#

print('starting DE optimization for MLE PP-GPD parameters for all stations in set...')

# need the likelihood function
source('likelihood_ppgpd.R')

deoptim.all <- vector('list', length(names_covariates))
names(deoptim.all) <- names_covariates
for (cc in names_covariates) {
  deoptim.all[[cc]] <- vector('list', nmodel); names(deoptim.all[[cc]]) <- types.of.model
  for (i in 1:nmodel) {
    deoptim.all[[cc]][[types.of.model[i]]] <- mat.or.vec(length(data_all), length(parnames_all[[types.of.model[i]]]))
    rownames(deoptim.all[[cc]][[types.of.model[i]]]) <- names(data_all)
    colnames(deoptim.all[[cc]][[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
  }
}

for (dd in 1:(length(data_all)-1)) {

  print(paste('starting to calculate MLE PP-GPD parameters for tide gauge data set ',dd,' / ',length(data_all),sep=''))
  tbeg0 <- proc.time()
  data_all[[dd]]$deoptim <- vector('list', nmodel)
  names(data_all[[dd]]$deoptim) <- types.of.model

  for (gpd.type in types.of.gpd) {

    print(paste('  - starting DE optimization for model',gpd.type,'...'))
    tbeg <- proc.time()

    irem_aux <- NULL

    # if tide gauge record starts before auxiliary forcing, clip it
    if(data_all[[dd]]$gpd$year[1] < covariates[1,'year']) {
      irem <- which(data_all[[dd]]$gpd$year < covariates[,'year'][1])
      data_all[[dd]]$gev_year$year <- data_all[[dd]]$gev_year$year[-irem]
      data_all[[dd]]$gev_year$lsl_max <- data_all[[dd]]$gev_year$lsl_max[-irem]
      data_all[[dd]]$gpd$year <- data_all[[dd]]$gpd$year[-irem]
      data_all[[dd]]$gpd$counts <- data_all[[dd]]$gpd$counts[-irem]
      data_all[[dd]]$gpd$excesses <- data_all[[dd]]$gpd$excesses[-irem]
      data_all[[dd]]$gpd$time_length <- data_all[[dd]]$gpd$time_length[-irem]
      data_all[[dd]]$gpd$time_length_all <- sum(data_all[[dd]]$gpd$time_length)
      data_all[[dd]]$gpd$counts_all <- sum(unlist(data_all[[dd]]$gpd$counts), na.rm=TRUE)
      data_all[[dd]]$gpd$excesses_all <- unlist(data_all[[dd]]$gpd$excesses)[!is.na(unlist(data_all[[dd]]$gpd$excesses))]
    } else if(data_all[[dd]]$gpd$year[1] > covariates[1,'year']) {
    # if begins after the forcing, clip the forcing
      irem_aux <- c(irem_aux, which(covariates[,'year'] < data_all[[dd]]$gpd$year[1]))
    }

    # if tide gauge record ends after auxiliary forcing, clip it
    if(max(data_all[[dd]]$gpd$year) > max(covariates[,'year'])) {
      irem <- which(data_all[[dd]]$gpd$year > max(covariates[,'year']))
      data_all[[dd]]$gev_year$year <- data_all[[dd]]$gev_year$year[-irem]
      data_all[[dd]]$gev_year$lsl_max <- data_all[[dd]]$gev_year$lsl_max[-irem]
      data_all[[dd]]$gpd$year <- data_all[[dd]]$gpd$year[-irem]
      data_all[[dd]]$gpd$counts <- data_all[[dd]]$gpd$counts[-irem]
      data_all[[dd]]$gpd$excesses <- data_all[[dd]]$gpd$excesses[-irem]
      data_all[[dd]]$gpd$time_length <- data_all[[dd]]$gpd$time_length[-irem]
      data_all[[dd]]$gpd$time_length_all <- sum(data_all[[dd]]$gpd$time_length)
      data_all[[dd]]$gpd$counts_all <- sum(unlist(data_all[[dd]]$gpd$counts), na.rm=TRUE)
      data_all[[dd]]$gpd$excesses_all <- unlist(data_all[[dd]]$gpd$excesses)[!is.na(unlist(data_all[[dd]]$gpd$excesses))]
    } else if(max(data_all[[dd]]$gpd$year) < max(covariates[,'year'])) {
    # if ends before forcing, clip the forcing
      irem_aux <- c(irem_aux, which(covariates[,'year'] > max(data_all[[dd]]$gpd$year)))
    }

    covariates_trimmed <- covariates
    if(length(irem_aux) > 0) {covariates_trimmed <- covariates_trimmed[-irem_aux,]}

    for (cc in names_covariates) {

      auxiliary_in <- covariates_trimmed[,cc]
      forc_max <- max(auxiliary_in)

      if(gpd.type=='gpd3') {
        auxiliary <- NULL
      } else {
        auxiliary <- auxiliary_in
      }

      out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=bound_lower_set[[gpd.type]], upper=bound_upper_set[[gpd.type]],
                           DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                           parnames=parnames_all[[gpd.type]], data_calib=data_all[[dd]]$gpd, auxiliary=auxiliary)
      deoptim.all[[cc]][[gpd.type]][dd,] <- out.deoptim$optim$bestmem
      colnames(deoptim.all[[cc]][[gpd.type]]) <- parnames_all[[gpd.type]]

    }
    tend <- proc.time()
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
  tend0 <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend0-tbeg0)[3]/60,2),' minutes', sep=''))
}

dd <- length(data_all) # Norfolk is set up differently

print(paste('starting to calculate MLE PP-GPD parameters for tide gauge data set ',dd,' / ',length(data_all),sep=''))
tbeg0 <- proc.time()
data_all[[dd]]$deoptim <- vector('list', nmodel)
names(data_all[[dd]]$deoptim) <- types.of.model

for (gpd.type in types.of.gpd) {

  print(paste('  - starting DE optimization for model',gpd.type,'...'))
  tbeg <- proc.time()

  irem_aux <- NULL

  # if tide gauge record starts before auxiliary forcing, clip it
  if(data_all[[dd]]$year[1] < covariates[1,'year']) {
    irem <- which(data_all[[dd]]$year < covariates[,'year'][1])
    data_all[[dd]]$year <- data_all[[dd]]$year[-irem]
    data_all[[dd]]$counts <- data_all[[dd]]$counts[-irem]
    data_all[[dd]]$excesses <- data_all[[dd]]$excesses[-irem]
    data_all[[dd]]$time_length <- data_all[[dd]]$time_length[-irem]
    data_all[[dd]]$time_length_all <- sum(data_all[[dd]]$time_length)
    data_all[[dd]]$counts_all <- sum(unlist(data_all[[dd]]$counts), na.rm=TRUE)
    data_all[[dd]]$excesses_all <- unlist(data_all[[dd]]$excesses)[!is.na(unlist(data_all[[dd]]$excesses))]
  } else if(data_all[[dd]]$year[1] > covariates[1,'year']) {
  # if begins after the forcing, clip the forcing
    irem_aux <- c(irem_aux, which(covariates[,'year'] < data_all[[dd]]$year[1]))
  }

  # if tide gauge record ends after auxiliary forcing, clip it
  if(max(data_all[[dd]]$year) > max(covariates[,'year'])) {
    irem <- which(data_all[[dd]]$year > max(covariates[,'year']))
    data_all[[dd]]$year <- data_all[[dd]]$year[-irem]
    data_all[[dd]]$counts <- data_all[[dd]]$counts[-irem]
    data_all[[dd]]$excesses <- data_all[[dd]]$excesses[-irem]
    data_all[[dd]]$time_length <- data_all[[dd]]$time_length[-irem]
    data_all[[dd]]$time_length_all <- sum(data_all[[dd]]$time_length)
    data_all[[dd]]$counts_all <- sum(unlist(data_all[[dd]]$counts), na.rm=TRUE)
    data_all[[dd]]$excesses_all <- unlist(data_all[[dd]]$excesses)[!is.na(unlist(data_all[[dd]]$excesses))]
  } else if(max(data_all[[dd]]$year) < max(covariates[,'year'])) {
  # if ends before forcing, clip the forcing
    irem_aux <- c(irem_aux, which(covariates[,'year'] > max(data_all[[dd]]$year)))
  }

  covariates_trimmed <- covariates
  if(length(irem_aux) > 0) {covariates_trimmed <- covariates_trimmed[-irem_aux,]}

  for (cc in names_covariates) {

    auxiliary_in <- covariates_trimmed[,cc]
    forc_max <- max(auxiliary_in)

    if(gpd.type=='gpd3') {
      auxiliary <- NULL
    } else {
      auxiliary <- auxiliary_in
    }

    out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=bound_lower_set[[gpd.type]], upper=bound_upper_set[[gpd.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[gpd.type]], data_calib=data_all[[dd]], auxiliary=auxiliary)
    deoptim.all[[cc]][[gpd.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.all[[cc]][[gpd.type]]) <- parnames_all[[gpd.type]]

  }
  tend <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
}
tend0 <- proc.time()
print(paste('... done. Took ',round(as.numeric(tend0-tbeg0)[3]/60,2),' minutes', sep=''))

print('...done.')

#
#===============================================================================
# check distributions, fit priors
#===============================================================================
#

print('fitting prior distributions to the MLE parameters...')

# fit gamma and normal priors
# -> centered at the medians
# -> with standard deviation equal to half the max-min range
#    (or do empirical sd? might underestimate though - take wider)

# assign which parameters have which priors
if(exists('gamma.priors')) {rm(list=c('gamma.priors','normal.priors','uniform.priors'))}
gamma.priors <- c('lambda','lambda0','sigma','sigma0')
normal.priors <- c('lambda1','sigma1','xi','xi0','xi1')
uniform.priors <- NULL

priors_normalgamma <- vector('list', length(names_covariates))
names(priors_normalgamma) <- names_covariates
for (cc in names_covariates) {
  priors_normalgamma[[cc]] <- vector('list', nmodel); names(priors_normalgamma[[cc]]) <- types.of.model
  for (model in types.of.model) {
    priors_normalgamma[[cc]][[model]] <- vector('list', length(parnames_all[[model]]))
    names(priors_normalgamma[[cc]][[model]]) <- parnames_all[[model]]
    for (par in parnames_all[[model]]) {
      priors_normalgamma[[cc]][[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
      if(!is.na(match(par, uniform.priors))) {
         names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','shape','rate')
         priors_normalgamma[[cc]][[model]][[par]]$type <- 'uniform'
         priors_normalgamma[[cc]][[model]][[par]]$lower <- bound_lower_set[[model]][match(par,parnames_all[[model]])]
         priors_normalgamma[[cc]][[model]][[par]]$upper <- bound_upper_set[[model]][match(par,parnames_all[[model]])]
      } else if(!is.na(match(par, gamma.priors))) { # shape=alpha, rate=beta, mean=shape/rate, var=shape/rate^2
        names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','shape','rate')
        priors_normalgamma[[cc]][[model]][[par]]$type <- 'gamma'
##        priors_normalgamma[[cc]][[model]][[par]]$rate <- median(deoptim.all[[cc]][[model]][,par]) / (0.5*(max(deoptim.all[[cc]][[model]][,par])-min(deoptim.all[[cc]][[model]][,par])))^2
##        priors_normalgamma[[cc]][[model]][[par]]$shape <- median(deoptim.all[[cc]][[model]][,par]) * priors_normalgamma[[cc]][[model]][[par]]$rate
        priors_normalgamma[[cc]][[model]][[par]]$rate <- mean(deoptim.all[[cc]][[model]][,par]) / var(deoptim.all[[cc]][[model]][,par])
        priors_normalgamma[[cc]][[model]][[par]]$shape <- mean(deoptim.all[[cc]][[model]][,par]) * priors_normalgamma[[cc]][[model]][[par]]$rate
      } else if(!is.na(match(par, normal.priors))) {
        names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','mean','sd')
        priors_normalgamma[[cc]][[model]][[par]]$type <- 'normal'
##        priors_normalgamma[[cc]][[model]][[par]]$mean <- median(deoptim.all[[cc]][[model]][,par])
##        priors_normalgamma[[cc]][[model]][[par]]$sd   <- 0.5*(max(deoptim.all[[cc]][[model]][,par])-min(deoptim.all[[cc]][[model]][,par]))
        priors_normalgamma[[cc]][[model]][[par]]$mean <- mean(deoptim.all[[cc]][[model]][,par])
        priors_normalgamma[[cc]][[model]][[par]]$sd   <- sd(deoptim.all[[cc]][[model]][,par])
      }
    }
  }
}

print('...done.')


#
#===============================================================================
# "fit" wide uniform priors (just using the bounds for the DE optim search)
#===============================================================================
#

print('fitting prior distributions to the uniform bounds for MLE search...')

# all parameters have uniform bounds, given by bound_lower_set and bound_upper_set
priors_uniform <- vector('list', length(names_covariates)); names(priors_uniform) <- names_covariates
for (cc in names_covariates) {
  priors_uniform[[cc]] <- vector('list', nmodel); names(priors_uniform[[cc]]) <- types.of.model
  for (model in types.of.model) {
    priors_uniform[[cc]][[model]] <- vector('list', length(parnames_all[[model]])); names(priors_uniform[[cc]][[model]]) <- parnames_all[[model]]
    for (par in parnames_all[[model]]) {
      priors_uniform[[cc]][[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
      names(priors_uniform[[cc]][[model]][[par]]) <- c('type','lower','upper'); priors_uniform[[cc]][[model]][[par]]$type <- 'uniform'
      priors_uniform[[cc]][[model]][[par]]$lower <- bound_lower_set[[model]][match(par,parnames_all[[model]])]
      priors_uniform[[cc]][[model]][[par]]$upper <- bound_upper_set[[model]][match(par,parnames_all[[model]])]
    }
  }
}

print('...done.')

#
#===============================================================================
# save priors and initial values (from DE optim for Delfzijl, Balboa, and Norfolk)
# and read later when calibrating with MCMC
#===============================================================================
#

# rds -> save single object; the only one we need is 'priors'
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.priors.normalgamma <- paste(output.dir,'surge_priors_normalgamma_',appen,'_',today,'.rds', sep='')
filename.priors.uniform <- paste(output.dir,'surge_priors_uniform_',appen,'_',today,'.rds', sep='')
filename.mles <- paste(output.dir,'surge_MLEs_',appen,'_',today,'.rds', sep='')

print(paste('saving priors and DE optim output as .rds files to read and use later...',sep=''))

saveRDS(priors_normalgamma, file=filename.priors.normalgamma)
saveRDS(priors_uniform, file=filename.priors.uniform)
saveRDS(deoptim.all, file=filename.mles)


print('...done.')

#
#===============================================================================
# End
#===============================================================================
#

if(FALSE) {

deoptim.all <- readRDS('../data/surge_MLEs_ppgpd_decl3-pot99_20Aug2018.rds')
priors_normalgamma <- readRDS('../data/surge_priors_normalgamma_ppgpd_decl3-pot99_20Aug2018.rds')

plot_priors <- function(cc, appen='') {

  nbins <- 12 # note that there are only 30 sites...
  frac.ran <- 0.35 # extend axes above/below max/min range

  # get standard limits for each of the 6 parameters
  lims <- mat.or.vec(6,2)
  pp <- 1 # lambda0
    parameters.pooled <- c(deoptim.all[[cc]]$gpd3[,1], deoptim.all[[cc]]$gpd4[,1],
                           deoptim.all[[cc]]$gpd5[,1], deoptim.all[[cc]]$gpd6[,1])
    ran <- diff(quantile(parameters.pooled, c(0,1)))
    lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
  pp <- 2 # lambda1
    parameters.pooled <- c(deoptim.all[[cc]]$gpd4[,2], deoptim.all[[cc]]$gpd5[,2], deoptim.all[[cc]]$gpd6[,2])
    ran <- diff(quantile(parameters.pooled, c(0,1)))
    lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
  pp <- 3 # sigma0
    parameters.pooled <- c(log(deoptim.all[[cc]]$gpd3[,2]), log(deoptim.all[[cc]]$gpd4[,3]),
                            deoptim.all[[cc]]$gpd5[,3], deoptim.all[[cc]]$gpd6[,3])
    ran <- diff(quantile(parameters.pooled, c(0,1)))
    lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
  pp <- 4 # sigma1
    parameters.pooled <- c(deoptim.all[[cc]]$gpd5[,4], deoptim.all[[cc]]$gpd6[,4])
    ran <- diff(quantile(parameters.pooled, c(0,1)))
    lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
  pp <- 5 # xi0
    parameters.pooled <- c(deoptim.all[[cc]]$gpd5[,5], deoptim.all[[cc]]$gpd6[,5])
    ran <- diff(quantile(parameters.pooled, c(0,1)))
    lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
  pp <- 6 # xi1
    parameters.pooled <- c(deoptim.all[[cc]]$gpd6[,6])
    ran <- diff(quantile(parameters.pooled, c(0,1)))
    lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)

  if (nchar(appen)==0) {
    plotname <- paste('../figures/priors_normalgamma_',cc,'.pdf',sep='')
  } else {
    plotname <- paste('../figures/priors_normalgamma_',cc,'_',appen,'.pdf',sep='')
  }
  pdf(plotname, height=7, width=10, colormodel='cmyk')
  par(mfrow=c(4,6))
  ##=============================
  model <- 'gpd3'
  pp <- 1; par(mai=c(.25,.59,.25,.01))
  box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000)
  lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[cc]][[model]]$lambda$shape,
                      rate=priors_normalgamma[[cc]][[model]]$lambda$rate), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext('Probability density', side=2, line=1.2, cex=1)
  mtext('ST', side=2, line=3, cex=1)
  plot.new()
  pp <- 2; par(mai=c(.25,.3,.25,.3))
  box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
  hist(log(deoptim.all[[cc]][[model]][,pp]), xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  rate.tmp <- median(log(deoptim.all[[cc]][[model]][,2])) / (0.5*(max(log(deoptim.all[[cc]][[model]][,2]))-min(log(deoptim.all[[cc]][[model]][,2]))))^2
  shape.tmp <- median(log(deoptim.all[[cc]][[model]][,2])) * rate.tmp
  x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=shape.tmp, rate=rate.tmp), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  plot.new()
  pp <- 3; par(mai=c(.25,.3,.25,.3))
  box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$xi$mean, sd=priors_normalgamma[[cc]][[model]]$xi$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  plot.new()
  mtext(cc, side=3)
  ##=============================
  model <- 'gpd4'
  pp <- 1; par(mai=c(.3,.59,.2,.01))
  box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000)
  lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[cc]][[model]]$lambda0$shape,
                               rate=priors_normalgamma[[cc]][[model]]$lambda0$rate), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext('Probability density', side=2, line=1.2, cex=1)
  mtext('NS1', side=2, line=3, cex=1)
  pp <- 2; par(mai=c(.3,.35,.2,.25))
  box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$lambda1$mean,
                     sd=priors_normalgamma[[cc]][[model]]$lambda1$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  pp <- 3; par(mai=c(.3,.3,.2,.3))
  box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
  hist(log(deoptim.all[[cc]][[model]][,pp]), xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  rate.tmp <- median(log(deoptim.all[[cc]][[model]][,3])) / (0.5*(max(log(deoptim.all[[cc]][[model]][,3]))-min(log(deoptim.all[[cc]][[model]][,3]))))^2
  shape.tmp <- median(log(deoptim.all[[cc]][[model]][,3])) * rate.tmp
  x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=shape.tmp, rate=rate.tmp), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  plot.new()
  pp <- 4; par(mai=c(.3,.3,.2,.3))
  box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$xi$mean,
                     sd=priors_normalgamma[[cc]][[model]]$xi$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  plot.new()
  ##=============================
  model <- 'gpd5'
  pp <- 1; par(mai=c(.35,.59,.15,.01))
  box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000)
  lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[cc]][[model]]$lambda0$shape,
                      rate=priors_normalgamma[[cc]][[model]]$lambda0$rate), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext('Probability density', side=2, line=1.2, cex=1)
  mtext('NS2', side=2, line=3, cex=1)
  pp <- 2; par(mai=c(.35,.35,.15,.25))
  box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$lambda1$mean,
                     sd=priors_normalgamma[[cc]][[model]]$lambda1$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  pp <- 3; par(mai=c(.35,.3,.15,.3))
  box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000)
  lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[cc]][[model]]$sigma0$shape,
                      rate=priors_normalgamma[[cc]][[model]]$sigma0$rate), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  pp <- 4; par(mai=c(.35,.3,.15,.3))
  box.width <- diff(lims[4,])/nbins; box.edges <- seq(from=lims[4,1], to=lims[4,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[4,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[4,1], to=lims[4,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$sigma1$mean,
                     sd=priors_normalgamma[[cc]][[model]]$sigma1$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  pp <- 5; par(mai=c(.35,.3,.15,.3))
  box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$xi$mean,
                     sd=priors_normalgamma[[cc]][[model]]$xi$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  plot.new()
  ##=============================
  model <- 'gpd6'
  pp <- 1; par(mai=c(.5,.59,.01,.01))
  box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000)
  lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[cc]][[model]]$lambda0$shape,
                      rate=priors_normalgamma[[cc]][[model]]$lambda0$rate), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext('Probability density', side=2, line=1.2, cex=1)
  mtext(expression(lambda[0]), side=1, line=2.7, cex=1)
  mtext('NS3', side=2, line=3, cex=1)
  pp <- 2; par(mai=c(.5,.35,.01,.25))
  box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$lambda1$mean,
                     sd=priors_normalgamma[[cc]][[model]]$lambda1$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext(expression(lambda[1]), side=1, line=2.7, cex=1)
  pp <- 3; par(mai=c(.5,.3,.01,.3))
  box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000)
  lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[cc]][[model]]$sigma0$shape,
                      rate=priors_normalgamma[[cc]][[model]]$sigma0$rate), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext(expression(sigma[0]), side=1, line=2.7, cex=1)
  pp <- 4; par(mai=c(.5,.3,.01,.3))
  box.width <- diff(lims[4,])/nbins; box.edges <- seq(from=lims[4,1], to=lims[4,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[4,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[4,1], to=lims[4,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$sigma1$mean,
                     sd=priors_normalgamma[[cc]][[model]]$sigma1$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext(expression(sigma[1]), side=1, line=2.7, cex=1)
  pp <- 5; par(mai=c(.5,.3,.01,.3))
  box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$xi0$mean,
                     sd=priors_normalgamma[[cc]][[model]]$xi0$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext(expression(xi[0]), side=1, line=2.7, cex=1)
  pp <- 6; par(mai=c(.5,.3,.01,.3))
  box.width <- diff(lims[6,])/nbins; box.edges <- seq(from=lims[6,1], to=lims[6,2], by=box.width)
  hist(deoptim.all[[cc]][[model]][,pp], xlim=lims[6,], freq=FALSE, main='', xlab='', ylab='',
       breaks=box.edges, yaxt='n', yaxs='i')
  x.tmp <- seq(from=lims[6,1], to=lims[6,2], length.out=1000)
  lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[cc]][[model]]$xi1$mean,
                     sd=priors_normalgamma[[cc]][[model]]$xi1$sd), col='red', lwd=2)
  u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
  mtext(expression(xi[1]), side=1, line=2.7, cex=1)
  ##=============================
  dev.off()
}

# make plots of the prior distirbutions for each covariate
for (cc in names_covariates) {plot_priors(cc)}

}
