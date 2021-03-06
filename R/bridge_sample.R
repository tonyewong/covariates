#===============================================================================
# This file contains the working script for estimating the          #
# marginal likelihood of different GPD surge models and data        #
# experiments.                                                      #
#                                                                   #
# This revised version calculates all likelihood estimates as a parallel
# 'foreach' loop.
#
# Original code: Vivek Srikrishnan (Penn State) 2017
# Modified code: Tony Wong (CU Boulder) 2018
#===============================================================================

# import libraries
library(mvtnorm)  # used for multivariate normal samples and densities
#library(coda)     # used to compute spectral density of MCMC draws at frequency 0
library(extRemes) # used to compute GEV densities within posterior likelihood function
library(foreach)
library(doParallel)
library(ncdf4)

appen <- '_threshold99'
calib_date <- '24Aug2018'
dt.decluster <- 3
pot.threshold <- 0.99
type.of.priors <- 'normalgamma'     # can be 'uniform' or 'normalgamma'

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/covariates/R')
  path.R <- '/Users/tony/codes/covariates/R'
  # set data and save directories
  path.data <- '/Users/tony/codes/covariates/data'
  path.output <- '/Users/tony/codes/covariates/output'
  path.save <- paste('/Users/tony/codes/covariates/output/bma',appen,sep='')
  nnode <- 4          # number of CPUs to use
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/covariates/R')
  path.R <- '/home/scrim/axw322/codes/covariates/R'
  path.data <- '/home/scrim/axw322/codes/covariates/data'
  path.output <- '/home/scrim/axw322/codes/covariates/output'
  path.save <- '/home/scrim/axw322/codes/covariates/output/bma'
  path.save <- paste('/home/scrim/axw322/codes/covariates/output/bma',appen,sep='')
  nnode <- 10          # number of CPUs to use
}

# import file containing the log likelihood calculations
source(paste(path.R,'likelihood_ppgpd.R',sep='/'))

source(paste(path.R,'get_timeseries_covariates.R',sep='/'))

# set up table of experiment parameters
# some of these combinations (in terms of the data length and sites) don't exist, but they
# just throw errors anyway. There is probably a better way to construct this data frame.
# could include the type of prior in this as well

#sites <- c('delfzijl','norfolk')
sites <- c('norfolk')
n_sites <- length(sites)
n_covar <- length(names_covariates)
gpd.models <- c('gpd3','gpd4','gpd5','gpd6')
n_model <- length(gpd.models)

# list to store the actual output
output <- vector('list', n_covar); names(output) <- names_covariates
for (cc in names_covariates) {
  output[[cc]] <- vector('list', n_model)
  names(output[[cc]]) <- gpd.models
}

# data frame to store the experimental details
experiments <- expand.grid(station  =sites,
                           covariate=names_covariates,
                           gpd.model=gpd.models)
#irem <- NULL
#irem <- c(irem, which((experiments[,'data.length']==89) & experiments[,'station']=='delfzijl'))
#irem <- c(irem, which((experiments[,'data.length']==137 | experiments[,'data.length']==110 | experiments[,'data.length']==90) & experiments[,'station']=='norfolk'))
#experiments <- experiments[-irem,]
n_experiments <- nrow(experiments)
output <- vector('list', n_experiments)

#cores = detectCores()
#cl <- makeCluster(cores[1]-1) #not to overload your computer
cl <- makeCluster(nnode)
print(paste('Starting cluster with ',nnode,' cores', sep=''))
registerDoParallel(cl)

source('bridge_sample_functions.R')
export.names <- c('bridge.samp.rel.err','bridge.samp.iter','recip.imp.samp','experiments','log_post_ppgpd','log_like_ppgpd','log_prior_ppgpd','path.R','calib_date','dt.decluster','pot.threshold','type.of.priors')

finalOutput <- foreach(ee=1:n_experiments,
                            .packages=c('mvtnorm','extRemes','ncdf4'),
                            .export=export.names,
                            .inorder=FALSE) %dopar% {

  setwd(path.R)
  source(paste(path.R,'likelihood_ppgpd.R',sep='/'))
  source(paste(path.R,'get_timeseries_covariates.R',sep='/'))
  # get parameters for this particular experiment
  print(experiments[ee,])
  station <- experiments[ee,'station']
  gpd.model <- experiments[ee,'gpd.model']
  covariate <- experiments[ee,'covariate']

  # set output (saved as .RData; to be collected into a single output file later) file name
  filename.out <- paste('ml_',station,'_',covariate,'_',gpd.model,'.RData',sep='')
  if (file.exists(paste(path.save, filename.out, sep='/'))) {
     #stop('Output file already exists!')
     print('Output file already exists!')
     output[[ee]] <- 'done!'
 } else {

  # read in calibration output file
  print('loading calibration file...')

  setwd(path.data)

  ##type.of.priors <- 'normalgamma'     # can be 'uniform' or 'normalgamma'
  filename.priors <- Sys.glob(paste('surge_priors_',type.of.priors,'_ppgpd_decl',dt.decluster,'-pot',pot.threshold*100,'_*','.rds',sep='')) # is in the output directory

  # use this if multiple files exist for the same location and prior
  #calib_date <- '31Dec2017'
  #pot.threshold <- '_pot95' # for main text results (_pot99)
  priors <- readRDS(filename.priors)
  if (exists('calib_date')) {
    setwd(path.output)
    filename.calib <- paste('mcmc_ppgpd-experiments_',station,'_',type.of.priors,'_decl',dt.decluster,'-pot',pot.threshold*100,'_',calib_date,'.RData',sep='')
  } else {
    setwd(path.output)
    filename.calib <- Sys.glob(paste('mcmc_ppgpd-experiments_',station,'_',type.of.priors,'_decl',dt.decluster,'-pot',pot.threshold*100,'_*','.RData',sep=''))
  }
  load(filename.calib)

  print('done!')

  # set name for data experiment from data.length variable
  cc <- paste(covariate)

  # compute the burn-in length for sampling based on the ifirst object
  full.length <- nrow(amcmc_out[[cc]][[gpd.model]][[1]]$samples)
  ##burned.length <- nrow(chains_burned[[gpd.exp]][[gpd.model]][[1]])
  burn.in <- ifirst[[cc]]

  # set number of samples to use for estimate
  post.samp.num <- 50000
  imp.samp.num <- 50000

  # burn in samples and log.p values
  post.samples <- amcmc_out[[cc]][[gpd.model]][[1]]$samples
  post.samples <- post.samples[(burn.in+1):full.length,]
  post.ll <- amcmc_out[[cc]][[gpd.model]][[1]]$log.p
  post.ll <- post.ll[(burn.in+1):full.length]

  # fit normal approximation to the posteriors
  post.mean <- colMeans(post.samples)
  post.cov <- cov(post.samples)

  # get posterior samples
  print('sampling from posterior distribution...')

  samp.names <- c('samples','log.imp','log.p')
  post.samp <- setNames(vector("list",length(samp.names)),samp.names)
  samp.idx <- sample(x=nrow(post.samples), size=post.samp.num, replace=TRUE)
  post.samp$samples <- post.samples
  post.samp$samples <- post.samples[samp.idx,]
  # get posterior log-likelihood of sampled posterior values
  post.samp$log.p <- post.ll
  post.samp$log.p <- post.ll[samp.idx]
  # get importance log-likelhood of posterior samples
  post.samp$log.imp <- dmvnorm(x=post.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)

  print('done!')

  # get importance samples and likelihood
  print('sampling from importance distribution...')

  imp.samp <- setNames(vector("list",length(samp.names)),samp.names)
  imp.samp$samples <- rmvnorm(n=imp.samp.num, mean=post.mean, sigma=post.cov)
  imp.samp$log.imp <- dmvnorm(x=imp.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)
  colnames(imp.samp$samples) <- colnames(post.samp$samples)
  # compute posterior log-likelihood of importance samples

  # set auxiliary parameters as trimmed_forcing for the relevant model
  if (gpd.model == 'gpd3') {aux <- NULL
  } else {aux <- covariates[,cc]; forc_max <- max(aux)}

  imp.samp$log.p <- apply(imp.samp$samples, 1, log_post_ppgpd,
                          parnames=colnames(post.samp$samples),
                          priors=priors[[cc]],
                          data_calib=data_calib,
                          model = gpd.model,
                          auxiliary=aux)

  print('done!')

  print('beginning bridge sampling recursion...')

  # set tolerance for halting of iteration
  TOL <- 1e-10

  # initialize storage for estimates
  ml <- mat.or.vec(nr=1,nc=1)

  # initialize with starting value
  # we can't quite start with the reciprocal importance sampling estimate from
  # Gelfand and Dey (1994) due to numerics (we get 0 values when we exponentiate
  # the difference of the importance log-densities and posterior log-likelihoods), so we just
  # average the ratios on a log scale.
  ml[1] <- -mean(post.samp$log.imp - post.samp$log.p)
  ml[2] <- bridge.samp.iter(ml[1], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

  # iterate until within tolerance.
  t <- 2
  while (abs(ml[t] - ml[t-1]) >= TOL) {
    ml[t+1] <- bridge.samp.iter(ml[t], post.samp[c('log.p', 'log.imp')], imp.samp[c('log.p', 'log.imp')])
    t <- t+1
  }

  print('done!')

  print('computing relative standard error of estimate')

  # compute the relative standard error of the bridge sampling estimator
  # we can treat the posterior samples as iid due to re-sampling from the posterior,
  # so we use the error formula from Fruhwirth-Schnatter (2004) with the spectral density
  # at frequency 0 set equal to 1.

  re.sq <- bridge.samp.rel.err(ml[length(ml)], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

  # save result of run
  # if save directory doesn't exist, create it
  #ifelse(!dir.exists(path.save), dir.create(path.save), FALSE)
  setwd(path.save)

  save(list=c('post.samp','imp.samp', 'ml', 're.sq', 'station', 'covariate', 'gpd.model'), file=filename.out)
  output[[ee]] <- 'done!'
 }
}
stopCluster(cl)
#data_many <- finalOutput
#names(data_many) <- names(data_set)


#===============================================================================
# end
#===============================================================================
