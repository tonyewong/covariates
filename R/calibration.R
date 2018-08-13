#===============================================================================
# Calibration of PP-GPD model(s) for storm surge at Norfolk, Virgina, USA.
# This version only does the four ST/NS1/NS2/NS3 experiments using ALL of the
# data
#
# Questions? Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================

rm(list=ls())

# vvv IMPORTANT SETTINGS YOU SHOULD MODIFY, DEPENDING ON THE EXPERIMENT vvv

station <- 'norfolk'             # can be 'norfolk' only for this study
type.of.priors <- 'normalgamma'      # can be only 'normalgamma' for this study
pot.threshold <- 0.99            # GPD threshold (percentile, 0-1)
dt.decluster <- 3                # declustering time-scale (days)

niter_mcmc_prelim000 <- 5e3      # number of MCMC iterations (PRELIMINARY chains)
nnode_mcmc_prelim000 <- 1        # number of CPUs to use (PRELIMINARY chains)
niter_mcmc_prod000 <- 1e5        # number of MCMC iterations (PRODUCTION chains)
#nnode_mcmc_prod000 <- 10        # number of CPUs to use (PRODUCTION chains)
gamma_mcmc000 <- 0.66            # speed of adaptation (0.5=faster, 1=slowest)
                                 # 0.5 is too fast for the time covariate.

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/covariates/R')
  nnode_mcmc_prod000 <- 1          # number of CPUs to use (PRODUCTION chains)
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/covariates/R')
  nnode_mcmc_prod000 <- 10          # number of CPUs to use (PRODUCTION chains)
}

#
#===============================================================================
# Here, set the file names for the prior distribution RDS file from the 'fit
# priors' script and the calibration data files from the processing scripts.
# This shoudl work, as long as you have not repeated experiments.
#===============================================================================
#

output.dir <- '../output/'
dat.dir <- '../data/'

filename.mles <- Sys.glob(paste(dat.dir,'surge_MLEs_ppgpd_decl',dt.decluster,'-pot',100*pot.threshold,'_*','.rds',sep=''))
filename.priors <- Sys.glob(paste(dat.dir,'surge_priors_',type.of.priors,'_ppgpd_decl',dt.decluster,'-pot',100*pot.threshold,'_*','.rds',sep=''))

appen <- paste('ppgpd-experiments_',station,'_',type.of.priors,'_decl',dt.decluster,'-pot',pot.threshold*100,sep='')

if (station=='delfzijl') {
  ind.in.mles <- 29
} else if (station=='norfolk') {
  ind.in.mles <- 30
}

# ^^^ IMPORTANT SETTINGS YOU SHOULD MODIFY, DEPENDING ON THE EXPERIMENT ^^^



# Name the calibrated parameters output file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
filename.parameters <- paste(output.dir,'calibratedParameters_',appen,'_',today,'.nc',sep='')

# Name the saved progress RData workspace image file
filename.mcmc <- paste(output.dir,'mcmc_',appen,'_',today,'.RData', sep='')

#
#===============================================================================
# Nothing below here should need to be modified. And really, all you ought to
# need to do is adjust the name of the station to calibrate, or the set of prior
# distributions for a supplementary experiment.
#===============================================================================
#

# On with the show!

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
# parameters). yields: nao, time, sealevel, temperature all Nyear x 2 with year
# as first column and data as the second. And data_calib calibration data for
# PP/GPD trimmed to match the time periods of overlap for these auxiliary fields
#===============================================================================
#

print('getting forcing data and trimming data_calib...')

source('get_timeseries_covariates.R')

print('...done.')

#
#===============================================================================
# set up PP-GPD model parameters
#===============================================================================
#

print('setting up PP-GPD model parameters from DE optimization...')

source('parameter_setup.R')
source('likelihood_ppgpd.R')

priors <- readRDS(filename.priors)
mle.fits <- readRDS(filename.mles)

initial.values <- vector('list', nmodel); names(initial.values) <- types.of.gpd
for (model in types.of.gpd) {initial.values[[model]] <- mle.fits[[model]][ind.in.mles,]}

print('...done.')

#
#===============================================================================
# set up and run PRELIMINARY Markov chain Monte Carlo (MCMC) calibration
#===============================================================================
#

# first, do a set of single-chain preliminary calibrations to get estimates of
# the jump covariance matrix

nnode_mcmc <- nnode_mcmc_prelim000
niter_mcmc <- niter_mcmc_prelim000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_prelim <- vector('list', length(names_covariates)); names(amcmc_prelim) <- names_covariates
for (cc in names_covariates) {amcmc_prelim[[cc]] <- vector('list', nmodel); names(amcmc_prelim[[cc]]) <- types.of.gpd}

for (cc in names_covariates) {
  print(paste('Starting preliminary calibration for covariate ',cc,'...', sep=''))

  for (model in types.of.gpd) {
      print(paste('Starting preliminary calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

      if(model=='gpd3') {auxiliary <- NULL
      } else {auxiliary <- covariates[,cc]; forc_max <- max(auxiliary)}

      accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
      step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
      tbeg=proc.time()
      amcmc_prelim[[cc]][[model]] = MCMC(log_post_ppgpd, niter_mcmc, initial.values[[model]],
                                adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                                gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                                parnames=parnames_all[[model]], data_calib=data_calib,
                                priors=priors, auxiliary=auxiliary, model=model)
      tend=proc.time()

      print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }

}


# save progress
print(paste('saving preliminary results as .RData file (',filename.mcmc,') to read and use later...',sep=''))
save.image(file=filename.mcmc)
print('...done.')


#
#===============================================================================
# set up and run PRODUCTION MCMC calibration
#===============================================================================
#

# then use these initial estimates of step_mcmc to launch the parallel chains
# (from amcmc_prelim[[cc]][[model]]$cov.jump)

nnode_mcmc <- nnode_mcmc_prod000
niter_mcmc <- niter_mcmc_prod000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_out <- vector('list', length(names_covariates)); names(amcmc_out) <- names_covariates
for (cc in names_covariates) {amcmc_out[[cc]] <- vector('list', nmodel); names(amcmc_out[[cc]]) <- types.of.gpd}

for (cc in names_covariates) {
  print(paste('Starting production calibration for covariates ',cc,'...', sep=''))

  for (model in types.of.gpd) {
    print(paste('Starting production calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

    initial_parameters <- amcmc_prelim[[cc]][[model]]$samples[amcmc_prelim[[cc]][[model]]$n.sample,]
    if(model=='gpd3') {auxiliary <- NULL
    } else {auxiliary <- covariates[,cc]; forc_max <- max(auxiliary)}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- amcmc_prelim[[cc]][[model]]$cov.jump
    if(nnode_mcmc==1) {
      # do single chain
      tbeg=proc.time()
      amcmc_out[[cc]][[model]] = MCMC(log_post_ppgpd, niter_mcmc, initial_parameters,
                             adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                             gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                             parnames=parnames_all[[model]], data_calib=data_calib,
                             priors=priors, auxiliary=auxiliary, model=model)
      tend=proc.time()
    } else if(nnode_mcmc > 1) {
      # do parallel chains
      tbeg <- proc.time()
      amcmc_out[[cc]][[model]] <- MCMC.parallel(log_post_ppgpd, niter_mcmc, initial_parameters,
                             n.chain=nnode_mcmc, n.cpu=nnode_mcmc, packages='extRemes',
                             scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                             gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                             parnames=parnames_all[[model]], data_calib=data_calib,
                             priors=priors, auxiliary=auxiliary, model=model)
      tend <- proc.time()
    }
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
  print(paste('saving production results (so far) as .RData file (',filename.mcmc,') to read and use later...',sep=''))
  save.image(file=filename.mcmc)
  print('...done.')
}


#
#===============================================================================
# convergence diagnostics
#===============================================================================
#

# Gelman and Rubin diagnostics - determine and chop off for burn-in
niter.test <- seq(from=round(0.1*niter_mcmc), to=niter_mcmc, by=round(0.05*niter_mcmc))
gr.test <- vector('list', length(names_covariates)); names(gr.test) <- names_covariates
for (cc in names_covariates) {
  gr.test[[cc]] <- mat.or.vec(length(niter.test), nmodel)
  colnames(gr.test[[cc]]) <- types.of.gpd
}
gr.tmp <- rep(NA, length(niter.test))

for (cc in names_covariates) {
  for (model in types.of.gpd) {
    if(nnode_mcmc == 1) {
      # don't do GR stats, just cut off first half of chains
      print('only one chain; will lop off first half for burn-in instead of doing GR diagnostics')
    } else if(nnode_mcmc > 1) {
      # this case is FAR more fun
      # accumulate the names of the soon-to-be mcmc objects
      string.mcmc.list <- 'mcmc1'
      for (m in 2:nnode_mcmc) {
        string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
      }
      for (i in 1:length(niter.test)) {
        for (m in 1:nnode_mcmc) {
          # convert each of the chains into mcmc object
          eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_out[[cc]][[model]][[m]]$samples[1:niter.test[i],])', sep='')))
        }
        eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

        gr.test[[cc]][i,model] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
      }
    } else {print('error - nnode_mcmc < 1 makes no sense')}
  }
}

# Monitor posterior 5, 50 and 95% quantiles for drift
# Only checking for one of the chains
quant <- vector('list', length(names_covariates)); names(quant) <- names_covariates
for (cc in names_covariates) {quant[[cc]] <- vector('list', nmodel); names(quant[[cc]]) <- types.of.gpd}
names.monitor <- c('q05', 'q50', 'q95')
for (cc in names_covariates) {
  for (model in types.of.gpd) {
    quant[[cc]][[model]] <- vector('list', 3); names(quant[[cc]][[model]]) <- names.monitor
    for (q in names.monitor) {
      quant[[cc]][[model]][[q]] <- mat.or.vec(length(niter.test)-1, length(parnames_all[[model]]))
      for (i in 1:(length(niter.test)-1)) {
        if(nnode_mcmc==1) {
          quant[[cc]][[model]][[q]][i,] <- apply(X=amcmc_out[[cc]][[model]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
        } else {
          quant[[cc]][[model]][[q]][i,] <- apply(X=amcmc_out[[cc]][[model]][[1]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
        }
      }
    }
  }
}


# Heidelberger and Welch diagnostics?
hw.diag <- vector('list', length(names_covariates)); names(hw.diag) <- names_covariates
for (cc in names_covariates) {
  hw.diag[[cc]] <- vector('list', nmodel); names(hw.diag[[cc]]) <- types.of.gpd
  for (model in types.of.gpd) {
    hw.diag[[cc]][[model]] <- heidel.diag(as.mcmc(amcmc_out[[cc]][[model]][[1]]$samples), eps=0.1, pvalue=0.05)
  }
}

#
#===============================================================================
# Chop off burn-in
#===============================================================================
#

# Note: here, we are only using the Gelman and Rubin diagnostic. But this is
# only after looking at the quantile stability as iterations increase, as well
# as the Heidelberger and Welch diagnostics, which suggest the chains are okay.
# 'ifirst' is the first spot where the GR stat gets to and stays below gr.max
# for all of the models.
# save a separate ifirst for each experiment
ifirst <- rep(NA, length(names_covariates)); names(ifirst) <- names_covariates
for (cc in names_covariates) {
  if(nnode_mcmc==1) {
    ifirst[[cc]] <- round(0.5*niter_mcmc)
  } else {
    gr.max <- 1.1
    lgr <- rep(NA, length(niter.test))
    for (i in 1:length(niter.test)) {lgr[i] <- all(gr.test[[cc]][i,] < gr.max)}
    for (i in seq(from=length(niter.test), to=1, by=-1)) {
      if( all(lgr[i:length(lgr)]) ) {ifirst[[cc]] <- niter.test[i]}
    }
  }
}

chains_burned <- vector('list', length(names_covariates)); names(chains_burned) <- names_covariates
for (cc in names_covariates) {
  chains_burned[[cc]] <- vector('list', nmodel); names(chains_burned[[cc]]) <- types.of.gpd
  for (model in types.of.gpd) {
    if(nnode_mcmc > 1) {
      chains_burned[[cc]][[model]] <- vector('list', nnode_mcmc)
      for (m in 1:nnode_mcmc) {
        chains_burned[[cc]][[model]][[m]] <- amcmc_out[[cc]][[model]][[m]]$samples[(ifirst[[cc]]+1):niter_mcmc,]
      }
    } else {
      chains_burned[[cc]][[model]] <- amcmc_out[[cc]][[model]]$samples[(ifirst[[cc]]+1):niter_mcmc,]
    }
  }
}

#
#===============================================================================
# possible thinning?
#===============================================================================
#

# If no thinning, then this initialization will remain
chains_burned_thinned <- chains_burned

if(FALSE) {#==========================

acf_cutoff <- 0.05
lag_max <- 0.01*niter_mcmc # if we end up with fewer than 100 samples, what are we even doing?
niter_thin <- rep(0, nmodel); names(niter_thin) <- types.of.model
for (model in types.of.model) {
  for (p in 1:length(parnames_all[[model]])) {
    if(nnode_mcmc > 1) {acf_tmp <- acf(x=chains_burned[[model]][[1]][,p], plot=FALSE, lag.max=lag_max)}
    else {acf_tmp <- acf(x=chains_burned[[model]][,p], plot=FALSE, lag.max=lag_max)}
    niter_thin[[model]] <- max(niter_thin[[model]], acf_tmp$lag[which(acf_tmp$acf < acf_cutoff)[1]])
  }
  nthin <- max(niter_thin, na.rm=TRUE)
  if(nnode_mcmc > 1) {
    for (m in 1:nnode_mcmc) {
      chains_burned_thinned[[model]][[m]] <- chains_burned[[model]][[m]][seq(from=1, to=nrow(chains_burned[[model]][[m]]), by=nthin),]
    }
  } else {
    chains_burned_thinned[[model]] <- chains_burned[[model]][seq(from=1, to=nrow(chains_burned[[model]]), by=nthin),]
  }
}

}#====================================


# thin to a target number of samples?
if(TRUE) {#===========================

n.sample <- 10000

for (cc in names_covariates) {
  for (model in types.of.gpd) {
    if(nnode_mcmc == 1) {
      ind.sample <- sample(x=1:nrow(chains_burned[[cc]][[model]]), size=n.sample, replace=FALSE)
      chains_burned_thinned[[cc]][[model]] <- chains_burned[[cc]][[model]][ind.sample,]
    } else {
      n.sample.sub <- rep(NA, nnode_mcmc)
      # the case where desired sample size is divisible by the number of chains
      if(round(n.sample/nnode_mcmc) == n.sample/nnode_mcmc) {
        n.sample.sub[1:nnode_mcmc] <- n.sample/nnode_mcmc
      } else {
      # the case where it is not
        n.sample.sub[2:nnode_mcmc] <- round(n.sample/nnode_mcmc)
        n.sample.sub[1] <- n.sample - sum(n.sample.sub[2:nnode_mcmc])
      }
      for (m in 1:nnode_mcmc) {
        ind.sample <- sample(x=1:nrow(chains_burned[[cc]][[model]][[m]]), size=n.sample.sub[m], replace=FALSE)
        chains_burned_thinned[[cc]][[model]][[m]] <- chains_burned[[cc]][[model]][[m]][ind.sample,]
      }
    }
  }
}

}#====================================


# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).
parameters.posterior <- vector('list', length(names_covariates)); names(parameters.posterior) <- names_covariates
covjump.posterior <- vector('list', length(names_covariates)); names(covjump.posterior) <- names_covariates
for (cc in names_covariates) {
  parameters.posterior[[cc]] <- vector('list', nmodel); names(parameters.posterior[[cc]]) <- types.of.gpd
  covjump.posterior[[cc]]    <- vector('list', nmodel); names(covjump.posterior[[cc]])    <- types.of.gpd
  for (model in types.of.gpd) {
    if(nnode_mcmc==1) {
      parameters.posterior[[cc]][[model]] <- chains_burned_thinned[[cc]][[model]]
      covjump.posterior[[cc]][[model]] <- amcmc_out[[cc]][[model]]$cov.jump
    } else {
      parameters.posterior[[cc]][[model]] <- chains_burned_thinned[[cc]][[model]][[1]]
      covjump.posterior[[cc]][[model]] <- amcmc_out[[cc]][[model]][[1]]$cov.jump
      for (m in 2:nnode_mcmc) {
        parameters.posterior[[cc]][[model]] <- rbind(parameters.posterior[[cc]][[model]], chains_burned_thinned[[cc]][[model]][[m]])
      }
    }
  }
}

# save results in case you need to revisit later
print(paste('saving MCMC results as .RData file (',filename.mcmc,') to read and use later...',sep=''))
save.image(file=filename.mcmc)
save(list=c('amcmc_out','ifirst'), file=filename.mcmc) # just do this later
print('...done.')

#
#===============================================================================
# write output file
#===============================================================================
#


## TESTING =====================
## TESTING =====================

if(FALSE) {

parameters.posterior <- vector('list', length(names_covariates))
names(parameters.posterior) <- names_covariates
covjump.posterior <- vector('list', length(names_covariates))
names(covjump.posterior) <- names_covariates
for (cc in names_covariates) {
  parameters.posterior[[cc]] <- vector('list', nmodel)
  names(parameters.posterior[[cc]]) <- types.of.gpd
  covjump.posterior[[cc]] <- vector('list', nmodel)
  names(covjump.posterior[[cc]]) <- types.of.gpd
  for (model in types.of.gpd) {
    parameters.posterior[[cc]][[model]] <- amcmc_prelim[[cc]][[model]]$samples[(0.5*niter_mcmc_prelim000+1):niter_mcmc_prelim000,]
    covjump.posterior[[cc]][[model]] <- amcmc_prelim[[cc]][[model]]$cov.jump
  }
}

amcmc_out <- amcmc_prelim
ifirst <- vector('list', length(names_covariates))
names(ifirst) <- names_covariates
for (cc in names_covariates) {
  ifirst[[cc]] <- 0
}

}

## TESTING =====================
## TESTING =====================


## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape we expect, just transpose it
lmax=0
for (model in types.of.gpd) {for (i in 1:length(parnames_all[[model]])){lmax=max(lmax,nchar(parnames_all[[model]][i]))}}

dim.parameters <- vector('list', nmodel); names(dim.parameters) <- types.of.gpd
dim.parnames   <- vector('list', nmodel); names(dim.parnames)   <- types.of.gpd
var.parnames   <- vector('list', nmodel); names(var.parnames)   <- types.of.gpd
var.parameters <- vector('list', length(types.of.gpd)); names(var.parameters) <- types.of.gpd
var.covjump    <- vector('list', length(types.of.gpd)); names(var.covjump)    <- types.of.gpd
for (cc in names_covariates) {
  var.parameters[[cc]] <- vector('list', nmodel); names(var.parameters[[cc]]) <- types.of.gpd
  var.covjump[[cc]]    <- vector('list', nmodel); names(var.covjump[[cc]])    <- types.of.gpd
}
dim.ensemble   <- vector('list', nmodel); names(dim.ensemble) <- types.of.gpd
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.time <- ncdim_def('ntime', '', 1:nrow(covariates), unlim=FALSE)
var.time <- ncvar_def('time', '', dim.time, -999)
for (model in types.of.gpd) {
  dim.parameters[[model]] <- ncdim_def(paste('n.parameters.',model,sep=''), '', 1:length(parnames_all[[model]]), unlim=FALSE)
  dim.ensemble[[model]]   <- ncdim_def(paste('n.ensemble.',model,sep=''), 'ensemble member', 1:nrow(parameters.posterior[[cc]][[model]]), unlim=FALSE)
  var.parnames[[model]]   <- ncvar_def(paste('parnames.',model,sep=''), '', list(dim.name,dim.parameters[[model]]), prec='char')
  for (cc in names_covariates) {
    var.parameters[[cc]][[model]] <- ncvar_def(paste('parameters.',cc,'.',model,sep=''), '', list(dim.parameters[[model]],dim.ensemble[[model]]), -999)
    var.covjump[[cc]][[model]] <- ncvar_def(paste('covjump.',cc,'.',model,sep=''), '', list(dim.parameters[[model]],dim.parameters[[model]]), -999)
  }
}
# length(types.of.gpd) * nmodel * 2 accounts for sizes of var.parmaeters and var.covjump
output.to.file <- vector('list', (length(types.of.gpd)*nmodel*2) + length(var.parnames) + 1)
output.to.file[[1]] <- var.time
# this counter will keep track of how many items are on the output.to.file list so far
cnt <- 2
for (model in types.of.gpd) {
  output.to.file[[cnt]] <- var.parnames[[model]]
  cnt <- cnt + 1
  for (cc in names_covariates) {
    output.to.file[[cnt]] <- var.parameters[[cc]][[model]]
    output.to.file[[cnt+1]] <- var.covjump[[cc]][[model]]
    cnt <- cnt + 2
  }
}
#outnc <- nc_create(filename.parameters, list(var.parameters, var.parnames, var.covjump))
outnc <- nc_create(filename.parameters, output.to.file)
ncvar_put(outnc, var.time, time[,1])
for (model in types.of.gpd) {
  ncvar_put(outnc, var.parnames[[model]], parnames_all[[model]])
  for (cc in names_covariates) {
    ncvar_put(outnc, var.parameters[[cc]][[model]], t(parameters.posterior[[cc]][[model]]))
    ncvar_put(outnc, var.covjump[[cc]][[model]], covjump.posterior[[cc]][[model]])
  }
}
nc_close(outnc)

# save results in case you need to revisit later
print(paste('saving MCMC results as .RData file (',filename.mcmc,') to read and use later...',sep=''))
##save.image(file=filename.mcmc)
save(list=c('amcmc_out','ifirst','data_calib','priors'), file=filename.mcmc)
print('...done.')

#
#===============================================================================
# End
#===============================================================================
#
