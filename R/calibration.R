#===============================================================================
# calibration.R
#
# Calibration of PP-GPD model for storm surge at Delfzijl, The Netherlands;
# and Norfolk, Virgina, USA.
# This version does the varying record length experiments and includes the
# experiment using all of the data for each site, eliminating the need for a
# separate script.
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

rm(list=ls())

# vvv IMPORTANT SETTINGS YOU SHOULD MODIFY, DEPENDING ON THE EXPERIMENT vvv

station <- 'delfzijl'             # can be 'delfzijl' or 'norfolk'
type.of.priors <- 'normalgamma'      # can be either 'uniform' or 'normalgamma'
pot.threshold <- 0.99            # GPD threshold (percentile, 0-1)
dt.decluster <- 3                # declustering time-scale (days)
filename.datacalib <- "../data/tidegauge_processed_norfolk-delfzijl_decl3-pot99-annual_07Jun2018.rds"

niter_mcmc_prelim000 <- 5e3      # number of MCMC iterations (PRELIMINARY chains)
nnode_mcmc_prelim000 <- 1        # number of CPUs to use (PRELIMINARY chains)
niter_mcmc_prod000 <- 5e3        # number of MCMC iterations (PRODUCTION chains)
#nnode_mcmc_prod000 <- 10          # number of CPUs to use (PRODUCTION chains)
gamma_mcmc000 <- 0.5             # speed of adaptation (0.5=faster, 1=slowest)

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/datalengths/R')
  nnode_mcmc_prod000 <- 1          # number of CPUs to use (PRODUCTION chains)
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/datalengths/R')
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

appen <- paste('datalengths_',station,'_',type.of.priors,'_decl',dt.decluster,'-pot',pot.threshold*100,sep='')

if (station=='delfzijl') {
  ind.in.mles <- 29
} else if (station=='norfolk') {
  ind.in.mles <- 30
}

print('')
print('-------------------------------------------------------------')
print(paste('Declustering time-scale:',dt.decluster,sep=' '))
print(paste('POT threshold:',pot.threshold,sep=' '))
print(paste('Prior MLEs filename:',filename.mles,sep=' '))
print(paste('Prior fits filename:',filename.priors,sep=' '))
print(paste('Calibration data filename:',filename.datacalib,sep=' '))
print('-------------------------------------------------------------')
print('')

# ^^^ IMPORTANT SETTINGS YOU SHOULD MODIFY, DEPENDING ON THE EXPERIMENT ^^^

# Name the calibrated parameters output file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
filename.parameters <- paste(output.dir,'calibratedParameters_',appen,'_',today,'.nc',sep='')

# Name the saved progress RData workspace image file
filename.prelim <- paste(output.dir,'mcmc-prelim_',appen,'_',today,'.RData', sep='')
filename.prod <- paste(output.dir,'mcmc-prod_',appen,'_',today,'.RData', sep='')
filename.results <- paste(output.dir,'mcmc-results_',appen,'_',today,'.RData', sep='')

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
# read processed data object
#===============================================================================
#

print('reading processed tide gauge data...')

data_calib <- readRDS(filename.datacalib)[[station]]
data.experiments <- names(data_calib)
n.experiments <- length(data.experiments)

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
initial.values <- mle.fits$gpd3[ind.in.mles,]

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

amcmc_prelim <- vector('list', n.experiments); names(amcmc_prelim) <- data.experiments
auxiliary <- NULL # left as an option if non-stationary models are desired
accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames)
step_mcmc <- as.numeric(0.05*apply(X=mle.fits$gpd3, MARGIN=2, FUN=sd))

for (data.exp in data.experiments) {
  print(paste('Starting preliminary calibration for experiment ',data.exp,'/',n.experiments,' -- ', nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))
  tbeg=proc.time()
  amcmc_prelim[[data.exp]] = MCMC(log_post_ppgpd, niter_mcmc, initial.values,
                              adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                              gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                              parnames=parnames, data_calib=data_calib[[data.exp]],
                              priors=priors, auxiliary=auxiliary, model='gpd3')
  tend=proc.time()
  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
}

# save progress
print(paste('saving preliminary results as .RData file (',filename.prelim,') to read and use later...',sep=''))
save(list=c('amcmc_prelim'), file=filename.prelim)
print('...done.')

#
#===============================================================================
# set up and run PRODUCTION MCMC calibration
#===============================================================================
#

# then use these initial estimates of step_mcmc to launch the parallel chains
# (from amcmc_prelim[[data.exp]]$cov.jump)

nnode_mcmc <- nnode_mcmc_prod000
niter_mcmc <- niter_mcmc_prod000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames)

amcmc_out <- vector('list', n.experiments); names(amcmc_out) <- data.experiments

for (data.exp in data.experiments) {
  print(paste('Starting production calibration for experiment ',data.exp,'/',n.experiments,' -- ', nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))
  initial_parameters <- amcmc_prelim[[data.exp]]$samples[amcmc_prelim[[data.exp]]$n.sample,]
  step_mcmc <- amcmc_prelim[[data.exp]]$cov.jump
  tbeg <- proc.time()
  if(nnode_mcmc==1) {
    # do single chain

    amcmc_out[[data.exp]] <- MCMC(log_post_ppgpd, niter_mcmc, initial_parameters,
                           adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                           gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                           parnames=parnames, data_calib=data_calib[[data.exp]],
                           priors=priors, auxiliary=auxiliary, model='gpd3')
  } else if(nnode_mcmc > 1) {
    # do parallel chains
    amcmc_out[[data.exp]] <- MCMC.parallel(log_post_ppgpd, niter_mcmc, initial_parameters,
                           n.chain=nnode_mcmc, n.cpu=nnode_mcmc, packages='extRemes',
                           scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                           gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                           parnames=parnames, data_calib=data_calib[[data.exp]],
                           priors=priors, auxiliary=auxiliary, model='gpd3')

  }
  tend <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  print(paste('saving production results (so far) as .RData file (',filename.prod,') to read and use later...',sep=''))
  save(list=c('amcmc_out'), file=filename.prod)
  print('...done.')
}

#
#===============================================================================
# convergence diagnostics
#===============================================================================
#

# Gelman and Rubin diagnostics - determine and chop off for burn-in
niter.test <- seq(from=round(0.1*niter_mcmc), to=niter_mcmc, by=round(0.05*niter_mcmc))
gr.test <- vector('list', n.experiments); names(gr.test) <- data.experiments
for (data.exp in data.experiments) {
  gr.test[[data.exp]] <- rep(0, length(niter.test))
}
gr.tmp <- rep(NA, length(niter.test))

for (data.exp in data.experiments) {
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
        eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_out[[data.exp]][[m]]$samples[1:niter.test[i],])', sep='')))
      }
      eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

      gr.test[[data.exp]][i] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
    }
  } else {print('error - nnode_mcmc < 1 makes no sense')}
}

# Monitor posterior 5, 50 and 95% quantiles for drift
# Only checking for one of the chains
quant <- vector('list', length(data.experiments)); names(quant) <- data.experiments
names.monitor <- c('q05', 'q50', 'q95')
for (data.exp in data.experiments) {
  quant[[data.exp]] <- vector('list', 3); names(quant[[data.exp]]) <- names.monitor
  for (q in names.monitor) {
    quant[[data.exp]][[q]] <- mat.or.vec(length(niter.test)-1, length(parnames))
    for (i in 1:(length(niter.test)-1)) {
      if(nnode_mcmc==1) {
        quant[[data.exp]][[q]][i,] <- apply(X=amcmc_out[[data.exp]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
      } else {
        quant[[data.exp]][[q]][i,] <- apply(X=amcmc_out[[data.exp]][[1]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
      }
    }
  }
}

if(FALSE) {
# examples monitoring of stability of quantiles:
data.exp <- 'y14'
par(mfrow=c(3,1))
for (p in 1:length(parnames)) {
  ran <- max(quant[[data.exp]]$q95[,p])-min(quant[[data.exp]]$q05[,p])
  lb <- min(quant[[data.exp]]$q05[,p]) - 0.05*ran; ub <- max(quant[[data.exp]]$q95[,p]) + 0.05*ran
  plot(niter.test[1:(length(niter.test)-1)], quant[[data.exp]]$q50[,p], type='l',
    ylim=c(lb,ub), ylab=parnames[p], xlab='From HERE to end of chain')
  lines(niter.test[1:(length(niter.test)-1)], quant[[data.exp]]$q05[,p], lty=2); lines(niter.test[1:(length(niter.test)-1)], quant[[data.exp]]$q95[,p], lty=2);
}
# the thing to note is that these stabilize as you include more members (i.e.,
# as you move from right to left)
}


# Heidelberger and Welch diagnostics?
hw.diag <- vector('list', n.experiments); names(hw.diag) <- data.experiments
for (data.exp in data.experiments) {
  hw.diag[[data.exp]] <- heidel.diag(as.mcmc(amcmc_out[[data.exp]][[1]]$samples), eps=0.1, pvalue=0.05)
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
gr.max <- 1.1
ifirst <- rep(NA, n.experiments); names(ifirst) <- data.experiments
for (data.exp in data.experiments) {
  if(nnode_mcmc==1) {
    ifirst[[data.exp]] <- round(0.5*niter_mcmc)
  } else {
    for (i in seq(from=length(niter.test), to=1, by=-1)) {
      if( all(gr.test[[data.exp]][i:length(lgr)] < gr.max) ) {ifirst[[data.exp]] <- niter.test[i]}
    }
  }
}

chains_burned <- vector('list', n.experiments); names(chains_burned) <- data.experiments
for (data.exp in data.experiments) {
  if(nnode_mcmc > 1) {
    chains_burned[[data.exp]] <- vector('list', nnode_mcmc)
    for (m in 1:nnode_mcmc) {
      chains_burned[[data.exp]][[m]] <- amcmc_out[[data.exp]][[m]]$samples[(ifirst[[data.exp]]+1):niter_mcmc,]
    }
  } else {
    chains_burned[[data.exp]] <- amcmc_out[[data.exp]]$samples[(ifirst[[data.exp]]+1):niter_mcmc,]
  }
}

# thin to a target number of samples?

# If no thinning, then this initialization will remain
chains_burned_thinned <- chains_burned

if(TRUE) {#===========================

n.sample <- 10000

for (data.exp in data.experiments) {
  if(nnode_mcmc == 1) {
    ind.sample <- sample(x=1:nrow(chains_burned[[data.exp]]), size=n.sample, replace=FALSE)
    chains_burned_thinned[[data.exp]] <- chains_burned[[data.exp]][ind.sample,]
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
      ind.sample <- sample(x=1:nrow(chains_burned[[data.exp]][[m]]), size=n.sample.sub[m], replace=FALSE)
      chains_burned_thinned[[data.exp]][[m]] <- chains_burned[[data.exp]][[m]][ind.sample,]
    }
  }
}

}#====================================


# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).
parameters.posterior <- vector('list', n.experiments); names(parameters.posterior) <- data.experiments
covjump.posterior <- vector('list', n.experiments); names(covjump.posterior) <- data.experiments
for (data.exp in data.experiments) {
  if(nnode_mcmc==1) {
    parameters.posterior[[data.exp]] <- chains_burned_thinned[[data.exp]]
    covjump.posterior[[data.exp]] <- amcmc_out[[data.exp]]$cov.jump
  } else {
    parameters.posterior[[data.exp]] <- chains_burned_thinned[[data.exp]][[1]]
    covjump.posterior[[data.exp]] <- amcmc_out[[data.exp]][[1]]$cov.jump
    for (m in 2:nnode_mcmc) {
      parameters.posterior[[data.exp]] <- rbind(parameters.posterior[[data.exp]], chains_burned_thinned[[data.exp]][[m]])
    }
  }
}

# save results in case you need to revisit later
print(paste('saving MCMC results as .RData file (',filename.prod,') to read and use later...',sep=''))
save(list=c('amcmc_out','ifirst'), file=filename.prod)
save(list=c('parameters.posterior'), file=filename.results)
print('...done.')

#
#===============================================================================
# write output file
#===============================================================================
#

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape we expect, just transpose it
lmax=0
for (i in 1:length(parnames)) {lmax=max(lmax,nchar(parnames[i]))}

dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.parameters <- ncdim_def('n.parameters', '', 1:length(parnames), unlim=FALSE)
dim.ensemble   <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior[[data.exp]]), unlim=FALSE)

var.parnames   <- ncvar_def('parnames', '', list(dim.name, dim.parameters), prec='char')
var.parameters <- vector('list', n.experiments); names(var.parameters) <- data.experiments
var.covjump    <- vector('list', n.experiments); names(var.covjump)    <- data.experiments
for (data.exp in data.experiments) {
  var.parameters[[data.exp]] <- ncvar_def(paste('parameters.',data.exp,sep=''), '', list(dim.parameters, dim.ensemble), -999)
  var.covjump[[data.exp]] <- ncvar_def(paste('covjump.',data.exp,sep=''), '', list(dim.parameters, dim.parameters), -999)
}

# n.experiments * 2 accounts for sizes of var.parmaeters and var.covjump
output.to.file <- vector('list', (n.experiments*2) + 1)
# this counter will keep track of how many items are on the output.to.file list so far
cnt <- 1
output.to.file[[cnt]] <- var.parnames
cnt <- cnt + 1
for (data.exp in data.experiments) {
  output.to.file[[cnt]] <- var.parameters[[data.exp]]
  output.to.file[[cnt+1]] <- var.covjump[[data.exp]]
  cnt <- cnt + 2
}
#outnc <- nc_create(filename.parameters, list(var.parameters, var.parnames, var.covjump))
outnc <- nc_create(filename.parameters, output.to.file)
ncvar_put(outnc, var.parnames, parnames)
for (data.exp in data.experiments) {
  ncvar_put(outnc, var.parameters[[data.exp]], t(parameters.posterior[[data.exp]]))
  ncvar_put(outnc, var.covjump[[data.exp]], covjump.posterior[[data.exp]])
}
nc_close(outnc)


#
#===============================================================================
# End
#===============================================================================
#
