#===============================================================================
# analysis.R
#
# 1. Make projections
#
# Questions?  Tony Wong (<anthony.e.wong@colorado.edu>)
#===============================================================================

rm(list=ls())
setwd('Users/tony/codes/covariates/R')
plot.dir <- '/Users/tony/codes/covariates/figures'

rl_years <- c(20, 100) # in years return period, which return levels do you want?
proj_years <- c(2065)  # what years do we want return level projections for?

# Read in BMA results ==========================================================
bma.weights.all <- readRDS('../output/bma_weights_all_threshold99.rds')
bma.weights     <- readRDS('../output/bma_weights_threshold99.rds')


# Load the covariate names and values ==========================================
source('get_timeseries_covariates.R')


# Load parameter information ===================================================
source('parameter_setup.R')


# Load the calibrated parameters ===============================================
filename.parameters <- '../output/calibratedParameters_ppgpd-experiments_norfolk_normalgamma_decl3-pot99_12Aug2018.nc'

gpd.parameters <- vector('list', length(names_covariates))
names(gpd.parameters) <- names_covariates

for (cc in names_covariates) {
  gpd.parameters[[cc]] <- vector('list', nmodel)
  names(gpd.parameters[[cc]]) <- types.of.gpd
}

ncdata <- nc_open(filename.parameters)

for (cc in names_covariates) {
  for (model in types.of.gpd) {
    gpd.parameters[[cc]][[model]] <- t(ncvar_get(ncdata, paste('parameters',cc,model,sep='.')))
  }
}
nc_close(ncdata)
n.ensemble <- nrow(gpd.parameters[[cc]][[model]])


# Calculate return levels ======================================================
# First two things are set at the beginning so it is easier to modify and see
# what is being used
##rl_years <- c(20, 100) # in years return period, which return levels do you want?
##proj_years <- c(2065)
names_rl <- NULL
for (ii in 1:length(rl_years)) {
  names_rl <- c(names_rl, paste('y',rl_years[ii], sep=''))
}
n.rl <- length(rl_years)
names(rl_years) <- names_rl

rl <- vector('list', length(names_covariates)+1)
names(rl) <- c(names_covariates, 'bma')
for (cc in names_covariates) {
  rl[[cc]] <- vector('list', nmodel+1)
  names(rl[[cc]]) <- c(types.of.gpd, 'bma')
  for (model in names(rl[[cc]])) {
    rl[[cc]][[model]] <- array(NA, dim=c(n.ensemble, n.rl))
    colnames(rl[[cc]][[model]]) <- names_rl
  }
}

for (year in proj_years) {
  ind_forc <- which(covariates_proj[,'year']==year)
  forcing <- covariates_proj[ind_forc,]
  for (cc in names_covariates) {
    for (model in types.of.gpd) {
      print(paste('Covariate',cc,'and model',model,'...'))

      forcing[[cc]]
    }
  }

}


pb <- txtProgressBar(min=0,max=n.ensemble, initial=0,style=3)
for (sow in 1:n.ensemble) {
  if (length(parnames_all[[model]])==3) {
    lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda', parnames_all[[model]])]
    sigma <- gpd.parameters[[cc]][[model]][sow,match('sigma', parnames_all[[model]])]
    xi <- gpd.parameters[[cc]][[model]][sow,match('xi', parnames_all[[model]])]
  } else if(length(parnames_all[[model]])==4) {
    lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda0', parnames_all[[model]])] +
              gpd.parameters[[cc]][[model]][sow,match('lambda1', parnames_all[[model]])]*nao.years[[year]]
    sigma <- gpd.parameters[[cc]][[model]][sow,match('sigma', parnames_all[[model]])]
    xi <- gpd.parameters[[cc]][[model]][sow,match('xi', parnames_all[[model]])]
  } else if(length(parnames_all[[model]])==5) {
    lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda0', parnames_all[[model]])] +
              gpd.parameters[[cc]][[model]][sow,match('lambda1', parnames_all[[model]])]*nao.years[[year]]
    sigma <- exp(gpd.parameters[[cc]][[model]][sow,match('sigma0', parnames_all[[model]])] +
                 gpd.parameters[[cc]][[model]][sow,match('sigma1', parnames_all[[model]])]*nao.years[[year]])
    xi <- gpd.parameters[[cc]][[model]][sow,match('xi', parnames_all[[model]])]
  } else if(length(parnames_all[[model]])==6) {
    lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda0', parnames_all[[model]])] +
              gpd.parameters[[cc]][[model]][sow,match('lambda1', parnames_all[[model]])]*nao.years[[year]]
    sigma <- exp(gpd.parameters[[cc]][[model]][sow,match('sigma0', parnames_all[[model]])] +
                 gpd.parameters[[cc]][[model]][sow,match('sigma1', parnames_all[[model]])]*nao.years[[year]])
    xi <- gpd.parameters[[cc]][[model]][sow,match('xi0', parnames_all[[model]])] +
          gpd.parameters[[cc]][[model]][sow,match('xi1', parnames_all[[model]])]*nao.years[[year]]
  }
  rl100[[site]][[data.len]][[model]][[year]][sow] <- rlevd(100, scale=sigma, shape=xi,
                                                           threshold=data_calib$threshold,
                                                           type='GP',
                                                           npy=365.25,
                                                           rate=lambda)
  setTxtProgressBar(pb, sow)
}
close(pb)


#===============================================================================


#===============================================================================
# End
#===============================================================================
