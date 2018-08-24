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
filename.parameters <- '../output/calibratedParameters_ppgpd-experiments_norfolk_normalgamma_decl3-pot99_21Aug2018.nc'

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
rl$bma <- array(NA, dim=c(n.ensemble, n.rl))
colnames(rl$bma) <- names_rl

for (year in proj_years) {
#### WARNING
#### WARNING -- WOULD NEED TO ADD ANOTHER LAYER TO rl IF USING MORE proj_years
#### WARNING
  ind_forc <- which(covariates_proj[,'year']==year)
  forcing <- covariates_proj[ind_forc,]
  for (cc in names_covariates) {
    for (model in types.of.gpd) {
      print(paste('Covariate',cc,'and model',model,'...'))
      pb <- txtProgressBar(min=0,max=n.ensemble, initial=0,style=3)
      for (sow in 1:n.ensemble) {
        for (yy in 1:length(rl_years)) {
          if (length(parnames_all[[model]])==3) {
            lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda', parnames_all[[model]])]
            sigma <- gpd.parameters[[cc]][[model]][sow,match('sigma', parnames_all[[model]])]
            xi <- gpd.parameters[[cc]][[model]][sow,match('xi', parnames_all[[model]])]
          } else if(length(parnames_all[[model]])==4) {
            lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda0', parnames_all[[model]])] +
                      gpd.parameters[[cc]][[model]][sow,match('lambda1', parnames_all[[model]])]*forcing[[cc]]
            sigma <- gpd.parameters[[cc]][[model]][sow,match('sigma', parnames_all[[model]])]
            xi <- gpd.parameters[[cc]][[model]][sow,match('xi', parnames_all[[model]])]
          } else if(length(parnames_all[[model]])==5) {
            lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda0', parnames_all[[model]])] +
                      gpd.parameters[[cc]][[model]][sow,match('lambda1', parnames_all[[model]])]*forcing[[cc]]
            sigma <- exp(gpd.parameters[[cc]][[model]][sow,match('sigma0', parnames_all[[model]])] +
                         gpd.parameters[[cc]][[model]][sow,match('sigma1', parnames_all[[model]])]*forcing[[cc]])
            xi <- gpd.parameters[[cc]][[model]][sow,match('xi', parnames_all[[model]])]
          } else if(length(parnames_all[[model]])==6) {
            lambda <- gpd.parameters[[cc]][[model]][sow,match('lambda0', parnames_all[[model]])] +
                      gpd.parameters[[cc]][[model]][sow,match('lambda1', parnames_all[[model]])]*forcing[[cc]]
            sigma <- exp(gpd.parameters[[cc]][[model]][sow,match('sigma0', parnames_all[[model]])] +
                         gpd.parameters[[cc]][[model]][sow,match('sigma1', parnames_all[[model]])]*forcing[[cc]])
            xi <- gpd.parameters[[cc]][[model]][sow,match('xi0', parnames_all[[model]])] +
                  gpd.parameters[[cc]][[model]][sow,match('xi1', parnames_all[[model]])]*forcing[[cc]]
          }
          rl[[cc]][[model]][sow,yy] <- rlevd(rl_years[yy], scale=sigma,
                                             shape=xi, threshold=data_calib$threshold,
                                             type='GP', npy=365.25, rate=lambda)
          setTxtProgressBar(pb, sow)
        }
      }
      close(pb)
    }
  }
}

# Calculate the BMA ensembles for each of the covariates individually
for (year in proj_years) {
  for (cc in names_covariates) {
    for (yy in 1:length(rl_years)) {
      rl[[cc]]$bma[,yy] <- bma.weights$Norfolk[[cc]][['gpd3']] * rl[[cc]]$gpd3[,yy] +
                           bma.weights$Norfolk[[cc]][['gpd4']] * rl[[cc]]$gpd4[,yy] +
                           bma.weights$Norfolk[[cc]][['gpd5']] * rl[[cc]]$gpd5[,yy] +
                           bma.weights$Norfolk[[cc]][['gpd6']] * rl[[cc]]$gpd6[,yy]
    }
  }
}

# Calculate all-together BMA ensemble
for (year in proj_years) {
  for (yy in 1:length(rl_years)) {
    rl$bma[,yy] <- bma.weights.all$Norfolk$time_gpd3 * rl$time$gpd3[,yy] +
                   bma.weights.all$Norfolk$time_gpd4 * rl$time$gpd4[,yy] +
                   bma.weights.all$Norfolk$time_gpd5 * rl$time$gpd5[,yy] +
                   bma.weights.all$Norfolk$time_gpd6 * rl$time$gpd6[,yy] +
                   bma.weights.all$Norfolk$temp_gpd4 * rl$temp$gpd4[,yy] +
                   bma.weights.all$Norfolk$temp_gpd5 * rl$temp$gpd5[,yy] +
                   bma.weights.all$Norfolk$temp_gpd6 * rl$temp$gpd6[,yy] +
                   bma.weights.all$Norfolk$sealevel_gpd4 * rl$sealevel$gpd4[,yy] +
                   bma.weights.all$Norfolk$sealevel_gpd5 * rl$sealevel$gpd5[,yy] +
                   bma.weights.all$Norfolk$sealevel_gpd6 * rl$sealevel$gpd6[,yy] +
                   bma.weights.all$Norfolk$nao_gpd4 * rl$nao$gpd4[,yy] +
                   bma.weights.all$Norfolk$nao_gpd5 * rl$nao$gpd5[,yy] +
                   bma.weights.all$Norfolk$nao_gpd6 * rl$nao$gpd6[,yy]
  }
}

# get rid of NAs
for (cc in names_covariates) {
  for (model in types.of.gpd) {
    ind_na <- NULL
    for (yy in 1:length(rl_years)) {
      ind_na <- c(ind_na, which(is.na(rl[[cc]][[model]][,yy])))
    }
    if(length(ind_na) > 0) {rl[[cc]][[model]] <- rl[[cc]][[model]][-ind_na,]}
  }
  ind_na <- NULL
  for (yy in 1:length(rl_years)) {
    ind_na <- c(ind_na, which(is.na(rl[[cc]]$bma[,yy])))
  }
  if(length(ind_na) > 0) {rl[[cc]]$bma <- rl[[cc]]$bma[-ind_na,]}
}

ind_na <- NULL
for (yy in 1:length(rl_years)) {
  ind_na <- c(ind_na, which(is.na(rl$bma[,yy])))
  if(length(ind_na) > 0) {rl$bma <- rl$bma[-ind_na,]}
}


# Calculate pdfs and quantiles =================================================
# for plotting of return levels

quantiles.to.grab <- c(.025, .25, .5, .75, .975)

pdf.rl <- q.rl <- vector('list', length(names_covariates)+1)
names(pdf.rl) <- names(q.rl) <- c(names_covariates, 'bma')
for (cc in names_covariates) {
  pdf.rl[[cc]] <- q.rl[[cc]] <- vector('list', nmodel+1)
  names(pdf.rl[[cc]]) <- names(q.rl[[cc]]) <- c(types.of.gpd, 'bma')
  for (model in names(pdf.rl[[cc]])) {
    pdf.rl[[cc]][[model]] <- q.rl[[cc]][[model]] <- vector('list', n.rl)
    names(pdf.rl[[cc]][[model]]) <- names(q.rl[[cc]][[model]]) <- names_rl
    for (yy in names_rl) {
      # converting to m from mm here
      pdf.rl[[cc]][[model]][[yy]] <- density(rl[[cc]][[model]][,yy]/1000, from=0, to=10)
      q.rl[[cc]][[model]][[yy]] <- quantile(rl[[cc]][[model]][,yy]/1000, quantiles.to.grab)
    }
  }
}
pdf.rl$bma <- q.rl$bma <- vector('list', n.rl)
names(pdf.rl$bma) <- names(q.rl$bma) <- names_rl
for (yy in names_rl) {
  pdf.rl$bma[[yy]] <- density(rl$bma[,yy]/1000, from=0, to=10)
  q.rl$bma[[yy]] <- quantile(rl$bma[,yy]/1000, quantiles.to.grab)
}


# Save image to use for plotting ===============================================
save.image('../output/analysis.RData')


#===============================================================================
# End
#===============================================================================
