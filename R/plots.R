#===============================================================================
# plots.R
#
# Questions?  Tony Wong (<anthony.e.wong@colorado.edu>)
#===============================================================================

rm(list=ls())
setwd('Users/tony/codes/covariates/R')
plot.dir <- '/Users/tony/codes/covariates/figures'


# Read in BMA results ==========================================================
bma.weights.all <- readRDS('../output/bma_weights_all_threshold99.rds')
bma.weights     <- readRDS('../output/bma_weights_all_threshold99.rds')


# Load the covariate names and values ==========================================
source('get_timeseries_covariates.R')


bw <- rev(sort(unlist(bma.weights.all)))
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
  } else if (grepl('Sea', cc)) {bw_totals['sealevel'] <- bw_totals['sealevel'] + bw[cc]
  } else if (grepl('Temp', cc)) {bw_totals['temp'] <- bw_totals['temp'] + bw[cc]
  } else if (grepl('NAO', cc)) {bw_totals['nao'] <- bw_totals['nao'] + bw[cc]
  }
}
#> bw_totals
#      time       temp   sealevel        nao
#0.09128698 0.29934667 0.43283316 0.17653319


pdf(paste(plot.dir,'bma_weights_all.pdf',sep='/'),width=5,height=4,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(0.65,1.7,.1,.5), las=1)
barplot(bw, horiz=TRUE, names.arg='', xlab='', ylab='', space=1, xlim=c(0,0.18))
mtext(side=1, text='BMA weight', line=2.2)
axis(2, at=seq(1.5,2*length(bw),2), labels=new_names)
#for (cc in 1:length(bw)) {text(bw[cc]+0.05, 1.5+2*(cc-1), paste(signif(bw[cc],2)), pos=2)}
for (cc in 1:length(bw)) {text(bw[cc]+0.04, 1.5+2*(cc-1), paste(round(bw[cc],digits=3)), pos=2)}
dev.off()



#==========
# individual covariates

site <- 'Norfolk'
better_names <- c('ST','NS1','NS2','NS3')
covar_names <- c('Time', 'Temperature','Sea level', 'NAO'); names(covar_names) <- names_covariates

bw_cov <- vector('list', length(names_covariates)); names(bw_cov) <- names_covariates
for (cc in names_covariates) {
  bw_cov[[cc]] <- bma.weights[[site]][[cc]]
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
