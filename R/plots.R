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
bma.weights     <- readRDS('../output/bma_weights_threshold99.rds')


# Load the covariate names and values ==========================================
source('get_timeseries_covariates.R')


# FIGURE ??? ===================================================================
# barplot of BMA weights all together

bw <- rev(sort(unlist(bma.weights.all)))
new_names <- NULL
for (i in 1:length(bw)) {
  cc_mod <- unlist(strsplit(names(bw)[i], split="[.]"))[3]
  cc_mod <- unlist(strsplit(cc_mod, split="[_]"))
  cc <- capitalize(cc_mod[1])
  mod <- cc_mod[2]
  if (cc=='Time' & mod=='gpd3') {cc='None'
  } else if (cc=='Nao') {cc <- 'NAO'} else if (cc=='Sealevel') {cc <- 'Sea level'} else if (cc=='Temp') {cc <- 'Temperature'}
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
#       time        temp    sealevel         nao
#0.001693089 0.286718417 0.469655883 0.118538700
## A check things work as expected:
#> as.numeric(sum(bw_totals) + bw['None, ST'])
#[1] 1
## And normalized to the non-stationary models only:
#> bw_totals/sum(bw_totals)
#       time        temp    sealevel         nao
#0.001931414 0.327077830 0.535766165 0.135224591



pdf(paste(plot.dir,'bma_weights_all.pdf',sep='/'),width=5,height=4,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(0.68,1.7,.05,.5), las=1)
barplot(bw, horiz=TRUE, names.arg='', xlab='', ylab='', space=1, xlim=c(0,0.21))
mtext(side=1, text='BMA weight', line=2.2)
axis(2, at=seq(1.5,2*length(bw),2), labels=new_names)
for (cc in 1:length(bw)) {text(bw[cc]+0.041, 1.5+2*(cc-1), paste(round(bw[cc],digits=3)), pos=2)}
dev.off()

#===============================================================================


# FIGURE ??? ===================================================================
# barplot of BMA weights for each individual covariate model

site <- 'Norfolk'
better_names <- c('ST','NS1','NS2','NS3')
covar_names <- c('Time', 'Temperature','Sea level', 'NAO'); names(covar_names) <- names_covariates

bw_cov <- vector('list', length(names_covariates)); names(bw_cov) <- names_covariates
for (cc in names_covariates) {
  bw_cov[[cc]] <- bma.weights[[site]][[cc]]
  names(bw_cov[[cc]]) <- better_names
}


pdf(paste(plot.dir,'bma_weights_covariates.pdf',sep='/'),width=5,height=5,colormodel='cmyk')
par(mfrow=c(2,2), mai=c(.6,.6,.3,.2), mgp=c(3,.5,0))
# Time
barplot(bw_cov$time, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0,1.18), xlim=c(0.5,8.5), main='Time covariate')
axis(2, at=c(0,.2,.4,.6,.8,1,1.2), labels=c('','0.2','','0.6','','1',''))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext(side=3, text=expression(a.), adj=-0.2, line=0.6)
for (cc in 1:length(bw_cov$time)) {text(1.5+2*(cc-1), bw_cov$time[cc]+0.035, paste(round(bw_cov$time[cc],digits=3)), pos=3)}
# Temperature
barplot(bw_cov$temp, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0, 0.39), xlim=c(0.5,8.5), main='Temperature covariate')
axis(2, at=seq(0,0.4,0.1), labels=c('0','0.1','0.2','0.3','0.4'))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext(side=3, text=expression(b.), adj=-0.2, line=0.6)
for (cc in 1:length(bw_cov$temp)) {text(1.5+2*(cc-1), bw_cov$temp[cc]+0.01, paste(round(bw_cov$temp[cc],digits=3)), pos=3)}
# Sea level
barplot(bw_cov$sealevel, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0, 0.38), xlim=c(0.5,8.5), main='Sea level covariate')
axis(2, at=seq(0,0.4,0.1), labels=c('0','0.1','0.2','0.3','0.4'))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext(side=3, text=expression(c.), adj=-0.2, line=0.6)
for (cc in 1:length(bw_cov$sealevel)) {text(1.5+2*(cc-1), bw_cov$sealevel[cc]+0.01, paste(round(bw_cov$sealevel[cc],digits=3)), pos=3)}
# NAO index
barplot(bw_cov$nao, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0, 0.62), xlim=c(0.5,8.5), main='NAO index covariate')
axis(2, at=seq(0,0.6,0.1), labels=c('0','','0.2','','0.4','','0.6'))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext(side=3, text=expression(d.), adj=-0.2, line=0.6)
for (cc in 1:length(bw_cov$nao)) {text(1.5+2*(cc-1), bw_cov$nao[cc]+0.015, paste(round(bw_cov$nao[cc],digits=3)), pos=3)}
dev.off()

#===============================================================================


#===============================================================================
# End
#===============================================================================
