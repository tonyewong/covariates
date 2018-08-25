#===============================================================================
# plots.R
#
# Questions?  Tony Wong (<anthony.e.wong@colorado.edu>)
#===============================================================================

rm(list=ls())
setwd('Users/tony/codes/covariates/R')
plot.dir <- '/Users/tony/codes/covariates/figures'


# Read in previous analysis work ===============================================
load('../output/analysis.RData')


# Read in BMA results ==========================================================
bma.weights.all <- readRDS('../output/bma_weights_all_threshold99.rds')
bma.weights     <- readRDS('../output/bma_weights_threshold99.rds')


# Load the covariate names and values ==========================================
source('get_timeseries_covariates.R')



#===============================================================================
# FIGURES
#===============================================================================



#===============================================================================
# FIGURE ???
# barplot of BMA weights all together
#===============================================================================

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
#0.2099934 0.1867303 0.1873299 0.1880989
## A check things work as expected:
#> as.numeric(sum(bw_totals) + bw['None, ST'])
#[1] 1
## And normalized to the non-stationary models only:
#> bw_totals/sum(bw_totals)
#       time        temp    sealevel         nao
#0.2719585 0.2418309 0.2426074 0.2436033


pdf(paste(plot.dir,'bma_weights_all.pdf',sep='/'),width=5,height=4,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(0.68,1.7,.05,.5), las=1)
barplot(bw, horiz=TRUE, names.arg='', xlab='', ylab='', space=1, xlim=c(0,0.52))
mtext(side=1, text='BMA weight', line=2.2)
axis(2, at=seq(1.5,2*length(bw),2), labels=new_names)
for (cc in 1:length(bw)) {text(bw[cc]+0.12, 1.5+2*(cc-1), paste(round(bw[cc],digits=3)), pos=2)}
dev.off()

#===============================================================================



#===============================================================================
# FIGURE ???
# barplot of BMA weights for each individual covariate model
#===============================================================================

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
barplot(bw_cov$time, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0,1), xlim=c(0.5,8.5))
axis(2, at=c(0,.2,.4,.6,.8,1), labels=c('','0.2','','0.6','','1'))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext('Time covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' a.')), line=.6, cex=1, adj=-0.25)
for (cc in 1:length(bw_cov$time)) {text(1.5+2*(cc-1), bw_cov$time[cc], paste(round(bw_cov$time[cc],digits=3)), pos=3)}
# Temperature
barplot(bw_cov$temp, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0,1), xlim=c(0.5,8.5))
axis(2, at=c(0,.2,.4,.6,.8,1), labels=c('','0.2','','0.6','','1'))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext('Temperature covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' b.')), line=.6, cex=1, adj=-0.25)
for (cc in 1:length(bw_cov$temp)) {text(1.5+2*(cc-1), bw_cov$temp[cc], paste(round(bw_cov$temp[cc],digits=3)), pos=3)}
# Sea level
barplot(bw_cov$sealevel, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0,1), xlim=c(0.5,8.5))
axis(2, at=c(0,.2,.4,.6,.8,1), labels=c('','0.2','','0.6','','1'))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext('Sea level covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' c.')), line=.6, cex=1, adj=-0.25)
for (cc in 1:length(bw_cov$sealevel)) {text(1.5+2*(cc-1), bw_cov$sealevel[cc], paste(round(bw_cov$sealevel[cc],digits=3)), pos=3)}
# NAO index
barplot(bw_cov$nao, names.arg=better_names, xlab='', ylab='', space=1, yaxt='n', ylim=c(0,1), xlim=c(0.5,8.5))
axis(2, at=c(0,.2,.4,.6,.8,1), labels=c('','0.2','','0.6','','1'))
mtext(side=1, text='Model', line=2.2)
mtext(side=2, text='BMA weight', line=2.2)
mtext('NAO index covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' d.')), line=.6, cex=1, adj=-0.25)
for (cc in 1:length(bw_cov$nao)) {text(1.5+2*(cc-1), bw_cov$nao[cc], paste(round(bw_cov$nao[cc],digits=3)), pos=3)}
dev.off()

#===============================================================================



#===============================================================================
# FIGURE ???
# Panel for each covariate, return level pdfs and boxplot beneath with quantiles
#===============================================================================

pdf(paste(plot.dir,'returnlevels_pdf_bar.pdf',sep='/'),width=7,height=5.5,colormodel='cmyk')
par(mfrow=c(2,2), mai=c(.6,.50,.4,.1))
yy <- 'y100'
offset <- 0.5
ytop <- 3.5
cc <- 'time' # =======================
plot(pdf.rl[[cc]]$bma[[yy]]$x, offset+pdf.rl[[cc]]$bma[[yy]]$y, type='l', xlim=c(1, 4), ylim=c(0, ytop+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl[[cc]]$gpd4[[yy]]$x, offset+pdf.rl[[cc]]$gpd4[[yy]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd5[[yy]]$x, offset+pdf.rl[[cc]]$gpd5[[yy]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd6[[yy]]$x, offset+pdf.rl[[cc]]$gpd6[[yy]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd3[[yy]]$x, offset+pdf.rl[[cc]]$gpd3[[yy]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(1,1)], q.rl[[cc]]$gpd3[[yy]][c(5,5)]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(2,2)], q.rl[[cc]]$gpd3[[yy]][c(4,4)]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl[[cc]]$gpd3[[yy]][3], q.rl[[cc]]$gpd3[[yy]][3]), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl[[cc]]$bma[[yy]][c(1,1)], q.rl[[cc]]$bma[[yy]][c(5,5)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl[[cc]]$bma[[yy]][c(2,2)], q.rl[[cc]]$bma[[yy]][c(4,4)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl[[cc]]$bma[[yy]][3], q.rl[[cc]]$bma[[yy]][3]), y2, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('Time covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' a.')), line=.6, cex=1, adj=0)
mtext(paste('100-year return level in',proj_years[1],'[m]'), side=1, line=2.4, cex=0.9);
axis(1, at=seq(1, 4, .25), labels=c('1','','1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.05)
legend(2.75, 3, c('ST','NS1','NS2','NS3','BMA'), lty=c(4,2,2,2,1), cex=1, bty='n', lwd=2,
       col=c('seagreen','darkorange3','mediumslateblue','mediumvioletred', 'black'))
#
cc <- 'temp' # =======================
plot(pdf.rl[[cc]]$bma[[yy]]$x, offset+pdf.rl[[cc]]$bma[[yy]]$y, type='l', xlim=c(1, 4), ylim=c(0, ytop+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl[[cc]]$gpd4[[yy]]$x, offset+pdf.rl[[cc]]$gpd4[[yy]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd5[[yy]]$x, offset+pdf.rl[[cc]]$gpd5[[yy]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd6[[yy]]$x, offset+pdf.rl[[cc]]$gpd6[[yy]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd3[[yy]]$x, offset+pdf.rl[[cc]]$gpd3[[yy]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(1,1)], q.rl[[cc]]$gpd3[[yy]][c(5,5)]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(2,2)], q.rl[[cc]]$gpd3[[yy]][c(4,4)]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl[[cc]]$gpd3[[yy]][3], q.rl[[cc]]$gpd3[[yy]][3]), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl[[cc]]$bma[[yy]][c(1,1)], q.rl[[cc]]$bma[[yy]][c(5,5)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl[[cc]]$bma[[yy]][c(2,2)], q.rl[[cc]]$bma[[yy]][c(4,4)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl[[cc]]$bma[[yy]][3], q.rl[[cc]]$bma[[yy]][3]), y2, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('Temperature covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' b.')), line=.6, cex=1, adj=0)
mtext(paste('100-year return level in',proj_years[1],'[m]'), side=1, line=2.4, cex=0.9);
axis(1, at=seq(1, 4, .25), labels=c('1','','1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.05)
#
cc <- 'sealevel' # =======================
plot(pdf.rl[[cc]]$bma[[yy]]$x, offset+pdf.rl[[cc]]$bma[[yy]]$y, type='l', xlim=c(1, 4), ylim=c(0, ytop+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl[[cc]]$gpd4[[yy]]$x, offset+pdf.rl[[cc]]$gpd4[[yy]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd5[[yy]]$x, offset+pdf.rl[[cc]]$gpd5[[yy]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd6[[yy]]$x, offset+pdf.rl[[cc]]$gpd6[[yy]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd3[[yy]]$x, offset+pdf.rl[[cc]]$gpd3[[yy]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(1,1)], q.rl[[cc]]$gpd3[[yy]][c(5,5)]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(2,2)], q.rl[[cc]]$gpd3[[yy]][c(4,4)]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl[[cc]]$gpd3[[yy]][3], q.rl[[cc]]$gpd3[[yy]][3]), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl[[cc]]$bma[[yy]][c(1,1)], q.rl[[cc]]$bma[[yy]][c(5,5)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl[[cc]]$bma[[yy]][c(2,2)], q.rl[[cc]]$bma[[yy]][c(4,4)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl[[cc]]$bma[[yy]][3], q.rl[[cc]]$bma[[yy]][3]), y2, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('Sea level covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' c.')), line=.6, cex=1, adj=0)
mtext(paste('100-year return level in',proj_years[1],'[m]'), side=1, line=2.4, cex=0.9);
axis(1, at=seq(1, 4, .25), labels=c('1','','1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.05)
#
cc <- 'nao' # =======================
plot(pdf.rl[[cc]]$bma[[yy]]$x, offset+pdf.rl[[cc]]$bma[[yy]]$y, type='l', xlim=c(1, 4), ylim=c(0, ytop+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl[[cc]]$gpd4[[yy]]$x, offset+pdf.rl[[cc]]$gpd4[[yy]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd5[[yy]]$x, offset+pdf.rl[[cc]]$gpd5[[yy]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd6[[yy]]$x, offset+pdf.rl[[cc]]$gpd6[[yy]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl[[cc]]$gpd3[[yy]]$x, offset+pdf.rl[[cc]]$gpd3[[yy]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(1,1)], q.rl[[cc]]$gpd3[[yy]][c(5,5)]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl[[cc]]$gpd3[[yy]][c(2,2)], q.rl[[cc]]$gpd3[[yy]][c(4,4)]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl[[cc]]$gpd3[[yy]][3], q.rl[[cc]]$gpd3[[yy]][3]), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl[[cc]]$bma[[yy]][c(1,1)], q.rl[[cc]]$bma[[yy]][c(5,5)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl[[cc]]$bma[[yy]][c(2,2)], q.rl[[cc]]$bma[[yy]][c(4,4)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl[[cc]]$bma[[yy]][3], q.rl[[cc]]$bma[[yy]][3]), y2, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('NAO index covariate', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' d.')), line=.6, cex=1, adj=0)
mtext(paste('100-year return level in',proj_years[1],'[m]'), side=1, line=2.4, cex=0.9);
axis(1, at=seq(1, 4, .25), labels=c('1','','1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.05)
dev.off()

#===============================================================================



#===============================================================================
# FIGURE ???
# Panel for each covariate, return level pdfs and boxplot beneath with quantiles
#===============================================================================

pdf(paste(plot.dir,'returnlevels_bma_pdf_bar.pdf',sep='/'),width=4,height=4,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.50,.25,.25))
yy <- 'y100'
offset <- 0.48
plot(pdf.rl$bma[[yy]]$x, offset+pdf.rl$bma[[yy]]$y, type='l', xlim=c(1.5, 4), ylim=c(0, 4+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl$time$bma[[yy]]$x, offset+pdf.rl$time$bma[[yy]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl$temperature$bma[[yy]]$x, offset+pdf.rl$temperature$bma[[yy]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl$sealevel$bma[[yy]]$x, offset+pdf.rl$sealevel$bma[[yy]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl$nao$bma[[yy]]$x, offset+pdf.rl$nao$bma[[yy]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.41, .47)
polygon(c(q.rl$time$bma[[yy]][c(1,1)], q.rl$time$bma[[yy]][c(5,5)]), c(y1,rev(y1)), col='orange1', border=NA)
polygon(c(q.rl$time$bma[[yy]][c(2,2)], q.rl$time$bma[[yy]][c(4,4)]), c(y1,rev(y1)), col='darkorange2', border=NA)
lines(c(q.rl$time$bma[[yy]][3], q.rl$time$bma[[yy]][3]), y1, lwd=2.5, col='darkorange4')
y2 = c(.32, .38)
polygon(c(q.rl$temp$bma[[yy]][c(1,1)], q.rl$temp$bma[[yy]][c(5,5)]), c(y2,rev(y2)), col='plum2', border=NA)
polygon(c(q.rl$temp$bma[[yy]][c(2,2)], q.rl$temp$bma[[yy]][c(4,4)]), c(y2,rev(y2)), col='slateblue3', border=NA)
lines(c(q.rl$temp$bma[[yy]][3], q.rl$temp$bma[[yy]][3]), y2, lwd=2.5, col='slateblue4')
y3 = c(.23, .29)
polygon(c(q.rl$sealevel$bma[[yy]][c(1,1)], q.rl$sealevel$bma[[yy]][c(5,5)]), c(y3,rev(y3)), col='salmon', border=NA)
polygon(c(q.rl$sealevel$bma[[yy]][c(2,2)], q.rl$sealevel$bma[[yy]][c(4,4)]), c(y3,rev(y3)), col='violetred', border=NA)
lines(c(q.rl$sealevel$bma[[yy]][3], q.rl$sealevel$bma[[yy]][3]), y3, lwd=2.5, col='red4')
y4 = c(.14, .2)
polygon(c(q.rl$nao$bma[[yy]][c(1,1)], q.rl$nao$bma[[yy]][c(5,5)]), c(y4,rev(y4)), col='palegreen1', border=NA)
polygon(c(q.rl$nao$bma[[yy]][c(2,2)], q.rl$nao$bma[[yy]][c(4,4)]), c(y4,rev(y4)), col='seagreen3', border=NA)
lines(c(q.rl$nao$bma[[yy]][3], q.rl$nao$bma[[yy]][3]), y4, lwd=2.5, col='darkgreen')
y5 = c(.05, .11)
polygon(c(q.rl$bma[[yy]][c(1,1)], q.rl$bma[[yy]][c(5,5)]), c(y5,rev(y5)), col='gray60', border=NA)
polygon(c(q.rl$bma[[yy]][c(2,2)], q.rl$bma[[yy]][c(4,4)]), c(y5,rev(y5)), col='gray30', border=NA)
lines(c(q.rl$bma[[yy]][3], q.rl$bma[[yy]][3]), y5, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext(paste('100-year return level in',proj_years[1],'[m]'), side=1, line=2.4, cex=0.9);
axis(1, at=seq(1, 4, .25), labels=c('1','','1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.05)
legend(2.8, 3.7, c('Time','Temperature','Sea level','NAO','BMA'), lty=c(4,2,2,2,1), cex=1, bty='n', lwd=2,
       col=c('darkorange3','mediumslateblue','mediumvioletred','seagreen','black'))
dev.off()

#===============================================================================




#===============================================================================
# End
#===============================================================================
