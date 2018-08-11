#===============================================================================
# Processing of Norfolk, Virginia, USA, and Delfzijl, Netherlands data to yield
# a time series of 99th percentiles of declustered daily maximum sea levels.
#
# Implicit inputs:
#   days.daily.max
#   years.daily.max
#   sl.daily.max
#
#===============================================================================

load('processing_inprog_norfolk.RData')

ts.percentile <- 0.99

ts.years <- unique(years.daily.max)

min.data <- 0.9*365.25 # only take years that have enough data points

ts.slmax <- rep(NA, length(ts.years))

for (y in ts.years) {
  ind.thisyear <- which(years.daily.max==y)
  if (length(ind.thisyear) >= min.data) {
    ts.slmax[match(y,ts.years)] <- quantile(sl.daily.max[ind.thisyear], ts.percentile)
  }
}


#===============================================================================
## Get other time series and look at these as predictors for ts.slmax in a
## multivariate regression sense (or BSTS?)


# get NAO index ================================================================
nao_dat <- read.table('../data/nao_3dp.dat')
colnames(nao_dat) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','ann')
ibeg <- which(nao_dat['year']==1850)
iend <- max(which(nao_dat[,'ann']!=-99.99))
time_hist <- nao_dat[ibeg:iend, 'year']

# get DJF means
nao_hist <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao_hist[y] <- mean( c(nao_dat$dec[ibeg+y-1], nao_dat$jan[ibeg+y], nao_dat$feb[ibeg+y]) )
}
time_hist <- nao_dat[ibeg:iend, 'year']

# re-normalize (mean and stdev) relative to 2001-2016 mean/stdev, so it is
# consistent with the projections
ind_norm <- which(time_hist==2001):which(time_hist==2016)
nao_hist <- (nao_hist - mean(nao_hist[ind_norm]))/sd(nao_hist[ind_norm])

nao <- cbind(time_hist, nao_hist)
colnames(nao) <- c('year','nao')


# get GMST =====================================================================
data.tmp <- read.table('../data/noaa_temperature_1880-2017.csv', header = TRUE, sep=',')
time_hist <- data.tmp$Year
temperature_hist <- data.tmp$Value
temperature <- cbind(time_hist, temperature_hist)
colnames(temperature) <- c('year','temperature')


# get GMSL =====================================================================
sl.dat.new = read.table("../data/GMSL_ChurchWhite2011_yr_2015.txt")     #Reconstructed sea-level
SL.time= sl.dat.new[,1]-0.5     # times are in half-year
SL.new = sl.dat.new[,2]/1000    # data are in mm
SL.err = sl.dat.new[,3]/1000
ibeg=which(SL.time==1961); iend=which(SL.time==1990);
SL.new = SL.new - mean(SL.new[ibeg:iend])               # make sure SL data are relative to 1961-1990
sealevel <- cbind(SL.time, SL.new)
colnames(sealevel) <- c('year','sealevel')


# just use Time (year relative to first year?) =================================
time <- cbind(ts.years, ts.years)
colnames(time) <- c('year','year')


# trim all potential features to ts.years
min_year <- max( min(nao[,'year']),      min(temperature[,'year']),
                 min(sealevel[,'year']), min(time[,'year'])        )
max_year <- min( max(nao[,'year']),      max(temperature[,'year']),
                 max(sealevel[,'year']), max(time[,'year'])        )

nao <- nao[which(nao[,1]==min_year):which(nao[,1]==max_year), ]
temperature <- temperature[which(temperature[,1]==min_year):which(temperature[,1]==max_year), ]
sealevel <- sealevel[which(sealevel[,1]==min_year):which(sealevel[,1]==max_year), ]
time <- time[which(time[,1]==min_year):which(time[,1]==max_year), ]
ts.slmax <- ts.slmax[which(ts.years==min_year):which(ts.years==max_year)]
ts.years <- ts.years[which(ts.years==min_year):which(ts.years==max_year)]


## Fit multivariate linear regression
data_calib <- data.frame(cbind(ts.slmax, nao[,2], temperature[,2], sealevel[,2], time[,2]))
colnames(data_calib) <- c('slmax','nao','temperature','sealevel','year')
fit <- lm(ts.slmax ~ nao + temperature + sealevel + year, data=data_calib)

#===============================================================================
# Okay, so multivariate linear regression does not work so well.
# Instead, fit BMA ensemble for each of the potential covariates and get the
# AIC/BIC/etc for each of the BMA ensemble max posterior score
# --> using PP/GPD and full data
# --> fit ST, NS1, NS2, NS3 for each of the covariates
# -->



#===============================================================================
# End
#===============================================================================
