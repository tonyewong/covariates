#===============================================================================
# get_timeseries_covariates.R
#
# First, everything is normalized to first 20 years (1928-1937) of the tide
# gauge data. Then, the covariates are each normalized so that the min-max
# range for the hindcast period is 0-1.
#
# Questions? Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================

data_calib <- readRDS('../data/tidegauge_processed_norfolk-delfzijl_decl3-pot99-annual_10Jul2018.rds')

years <- data_calib$norfolk$y89$year


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
SL.new = SL.new - mean(SL.new[ibeg:iend])   # SL data are relative to 1961-1990
sealevel <- cbind(SL.time, SL.new)
colnames(sealevel) <- c('year','sealevel')


# just use Time (year relative to first year?) =================================
time <- cbind(years, years)
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

# trim the data_calib object
data_calib <- data_calib$norfolk$y89
ind_keep <- which(data_calib$year >= min_year & data_calib$year <= max_year)
data_calib$year <- data_calib$year[ind_keep]
data_calib$counts <- data_calib$counts[ind_keep]
data_calib$excesses[-ind_keep] <- NULL
data_calib$excesses_all <- NULL # make sure this doesn't ruin anything
data_calib$time_length <- data_calib$time_length[ind_keep]
data_calib$time_length_all <- sum(data_calib$time_length)


#===============================================================================
# create a single array to hold all the possible covariates

names_covariates <- c('time','temp','sealevel','nao')
covariates <- cbind(time[,2], temperature[,2], sealevel[,2], nao[,2])
colnames(covariates) <- names_covariates


#===============================================================================
#===============================================================================
# get time series for forcing/making return level projections


# CNRM temperature projections =================================================
ncdata <- nc_open('../data/global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc')
   temperature_proj <- ncvar_get(ncdata, 'tas')
   time_proj <- ncvar_get(ncdata, 'time')
nc_close(ncdata)

# normalize like covariates for calibration
year_proj <- floor(time_proj/10000)
ibeg <- which(year_proj==data_calib$year[1])
iend <- which(year_proj==data_calib$year[20])
temperature_proj <- temperature_proj - mean(temperature_proj[ibeg:iend])
temperature_proj <- cbind(year_proj, temperature_proj)
colnames(temperature_proj) <- c('year','temp')


# Time =========================================================================
timec_proj <- temperature_proj[,'year'] - mean(temperature_proj[ibeg:iend,'year'])
timec_proj <- cbind(temperature_proj[,'year'], timec_proj)
colnames(timec_proj) <- c('year','time')


# BRICK sea level projection (from Wong and Keller 2017) =======================
ncdata <- nc_open("../data/BRICK_GMSL_WongKeller2017.nc")
gmsl_rcp85 <- t(ncvar_get(ncdata, 'GlobalSeaLevel_RCP85'))
time_proj <- ncvar_get(ncdata, 'time_proj')
nc_close(ncdata)

# get projection that yields median in 2100 under RCP8.5
ibeg <- 1
iend <- nrow(gmsl_rcp85)
# need only odd number so median = element
if (iend%%2 == 0) {iend <- iend-1}
ind_median <- which(gmsl_rcp85[,length(time_proj)]==median(gmsl_rcp85[ibeg:iend,length(time_proj)]))
gmsl_proj <- cbind(time_proj, gmsl_rcp85[ind_median,])
colnames(gmsl_proj) <- c('year','sealevel')

# normalize as in hindcast
ibeg <- which(time_proj==data_calib$year[1])
iend <- which(time_proj==data_calib$year[20])
gmsl_proj[,2] <- gmsl_proj[,2] - mean(gmsl_proj[ibeg:iend,2])


# NAO index ====================================================================
# (using method of Stephenson et al 2006, as in Wong et al 2018)
file.in <- '../data/DMIEH5_SRA1B_4_MM_psl.1-1200.nc'
ncdata <- nc_open(file.in)
  psl <- ncvar_get(ncdata, 'psl')
  lon <- ncvar_get(ncdata, 'lon')
  lat <- ncvar_get(ncdata, 'lat')
  time <- ncvar_get(ncdata, 'time')  # hours after 2001-01-31 (2001-2100 data)
nc_close(ncdata)

n_month <- length(time)
n_year <- n_month/12 -1

# As in Stephenson et al, 2006 (doi:  10.1007/s00382-006-0140-x)
# use regional, since there is disagreement over how exactly to do the PCA
# (which EOFs to rotate, e.g.) and a physical interpretation is direct this way
ilat_azores <- which(lat >= 20 & lat <= 55)
ilon_azores <- which(lon >=(360-90) | lon <= 60)

ilat_iceland <- which(lat >= 55 & lat <= 90)
ilon_iceland <- which(lon >=(360-90) | lon <= 60)

psl_azores <- psl[ilon_azores, ilat_azores, ]
psl_iceland <- psl[ilon_iceland, ilat_iceland, ]

# Don't do area-weighting as per explanation from Jesse Nusbaumer 27 March 2018
# email. Jesse says:
#   Don't area-weight the results.  Instead, just take the average of the SLP
#   over the specified region, treating every grid box the same.  The reason is
#   because in GCM papers they will usually specifically state "area-weighted"
#   if the actual meters-squared area is being taken into account, otherwise it
#   is assumed to just be a non-weighted average.  Also when it comes to a lot
#   of these indices the actual physical values (e.g. the total amount of
#   atmospheric mass producing the surface pressure) isn't really important,
#   they are just using the values to produce a strong statistical relationship.

# take the bulk average over each area, for each month
psl_azores_spavg <- apply(psl_azores, 3, mean)
psl_iceland_spavg <- apply(psl_iceland, 3, mean)

# Standardize each site separately (as discussed in Jones et al 1997). Do relative
# to 2001-2016 mean/stdev so we can be consistent between projections and the
# historical record
psl_azores_std <- rep(-999, n_month)
psl_iceland_std <- rep(-999, n_month)
for (m in 1:12) {
  ind_this_month <- seq(from=m, to=n_month, by=12)
  # first 16 are 2001-2016
  psl_tmp <- psl_azores_spavg[ind_this_month]
  psl_azores_std[ind_this_month] <- (psl_tmp - mean(psl_tmp[1:16]))/sd(psl_tmp[1:16])
  psl_tmp <- psl_iceland_spavg[ind_this_month]
  psl_iceland_std[ind_this_month] <- (psl_tmp - mean(psl_tmp[1:16]))/sd(psl_tmp[1:16])
}

# get SLP difference
nao_monthly <- psl_azores_std - psl_iceland_std

# get winter mean
nao_proj <- rep(-999, n_year)
for (y in 1:n_year) {
  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:14])  # DJF
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:15])  # DJFM
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 1:12])   # annual
}

time_proj <- 2001:2099

# paste together with the nao (historical)
ibeg <- which(time_proj==(max(nao[,1])+1))
iend <- length(time_proj)
nao_proj <- cbind(c(nao[,1], time_proj[ibeg:iend]),
                  c(nao[,2], nao_proj[ibeg:iend]))

# and normalize all together
ibeg <- which(nao_proj[,1]==data_calib$year[1])
iend <- which(nao_proj[,1]==data_calib$year[20])
nao_proj[,2] <- nao_proj[,2] - mean(nao_proj[ibeg:iend,2])
colnames(nao_proj) <- c('year','nao')


# put all together as in covariates ============================================

# trim to only the period we have all of them
min_year <- max( min(nao_proj[,'year']),  min(temperature_proj[,'year']),
                 min(gmsl_proj[,'year']), min(timec_proj[,'year'])        )
max_year <- min( max(nao_proj[,'year']),  max(temperature_proj[,'year']),
                 max(gmsl_proj[,'year']), max(timec_proj[,'year'])        )

year_proj <- min_year:max_year

ibeg <- which(timec_proj[,'year']==min_year)
iend <- which(timec_proj[,'year']==max_year)
timec_proj <- timec_proj[ibeg:iend,]
ibeg <- which(nao_proj[,'year']==min_year)
iend <- which(nao_proj[,'year']==max_year)
nao_proj <- nao_proj[ibeg:iend,]
ibeg <- which(temperature_proj[,'year']==min_year)
iend <- which(temperature_proj[,'year']==max_year)
temperature_proj <- temperature_proj[ibeg:iend,]
ibeg <- which(gmsl_proj[,'year']==min_year)
iend <- which(gmsl_proj[,'year']==max_year)
gmsl_proj <- gmsl_proj[ibeg:iend,]

covariates_proj <- cbind(year_proj, timec_proj[,2], temperature_proj[,2], gmsl_proj[,2], nao_proj[,2])
colnames(covariates_proj) <- c('year',names_covariates) # time temp sealevel nao


# first, normalize the covariates also relative to the first 20 years
# of available tide gauge data
for (cc in names_covariates) {
  covariates[,cc] <- covariates[,cc] - mean(covariates[1:20, cc])
  covariates_proj[,cc] <- covariates_proj[,cc] - mean(covariates_proj[1:20, cc])
}

# now paste in the historical data for the period it's available for all time
# series (1928-2013)
ibeg <- 1
iend <- nrow(covariates)
covariates_proj[ibeg:iend, names_covariates] <- covariates

# finally, normalize so the hindcast period (covariates) are on 0-1
norms <- array(NA, c(length(names_covariates),2))
rownames(norms) <- names_covariates
colnames(norms) <- c('min','max')

for (cc in names_covariates) {
  # then, between 0 and 1
  norms[cc,'min'] <- min(covariates[,cc])
  norms[cc,'max'] <- max(covariates[,cc])
  covariates[,cc] <- (covariates[,cc] - norms[cc,'min'])/(norms[cc,'max'] - norms[cc,'min'])
  covariates_proj[,cc] <- (covariates_proj[,cc] - norms[cc,'min'])/(norms[cc,'max'] - norms[cc,'min'])
}


#===============================================================================
# End
#===============================================================================
