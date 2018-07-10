#===============================================================================
# processing of Delfzijl, Netherlands data.
#
#
# Inputs:
#   dt.decluster      declustering time scale (in days)
#   detrend.method    either 'linear' or 'annual' [block means]
#   pot.threshold     percentile (0-1) for GPD peaks-over-thresholds cutoff
#
#===============================================================================

print('starting to process Delfzijl tide gauge data')

filename.saveprogress <- '../output/processing_delfzijl.RData'
print(paste('saving intermediate progress into file:',filename.saveprogress))

#===
#=== read in tide gauge data
#===

data <- read.table('../data/id1-DELFZL-187901010000-201607072359_reduced.csv', header = TRUE, sep=',')
names(data)[1:3] <- c('date','time','sl')
data$sl <- data$sl*10 # convert to mm from cm

# separate date-stamp into year / month / day / hour / minute
times.char  <- as.character(data$time)
data$hour   <- as.numeric(unlist(strsplit(times.char, split=":"))[2*seq(from=1, to=length(times.char), by=1)-1])
data$minute <- as.numeric(unlist(strsplit(times.char, split=":"))[2*seq(from=1, to=length(times.char), by=1)])

# time in days since 01 January 1960
data$time.days <- as.numeric(mdy.date(month=data$month, day=data$day, year=data$year)) + data$hour/24 + data$minute/(24*60)

#
#===============================================================================
# daily POT series, for PP-GPD model
#===============================================================================
#

tbeg <- proc.time()

# difference between time stamps (in units of days)
time.diff <- diff(data$time.days)

# check that they are in the proper order, ascending
print(paste('Are there any times out of order? ',any(time.diff < 0), sep=''))

# put everything back in order - make sure you do this for the sea levels and
# other fields too, so do as a matrix. also recalculate the time differences,
# which you will need for averaging
data <- data[order(data$time.days),]
time.diff <- diff(data$time.days)

#===
#=== first, need to *subsample* up to 3-hourly time series (c.f. email with Rijkwaterstaat, Alexander Bakker)
#===

# what is three hours? in units of days
three.hours <- 3/24

# where are there gaps longer than three hours? (+10sec for precision)
igap <- which(time.diff > (three.hours+10/(24*60*60)))

# where are there gaps shorter than three hours? (-10sec)
# should find lots of places. need to average these up to 3 hours.
iokay <- which( (time.diff < (three.hours+10/(24*60*60))) & (time.diff > (three.hours-10/(24*60*60))) )
ishort <- which(time.diff < (three.hours-10/(24*60*60)))

# turns out with the Rijkwaterstaat data set, have a lot of different time
# intervals, but no gaps longer than 3 hours.  upscale the less-than-3-hourly
# data to 3-hourly

# in a perfect world, these are the times at which we have obs, and we just
# need to pluck out the observations at these time points
time_beg <- min(data$time.days)
time_end <- max(data$time.days)
delta_t <- 3/24   # time step of 3 hours, units of days
time_3hour <- seq(from=time_beg, to=time_end, by=delta_t)

# go through the entire time.days vector; you know there are no obs gaps longer
# than 3 hours, so each of the 3-hour intervals has at least one observation in
# it. grab the one closest to the mark (without going over?)
# warning: this will take a while (~1 hour, maybe a bit more, on laptop)

ind_map3hour <- sapply(1:length(time_3hour), function(t) {which.min(abs(time_3hour[t]-data$time.days))})
sl_3hour <- data$sl[ind_map3hour]
year_3hour <- data$year[ind_map3hour]

# save progress
save.image(file=filename.saveprogress)

#===
#=== Detrend by either subtracting linear sea-level trend (fit to monthly
#=== means) or by subtracting annual means (moving 1-year window)
#===

print(paste('Detrending using method `',detrend.method,'` ...', sep=''))

if (detrend.method=='linear') {

  # calculate monthly means

  dates.new <- date.mdy(time_3hour)
  date.beg <- date.mdy(min(time_3hour))
  date.end <- date.mdy(max(time_3hour))

  # what the years in which we have data?
  years.unique <- unique(dates.new$year)

  # in each year, what are the months with at least 90% of the data?
  months.this.year <- vector('list', length(years.unique))
  names(months.this.year) <- years.unique
  years.to.remove <- NULL
  for (year in years.unique) {
    ind.this.year <- which(dates.new$year == year)
    months.to.keep <- NULL
    for (month in 1:12) {
      ind.this.month <- which(dates.new$month[ind.this.year] == month)
      days.this.month <- monthDays(paste(year,'-',month,'-','1', sep=''))
      hours.this.month <- days.this.month * 24
      # *3 because these are 3-hourly blocks
      perc.data.this.month <- 3*length(ind.this.month)/hours.this.month
      if (perc.data.this.month >= 0.9) {months.to.keep <- c(months.to.keep, month)}
    }
    if(length(months.to.keep)>0) {months.this.year[[year]] <- months.to.keep }
    else                         {years.to.remove <- c(years.to.remove, year)}
  }
  if(length(years.to.remove)>0) {years.unique <- years.unique[-match(years.to.remove, years.unique)]}

  # get the mean time (in days releative to 1 Jan 1960) of the observations of
  # each month we are using to fit the trend for SLR
  times.month <- rep(NA, length(unlist(months.this.year)))
  sl.month    <- rep(NA, length(unlist(months.this.year)))
  cnt <- 1
  for (year in years.unique) {
    ind.this.year <- which(dates.new$year == year)
    for (month in months.this.year[[year]]) {
      ind.this.month <- which(dates.new$month[ind.this.year] == month)
      times.month[cnt] <- mean(data$time.days[ind.this.year[ind.this.month]])
      sl.month[cnt]    <- mean(data$sl[ind.this.year[ind.this.month]])
      cnt <- cnt + 1
    }
  }

  # fit trend
  slr.trend <- lm(sl.month ~ times.month)
  slr.trend.3hour <- slr.trend$coefficients[1] + slr.trend$coefficients[2]*data$time.days

  # subtract off from the 1-hourly data
  data$sl.detrended <- data$sl - slr.trend.3hour

} else if(detrend.method=='annual') {

  # what the years in which we have data?
  dates.new <- date.mdy(time_3hour)
  years.unique <- unique(dates.new$year)

  # get a placeholder -- want to be using the 3-hourly time series
  sl_3hour_detrended <- sl_3hour
  time.days.beg <- min(time_3hour)
  time.days.end <- max(time_3hour)

  pb <- txtProgressBar(min=0,max=length(time_3hour),initial=0,style=3)
  for (tt in 1:length(time_3hour)) {
    # if within half a year of either end of the time series, include either the
    # entire first year or entire last year to get a full year's worth of data in
    # the subtracted mean
    if (time_3hour[tt] - time.days.beg < (365.25*0.5)) {
      ind.close <- which(time_3hour - time.days.beg <= 365.25)
    } else if(time.days.end - time_3hour[tt] < (365.25*0.5)) {
      ind.close <- which(time.days.end - time_3hour <= 365.25)
    } else {
      ind.close <- which(abs(time_3hour-time_3hour[tt]) <= (365.25*0.5) )
    }
    sl_3hour_detrended[tt] <- sl_3hour[tt] - mean(sl_3hour[ind.close])
    setTxtProgressBar(pb, tt)
  }
  close(pb)
} else {
  print('ERROR: unknown detrend.method value')
}

print('  ... done.')

# save progress
save.image(file=filename.saveprogress)

#===
#=== daily block maxima; calculate 99% quantile as GPD threshold
#===

# how many days in each year have at least 90% of their values?
days.all <- floor(data$time.days)
days.unique <- unique(days.all)
ind.days.to.remove <- NULL
print('... filtering down to do a daily maxima time series of only the days with at least 90% of data ...')
pb <- txtProgressBar(min=min(days.unique),max=max(days.unique),initial=0,style=3)
for (day in days.unique) {
  ind.today <- which(floor(data$time.days) == day)
  perc.data.today <- length(ind.today)/24
  if(perc.data.today < 0.9) {ind.days.to.remove <- c(ind.days.to.remove, match(day, days.unique))}
  setTxtProgressBar(pb, day)
}
close(pb)
days.daily.max <- days.unique[-ind.days.to.remove]
n.days <- length(days.daily.max)

# calculate the daily maximum sea levels on the days of 'days.daily.max'
sl.daily.max <- rep(NA, n.days)
years.daily.max <- rep(NA, n.days)
print('... calculating time series of daily maxima ...')
pb <- txtProgressBar(min=0,max=n.days,initial=0,style=3)
for (day in days.daily.max) {
  cnt <- match(day,days.daily.max)
  ind.today <- which(days.all == day)
  sl.daily.max[cnt] <- max(data$sl.detrended[ind.today])
  years.daily.max[cnt] <- data$year[ind.today][1]
  setTxtProgressBar(pb, cnt)
}
close(pb)

# save progress
save.image(file=filename.saveprogress)

#===
#=== find all the excesses, "declustering" = if two are within dt.decluster of
#=== each other, take only the maximum of the two (so make sure you save the
#=== times of each excess)
#===

# add in data dt years at a time
dt <- 1
year.beg <- min(years.daily.max)
year.end <- max(years.daily.max)
delfzijl_total <- year.end - year.beg + 1
delfzijl_years <- rev(seq(from=delfzijl_total, to=10, by=-dt))
# the rev() accounts for possibly not ending at 10 years; we want to include
# the longest data set possible, and would be fine to leave out the shorter
# subset experiment

len_delfzijl <- NULL
for (yy in delfzijl_years) {
  len_delfzijl <- c(len_delfzijl, paste('y', yy,sep=''))
}
names(delfzijl_years) <- len_delfzijl

data_delfzijl <- vector('list', length(delfzijl_years))
names(data_delfzijl) <- len_delfzijl

print('... getting threshold excesses ...')

for (ll in 1:length(delfzijl_years)) {
  data_delfzijl[[ll]]$dt.decluster <- dt.decluster

  # find all the daily maxima from this time subset
  ind.subset <- which((years.daily.max > (year.end - delfzijl_years[ll])) & (years.daily.max <= year.end))
  sl.daily.max.subset <- sl.daily.max[ind.subset]
  days.daily.max.subset <- days.daily.max[ind.subset]
  years.daily.max.subset <- years.daily.max[ind.subset]
  years.unique <- unique(years.daily.max.subset)
  data_delfzijl[[ll]]$p.threshold <- pot.threshold
  data_delfzijl[[ll]]$threshold <- as.numeric(quantile(sl.daily.max.subset, pot.threshold))

  # get the exceedances data
  ind.exceed <- which(sl.daily.max.subset >= data_delfzijl[[ll]]$threshold)
  days.exceed <- days.daily.max.subset[ind.exceed]
  sl.exceed <- sl.daily.max.subset[ind.exceed]
  years.exceed <- years.daily.max.subset[ind.exceed]

  # declustering
  declustered.exceed <- decluster_timeseries(time=days.exceed, year=years.exceed,
                                             time.series=sl.exceed, min.dt=data_delfzijl[[ll]]$dt.decluster)
  days.exceed.decl <- declustered.exceed$time
  years.exceed.decl <- declustered.exceed$year
  sl.exceed.decl <- declustered.exceed$time.series

  # data to list

  data_delfzijl[[ll]]$counts <- data_delfzijl[[ll]]$year <- data_delfzijl[[ll]]$time_length <- rep(NA, length(years.unique))
  data_delfzijl[[ll]]$excesses <- vector('list', length(years.unique))

  for (yy in 1:length(years.unique)) {
    ind.hits.this.year <- which(years.exceed.decl == years.unique[yy])
    data_delfzijl[[ll]]$counts[yy] <- length(ind.hits.this.year)
    data_delfzijl[[ll]]$year[yy] <- years.unique[yy]
    data_delfzijl[[ll]]$time_length[yy] <- length(which(years.daily.max == years.unique[yy]))
    if(length(ind.hits.this.year) > 0) {data_delfzijl[[ll]]$excesses[[yy]] <- sl.exceed.decl[ind.hits.this.year]
    } else                             {data_delfzijl[[ll]]$excesses[[yy]] <- NA}
  }

  # alternatively, could bin em all together. but this won't allow for potential
  # non-stationary behavior in the poisson process
  data_delfzijl[[ll]]$excesses_all <- sl.exceed.decl
  data_delfzijl[[ll]]$counts_all <- length(sl.exceed.decl)
  data_delfzijl[[ll]]$time_length_all <- length(days.daily.max.subset)

}

tend <- proc.time()
print(paste('  ... done. Took ', (tend[3]-tbeg[3])/60, ' minutes.',sep=''))

# save final 'data_delfzijl' object to RDS to use later
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_delfzijl_decl',dt.decluster,'-pot',pot.threshold*100,'-',detrend.method,'_',today,'.rds', sep='')
print(paste('Saving processed data to file:',filename.output))
saveRDS(data_delfzijl, file=filename.output)

#===============================================================================

print('done processing the Sewells Point/delfzijl tide gauge data set')

#===============================================================================
# End
#===============================================================================
