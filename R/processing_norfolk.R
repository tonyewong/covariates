#===============================================================================
# processing of Norfolk, Virginia, USA, data.
#
#
# Inputs:
#   dt.decluster      declustering time scale (in days)
#   detrend.method    either 'linear' or 'annual' [block means]
#   pot.threshold     percentile (0-1) for GPD peaks-over-thresholds cutoff
#
#===============================================================================

print('starting to process Norfolk (Sewells Point) tide gauge data')

#===
#=== read in tide gauge data, they're all already hourly series
#===

dat.dir <- '../data/tide_gauge_SewellsPoint/'
filetype <- 'txt'
septype <- '\t'

files.tg <- list.files(path=dat.dir,pattern=filetype)

data <- read.table(paste(dat.dir,files.tg[1],sep=''), header = TRUE, sep=septype)
if(length(files.tg) > 1) {
    for (ff in 2:length(files.tg)) {
        data <- rbind(data, read.table(paste(dat.dir,files.tg[ff],sep=''), header = TRUE, sep=septype))
    }
}

# convert sea levels from m to mm, consistent with the other data
data$sl <- 1000* data$sl

data$year   <- as.numeric(substr(as.character(data$date), start=1, stop=4))
data$month  <- as.numeric(substr(as.character(data$date), start=5, stop=6))
data$day    <- as.numeric(substr(as.character(data$date), start=7, stop=8))
# where is the ':' in the times? tells us where to separate off hte substring
# giving the hour of the observations.
colon.location <- regexpr(':', data$time)
ind.one.digit <- which(colon.location==2)
ind.two.digit <- which(colon.location==3)
data$hour <- rep(NA, length(data$day))
data$hour[ind.one.digit] <- as.numeric(substr(data$time[ind.one.digit], 1,1))
data$hour[ind.two.digit] <- as.numeric(substr(data$time[ind.two.digit], 1,2))

# time in days since 01 January 1960
data$time.days <- as.numeric(mdy.date(month=data$month, day=data$day, year=data$year)) + data$hour/24

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
# Note - none for Sewells Point/Norfolk data
#data <- data[order(data$time.days),]
#time.diff <- diff(data$time.days)

# what is one hour? in units of days
one.hours <- 1/24

# where are there gaps longer than one hours? (+10sec for precision)
igap <- which(time.diff > (one.hours+10/(24*60*60)))

#===
#=== Detrend by either subtracting linear sea-level trend (fit to monthly
#=== means) or by subtracting annual means (moving 1-year window)
#===

print(paste('Detrending using method `',detrend.method,'` ...', sep=''))

# what the years in which we have data?
dates.new <- date.mdy(data$time.days)
years.unique <- unique(dates.new$year)

# get a placeholder
data$sl.detrended <- data$sl
time.days.beg <- min(data$time.days)
time.days.end <- max(data$time.days)

pb <- txtProgressBar(min=0,max=length(data$time.days),initial=0,style=3)
for (tt in 1:length(data$time.days)) {
  # if within half a year of either end of the time series, include either the
  # entire first year or entire last year to get a full year's worth of data in
  # the subtracted mean
  if (data$time.days[tt] - time.days.beg < (365.25*0.5)) {
    ind.close <- which(data$time.days - time.days.beg <= 365.25)
  } else if(time.days.end - data$time.days[tt] < (365.25*0.5)) {
    ind.close <- which(time.days.end - data$time.days <= 365.25)
  } else {
    ind.close <- which(abs(data$time.days-data$time.days[tt]) <= (365.25*0.5) )
  }
  data$sl.detrended[tt] <- data$sl[tt] - mean(data$sl[ind.close])
  setTxtProgressBar(pb, tt)
}
close(pb)

print('  ... done.')

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

#===
#=== find all the excesses, "declustering" = if two are within dt.decluster of
#=== each other, take only the maximum of the two (so make sure you save the
#=== times of each excess)
#===

# add in data dt years at a time
dt <- 1
year.beg <- min(years.daily.max)
year.end <- max(years.daily.max)
norfolk_total <- year.end - year.beg + 1
norfolk_years <- rev(seq(from=norfolk_total, to=10, by=-dt))
# the rev() accounts for possibly not ending at 10 years; we want to include
# the longest data set possible, and would be fine to leave out the shorter
# subset experiment

len_norfolk <- NULL
for (yy in norfolk_years) {
  len_norfolk <- c(len_norfolk, paste('y', yy,sep=''))
}
names(norfolk_years) <- len_norfolk

data_norfolk <- vector('list', length(norfolk_years))
names(data_norfolk) <- len_norfolk

print('... getting threshold excesses ...')

for (ll in 1:length(norfolk_years)) {
  data_norfolk[[ll]]$dt.decluster <- dt.decluster

  # find all the daily maxima from this time subset
  ind.subset <- which((years.daily.max > (year.end - norfolk_years[ll])) & (years.daily.max <= year.end))
  sl.daily.max.subset <- sl.daily.max[ind.subset]
  days.daily.max.subset <- days.daily.max[ind.subset]
  years.daily.max.subset <- years.daily.max[ind.subset]
  years.unique <- unique(years.daily.max.subset)
  data_norfolk[[ll]]$p.threshold <- pot.threshold
  data_norfolk[[ll]]$threshold <- as.numeric(quantile(sl.daily.max.subset, pot.threshold))

  # get the exceedances data
  ind.exceed <- which(sl.daily.max.subset >= data_norfolk[[ll]]$threshold)
  days.exceed <- days.daily.max.subset[ind.exceed]
  sl.exceed <- sl.daily.max.subset[ind.exceed]
  years.exceed <- years.daily.max.subset[ind.exceed]

  # declustering
  declustered.exceed <- decluster_timeseries(time=days.exceed, year=years.exceed,
                                             time.series=sl.exceed, min.dt=data_norfolk[[ll]]$dt.decluster)
  days.exceed.decl <- declustered.exceed$time
  years.exceed.decl <- declustered.exceed$year
  sl.exceed.decl <- declustered.exceed$time.series

  # data to list

  data_norfolk[[ll]]$counts <- data_norfolk[[ll]]$year <- data_norfolk[[ll]]$time_length <- rep(NA, length(years.unique))
  data_norfolk[[ll]]$excesses <- vector('list', length(years.unique))

  for (yy in 1:length(years.unique)) {
    ind.hits.this.year <- which(years.exceed.decl == years.unique[yy])
    data_norfolk[[ll]]$counts[yy] <- length(ind.hits.this.year)
    data_norfolk[[ll]]$year[yy] <- years.unique[yy]
    data_norfolk[[ll]]$time_length[yy] <- length(which(years.daily.max == years.unique[yy]))
    if(length(ind.hits.this.year) > 0) {data_norfolk[[ll]]$excesses[[yy]] <- sl.exceed.decl[ind.hits.this.year]
    } else                             {data_norfolk[[ll]]$excesses[[yy]] <- NA}
  }

  # alternatively, could bin em all together. but this won't allow for potential
  # non-stationary behavior in the poisson process
  data_norfolk[[ll]]$excesses_all <- sl.exceed.decl
  data_norfolk[[ll]]$counts_all <- length(sl.exceed.decl)
  data_norfolk[[ll]]$time_length_all <- length(days.daily.max.subset)

}

tend <- proc.time()
print(paste('  ... done. Took ', (tend[3]-tbeg[3])/60, ' minutes.',sep=''))

# save final 'data_norfolk' object to RDS to use later
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_norfolk_decl',dt.decluster,'-pot',pot.threshold*100,'-',detrend.method,'_',today,'.rds', sep='')
print(paste('Saving processed data to file:',filename.output))
saveRDS(data_norfolk, file=filename.output)

#===============================================================================

print('done processing the Sewells Point/Norfolk tide gauge data set')

#===============================================================================
# End
#===============================================================================
