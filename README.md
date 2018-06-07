---
# datalengths

Questions? Tony Wong <anthony.e.wong@colorado.edu>

---
## Pipeline

### 1 Processing

Subtract moving window 1-year mean, etc TODO

Threshold selection, etc

New data set for each interval `dyear` between 30 years and full data.  Do not use each year because some years will have 0 exceedances.  Instead, calculate the minimum number of years (integer) so that each added block will add at least one data point.

### 2 Prior distributions

Using only stationary GPD, priors by 

### 3 Calibration

MCMC...

Key figure:  box-whisker plots of distributions across data lengths.  do as a surface, looks cool!

### 4 Kolmogorov-Smirnov test 

Asks the question:  are the resulting posteriors for a given return level but different data length experiment from the same distribution?

D = max |E1(i) - E2(i)|

empirical CDF

Key figure:  test statistic D as number of years of data increases from 30 -> Full (dashed line for accept/reject boundary)

### 5 Kullback-Leibler divergence

Asks the question:  what is the value of the additional information as new data becomes available?

Only assess this after the KS test shows we are sampling from same distribution

Key figure:  plot of marginal value of information.  How this is presented depends on what form this MVOI takes; if it's constant, a histogram works, but if it's linear or anything else, a plot of this functional form is awesome.

### 6 Conclusions

Provides guidance on two key questions (1) how much data do we need to constrain a stationary PP/GPD model? and (2) what is the value of additional data?
