---
# Extreme value models - covariates

Questions? Tony Wong <anthony.e.wong@colorado.edu>

---
## Pipeline

### 1 Processing

Same as in Wong et al 2018. Exact same data_calib$norfolk object as used there, with full-data form. To run, use the `process_data.R` script in the `R` directory.

In a nutshell, I take raw hourly sea level data from the University of Hawaii Sea Level Center (UHSLC) and (1) detrend using a moving window 1-year mean, (2) take daily maxima, (3) use the 99th percentile as the threshold for extreme events and take only exceedances of this threshold as "extremes", and (4) decluster using a 3-day time-scale.

### 2 Prior distributions

Priors also same as in Wong et al 2018. Using normal and gamma distributions, fit to the distribution of about 30 long tide gauge records from UHSLC, plus a couple other long records. This is done using the `fit_priors.R` routine.

### 3 Calibration

I use a robust adaptive Markov chain Monte Carlo approach for parameter calibration. The paper associated with this work has some references that are good additional reading. The essence of this approach is to create a Markov chain of Poisson process/generalized Pareto distribution model parameters, whose stationary distribution is the posterior distribution of these parameters, given the processed set of peaks-over-thresholds data.

The robust adaptive approach is used because it is nice to be able to use a target acceptance rate (23% for a large number of parameters, or 44% for a single parameter, and somewhere in between for intermediate numbers of parameters). The target acceptance rate is achieved by adapting the proposal covariance matrix (where multivariate normal proposals are used), based on previous samples. The rate of this adaptation decreases as sampling continues, so that once the stationary distribution is reached, the tuned covariance matrix is not screwed up. Gelman and Rubin diagnostics and visual inspection are used to evaluate convergence of the Markov chains to their stationary distribution.

This is all done using the `calibration.R` routine.

### 4 Bayesian model averaging

First, do for each covariate individually to see how the BMA weights (model choice) is affected by which covariate you use. Is a particular model structured favored across the board?

Second, do using all 16 candidate models (four covariates, times four model structures) to examine which covariates/structures are favored.

This is done using the `bridge_sample.R` and then `compute_weights.R` routines.

### 5 Analysis and plots

Run `analysis.R` and `plots.R`.

If any results or analysis files are missing, or if you are trying to replicate or extend on this work but getting stuck, do not hesitate to send me an email. I am always happy to chat!

-Tony
