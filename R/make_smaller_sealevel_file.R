# make_smaller_sealevel_file.R
#
# Read in the data from Wong and Keller 2017, gamma priors, sea level
# projections. Then, spit out only the GMSL projections. Otherwise, the whole
# results file is going to be big and we do not want to tote that around.
#===============================================================================

# read
ncdata <- nc_open("../data/BRICK_physical_fd-gamma_08May2017.nc")
gmsl_rcp26 <- t(ncvar_get(ncdata, 'GlobalSeaLevel_RCP26'))
gmsl_rcp45 <- t(ncvar_get(ncdata, 'GlobalSeaLevel_RCP45'))
gmsl_rcp60 <- t(ncvar_get(ncdata, 'GlobalSeaLevel_RCP60'))
gmsl_rcp85 <- t(ncvar_get(ncdata, 'GlobalSeaLevel_RCP85'))
time_proj <- t(ncvar_get(ncdata, 'time_proj'))
nc_close(ncdata)

# write
dim.tproj <- ncdim_def('time_proj', 'years', as.double(time_proj))
dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:nrow(gmsl_rcp26)), unlim=TRUE)
gsl.rcp26 <- ncvar_def('GlobalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                longname = 'Global sea level (RCP26)')
gsl.rcp45 <- ncvar_def('GlobalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                longname = 'Global sea level (RCP45)')
gsl.rcp60 <- ncvar_def('GlobalSeaLevel_RCP60', 'meters', list(dim.tproj, dim.ensemble), -999,
                longname = 'Global sea level (RCP60)')
gsl.rcp85 <- ncvar_def('GlobalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                longname = 'Global sea level (RCP85)')
outnc <- nc_create('../data/BRICK_GMSL_WongKeller2017.nc',
                   list(gsl.rcp26, gsl.rcp45, gsl.rcp60, gsl.rcp85), force_v4 = TRUE)
ncvar_put(outnc, gsl.rcp26, t(gmsl_rcp26))
ncvar_put(outnc, gsl.rcp45, t(gmsl_rcp45))
ncvar_put(outnc, gsl.rcp60, t(gmsl_rcp60))
ncvar_put(outnc, gsl.rcp85, t(gmsl_rcp85))
nc_close(outnc)

#===============================================================================
# End
#===============================================================================
