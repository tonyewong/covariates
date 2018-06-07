#===============================================================================
# parameter_setup.R
#
# set up parameters, parameter names, and ranges for daily peaks-over-threshold
# model
#
# Questions? Tony Wong (<anthony.e.wong@colorado.edu>)
#===============================================================================
# Copyright 2018 Tony Wong
#
# MESS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MESS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MESS.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

parnames <- c('lambda','sigma','xi')

# set up parameter bounds
bound_lower <- c(0,0,-2)
bound_upper <- c(1,800,2)

#===============================================================================
# End
#===============================================================================
