# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Run the code for the CPR data
#####################################

# ----------------- IMPORT MODULES -------------------------

from CPR_spatial import*
from CPR_MF import*

# ----------------------------------------------------------

# calculate the empirical disconnected SAR
SAC_emp()

# calculate the negative log-likelihood for a combination of spatial parameters
outputLL()

# fit the Mean Field parameters
calculations()

# calculate the Mean Field disconnected SAR
dSAR_MF()