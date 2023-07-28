# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# MAIN execution of NABB analysis
#####################################


# ------------ IMPORT FROM OTHER FILES ---------------------

from NABB_PCF import*
import NABB_discrete
import NABB_continuous

# ----------------------------------------------------------


# ----------------- IMPORT MODULES -------------------------

import os
import shutil
import time

# ----------------------------------------------------------

# calculate the empirical disconnected SAR
if(not os.path.exists("NABB_dSAR_emp.csv")):
	SAC_emp()

# create the 'PCF' directory if it does not exist
if(not os.path.exists("PCF")):
	os.mkdir("PCF")

# create the output directory and remove if it existed already
if(os.path.exists("PCF/PCF_out")):
	shutil.rmtree("PCF/PCF_out")
os.mkdir("PCF/PCF_out")

############################
# PARAMETERS
n_sub = 1000# number of samples in the subset
n_rep = 20	# number of replicates
n_cpu = 4	# number of CPU cores
############################

if(not os.path.exists("PCF/Output_"+str(n_sub))):
	os.mkdir("PCF/Output_"+str(n_sub))

# calculate the empirical PCF for a number of replicates 
# write the results to the csv file
for i in range(0,n_rep):
	print('\nreplicate ' + str(i+1) + ' of ' + str(n_rep))
	replicateRun(n_sub,i,n_cpu)

# do the summary calculations
NABB_discrete.calcReplicates(n_sub)
NABB_continuous.calcReplicates(n_sub)

# output the n_sub parameter for the plots
with open('PCF/PARAM_file.txt', 'w') as f:
	f.write(str(n_sub))
	f.close()