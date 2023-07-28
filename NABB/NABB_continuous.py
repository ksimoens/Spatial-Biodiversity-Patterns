# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Calculations for the CONTINUOUS Analytical Model and the NABB data set
#####################################


# ----------------- IMPORT MODULES -------------------------

import numpy as np
from scipy import special
from scipy import stats
import pandas as pd

# ----------------------------------------------------------


# calculate the spatial α parameter at a scale R
# requires: R = radius of the scale (km)
#			rho = ρ parameter of the PCF (km)
# 			lambd = λ parameter of the PCF (km)
# https://doi.org/10.1111/2041-210X.12319
def alpha_function(R,rho,lambd):

	# distances relative to the correlation length
	x = R/lambd
	a = np.pi * pow((R/rho),2) / ( 1. - 2./x * special.k1(x) * special.i1(x) / (special.i0(x)*special.k1(x) + special.i1(x)*special.k0(x)) )

	return(a)


# calculate the spatial β parameter at a scale R
# requires: R = radius of the scale (km)
#			rho = ρ parameter of the PCF (km)
# 			lambd = λ parameter of the PCF (km)
# 			n = the mean density of individuals per species in the total area (km^-2)
# https://doi.org/10.1111/2041-210X.12319
def beta_function(R,rho,lambd,n):
	x = R/lambd
	b = n * pow(rho,2) * ( 1. - 2./x * special.k1(x) * special.i1(x) / (special.i0(x)*special.k1(x) + special.i1(x)*special.k0(x)) )
	return(b)


# the limit of of alpha_function() for x -> 0 to avoid computing errors when R << λ 
def alpha_limit(R,rho,lambd):

	x = R/lambd
	a = np.pi * pow((R/rho),2) / (-1./8.*pow(x,2)*4.*np.log(x))

	return(a)


# the limit of of alpha_function() for x -> 0 to avoid computing errors when R << λ 
def beta_limit(R,rho,lambd,n):

	x = R/lambd
	b = n * pow(rho,2) * (-1./8.*pow(x,2)*4.*np.log(x))
	return(b)


# the integral in the calculations for the number of species in an area
# the integral goes from 1 to infinity
# requires: alpha = α parameter at the scale of the area
#			beta = β parameter at the scale of the area
# https://doi.org/10.1111/2041-210X.12319
def integral_R(alpha,beta):

	I_num = (1.-stats.gamma.cdf(x=1,a=alpha,scale=beta)) 
	return(I_num)


# the integral in the calculations for the number of species in an area
# the integral goes from low to up and returns the number of species with abundances between low and up
# requires: low = the lower limit for the abundance
#			up = the upper limit for the abundance
#			alpha = α parameter at the scale of the area
#			beta = β parameter at the scale of the area
# https://doi.org/10.1111/2041-210X.12319
def integral_R_bound(low,up,alpha,beta):

	I = stats.gamma.cdf(x=up,a=alpha,scale=beta) - stats.gamma.cdf(x=low,a=alpha,scale=beta)

	return(I)


# calculate the number of species in an area with radius R
# given the number of species Ss in an area with radius Rs
# can be used to upscale and downscale the number of species
# requires: Is = the normalised integral at the scale of the known area Rs
#			rho = ρ parameter of the PCF (km)
# 			lambd = λ parameter of the PCF (km)
# 			n = the mean density of individuals per species in the contiguous USA (km^-2)
#			R = the radius of the area to be estimated (km)
# https://doi.org/10.1111/2041-210X.12319
def S(Is,rho,lambd,n,R):

	# get the α and β parameters at the scale to be estimated Rs
	alpha_R = alpha_function(R,rho,lambd)
	beta_R = beta_function(R,rho,lambd,n)

	# calculate the number of species at scale R
	S = integral_R(alpha_R,beta_R) / Is

	return(S)


# calculate the spatial SAD at scale R in Preston plot and output dataframe
# requires:	R0 = the radius of a circle with the same surface area as the contiguous USA (km)
#			S0 = the (estimated) total number of indiviuals in the NABB data set
#			rho = ρ parameter of the PCF (km)
# 			lambd = λ parameter of the PCF (km)
# 			n = the mean density of individuals per species in the contiguous USA (km^-2)
#			R = the radius of the scale at which to estimate the SAD (km)
#			Nmax = the maximum number of individuals for which to calculate the probability
# https://doi.org/10.1111/2041-210X.12319
def sSAD(R0,S0,rho,lambd,n,R,N_max):

	# calculate the α and β parameter for the total area
	alpha_0 = alpha_function(R0,rho,lambd)
	beta_0 = beta_function(R0,rho,lambd,n)
	# calculate the normalised integral for the total area
	I_0 = integral_R(alpha_0,beta_0) / S0

	# create the containers for the number of individuals and the number of species
	S_list = np.asarray([])
	N_list = np.asarray([])

	# calculate the α and β parameter for the scale R to estimate
	alpha_R = alpha_function(R,rho,lambd)
	beta_R = beta_function(R,rho,lambd,n)

	N = 0
	# for all logarithmic abundances below N_max
	while(N < N_max+1):
		# calculate the bounded integral at scale R
		I_num = integral_R_bound(low=pow(2,N),up=pow(2,N+1),alpha=alpha_R,beta=beta_R) 
		# add to the container
		S_list = np.append(S_list,I_num / I_0)
		N_list = np.append(N_list,N)

		# continuous distribution
		N += 0.1

	# create the dataframe and write out
	darray = np.concatenate((N_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['N','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	return(df_out)


# calculate the theoretical upscaled disconnected SAR
# requires: Ss = the number of species at known disconnected area of combined samples
#			Rs = radius of the combined disconnected area (km)
#			rho = ρ parameter of the PCF (km)
# 			lambd = λ parameter of the PCF (km)
# 			n = the mean density of individuals per species in the contiguous USA (km^-2)
#			R0 = the radius of the contiguous USA (km)
def dSAR(Ss,Rs,rho,lambd,n,R0):

	# calculate the α and β parameters at the level of the known disconnected area Rs
	alpha_s = alpha_function(Rs,rho,lambd)
	beta_s = beta_function(Rs,rho,lambd,n)
	# calculate the normalised integral at the known scale Rs
	I_s = integral_R(alpha_s,beta_s) / Ss

	# create the containers for the number of samples and the number of species
	S_list = np.array([])
	R_list = np.array([])
	# surface area of a single sample (km^2)
	Ai = 50.*np.pi*0.4*0.4
	# get the total number of samples
	Nsample = 2287
	
	# for each number of samples from 100
	for k in range(100,Nsample+1):
		# get the radius for the disconnected upscaled area
		R = np.sqrt(k*Ai/np.pi)
		
		# calculate the number of species in the upscaled area
		S_R = S(I_s,rho,lambd,n,R)

		# add to the container
		S_list = np.append(S_list,S_R)
		R_list = np.append(R_list,R)

		k += 1

	darray = np.concatenate((R_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['R','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	return(df_out)


# do the calculations for all replicates
# requires: n_sub = the number of samples in the subset
def calcReplicates(n_sub):

	# global parameters
	Ai = 50.*np.pi*0.4*0.4		# surface area of a single sample (km^2)
	A0 = 7663941.7	 			# surface area of the contiguous USA (km^2)
	R0 = np.sqrt(A0/np.pi) 		# radius for the total area (km)

	# subset parameters 
	As = n_sub*Ai 				# combined surface area of the subset (km^2)
	Rs = np.sqrt(As/np.pi)		# corresponding radius (km)

	# data set parameters
	At = 2287.*Ai 				# combined surface area of all samples (km^2)
	Rt = np.sqrt(At/np.pi)		# corresponding radius (km)

	# get the replicate list with the empirical data (created in NABB_PCF.py)
	rep_list = pd.read_csv('PCF/Output_' +str(n_sub)+ '/rep_list_' + str(n_sub) + '.csv')

	# create new columns in the replicate list
	rep_list['S0'] = np.nan 		# total number of species
	rep_list['np'] = np.nan 		# mean density of individuals per species from subset (m^-2)
	rep_list['n'] = np.nan 			# mean density of individuals per species from total estimate (m^-2)
	rep_list['St'] = np.nan 		# number of species in all samples
	rep_list['alpha0'] = np.nan		# α parameter for the total area
	rep_list['beta0'] = np.nan 		# β parameter for the total area

	# create containers for the SAD and SAR
	SAD0_list = pd.DataFrame(columns=['N', 'S'])
	tSAD_list = pd.DataFrame(columns=['N', 'S'])
	dSAR_list = pd.DataFrame(columns=['R', 'S'])

	# for each replicate
	for i in range(0,len(rep_list)):

		# calculate the density from the subset
		n = rep_list['Nsum'][i]/rep_list['S'][i]/As
		rep_list['np'][i] = n

		# calculate the total number of species in the contiguous USA and all samples
		alpha_s = alpha_function(Rs,rep_list['rho'][i],rep_list['lambd'][i])
		beta_s = beta_function(Rs,rep_list['rho'][i],rep_list['lambd'][i],n)
		I_s = integral_R(alpha_s,beta_s) / rep_list['S'][i]
		rep_list['S0'][i] = S(I_s,rep_list['rho'][i],rep_list['lambd'][i],n,R0)
		rep_list['St'][i] = S(I_s,rep_list['rho'][i],rep_list['lambd'][i],n,Rt)
		
		# calculate the new estimate for the mean density per species
		n_down = n*rep_list['S'][i]/rep_list['S0'][i]
		rep_list['n'][i] = n_down

		# calculate the α and β parameters for the contiguous USA
		rep_list['alpha0'][i] = alpha_function(R0,rep_list['rho'][i],rep_list['lambd'][i])
		rep_list['beta0'][i] = beta_function(R0,rep_list['rho'][i],rep_list['lambd'][i],n_down)

		# create the SAD for a connected subarea of the same size as all samples combined
		tSAD_list = pd.concat( [tSAD_list , sSAD(R0,rep_list['S0'][i],rep_list['rho'][i],rep_list['lambd'][i],n_down,Rt,18)] )
		# create the SAD for the contiguous USA and add to the container
		SAD0_list = pd.concat( [SAD0_list , sSAD(R0,rep_list['S0'][i],rep_list['rho'][i],rep_list['lambd'][i],n_down,R0,24)] )
		# create the disconnected SAR and add to the container
		dSAR_list = pd.concat( [dSAR_list , dSAR(rep_list['S'][i],Rs,rep_list['rho'][i],rep_list['lambd'][i],n,R0) ] )


	# write out the calculated quantities
	rep_list.to_csv('PCF/Output_' +str(n_sub)+ '/rep_list_' + str(n_sub) + '_continuous.csv')

	# summarise for each Preston class and calculate mean and standard deviation
	tSAD_sum = tSAD_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	tSAD_sum.columns = ['N','Smean','Ssd']
	tSAD_sum.to_csv('PCF/Output_' +str(n_sub)+ '/tSAD_' + str(n_sub) + '_continuous.csv')
	
	SAD0_sum = SAD0_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	SAD0_sum.columns = ['N','Smean','Ssd']
	SAD0_sum.to_csv('PCF/Output_' +str(n_sub) +'/SAD0_' + str(n_sub) + '_continuous.csv')

	# summarise for each radius R and calculate mean and standard deviation
	dSAR_sum = dSAR_list.groupby(['R'],as_index=False).agg({'S':['mean','std']})
	dSAR_sum.columns = ['R','Smean','Ssd']
	dSAR_sum.to_csv('PCF/Output_' +str(n_sub) +'/dSAR_' + str(n_sub) + '_continuous.csv')