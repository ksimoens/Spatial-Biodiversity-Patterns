# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Calculations for the DISCRETE Analytical Model and the NABB data set
#####################################


# ----------------- IMPORT MODULES -------------------------

import numpy as np
from scipy import special
from scipy import stats
import pandas as pd

# ----------------------------------------------------------


# calculate the spatial mean number of individuals per species at a scale R
# requires: R = radius of the scale (km)
#			n = spatial mean number of individuals per species per area in the contiguous USA (km^-2)
# https://doi.org/10.1111/2041-210X.12319
def mean(R,n):

	return(n*np.pi*R*R)


# calculate the spatial variance of the number of individuals per species at scale R
# requires: R = radius of the scale (km)
#			n = spatial mean number of individuals per species per area in the contiguous USA (km^-2)
#			rho = ρ parameter of the PCF (km)
# 			lambd = λ parameter of the PCF (km)
# https://doi.org/10.1111/2041-210X.12319
def sigma(R,n,rho,lambd):

	# get the pre-factor
	factor = n*n*rho*rho*np.pi*R*R
	# distances relative to correlation length
	x = R / lambd
	# get Bessel functions
	bessels = ( 1. - 2./x * special.k1(x) * special.i1(x) / (special.i0(x)*special.k1(x) + special.i1(x)*special.k0(x)) )

	return(factor*bessels)


# calculate the limit of sigma() for x -> 0 to avoid computational errors when λ is large compared to R
def sigma_limit(R,n,rho,lambd):

	factor = n*n*rho*rho*np.pi*R*R
	x = R / lambd
	bessels = (-1./8.*pow(x,2)*4.*np.log(x))
	return(factor*bessels)


# calculate the spatial ξ parameter at a scale defined by Nmean and Nsigma
# requires: Nmean = the mean number of individuals per species at scale R
#			Nsigma = the variance of the number of individuals per species at scale R
def xi_function(Nmean,Nsigma):

	return(1.-Nmean/Nsigma)


# calculate the spatial r parameter at a scale defined by Nmean and Nsigma
# requires: Nmean = the mean number of individuals per species at scale R
#			Nsigma = the variance of the number of individuals per species at scale R
def r_function(Nmean,Nsigma):

	return(Nmean/(Nsigma/Nmean - 1.))


# calculate the spatial ξp parameter for a fraction of the contiguous USA p
# requires: xi = the ξ parameter at the scale of the contiguous USA
#			p = the sampled fraction
# https://doi.org/10.1126/sciadv.1701438
def xip(xi,p):

	return(p*xi/(1.-xi*(1.-p)))


# calculate the number of species at a larger scale -> upscaling
# requires: xi = ξ parameter at the larger (connected) scale
#			r = r parameter
#			p = the sampled fraction = smaller scale
#			Sp = the number of species found in the smaller scale
# https://doi.org/10.1126/sciadv.1701438
def Sup(xi,r,p,Sp):

	return( (1.-pow(1.-xi,r)) / (1.-pow(1.-xip(xi,p),r) ) * Sp )


# inverse of the Sup() method
# calculate the number of species at a smaller scale -> downscaling
# requires:	xi = ξ parameter at the larger (connected) scale
#			r = r parameter
#			p = the sampled fraction = smaller scale
#			S = the number of species found in the larger scale
# https://doi.org/10.1126/sciadv.1701438
def Sdown(xi,r,p,S):

	return( (1.-pow(1.-xip(xi,p),r)) / (1.-pow(1.-xi,r))*S )


# calculate the number of species in the contiguous USA
#	starting from a randomly sampled subset of samples
# requires: S = the number of species in the subset of samples
#			R0 = the radius of a circle with the same surface area as the contiguous USA (km)
#			rho = ρ parameter of the PCF (km)
#			lambd = λ parameter of the PCF (km)
#			N = the total number of individuals found in the subset of samples
#			A = total surface area of the subset (km^2)
#			p = fraction of total surface area sampled
def S0(S,R0,rho,lambd,N,A,p):

	# initialise the S0 prediction
	S0 = 0
	# first estimate of the mean density
	n = N/A/S
	# iteratively solve for S0
	i = 0
	while(i < 100):
		# calculate the ξ and r parameters
		xi = xi_function(mean(R0,n),sigma(R0,n,rho,lambd))
		r = r_function(mean(R0,n),sigma(R0,n,rho,lambd))
		# update the total number of species
		S0 = Sup(xi,r,p,S)
		# update the mean density
		n = N/A/S0
		i += 1

	return(S0)


# get the normalisation factor for the RSA / SAD
# requires:	xi = ξ parameter at the scale of the disconnected area
#			r = r parameter
# https://doi.org/10.1126/sciadv.1701438
def c(xi,r):

	return(1./(1.-pow(1.-xi,r)))


# get the integral for the number of species with abundances between low and up
# requires:	xi = ξ parameter at the scale of the disconnected area
#			r = r parameter
#			low = lower limit of the abundance
#			up = upper limit of the abundance
# https://doi.org/10.1126/sciadv.1701438 
def SAD(xi,r,low,up):

	return( c(xi,r) *( stats.nbinom.cdf(k=up,n=r,p=1-xi) - stats.nbinom.cdf(k=low,n=r,p=1-xi) ) )


# calculate the spatial SAD at scale p in Preston plots and output dataframe for a disconnected (sub)area
# requires:	R0 = the radius of a circle with the same surface area as the contiguous USA (km)
#			S0 = the (estimated) total number of indiviuals in the contiguous USA
#			xi = the ξ parameter at the scale of the contiguous USA
#			r = the r parameter 
#			n = the mean number of individuals per species per area in the contiguous USA (km^-2)
#			p = the fraction of sampled area in the subset of samples
#			Nmax = the maximum number of individuals for which to calculate the probability
# https://doi.org/10.1126/sciadv.1701438 
def sSAD(R0,S0,xi,r,n,p,Nmax):

	# calculate the ξp parameter at the scale of the subset of samples
	# downscale it from the estimate for the contiguous USA
	xip0 = xip(xi,p)
	# calculate the number of species at scale p
	Sp = Sdown(xi,r,p,S0)

	# create containers for the number of individuals and species
	N_list = np.array([])
	S_list = np.array([])
	# for each logarithmic bin smaller than Nmax
	N = 0
	while(N < Nmax+1):
		# calculate the number of species with abundances in the bin
		Sn = SAD(xip0,r,low=pow(2,N),up=pow(2,N+1))*Sp
		# add to the container
		N_list = np.append(N_list,N)
		S_list = np.append(S_list,Sn)

		N += 1

	# create a final dataframe and return
	darray = np.concatenate((N_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['N','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	return(df_out)


# calculate the SAD at the scale of the contiguous USA
# requires: xi = the ξ parameter at the scale of the contiguous USA
#			r = the r parameter
#			S0 = the number of species in the contiguous USA
def SAD0(xi,r,S0,Nmax):

	# create containers for the number of individuals and species
	N_list = np.array([])
	S_list = np.array([])
	# for each logarithmic bin smaller than Nmax
	N = 0
	while(N < Nmax+1):
		# calculate the number of species with abundances in the bin
		Sn = SAD(xi,r,low=pow(2,N),up=pow(2,N+1))*S0
		# add to the container
		N_list = np.append(N_list,N)
		S_list = np.append(S_list,Sn)

		N += 1

	# create a final dataframe and return
	darray = np.concatenate((N_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['N','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	return(df_out)


# calculate the theoretical disconnected SAR 
def dSAR(xi,r,R0,S0):
	
	# surface area of a single sample (km^2)
	Ai = 50.*np.pi*0.4*0.4
	A0 = R0*R0*np.pi
	Nsample = 2287

	# create containers for the number of species and the radius
	R_list = np.array([])
	S_list = np.array([])
	
	# for each number of samples from 100
	for k in range(100,Nsample+1):

		# radius at the scale of the subset of samples (km)
		R = np.sqrt(k*Ai/np.pi)
		# fraction of the area
		p = k*Ai/A0
		# calculate the number of species in the subset
		Sp = Sdown(xi,r,p,S0)
		# add to container
		S_list = np.append(S_list,Sp)
		R_list = np.append(R_list,R)

		k += 1

	# create dataframe and return
	darray = np.concatenate((R_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['R','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	return(df_out)


# do the calculations for all replicates
# requires: n_sub = the number of samples in the subset
def calcReplicates(n_sub):

	# global parameters
	Ai = 50.*np.pi*0.4*0.4 		# surface area of a single sample (km^2)
	A0 = 7663941.7	 			# surface area of the contiguous USA (km^2)
	R0 = np.sqrt(A0/np.pi) 		# radius for the contiguous USA (km)

	# subset parameters 
	As = n_sub*Ai 				# combined surface area of the subset (km^2)
	Rs = np.sqrt(As/np.pi)		# corresponding radius (km)

	# data set parameters
	At = 2287.*Ai 				# surface area of all samples combined (km^2)
	Rt = np.sqrt(At/np.pi) 		# corresponding radius (km)

	# get downscaled radius
	# if n_sub = 2287; downscale to 500 samples
	if(n_sub == 2287):
		Ac = 500.*Ai
	else:
		Ac = n_sub*Ai


	# get the replicate list with the empirical data (created in NABB_PCF.py)
	rep_list = pd.read_csv('PCF/Output_' +str(n_sub)+ '/rep_list_' + str(n_sub) + '.csv')

	# create new columns in the replicate list
	rep_list['S0'] = np.nan 		# total number of species
	rep_list['np'] = np.nan 		# mean density of individuals per species from subset (km^-2)
	rep_list['n'] = np.nan 			# mean density of individuals per species from total estimate (km^-2)
	rep_list['St'] = np.nan  		# number of species in all samples
	rep_list['Nmean'] = np.nan 		# mean number of individuals per species in contiguous USA
	rep_list['Nsd'] = np.nan 		# variance in the number of individuals per species in contiguous USA
	rep_list['r'] = np.nan 			# r parameter
	rep_list['xi'] = np.nan 		# ξ parameter for the contiguous USA
	rep_list['xip'] = np.nan     	# ξp parameter for the subset

	# create containers for the SAD and SAR
	sSAD_list = pd.DataFrame(columns=['N', 'S'])	# the SAD of a disconnected subarea with n_sub samples
	tSAD_list = pd.DataFrame(columns=['N', 'S'])	# the SAD of a disconnected subarea with all samples
	SAD0_list = pd.DataFrame(columns=['N', 'S'])	# the SAD at the scale of the contiguous USA
	dSAR_list = pd.DataFrame(columns=['R', 'S'])	# the disconnected (downscaled) SAR

	# for each replicate
	for i in range(0,len(rep_list)):

		# calculate the density from the subset
		n = rep_list['Nsum'][i]/rep_list['S'][i]/As
		rep_list['np'][i] = n

		# calculate the total number of species
		rep_list['S0'][i] = S0(rep_list['S'][i],R0,rep_list['rho'][i],rep_list['lambd'][i],rep_list['Nsum'][i],As,As/A0)
		
		# calculate the new estimate for the mean density per species
		n_down = n*rep_list['S'][i]/rep_list['S0'][i]
		rep_list['n'][i] = n_down

		# calculate the new mean and variance
		rep_list['Nmean'][i] = mean(R0,n_down)
		rep_list['Nsd'][i] = sigma(R0,n_down,rep_list['rho'][i],rep_list['lambd'][i])

		# calculate the negative binomial parameters
		rep_list['r'][i] = r_function(rep_list['Nmean'][i],rep_list['Nsd'][i])
		rep_list['xi'][i] = xi_function(rep_list['Nmean'][i],rep_list['Nsd'][i])
		rep_list['xip'][i] = xip(rep_list['xi'][i], As/A0)

		# calculate the downscaled number of species in all samples
		rep_list['St'][i] = Sdown(rep_list['xi'][i],rep_list['r'][i],At/A0,rep_list['S0'][i])

		# create the sSAD and add to the container
		sSAD_list = pd.concat( [sSAD_list , sSAD(R0,rep_list['S0'][i],rep_list['xi'][i],rep_list['r'][i],n_down,Ac/A0,np.log2(rep_list['Nsum'][i]))] )
		# create the SAD for the contiguous USA and add to the container
		SAD0_list = pd.concat( [SAD0_list , SAD0(rep_list['xi'][i],rep_list['r'][i],rep_list['S0'][i],24)] )
		# create the tSAD and add to the container
		tSAD_list = pd.concat( [tSAD_list , sSAD(R0,rep_list['S0'][i],rep_list['xi'][i],rep_list['r'][i],n_down,At/A0,18)] )
		# create the disconnected SAR and add to the container
		dSAR_list = pd.concat( [dSAR_list , dSAR(rep_list['xi'][i],rep_list['r'][i],R0,rep_list['S0'][i]) ] )
		

	# write out the calculated quantities
	rep_list.to_csv('PCF/Output_' +str(n_sub)+ '/rep_list_' + str(n_sub) + '_discrete.csv')

	# summarise for each Preston class and calculate mean and standard deviation
	sSAD_sum = sSAD_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	sSAD_sum.columns = ['N','Smean','Ssd']
	sSAD_sum.to_csv('PCF/Output_' +str(n_sub)+ '/sSAD_' + str(n_sub) + '_discrete.csv')

	SAD0_sum = SAD0_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	SAD0_sum.columns = ['N','Smean','Ssd']
	SAD0_sum.to_csv('PCF/Output_' +str(n_sub)+ '/SAD0_' + str(n_sub) + '_discrete.csv')

	tSAD_sum = tSAD_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	tSAD_sum.columns = ['N','Smean','Ssd']
	tSAD_sum.to_csv('PCF/Output_' +str(n_sub)+ '/tSAD_' + str(n_sub) + '_discrete.csv')

	# summarise for each radius R and calculate mean and standard deviation
	dSAR_sum = dSAR_list.groupby(['R'],as_index=False).agg({'S':['mean','std']})
	dSAR_sum.columns = ['R','Smean','Ssd']
	dSAR_sum.to_csv('PCF/Output_' +str(n_sub)+ '/dSAR_' + str(n_sub) + '_discrete.csv')