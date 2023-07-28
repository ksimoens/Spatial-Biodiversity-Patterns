# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Calculations for the DISCRETE Analytical Model and BCI data set
#####################################


# ----------------- IMPORT MODULES -------------------------

import numpy as np
from scipy import special
from scipy import stats
import pandas as pd

# ----------------------------------------------------------


# calculate the spatial mean number of individuals per species at a scale R
# requires: R = radius of the scale (m)
#			n = spatial mean number of individuals per species per area in the total area (m^-2)
# https://doi.org/10.1111/2041-210X.12319
def mean(R,n):

	return(n*np.pi*R*R)


# calculate the spatial variance of the number of individuals per species at scale R
# requires: R = radius of the scale (m)
#			n = spatial mean number of individuals per species per area in the total area (m^-2)
#			rho = ρ parameter of the PCF (m)
# 			lambd = λ parameter of the PCF (m)
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


# calculate the spatial ξp parameter for a fraction of the total area p
# requires: xi = the ξ parameter at the scale of the total area
#			p = the sampled fraction
# https://doi.org/10.1126/sciadv.1701438
def xip(xi,p):

	return(p*xi/(1.-xi*(1.-p)))


# calculate the number of species at a larger scale -> upscaling
# requires: xi = ξ parameter at the larger (continuous) scale
#			r = r parameter
#			p = the sampled fraction = smaller scale
#			Sp = the number of species found in the smaller scale
# https://doi.org/10.1126/sciadv.1701438
def Sup(xi,r,p,Sp):

	return( (1.-pow(1.-xi,r)) / (1.-pow(1.-xip(xi,p),r) ) * Sp )


# inverse of the Sup() method
# calculate the number of species at a smaller scale -> downscaling
# requires:	xi = ξ parameter at the larger (continuous) scale
#			r = r parameter
#			p = the sampled fraction = smaller scale
#			S = the number of species found in the larger scale
# https://doi.org/10.1126/sciadv.1701438
def Sdown(xi,r,p,S):

	return( (1.-pow(1.-xip(xi,p),r)) / (1.-pow(1.-xi,r))*S )


# calculate the number of species in the entire BCI area
#	starting from a randomly sampled subset of samples
# requires: S = the number of species in the subset of samples
#			R0 = the radius of a circle with the same surface area as the total BCI surface area (m)
#			rho = ρ parameter of the PCF (m)
#			lambd = λ parameter of the PCF (m)
#			N = the total number of individuals found in the subset of samples
#			A = total surface area of the subset (m^2)
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
		xi = xi_function(mean(R0,n),sigma_limit(R0,n,rho,lambd))
		r = r_function(mean(R0,n),sigma_limit(R0,n,rho,lambd))
		# update the total number of species
		S0 = Sup(xi,r,p,S)
		# update the mean density
		n = N/A/S0
		i += 1

	return(S0)


# get the normalisation factor for the RSA / SAD
# requires:	xi = ξ parameter at the scale of the total area
#			r = r parameter
# https://doi.org/10.1126/sciadv.1701438
def c(xi,r):

	return(1./(1.-pow(1.-xi,r)))


# get the integral for the number of species with abundances between low and up
# requires:	xi = ξ parameter at the scale of the total area
#			r = r parameter
#			low = lower limit of the abundance
#			up = upper limit of the abundance
# https://doi.org/10.1126/sciadv.1701438 
def SAD(xi,r,low,up):

	return( c(xi,r) *( stats.nbinom.cdf(k=up,n=r,p=1-xi) - stats.nbinom.cdf(k=low,n=r,p=1-xi) ) )


# get the integral for the number of species in a connected area with abundances between low and up
# requires:	xi = ξ parameter at the scale of the connected (sub)area
#			r = r parameter
#			low = lower limit of the abundance
#			up = upper limit of the abundance
def connSAD(xi,r,low,up):

	return(stats.nbinom.cdf(k=up,n=r,p=1.-xi) - stats.nbinom.cdf(k=low,n=r,p=1.-xi))



# calculate the spatial SAD at scale p in Preston plots and output dataframe for a disconnected (sub)area
# requires:	R0 = the radius of a circle with the same surface area as the total BCI surface area (m)
#			S0 = the (estimated) total number of indiviuals in the BCI data set
#			xi = the ξ parameter at the scale of the total area
#			r = the r parameter 
#			n = the mean number of individuals per species per area in the total area (m^-2)
#			p = the fraction of sampled area in the subset of samples
#			Nmax = the maximum number of individuals for which to calculate the probability
# https://doi.org/10.1126/sciadv.1701438 
def sSAD(R0,S0,xi,r,n,p,Nmax):

	# calculate the ξp parameter at the scale of the subset of samples
	# downscale it from the estimate for the total area
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


# summation of the SAD for the calculations of the number of species in a connected subarea
# requires:	xi = the ξ parameter at the scale of the connected subarea
# 			r = the r parameter
def summation(xi,r):

	return(1-stats.nbinom.cdf(k=1,n=r,p=1-xi))


# caluclate the SAD for a connected subarea
# requires: xi0 = the ξ parameter at the scale of the total area
#			r0 = the r parameter at the scale of the total area
#			S0 = the number of species in the total area
#			rho = the ρ parameter of the PCF (m)
#			lambd = the λ parameter of the PCF (m)
# 			n = the mean density of individuals per species for the total area (m^-2)
# 			R = the radius of the connected (sub)area (m)
#			Nmax = the maximum number of individuals for which to calculate the probability
def connsSAD(xi0,r0,S0,rho,lambd,n,R,Nmax):

	# calculate the summation for the entire area	
	sum0 = summation(xi0,r0)/S0

	# calculate the mean and variance at scale R
	Nmean_R = mean(R,n)
	Nsigma_R = sigma_limit(R,n,rho,lambd)
	
	# calculate the r,ξ parameters at scale R
	r_R = r_function(Nmean_R,Nsigma_R)
	xi_R = xi_function(Nmean_R,Nsigma_R)

	# create containers for the number of individuals and species
	N_list = np.array([])
	S_list = np.array([])
	# for each logarithmic bin smaller than Nmax
	N = 0
	while(N < Nmax+1):
		# calculate the number of species with abundances in the bin
		Sn = connSAD(xi_R,r_R,low=pow(2,N),up=pow(2,N+1))/sum0
		# add to the container
		N_list = np.append(N_list,N)
		S_list = np.append(S_list,Sn)

		N += 1

	# create a final dataframe and return
	darray = np.concatenate((N_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['N','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	return(df_out)


# calculate the SAD at the scale of the total area
# requires: xi = the ξ parameter at the scale of the total area
#			r = the r parameter
#			S0 = the number of species in the total area
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
	
	# surface area of a single sample (m^2)
	Ai = 25*25
	A0 = R0*R0*np.pi
	Nsample = int(A0/Ai)

	# create containers for the number of species and the radius
	R_list = np.array([])
	S_list = np.array([])
	
	# for each number of samples
	for k in range(20,Nsample+1):

		# radius at the scale of the subset of samples (m)
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


# calculate the theoretical connected SAR
# requires:	xi0 = the ξ parameter at the scale of the total area
#			r0 = the r parameter at the scale of the local area
#			rho = the ρ parameter of the PCF
#			lambd = the λ parameter of the PCF
# 			n = the mean density of individuals per species at the total area
#			R0 = the radius of the total area
#			S0 = the number of species in the total area
def cSAR(xi0,r0,rho,lambd,n,R0,S0):

	# calculate the summation for the entire area
	sum0 = summation(xi0,r0)/S0

	# create containers for the number of species and the radius
	R_list = np.array([])
	S_list = np.array([])

	# radius for a single sample
	Ai = 25*25 										# m^2
	Rmin = np.floor(np.sqrt(Ai/np.pi)).astype(int)	# m
	R0_int = np.ceil(R0).astype(int)				# m

	for R in range(Rmin,R0_int+1):

		# calculate the mean and variance at scale R
		Nmean_R = mean(R,n)
		Nsigma_R = sigma_limit(R,n,rho,lambd)
	
		# calculate the r,ξ parameters at scale R
		r_R = r_function(Nmean_R,Nsigma_R)
		xi_R = xi_function(Nmean_R,Nsigma_R)

		# calculate the number of species at scale R
		S_R = summation(xi_R,r_R) / sum0

		# add to the containers
		R_list = np.append(R_list,R)
		S_list = np.append(S_list,S_R)

	# create the output and write out
	darray = np.concatenate((R_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['R','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	return(df_out)


# do the calculations for all replicates
# requires: n_sub = the number of samples in the subset
def calcReplicates(n_sub):

	# global parameters
	Ai = 25.*25. 				# surface area of a single sample (m^2)
	A0 = 1000.*500. 			# surface area of the entire BCI data set (m^2)
	R0 = np.sqrt(A0/np.pi) 		# radius for the total area (m)

	# subset parameters 
	As = n_sub*Ai 				# combined surface area of the subset (m^2)
	Rs = np.sqrt(As/np.pi)		# corresponding radius (m)

	# R at which to calculate the SAD connected area 
	# if all samples are included, calculate SAD for 100 samples
	if(n_sub == 800):
		Ac = 100.*Ai
	else:
		Ac = As 
	Rc = np.sqrt(Ac/np.pi)

	# get the replicate list with the empirical data (created in BCI_PCF.py)
	rep_list = pd.read_csv('PCF/Output_' +str(n_sub)+ '/rep_list_' + str(n_sub) + '.csv')

	# create new columns in the replicate list
	rep_list['S0'] = np.nan 		# total number of species
	rep_list['np'] = np.nan 		# mean density of individuals per species from subset (m^-2)
	rep_list['n'] = np.nan 			# mean density of individuals per species from total estimate (m^-2)
	rep_list['Ss'] = np.nan 		# number of species in a connected subarea of size As
	rep_list['Nmean'] = np.nan 		# mean number of individuals per species in total area
	rep_list['Nsd'] = np.nan 		# variance in the number of individuals per species in total area
	rep_list['r'] = np.nan 			# r parameter
	rep_list['xi'] = np.nan 		# ξ parameter for the total area
	rep_list['xip'] = np.nan     	# ξp parameter for the subset

	# create containers for the SAD and SAR
	sSAD_list = pd.DataFrame(columns=['N', 'S'])	# the SAD of a disconnected subarea with n_sub samples
	SAD0_list = pd.DataFrame(columns=['N', 'S'])	# the SAD at the scale of the total area
	cSAD_list = pd.DataFrame(columns=['N', 'S'])	# the SAD of a connected subarea with the same surface area as the samples
	dSAR_list = pd.DataFrame(columns=['R', 'S'])	# the disconnected (downscaled) SAR
	cSAR_list = pd.DataFrame(columns=['R', 'S'])	# the connected (downscaled) SAR

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
		rep_list['Nsd'][i] = sigma_limit(R0,n_down,rep_list['rho'][i],rep_list['lambd'][i])

		# calculate the negative binomial parameters
		rep_list['r'][i] = r_function(rep_list['Nmean'][i],rep_list['Nsd'][i])
		rep_list['xi'][i] = xi_function(rep_list['Nmean'][i],rep_list['Nsd'][i])
		rep_list['xip'][i] = xip(rep_list['xi'][i], As/A0)

		sum0 = summation(rep_list['xi'][i],rep_list['r'][i])/rep_list['S0'][i]
		Nmean_R = mean(Rs,n_down)
		Nsigma_R = sigma_limit(Rs,n_down,rep_list['rho'][i],rep_list['lambd'][i])
		r_R = r_function(Nmean_R,Nsigma_R)
		xi_R = xi_function(Nmean_R,Nsigma_R)
		rep_list['Ss'][i] = summation(xi_R,r_R) / sum0

		# create the sSAD and add to the container
		sSAD_list = pd.concat( [sSAD_list , sSAD(R0,rep_list['S0'][i],rep_list['xi'][i],rep_list['r'][i],n_down,Ac/A0,np.log2(rep_list['Nsum'][i]))] )
		# create the SAD for the entire area and add to the container
		SAD0_list = pd.concat( [SAD0_list , SAD0(rep_list['xi'][i],rep_list['r'][i],rep_list['S0'][i],18)] )
		# create the SAD for a connected subarea of equal size as the combined samples
		cSAD_list = pd.concat( [cSAD_list , connsSAD(rep_list['xi'][i],rep_list['r'][i],rep_list['S0'][i],rep_list['rho'][i],rep_list['lambd'][i],rep_list['n'][i],Rc,np.log2(rep_list['Nsum'][i]))] )
		# create the disconnected SAR and add to the container
		dSAR_list = pd.concat( [dSAR_list , dSAR(rep_list['xi'][i],rep_list['r'][i],R0,rep_list['S0'][i]) ] )
		# create the connected SAR and add to the container
		cSAR_list = pd.concat( [cSAR_list , cSAR(rep_list['xi'][i],rep_list['r'][i],rep_list['rho'][i],rep_list['lambd'][i],rep_list['n'][i],R0,rep_list['S0'][i])] )


	# write out the calculated quantities
	rep_list.to_csv('PCF/Output_' +str(n_sub)+ '/rep_list_' + str(n_sub) + '_discrete.csv')

	# summarise for each Preston class and calculate mean and standard deviation
	sSAD_sum = sSAD_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	sSAD_sum.columns = ['N','Smean','Ssd']
	sSAD_sum.to_csv('PCF/Output_' +str(n_sub)+ '/sSAD_' + str(n_sub) + '_discrete.csv')

	SAD0_sum = SAD0_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	SAD0_sum.columns = ['N','Smean','Ssd']
	SAD0_sum.to_csv('PCF/Output_' +str(n_sub)+ '/SAD0_' + str(n_sub) + '_discrete.csv')

	cSAD_sum = cSAD_list.groupby(['N'],as_index=False).agg({'S':['mean','std']})
	cSAD_sum.columns = ['N','Smean','Ssd']
	cSAD_sum.to_csv('PCF/Output_' +str(n_sub)+ '/cSAD_' + str(n_sub) + '_discrete.csv')

	# summarise for each radius R and calculate mean and standard deviation
	dSAR_sum = dSAR_list.groupby(['R'],as_index=False).agg({'S':['mean','std']})
	dSAR_sum.columns = ['R','Smean','Ssd']
	dSAR_sum.to_csv('PCF/Output_' +str(n_sub)+ '/dSAR_' + str(n_sub) + '_discrete.csv')

	cSAR_sum = cSAR_list.groupby(['R'],as_index=False).agg({'S':['mean','std']})
	cSAR_sum.columns = ['R','Smean','Ssd']
	cSAR_sum.to_csv('PCF/Output_' +str(n_sub)+ '/cSAR_' + str(n_sub) + '_discrete.csv')