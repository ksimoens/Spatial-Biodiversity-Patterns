# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Calculations for the MEAN FIELD Model
# 	on the CPR data
#####################################

# ----------------- IMPORT MODULES -------------------------

import pandas as pd
import numpy as np
from scipy import special
from scipy import stats
from scipy.optimize import minimize

# ----------------------------------------------------------


# calculate ξ for fewer samples p than a reference number of samples p0
# 	given ξp0 (https://doi.org/10.1111/oik.06754)
def U(p,p0,xip0):

	return(p*xip0/(p0+xip0*(p-p0)))


# calculate the number of species in fewer samples pk than a reference number of samples p
# 	given r, ξp and Sp (https://doi.org/10.1111/oik.06754)
def Spk(pk,p,r,xip,Sp):

	# renormalise the fraction so that p = 1
	pkN = pk/p

	# the predicted number of species in k samples
	Spk = (1.-pow(1.-U(pkN,1.,xip),r))/(1.-pow(1.-xip,r))*Sp

	return(Spk)


# calculate the log likelihood of the parameter values given  
# 	the empirical data  
def SAR_MF(paramFIT):

	# r and ξp = parameters of the negative binomial
	r = paramFIT[0]
	xip = paramFIT[1]
	sd = paramFIT[2]

	# get the empirical disconnected SAR
	# calculated in CPR_spatial.py
	df_data = pd.read_csv('CPR_dSAR_emp.csv')

	# get the number of included samples
	xdata = df_data['k'].to_numpy()
	# the total number of samples in the data set
	kmax = np.max(xdata)
	# get the number of species in the samples
	ydata = df_data['S'].to_numpy()
	# the total number of species
	Sp = np.max(ydata)

	# container for the mean field predictions
	ypred = np.array([None]*len(xdata))
	
	# for each number of samples
	i = 0
	while(i < len(ypred)):
		k = xdata[i]
		# the fraction of the total number of samples
		pk = k / kmax
		# the predicted number of species in k samples
		Sk = Spk(pk,1.,r,xip,Sp)		
		# add to the container
		ypred[i] = Sk

		i += 1	

	# calculate the log-likelihood for the current combination of parameters
	LL = -np.sum( stats.norm.logpdf(ydata, loc=ypred, scale=sd ) )

	return(LL)


# fit the r and ξp parameters to the empirical data
def fitMF():

	# containers for the end products
	LL = np.array([]) # log-likelihood
	r = np.array([]) # r values
	xi = np.array([]) # ξp values

	# use different initial values to avoid getting stuck in local minima
	# r0 from 1e-10 to 1e-5
	r0 = 1e-10
	while(r0 < 1e-5):
		# ξp0 from 0.1 to 1
		xi0 = 0.1
		while(xi0 < 1.):

			# define the initial parameter values
			# sd0 = 1
			initParamFIT = [r0,xi0,1.]
			# the optimisation algorithm
			# bound r and sd at 0 at the lower end
			# bound ξp between 0 and 1
			results = minimize(fun=SAR_MF, x0=initParamFIT, bounds=((0.,None),(0.,1.),(0,None)),method='Nelder-Mead')

			LL = np.append(LL,results.fun)
			r = np.append(r,results.x[0])
			xi = np.append(xi,results.x[1])

			# in jumps of 0.1
			xi0 += 0.1
		# in multiplicative jumps of 5
		r0 *= 5

	# get the minimal negative log-likelihood
	i_min = np.nanargmin(LL)
	# and the corresponding parameters
	res = np.array([r[i_min],xi[i_min]])
	
	return(res)


# calculate the ξ parameter at the scale of the entire CPR grid (https://doi.org/10.1111/oik.06754)
# requires	p: the fraction of surface area covered by all the samples
#			ξp: negative binomial parameter for all the samples
def xi0(p,xip):

	return(xip/(p+xip*(1.-p)))


# calculate the number of species at the scale of the entire CPR grid (https://doi.org/10.1111/oik.06754)
# requires 	p: the fraction of surface area covered by all the samples
#			r: negative binomial parameter
#			ξ: negative binomial parameter at the scale of the entire CPR grid
#			ξp: negative binomial parameter for all the samples
#			Sp: number of species in the samples
def S0(p,r,xi,xip,Sp):

	return(Sp*(1.-pow(1.-xi,r))/(1.-pow(1.-xip,r)))


# calculate the number of individuals at the scale of the entire CPR grid (https://doi.org/10.1111/oik.06754)
# requires:	r: negative binomial parameter
#			ξ: negative binomial parameter at the scale of the entire CPR grid
#			S0: number of species at the scale of the entire CPR grid
def N0(r,xi,S0):

	return(S0*1./(1.-pow(1.-xi,r))*r*xi/(1.-xi))


# mean number of individuals per species at the scale of R (https://doi.org/10.1111/2041-210X.12319)
# requires	R: the radius of a circular area with the same surface area as the desired scale area
#			n: the mean number of individuals per species per area in the entire CPR grid
def Nmean(R,n):

	return(n*np.pi*R*R)


# variance of individuals per species at the scale of R (https://doi.org/10.1111/2041-210X.12319)
# requires	R: the radius of a circular area with the same surface area as the desired scale area
#			n: the mean number of individuals per species per area in the entire CPR grid
#			ρ: parameter of the PCF
#			λ: parameter of the PCF
def Nsigma(R,n,rho,lambd):

	# pre-factor
	factor = n*n*rho*rho*np.pi*R*R
	# distances = relative to correlation length
	x = R / lambd
	bessels = ( 1. - 2./x * special.k1(x) * special.i1(x) / (special.i0(x)*special.k1(x) + special.i1(x)*special.k0(x)) )

	return(factor*bessels)


# spatial expression of the r parameter
def r_function(N_mean,N_sigma):

	return( N_mean*N_mean/(N_sigma - N_mean) )


# spatial expression of the ξ parameter
def xi_function(N_mean,N_sigma):

	return( 1. - N_mean/N_sigma )


# do the calculations and write out results
def calculations():

	# fit the Mean Field parameters
	res_MF = fitMF()
	r_MF = res_MF[0]
	xip_MF = res_MF[1]

	# global parameters
	A0 = 6885083 # surface area of the entire CPR grid (water)
	R0 = np.sqrt(A0/np.pi) # corresponding radius
	Ai = 19.*0.0127/1000. # surface area of one sample
	k_max = pd.read_csv('CPR_dSAR_emp.csv')['k'].max() # total number of samples
	S_max = pd.read_csv('CPR_dSAR_emp.csv')['S'].max() # total number of species in the samples

	# calculated parameters
	p = k_max*Ai/A0 # total fraction of the total area that is sampled

	# Mean Field ξ parameter at the scale of the total area
	xi0_MF = xi0(p,xip_MF)

	# Mean Field number of species at the scale of the total area
	S0_MF = S0(p,r_MF,xi0_MF,xip_MF,S_max)

	# Mean Field number of individuals at the scale of the total area
	N0_MF = N0(r_MF,xi0_MF,S0_MF)

	# Mean Field mean number of individuals per species per area in the total area
	n_MF = N0_MF / S0_MF / A0

	# spatial parameters
	# typical values based on other data sets
	rho = 10000.
	lambd = 1000.

	# Mean Field mean number of individuals per species at the scale of the total area
	N_mean = Nmean(R0,n_MF)
	# Spatial variance of the number of individuals per species at the scale of the total area
	N_sigma = Nsigma(R0,n_MF,rho,lambd)

	# Spatial r parameter at the scale of the total area
	r_SP = r_function(N_mean,N_sigma)
	# Spatial ξ parameter at the scale of the total area
	xi_SP = xi_function(N_mean,N_sigma)

	darray = np.column_stack((r_MF,xip_MF,xi0_MF,S0_MF,N0_MF,n_MF,rho,lambd,N_mean,N_sigma,r_SP,xi_SP))
	names_out = ['r_MF','ξp_MF','ξ0_MF','S0_MF','N0_MF','n_MF','ρ','λ','N_mean','N_sigma','r_SP','ξ_SP']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	df_out.to_csv('CPR_MF.csv')


# write out the Mean Field disconnected SAR
def dSAR_MF():

	# fit the Mean Field parameters
	res_MF = fitMF()
	r_MF = res_MF[0]
	xip_MF = res_MF[1]

	# get the total number of samples
	k_max = pd.read_csv('CPR_dSAR_emp.csv')['k'].max()
	# get the total number of species in all the samples
	S_max = pd.read_csv('CPR_dSAR_emp.csv')['S'].max()

	# create a container for the fraction of the number of samples
	pk_list = np.array(range(500,k_max+1)) / k_max

	# calculate the disconnected SAR
	Spk_list = Spk(pk_list,1.,r_MF,xip_MF,S_max)

	# create dataframe and write out
	darray = np.column_stack((pk_list,Spk_list))
	names_out = ['pk','Spk']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	df_out.to_csv('CPR_dSAR_MF.csv')