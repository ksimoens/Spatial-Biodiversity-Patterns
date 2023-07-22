# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Calculations for the SPATIAL Continuous Model
# 	on the CPR data
#####################################

# ----------------- IMPORT MODULES -------------------------

import pandas as pd
import numpy as np
from scipy import special
from scipy import stats

# ----------------------------------------------------------


# calculate the number of species in k samples of dat_sample
def S_emp(dat_sample,k):
	# get subsample of total dat_sample
	# remove coordinates
	dat_sub = dat_sample.sample(n=k,replace=False).drop(['x','y'],axis=1)

	# count the number of species
	dat_sum = dat_sub.sum(axis=0)
	S = len(dat_sum[dat_sum != 0])

	return(S)


# calculate the empirical disconnected SAR
def SAC_emp():
	# read in the diversity matrix for samples
	dat_sample = pd.read_csv('CPR_samples_div.csv',index_col=0)

	# get species and samples containers
	S_list = np.array([])
	k_list = np.array([])

	# for each 1000 samples
	for k in range(1000,len(dat_sample)+1,1000):
		print(str(k) + ' samples of ' + str(len(dat_sample)))
		# replicate 100 times and add to container
		S_sub_list = np.empty(100)
		for i in range(0,100):
			S_sub_list[i] = S_emp(dat_sample,k)

		# result is the mean of all replicates
		k_list = np.append(k_list,k)
		S_list = np.append(S_list,np.mean(S_sub_list))

	# create output dataframe
	darray = np.concatenate((k_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['k','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	df_out.to_csv('CPR_dSAR_emp.csv')


# α function (https://doi.org/10.1111/2041-210X.12319)
# parameter of the spatial Γ-distribution
# requires the PCF parameters ρ and λ
def alpha_function(R,rho,lambd):
	# distances are scaled with the 'correlation length'
	x = R/lambd
	a = np.pi * pow((R/rho),2) / ( 1. - 2./x * special.k1(x) * special.i1(x) / (special.i0(x)*special.k1(x) + special.i1(x)*special.k0(x)) )
	return(a)

# β function (https://doi.org/10.1111/2041-210X.12319)
# parameter of the spatial Γ-distribution
# requires the PCF parameters ρ and λ
def beta_function(R,rho,lambd,n):
	# distances are scaled with the 'correlation length'
	x = R/lambd
	b = n * pow(rho,2) * ( 1. - 2./x * special.k1(x) * special.i1(x) / (special.i0(x)*special.k1(x) + special.i1(x)*special.k0(x)) )
	return(b)

# integral in SAR / SAD formulae (https://doi.org/10.1111/2041-210X.12319)
# requires the SAD parameters α and β
def integral_R(alpha,beta):
	# integral of the Γ-distribution from 1 to infinity
	# = integral from 0 to infinity minus integral from 0 to 1 
	I_num = (1.-stats.gamma.cdf(x=1.,a=alpha,scale=beta))
	return(I_num)

# calculate the number of species at scale R
# 	give the number of species at scale R0
# formula for SAR (https://doi.org/10.1111/2041-210X.12319)
# requires the PCF parameters ρ and λ
def S(R,rho,lambd,n,R0,S0):

	# α parameter at scale R0
	alpha_0 = alpha_function(R0,rho,lambd)
	# β parameter at scale R0
	beta_0 = beta_function(R0,rho,lambd,n)
	# intagral at scale R0
	# corrected with the number of species S0 at scale R0
	I_0 = integral_R(alpha_0,beta_0) / S0

	# α parameter at scale R
	alpha_R = alpha_function(R,rho,lambd)
	# β parameter at scale R
	beta_R = beta_function(R,rho,lambd,n)
	# number of species at scale R
	S_R = integral_R(alpha_R,beta_R) / I_0

	return(S_R)

# calculate the log-likelihood for a combination of input parameter values
def SARfit(paramFIT):

	# the parameters to be fitted
	# the PCF parameters ρ and λ
	# the total density per species
	rho = paramFIT[0]
	lambd = paramFIT[1]
	n = paramFIT[2]

	# the standard deviation of the statistical distribution
	# fixed as 1 for the plots
	sd = 1

	# get the empirical disconnected SAR
	# created on line 37
	df_data = pd.read_csv('CPR_dSAR_emp.csv')

	# only use points starting from half the total number of samples
	# 	in order to avoid over-stretching the model
	df_data = df_data[df_data['k'] > np.max(df_data['k'])/2]

	# the number of samples
	xdata = df_data['k'].to_numpy()
	# the empirical number of species
	ydata = df_data['S'].to_numpy()

	# container for the estimated number of species
	ypred = np.array([None]*len(xdata))

	# surface area of a single sample (CPR specifications)
	# 	= (tow length) x (horizontal aperture)
	Ai = 19.*0.0127/1000.

	# the fixed parameters
	# S0 is the number of species at the scale R0 of all samples
	R0 = np.sqrt(np.max(xdata)*Ai/np.pi)
	S0 = np.max(ydata)
	
	# for each number of samples in the data
	i = 0
	while(i < len(ypred)):
		k = xdata[i]
		# calculate the surfacte area of the subsample (circle)
		R = np.sqrt(k*Ai/np.pi)
		# estimated number of species
		S_R = S(R,rho,lambd,n,R0,S0)
		ypred[i] = S_R

		i += 1

	# calculate the negative log-likelihood
	LL = -np.sum( stats.norm.logpdf(ydata, loc=ypred, scale=sd ) )

	return(LL)

# calculate the likelihood surface
# takes a couple of hours to run (1000000 combinations)
def outputLL():

	# get containers for the quantities
	rList = np.array([]) # ρ
	lList = np.array([]) # λ
	nList = np.array([]) # n
	LLlist = np.array([]) # log likelihood

	# counter
	c = 0

	# from ρ = 50 to 4950
	r = 50.
	while(r < 5000.):
		# from λ = 50 to 4950
		l = 50.
		while(l < 5000.):
			# from n = 1 to 99
			n = 1.
			while(n < 100.):

				print(f"{c/1e6*100.:.4f}" + ' %',end='\r')

				# combine the parameters
				paramFIT = [r,l,n]
				# calculate the log-likelihood
				LL = SARfit(paramFIT)
				# add to the containers
				rList = np.append(rList,r)
				lList = np.append(lList,l)
				nList = np.append(nList,n)
				LLlist = np.append(LLlist,LL)

				# n in steps of 1
				n += 1.
				c += 1

			# λ in steps of 50
			l += 50.

		# ρ in steps of 50
		r += 50.

	# create dataframe and write out
	darray = np.column_stack((rList,lList,nList,LLlist))
	names_out = ['r','l','n','LL']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	df_out.to_csv('CPR_likelihood.csv')
