# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Calculate the empirical PCF for the BCI data set
#####################################


# ----------------- IMPORT MODULES -------------------------

import pandas as pd
import numpy as np 
import multiprocessing as mp
import os
import shutil
import re
import random
from scipy.optimize import minimize
from scipy import special
from scipy import stats
from csv import writer
from functools import partial
import matplotlib.pyplot as plt

# ----------------------------------------------------------


# calculate the number of species in connected area of size (Lx/kx) x (Ly/ky)
# requires: dat_sample = the total list of sampled to be divided
#			kx = the division in the horizontal direction (= 2 -> divide in half)
#			ky = the division in the vertical direction (= 2 -> divide in half)
def S_conn(dat_sample,kx,ky):

	# get size of subareas
	Lkx = 1000. / kx
	Lky = 500. / ky

	# subdivede the total area in equal connected subareas
	dat_sample['cellx'] = np.floor(dat_sample['x'] / Lkx) 
	dat_sample['celly'] = np.floor(dat_sample['y'] / Lky) 

	# create a unique identifier for each subarea
	dat_sample = dat_sample.astype({'cellx':'int','celly':'int'})
	dat_sample['cell'] = dat_sample['cellx'].astype(str) + '_' + dat_sample['celly'].astype(str)
	# rearrange the dataframe
	dat_sample = dat_sample.drop(['x','y','cellx','celly'],axis=1)
	cols = list(dat_sample.columns)
	cols = [cols[-1]] + cols[:-1]
	dat_sample = dat_sample[cols]

	# combine the number of individuals for each species and cell
	dat_sum = dat_sample.groupby(['cell']).sum()

	# return the mean number of species in a cell
	return(np.mean(dat_sum.astype(bool).sum(axis=1)))


# calculate the empirical connected SAR
def SAR_emp():

	# read in the diversity matrix for samples
	dat_sample = pd.read_csv('BCI_grid_25_25.csv',index_col=0)

	# transform the number of samples to a radius in metres
	Lx = 1000. 			# horizontal size of the total area (m)
	Ly = 500. 			# vertical size of the total area (m)

	# get the total number of species and initialise containers
	S_list = np.array([len(dat_sample.columns)-2])
	R_list = np.array([np.sqrt(Lx*Ly/np.pi)])

	# manually create subsets

	# area in half (two ways to do that)
	S = np.mean([S_conn(dat_sample,2,1),S_conn(dat_sample,1,2)])
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/2.*Ly/np.pi))

	# area in four (two ways considered)
	S = np.mean([S_conn(dat_sample,2,2),S_conn(dat_sample,4,1)])
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/4.*Ly/np.pi))

	# area in eight (one way considered)
	S = S_conn(dat_sample,4,2)
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/4.*Ly/2./np.pi))

	# area in ten (one way considered)
	S = S_conn(dat_sample,5,2)
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/5.*Ly/2./np.pi))

	# area in sixteen (two ways considered)
	S = np.mean([S_conn(dat_sample,4,4),S_conn(dat_sample,8,2)])
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/4.*Ly/4./np.pi))

	# area in 20 (two ways considered)
	S = np.mean([S_conn(dat_sample,5,4),S_conn(dat_sample,4,5)])
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/5.*Ly/4./np.pi))

	# area in 25 (one way)
	S = np.mean([S_conn(dat_sample,5,5),S_conn(dat_sample,5,5)])
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/5.*Ly/5./np.pi))

	# area in 32 (one way considered)
	S = S_conn(dat_sample,8,4)
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/8.*Ly/4./np.pi))

	# area in 50 (one way considered)
	S = S_conn(dat_sample,10,5)
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/10.*Ly/5./np.pi))

	# area in 100 (two ways considered)
	S = np.mean([S_conn(dat_sample,20,5),S_conn(dat_sample,10,10)])
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/20.*Ly/5./np.pi))

	# area in 200 (one way considered)
	S = S_conn(dat_sample,20,10)
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/20.*Ly/10./np.pi))

	# area in 400 (one way considered)
	S = S_conn(dat_sample,20,20)
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/20.*Ly/20./np.pi))

	# area in 800 = single samples (one way)
	S = S_conn(dat_sample,40,20)
	S_list = np.append(S_list,S)
	R_list = np.append(R_list,np.sqrt(Lx/40.*Ly/20./np.pi))

	# create output dataframe
	darray = np.concatenate((R_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['R','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	df_out.to_csv('BCI_cSAR_emp.csv')


# calculate the number of species in k samples of dat_sample
# requires: dat_sample = the total list of samples to be subsetted
#			k = the number of samples in the subset
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
	dat_sample = pd.read_csv('BCI_grid_25_25.csv',index_col=0)

	# get species and samples containers
	S_list = np.array([])
	k_list = np.array([])

	# for each 20 samples
	for k in range(20,len(dat_sample)+1,20):
		print(str(k) + ' samples of ' + str(len(dat_sample)))
		# replicate 100 times and add to container
		S_sub_list = np.empty(100)
		for i in range(0,100):
			S_sub_list[i] = S_emp(dat_sample,k)

		# result is the mean of all replicates
		k_list = np.append(k_list,k)
		S_list = np.append(S_list,np.mean(S_sub_list))

	# transform the number of samples to a radius in metres
	Ai = 25*25 						# surface area of one sample (m^2)
	Ak = k_list*Ai 					# combined surface area of samples (m^2)
	R_list = np.sqrt(Ak/np.pi)		# corresponding radius (m)

	# create output dataframe
	darray = np.concatenate((R_list,S_list)).reshape((-1, 2), order='F')
	names_out = ['R','S']
	df_out = pd.DataFrame(data=darray, columns=names_out)

	df_out.to_csv('BCI_dSAR_emp.csv')


# randomly select n samples from the total number of samples
def subsetSamples(n):

	# get the diversity matrix for the samples
	df_all = pd.read_csv('BCI_grid_25_25.csv',index_col=0)

	# select the samples
	subset = random.sample(range(1,len(df_all)+1),n)
	df_sub = df_all[df_all.index.isin(subset)]
	# remove species that are absent from the subset
	df_sub = df_sub.loc[:,(df_sub.sum(axis=0) != 0)]

	return(df_sub)


# get a list of all the unique sample pairs in the subset
def calcPairs(df_sub):

	# create the containers for the pair indices
	L = len(df_sub)
	# list of indices in the subset
	index_list = df_sub.index
	pair_list_A = np.empty(int(L*(L-1)/2))
	pair_list_B = np.empty(int(L*(L-1)/2))

	# for all unique pairs
	k = 0
	for i in range(0,L-1):
		for j in range(i+1,L):
			# get indices of both samples
			pair_list_A[k] = int(index_list[i])
			pair_list_B[k] = int(index_list[j])

			k += 1

	# create a dataframe with the indices of the pairs
	pair_list = pd.DataFrame({'A':pair_list_A,'B':pair_list_B},columns=['A','B'])
	pair_list = pair_list.astype({'A':'int','B':'int'})

	return(pair_list)


# calculate the Euclidean distance between pairs
def calcDist(df_dist):

	df_dist_T = df_dist.transpose()
	df_dist_T['square'] = pow(df_dist_T.iloc[:,0] - df_dist_T.iloc[:,1],2)

	dist_pair = np.sqrt(np.sum(df_dist_T['square']))

	return(dist_pair)


# calculate the empirical PCF between pairs
# https://doi.org/10.1111/2041-210X.12319
def calcPCF(df_PCF):

	df_PCF_T = df_PCF.transpose()

	# calculate the correlation factor
	df_PCF_T['correl'] = df_PCF_T.iloc[:,0]*df_PCF_T.iloc[:,1]
	# normalise with the mean number of individuals in both samples
	PCF_pair = np.mean(df_PCF_T.iloc[:,2])/np.mean(df_PCF_T.iloc[:,0])/np.mean(df_PCF_T.iloc[:,1])

	return(PCF_pair)


# calculate additional quantities:
# 	the number of species
# 	the total number of individuals
# 	the mean number of individuals per species
# 	the variance of the total number of individuals per species
def calcMeanSumSdevS(df_sub):

	S = len(df_sub.columns) - 2

	# calculate all the individuals per species
	colSums = df_sub.sum().to_frame()
	colSums = colSums.iloc[2:len(colSums),]
	colSums = colSums.rename(columns={colSums.columns[0]:'count'},errors='raise')
	colSums.columns = colSums.columns.astype(str)
	# calculate the squared difference with the mean
	colSums['var'] = pow(colSums['count'] - np.mean(colSums['count']),2)

	Nsum = int(np.sum(colSums['count']))
	Nmean = np.mean(colSums['count'])
	Nsd = np.mean(colSums['var'])

	return(np.array([S,Nsum,Nmean,Nsd]))


# create the subset and subsequent lists
# requires:	n_sub = number of samples in the subset
#			n_cpu = number of computer cores
def createSubsetFrames(n_sub,n_cpu):

	# create the 'PCF'
	if(not os.path.exists("PCF")):
		os.mkdir("PCF")

	# create the output directory and remove if it existed already
	if(os.path.exists("PCF/PCF_out")):
		shutil.rmtree("PCF/PCF_out")
	os.mkdir("PCF/PCF_out")
	
	# create the subset
	df_sub = subsetSamples(n_sub)
	# create a list of the sample pairs in the subset
	pair_list = calcPairs(df_sub)
	
	# calculate an indication of the the number of pairs per core
	n_node = int(np.floor(len(pair_list) / n_cpu))

	return([df_sub,pair_list,n_node])


# create a subset of the pairs for a single computer node 
# requires:	start = lowest index of subset in the total list
#			pair_list = list of sample pairs in the subset
def makeNodeSubset(start,pair_list,n_node):

	# if n_node cannot fit anymore
	if (len(pair_list) - 2*n_node) < start:
		# the node receives the remaining pairs in the list
		pair_node = pair_list[pair_list.index.isin(range(start,len(pair_list)))]

	else:
		# the node receives n_node pairs
		end = start + (n_node-1)
		pair_node = pair_list[pair_list.index.isin(range(start,end+1))]

	return(pair_node)


# calculate the PCF and distance of all pairs in the node
# requires: start = lowest index of subset in the total list 
#			pair_list = list of sample pairs in the subset
#			df_sub = subset of samples
#			n_node = number of pairs for each computer core
def calcPCFnode(start,pair_list,df_sub,n_node):

	# get the pairs to calculate
	pair_node = makeNodeSubset(start,pair_list,n_node)

	# create containers for the distances and PCF 
	dist_node = np.empty(len(pair_node))
	PCF_node = np.empty(len(pair_node))

	# for each pair
	for i in range(0,len(pair_node)):

		# calculate PCF and distance of the pair
		df_pair = df_sub[df_sub.index.isin([pair_node.iloc[i]['A'],pair_node.iloc[i]['B']])]
			
		dist_node[i] = calcDist(df_pair[['x','y']])

		PCF_node[i] = calcPCF(df_pair.iloc[:,2:])

		print('pair ' + str(i+1) + ' of ' + str(len(pair_node)),end='\r')
		
	# create a dataframe for the node and write out
	# tagged with unique 'start' identifier
	df_out = pd.DataFrame({'distance':dist_node,'PCF':PCF_node},columns=['distance','PCF'])
	df_out.to_csv('PCF/PCF_out/PCF_' + str(int(start/n_node)) + '.csv')


# calculate the PCF for all pairs
# distribute the pairs over the computer nodes
def calcPCFall(n_cpu,n_node,pair_list,df_sub):

	# create a list for all the 'start' values
	start_list = np.empty(n_cpu)
	for i in range(0,n_cpu):
		start_list[i] = int(i*n_node)

	start_list = start_list.astype(int)

	calcPCFnode.cte = partial(calcPCFnode,pair_list=pair_list,df_sub=df_sub,n_node=n_node)
		
	# create the parallel computer
	pool = mp.Pool(n_cpu)

	pool.map(calcPCFnode.cte, [s for s in start_list])

	pool.close()


# combine the results of all the computer nodes
# requires: n_sub = number of samples in the subset
#			rep = replicate number
def combineNodes(n_sub,rep):

	# create a dataframe container for the results
	df_final = pd.DataFrame({'distance':[],'PCF':[]},columns=['distance','PCF'])

	# for each output file
	with os.scandir('PCF/PCF_out/') as l:
		for file in l:
			# read the results
			df_node = pd.read_csv(file,index_col=0)
			# add them to the container
			df_final = pd.concat([df_final, df_node], ignore_index=True)

	# remove the 'PCF_out' directory
	shutil.rmtree("PCF/PCF_out")

	# create the n_sub directory if it does not exist
	if(not os.path.exists("PCF/PCF_" + str(n_sub))):
		os.mkdir("PCF/PCF_" + str(n_sub))
	# write out the replicate PCF calculations
	df_final.to_csv('PCF/PCF_' + str(n_sub) + '/PCF_' + str(n_sub) + '_' + str(rep) + '.csv')
	

# calculate the negative log-likelihood of a combination of parameters 
# theoretical PCF: https://doi.org/10.1111/2041-210X.12319
def bessel(params,params_rep):

	# parameters to be fitted
	rho = params[0]
	lambd = params[1]   
	sd = params[2]

	# fixed parameters
	n_sub = params_rep[0]
	rep = params_rep[1]

	# get the replicate results
	df_data = pd.read_csv('PCF/PCF_' + str(n_sub) + '/PCF_' + str(n_sub) + '_' + str(rep) + '.csv')
		
	# calculate the mean PCF for each unique distance
	df_data = df_data.groupby('distance')['PCF'].mean().reset_index()
	# only include distances smaller than 500 m to avoid spurious correlations
	df_data = df_data[(df_data['distance']<500)]

	# the empirical distance
	xdata = df_data['distance'].to_numpy()
	# the empirical PCF
	ydata = df_data['PCF'].to_numpy()
	# the theoretical PCF
	yPred = 1. + 1./2./np.pi*pow(rho/lambd,2)*special.k0(xdata/lambd)

	# calculate negative log likelihood
	LL = -np.sum( stats.norm.logpdf(ydata, loc=yPred, scale=sd ) )

	return(LL)


# fit the theoretical PCF to the empirical PCF
# requires:	n_sub = number of samples in the subset
#			rep = replicate number
def fitDistribution(n_sub,rep):

	# initial parameters to be fitted
	initParams = [10e6, 3e6, 1.]
	# fixed parameters
	params_rep = [n_sub,rep]

	# minimise the negative log-likelihood
	# ρ,λ and sd bounded at zero (non-negative)
	results = minimize(bessel, initParams, bounds=((0,None),(0,None),(0,None)), method='Nelder-Mead',args=params_rep)

	return(results.x[[0,1]])

def replicateRun(n_sub,rep,n_cpu):

	# create the necessary data frames and lists
	sub_list = createSubsetFrames(n_sub,n_cpu)
	df_sub = sub_list[0]
	pair_list = sub_list[1]
	n_node = sub_list[2]

	# calculate the PCF for all pairs -> n_cpu csv files
	calcPCFall(n_cpu,n_node,pair_list,df_sub)
	# combine the files to one csv file
	combineNodes(n_sub,rep)

	# fit ρ and λ of the theoretical PCF
	rhoLambda = np.array(fitDistribution(n_sub,rep))
	# get the other quantities
	meanSum = np.array(calcMeanSumSdevS(df_sub))
	# combine the output in one line
	output_line = np.concatenate((rhoLambda,meanSum))

	# if the output csv file does not exist, create it and write the column names
	if(not os.path.exists('PCF/Output_' +str(n_sub) +'/rep_list_' + str(n_sub) + '.csv')):
		df_names = pd.DataFrame({'rho':[],'lambd':[],'S':[],'Nsum':[],'Nmean':[],'Nsd':[]},columns=['rho','lambd','S','Nsum','Nmean','Nsd'])
		df_names.to_csv('PCF/Output_' +str(n_sub) +'/rep_list_' + str(n_sub) + '.csv',index=False)

	# write a new line to the output file
	with open('PCF/Output_' +str(n_sub)+ '/rep_list_' + str(n_sub) + '.csv','a') as f:
		writer_object = writer(f)
		writer_object.writerow(output_line)
		f.close()

	#os.remove('PCF_' + str(n_sub) + '/PCF_' + str(n_sub) + '_' + str(rep) + '.csv')