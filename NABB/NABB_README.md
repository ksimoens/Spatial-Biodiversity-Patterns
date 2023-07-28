# Spatial Biodiversity Patterns
## Master Thesis IMBRSea

The Physics of Biodiversity:
exploring the dynamics behind spatial biodiversity patterns  
<br/>
contact: kobe.simoens@imbrsea.eu  
date: 01/08/2023  

---

## Analytical calculations for the NABB data set

---

### Python scripts

---

#### NABB_PCF.py

Python code to calculate the empirical patterns: PCF and disconnected SAR.  
The code contains the following methods:   

- **S_emp()** : calculates the number of species in a subset of samples

- **SAC_emp()** : calculates the empirical disconnected SAR

- **subsetSamples()** : returns a random subset of the total number of samples

- **calcPairs()** : returns a list of all the pairs in a subset of samples

- **calcDist()** : calculates the (Euclidean) distance between two samples

- **calcPCF()** : calculates the correlation between two samples -> PCF

- **calcMeanSumSdevS()** : calculates extra quantities of the subset of samples:  
	
	- the mean number of individuals per species
	- the total number of individuals in the subset
	- the variance of the number of individuals per species
	- the total number of species in the subset 

- **createSubsetFrames()** : creates the lists with the subset of samples

- **makeNodeSubset()** : allocates a subset of pairs to a computer node

- **calcPCFnode()** : calculates the PCF for all pairs in the node subset

- **calcPCFall()** : calculates the PCF for all the computer nodes

- **combineNodes()** : combines the output of all the computer nodes

- **bessel()** : calculates the negative log-likelihood for a combination of parameters of the theoritical PCF 

- **fitDistribution()** : fits the parameters of the theoretical PCF to the empirical PCF

- **replicateRun()** : does all the calculations for a given subset of samples and writes out the results

---

#### NABB_discrete.py

Python code to do all the calculations of the Discrete Analytical Model.  
The code upscales the patterns from a random subset of samples and downscales again from upscaled predictions.  
The model combines https://doi.org/10.1111/2041-210X.12319 with https://doi.org/10.1126/sciadv.1701438 .  
The code contains the following methods:  

- **mean()** : calculates the mean number of individuals per species at a certain connected scale  

- **sigma()** : calculates the variance of the number of individuals per species at a certain connected scale

- **sigma_limit()** : as for **sigma()**, but in the mathematical limit for R << λ  

- **xi_function()** : calculates the ξ parameter of the negative binomial SAD at a certain connected scale

- **r_function()** : calculates the r parameter of the negative binomial SAD at a certain connected scale

- **xip()** : calculates the ξp parameter of the negative binomial SAD at a certain disconnected scale

- **Sup()** : calculates the number of species at a larger connected scale 

- **Sdown()** : calculates the number of species at a disconnected smaller scale 

- **S0()** : calculates the upscaled number of species at the total connected scale iteratively

- **c()** : calculates the normalisation factor of the disconnected negative binomial SAD

- **SAD()** : calculates the number of species with abundances between limits at a disconnected scale

- **sSAD()** : calculates the spatial SAD at a disconnected scale

- **SAD0()** : calculates the spatial SAD at the total (connected) scale

- **dSAR()** : calculates the downscaled disconnected SAR

- **calcReplicates()** : does all the calculations for all the replicates of the empirical PCF and combines and writes out the results

---

#### NABB_continuous.py

Python code to do all the calculations of the Continuous Analytical Model.  
The code upscales the patterns from a random subset of samples and downscales again from upscaled predictions.  
The model implements https://doi.org/10.1111/2041-210X.12319 .    
The code contains the following methods:  

- **alpha_function()** : calculates the α parameter of the Γ distribution SAD at a connected scale

- **beta_function()** : calculates the β parameter of the Γ distribution SAD at a connected scale

- **alpha_limit()** : as for **alpha_function()**, but in the mathematical limit for R << λ  

- **beta_limit()** : as for **beta_function()**, but in the mathematical limit for R << λ  

- **integral_R()** : calculates the normalisation of the number of species at a connected scale

- **integral_R_bound()** : calculates the number of species with abundances between limits at a connected scale

- **S()** : calculates the number of species at a connected scale

- **sSAD()** : calculates the spatial SAD at a connected scale

- **dSAR()** : calculates the upscaled disconnected SAR

- **calcReplicates()** : does all the calculations for all the replicates of the empirical PCF and combines and writes out the results

---

#### NABB_main.py

Executes all the functions in the Python scripts and manages the filing system.  
Change the parameters here.  
Default:  

- randomly sample 1000 of the 2287 samples;
- replicate the sampling 20 times;
- use 4 computer nodes.

!!! Executing everything will take a couple of hours depending on the parameters.

---

### R scripts

---

#### NABB_transform.R

Reprojects the WGS 84 coordinates of **NABB_routes_div_2021.csv** to  
Albers projection CRS 5070 coordinates in **NABB_routes_div_2021_proj.csv**.

---

#### NABB_plot.R

R script to generate the plots.  
The script requires csv files generated in the Python scripts.  
The code contains the following functions:

- **combinePCF()** : collects and summarises the empirical PCF from the replicate files

- **bessel()** : the theoretical PCF equation

- **calcPCF()** : calculates and summarises the theoretical PCF from the replicate parameters

- **plotPCF()** : plots the empirical and theoretical PCF

- **calcsSADdata()** : calculates the empirical SAD of a disconnected subset of samples  

- **plotsSAD_discrete()** : plots the spatial SAD at a disconnected scale (Discrete Analytical Model)

- **plotSAD0_discrete()** : plots the SAD at the scale of the total area (Discrete Analytical Model)

- **plottSAD_discrete()** : plots the SAD at the disconnected scale of all the samples in the data set (Discrete Analytical Model)

- **plotdSAR_discrete()** : plots the disconnected SAR (Discrete Analytical Model)

- **plotSAD0_continuous()** : plots the SAD at the scale of the total area (Continuous Analytical Model)

- **plottSAD_continuous()** : plots the spatial SAD at the connected scale of all the samples (Continuous Analytical Model)

- **plotdSAR_continuous()** : plots the disconnected SAR (Continuous Analytical Model)

--- 

### csv files

---

#### NABB_routes_div_2021.csv

This file contains the diversity matrix for the individual routes.  
The total data set contains 2287 routes with approximate surface area of 50 x 0.4 x 0.4 x π km<sup>2</sup>.  
The csv file is generated in **NABB_extractData.R** of https://github.com/ksimoens/Thesis-Data-Analysis.git .  
Columns are:

- **Longitude** : longitudinal coordinate of the sample  
Projection: WGS 84

- **Latitude** : latitudinal coordinate of the sample  
Projection: WGS 84

- **#key#**:  
All remaining columns represent a unique species.  
Column names are the unique AOU identifiers for the species.  
Values are the total number of individuals of a particular species counted along a particular route.

---

#### NABB_routes_div_2021_proj.csv

As for **NABB_routes_div_2021.csv** but with the coordinates in the Albers projection CRS 5070.  
Generated in **NABB_transform.R**. Columns are: 

- **x** : longitudinal coordinate of the sample  
Projection: Albers projection CRS 5070

- **y** : latitudinal coordinate of the sample  
Projection: Albers projection CRS 5070

- **#key#**:  
All remaining columns represent a unique species.  
Column names are the unique AOU identifiers for the species.  
Values are the total number of individuals of a particular species counted along a particular route.

---

#### NABB_dSAR_emp.csv

This file contains the empirical disconnected SAR as calculated in **NABB_PCF()**.  
Columns are:

- **R**: the radius of the combined disconnected area  
Unit: kilokilometres (km)

- **S**: the mean number of species found in the disconnected area

---

### Directories

---

#### PCF

This directory contains all the output. The results are automatically written to subdirectories of this directory.

---

#### PCF/PCF_X

This directory contains the replicate runs of the calculation of the empirical PCF.  
'X' denotes the number of samples from which the PCF is calculated.  
The replicate runs are stored in:  
**PCF_X_Y.csv**, in which 'Y' denotes the replicate number.  
Columns are:

- **distance**: inbetween distance of a pair of samples  
Unit: kilokilometres (km)

- **PCF** : mean inbetween correlation of a pair of samples at a certain inbetween distance

---

#### PCF/Output_X

This directory contains all the output of the calculations in the Python code and the plots from the R scripts.  
'X' denotes the number of samples from which the PCF is calculated.  
The directory contains the following csv files:

- **rep_list_X.csv** : the empirical parameters of the PCF calculations for all replicates = rows (**NABB_PCF.py**)  
Columns are:

	- *rho* : the fitted ρ parameter of the PCF  
	Unit: kilokilometres (km)
	- *lambd* : the fitted λ parameter of the PCF  
	Unit: kilokilometres (km)
	- *S* : the number of species in the subset of samples
	- *Nsum* : the number of individuals in the subset of samples
	- *Nmean* : the mean number of individuals per species in the subset of samples
	- *Nsd* : the variance of the number of individuals per species in the subset of samples

- **rep_list_X_discrete.csv** : the updated replicate list with the theoretical quantities added (**NABB_discrete.py**)  
Columns are:  

	- *rho* : the fitted ρ parameter of the PCF  
	Unit: kilokilometres (km)
	- *lambd* : the fitted λ parameter of the PCF  
	Unit: kilokilometres (km)
	- *S* : the number of species in the subset of samples
	- *Nsum* : the number of individuals in the subset of samples
	- *Nmean* : the mean number of individuals per species in the subset of samples
	- *Nsd* : the variance of the number of individuals per species in the subset of samples
	- *S0* : the number of species at the scale of the contiguous USA
	- *np* : the mean density per species in the subset of samples  
	Unit: per square kilokilometres (km<sup>-2</sup>)
	- *n* : the mean density per species at the scale of the total area   
	Unit: per square kilokilometres (km<sup>-2</sup>)
	- *St* : the number of species in all samples of the data set
	- *r* : the r parameter of the negative binomial SAD at the connected scale of the total area  
	- *xi* : the ξ parameter of the negative binomial SAD at the connected scale of the total area
	- *xip* : the ξp parameter of the negative binomial SAD at the disconnected scale of the subset of samples

- **rep_list_X_continuous.csv** : the updated replicate list with the theoretical quantities added (**NABB_continuous.py**)  
Columns are:  

	- *rho* : the fitted ρ parameter of the PCF  
	Unit: kilokilometres (km)
	- *lambd* : the fitted λ parameter of the PCF  
	Unit: kilometres (km)
	- *S* : the number of species in the subset of samples
	- *Nsum* : the number of individuals in the subset of samples
	- *Nmean* : the mean number of individuals per species in the subset of samples
	- *Nsd* : the variance of the number of individuals per species in the subset of samples
	- *S0* : the number of species at the scale of the contiguous USA
	- *np* : the mean density per species in the subset of samples  
	Unit: per square kilometres (km<sup>-2</sup>)
	- *n* : the mean density per species at the scale of the total area   
	Unit: per square kilometres (km<sup>-2</sup>)
	- *St* : the number of species in all samples of the data set
	- *alpha0* : the α parameter of the Γ distribution SAD at the scale of the total area
	- *beta0* : the β parameter of the Γ distribution SAD at the scale of the total area
	- *Sdown* : the number of species at the connected scale of X samples

- **SAD0_X_discrete.csv** : the SAD at the scale of the contiguous USA (**NABB_discrete.py**)  
Columns are:

	- *N* : Preston classes of the abundances  
	- *Smean* : mean number of species with abundances in the Preston bin
	- *Ssd* : standard deviation of the number of species with abundances in the Preston bin   

- **sSAD_X_discrete.csv** : the spatial SAD at the disconnected scale of X samples (**NABB_discrete.py**)  
Columns are:

	- *N* : Preston classes of the abundances  
	- *Smean* : mean number of species with abundances in the Preston bin
	- *Ssd* : standard deviation of the number of species with abundances in the Preston bin  

- **tSAD_X_discrete.csv** : the spatial SAD at the disconnected scale of all samples (**NABB_discrete.py**)  
Columns are:

	- *N* : Preston classes of the abundances  
	- *Smean* : mean number of species with abundances in the Preston bin
	- *Ssd* : standard deviation of the number of species with abundances in the Preston bin  

- **dSAR_X_discrete.csv** : the disconnected SAR (**NABB_discrete.py**)  
Columns are: 

	- *R* : radius of a disconnected area  
	Unit: kilometres (km)
	- *Smean* : mean number of species in a disconnected area with radius *R*
	- *Ssd* : standard deviation of the number of species in a disconnected area with radius *R*    

- **SAD0_X_continuous.csv** : the SAD at the scale of the contiguous USA (**NABB_continuous.py**)  
Columns are:

	- *N* : Preston classes of the abundances  
	- *Smean* : mean number of species with abundances in the Preston bin
	- *Ssd* : standard deviation of the number of species with abundances in the Preston bin 

- **tSAD_X_continuous.csv** : the spatial SAD at the connected scale of all samples (**NABB_continuous.py**)  
Columns are:

	- *N* : Preston classes of the abundances  
	- *Smean* : mean number of species with abundances in the Preston bin
	- *Ssd* : standard deviation of the number of species with abundances in the Preston bin  

- **dSAR_X_continuous.csv** : the upscaled disconnected SAR (**NABB_continuous.py**)  
Columns are: 

	- *R* : radius of a disconnected area  
	Unit: kilometres (km)
	- *Smean* : mean number of species in a disconnected area with radius *R*
	- *Ssd* : standard deviation of the number of species in a disconnected area with radius *R*   

---

#### PCF/Output_X/plots

This directory contains all the plots generated in **NABB_plot.R**.  
'X' denotes the number of samples from which the PCF is calculated.  
The directory contains the png files corresponding with the csv files in the output.  
The names of the png files are identical to the names of the csv files.  
Additionally:

- **PCF_X.png** : plot of the summarised empirical and theoretical PCF

---

#### Output_1000 + PCF_1000

The project shows an example of a run for a subset of 1000 samples.  
The sampling is replicated 20 times. Only one PCF file is present due to the size of these files.  
The directories contain all files listed above.

#### Output_500 + Output_1500 + Output_2287

The results of the analysis with 500, 1500 and 2287 samples as starting point.  
The sampling is replicated 20 times for 500 and 1500.  
The files only contain the rep_list_X.csv files with the parameter results. 




## Instructions

1. Download the files and put them in a local directory
2. Execute the 'runCode.sh' bash file: $ ./runCode.sh  
!!! The full run can take a couple of hours depending on the parameters.  
!!! The programme uses 4 computer nodes. Change this in **NABB_main.py** if necessary.  
3. or run the methods in **NABB_main.py** and **NABB_plot.R*** seperately  
4. Change the parameters in **NABB_main.py** as desired.  
The code should work for any number of starting samples *n_sub*.    
For now, the code only runs for the 2287 routes of year 2021 in **NABB_routes_div_2021.csv**.
