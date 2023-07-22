# Spatial Biodiversity Patterns
## Master Thesis IMBRSea

The Physics of Biodiversity:
exploring the dynamics behind spatial biodiversity patterns  
<br/>
contact: kobe.simoens@imbrsea.eu  
date: 01/08/2023  

---

## Analytical calculations for the CPR data set

---

### Python scripts

---

#### CPR_spatial.py

Python code to fit the Continuous Analytical Model directly to the disconnected SAR.  
The code contains the following methods:  

- **S_emp()** : calculates the empirical number of species in a subset of samples  

- **SAC_emp()** : calculates the empirical disconnected SAR or SAC  
!!! This method takes a few minutes to run.

- **alpha_function()** : calculates the α parameter of the Γ distribution

- **beta_function()** : calculates the β parameter of the Γ distribution

- **integral_R()** : calculates the integral in the SAR / SAD expression

- **S()** : calculates the predicted number of species at a certain scale

- **SARfit()** : calculates the negative log-likelihood for a given combination of parameters

- **outputLL()** : combines the negative log-likelihood values for a range of parameter values  
!!! This method takes a few hours to run.

---

#### CPR_MF.py

Python code to fit the Mean Field Model to the disconnected SAR.  
Derived quantities are calculated and an attempt is made to use these predictions in the fitting of the Continuous Analytical Model.  
The code contains the following methods:  

- **U()** : calculates the ξp parameter for a number of samples smaller than the total number of samples

- **Spk()** : calculates the number of species in a number of samples smaller than the total number of samples

- **SAR_MF()** : calculates the negative log-likelihood for a given combination of parameters

- **fitMF()** : fits the Mean Field parameters to the empirical disconnected SAR

- **xi0()** : calculates the ξ parameter at the scale of the entire CPR grid

- **S0()** : calculates the number of species at the scale of the entire CPR grid

- **N0()** : calculates the number of individuals at the scale of the entire CPR grid

- **Nmean()** : calculates the Mean Field mean number of individuals per species at the scale of the entire CPR grid

- **Nsigma()** : calculates the spatial variance of the number of individuals per species at the scale of the entire CPR grid

- **r_function()** : calculates the spatial r parameter at the scale of the entire CPR grid

- **xi_function()** : calculates the spatial ξ parameter at the scale of the entire CPR grid

- **calculations()** : performs all the calculations and write out to a csv file

- **dSAR_MF()** : calculates the Mean Field prediction for the disconnected SAR

---

#### CPR_main.py

Executes all the functions in the Python scripts.  
!!! Executing everything will take several hours.

---

### R script

---

#### CPR_plot.R

R script to generate the plots.  
The script requires csv files generated in the Python scripts.  
The code contains the following functions:

- **plotLL()** : plot the likelihood surface for a combination of parameter values  

- **plotdSAR()** : plot the disconnected SAR, both empirical and the Mean Field fit

--- 

### csv files

---

#### CPR_samples_div.csv

This file contains the diversity matrix for the individual samples.  
The csv file is generated in **CPR_extractData.R** of https://github.com/ksimoens/Thesis-Data-Analysis.git.  
Columns are:

- **x** : longitudinal coordinate of the sample  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **y** : latitudinal coordinate of the sample  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **#key#**:  
All remaining columns represent a unique species.  
Column names are the unique GBIF identifiers for the species.  
Values are the total number of entries of a particular species found in a particular sample.  
This number is only indicative and cannot be used as a real abundance.

---

#### CPR_dSAR_emp.csv

This file contains the empirical disconnected SAR as calculated in 'CPR_spatial.SAC_emp()'.  
Columns are:

- **k**: the number of samples included in the subset of samples

- **S**: the mean number of species found in the subset of samples

---

#### CPR_likelihood.csv

This file contains the values of the negative log-likelihood for around 1 million combinations of ρ, λ and n parameters as compiled in 'CPR_spatial.outputLL()'.  
Columns are:

- **r** : ρ parameter value  
Unit: kilometres (km)

- **l** : λ parameter value  
Unit: kilometres (km)

- **n** : n parameter value  
Unit: per square kilometres (km<sup>-2</sup>)

- **LL** : negative log-likelihood value

---

#### CPR_MF.csv

This file contains the output of the Mean Field calculations.  
Columns are:

- **r_MF** : Mean Field fit of the r parameter ('CPR_MF.fitMF()')

- **ξp_MF** : Mean Field fit of the ξp parameter ('CPR_MF.fitMF()')

- **ξ0_MF** : Mean Field prediction of the ξ parameter at the scale of the entire CPR grid ('CPR_MF.xi0()')

- **S0_MF** : Mean Field prediction of the total number of species in the entire CPR grid ('CPR_MF.S0()')

- **N0_MF** : Mean Field prediction of the total number of individuals in the entire CPR grid ('CPR_MF.N0()')

- **n_MF** : Mean Field prediction of the mean number of individuals per species per area in the entire CPR grid  
Unit: per square kilometres (km<sup>-2</sup>)

- **ρ** : Fixed arbitrary choice for the spatial ρ parameter  
Unit: kilometres (km)

- **λ** : Fixed arbitrary choice for the spatial λ parameter  
Unit: kilometres (km)

- **N_mean** : Mean Field prediction of the mean number of individuals per species in the entire CPR grid ('CPR_MF.Nmean()')

- **N_sigma** : Spatial variance of the number of individuals per species in the entire CPR grid ('CPR_MF.Nsigma()')

- **r_SP** : Spatial prediction of the r parameter ('CPR_MF.r_function()')

- **ξ_SP** : Spatial prediction of the ξ parameter in the entire CPR grid ('CPR_MF.xi_function()')

---

#### CPR_dSAR_MF.csv

This file contains the theoretical Mean Field prediction for the disconnected SAR as calculated in 'CPR_MF.dSAR_MF()'.  
Columns are:

- **pk** : fraction of the total number of samples included in the subset of samples

- **Spk** : number of species found in the subset of samples

---

## Instructions

1. Download the files and put them in a local directory
2. Execute the 'runCode.sh' bash file: $ ./runCode.sh  
!!! The full run will take a couple of hours.  
3. or run the methods in CPR_main.py and the CPR_plot.R seperately 
