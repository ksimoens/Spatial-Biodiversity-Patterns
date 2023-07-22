# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Make the PLOTs for the CPR data
#####################################


# --------------------- LOAD PACKAGES ----------------------

library(tidyverse)
library(ggpubr)

# ----------------------------------------------------------


# plot the log-likelihood for combinations of ρ,λ and n
plotLL <- function(){

	# get the likelihood matrix
	df_LL <- read.csv('CPR_likelihood.csv',header=TRUE,row.names=1)

	#############################
	# subset for fixed n
	#############################
	df_n_1_99 <- df_LL %>% dplyr::filter(n %in% c(1,99)) 	
	# get n as factor
	df_n_1_99$n <- factor(df_n_1_99$n)
	levels(df_n_1_99$n) <- c(expression(n == 1~km^-~2), expression(n == 99~km^-~2))
	# make the plot
	pn <- df_n_1_99 %>% ggplot() + geom_tile(aes(x=r,y=l,fill=LL)) + theme_bw() + facet_grid(~n,labeller=label_parsed) +
					scale_fill_viridis_c(name='-log(L)',option='magma',direction=-1,na.value=rgb(0,0,0,0)) +
					xlab('ρ (km)') + ylab('λ (km)') 

	#############################
	# subset for fixed ρ
	#############################
	df_r_1000_4950 <- df_LL %>% dplyr::filter(r %in% c(1000,4950))
	df_r_1000_4950$r <- factor(df_r_1000_4950$r)
	levels(df_r_1000_4950$r) <- c('ρ = 1000 km','ρ = 4950 km')

	pr <- df_r_1000_4950 %>% ggplot() + geom_tile(aes(x=n,y=l,fill=LL)) + theme_bw() + facet_grid(~r) +
					scale_fill_viridis_c(name='-log(L)',option='magma',direction=-1,na.value=rgb(0,0,0,0)) +
					xlab(expression(n~(km^-~2))) + ylab('λ (km)') 

	#############################
	# subset for fixed λ
	#############################
	df_l_1000_2500 <- df_LL %>% dplyr::filter(l %in% c(1000,2500))
	df_l_1000_2500$l <- factor(df_l_1000_2500$l)
	levels(df_l_1000_2500$l) <- c('λ = 1000 km','λ = 2500 km')

	pl <- df_l_1000_2500 %>% ggplot() + geom_tile(aes(x=n,y=r,fill=LL)) + theme_bw() + facet_grid(~l) +
					scale_fill_viridis_c(name='-log(L)',option='magma',direction=-1,na.value=rgb(0,0,0,0)) +
					xlab(expression(n~(km^-~2))) + ylab('ρ (km)') 

	# combine plots
	q <- ggarrange(pn,pr,pl,ncol=1)
	# uncomment to get the png file
	q %>% ggsave('CPR_likelihood.png',.,device='png',width=23,height=30,units='cm')	

}


# plot the disconnected SAR
plotdSAR <- function(){

	# get empirical SAR
	dSAR_emp <- read.csv('CPR_dSAR_emp.csv',header=TRUE,row.names=1) %>% 
					dplyr::mutate(pk = k/max(k))

	# get mean field predicted SAR
	dSAR_MF <- read.csv('CPR_dSAR_MF.csv',header=TRUE,row.names=1)

	# make the plot
	p <- ggplot() + geom_point(data=dSAR_emp,aes(x=pk,y=S)) + 
					geom_line(data=dSAR_MF,aes(x=pk,y=Spk),col='red') +
					xlab('fraction of samples') + ylab('number of species') +
					theme_bw()

	p %>% ggsave('CPR_dSAR.png',.,device='png',width=15,height=10,units='cm',bg='white')

}


# execute the functions
plotLL()
plotdSAR()