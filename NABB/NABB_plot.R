# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Make the PLOTs for the NABB data
#####################################


# --------------------- LOAD PACKAGES ----------------------

library(tidyverse)

# ----------------------------------------------------------


# get the empirical PCF
combinePCF <- function(n_sub){

	# list all replicate files
	file_list <- list.files(path=paste0('PCF/PCF_',n_sub),full.names=TRUE)

	# final container
	df_final <- matrix(nrow=0,ncol=2) %>% as.data.frame()
	names(df_final) <- c('distance','PCF')

	# for each replicte
	for(file in file_list){

		rep <- read.csv(file,header=TRUE,row.names=1) %>% 
				dplyr::mutate(distance_red=round(distance/50)*50) %>%
				dplyr::filter(distance_red < 2050 & distance_red > 0) %>%
				dplyr::select(c(distance_red,PCF))

		df_final <- rbind(df_final,rep)

	}

	# average over samples and summarise for distance_red
	df_final <- df_final %>% dplyr::group_by(distance_red) %>%
								dplyr::summarise(PCF=mean(PCF)) %>%
								dplyr::rename(distance=distance_red)

	return(df_final)

}

# theoretical PCF
bessel <- function(R,rho,lambd){

	return(1 + 1/2/pi*(rho/lambd)^2*besselK(R/lambd,0))

}


# get the theoretical PCF
calcPCF <- function(n_sub){

	# get the list with the results of the empirical PCF
	rep_list <- read.csv(paste0('PCF/Output_',n_sub,'/rep_list_',n_sub,'.csv'),header=TRUE,row.names=NULL)
	
	# get a list of distance values
	R_list <- 40:2010

	# create container for final results
	df_final <- matrix(nrow=0,ncol=2) %>% as.data.frame()
	names(df_final) <- c('distance','PCF')

	for(i in 1:nrow(rep_list)){
		df_final <- rbind(df_final,
							data.frame(distance=R_list,PCF=bessel(R_list,rep_list[i,1],rep_list[i,2])) )
	}

	df_final <- df_final %>% dplyr::group_by(distance) %>%
								dplyr::summarise(PCFmean=mean(PCF),PCFsd=sd(PCF))
	
	return(df_final)

}


# plot the empirical and fitted PCF
plotPCF <- function(n_sub){

	# get the empirical PCF
	PCF_emp <- combinePCF(n_sub)
	# get the theoretical PCF
	PCF_the <- calcPCF(n_sub)

	# create a dataframe for the polygon
	PCF_pol <- data.frame(distance=c(PCF_the$distance,rev(PCF_the$distance)),
							rib=c(PCF_the$PCFmean-3*PCF_the$PCFsd,rev(PCF_the$PCFmean+3*PCF_the$PCFsd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=PCF_pol,aes(x=distance,y=rib),fill='red',alpha=0.4) +
					geom_line(data=PCF_the,aes(x=distance,y=PCFmean)) + 
					geom_point(data=PCF_emp,aes(x=distance,y=PCF)) +
					theme_bw() + xlab('distance (km)') + ylab('pair correlation function') + 
					ylim(0,NA)

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/PCF_',n_sub,'.png'),.,device='png',width=15,height=10,units='cm')

}


# calculate the empirical spatial SAD for a subset of samples 
# requires: n_sub = the number of samples in the subset
#			n_rep = the number of replicate sampling for the empirical data
calcsSADdata <- function(n_sub,n_rep){

	# read in the empirical diversity matrix
	dat <- read.csv('NABB_routes_div_2021.csv',header=TRUE,row.names=1) %>% dplyr::select(-c(Longitude,Latitude))

	# create the container for the output
	df_final <- matrix(ncol=2,nrow=0) %>% as.data.frame()
	names(df_final) <- c('logCount','S')

	# replicate the sampling n_rep times
	for(i in 1:n_rep){

		# make the subset of samples
		dat_sub <- dat[sample(1:nrow(dat),n_sub),]

		# combine all counts per species
		df_sum <- dat_sub %>% tidyr::pivot_longer(1:ncol(.),values_to='count',names_to='sp') %>% dplyr::group_by(sp) %>%
							dplyr::summarise(count=sum(count)) %>%	dplyr::mutate(logCount=log2(count))

		# the maximum number of individuals to be considered (log transformed)
		Nmax <- max(df_sum$logCount) %>% ceiling()
		# start from 1 individual
		N <- 0

		# create the output containers
		N_list <- rep(NA,Nmax+1)
		S_list <- rep(NA,Nmax+1)

		# create the first entry (Preston binning)
		# https://doi.org/10.1038/nature04030
		N_list[1] <- N
		S0 <- (df_sum %>% dplyr::filter(logCount==0) %>% nrow())/2
		S_list[1] <- S0

		# for each abundance bin
		for(N in 1:Nmax){
			# calculate the Preston bins
			N_list[N+1] <- N
			Sn <- (df_sum %>% dplyr::filter(logCount==(N-1) | logCount==N) %>% nrow())/2
			Sn <- Sn + (df_sum %>% dplyr::filter(logCount > N-1 & logCount < N) %>% nrow())
			S_list[N+1] <- Sn
		}

		# create the dataframe and add to container
		df <- data.frame(logClass=N_list,S=S_list)
		df_final <- rbind(df_final,df)

	}

	# combine the replicates
	df_plot <- df_final %>% dplyr::group_by(logClass) %>% dplyr::summarise(Smean=mean(S),Ssd=sd(S))
	names(df_plot) <- c('N','Smean','Ssd')

	return(df_plot)

}


# plot the empirical and theoretical sSAD for a subset of samples
# requires: n_sub = the number of samples in the subset
#			n_rep = the number of replicate sampling for the empirical data
plotsSAD_discrete <- function(n_sub){

	# create the empirical sSAD 
	# if all samples are included, downscale to 500 samples
	if(n_sub == 2287){
		df_emp <- calcsSADdata(500,100)
	} else {
		df_emp <- calcsSADdata(n_sub,100)
	}

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/sSAD_',n_sub,'_discrete.csv'),header=TRUE,row.names=1) %>%
				dplyr::filter(N <= max(df_emp$N))

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_bar(data=df_emp,aes(x=N,y=Smean),stat='identity') +
					geom_errorbar(data=df_emp,aes(x=N,ymin=Smean-Ssd,ymax=Smean+Ssd),width=0.1) +
					geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=0:max(df_emp$N),breaks=0:max(df_emp$N)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/sSAD_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the theoretical SAD for the contiguous USA
# requires: n_sub = the number of samples in the subset
plotSAD0_discrete <- function(n_sub){

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/SAD0_',n_sub,'_discrete.csv'),header=TRUE,row.names=1)

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=seq(0,max(df_the$N),2),breaks=seq(0,max(df_the$N),2)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/SAD0_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical sSAD for all samples
# requires: n_sub = the number of samples in the subset
plottSAD_discrete <- function(n_sub){

	# create the empirical sSAD 
	df_emp <- calcsSADdata(2287,1)

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/tSAD_',n_sub,'_discrete.csv'),header=TRUE,row.names=1) %>%
				dplyr::filter(N <= max(df_emp$N))

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_bar(data=df_emp,aes(x=N,y=Smean),stat='identity') +
					geom_errorbar(data=df_emp,aes(x=N,ymin=Smean-Ssd,ymax=Smean+Ssd),width=0.1) +
					geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=0:max(df_emp$N),breaks=0:max(df_emp$N)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/tSAD_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical disconnected SAR
# requires: n_sub = the number of samples in the subset used as starting point
plotdSAR_discrete <- function(n_sub){

	# create the empirical disconnected SAR
	dSAR_emp <- read.csv('NABB_dSAR_emp.csv',header=TRUE,row.names=1)

	# get the theoretical disconnected SAR
	dSAR_the <- read.csv(paste0('PCF/Output_',n_sub,'/dSAR_',n_sub,'_discrete.csv'),header=TRUE,row.names=1)

	# create dataframe for the theory errorbars
	dSAR_pol <- data.frame(R=c(dSAR_the$R,rev(dSAR_the$R)),rib=c(dSAR_the$Smean-3*dSAR_the$Ssd,rev(dSAR_the$Smean+3*dSAR_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=dSAR_pol,aes(x=R,y=rib),fill='red',alpha=0.4) +
					geom_point(data=dSAR_emp,aes(x=R,y=S)) + 
					geom_line(data=dSAR_the,aes(x=R,y=Smean)) + 
					theme_bw() + xlab('R (km)') + ylab('number of species') +
					scale_x_continuous(limits=c(0,NA),labels=seq(0,125,25),breaks=seq(0,125,25)) +
					scale_y_continuous(limits=c(0,NA),labels=seq(0,600,100),breaks=seq(0,600,100))

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/dSAR_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm')

}


# plot the theoretical SAD for the contiguous USA
# requires: n_sub = the number of samples in the subset
plotSAD0_continuous <- function(n_sub){

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/SAD0_',n_sub,'_continuous.csv'),header=TRUE,row.names=1) 

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() +	geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=seq(0,max(df_the$N),2),breaks=seq(0,max(df_the$N),2)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/SAD0_',n_sub,'_continuous.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical SAD for a connected subarea equivalent to all the samples
# requires: n_sub = the number of samples in the subset
plottSAD_continuous <- function(n_sub){

	# create the empirical SAD
	df_emp <- calcsSADdata(2287,1)

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/tSAD_',n_sub,'_continuous.csv'),header=TRUE,row.names=1) %>%
				dplyr::filter(N <= max(df_emp$N))

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_bar(data=df_emp,aes(x=N,y=Smean),stat='identity') +
					geom_errorbar(data=df_emp,aes(x=N,ymin=Smean-Ssd,ymax=Smean+Ssd),width=0.1) +
					geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=0:max(df_emp$N),breaks=0:max(df_emp$N)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/tSAD_',n_sub,'_continuous.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical disconnected SAR
# requires: n_sub = the number of samples in the subset used as starting point
plotdSAR_continuous <- function(n_sub){

	# create the empirical disconnected SAR
	dSAR_emp <- read.csv('NABB_dSAR_emp.csv',header=TRUE,row.names=1)

	# get the theoretical disconnected SAR
	dSAR_the <- read.csv(paste0('PCF/Output_',n_sub,'/dSAR_',n_sub,'_continuous.csv'),header=TRUE,row.names=1)

	# create dataframe for the theory errorbars
	dSAR_pol <- data.frame(R=c(dSAR_the$R,rev(dSAR_the$R)),rib=c(dSAR_the$Smean-3*dSAR_the$Ssd,rev(dSAR_the$Smean+3*dSAR_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=dSAR_pol,aes(x=R,y=rib),fill='red',alpha=0.4) +
					geom_point(data=dSAR_emp,aes(x=R,y=S)) + 
					geom_line(data=dSAR_the,aes(x=R,y=Smean)) + 
					theme_bw() + xlab('R (km)') + ylab('number of species') +
					scale_x_continuous(limits=c(0,NA),labels=seq(0,125,25),breaks=seq(0,125,25)) +
					scale_y_continuous(limits=c(0,NA),labels=seq(0,600,100),breaks=seq(0,600,100))

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/dSAR_',n_sub,'_continuous.png'),.,device='png',width=15,height=10,units='cm')

}

# make the plots
# get n_sub = the number of samples in the subsample to start with
# uncomment to manually run these plots and set n_sub
PAR_file <- file('PCF/PARAM_file.txt',open='r')
on.exit(close(PAR_file))
PAR_file_lines <- readLines(PAR_file)
n_sub <- strtoi(PAR_file_lines)

# create the directory for the plots
dir.create(paste0('PCF/Output_',n_sub,'/plots'))

# plot the PCF only if the empirical PCF are available
if(file.exists(paste0('PCF/PCF_',n_sub,'/PCF_',n_sub,'_0.csv'))){
	plotPCF(n_sub)
}

plotsSAD_discrete(n_sub)
plotSAD0_discrete(n_sub)
plottSAD_discrete(n_sub)
plotdSAR_discrete(n_sub)

plottSAD_continuous(n_sub)
plotSAD0_continuous(n_sub)
plotdSAR_continuous(n_sub)

unlink('PCF/PARAM_file.txt')