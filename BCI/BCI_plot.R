# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# Make the PLOTs for the BCI data
#####################################


# --------------------- LOAD PACKAGES ----------------------

library(tidyverse)

# ----------------------------------------------------------


# calculate the empirical spatial SAD for a subset of samples 
# requires: n_sub = the number of samples in the subset
#			n_rep = the number of replicate sampling for the empirical data
calcsSADdata <- function(n_sub,n_rep){

	# read in the empirical diversity matrix
	dat <- read.csv('BCI_grid_25_25.csv',header=TRUE,row.names=1) %>% dplyr::select(-c(x,y))

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
plotsSAD_discrete <- function(n_sub,n_rep){

	# create the empirical sSAD 
	df_emp <- calcsSADdata(n_sub,n_rep)

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/sSAD_',n_sub,'_discrete.csv'),header=TRUE,row.names=1) %>%
				dplyr::filter(N <= max(df_emp$N))

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_bar(data=df_emp,aes(x=N,y=Smean),stat='identity') +
					geom_errorbar(data=df_emp,aes(x=N,ymin=Smean-3*Ssd,ymax=Smean+3*Ssd),width=0.1) +
					geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=0:max(df_emp$N),breaks=0:max(df_emp$N)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/sSAD_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical SAD for the total area
# requires: n_sub = the number of samples in the subset
plotSAD0_discrete <- function(n_sub){

	# create the empirical SAD
	df_emp <- calcsSADdata(800,1)

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/SAD0_',n_sub,'_discrete.csv'),header=TRUE,row.names=1) %>%
				dplyr::filter(N <= max(df_emp$N))

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_bar(data=df_emp,aes(x=N,y=Smean),stat='identity') +
					geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=0:max(df_emp$N),breaks=0:max(df_emp$N)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/SAD0_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# extract connected areas from the data
# requires: n_sub = the number of subsamples in an equivalent disconnected area
extractConnAreas <- function(n_sub){

	# the size of the entire area (m)
	Lx <- 1000
	Ly <- 500
	# number of connected cells in the entire area
	Ncell <- 800 / n_sub
	# get the division of the area in the two directions
	kx <- ceiling(sqrt(Ncell))
	ky <- Ncell/kx
	while(ky%%1!=0){
		kx <- kx + 1
		ky <- Ncell/kx
	}

	# get all the samples and create the cells
	dat_samples <- read.csv('BCI_grid_25_25.csv',header=TRUE,row.names=1) %>%
					dplyr::mutate(cellx=floor(x/(Lx/kx)),celly=floor(y/(Ly/ky))) %>%
					dplyr::mutate(cell=paste0(cellx,'_',celly)) %>%
					dplyr::select(-c(x,y,cellx,celly))
	# combine the counts for each cell
	dat_cells <- dat_samples %>% dplyr::group_by(cell) %>% dplyr::summarise_all(sum) %>%
					dplyr::select(-c(cell))

	return(dat_cells)

}


# calculate the empirical spatial SAD for a connected subarea 
# requires: n_sub = the number of samples in the connected subarea
calccSADdata <- function(n_sub,n_samples){

	# get the counts per cell
	if(n_sub == 800){
		dat_cells <- extractConnAreas(100)
	} else {
		dat_cells <- extractConnAreas(n_sub)
	}

	# create the container for the output
	df_final <- matrix(ncol=2,nrow=0) %>% as.data.frame()
	names(df_final) <- c('logCount','S')

	# replicate the sampling n_rep times
	for(i in 1:nrow(dat_cells)){

		# make the subset of the cell
		dat_sub <- dat_cells[i,]

		# combine all counts per species
		df_sum <- dat_sub %>% tidyr::pivot_longer(1:ncol(.),values_to='count',names_to='sp') %>% 
						dplyr::mutate(logCount=log2(count))

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


# plot the empirical and theoretical SAD for a connected subarea with the same area as the samples
# requires: n_sub = the number of samples in the subset
plotcSAD_discrete <- function(n_sub){

	# create the empirical SAD
	df_emp <- calccSADdata(n_sub)

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/cSAD_',n_sub,'_discrete.csv'),header=TRUE,row.names=1) %>%
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

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/cSAD_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical disconnected SAR
# requires: n_sub = the number of samples in the subset used as starting point
plotdSAR_discrete <- function(n_sub){

	# create the empirical disconnected SAR
	dSAR_emp <- read.csv('BCI_dSAR_emp.csv',header=TRUE,row.names=1)

	# get the theoretical disconnected SAR
	dSAR_the <- read.csv(paste0('PCF/Output_',n_sub,'/dSAR_',n_sub,'_discrete.csv'),header=TRUE,row.names=1)

	# create dataframe for the theory errorbars
	dSAR_pol <- data.frame(R=c(dSAR_the$R,rev(dSAR_the$R)),rib=c(dSAR_the$Smean-3*dSAR_the$Ssd,rev(dSAR_the$Smean+3*dSAR_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=dSAR_pol,aes(x=R,y=rib),fill='red',alpha=0.4) +
					geom_point(data=dSAR_emp,aes(x=R,y=S)) + 
					geom_line(data=dSAR_the,aes(x=R,y=Smean)) + 
					theme_bw() + xlab('radius (m)') + ylab('number of species')

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/dSAR_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm')

}


# plot the empirical and theoretical connected SAR
# requires: n_sub = the number of samples in the subset used as starting point
plotcSAR_discrete <- function(n_sub){

	# create the empirical connected SAR
	cSAR_emp <- read.csv('BCI_cSAR_emp.csv',header=TRUE,row.names=1)

	# get the theoretical connected SAR
	cSAR_the <- read.csv(paste0('PCF/Output_',n_sub,'/cSAR_',n_sub,'_discrete.csv'),header=TRUE,row.names=1)

	# create dataframe for the theory errorbars
	cSAR_pol <- data.frame(R=c(cSAR_the$R,rev(cSAR_the$R)),rib=c(cSAR_the$Smean-3*cSAR_the$Ssd,rev(cSAR_the$Smean+3*cSAR_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=cSAR_pol,aes(x=R,y=rib),fill='red',alpha=0.4) +
					geom_point(data=cSAR_emp,aes(x=R,y=S)) + 
					geom_line(data=cSAR_the,aes(x=R,y=Smean)) + 
					theme_bw() + xlab('radius (m)') + ylab('number of species')

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/cSAR_',n_sub,'_discrete.png'),.,device='png',width=15,height=10,units='cm')

}


# plot the empirical and theoretical SAD for the total area
# requires: n_sub = the number of samples in the subset
#			n_samples = the total number of samples
plotSAD0_continuous <- function(n_sub){

	# create the empirical SAD
	df_emp <- calcsSADdata(800,1)

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/SAD0_',n_sub,'_continuous.csv'),header=TRUE,row.names=1) %>%
				dplyr::filter(N <= max(df_emp$N))

	# create dataframe for the theory errorbars
	df_pol <- data.frame(N=c(df_the$N,rev(df_the$N)),rib=c(df_the$Smean-3*df_the$Ssd,rev(df_the$Smean+3*df_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_bar(data=df_emp,aes(x=N,y=Smean),stat='identity') +
					geom_polygon(data=df_pol,aes(x=N,y=rib),fill='red',alpha=0.4) +
					geom_line(data=df_the,aes(x=N,y=Smean)) + 
					scale_x_continuous(labels=0:max(df_emp$N),breaks=0:max(df_emp$N)) +
					theme_bw() + xlab('log2(upper abundance)') + ylab('number of species') + 
					theme(panel.grid.minor.x=element_blank()) 

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/SAD0_',n_sub,'_continuous.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical SAD for a connected subarea
# requires: n_sub = the number of samples in the subset
plotcSAD_continuous <- function(n_sub){

	# create the empirical SAD
	df_emp <- calccSADdata(n_sub)

	# get the theoretical sSAD
	df_the <- read.csv(paste0('PCF/Output_',n_sub,'/cSAD_',n_sub,'_continuous.csv'),header=TRUE,row.names=1) %>%
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

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/cSAD_',n_sub,'_continuous.png'),.,device='png',width=15,height=10,units='cm',bg='white')

}


# plot the empirical and theoretical disconnected SAR
# requires: n_sub = the number of samples in the subset used as starting point
plotdSAR_continuous <- function(n_sub){

	# create the empirical disconnected SAR
	dSAR_emp <- read.csv('BCI_dSAR_emp.csv',header=TRUE,row.names=1)

	# get the theoretical disconnected SAR
	dSAR_the <- read.csv(paste0('PCF/Output_',n_sub,'/dSAR_',n_sub,'_continuous.csv'),header=TRUE,row.names=1)

	# create dataframe for the theory errorbars
	dSAR_pol <- data.frame(R=c(dSAR_the$R,rev(dSAR_the$R)),rib=c(dSAR_the$Smean-3*dSAR_the$Ssd,rev(dSAR_the$Smean+3*dSAR_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=dSAR_pol,aes(x=R,y=rib),fill='red',alpha=0.4) +
					geom_point(data=dSAR_emp,aes(x=R,y=S)) + 
					geom_line(data=dSAR_the,aes(x=R,y=Smean)) + 
					theme_bw() + xlab('radius (m)') + ylab('number of species')

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/dSAR_',n_sub,'_continuous.png'),.,device='png',width=15,height=10,units='cm')

}


# plot the empirical and theoretical connected SAR
# requires: n_sub = the number of samples in the subset used as starting point
plotcSAR_continuous <- function(n_sub){

	# create the empirical connected SAR
	cSAR_emp <- read.csv('BCI_cSAR_emp.csv',header=TRUE,row.names=1)

	# get the theoretical connected SAR
	cSAR_the <- read.csv(paste0('PCF/Output_',n_sub,'/cSAR_',n_sub,'_continuous.csv'),header=TRUE,row.names=1)

	# create dataframe for the theory errorbars
	cSAR_pol <- data.frame(R=c(cSAR_the$R,rev(cSAR_the$R)),rib=c(cSAR_the$Smean-3*cSAR_the$Ssd,rev(cSAR_the$Smean+3*cSAR_the$Ssd)))

	# make the plot
	p <- ggplot() + geom_polygon(data=cSAR_pol,aes(x=R,y=rib),fill='red',alpha=0.4) +
					geom_point(data=cSAR_emp,aes(x=R,y=S)) + 
					geom_line(data=cSAR_the,aes(x=R,y=Smean)) + 
					theme_bw() + xlab('radius (m)') + ylab('number of species')

	p %>% ggsave(paste0('PCF/Output_',n_sub,'/plots/cSAR_',n_sub,'_continuous.png'),.,device='png',width=15,height=10,units='cm')

}

# make the plots
# get n_sub = the number of samples in the subsample to start with
PAR_file <- file('PCF/PARAM_file.txt',open='r')
on.exit(close(PAR_file))
PAR_file_lines <- readLines(PAR_file)
n_sub <- strtoi(PAR_file_lines)

if(n_sub == 800){
	n_rep <- 1
} else {
	n_rep <- 100
} 

# create the directory for the plots
dir.create(paste0('PCF/Output_',n_sub,'/plots'))

plotsSAD_discrete(n_sub,n_rep)
plotSAD0_discrete(n_sub)
plotcSAD_discrete(n_sub)
plotdSAR_discrete(n_sub)
plotcSAR_discrete(n_sub)

plotcSAD_continuous(n_sub)
plotSAD0_continuous(n_sub)
plotdSAR_continuous(n_sub)
plotcSAR_continuous(n_sub)

unlink('PCF/PARAM_file.txt')