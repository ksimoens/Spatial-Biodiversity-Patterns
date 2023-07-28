# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Calculatons in the Analytical Models
#####################################
# TRANSFORM the route coordinates
#####################################

# ------------------ LOAD PACKAGES -------------------------

library(tidyverse)
library(sf)

# ----------------------------------------------------------

routes <- read.csv('NABB_routes_div_2021.csv',header=TRUE,row.names=1) 
div <- routes %>% dplyr::select(-c(Longitude,Latitude))

coords <- routes %>% dplyr::select(c(Longitude,Latitude)) %>% 
						sf::st_as_sf(coords=c('Longitude','Latitude'),crs=4326) %>%
						sf::st_transform(crs=5070) %>% sf::st_coordinates() %>%
						as.data.frame() %>% dplyr::rename(x=X,y=Y)

routes <- cbind(coords,div)
write.csv(routes,'NABB_routes_div_2021_proj.csv')