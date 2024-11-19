#Written by Cindy Norton 2022
#Script segments a point cloud into clusters using watershed segmentation on rasterized point cloud
#install.packages('VoxR')
#remotes::install_github('bi0m3trics/spanner')
#install.packages("remotes")
#install.packages(c('spanner','readxl',"dplR","dplyr","tidyr", "geosphere", "lidR", "raster", "TreeLS", "rgdal", "rgeos", "sp","sf","tibble","ggplot2","tidyverse"),dependencies=TRUE)
#remotes::install_github("Jean-Romain/lidR", dependencies=TRUE)
#remotes::install_github('tiagodc/TreeLS',dependencies=TRUE)
#remotes::install_github('cszang/treeclim',dependencies=TRUE)
#install.packages("janitor")

#Packages <- c("dplyr", "ggplot2",  "httr", "rjson", "splitstackshape", "jsonlite", "curl","purrr","reshape2","tidyr","stringr","broom","modelr","lubridate")
#lapply(Packages, library, character.only = TRUE)
#library("pracma")
Packages <- c('cowplot','matrixStats','patchwork','treeclim','remotes','readxl',"dplR","dplyr","tidyr", "geosphere", "lidR", 'terra',"sf","tibble","ggplot2","tidyverse","lubridate","janitor")
lapply(Packages, library, character.only = TRUE)


setwd("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/")

options(scipen=999)

############################################ Tree Ring and Ground Fusion ####################################################
#SET the path for folder where you have the file stoblue, each function will load the data format
groundA =   readxl::read_excel('//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=1)
groundC =   readxl::read_excel('//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=2)
groundB =   readxl::read_excel('//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=3)
MBA<- "//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Treerings/MtBigelow_Treerings/MBA_revised.rwl"
MBB<- "//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Treerings/MtBigelow_Treerings/MBB.rwl"
MBC<- "//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Treerings/MtBigelow_Treerings/MBC.rwl"

mtB_A_xyz <- "//snow/projects/Babst_Lidar_treering_CLNF/CLNorton/outputs/TLS_Bigelow_A_PCGR.xyz"
mtB_C_xyz <- "//snow/projects/Babst_Lidar_treering_CLN/CLNorton/outputs/TLS_Bigelow_C_GRPC.xyz"

drone <-  "//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/outputs/Mt.Bigelow_lidarFlight/Mt_Bigelow_2022_Final.las"
climate <- "//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/outputs/Mt.bigelow_prism_metereological.xlsx"
#noaa_climate_csv <- read.csv("https://www.ncei.noaa.gov/access/services/data/v1?dataset=daily-summaries&dataTypes=PRCP,TMAX,TMIN&stations=USC00025737,USC00025734,USC00025735,USC00025733,USC00025732,USC00026202&startDate=1950-01-01&endDate=2022-12-19&format=csv&options=includeStationName:1&includeStationName=true&units=metric")
#rwl<-dplR::read.rwl(tree) #read file path of rwl


###FUNCTION 1 - tree ring data parse and formatting###
#raw tree ring data parsing and mean tree cores
treering_parse <- function(tree) {
  rwl<-dplR::read.rwl(tree) #read file path of rwl
  detrend.rwi <- dplR::detrend(rwl = rwl, method = "Spline") %>% #detrend rwl
    t() %>% #transpose years as columns and row as each core sample 
    as.data.frame() %>% #as data frame
    dplyr::rename_with( ~ paste0(.x)) %>% #adding "year_" to each year column name for unique identifiers
    tibble::rownames_to_column()%>% #converts sample row name into a column to be parsed
    tidyr::separate(rowname,
                    into = c("site", "coreID"),
                    sep = "(?<=[A-ZA-Z])(?=[0-9])") %>% ## separates column of site and treeID
    tidyr::separate(coreID,
                    into = c("treeID", "core"),
                    sep = "(?<=[0-9])(?=[A-ZA-Za-z])")
  return(detrend.rwi)}
MBA_parse <- treering_parse(MBA)
MBB_parse <- treering_parse(MBB)
MBC_parse <- treering_parse(MBC)

#write.csv(data.transpose, file = "datatranspose.csv")
















##RAW GROUND DATA
MBA_coor <-c(-110.72624,32.41585) 	
MBB_coor <-c(-110.72578,32.41596) 	
MBC_coor <-c(-110.72517,32.41548) 	


###FUNCTION 2 - ground data azimuth and distance to coordinates###
ground_parse <- function(coor, ground) { #read xlsx ground metadata file
  dist_azim <- ground
  coor_list <- list() #creates empty list for new coordinates
  
  for (i in 1:nrow(dist_azim)) { #for loop to create a lat long for each azimuth and distance from plot center lat and lon
    stemAzimuth <- dist_azim$Azimuth[i] #grab azimuth
    stemDistance <- dist_azim$Distance[i] #grab distance
    
    new_coordinates <- geosphere::destPoint(coor, stemAzimuth,  stemDistance) %>% #estimates new lat lon coordinate
      as.data.frame() %>% #make as data frame
      add_column(treeID = dist_azim$treeID[i]) %>% #add treeID column of that tree
      add_column(stemAzimuth = stemAzimuth) %>% #add Azimuth column of that tree
      add_column(stemDistance = stemDistance) %>% #add Distance column of that tree
      add_column(species = dist_azim$Species[i]) %>% #add Species column of that tree
      add_column(Height = dist_azim$Height[i])%>% #add Height column of that tree
      add_column(DBH = dist_azim$DBH[i]) #add DBH column of that tree
    
    
    coor_list[[i]] <- new_coordinates #put result of loop in list
    
  }
  
  coor_df<- coor_list %>% bind_rows() %>% na.omit() #remove NA rows
  
  
  return(coor_df)
}

MBA_ground <- ground_parse(MBA_coor, groundA) 
MBB_ground <- ground_parse(MBB_coor, groundB) 
MBC_ground <- ground_parse(MBC_coor, groundC) 











###FUNCTION 3 -  tree ring and ground data merge###
###TREE rings and ground merge#####
treering_ground <- function (treering_parse, ground_parse, plot) {
  
  treering_ground_join <- treering_parse %>% 
    subset(site == plot) %>% #subsets to input plot
    lapply(as.numeric)%>% #makes numeric
    as.data.frame() %>% #makes into df
    group_by(treeID) %>% #groups data by the treeID
    summarise(across(-c(site, core), mean, na.rm = TRUE))%>% #a mean summary of multiple cores per tree
    left_join(ground_parse, by = "treeID") 
  #joins tree ring data with ground metadat using treeID
  colnames(treering_ground_join)<-gsub("X","",colnames(treering_ground_join))

  return(treering_ground_join)
}

MBA_treering_ground <- treering_ground(MBA_parse, MBA_ground, "MBA")
MBB_treering_ground <- treering_ground(MBB_parse, MBB_ground, "MBB")
MBC_treering_ground <- treering_ground(MBC_parse, MBC_ground, "MBC")


























###FUNCTION 4 - DBH and Percentage reconstruction###
###DBH RECONSTRUCTIONS
DBH_recon <- function(tree,ground_parse){
  
  rwl<-dplR::read.rwl(tree)%>%#read file path of rwl
    t() %>% #transpose years as columns and row as each core sample 
    as.data.frame() %>% #as data frame
    dplyr::rename_with( ~ paste0("year_",.x)) %>% #adding "year_" to each year column name for unique identifiers
    tibble::rownames_to_column()%>% #converts sample row name into a column to be parsed
    tidyr::separate(rowname,
                    into = c("site", "coreID"),
                    sep = "(?<=[A-ZA-Z])(?=[0-9])") %>% ## separates column of site and treeID
    tidyr::separate(coreID,
                    into = c("treeID", "core"),
                    sep = "(?<=[0-9])(?=[A-ZA-Za-z])")%>%
    rev()
  
  treering_tree <- rwl %>% 
    lapply(as.numeric)%>% #makes numeric
    as.data.frame() %>% #makes into df
    group_by(treeID) %>% #groups data by the treeID
    summarise(across(-c(site, core), mean, na.rm = TRUE))%>%
    left_join(ground_parse, by = "treeID") %>%
    as.data.frame()
  
  treering_tree[is.na(treering_tree)] <- 0  #zero all NA
  year_cols <- which(substr(x = colnames(treering_tree), start = 1, stop = 4) == "year")  # Identify the columns that have data of interest which are year columns
  treering_tree$sum <- rowSums(x = treering_tree[, year_cols], na.rm = TRUE)  #row sums all the years into column, removing NA 
  csum <- treering_tree  # Make a copy of that data frame that we will use for difference
  csum[, year_cols] <- sapply(1:ncol(csum[, year_cols]), function(col){rowSums(csum[, year_cols][1:col], na.rm = TRUE)}) #convert years to cumulative sum of year columns
  dif <- csum # Make a copy of that data frame that we will use for difference
  dif[, year_cols] <- (dif$sum) - (dif[, year_cols]) #calculate difference
  frac <- csum # Make a copy of that data frame that we will use for difference
  frac[, year_cols] <- dif[, year_cols]/dif$sum #calculate difference
  dbh <- csum # Make a copy of that data frame that we will use for difference
  dbh[, year_cols] <- as.numeric(frac$DBH)*frac[, year_cols] #calculate difference

  # Reality check to see that percentages add up to 100
  
  return(dbh)
}



MBA_DBH_recon_results <- DBH_recon(MBA, MBA_ground)
MBB_DBH_recon_results <- DBH_recon(MBB, MBB_ground)
MBC_DBH_recon_results <- DBH_recon(MBC, MBC_ground)


###DBH RECONSTRUCITONS


####MERGE PLOTS AND SUBSET BY SPECIES### format for treeclim ###
MBA_ponderosa_growth = MBA_DBH_recon_results %>% filter(species == "ponderosa") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBA"))%>%
  rename(year = 1)
MBB_ponderosa_growth = MBB_DBH_recon_results %>% filter(species == "ponderosa") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBB"))%>%
  rename(year = 1)
MBC_ponderosa_growth = MBC_DBH_recon_results %>% filter(species == "ponderosa") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBC"))%>%
  rename(year = 1)


MB_all_ponderosa_growth = MBB_ponderosa_growth %>% left_join(MBA_ponderosa_growth, by = "year") %>% 
  left_join(MBC_ponderosa_growth, by = "year")%>%
  slice(-c(217:224))%>%
  mutate_at(c(2:46), as.numeric)%>%
  mutate(ponderosa_mean = rowMeans(across(2:46), na.rm = TRUE))%>%
  apply(2,rev) %>%
  as.data.frame()%>%
  mutate_at(-1,as.numeric) %>%
  mutate_if(is.numeric, ~na_if(., 0))
  
#write.csv(MB_all_ponderosa_growth, file = "MB_all_ponderosa_growth_100323.csv")

#windows()
matplot(rownames(MB_all_ponderosa_growth), MB_all_ponderosa_growth, type='l', xlab='year', ylab='DBH', col=1:5)



MBA_menziesii_growth = MBA_DBH_recon_results %>% filter(species == "menziesii") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBA"))%>%
  rename(year = 1)
MBB_menziesii_growth = MBB_DBH_recon_results %>% filter(species == "menziesii") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBB"))%>%
  rename(year = 1)
MBC_menziesii_growth = MBC_DBH_recon_results %>% filter(species == "menziesii") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBC"))%>%
  rename(year = 1)



MB_all_menziesii_growth = MBB_menziesii_growth %>% left_join(MBA_menziesii_growth, by = "year") %>% 
  left_join(MBC_menziesii_growth, by = "year") %>%
  slice(-c(217:224))%>%                                                                 
  mutate_at(c(2:17), as.numeric) %>%
  mutate(menziesii_mean = rowMeans(across(2:17), na.rm = TRUE))%>%
  apply(2,rev) %>%
  as.data.frame()%>%
  mutate_at(-1,as.numeric) %>%
  mutate_if(is.numeric, ~na_if(., 0)) 

#write.csv(MB_all_menziesii_growth, file = "MB_all_menziesii_growth_100323.csv")

windows()
matplot(rownames(MB_all_menziesii_growth), MB_all_menziesii_growth, type='l', xlab='year', ylab='DBH', col=1:5)


MBA_strobiformis_growth = MBA_DBH_recon_results %>% filter(species == "strobiformis") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBA"))%>%
  rename(year = 1)
MBB_strobiformis_growth = MBB_DBH_recon_results %>% filter(species == "strobiformis") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBB"))%>%
  rename(year = 1)
MBC_strobiformis_growth = MBC_DBH_recon_results %>% filter(species == "strobiformis") %>% t() %>% as.data.frame() %>% rownames_to_column() %>% rename(year = rowname)%>%
  row_to_names(1) %>%
  rename_at(vars(-treeID),function(x) paste0(x,"_MBC"))%>%
  rename(year = 1)


MB_all_strobiformis_growth = MBB_strobiformis_growth %>% left_join(MBA_strobiformis_growth, by = "year") %>% 
  left_join(MBC_strobiformis_growth, by = "year") %>%
  slice(-c(217:224)) %>%
  mutate_at(-c(1), as.numeric) %>%
  mutate(strobiformis_mean = rowMeans(across(2:42), na.rm = TRUE))%>%
  apply(2,rev) %>%
  as.data.frame()%>%
  mutate_at(-1,as.numeric)%>% 
  mutate_if(is.numeric, ~na_if(., 0))
  

#write.csv(MB_all_strobiformis_growth, file = "MB_all_strobiformis_growth_100323.csv")


windows()
matplot(rownames(MB_all_strobiformis_growth), MB_all_strobiformis_growth, type='l', xlab='year', ylab='DBH', col=1:5)


MB_all_growth_recon <- MB_all_ponderosa_growth %>% 
  left_join(MB_all_menziesii_growth ,by = "year") %>% 
  left_join(MB_all_strobiformis_growth ,by = "year") %>% na.omit() %>%
  t()%>%
  as.data.frame()%>%
  purrr::set_names(as.character(slice(., 1))) %>%
  slice(-1) %>%
  mutate_all(as.numeric) %>% rev() %>%t()%>%
  as.data.frame()%>%
  slice(-c(1:7)) 
#write.csv(MB_all_growth_recon,'MB_all_growth_recon.csv')

windows()
plot(MB_all_ponderosa_growth$ponderosa_mean, type = 'l', lwd = 2, col = 'red')
lines(MB_all_menziesii_growth$menziesii_mean, lwd = 2, col = 'green')
lines(MB_all_strobiformis_growth$strobiformis_mean, lwd = 2, col = 'blue')
legend("bottomright",legend=c("PIPO","PM","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n", lwd = 2)

#P. ponderosa 	Pine, ponderosa 	122 	0.38 	9 
#Pseudotsuga menziesii 	Douglas-fir 	202 	0.45 	23 
#P. strobus 	Pine, eastern white 	129 	0.34 	14 

#Chojnacky et al. 2014 blanket group
#Group	Taxa	Median specific gravity	β0	β1	
#Conifer	Pseudotsuga	0.45	−2.4623	2.4852	
#Conifer	Pinus ≥ 0.45 spg	0.47	−30.1006	2.6465
#Conifer	Pinus < 0.45 spg	0.39	−2.6177	2.4638	

##PONDEROSA
#Callaway et al. 1994 desert
#LM = -1.398 + 1.805(log dbh),
#ln(biomass) = β0 + β1 ln(dbh)

#kaye et al. 2005 flagstaff
#biomass = a × e^[b+ln(dbh)×c]
#1.0469e^(-4.1279+ln(DBH)2.2391)
#write.csv(MB_all,'MB_all.csv')
##biomass
#years
#3:48, 3:40

MB_all_PIPO_biomass = read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_ponderosa_growth_100323_subset.csv")%>%
  mutate(across(3:40, ~ exp((-2.6177+2.4638 * log(.)))))%>%
  slice(-(1:144))%>%
  mutate(PIPO_biomass = rowMeans(dplyr::select(., 3:40),na.rm=TRUE), 
         PIPO_biomass_sdev = rowSds(as.matrix(.[3:40]),na.rm=TRUE))%>%
  mutate(across(3:40, ~ (((.))-lag(.))))%>%
  mutate(PIPO_increment = rowMeans(dplyr::select(., 3:40),na.rm=TRUE), 
         PIPO_increment_sdev = rowSds(as.matrix(.[3:40]),na.rm=TRUE))
#3:43
MB_all_PISF_biomass =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_strobiformis_growth_100323_subset.csv")%>%
  mutate(across(3:35, ~ exp((-2.6177+2.4638 * log(.)))))%>%
  slice(-(1:144))%>%
  mutate(PISF_biomass = rowMeans(dplyr::select(., 3:35),na.rm=TRUE), 
         PISF_biomass_sdev = rowSds(as.matrix(.[3:35]),na.rm=TRUE))%>%
  mutate(across(3:35, ~ (((.))-lag(.))))%>%
  mutate(PISF_increment = rowMeans(dplyr::select(., 3:35),na.rm=TRUE), 
         PISF_increment_sdev = rowSds(as.matrix(.[3:35]),na.rm=TRUE))
#3:18
MB_all_PM_biomass =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_menziesii_growth_100323_subset.csv")%>%
  mutate(across(3:15, ~ exp((-2.4623+2.4852 * log(.)))))%>%
  slice(-(1:144))%>%
  mutate(PM_biomass = rowMeans(dplyr::select(., 3:15),na.rm=TRUE), 
         PM_biomass_sdev = rowSds(as.matrix(.[3:15]),na.rm=TRUE))%>%
  mutate(across(3:15, ~ (((.))-lag(.))))%>%
  mutate(PM_increment = rowMeans(dplyr::select(., 3:15),na.rm=TRUE), 
         PM_increment_sdev = rowSds(as.matrix(.[3:15]),na.rm=TRUE))

##basal increment
#years
#MB_all_PIPO_basal = read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_ponderosa_growth_100323_subset.csv")%>% 
#  mutate(across(3:48, ~ pi*((./2)^2)))%>%
#  #slice(-(1:144))%>%
#  mutate(PIPO_basal = rowMeans(dplyr::select(., 3:48),na.rm=TRUE), 
#         PIPO_basal_sdev = rowSds(as.matrix(.[3:48]),na.rm=TRUE))%>%
#  mutate(across(3:48, ~ (((.))-lag(.))))%>%
#  mutate(PIPO_increment = rowMeans(dplyr::select(., 3:48),na.rm=TRUE), 
#         PIPO_increment_sdev = rowSds(as.matrix(.[3:48]),na.rm=TRUE))

#MB_all_PISF_basal =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_strobiformis_growth_100323_subset.csv")%>%
#  mutate(across(3:43, ~ pi*((./2)^2)))%>%
#  #slice(-(1:144))%>%
#  mutate(PISF_basal = rowMeans(dplyr::select(., 3:43),na.rm=TRUE), 
#         PISF_basal_sdev = rowSds(as.matrix(.[3:43]),na.rm=TRUE))%>%
#  mutate(across(3:43, ~ (((.))-lag(.))))%>%
#  mutate(PISF_increment = rowMeans(dplyr::select(., 3:43),na.rm=TRUE), 
#         PISF_increment_sdev = rowSds(as.matrix(.[3:43]),na.rm=TRUE))

#MB_all_PM_basal =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_menziesii_growth_100323.csv")%>%
#  mutate(across(3:18, ~ pi*((./2)^2)))%>%
#  #slice(-(1:144))%>%
#  mutate(PM_basal = rowMeans(dplyr::select(., 3:18),na.rm=TRUE), 
#         
#         PM_basal_sdev = rowSds(as.matrix(.[3:18]),na.rm=TRUE))%>%
#  mutate(across(3:18, ~ (((.))-lag(.))))%>%
#  mutate(PM_increment = rowMeans(dplyr::select(., 3:18),na.rm=TRUE), PM_increment_sdev = rowSds(as.matrix(.[3:18]),na.rm=TRUE))

#AGE
#BIOMASS
#age
MB_all_PIPO_biomass_age = read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_ponderosa_growth_100323_narm.csv")%>% 
  mutate(ponderosa_biomass = exp(-2.6177+2.4638 * log(as.numeric(PIPO))))%>%
  mutate(across(3:48, ~ exp((-2.6177+2.4638 * log(.)))))%>%
  #slice(-(50:216))%>%
  mutate(PIPO_biomass_age = rowMeans(dplyr::select(., 3:48),na.rm=TRUE), PIPO_biomass_age_sdev = rowSds(as.matrix(.[3:48]),na.rm=TRUE))%>%
  mutate(across(3:48, ~ (((.))-lag(.))))%>%
  mutate(PIPO_increment_age = rowMeans(dplyr::select(., 3:48),na.rm=TRUE), PIPO_increment_sdev_age = rowSds(as.matrix(.[3:48]),na.rm=TRUE))

MB_all_PISF_biomass_age =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_strobiformis_growth_100323_narm.csv")%>% 
  mutate(strobiformis_biomass = exp(-2.6177+2.4638 * log(as.numeric(PISF)))) %>%
  mutate(across(3:43, ~ exp((-2.6177+2.4638 * log(.)))))%>%  
  #slice(-(50:216))%>%
  mutate(PISF_biomass_age = rowMeans(dplyr::select(., 3:43),na.rm=TRUE), PISF_biomass_age_sdev = rowSds(as.matrix(.[3:43]),na.rm=TRUE))%>%
  mutate(across(3:43, ~ (((.))-lag(.))))%>%
  mutate(PISF_increment_age = rowMeans(dplyr::select(., 3:43),na.rm=TRUE), PISF_increment_sdev_age = rowSds(as.matrix(.[3:43]),na.rm=TRUE))

MB_all_PM_biomass_age =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_menziesii_growth_100323_narm.csv")%>% 
  mutate(menziesii_biomass = exp(-2.4623 + 2.4852 * log(as.numeric(PSME)))) %>%
  mutate(across(3:18, ~ exp((-2.4623+2.4852 * log(.)))))%>%
  #slice(-(50:216))%>%
  mutate(PM_biomass_age = rowMeans(dplyr::select(., 3:18),na.rm=TRUE), PM_biomass_age_sdev = rowSds(as.matrix(.[3:18]),na.rm=TRUE))%>%
  mutate(across(3:18, ~ (((.))-lag(.))))%>%
  mutate(PM_increment_age = rowMeans(dplyr::select(., 3:18),na.rm=TRUE), PM_increment_sdev_age = rowSds(as.matrix(.[3:18]),na.rm=TRUE))
#BASAL AREA
#age
#MB_all_PIPO_basal_age = read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_ponderosa_growth_100323_narm.csv")%>% 
#  mutate(across(3:48, ~ pi*((./2)^2)))%>%
#  #slice(-(50:216))%>%
#  mutate(PIPO_basal_age = rowMeans(dplyr::select(., 3:48),na.rm=TRUE), PIPO_basal_age_sdev = rowSds(as.matrix(.[3:48]),na.rm=TRUE))%>%
#  mutate(across(3:48, ~ (((.))-lag(.))))%>%
#  mutate(PIPO_increment_age = rowMeans(dplyr::select(., 3:48),na.rm=TRUE), PIPO_increment_sdev_age = rowSds(as.matrix(.[3:48]),na.rm=TRUE))

#MB_all_PISF_basal_age =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_strobiformis_growth_100323_narm.#csv")%>%
#  mutate(across(3:43, ~ pi*((./2)^2)))%>%
#  #slice(-(50:216))%>%
#  mutate(PISF_basal_age = rowMeans(dplyr::select(., 3:43),na.rm=TRUE), PISF_basal_age_sdev = rowSds(as.matrix(.[3:43]),na.rm=TRUE))%>%
#  mutate(across(3:43, ~ (((.))-lag(.))))%>%
#  mutate(PISF_increment_age = rowMeans(dplyr::select(., 3:43),na.rm=TRUE), PISF_increment_sdev_age = rowSds(as.matrix(.[3:43]),na.rm=TRUE))

#MB_all_PM_basal_age =  read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_menziesii_growth_100323_narm.csv")%>%
#  mutate(across(3:18, ~ pi*((./2)^2)))%>%
#  #slice(-(50:216))%>%
#  mutate(PM_basal_age = rowMeans(dplyr::select(., 3:18),na.rm=TRUE), PM_basal_age_sdev = rowSds(as.matrix(.[3:18]),na.rm=TRUE))%>%
#  mutate(across(3:18, ~ (((.))-lag(.))))%>%
#  mutate(PM_increment_age = rowMeans(dplyr::select(., 3:18),na.rm=TRUE), PM_increment_sdev_age = rowSds(as.matrix(.[3:18]),na.rm=TRUE))



MB_all_biomass_years = cbind(MB_all_PIPO_biomass,
                    MB_all_PISF_biomass,
                    MB_all_PM_biomass) %>% as.data.frame()
#MB_all_basal_years = cbind(MB_all_PIPO_basal,
#                            MB_all_PISF_basal,

MB_all_biomass_age <- MB_all_PIPO_biomass_age %>%
  select(ponderosa_biomass = ponderosa_biomass) %>%
  mutate(strobiformis_biomass = MB_all_PISF_biomass_age$strobiformis_biomass,
         menziesii_biomass = MB_all_PM_biomass_age$menziesii_biomass)%>% rownames_to_column()%>%
  filter(as.numeric(rowname) <= 75)


#MB_all_basal_age = cbind(MB_all_PIPO_basal_age,
#                           MB_all_PISF_basal_age,
#                           MB_all_PM_basal_age) %>% as.data.frame()


#basal_PIPO = mean(MB_all_basal_years$PIPO_increment, na.rm =TRUE)
#basal_PISF = mean(MB_all_basal_years$PISF_increment, na.rm =TRUE)
#basal_PM = mean(MB_all_basal_years$PM_increment, na.rm =TRUE)
biomass_PIPO = mean(MB_all_biomass_years$PIPO_increment, na.rm=TRUE)
biomass_PISF = mean(MB_all_biomass_years$PISF_increment, na.rm=TRUE)
biomass_PM = mean(MB_all_biomass_years$PM_increment, na.rm=TRUE)

#basal_PIPO_age = mean(MB_all_basal_age$PIPO_increment_age, na.rm =TRUE)
#basal_PISF_age = mean(MB_all_basal_age$PISF_increment_age, na.rm =TRUE)
#basal_PM_age = mean(MB_all_basal_age$PM_increment_age, na.rm =TRUE)
biomass_PIPO_age = mean(MB_all_biomass_age$ponderosa_biomass, na.rm=TRUE)
biomass_PISF_age = mean(MB_all_biomass_age$strobiformis_biomass, na.rm=TRUE)
biomass_PM_age = mean(MB_all_biomass_age$menziesii_biomass, na.rm=TRUE)

#write.csv(MB_all_biomass_years,'MB_all_100323_biomass_years.csv')
#write.csv(MB_all_basal_years,'MB_all_100323_basal_years.csv')
#write.csv(MB_all_biomass_age,'MB_all_06242024_biomass_age.csv')
#write.csv(MB_all_basal_age,'MB_all_100323_basal_age.csv')


#ba <- ((diam/2)^2)*pi


#P. ponderosa 	Pine, ponderosa 	122 	0.38 	9 
#Pseudotsuga menziesii 	Douglas-fir 	202 	0.45 	23 
#P. strobus 	Pine, eastern white 	129 	0.34 	14 

#Chojnacky et al. 2014 blanket group
#Group	Taxa	Median specific gravity	β0	β1	
#Conifer	Pseudotsuga	0.45	−2.4623	2.4852	
#Conifer	Pinus ≥ 0.45 spg	0.47	−30.1006	2.6465
#Conifer	Pinus < 0.45 spg	0.39	−2.6177	2.4638	

##PIPO
#Callaway et al. 1994 desert
#LM = -1.398 + 1.805(log dbh),
#ln(biomass) = β0 + β1 ln(dbh)

#kaye et al. 2005 flagstaff
#biomass = a × e^[b+ln(dbh)×c]
#1.0469e^(-4.1279+ln(DBH)2.2391)

#PISF Hernandez et al . 2020
#Chapman-Richards	y=a(1−(e^−bt))^c)+e
#Schumacher	y=ae^(−b(1t)))+e


#PM
#write.csv(MB_all,'MB_all.csv')

MB_all_growth_years_melt = read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_years_subsetmelt.csv")
MB_all_growth_years = read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_100323_increment_update_subset.csv")
MB_all_basal_years = read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/MB_all_100323_basal_years.csv")


windows()
plot(MB_all_biomass_age$ponderosa_biomass, type = 'l', lwd = 2, col = 'red', ylab = "Biomass", xlab = "Age") #xlim=c(0,80), ylim = c(0,100))
lines(MB_all_biomass_age$strobiformis_biomass, lwd = 2, col = 'blue')
lines(MB_all_biomass_age$menziesii_biomass, lwd = 2, col = 'green')
legend("topleft",legend=c("PIPO","PISF","PSME"),lty=1, col = c('red','blue','green'),cex = 1,bty="n", lwd = 2)

MB_all_biomass_age_long <- MB_all_biomass_age %>%
  pivot_longer(cols = c(ponderosa_biomass, strobiformis_biomass, menziesii_biomass),
               names_to = "species", values_to = "biomass")


# Set custom colors and labels for species
species_colors <- c("ponderosa_biomass" = 'red', "strobiformis_biomass" = 'blue', "menziesii_biomass" = 'green')
species_labels <- c("ponderosa_biomass" = "PIPO", "strobiformis_biomass" = "PISF", "menziesii_biomass" = "PSME")

windows()
ggplot(MB_all_biomass_age_long, aes(x = as.numeric(rowname), y = biomass, color = species)) +
  geom_smooth(method = "loess", se = FALSE, linetype = "solid", aes(color = species)) +  # Loess smoothing
  scale_color_manual(values = species_colors, labels = species_labels) +  # Match the base plot colors and labels
  labs(y = "Biomass", x = "Age", color = "Species") +
  theme_minimal() +
  theme(
    legend.position = "left",  # Position the legend on the left
    axis.title = element_text(size = 14, face = "bold"),  # Bold and increase axis titles
    axis.text = element_text(size = 12, face = "bold"),   # Bold and increase axis text
    axis.ticks = element_line(size = 1)  # Make axis ticks thicker
  )

# Calculate growth rate for each species in MB_all_biomass_age_long
MB_all_biomass_growthrate_long <- MB_all_biomass_age_long %>%
  group_by(species) %>%
  mutate(
    growth_rate = (biomass - lag(biomass)) / (as.numeric(rowname) - lag(as.numeric(rowname)))
  )
library("broom")

# Calculate slope for each species using linear regression
slope_table <- MB_all_biomass_growthrate_long %>%
  group_by(species) %>%
  do(tidy(lm(biomass ~ as.numeric(rowname), data = .))) %>%
  filter(term == "as.numeric(rowname)") %>%
  select(species, slope = estimate) %>%
  mutate(slope = round(slope, 2))

# Summarize growth rate for each species (mean growth rate)
growth_rate_table <- MB_all_biomass_growthrate_long %>%
  group_by(species) %>%
  summarize(mean_growth_rate = round(mean(growth_rate, na.rm = TRUE), 2))

#windows()
#plot(biomass_PIPO$PIPO_biomass, type = 'l', lwd = 2, col = 'red')
#lines(biomass_PISF$PISF_biomass, lwd = 2, col = 'blue')
#lines(biomass_PSME$PM_biomass, lwd = 2, col = 'green')
#legend("topright",legend=c("PIPO","PISF","PSME"),lty=1, col = c('red','blue','green'),cex = 1,bty="n", lwd = 2)

#basal_PIPO = mean(MB_all_basal_years$PIPO_increment, na.rm =TRUE)
#basal_PISF = mean(MB_all_basal_years$PISF_increment, na.rm =TRUE)
#basal_PM = mean(MB_all_basal_years$PM_increment, na.rm =TRUE)
biomass_PIPO = mean(MB_all_biomass_years$PIPO_increment, na.rm=TRUE)
biomass_PISF = mean(MB_all_biomass_years$PISF_increment, na.rm=TRUE)
biomass_PM = mean(MB_all_biomass_years$PM_increment, na.rm=TRUE)

#basal_PIPO_age = mean(MB_all_basal_age$PIPO_increment_age, na.rm =TRUE)
#basal_PISF_age = mean(MB_all_basal_age$PISF_increment_age, na.rm =TRUE)
#basal_PM_age = mean(MB_all_basal_age$PM_increment_age, na.rm =TRUE)
biomass_PIPO_age = mean(MB_all_biomass_age$PIPO_increment_age, na.rm=TRUE)
biomass_PISF_age = mean(MB_all_biomass_age$PISF_increment_age, na.rm=TRUE)
biomass_PM_age = mean(MB_all_biomass_age$PM_increment_age, na.rm=TRUE)

# Make the plot
#ggplot(data=MB_all_growth_years_melt, aes(x=year, y=biomass_allometry_ground, ymin=mean - stdev, ymax=mean + stdev, fill=species, linetype=species)) + 
#  geom_line() + 
#  geom_ribbon(alpha=0.3) +
#  xlab(as.expression(expression( paste("Radius (", R[500], ")") ))) + 
#  ylab("DBH Reconstruction")

windows()
plot(MB_all_growth_years$PIPO_increment, type = 'l', lwd = 2, col = 'red', xlim=c(0,280), ylim = c(0,25))
lines(MB_all_growth_years$PISF_increment, lwd = 2, col = 'blue')
lines(MB_all_growth_years$PM_increment, lwd = 2, col = 'green')
#abline(0, 0.27, col = 'red', lty = "dashed",lwd = 1.5)
legend("topright",legend=c("PIPO","PM","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n", lwd = 2)

# Add shaded regions for standard deviations
#polygon(c(MB_all_growth_years$PIPO_sdev, rev(MB_all_growth_years$PIPO_sdev)), c(MB_all_growth_years$PIPO_mean - MB_all_growth_years$PIPO_sdev, rev(MB_all_growth_years$PIPO_mean + MB_all_growth_years$PIPO_sdev)), col = alpha('green', 0.3), border = NA)
#polygon(c(x, rev(x)), c(MB_all_growth_years$PISF_mean - MB_all_growth_years$PISF_sd, rev(MB_all_growth_years$PISF_mean + MB_all_growth_years$PISF_sd)), col = alpha('green', 0.3), border = NA)
#polygon(c(x, rev(x)), c(MB_all_growth_years$PM_mean - MB_all_growth_years$PM_sd, rev(MB_all_growth_years$PM_mean + MB_all_growth_years$PM_sd)), col = alpha('blue', 0.3), border = NA)

mod1 = lm(MB_all_growth_years$PIPO_increment~MB_all_growth_years$X, rm.na=TRUE)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(MB_all_growth_years$PISF_increment~MB_all_growth_years$X)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = "green")


mod1 = lm(MB_all_growth_years$PM_increment~MB_all_growth_years$X)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'green')







# Identify the columns that have data of interest


##PERCENTAGE DBH####
DBH_perc <- function(diameter){
  year_cols <- which(substr(x = colnames(diameter), start = 1, stop = 4) == "year")# Identify the columns that have data of interest
  DBH_perc <- diameter # Make a copy of that data frame that we will use for percentages
  DBH_perc[, year_cols]<- (DBH_perc[, year_cols]/as.numeric(DBH_perc$DBH))* 100 # Do percentage calculations for only those columns with data of interest
  return(DBH_perc)
}


#DBH_perc_results <- DBH_perc(DBH_recon_results)
##PERCENTAGE DBH####

  
#ba <- ((diam/2)^2)*pi
#bai <- ba
#for(i in 2:nrow(ba)){bai[i,] <- ba[i,]-ba[i-1,]}












####MERGE PLOTS AND SUBSET BY SPECIES### format for treeclim ###
MBA_ponderosa = MBA_treering_ground %>% filter(species == "ponderosa") %>% t() %>% as.data.frame() %>% rownames_to_column() 
MBB_ponderosa = MBB_treering_ground %>% filter(species == "ponderosa") %>% t() %>% as.data.frame() %>% rownames_to_column() 
MBC_ponderosa = MBC_treering_ground %>% filter(species == "ponderosa") %>% t() %>% as.data.frame() %>% rownames_to_column() 
MB_all_ponderosa = MBA_ponderosa %>% left_join(MBB_ponderosa, by = "rowname") %>% 
  left_join(MBC_ponderosa, by = "rowname") %>% 
  slice(-c(1,156:162)) %>% 
  `rownames<-`(.[,1]) %>% 
  subset(rowname >= 1930 & rowname < 2021)%>%
  dplyr::select(-rowname) %>%
  mutate_all(function(x) as.numeric(as.character(x))) 


MB_all_ponderosa_detrend.rwi.crn <- chron(MB_all_ponderosa)
#plot(MB_all_ponderosa_detrend.rwi.crn,xlab="Year",ylab="RWI")

MBA_menziesii = MBA_treering_ground %>% filter(species == "menziesii")%>% t() %>% as.data.frame() %>% rownames_to_column() 
MBB_menziesii = MBB_treering_ground %>% filter(species == "menziesii")%>% t() %>% as.data.frame() %>% rownames_to_column() 
MBC_menziesii = MBC_treering_ground %>% filter(species == "menziesii")%>% t() %>% as.data.frame() %>% rownames_to_column() 
MB_all_menziesii = MBA_menziesii %>% left_join(MBB_menziesii, by = "rowname") %>%
  left_join(MBC_menziesii, by = "rowname") %>% 
  slice(-c(1,156:162)) %>% 
  `rownames<-`(.[,1]) %>% 
  subset(rowname >= 1930 & rowname < 2021)%>%
  dplyr::select(-rowname)


MB_all_menziesii_detrend.rwi.crn <- chron(MB_all_menziesii)
#plot(MB_all_menziesii_detrend.rwi.crn,xlab="Year",ylab="RWI")


MBA_strobiformis = MBA_treering_ground %>% filter(species == "strobiformis")%>% t() %>% as.data.frame() %>% rownames_to_column() 
MBB_strobiformis = MBB_treering_ground %>% filter(species == "strobiformis")%>% t() %>% as.data.frame() %>% rownames_to_column() 
MBC_strobiformis = MBC_treering_ground %>% filter(species == "strobiformis")%>% t() %>% as.data.frame() %>% rownames_to_column() 
MB_all_strobiformis = MBA_strobiformis %>% left_join(MBB_strobiformis, by = "rowname") %>% 
  left_join(MBC_strobiformis, by = "rowname") %>% 
  slice(-c(1,156:162)) %>% 
  `rownames<-`(.[,1])%>% 
  subset(rowname >= 1930 & rowname < 2021)%>%
  dplyr::select(-rowname)

MB_all_strobiformis_detrend.rwi.crn <- chron(MB_all_strobiformis)
plot(MB_all_strobiformis_detrend.rwi.crn,xlab="Year",ylab="RWI")
####MERGE PLOTS AND SUBSET BY SPECIES###



























#######MEAN BY species BY PLOT --- mean of mean species by plot######

MBA_mean_group <- MBA_treering_ground %>% as.data.frame() %>% dplyr::select(-c(treeID,lon,lat,stemAzimuth,stemDistance,DBH,Height)) %>% 
group_by(species)  %>% summarise(across(everything(), list(mean), na.rm = TRUE))%>% filter(species == "menziesii" | species == "ponderosa" | species == "strobiformis")

MBA_mean_group_t <- t(MBA_mean_group) %>% as.data.frame() %>% rownames_to_column() 
names(MBA_mean_group_t) <- MBA_mean_group_t %>% slice(1) %>% unlist()
MBA_mean_group_t <- MBA_mean_group_t %>% slice(-1)
MBA_mean_group_t <- MBA_mean_group_t %>% rename(date = species)


MBB_mean_group <- MBB_treering_ground %>% as.data.frame() %>% dplyr::select(-c(treeID,lon,lat,stemAzimuth,stemDistance,DBH,Height)) %>% 
  group_by(species)  %>% summarise(across(everything(), list(mean), na.rm = TRUE))%>% filter(species == "menziesii" | species == "ponderosa" | species == "strobiformis")

MBB_mean_group_t <- t(MBB_mean_group) %>% as.data.frame() %>% rownames_to_column() 
names(MBB_mean_group_t) <- MBB_mean_group_t %>% slice(1) %>% unlist()
MBB_mean_group_t <- MBB_mean_group_t %>% slice(-1)
MBB_mean_group_t <- MBB_mean_group_t %>% rename(date = species)

MBC_mean_group <- MBC_treering_ground %>% as.data.frame() %>% dplyr::select(-c(treeID,lon,lat,stemAzimuth,stemDistance,DBH,Height)) %>%
  group_by(species) %>% summarise(across(everything(), list(mean), na.rm = TRUE)) %>% filter(species == "menziesii" | species == "ponderosa" | species == "strobiformis")

MBC_mean_group_t <- t(MBC_mean_group) %>% as.data.frame() %>% rownames_to_column() 
names(MBC_mean_group_t) <- MBC_mean_group_t %>% slice(1) %>% unlist()
MBC_mean_group_t <- MBC_mean_group_t %>% slice(-1)
MBC_mean_group_t <- MBC_mean_group_t %>% rename(date = species)


MB_all_rbind = rbind(MBA_mean_group_t,MBB_mean_group_t) %>% rbind(MBC_mean_group_t) %>%
  group_by(date)%>%
  mutate_at(c(2:4), as.numeric)%>%
  summarise_at(c("menziesii","ponderosa", "strobiformis"), mean, na.rm=TRUE)%>%
  rename("year" = "date" )%>%
  subset(year > 1901 & year < 2021)

  

MB_all_rbind$year<-gsub("_1","",as.character(MB_all_rbind$year))
MB_all_rbind$year<-gsub("X","",as.character(MB_all_rbind$year))



#######MEAN BY GROUP BY PLOT --- mean of mean group by plot######































#########CLIMATE##########
prism_met <- readxl::read_xlsx(climate)  %>% mutate(date = as.Date(date, format="%Y-%m-%d"))%>%
  mutate(month = month(date), year = year(date)) %>%
  relocate(year, .before = tmean) %>%
  relocate(month, .after = year) %>%
  subset(year >= 1930 & year < 2022) %>%
  rename(prec = ppt, temp = tmean) 


# converts to date format
prism_met$date <- date(prism_met$date)

# add in columns
# add in columns
prism_met <- mutate(prism_met,
                    season = case_when(
                      month(date) %in% c(7,8,9) ~ "prev_mons",
                      month(date) %in% c(10) ~ "prev_fall",
                      month(date) %in% c(11, 12,1, 2,3) ~ "winter",
                      month(date) %in% c(4,5,6) ~ "pre_mons",
                      T ~ NA_character_))
# add in columns
#WATERYEAR
wateryear_breaks <- seq(as.Date("1900-10-01"), length=122, by="year")  
years.wateryear.breaks = as.numeric(format(wateryear_breaks,"%Y"))
labels.wateryear = years.wateryear.breaks[2:length(wateryear_breaks)]  # Why are we leaving off the first year in the water years label?
prism_met$wateryear <- cut(prism_met$date, wateryear_breaks,labels=labels.wateryear)

wateryear <- prism_met %>% group_by(wateryear) %>% 
  summarise(wateryear_ppt=sum(prec), wateryear_tmean = mean(temp), wateryear_vpd = mean(vpdmax)) 


#seasons
#LAG weather
prev_mons <- prism_met %>% filter(season == 'prev_mons') %>% group_by(year) %>% 
  summarise(prev_mons_ppt=sum(prec), prev_mons_tmean = mean(temp), prev_mons_vpd = mean(vpdmax)) %>%
  mutate(prev_mons_ppt = lag(prev_mons_ppt, default = NA),
         prev_mons_tmean = lag(prev_mons_tmean, default = NA),
         prev_mons_vpd = lag(prev_mons_vpd, default = NA))%>%
  mutate(year = factor(year))
prev_fall <- prism_met %>% filter(season == 'prev_fall')  %>% group_by(year)%>% 
  summarise(prev_fall_ppt=sum(prec), prev_fall_tmean = mean(temp), prev_fall_vpd = mean(vpdmax)) %>%
  mutate(prev_fall_ppt = lag(prev_fall_ppt, default = NA),
         prev_fall_tmean = lag(prev_fall_tmean, default = NA),
         prev_fall_vpd = lag(prev_fall_vpd, default = NA))%>%
  mutate(year = factor(year))

#WINTER
winter <- prism_met %>% filter(season == 'winter')
winter_breaks <- seq(as.Date("1900-11-01"), length=122, by="year")  
years.winter.breaks = as.numeric(format(winter_breaks,"%Y"))
labels.winter = years.winter.breaks[2:length(winter_breaks)]  # Why are we leaving off the first year in the water years label?
winter$winter <- cut(winter$date, winter_breaks,labels=labels.winter)

winter <- winter %>% group_by(winter) %>% 
  summarise(winter_ppt=sum(prec), winter_tmean = mean(temp), winter_vpd = mean(vpdmax)) 


#CURRENT weather
pre_mons <- prism_met %>% filter(season == 'pre_mons') %>% group_by(year) %>% 
  summarise(pre_mons_ppt=sum(prec), pre_mons_tmean = mean(temp), pre_mons_vpd = mean(vpdmax)) %>%
  mutate(year = factor(year))

monsoon <- prism_met %>% filter(month >= 7 & month <= 9)  %>% group_by(year)%>% 
  summarise(monsoon_ppt=sum(prec), monsoon_tmean = mean(temp), monsoon_vpd = mean(vpdmax)) %>%
  mutate(year = factor(year))


#SUMMARISE PRISM
#coerce all as double
summarise_prism_met <- prev_mons %>%
  left_join(prev_fall)%>%
  left_join(winter,by = c("year" = "winter")) %>% 
  left_join(pre_mons) %>%
  left_join(monsoon) %>%
  left_join(prev_fall)%>%
  left_join(wateryear,by = c("year" = "wateryear"))%>%
  left_join(MB_all_rbind)%>%
  mutate(year = as.numeric(year)) %>%  # Ensure the year column is numeric
  filter(year > 1930)







top_10_percent_lowest_PPT_wateryear <- summarise_prism_met %>% 
  select(year,wateryear_ppt)%>%
  arrange(wateryear_ppt) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))


top_10_percent_lowest_PPT_winter <- summarise_prism_met %>%
  select(year,winter_ppt)%>%#select columns
  arrange(winter_ppt) %>%#arrange by smallest to largest
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_PPT_monsoon <- summarise_prism_met %>%
  select(year,monsoon_ppt)%>%
  arrange(monsoon_ppt) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_PPT_pre_monsoon <- summarise_prism_met %>% 
  select(year,pre_mons_ppt)%>%
  arrange(pre_mons_ppt) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))


top_10_percent_lowest_PPT_prev_fall <- summarise_prism_met %>% 
  select(year,prev_fall_ppt)%>%
  arrange(prev_fall_ppt) %>%
  slice(1:(n() * 0.10)) %>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))










top_10_percent_lowest_vpd_wateryear <- summarise_prism_met %>% 
  select(year,wateryear_vpd)%>%
  arrange(wateryear_vpd) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_vpd_winter <- summarise_prism_met %>%
  select(year,winter_vpd)%>%
  arrange(winter_vpd) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_vpd_monsoon <- summarise_prism_met %>%
  select(year,monsoon_vpd)%>%
  arrange(monsoon_vpd) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_vpd_pre_monsoon <- summarise_prism_met %>% 
  select(year,pre_mons_vpd)%>%
  arrange(pre_mons_vpd) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_vpd_prev_fall <- summarise_prism_met %>% 
  select(year,prev_fall_vpd)%>%
  arrange(prev_fall_vpd) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))
















top_10_percent_lowest_tmean_wateryear <- summarise_prism_met %>% 
  select(year,wateryear_tmean)%>%
  arrange(wateryear_tmean) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_tmean_winter <- summarise_prism_met %>%
  select(year,winter_tmean)%>%
  arrange(winter_tmean) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_tmean_monsoon <- summarise_prism_met %>%
  select(year,monsoon_tmean)%>%
  arrange(monsoon_tmean) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_tmean_pre_monsoon <- summarise_prism_met %>% 
  select(year,pre_mons_tmean)%>%
  arrange(pre_mons_tmean) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))

top_10_percent_lowest_tmean_prev_fall <- summarise_prism_met %>% 
  select(year,prev_fall_tmean)%>%
  arrange(prev_fall_tmean) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))







#SUMMARISE PRISM
#coerce all as double
#summarise_prism_met_15 <- top_10_percent_lowest_PPT_wateryear %>%
#  cbind(top_10_percent_lowest_PPT_winter)%>%
#  cbind(top_10_percent_lowest_PPT_monsoon)%>%
#  cbind(top_10_percent_lowest_PPT_pre_monsoon)%>%
#  cbind(top_10_percent_lowest_PPT_prev_fall)%>%
#  cbind(top_10_percent_lowest_vpd_wateryear) %>%
#  cbind(top_10_percent_lowest_vpd_winter)%>%
#  cbind(top_10_percent_lowest_vpd_monsoon)%>%
#  cbind(top_10_percent_lowest_vpd_pre_monsoon)%>%
#  cbind(top_10_percent_lowest_vpd_prev_fall)%>%
# cbind(top_10_percent_lowest_tmean_wateryear) %>%
#  cbind(top_10_percent_lowest_tmean_winter)%>%
#  cbind(top_10_percent_lowest_tmean_monsoon)%>%
#  cbind(top_10_percent_lowest_tmean_pre_monsoon)%>%
#  cbind(top_10_percent_lowest_tmean_prev_fall) 




prism_year <- prism_met %>% group_by(year) %>% 
  summarise(year_ppt=sum(prec), year_tmean = mean(temp), year_vpd = mean(vpdmax)) %>% 
  mutate(year = factor(year))

prism_month <- prism_met %>% group_by(month) %>% 
  summarise(month_ppt=mean(prec), month_tmean = mean(temp), month_vpd = mean(vpdmax))  

#seasons

menziesii = MB_all_rbind %>% dplyr::select(year,menziesii) %>% na.omit()
ponderosa = MB_all_rbind %>% dplyr::select(year,ponderosa) %>% na.omit()
strobiformis = MB_all_rbind %>% dplyr::select(year,strobiformis) %>% na.omit()
  
#mean stations
#noaa_summary = noaa_climate_csv %>%
#  group_by(DATE) %>%
#  mutate_at(c(4:6), as.numeric)%>%
#  summarize(avg_prcp = mean(PRCP, na.rm = TRUE),avg_tmax = mean(TMAX, na.rm = TRUE))

#noaa_summary_month = noaa_summary %>% na.omit() %>% mutate(DATE = as.Date(DATE, format="%Y-%m-%d"))%>%
#  mutate(month = month(DATE), year = year(DATE))%>%
#  group_by(month,year) %>% # group by the day column
#  summarise(total_precip=sum(avg_prcp), mean_tmax =mean(avg_tmax)) %>% 
#  mutate(year = factor(year))



MB_all_rbind_final = MB_all_rbind %>% mutate_at(c(1), as.double())
yearly_prism_treering <- MB_all_rbind_final%>% left_join(summarise_prism_met)

#monthly climate yearly ring width
#prism_noaa_treering = prism_year %>% left_join(noaa_summary_month)%>%
#  left_join(MB_all_rbind_final) 




#for (i in 1:nrow(prism_noaa_treering)){prism_noaa_treering$ppt_corr[i] = cor(prism_noaa_treering$ppt[1:i], prism_noaa_treering$ponderosa[1:i])
#}
#####CLIMATE#####

#write.csv(MB_all_ponderosa, file="MB_all_ponderosa.csv")
#write.csv(prism_met, file="prism_met.csv")
#write.csv(MB_all_ponderosa_detrend.rwi.crn, file="MB_all_ponderosa_detrend.rwi.crn.csv")








#######EDITING FOR PAPER ANLAYSIS and FIGURES
#####years
# Load necessary packages
#install.packages("SPEI")
library(SPEI)


# Latitude of your location
latitude <- 32.415476 # Replace with your actual latitude

# Adjust prism_met to water year (October to September)
prism_met$water_year <- ifelse(month(prism_met$date) >= 10, year(prism_met$date) + 1, year(prism_met$date))
prism_met$water_month <- ifelse(month(prism_met$date) >= 10, month(prism_met$date) - 9, month(prism_met$date) + 3)

# Calculate monthly average temperature and total precipitation for the water year
prism_met_spei <- prism_met %>%
  mutate(pet =thornthwaite(temp,latitude))%>%
  mutate(wb = prec - pet)%>%
  mutate(year = as.numeric(year)) %>%  # Ensure the year column is numeric
  filter(year > 1930)

# Calculate SPEI for a yearly timescale (12 months)
spei_wateryear_pet <- prism_met_spei %>% select(water_year,wb)%>%
  na.omit()

spei_wateryear <- spei(spei_wateryear_pet$wb, scale=12)

# Calculate SPEI for a yearly timescale (12 months)
spei_winter_pet <- prism_met_spei %>% select(year,season,wb)%>%
  filter(season =="winter")%>%
  na.omit()

spei_winter <- spei(spei_winter_pet$wb, scale=5)

# Calculate SPEI for a yearly timescale (12 months)
spei_monsoon_pet <- prism_met_spei %>% select(year,month,wb)%>%
  filter(month >= 7 & month <= 9)%>%
  na.omit()

spei_monsoon <- spei(spei_monsoon_pet$wb, scale=3)



# Calculate SPEI for a yearly timescale (12 months)
spei_pre_monsoon_pet <- prism_met_spei %>% select(year,month,wb)%>%
  filter(month >= 4 & month <= 6)%>%
  na.omit()

spei_pre_monsoon <- spei(spei_pre_monsoon_pet$wb, scale=3)

# Calculate SPEI for a yearly timescale (12 months)
spei_prev_fall_pet <- prism_met_spei %>% select(year,month,wb)%>%
  filter(month == 10)%>%
  na.omit()

spei_prev_fall <- spei(spei_prev_fall_pet$wb, scale=3)











# Extract yearly SPEI values
spei_yearly_mean <- as.data.frame(spei_wateryear$fitted) %>% cbind(spei_wateryear_pet) %>% rename(SPEI=1) %>%
  group_by(water_year) %>%
  summarise(SPEI = mean(SPEI, na.rm = TRUE))%>%na.omit()

  # Extract yearly SPEI values
spei_winter_mean <- as.data.frame(spei_winter$fitted) %>% cbind(spei_winter_pet) %>% rename(SPEI=1) %>%
  group_by(year)%>%
  summarise(SPEI = mean(SPEI, na.rm = TRUE))%>%na.omit()

  # Extract yearly SPEI values
spei_monsoon_mean <- as.data.frame(spei_monsoon$fitted) %>% cbind(spei_monsoon_pet) %>% rename(SPEI=1) %>%
  group_by(year) %>%
  summarise(SPEI = mean(SPEI, na.rm = TRUE))%>%na.omit()

# Extract yearly SPEI values
spei_prev_fall_mean <- as.data.frame(spei_prev_fall$fitted) %>% cbind(spei_prev_fall_pet) %>% rename(SPEI=1) %>%
  group_by(year)%>%
  summarise(SPEI = mean(SPEI, na.rm = TRUE))%>%na.omit()

# Extract yearly SPEI values
spei_pre_monsoon_mean <- as.data.frame(spei_pre_monsoon$fitted) %>% cbind(spei_pre_monsoon_pet) %>% rename(SPEI=1) %>%
  group_by(year) %>%
  summarise(SPEI = mean(SPEI, na.rm = TRUE))%>%na.omit()



# Plot the yearly SPEI values
plot(spei_yearly_mean$water_year, spei_yearly_mean$SPEI, type='l', main="Yearly SPEI (Water Year)", xlab="Water Year", ylab="SPEI")


##SEA iterations
top_10_percent_lowest_SPEI_wateryear <- spei_yearly_mean %>%
  arrange(SPEI) %>%
  slice(1:(n() * 0.10))%>%
  mutate(water_year = as.numeric(water_year)) %>%
  # Filter out water_years that are consecutive or within two water_years
  filter(!((water_year - lag(water_year) <= 2 & water_year - lag(water_year) > 0) |
             (lead(water_year) - water_year <= 2 & lead(water_year) - water_year > 0))) %>%
  # Revert water_year back to character if needed
  mutate(water_year = as.numeric(water_year))


top_10_percent_lowest_SPEI_winter <- spei_winter_mean %>%
  arrange(SPEI) %>%
  slice(1:(n() * 0.10))%>%
mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))


top_10_percent_lowest_SPEI_monsoon <- spei_monsoon_mean %>%
  arrange(SPEI) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))


top_10_percent_lowest_SPEI_prev_fall <- spei_prev_fall_mean %>%
  arrange(SPEI) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))


top_10_percent_lowest_SPEI_pre_monsoon <- spei_pre_monsoon_mean %>%
  arrange(SPEI) %>%
  slice(1:(n() * 0.10))%>%
  mutate(year = as.numeric(year)) %>%
  # Filter out years that are consecutive or within two years
  filter(!((year - lag(year) <= 2 & year - lag(year) > 0) |
             (lead(year) - year <= 2 & lead(year) - year > 0))) %>%
  # Revert year back to character if needed
  mutate(year = as.numeric(year))















##EPOCH
library(graphics)
library(utils)





windows()
par(mfrow = c(1, 5))
wateryear_PPT <- top_10_percent_lowest_PPT_wateryear$year

# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_PPT)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_PPT)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_menziesii$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_PPT)

strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_strobiformis$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, wateryear_PPT,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, wateryear_PPT,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, wateryear_PPT,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0, alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="wateryear_PPT")















winter_PPT <-top_10_percent_lowest_PPT_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_PPT)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_PPT)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_PPT)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_PPT,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_PPT,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_PPT,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_PPT")
















monsoon_PPT <-top_10_percent_lowest_PPT_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_PPT)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_PPT)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_PPT)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, monsoon_PPT,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, monsoon_PPT,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, monsoon_PPT,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="monsoon_PPT")





















pre_monsoon_PPT <-top_10_percent_lowest_PPT_pre_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_PPT)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_PPT)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_PPT)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, pre_monsoon_PPT,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, pre_monsoon_PPT,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, pre_monsoon_PPT,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="pre_monsoon_PPT")














prev_fall_PPT <-top_10_percent_lowest_PPT_prev_fall$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_PPT)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_PPT)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_PPT)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, prev_fall_PPT,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, prev_fall_PPT,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, prev_fall_PPT,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="prev_fall_PPT")













































windows()
par(mfrow = c(1, 5))

wateryear_tmean <- top_10_percent_lowest_tmean_wateryear$year

# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_tmean)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_tmean)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_menziesii$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_tmean)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_strobiformis$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, wateryear_tmean,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, wateryear_tmean,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, wateryear_tmean,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0, alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="wateryear_tmean")















winter_tmean <-top_10_percent_lowest_tmean_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_tmean)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_tmean)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_tmean)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_tmean,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_tmean,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_tmean,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_tmean")
















monsoon_tmean <-top_10_percent_lowest_tmean_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_tmean)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_tmean)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_tmean)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, monsoon_tmean,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, monsoon_tmean,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, monsoon_tmean,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="monsoon_tmean")





















pre_monsoon_tmean <-top_10_percent_lowest_tmean_pre_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_tmean)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_tmean)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_tmean)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, pre_monsoon_tmean,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, pre_monsoon_tmean,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, pre_monsoon_tmean,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="pre_monsoon_tmean")














prev_fall_tmean <-top_10_percent_lowest_tmean_prev_fall$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_tmean)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_tmean)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_tmean)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, prev_fall_tmean,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, prev_fall_tmean,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, prev_fall_tmean,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="prev_fall_tmean")




























































windows()
par(mfrow = c(1, 5))
wateryear_vpd <- top_10_percent_lowest_vpd_wateryear$year

# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_vpd)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_vpd)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_menziesii$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_vpd)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_strobiformis$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, wateryear_vpd,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, wateryear_vpd,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, wateryear_vpd,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0, alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="wateryear_vpd")















winter_vpd <-top_10_percent_lowest_vpd_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_vpd)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_vpd)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_vpd)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_vpd,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_vpd,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_vpd,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_vpd")
















monsoon_vpd <-top_10_percent_lowest_vpd_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_vpd)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_vpd)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_vpd)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, monsoon_vpd,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, monsoon_vpd,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, monsoon_vpd,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="monsoon_vpd")





















pre_monsoon_vpd <-top_10_percent_lowest_vpd_pre_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_vpd)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_vpd)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_vpd)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, pre_monsoon_vpd,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, pre_monsoon_vpd,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, pre_monsoon_vpd,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="pre_monsoon_vpd")














prev_fall_vpd <-top_10_percent_lowest_vpd_prev_fall$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_vpd)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_vpd)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_vpd)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, prev_fall_vpd,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, prev_fall_vpd,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, prev_fall_vpd,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="prev_fall_vpd")










































windows()
par(mfrow = c(1, 5))



wateryear_SPEI <- top_10_percent_lowest_SPEI_wateryear$water_year

# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_SPEI)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_SPEI)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_menziesii$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% wateryear_SPEI)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_strobiformis$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, wateryear_SPEI,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, wateryear_SPEI,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, wateryear_SPEI,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0, alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="wateryear_SPEI")















winter_SPEI <-top_10_percent_lowest_SPEI_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_SPEI)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_SPEI)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_SPEI)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_SPEI,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_SPEI,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_SPEI,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_SPEI")
















monsoon_SPEI <-top_10_percent_lowest_SPEI_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_SPEI)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_SPEI)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% monsoon_SPEI)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, monsoon_SPEI,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, monsoon_SPEI,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, monsoon_SPEI,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="monsoon_SPEI")





















pre_monsoon_SPEI <-top_10_percent_lowest_SPEI_pre_monsoon$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_SPEI)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>% rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_SPEI)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% pre_monsoon_SPEI)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, pre_monsoon_SPEI,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, pre_monsoon_SPEI,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, pre_monsoon_SPEI,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="pre_monsoon_SPEI")














prev_fall_SPEI <-top_10_percent_lowest_SPEI_prev_fall$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_SPEI)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_SPEI)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% prev_fall_SPEI)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, prev_fall_SPEI,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, prev_fall_SPEI,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, prev_fall_SPEI,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n", 
     cex.lab = 2, font.lab = 2,  # Increase and bold axis labels
     cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

legend("bottomleft",legend=c("PIPO","PSME","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)
title(main="prev_fall_SPEI")
































































windows()
par(mfrow = c(1, 4))







winter_PPT <-top_10_percent_lowest_PPT_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_PPT)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_PPT)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_PPT)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_PPT,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_PPT,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_PPT,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_PPT")











winter_tmean <-top_10_percent_lowest_tmean_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_tmean)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_tmean)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_tmean)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_tmean,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_tmean,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_tmean,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_tmean")




















winter_vpd <-top_10_percent_lowest_vpd_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_vpd)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_vpd)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_vpd)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_vpd,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_vpd,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_vpd,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_vpd")















winter_SPEI <-top_10_percent_lowest_SPEI_winter$year

# Calculate normal mean and stdev growth rates for each species
# Ponderosa Pine
normal_ponderosa <- MB_all_ponderosa_detrend.rwi.crn %>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_SPEI)
ponderosa_mean <- mean(normal_ponderosa$std, na.rm = TRUE)
ponderosa_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Douglas-fir
normal_menziesii <- MB_all_menziesii_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_SPEI)
menziesii_mean <- mean(normal_menziesii$std, na.rm = TRUE)
menziesii_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Southwestern White Pine
normal_strobiformis <- MB_all_strobiformis_detrend.rwi.crn%>% 
  rownames_to_column()%>%   rename(year = rowname)%>%
  filter(!year %in% winter_SPEI)
strobiformis_mean <- mean(normal_strobiformis$std, na.rm = TRUE)
strobiformis_sd <- sd(normal_ponderosa$std,na.rm=TRUE)

# Calculate upper and lower bounds for the standard error shading
# Calculate upper and lower bounds for each series
upper_bound_ponderosa <- rep(ponderosa_mean + ponderosa_sd, 5)
lower_bound_ponderosa <- rep(ponderosa_mean - ponderosa_sd, 5)

upper_bound_menziesii <- rep(menziesii_mean + menziesii_sd, 5)
lower_bound_menziesii <- rep(menziesii_mean - menziesii_sd, 5)

upper_bound_strobiformis <- rep(strobiformis_mean + strobiformis_sd, 5)
lower_bound_strobiformis <- rep(strobiformis_mean - strobiformis_sd, 5)

ponderosa.sea <- sea(MB_all_ponderosa_detrend.rwi.crn, winter_SPEI,lag = 2, resample = 1000)
ponderosa.foo <- ponderosa.sea$se.unscaled
names(ponderosa.foo) <- ponderosa.sea$lag

menziesii.sea <- sea(MB_all_menziesii_detrend.rwi.crn, winter_SPEI,lag = 2, resample = 1000)
menziesii.foo <- menziesii.sea$se.unscaled
names(menziesii.foo) <- menziesii.sea$lag

strobiformis.sea <- sea(MB_all_strobiformis_detrend.rwi.crn, winter_SPEI,lag = 2, resample = 1000)
strobiformis.foo <- strobiformis.sea$se.unscaled
names(strobiformis.foo) <- strobiformis.sea$lag

#windows()
plot(ponderosa.foo, col = 'red',  type = "l", ylim = c(0.7,1.3),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.foo, col = ifelse(ponderosa.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(menziesii.foo, col = ifelse(menziesii.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))
points(strobiformis.foo, col = ifelse(strobiformis.sea$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.sea$p < 00.10,8, 16),cex=ifelse(ponderosa.sea$p < 00.10,1.5, 1))

# Add mean lines
abline(h = ponderosa_mean, col = 'red', lty = 2)
abline(h = menziesii_mean, col = 'green', lty = 2)
abline(h = strobiformis_mean, col = 'blue', lty = 2)

polygon(x = c(1:length(ponderosa.foo), rev(1:length(ponderosa.foo))),
        y = c(lower_bound_ponderosa, rev(upper_bound_ponderosa)),
        col = rgb(0, 0, 0,alpha = 00.10), border = 'red')

# Shading for Douglas-fir with transparency
polygon(x = c(1:length(menziesii.foo), rev(1:length(menziesii.foo))),
        y = c(lower_bound_menziesii, rev(upper_bound_menziesii)),
        col = rgb(0, 0, 1,alpha = 00.10), border = 'green')

# Shading for Southwestern White Pine with transparency
polygon(x = c(1:length(strobiformis.foo), rev(1:length(strobiformis.foo))),
        y = c(lower_bound_strobiformis, rev(upper_bound_strobiformis)),
        col = rgb(1, 0, 0,alpha = 00.10), border = 'blue')

# Add legend and axis
legend("bottomleft", legend = c("PIPO", "PSME", "PISF"), lty = 1, col = c('red', 'green', 'blue'), cex = 1, bty = "n")
axis(1, at = 1:5, labels = c("-2", "-1", "0", "1", "2"),cex.main = 2)
title(main = "winter_SPEI")













































































































































































sea_mod <- function (x, key, lag = 2, symmetric = 1, remove_evnts = 0, resample = 1000)
  
{
  
  if (!is.data.frame(x)) {
    
    stop("'x' must be a data.frame")
    
  }
  
  stopifnot(is.numeric(lag), length(lag) == 1, is.finite(lag),
            
            lag >= 0, round(lag) == lag)
  
  stopifnot(is.numeric(resample), length(resample) == 1, is.finite(resample),
            
            resample >= 1, round(resample) == resample)
  
  if (ncol(x) >= 1) {
    
    x.unscaled <- x[1]
    
  } else {
    
    stop("'x' must have at least one column")
    
  }
  
  rnames <- row.names(as.matrix(x.unscaled))
  
  if (is.null(rnames)) {
    
    stop("'x' must have non-automatic row.names")
    
  }
  
  rnames <- as.numeric(rnames)
  
  x.scaled <- data.frame(scale(x.unscaled))
  
  n <- length(key)
  
  if (symmetric==1) {
    
    m <- 2 * lag + 1
    
    se.table <- matrix(NA_real_, ncol = m, nrow = n)
    
    se.unscaled.table <- se.table
    
    yrs.base <- (-lag):(m - lag - 1)
    
  } else {
    
    m <- lag + 1
    
    se.table <- matrix(NA_real_, ncol = m, nrow = n)
    
    se.unscaled.table <- se.table
    
    yrs.base <- 0:lag
    
  }
  
  seq.n <- seq_len(n)
  
  for (i in seq.n) {
    
    yrs <- as.character(key[i] + yrs.base)
    
    se.table[i, ] <- x.scaled[yrs, ]
    
    se.unscaled.table[i, ] <- x.unscaled[yrs, ]
    
  }
  
  se <- colMeans(se.table, na.rm = TRUE)
  
  se.unscaled <- colMeans(se.unscaled.table, na.rm = TRUE)
  
  re.table <- matrix(NA_real_, ncol = m, nrow = resample)
  
  re.subtable <- matrix(NA_real_, ncol = m, nrow = n)
  
  l <- length(rnames)
  
  ifelse(symmetric==1, trim <- c(1:lag, (l - lag + 1):l),trim <- c((l - lag + 1):l))
  
  rnames.blue <- rnames[-trim]
  
  if (remove_evnts==1) {
    
    # F&P
    
    cnt <- 1
    
    n.blue <- length(blueuce(intersect, c(list(rnames.blue),list(key)))) # to account for years removed by trim
    
    yrs.base.pos <-yrs.base[yrs.base>(-1)]
    
    ev_and_lag_yrs <- rep(NA,length(yrs.base.pos)*n.blue)
    
    for (y in 1:length(yrs.base.pos)) {
      
      ev_and_lag_yrs[cnt:(cnt+n.blue-1)] <- which(rnames.blue %in% key)+yrs.base.pos[y]
      
      cnt <- cnt + n.blue
      
    }
    
    rnames.blue2 <- rnames.blue[-ev_and_lag_yrs]
    
  } else {
    
    rnames.blue2 <- rnames.blue
    
  }     
  
  for (k in seq_len(resample)) {
    
    rand.key <- sample(rnames.blue2, n, replace = TRUE)
    
    for (i in seq.n) {
      
      re.subtable[i, ] <- x.scaled[as.character(rand.key[i] + yrs.base), ]
      
    }
    
    re.table[k, ] <- colMeans(re.subtable, na.rm = TRUE)
    
  }
  
  p <- rep(NA_real_, m)
  
  re_ecdfs <- apply(re.table, 2, ecdf)
  
  ci <- apply(re.table, 2, quantile, probs = c(0.025, 0.975,
                                               
                                               0.005, 0.995))
  
  if (resample < 1000) {
    
    warning("'resample' is lower than 1000, potentially leading to less accuracy in estimating p-values and confidence bands.")
    
  }
  
  for (i in seq_len(m)) {
    
    if (is.na(se[i])) {
      
      warning(gettextf("NA result at position %d. ", i),
              
              "You could check whether 'key' years are in range.")
      
    }
    
    else {
      
      if (se[i] < 0) {
        
        p[i] <- re_ecdfs[[i]](se[i])
        
      }
      
      else {
        
        p[i] <- 1 - re_ecdfs[[i]](se[i])
        
      }
      
    }
    
  }
  
  
  
  if (symmetric==1) {
    
    data.frame(symmetric,remove_evnts,lag = c(-lag:lag), se, se.unscaled, p, ci.95.lower = ci[1,],
               
               ci.95.upper = ci[2, ], ci.99.lower = ci[3, ], ci.99.upper = ci[4,])
    
  } else {
    
    data.frame(symmetric,remove_evnts, lag = c(0:lag), se, se.unscaled, p, ci.95.lower = ci[1,],
               
               ci.95.upper = ci[2, ], ci.99.lower = ci[3, ], ci.99.upper = ci[4,])
    
  }
  
}
####EPOCH########
ponderosa.wateryear_seamod <- sea_mod(MB_all_ponderosa_detrend.rwi.crn, top_10_percent_lowest_PPT_wateryear,lag = 2,symmetric = 1, remove_evnts = 1, resample = 1000)
ponderosa.wateryear_foo <- ponderosa.wateryear_seamod$se.unscaled
names(ponderosa.wateryear_foo) <- ponderosa.wateryear_seamod$lag

menziesii.wateryear_seamod <- sea_mod(MB_all_menziesii_detrend.rwi.crn, wateryear_drought,lag = 2,symmetric = 1, remove_evnts = 1, resample = 1000)
menziesii.wateryear_foo <- menziesii.wateryear_seamod$se.unscaled
names(menziesii.wateryear_foo) <- menziesii.wateryear_seamod$lag

strobiformis.wateryear_seamod <- sea_mod(MB_all_strobiformis_detrend.rwi.crn, wateryear_drought,lag = 2,symmetric = 1, remove_evnts = 1, resample = 1000)
strobiformis.wateryear_foo <- strobiformis.wateryear_seamod$se.unscaled
names(strobiformis.wateryear_foo) <- strobiformis.wateryear_seamod$lag

windows()
plot(ponderosa.wateryear_foo, col = 'red',  type = "l", ylim = c(0.73,10.10),lwd=2,
     ylab = "RWI", xlab = "Superposed Epoch", xaxt ="n",  cex.lab = 2, font.lab = 2,  cex.axis = 2, font.axis = 2)
lines(menziesii.wateryear_foo, col = 'green',  type = "l",lwd =2)
lines(strobiformis.wateryear_foo, col = 'blue',  type = "l", lwd=2)
points(ponderosa.wateryear_foo, col = ifelse(ponderosa.wateryear_seamod$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.wateryear_seamod$p < 00.10,8, 16),cex=ifelse(ponderosa.wateryear_seamod$p < 00.10,1.5, 1))
points(menziesii.wateryear_foo, col = ifelse(menziesii.wateryear_seamod$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.wateryear_seamod$p < 00.10,8, 16),cex=ifelse(ponderosa.wateryear_seamod$p < 00.10,1.5, 1))
points(strobiformis.wateryear_foo, col = ifelse(strobiformis.wateryear_seamod$p < 00.10, "black", 'black'), 
       pch =ifelse(ponderosa.wateryear_seamod$p < 00.10,8, 16),cex=ifelse(ponderosa.wateryear_seamod$p < 00.10,1.5, 1))
legend("bottomleft",legend=c("PIPO","PM","PISF"),lty=1, col = c('red','green','blue'),cex = 1,bty="n")
axis(1, at=1:5, labels=c("-2","-1","0","1","2"),cex.main = 2)



########TREE CLIMATE ASSESSMENT##########
length(unique(rownames(MB_all_ponderosa_detrend.rwi.crn))) >= (ncol(prism_met)-1) 




tmp_pivot <- prism_met %>% dplyr::select(2:4) %>%  mutate(across(1:3, as.character)) %>% 
  pivot_wider(names_from = month, values_from = temp) %>% lapply(as.numeric)%>% as.data.frame() #makes numeric

ppt_pivot <- prism_met %>% dplyr::select(c(2:3,5)) %>%  mutate(across(1:3, as.character)) %>% 
  pivot_wider(names_from = month, values_from = prec) %>% lapply(as.numeric)%>% as.data.frame() #makes numeric

vpd_pivot <- prism_met %>% dplyr::select(c(2:3,6)) %>%  mutate(across(1:3, as.character)) %>% 
  pivot_wider(names_from = month, values_from = vpdmax) %>% lapply(as.numeric)%>% as.data.frame() #makes numeric




#ponderosa
ponderosa_response = dcc(MB_all_ponderosa_detrend.rwi.crn,
                         climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                         var_names = c("tmp", "ppt", "vpd"),
                         selection = 6:10, dynamic = "moving", sb = FALSE,
                         method = "response") 

ponderosa_response_species = ponderosa_response$coef %>% mutate(species = "ponderosa")

ponderosa_corr = dcc(MB_all_ponderosa_detrend.rwi.crn,
                            climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                            var_names = c("tmp", "ppt", "vpd"),
                            method = "correlation") 

ponderosa_corr_species = ponderosa_corr$coef %>% mutate(species = "ponderosa")

ponderosa_corr_window <- dcc(MB_all_ponderosa_detrend.rwi.crn,
                      climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                      selection = 6:10, dynamic = "moving", sb = FALSE,
                      var_names = c("tmp", "ppt", "vpd"))

ponderosa_corr_evolving <- dcc(MB_all_ponderosa_detrend.rwi.crn,
                             climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                             selection = 6:10, dynamic = "evolving", sb = FALSE,
                             var_names = c("tmp", "ppt", "vpd"))
ponderosa_corr_temporal<- dcc(MB_all_ponderosa_detrend.rwi.crn,
           climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
           var_names = c("tmp", "ppt", "vpd"),
           selection = .mean(c(1:3,1:12)) + .mean(6:9), method = "cor",
           dynamic = "moving", win_size = 30, sb = FALSE)


ponderosa_sc <- seascorr(MB_all_ponderosa_detrend.rwi.crn,
                climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                var_names = c("tmp", "ppt", "vpd"))

ponderosa_sc_V2 <- seascorr(MB_all_ponderosa_detrend.rwi.crn,
                climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                var_names = c("tmp", "ppt", "vpd"), primary = "ppt", secondary = "tmp")





#menziesii
menziesii_response = dcc(MB_all_menziesii_detrend.rwi.crn,
                         climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                         var_names = c("tmp", "ppt", "vpd"),
                         method = "response") 
menziesii_response_species = menziesii_response$coef %>% mutate(species = "menziesii")

menziesii_corr = dcc(MB_all_menziesii_detrend.rwi.crn,
                     climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                     var_names = c("tmp", "ppt", "vpd"),
                     method = "correlation")
menziesii_corr_species = menziesii_corr$coef %>% mutate(species = "menziesii")

menziesii_corr_window <- dcc(MB_all_menziesii_detrend.rwi.crn,
                             climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                             selection = 6:10, dynamic = "moving", sb = FALSE,
                             var_names = c("tmp", "ppt", "vpd"))
menziesii_corr_evolving <- dcc(MB_all_menziesii_detrend.rwi.crn,
                               climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                               selection = 6:10, dynamic = "evolving", sb = FALSE,
                               var_names = c("tmp", "ppt", "vpd"))
menziesii_corr_temporal<- dcc(MB_all_menziesii_detrend.rwi.crn,
                              climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                              var_names = c("tmp", "ppt", "vpd"),
                              selection = .mean(c(1:3,1:12)) + .mean(6:9), method = "cor",
                              dynamic = "moving", win_size = 30, sb = FALSE)

menziesii_sc <- seascorr(MB_all_menziesii_detrend.rwi.crn,
                         climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                         var_names = c("tmp", "ppt", "vpd"))

menziesii_sc_V2 <- seascorr(MB_all_menziesii_detrend.rwi.crn,
                            climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                            var_names = c("tmp", "ppt", "vpd"), primary = "ppt", secondary = "tmp")





tmp_pivot <- prism_met %>% dplyr::select(2:4) %>%  mutate(across(1:3, as.character)) %>% 
  pivot_wider(names_from = month, values_from = temp) %>% lapply(as.numeric)%>% as.data.frame() %>% slice(-c(1:7)) #makes numeric 

ppt_pivot <- prism_met %>% dplyr::select(c(2:3,5)) %>%  mutate(across(1:3, as.character)) %>% 
  pivot_wider(names_from = month, values_from = prec) %>% lapply(as.numeric)%>% as.data.frame() %>% slice(-c(1:7)) #makes numeric

vpd_pivot <- prism_met %>% dplyr::select(c(2:3,6)) %>%  mutate(across(1:3, as.character)) %>% 
  pivot_wider(names_from = month, values_from = vpdmax) %>% lapply(as.numeric)%>% as.data.frame() %>% slice(-c(1:7)) #makes numeric


#strobiformis
strobiformis_response = dcc(MB_all_strobiformis_detrend.rwi.crn,
                         climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                         var_names = c("tmp", "ppt", "vpd"),
                         method = "response") 
strobiformis_response_species = strobiformis_response$coef %>% mutate(species = "strobiformis")

strobiformis_corr = dcc(MB_all_strobiformis_detrend.rwi.crn,
                     climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                     var_names = c("tmp", "ppt", "vpd"),
                     method = "correlation")
strobiformis_corr_species = strobiformis_corr$coef %>% mutate(species = "strobiformis")

strobiformis_corr_window <- dcc(MB_all_strobiformis_detrend.rwi.crn,
                             climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                             selection = 6:10, dynamic = "moving", sb = FALSE,
                             var_names = c("tmp", "ppt", "vpd"))
strobiformis_corr_evolving <- dcc(MB_all_strobiformis_detrend.rwi.crn,
                               climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                               method = "correlation",
                               selection = 6:10, dynamic = "evolving", sb = FALSE,
                               var_names = c("tmp", "ppt", "vpd"))

strobiformis_corr_temporal<- dcc(MB_all_strobiformis_detrend.rwi.crn,
                              climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                              var_names = c("tmp", "ppt", "vpd"),
                              selection = .mean(c(1:3,1:12)) + .mean(6:9), method = "cor",
                              dynamic = "moving", win_size = 30, sb = FALSE)

strobiformis_sc <- seascorr(MB_all_strobiformis_detrend.rwi.crn,
                         climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                         var_names = c("tmp", "ppt", "vpd"))

strobiformis_sc_V2 <- seascorr(MB_all_strobiformis_detrend.rwi.crn,
                            climate = list(tmp_pivot, ppt_pivot, vpd_pivot),
                            var_names = c("tmp", "ppt", "vpd"), primary = "ppt", secondary = "tmp")



species_bind_corr = ponderosa_corr_species %>% rbind(menziesii_corr_species) %>% rbind(strobiformis_corr_species)
species_bind_resp = ponderosa_response_species %>% rbind(menziesii_response_species) %>% rbind(strobiformis_response_species)

tmp_cor = species_bind_corr %>% filter(varname == "tmp")
ppt_cor = species_bind_corr %>% filter(varname == "ppt")
vpd_cor = species_bind_corr %>% filter(varname == "vpd")

tmp_resp = species_bind_resp %>% filter(varname == "tmp")
ppt_resp = species_bind_resp %>% filter(varname == "ppt")
vpd_resp = species_bind_resp %>% filter(varname == "vpd")













































##################PLOTS##############
####TIMESERIES PLOTS####
###treering - species
windows()
plot(MB_all_rbind$menziesii~MB_all_rbind$year, type="l",xlab="Year",ylab="RWI", col = 'green')
lines(MB_all_rbind$ponderosa~MB_all_rbind$year, col = 'red')
lines(MB_all_rbind$strobiformis~MB_all_rbind$year, col = 'blue')

par(new=T)
plot(summarise_prism_met$wateryear_ppt ~summarise_prism_met$year, col = "orange",type="l", 
     axes=F, ylim=c(min(summarise_prism_met$wateryear_ppt, na.rm = TRUE),max(summarise_prism_met$wateryear_ppt, na.rm = TRUE)), xlab="", ylab="", main="")


par(new=T)
plot(summarise_prism_met$winter_ppt~summarise_prism_met$year, col = "black",type="l", 
     axes=F, ylim=c(min(summarise_prism_met$winter_ppt, na.rm = TRUE),max(summarise_prism_met$winter_ppt, na.rm = TRUE)), xlab="", ylab="", main="")

par(new=T)
plot(summarise_prism_met$pre_mons_vpd ~summarise_prism_met$year, col = "orange",type="l", 
     axes=F, ylim=c(min(summarise_prism_met$pre_mons_vpd, na.rm = TRUE),max(summarise_prism_met$pre_mons_vpd, na.rm = TRUE)), xlab="", ylab="", main="")




###mean longterm monthly timeseries cliamte
windows()
par(mar=c(5, 12, 4, 4) + 0.1)

plot(prism_month$month, prism_month$month_ppt, axes=F, ylim=c(0,max(prism_month$month_ppt)), xlab="", ylab="",type="l",col='black', main="")
points(prism_month$month, prism_month$month_ppt,pch=20,col='black')
axis(2, ylim=c(0,max(prism_month$month)),col='black',lwd=2)
mtext(2,text="Mean Precipitation",line=2)

par(new=T)
plot(prism_month$month, prism_month$month_vpd, axes=F, ylim=c(0,max(prism_month$month_vpd)), xlab="", ylab="",type="l",lty=2,col='black', main="")
points(prism_month$month, prism_month$month_vpd,pch=20,col='black')
axis(2, ylim=c(0,max(prism_month$month)),lwd=2,line=3.5)
mtext(2,text="Mean Max Vapor Pressure Deficit",line=5.5)

par(new=T)
plot(prism_month$month, prism_month$month_tmean, axes=F, ylim=c(0,max(prism_month$month_tmean)), xlab="", ylab="",type="l",lty=3,col='black', main="")
points(prism_month$month, prism_month$month_tmean,pch=20,col='black')
axis(2, ylim=c(0,max(prism_month$month)),lwd=2,line=7)
mtext(2,text="Mean Temperature",line=9)
axis(1,pretty(range(prism_month$month),10))

mtext("Month",side=1,col='black',line=2)
legend(x=6.5,y=2,legend=c("Mean Precipitation","Mean Max Vapor Pressure Deficit","Mean Temperature"),lty=c(1,2,3), col = c('black','black','black'),cex = 0.6,bty="n")

##########timeseries plots#########################




































########CLIMATE and TREE assessment PLOTS###############

species_temp_corr <- tmp_cor %>%
  ggplot(aes(x=factor(id), y=coef, group=species, color=species)) +
  geom_point(aes(fill = significant, pch = significant, size = significant))  +
  geom_line() + scale_color_manual(values = c('red','green', 'blue')) + ylab("Correlation Coefficients") + xlab("month")+ ggtitle("tmp") +theme_light()+
  scale_x_discrete( labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP"))

species_ppt_corr <- ppt_cor %>%
  ggplot(aes(x=factor(id), y=coef, group=species, color=species)) +
  geom_point(aes(fill = significant, pch = significant, size = significant)) +
  geom_line() + scale_color_manual(values = c('red','green', 'blue')) + ylab("Correlation Coefficients") + xlab("month")+ ggtitle("ppt") +theme_light()+
  scale_x_discrete( labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP"))

species_vpd_corr <- vpd_cor %>%
  ggplot(aes(x=factor(id), y=coef, group=species, color=species)) +
  geom_point(aes(fill = significant, pch = significant, size = significant)) +
  geom_line() + scale_color_manual(values = c('red','green', 'blue')) + ylab("Correlation Coefficients") + xlab("month")+ ggtitle("vpd") +theme_light()+
  scale_x_discrete( labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP"))

windows()
species_temp_corr/species_ppt_corr/species_vpd_corr








species_temp_resp <- tmp_resp %>%
  ggplot(aes(x=factor(id), y=coef, group=species, color=species)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper,linetype = significant), width = 0.5, linewidth = 1) +
  geom_point(aes(fill = significant, pch = significant, size = significant))  +
scale_color_manual(values = c('green','red', 'blue')) + ylab("Response Coefficients") + xlab("month")+ ggtitle("tmp") +theme_light()+
  scale_linetype_manual(values = c("FALSE" = "dotted", "TRUE" = "solid")) +  # Specify linetypes here
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold")) +scale_x_discrete( labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP"))

species_ppt_resp <- ppt_resp %>%
  ggplot(aes(x=factor(id), y=coef, group=species, color=species)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper,linetype = factor(significant)), width = 0.5, linewidth = 1) +
  geom_point(aes(fill = significant, pch = significant, size = significant)) +
  scale_color_manual(values = c('green','red', 'blue')) + ylab("Response Coefficients") + xlab("month")+ ggtitle("ppt") +theme_light()+
  scale_linetype_manual(values = c("FALSE" = "dotted", "TRUE" = "solid")) +  # Specify linetypes here
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold")) +scale_x_discrete( labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP"))

species_vpd_resp <- vpd_resp %>%
  ggplot(aes(x=factor(id), y=coef, group=species, color=species)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper,linetype = factor(significant)), width = 0.5, linewidth = 1) +
  geom_point(aes(fill = significant, pch = significant, size = significant)) +
  scale_linetype_manual(values = c("FALSE" = "dotted", "TRUE" = "solid")) +  # Specify linetypes here
  scale_color_manual(values = c('green','red', 'blue')) + ylab("Response Coefficients") + xlab("month")+ ggtitle("vpd") +theme_light()+
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold")) +scale_x_discrete( labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP"))

windows()
species_temp_resp/species_ppt_resp/species_vpd_resp















p_temporal = plot(ponderosa_corr_temporal)+  scale_fill_gradient(low = "yellow", high = 'blue', na.value = NA) 
m_temporal = plot(menziesii_corr_temporal)+  scale_fill_gradient(low = "yellow", high = 'blue', na.value = NA)
s_temporal = plot(strobiformis_corr_temporal)+  scale_fill_gradient(low = "yellow", high = 'blue', na.value = NA)
windows()
p_temporal/m_temporal/s_temporal

p_sc = plot(ponderosa_sc) + theme_bw() + scale_fill_brewer(palette = "Dark2") 
m_sc = plot(menziesii_sc) + theme_bw() + scale_fill_brewer(palette = "Dark2") 
s_sc = plot(strobiformis_sc) + theme_bw() + scale_fill_brewer(palette = "Dark2") 
windows()
p_sc/m_sc/s_sc

p_sc_V2 =plot(ponderosa_sc_V2)+ theme_bw() + scale_fill_brewer(palette = "Dark2") 
m_sc_V2 = plot(menziesii_sc_V2)+ theme_bw() + scale_fill_brewer(palette = "Dark2") 
s_sc_V2 = plot(strobiformis_sc_V2)+ theme_bw() + scale_fill_brewer(palette = "Dark2") 
windows()
p_sc_V2/m_sc_V2/s_sc_V2



p_corr = plot(ponderosa_corr)+ scale_color_manual(values = c('green', 'red', 'blue')) + ylab("Correlation Coefficients") + ggtitle("ponderosa") 
m_corr = plot(menziesii_corr)+ scale_color_manual(values = c('green', 'red', 'blue')) + ylab("Correlation Coefficients") + ggtitle("menziesii")
s_corr = plot(strobiformis_corr)+ scale_color_manual(values = c('green', 'red', 'blue')) + ylab("Correlation Coefficients") + ggtitle("strobiformis")
windows()
p_corr/m_corr/s_corr


p_response = plot(ponderosa_response)+ scale_color_manual(values = c('green', 'red', 'blue')) + ylab("Response Coefficients") + ggtitle("ponderosa")
m_response = plot(menziesii_response)+ scale_color_manual(values = c('green', 'red', 'blue')) + ylab("Response Coefficients") + ggtitle("menziesii")
s_response = plot(strobiformis_response)+ scale_color_manual(values = c('green', 'red', 'blue')) + ylab("Response Coefficients") + ggtitle("strobiformis")
windows()
p_response/m_response/s_response


















##SIGNIFICANT CLIMATE SEASONs###
windows()
plot(summarise_prism_met$ponderosa~summarise_prism_met$winter_ppt, xlab = "winter precipitation", ylab = "RWI",cex.lab = 1.5,cex.axis=1.5)
points(summarise_prism_met$menziesii~summarise_prism_met$winter_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$winter_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$winter_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_ppt), lwd=4, lty=2,col = 'blue')


mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$winter_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')








windows()
plot(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_vpd, xlab = "pre-monsoon VPD", ylab = "RWI",cex.lab = 1.5,cex.axis=1.5)
points(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_vpd, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_vpd, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



windows()
plot(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_ppt, xlab = "wateryear_ppt", ylab = "RWI",cex.lab = 1.5,cex.axis=1.5)
points(summarise_prism_met$menziesii~summarise_prism_met$wateryear_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_ppt), lwd=4, lty=2,col = 'blue')


mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')











windows()

plot(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_vpd, xlab = "wateryear_vpd", ylab = "RWI")
points(summarise_prism_met$menziesii~summarise_prism_met$wateryear_vpd, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_vpd, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')

















































































































###SEASONAL



#ppt
plot(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_ppt, xlab = "ppt", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_ppt), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')




plot(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_ppt, xlab = "ppt", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_ppt), lwd=4, lty=2,col = 'blue')


mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')






windows()
plot(summarise_prism_met$ponderosa~summarise_prism_met$winter_ppt, xlab = "winter precipitation", ylab = 'RWI')
points(summarise_prism_met$menziesii~summarise_prism_met$winter_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$winter_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$winter_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_ppt), lwd=4, lty=2,col = 'blue')


mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$winter_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')







plot(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_ppt, xlab = "ppt", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_ppt), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')







plot(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_ppt, xlab = "ppt", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$monsoon_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$monsoon_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_ppt), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$monsoon_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



windows()
plot(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_ppt, xlab = "wateryear_ppt", ylab = "RWI")
points(summarise_prism_met$menziesii~summarise_prism_met$wateryear_ppt, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_ppt, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_ppt), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_ppt), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_ppt), lwd=4, lty=2,col = 'blue')


mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_ppt)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



















#tmean
plot(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_tmean, xlab = "tmean", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_tmean, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_tmean, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_tmean), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_tmean), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_tmean), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')




plot(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_tmean, xlab = "tmean", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_tmean, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_tmean, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_tmean), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_tmean), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_tmean), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')





plot(summarise_prism_met$ponderosa~summarise_prism_met$winter_tmean, xlab = "tmean", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$winter_tmean, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$winter_tmean, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_tmean), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$winter_tmean), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_tmean), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$winter_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')




windos()
plot(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_tmean, xlab = "tmean", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_tmean, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_tmean, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_tmean), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_tmean), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_tmean), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



plot(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_tmean, xlab = "tmean", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$monsoon_tmean, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_tmean, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_tmean), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$monsoon_tmean), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_tmean), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$monsoon_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



plot(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_tmean, xlab = "tmean", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$wateryear_tmean, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_tmean, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_tmean), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_tmean), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_tmean), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_tmean)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')


















#vpd
plot(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_vpd, xlab = "vpd", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_vpd, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_vpd, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$prev_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



plot(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_vpd, xlab = "vpd", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_vpd, col = 'black')
points(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_vpd, col = 'black')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$prev_fall_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$prev_fall_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$prev_fall_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



plot(summarise_prism_met$ponderosa~summarise_prism_met$winter_vpd, xlab = "vpd", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$winter_vpd, col = 'black')
points(summarise_prism_met$strobiformis~summarise_prism_met$winter_vpd, col = 'black')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$winter_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$winter_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$winter_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$winter_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')


windows()
plot(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_vpd, xlab = "vpd", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_vpd, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_vpd, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$pre_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$pre_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$pre_mons_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')



plot(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_vpd, xlab = "vpd", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$monsoon_vpd, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_vpd, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$monsoon_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$monsoon_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$monsoon_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$monsoon_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')

windows()

plot(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_vpd, xlab = "wateryear_vpd", ylab = "menziesii")
points(summarise_prism_met$menziesii~summarise_prism_met$wateryear_vpd, col = 'green')
points(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_vpd, col = 'blue')
abline(lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_vpd), lwd=4, lty=2,col = 'red')
abline(lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_vpd), lwd=4, lty=2,col = 'green')
abline(lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_vpd), lwd=4, lty=2,col = 'blue')

mod1 = lm(summarise_prism_met$menziesii~summarise_prism_met$wateryear_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n', cex=1.5, text.col = 'green')


mod1 = lm(summarise_prism_met$ponderosa~summarise_prism_met$wateryear_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n', cex=1.5, text.col = 'red')


mod1 = lm(summarise_prism_met$strobiformis~summarise_prism_met$wateryear_vpd)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
my.p = modsum$coefficients[2,4]
my.slope = modsum$coefficients[2]
rp = vector('expression',3)
rp[1] = substitute(expression("" ~ italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression("" ~ italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = substitute(expression("" ~ italic(slope) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.slope, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n', cex=1.5, text.col = 'blue')




















