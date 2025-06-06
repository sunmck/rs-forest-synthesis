#########################################################
##
## Aggregate Data to Hexagons (parallel implementation)
##
#########################################################


library(terra)
library(sf)
library(plyr, lib.loc = "/usr/local/lib/R/site-library")
library(dplyr)
library(parallel)


# TODO: number of cores available for parallelization
num_cores <- 14


#### load rasterdata ####

# forest canopy cover loss (epsg3035)
fccl <- rast("./data/reproj/final/FCCL_fullseries_GER_epsg3035.tif") 

# fractional cover of standing deadwood (epsg3035)
fcsd_2018 <- rast("./data/reproj/final/DE-2018_epgs3035.tif") 
fcsd_2019 <- rast("./data/reproj/final/DE-2019_epgs3035.tif")
fcsd_2020 <- rast("./data/reproj/final/DE-2020_epgs3035.tif")
fcsd_2021 <- rast("./data/reproj/final/DE-2021_epgs3035.tif")
fcsd <- list(fcsd_2018,fcsd_2019,fcsd_2020,fcsd_2021)

# forest structure (epsg4326)
agbd_2017 <- rast("./data/reproj/final/2017_agbd_epgs3035.tif")
agbd_2018 <- rast("./data/reproj/final/2018_agbd_epgs3035.tif")
agbd_2019 <- rast("./data/reproj/final/2019_agbd_epgs3035.tif")
agbd_2020 <- rast("./data/reproj/final/2020_agbd_epgs3035.tif")
agbd_2021 <- rast("./data/reproj/final/2021_agbd_epgs3035.tif")
agbd_2022 <- rast("./data/reproj/final/2022_agbd_epgs3035.tif")
agbd_2023 <- rast("./data/reproj/final/2023_agbd_epgs3035.tif")
agbd <- list(agbd_2017,agbd_2018,agbd_2019,agbd_2020,agbd_2021,agbd_2022,agbd_2023)

cover_2017 <- rast("./data/reproj/final/2017_cove_epgs3035.tif")
cover_2018 <- rast("./data/reproj/final/2018_cove_epgs3035.tif")
cover_2019 <- rast("./data/reproj/final/2019_cove_epgs3035.tif")
cover_2020 <- rast("./data/reproj/final/2020_cove_epgs3035.tif")
cover_2021 <- rast("./data/reproj/final/2021_cove_epgs3035.tif")
cover_2022 <- rast("./data/reproj/final/2022_cove_epgs3035.tif")
cover_2023 <- rast("./data/reproj/final/2023_cove_epgs3035.tif")
cover <- list(cover_2017,cover_2018,cover_2019,cover_2020,cover_2021,cover_2022,cover_2023)

fdh_normal_2017 <- rast("./data/reproj/final/2017_fhd__epgs3035.tif")
fdh_normal_2018 <- rast("./data/reproj/final/2018_fhd__epgs3035.tif")
fdh_normal_2019 <- rast("./data/reproj/final/2019_fhd__epgs3035.tif")
fdh_normal_2020 <- rast("./data/reproj/final/2020_fhd__epgs3035.tif")
fdh_normal_2021 <- rast("./data/reproj/final/2021_fhd__epgs3035.tif")
fdh_normal_2022 <- rast("./data/reproj/final/2022_fhd__epgs3035.tif")
fdh_normal_2023 <- rast("./data/reproj/final/2023_fhd__epgs3035.tif")
fdh_normal <- list(fdh_normal_2017,fdh_normal_2018,fdh_normal_2019,fdh_normal_2020,fdh_normal_2021,fdh_normal_2022,fdh_normal_2023)

rh_95_2017 <- rast("./data/reproj/final/2017_rh_9_epgs3035.tif")
rh_95_2018 <- rast("./data/reproj/final/2018_rh_9_epgs3035.tif")
rh_95_2019 <- rast("./data/reproj/final/2019_rh_9_epgs3035.tif")
rh_95_2020 <- rast("./data/reproj/final/2020_rh_9_epgs3035.tif")
rh_95_2021 <- rast("./data/reproj/final/2021_rh_9_epgs3035.tif")
rh_95_2022 <- rast("./data/reproj/final/2022_rh_9_epgs3035.tif")
rh_95_2023 <- rast("./data/reproj/final/2023_rh_9_epgs3035.tif")
rh_95 <- list(rh_95_2017,rh_95_2018,rh_95_2019,rh_95_2020,rh_95_2021,rh_95_2022,rh_95_2023)

# tree species (epsg3035)
tree_species_filled <- rast("./data/reproj/final/Thuenen_Germany_stockedForest_TreeSpecies_combined_epsg3035.tif")
stocked_forest <- rast("./data/reproj/final/Thuenen_Germany_stockedForest_2018_epsg3035.tif")


#### load shapefiles ####

# federal states of germany (epsg4326)
federal_4326 <- st_read("/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_AdministrativeBoundary/gadm41_DEU_1.json") 
# districts of germany (epsg4326)
districts_4326 <- st_read("/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_AdministrativeBoundary/gadm41_DEU_2.json") 
# hexagons of germany (epsg4326)
#hex_4326 <- st_read("./data/shp/ger_hex_res6_epsg4326.gpkg") 
hex_3035 <- st_read("./data/shp/ger_hex_epsg3035.gpkg")
# forest landscapes
forest_landscapes <- st_read("./data/shp/forstl_gl_2011.shp") 
# forest growth districts
forest_district <- st_read("/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_SilviculturalGrowthDistricts/wgwb_wb_2020.shp") 


# reproject to epsg3035
federal_3035 <- st_transform(federal_4326, 3035)
districts_3035 <- st_transform(districts_4326, 3035)
#hex_3035 <- st_transform(hex_4326, 3035)
forest_landscape_3035 <- st_transform(forest_landscapes, 3035)
forest_district_3035 <- st_transform(forest_district, 3035)


#### create dataframe ####

hex_3035$id <- seq(1,nrow(hex_3035))
hex_names <- hex_3035$id
fulltab <- data.frame(matrix(NA, nrow = 0, ncol = 92))
colnames(fulltab) <- c("year", "id", "district", "federal", "forst_landscape", "forest_district", "forest_area_ha", "pine_perc", "spruce_perc", "beech_perc", "oak_perc", "other_treespecies_perc", "dec_perc", "con_perc",
                       "agbd_true_min", "agbd_true_max", "agbd_q1", "agbd_median", "agbd_mean", "agbd_q3", "cover_std_dev", "agbd_iqr", "agbd_outlier_min", "agbd_outlier_max", "agbd_true_range",
                       "cover_true_min", "cover_true_max", "cover_q1", "cover_median", "cover_mean", "cover_q3", "cover_std_dev", "cover_iqr", "cover_outlier_min", "cover_outlier_max", "cover_true_range",
                       "fdh_normal_true_min", "fdh_normal_true_max", "fdh_normal_q1", "fdh_normal_median", "fdh_normal_mean", "fdh_normal_q3", "cover_std_dev", "fdh_normal_iqr", "fdh_normal_outlier_min", "fdh_normal_outlier_max", "fdh_normal_true_range",
                       "rh_95_true_min", "rh_95_true_max", "rh_95_q1", "rh_95_median", "rh_95_mean", "rh_95_q3", "cover_std_dev", "rh_95_iqr", "rh_95_outlier_min", "rh_95_outlier_max", "rh_95_true_range",
                       "fcsd_true_min", "fcsd_true_max", "fcsd_q1", "fcsd_median", "fcsd_mean", "fcsd_q3", "cover_std_dev", "fcsd_iqr", "fcsd_outlier_min", "fcsd_outlier_max", "fcsd_true_range",
                       "fccl_ges_perc", "fccl_ges_ha", "fccl_ges_cum_perc", "fccl_ges_cum_ha",
                       "fccl_jan_perc", "fccl_feb_perc", "fccl_mar_perc", "fccl_apr_perc", "fccl_mai_perc", "fccl_jun_perc", "fccl_jul_perc", "fccl_aug_perc", "fccl_sep_perc", "fccl_okt_perc", "fccl_nov_perc", "fccl_dez_perc",
                       "pine_fccl", "spruce_fccl", "beech_fccl", "oak_fccl", "others_fccl", "dec_fccl", "con_fccl")



# Define the function to run the loop for each subset of hex_names
process_hexagon <- function(hex_name) {
  
  # get pixels per hexagon  
  sub_3035 <- subset(hex_3035,hex_3035$id==hex_name)

  ## get values for time-independent parameters
  sub_district <- head(st_intersection(districts_3035, sub_3035),1)
  district <- as.character(sub_district$NAME_2)
  federal <- as.character(sub_district$NAME_1)
  fl <- head(st_intersection(forest_landscape_3035, sub_3035),1) # take first entry of intersection
  fl <- fl$gl_name
  fl <- ifelse(length(fl) == 0, NA, fl) # ifempty character, fill with NA
  fd <- head(st_intersection(forest_district_3035, sub_3035),1)
  fd <- fd$bez_wg_bu
  fd <- ifelse(length(fd) == 0, NA, fd) 
  
  # area of stocked forest
  stocked_forest_values <- na.omit(values(mask(crop(stocked_forest, ext(sub_3035)), mask=sub_3035)))
  forest_px <- length(stocked_forest_values)   # number of stocked forest pixels
  forest_area_ha <- length(stocked_forest_values)/100   # forest area in ha as pixels are 10mx10m
  
  # fccl
  fccl_cropped <- mask(crop(fccl, ext(sub_3035)), mask=sub_3035)
  fccl_cropped_values <- na.omit(values(fccl_cropped))
  fccl_cropped_values <- fccl_cropped_values[fccl_cropped_values!=0] # remove non-forest pixels
  
  # tree species
  ts_cropped <- mask(crop(tree_species_filled, ext(sub_3035)), mask=sub_3035)
  unclass_cropped_values <- global(ts_cropped == 1, "sum", na.rm=T)$sum
  birch_cropped_values <- global(ts_cropped == 2, "sum", na.rm=T)$sum
  beech_cropped_values <- global(ts_cropped == 3, "sum", na.rm=T)$sum
  douglasfir_cropped_values <- global(ts_cropped == 4, "sum", na.rm=T)$sum
  oak_cropped_values <- global(ts_cropped == 5, "sum", na.rm=T)$sum
  alder_cropped_values <- global(ts_cropped == 6, "sum", na.rm=T)$sum
  spruce_cropped_values <- global(ts_cropped == 8, "sum", na.rm=T)$sum
  pine_cropped_values <- global(ts_cropped == 9, "sum", na.rm=T)$sum
  larch_cropped_values <- global(ts_cropped == 10, "sum", na.rm=T)$sum
  fir_cropped_values <- global(ts_cropped == 14, "sum", na.rm=T)$sum
  other_dec_cropped_values <- global(ts_cropped %in% c(16,17,20), "sum", na.rm=T)$sum
  other_con_cropped_values <- global(ts_cropped == 21, "sum", na.rm=T)$sum
  
  unclass_values = sum(unclass_cropped_values, birch_cropped_values, douglasfir_cropped_values, alder_cropped_values, larch_cropped_values, fir_cropped_values, other_dec_cropped_values, other_con_cropped_values)
  dec_cropped_values = sum(birch_cropped_values,beech_cropped_values,oak_cropped_values,alder_cropped_values,other_dec_cropped_values)
  con_cropped_values = sum(douglasfir_cropped_values,spruce_cropped_values,pine_cropped_values,larch_cropped_values,fir_cropped_values,other_con_cropped_values)
  
  unclass_per <- unclass_values / forest_px
  beech_per <- beech_cropped_values / forest_px
  oak_per <- oak_cropped_values / forest_px
  spruce_per <- spruce_cropped_values / forest_px
  pine_per <- pine_cropped_values / forest_px
  dec_per <- dec_cropped_values / forest_px
  con_per <- con_cropped_values / forest_px
  
  # Cross-tabulate fccl and tree species
  # Create a combined raster to analyze overlaps
  fccl_ts_combined <- c(fccl_cropped, ts_cropped)
  
  # Initialize result container for each hexagon
  hex_results <- data.frame(matrix(NA, nrow = 7, ncol = 92))
  colnames(hex_results) <- colnames(fulltab)
  
  fccl_ges_cum <- 0 # set cumulative couter to 0
  
  ## get values for time-dependent parameters
  for(j in 1:7){
    # agbd
    agbd_cropped_values <- na.omit(values(mask(crop(agbd[[j]], ext(sub_3035)), mask=sub_3035)))
    agbd_true_min <- min(agbd_cropped_values)
    agbd_true_max <- max(agbd_cropped_values)
    agbd_q1 <- quantile(agbd_cropped_values, 0.25)
    agbd_median <- median(agbd_cropped_values)
    agbd_mean <- mean(agbd_cropped_values)
    agbd_q3 <- quantile(agbd_cropped_values, 0.75)
    agbd_std_dev <- sd(agbd_cropped_values)
    agbd_iqr <- IQR(agbd_cropped_values)
    agbd_outlier_min <- agbd_q1 - 1.5 * agbd_iqr
    agbd_outlier_max <- agbd_q3 + 1.5 * agbd_iqr
    agbd_true_range <- agbd_true_max - agbd_true_min
    
    # cover
    cover_cropped_values <- na.omit(values(mask(crop(cover[[j]], ext(sub_3035)), mask=sub_3035)))
    cover_true_min <- min(cover_cropped_values)
    cover_true_max <- max(cover_cropped_values)
    cover_q1 <- quantile(cover_cropped_values, 0.25)
    cover_median <- median(cover_cropped_values)
    cover_mean <- mean(cover_cropped_values)
    cover_q3 <- quantile(cover_cropped_values, 0.75)
    cover_std_dev <- sd(cover_cropped_values)
    cover_iqr <- IQR(cover_cropped_values)
    cover_outlier_min <- cover_q1 - 1.5 * cover_iqr
    cover_outlier_max <- cover_q3 + 1.5 * cover_iqr
    cover_true_range <- cover_true_max - cover_true_min
    
    # fdh_normal
    fdh_normal_cropped_values <- na.omit(values(mask(crop(fdh_normal[[j]], ext(sub_3035)), mask=sub_3035)))
    fdh_normal_true_min <- min(fdh_normal_cropped_values)
    fdh_normal_true_max <- max(fdh_normal_cropped_values)
    fdh_normal_q1 <- quantile(fdh_normal_cropped_values, 0.25)
    fdh_normal_median <- median(fdh_normal_cropped_values)
    fdh_normal_mean <- mean(fdh_normal_cropped_values)
    fdh_normal_q3 <- quantile(fdh_normal_cropped_values, 0.75)
    fdh_normal_std_dev <- sd(fdh_normal_cropped_values)
    fdh_normal_iqr <- IQR(fdh_normal_cropped_values)
    fdh_normal_outlier_min <- fdh_normal_q1 - 1.5 * fdh_normal_iqr
    fdh_normal_outlier_max <- fdh_normal_q3 + 1.5 * fdh_normal_iqr
    fdh_normal_true_range <- fdh_normal_true_max - fdh_normal_true_min
    
    # rh_95
    rh_95_cropped_values <- na.omit(values(mask(crop(rh_95[[j]], ext(sub_3035)), mask=sub_3035)))
    rh_95_true_min <- min(rh_95_cropped_values)
    rh_95_true_max <- max(rh_95_cropped_values)
    rh_95_q1 <- quantile(rh_95_cropped_values, 0.25)
    rh_95_median <- median(rh_95_cropped_values)
    rh_95_mean <- mean(rh_95_cropped_values)
    rh_95_q3 <- quantile(rh_95_cropped_values, 0.75)
    rh_95_std_dev <- sd(rh_95_cropped_values)
    rh_95_iqr <- IQR(rh_95_cropped_values)
    rh_95_outlier_min <- rh_95_q1 - 1.5 * rh_95_iqr
    rh_95_outlier_max <- rh_95_q3 + 1.5 * rh_95_iqr
    rh_95_true_range <- rh_95_true_max - rh_95_true_min
    
    # fcsd
    if(j>=2 & j<=5){
      #fcsd
      fcsd_cropped_values <- na.omit(values(mask(crop(fcsd[[j-1]], ext(sub_3035)), mask=sub_3035)))
      fcsd_true_min <- min(fcsd_cropped_values)
      fcsd_true_max <- max(fcsd_cropped_values)
      fcsd_q1 <- quantile(fcsd_cropped_values, 0.25)
      fcsd_median <- median(fcsd_cropped_values)
      fcsd_mean <- mean(fcsd_cropped_values)
      fcsd_q3 <- quantile(fcsd_cropped_values, 0.75)
      fcsd_std_dev <- sd(fcsd_cropped_values)
      fcsd_iqr <- IQR(fcsd_cropped_values)
      fcsd_outlier_min <- fcsd_q1 - 1.5 * fcsd_iqr
      fcsd_outlier_max <- fcsd_q3 + 1.5 * fcsd_iqr
      fcsd_true_range <- fcsd_true_max - fcsd_true_min  
    } else{
      fcsd_true_min <- NA
      fcsd_true_max <- NA
      fcsd_q1 <- NA
      fcsd_median <- NA
      fcsd_mean <- NA
      fcsd_q3 <- NA
      fcsd_std_dev <- NA
      fcsd_iqr <- NA
      fcsd_outlier_min <- NA
      fcsd_outlier_max <- NA
      fcsd_true_range <- NA  
    }
    
    # fccl per forest area
    if(j==1){
      fccl_jan <- NA
      fccl_feb <- NA
      fccl_mar <- NA
      fccl_apr <- NA
      fccl_mai <- NA
      fccl_jun <- NA
      fccl_jul <- NA
      fccl_aug <- NA
      fccl_sep <- sum(fccl_cropped_values == 1) / forest_px
      fccl_okt <- sum(fccl_cropped_values == 2) / forest_px
      fccl_nov <- sum(fccl_cropped_values == 3) / forest_px
      fccl_dez <- sum(fccl_cropped_values == 4) / forest_px
      
      # Count pixels per tree species for each fccl class
      pine_fccl <- global((fccl_ts_combined[[1]] %in% c(1:4)) & (fccl_ts_combined[[2]] == 9), "sum", na.rm = TRUE)[1, 1] / forest_px
      spruce_fccl <- global((fccl_ts_combined[[1]] %in% c(1:4)) & (fccl_ts_combined[[2]] == 8), "sum", na.rm = TRUE)[1, 1] / forest_px
      beech_fccl <- global((fccl_ts_combined[[1]] %in% c(1:4)) & (fccl_ts_combined[[2]] == 3), "sum", na.rm = TRUE)[1, 1] / forest_px
      oak_fccl <- global((fccl_ts_combined[[1]] %in% c(1:4)) & (fccl_ts_combined[[2]] == 5), "sum", na.rm = TRUE)[1, 1] / forest_px
      others_fccl <- global((fccl_ts_combined[[1]] %in% c(1:4)) & !(fccl_ts_combined[[2]] %in% c(9, 8, 3, 5)), "sum", na.rm = TRUE)[1, 1] / forest_px
      dec_fccl <- global((fccl_ts_combined[[1]] %in% c(1:4)) & (fccl_ts_combined[[2]] %in% c(2, 3, 5, 6, 16, 17, 20)), "sum", na.rm = TRUE)[1, 1] / forest_px
      con_fccl <- global((fccl_ts_combined[[1]] %in% c(1:4)) & (fccl_ts_combined[[2]] %in% c(4, 8, 9, 10, 14, 21)), "sum", na.rm = TRUE)[1, 1] / forest_px
      
    } else if(j==2){
      fccl_jan <- sum(fccl_cropped_values == 5) / forest_px
      fccl_feb <- sum(fccl_cropped_values == 6) / forest_px
      fccl_mar <- sum(fccl_cropped_values == 7) / forest_px
      fccl_apr <- sum(fccl_cropped_values == 8) / forest_px
      fccl_mai <- sum(fccl_cropped_values == 9) / forest_px
      fccl_jun <- sum(fccl_cropped_values == 10) / forest_px
      fccl_jul <- sum(fccl_cropped_values == 11) / forest_px
      fccl_aug <- sum(fccl_cropped_values == 12) / forest_px
      fccl_sep <- sum(fccl_cropped_values == 13) / forest_px
      fccl_okt <- sum(fccl_cropped_values == 14) / forest_px
      fccl_nov <- sum(fccl_cropped_values == 15) / forest_px
      fccl_dez <- sum(fccl_cropped_values == 16) / forest_px
      
      # Count pixels per tree species for each fccl class
      pine_fccl <- global((fccl_ts_combined[[1]] %in% c(5:16)) & (fccl_ts_combined[[2]] == 9), "sum", na.rm = TRUE)[1, 1] / forest_px
      spruce_fccl <- global((fccl_ts_combined[[1]] %in% c(5:16)) & (fccl_ts_combined[[2]] == 8), "sum", na.rm = TRUE)[1, 1] / forest_px
      beech_fccl <- global((fccl_ts_combined[[1]] %in% c(5:16)) & (fccl_ts_combined[[2]] == 3), "sum", na.rm = TRUE)[1, 1] / forest_px
      oak_fccl <- global((fccl_ts_combined[[1]] %in% c(5:16)) & (fccl_ts_combined[[2]] == 5), "sum", na.rm = TRUE)[1, 1] / forest_px
      others_fccl <- global((fccl_ts_combined[[1]] %in% c(5:16)) & !(fccl_ts_combined[[2]] %in% c(9, 8, 3, 5)), "sum", na.rm = TRUE)[1, 1] / forest_px
      dec_fccl <- global((fccl_ts_combined[[1]] %in% c(5:16)) & (fccl_ts_combined[[2]] %in% c(2, 3, 5, 6, 16, 17, 20)), "sum", na.rm = TRUE)[1, 1] / forest_px
      con_fccl <- global((fccl_ts_combined[[1]] %in% c(5:16)) & (fccl_ts_combined[[2]] %in% c(4, 8, 9, 10, 14, 21)), "sum", na.rm = TRUE)[1, 1] / forest_px
      
    } else if(j==3){
      fccl_jan <- sum(fccl_cropped_values == 17) / forest_px
      fccl_feb <- sum(fccl_cropped_values == 18) / forest_px
      fccl_mar <- sum(fccl_cropped_values == 19) / forest_px
      fccl_apr <- sum(fccl_cropped_values == 20) / forest_px
      fccl_mai <- sum(fccl_cropped_values == 21) / forest_px
      fccl_jun <- sum(fccl_cropped_values == 22) / forest_px
      fccl_jul <- sum(fccl_cropped_values == 23) / forest_px
      fccl_aug <- sum(fccl_cropped_values == 24) / forest_px
      fccl_sep <- sum(fccl_cropped_values == 25) / forest_px
      fccl_okt <- sum(fccl_cropped_values == 26) / forest_px
      fccl_nov <- sum(fccl_cropped_values == 27) / forest_px
      fccl_dez <- sum(fccl_cropped_values == 28) / forest_px
      
      # Count pixels per tree species for each fccl class
      pine_fccl <- global((fccl_ts_combined[[1]] %in% c(17:28)) & (fccl_ts_combined[[2]] == 9), "sum", na.rm = TRUE)[1, 1] / forest_px
      spruce_fccl <- global((fccl_ts_combined[[1]] %in% c(17:28)) & (fccl_ts_combined[[2]] == 8), "sum", na.rm = TRUE)[1, 1] / forest_px
      beech_fccl <- global((fccl_ts_combined[[1]] %in% c(17:28)) & (fccl_ts_combined[[2]] == 3), "sum", na.rm = TRUE)[1, 1] / forest_px
      oak_fccl <- global((fccl_ts_combined[[1]] %in% c(17:28)) & (fccl_ts_combined[[2]] == 5), "sum", na.rm = TRUE)[1, 1] / forest_px
      others_fccl <- global((fccl_ts_combined[[1]] %in% c(17:28)) & !(fccl_ts_combined[[2]] %in% c(9, 8, 3, 5)), "sum", na.rm = TRUE)[1, 1] / forest_px
      dec_fccl <- global((fccl_ts_combined[[1]] %in% c(17:28)) & (fccl_ts_combined[[2]] %in% c(2, 3, 5, 6, 16, 17, 20)), "sum", na.rm = TRUE)[1, 1] / forest_px
      con_fccl <- global((fccl_ts_combined[[1]] %in% c(17:28)) & (fccl_ts_combined[[2]] %in% c(4, 8, 9, 10, 14, 21)), "sum", na.rm = TRUE)[1, 1] / forest_px
      
    } else if(j==4){
      fccl_jan <- sum(fccl_cropped_values == 29) / forest_px
      fccl_feb <- sum(fccl_cropped_values == 30) / forest_px
      fccl_mar <- sum(fccl_cropped_values == 31) / forest_px
      fccl_apr <- sum(fccl_cropped_values == 32) / forest_px
      fccl_mai <- sum(fccl_cropped_values == 33) / forest_px
      fccl_jun <- sum(fccl_cropped_values == 34) / forest_px
      fccl_jul <- sum(fccl_cropped_values == 35) / forest_px
      fccl_aug <- sum(fccl_cropped_values == 36) / forest_px
      fccl_sep <- sum(fccl_cropped_values == 37) / forest_px
      fccl_okt <- sum(fccl_cropped_values == 38) / forest_px
      fccl_nov <- sum(fccl_cropped_values == 39) / forest_px
      fccl_dez <- sum(fccl_cropped_values == 40) / forest_px
      
      # Count pixels per tree species for each fccl class
      pine_fccl <- global((fccl_ts_combined[[1]] %in% c(29:40)) & (fccl_ts_combined[[2]] == 9), "sum", na.rm = TRUE)[1, 1] / forest_px
      spruce_fccl <- global((fccl_ts_combined[[1]] %in% c(29:40)) & (fccl_ts_combined[[2]] == 8), "sum", na.rm = TRUE)[1, 1] / forest_px
      beech_fccl <- global((fccl_ts_combined[[1]] %in% c(29:40)) & (fccl_ts_combined[[2]] == 3), "sum", na.rm = TRUE)[1, 1] / forest_px
      oak_fccl <- global((fccl_ts_combined[[1]] %in% c(29:40)) & (fccl_ts_combined[[2]] == 5), "sum", na.rm = TRUE)[1, 1] / forest_px
      others_fccl <- global((fccl_ts_combined[[1]] %in% c(29:40)) & !(fccl_ts_combined[[2]] %in% c(9, 8, 3, 5)), "sum", na.rm = TRUE)[1, 1] / forest_px
      dec_fccl <- global((fccl_ts_combined[[1]] %in% c(29:40)) & (fccl_ts_combined[[2]] %in% c(2, 3, 5, 6, 16, 17, 20)), "sum", na.rm = TRUE)[1, 1] / forest_px
      con_fccl <- global((fccl_ts_combined[[1]] %in% c(29:40)) & (fccl_ts_combined[[2]] %in% c(4, 8, 9, 10, 14, 21)), "sum", na.rm = TRUE)[1, 1] / forest_px
      
    } else if(j==5){
      fccl_jan <- sum(fccl_cropped_values == 41) / forest_px
      fccl_feb <- sum(fccl_cropped_values == 42) / forest_px
      fccl_mar <- sum(fccl_cropped_values == 43) / forest_px
      fccl_apr <- sum(fccl_cropped_values == 44) / forest_px
      fccl_mai <- sum(fccl_cropped_values == 45) / forest_px
      fccl_jun <- sum(fccl_cropped_values == 46) / forest_px
      fccl_jul <- sum(fccl_cropped_values == 47) / forest_px
      fccl_aug <- sum(fccl_cropped_values == 48) / forest_px
      fccl_sep <- sum(fccl_cropped_values == 49) / forest_px
      fccl_okt <- sum(fccl_cropped_values == 50) / forest_px
      fccl_nov <- sum(fccl_cropped_values == 51) / forest_px
      fccl_dez <- sum(fccl_cropped_values == 52) / forest_px
      
      # Count pixels per tree species for each fccl class
      pine_fccl <- global((fccl_ts_combined[[1]] %in% c(41:52)) & (fccl_ts_combined[[2]] == 9), "sum", na.rm = TRUE)[1, 1] / forest_px
      spruce_fccl <- global((fccl_ts_combined[[1]] %in% c(41:52)) & (fccl_ts_combined[[2]] == 8), "sum", na.rm = TRUE)[1, 1] / forest_px
      beech_fccl <- global((fccl_ts_combined[[1]] %in% c(41:52)) & (fccl_ts_combined[[2]] == 3), "sum", na.rm = TRUE)[1, 1] / forest_px
      oak_fccl <- global((fccl_ts_combined[[1]] %in% c(41:52)) & (fccl_ts_combined[[2]] == 5), "sum", na.rm = TRUE)[1, 1] / forest_px
      others_fccl <- global((fccl_ts_combined[[1]] %in% c(41:52)) & !(fccl_ts_combined[[2]] %in% c(9, 8, 3, 5)), "sum", na.rm = TRUE)[1, 1] / forest_px
      dec_fccl <- global((fccl_ts_combined[[1]] %in% c(41:52)) & (fccl_ts_combined[[2]] %in% c(2, 3, 5, 6, 16, 17, 20)), "sum", na.rm = TRUE)[1, 1] / forest_px
      con_fccl <- global((fccl_ts_combined[[1]] %in% c(41:52)) & (fccl_ts_combined[[2]] %in% c(4, 8, 9, 10, 14, 21)), "sum", na.rm = TRUE)[1, 1] / forest_px
      
    } else if(j==6){
      fccl_jan <- sum(fccl_cropped_values == 53) / forest_px
      fccl_feb <- sum(fccl_cropped_values == 54) / forest_px
      fccl_mar <- sum(fccl_cropped_values == 55) / forest_px
      fccl_apr <- sum(fccl_cropped_values == 56) / forest_px
      fccl_mai <- sum(fccl_cropped_values == 57) / forest_px
      fccl_jun <- sum(fccl_cropped_values == 58) / forest_px
      fccl_jul <- sum(fccl_cropped_values == 59) / forest_px
      fccl_aug <- sum(fccl_cropped_values == 60) / forest_px
      fccl_sep <- sum(fccl_cropped_values == 61) / forest_px
      fccl_okt <- sum(fccl_cropped_values == 62) / forest_px
      fccl_nov <- sum(fccl_cropped_values == 63) / forest_px
      fccl_dez <- sum(fccl_cropped_values == 64) / forest_px
      
      # Count pixels per tree species for each fccl class
      pine_fccl <- global((fccl_ts_combined[[1]] %in% c(53:64)) & (fccl_ts_combined[[2]] == 9), "sum", na.rm = TRUE)[1, 1] / forest_px
      spruce_fccl <- global((fccl_ts_combined[[1]] %in% c(53:64)) & (fccl_ts_combined[[2]] == 8), "sum", na.rm = TRUE)[1, 1] / forest_px
      beech_fccl <- global((fccl_ts_combined[[1]] %in% c(53:64)) & (fccl_ts_combined[[2]] == 3), "sum", na.rm = TRUE)[1, 1] / forest_px
      oak_fccl <- global((fccl_ts_combined[[1]] %in% c(53:64)) & (fccl_ts_combined[[2]] == 5), "sum", na.rm = TRUE)[1, 1] / forest_px
      others_fccl <- global((fccl_ts_combined[[1]] %in% c(53:64)) & !(fccl_ts_combined[[2]] %in% c(9, 8, 3, 5)), "sum", na.rm = TRUE)[1, 1] / forest_px
      dec_fccl <- global((fccl_ts_combined[[1]] %in% c(53:64)) & (fccl_ts_combined[[2]] %in% c(2, 3, 5, 6, 16, 17, 20)), "sum", na.rm = TRUE)[1, 1] / forest_px
      con_fccl <- global((fccl_ts_combined[[1]] %in% c(53:64)) & (fccl_ts_combined[[2]] %in% c(4, 8, 9, 10, 14, 21)), "sum", na.rm = TRUE)[1, 1] / forest_px
      
    } else{
      fccl_jan <- sum(fccl_cropped_values == 65) / forest_px
      fccl_feb <- sum(fccl_cropped_values == 66) / forest_px
      fccl_mar <- sum(fccl_cropped_values == 67) / forest_px
      fccl_apr <- sum(fccl_cropped_values == 68) / forest_px
      fccl_mai <- sum(fccl_cropped_values == 69) / forest_px
      fccl_jun <- sum(fccl_cropped_values == 70) / forest_px
      fccl_jul <- sum(fccl_cropped_values == 71) / forest_px
      fccl_aug <- sum(fccl_cropped_values == 72) / forest_px
      fccl_sep <- sum(fccl_cropped_values == 73) / forest_px
      fccl_okt <- sum(fccl_cropped_values == 74) / forest_px
      fccl_nov <- sum(fccl_cropped_values == 75) / forest_px
      fccl_dez <- sum(fccl_cropped_values == 76) / forest_px
      
      # Count pixels per tree species for each fccl class
      pine_fccl <- global((fccl_ts_combined[[1]] %in% c(65:76)) & (fccl_ts_combined[[2]] == 9), "sum", na.rm = TRUE)[1, 1] / forest_px
      spruce_fccl <- global((fccl_ts_combined[[1]] %in% c(65:76)) & (fccl_ts_combined[[2]] == 8), "sum", na.rm = TRUE)[1, 1] / forest_px
      beech_fccl <- global((fccl_ts_combined[[1]] %in% c(65:76)) & (fccl_ts_combined[[2]] == 3), "sum", na.rm = TRUE)[1, 1] / forest_px
      oak_fccl <- global((fccl_ts_combined[[1]] %in% c(65:76)) & (fccl_ts_combined[[2]] == 5), "sum", na.rm = TRUE)[1, 1] / forest_px
      others_fccl <- global((fccl_ts_combined[[1]] %in% c(65:76)) & !(fccl_ts_combined[[2]] %in% c(9, 8, 3, 5)), "sum", na.rm = TRUE)[1, 1] / forest_px
      dec_fccl <- global((fccl_ts_combined[[1]] %in% c(65:76)) & (fccl_ts_combined[[2]] %in% c(2, 3, 5, 6, 16, 17, 20)), "sum", na.rm = TRUE)[1, 1] / forest_px
      con_fccl <- global((fccl_ts_combined[[1]] %in% c(65:76)) & (fccl_ts_combined[[2]] %in% c(4, 8, 9, 10, 14, 21)), "sum", na.rm = TRUE)[1, 1] / forest_px
    }
    
    # new fccl in % of stocked forest area in 2017
    fccl_ges <- sum(fccl_jan, fccl_feb, fccl_mar, fccl_apr, fccl_mai, fccl_jun, fccl_jul, fccl_aug, fccl_sep, fccl_okt, fccl_nov, fccl_dez, na.rm = TRUE)   
    
    # new fccl in ha
    fccl_ges_ha <- fccl_ges*forest_area_ha
    
    # cummulative fccl in % of stocked forest area in 2017
    fccl_ges_cum <- fccl_ges_cum + fccl_ges
    
    # cummulative fccl in ha of stocked forest area in 2017
    fccl_ges_cum_ha <- fccl_ges_cum*forest_area_ha
    
    
    # create table 
    # forest_area_ha: forest area in ha per hexagon based on stocked forest mask from 2017
    
    # Place computed values in row
    hex_results[j, ] <- c(
      2016 + j, hex_name, district, federal, fl, fd, forest_area_ha, pine_per, spruce_per, beech_per, oak_per, unclass_per, dec_per, con_per,
      agbd_true_min, agbd_true_max, agbd_q1, agbd_median, agbd_mean, agbd_q3, cover_std_dev,
      agbd_iqr, agbd_outlier_min, agbd_outlier_max, agbd_true_range,
      cover_true_min, cover_true_max, cover_q1, cover_median, cover_mean, cover_q3, cover_std_dev,
      cover_iqr, cover_outlier_min, cover_outlier_max, cover_true_range,
      fdh_normal_true_min, fdh_normal_true_max, fdh_normal_q1, fdh_normal_median, fdh_normal_mean,
      fdh_normal_q3, cover_std_dev, fdh_normal_iqr, fdh_normal_outlier_min, fdh_normal_outlier_max,
      fdh_normal_true_range, rh_95_true_min, rh_95_true_max, rh_95_q1, rh_95_median, rh_95_mean,
      rh_95_q3, cover_std_dev, rh_95_iqr, rh_95_outlier_min, rh_95_outlier_max, rh_95_true_range,
      fcsd_true_min, fcsd_true_max, fcsd_q1, fcsd_median, fcsd_mean, fcsd_q3, cover_std_dev,
      fcsd_iqr, fcsd_outlier_min, fcsd_outlier_max, fcsd_true_range, fccl_ges, fccl_ges_ha,
      fccl_ges_cum, fccl_ges_cum_ha, fccl_jan, fccl_feb, fccl_mar, fccl_apr, fccl_mai, fccl_jun,
      fccl_jul, fccl_aug, fccl_sep, fccl_okt, fccl_nov, fccl_dez,
      pine_fccl, spruce_fccl, beech_fccl, oak_fccl, others_fccl, dec_fccl, con_fccl
    )
  }
  return(hex_results)
}


results_list <- mclapply(hex_names, process_hexagon, mc.cores = num_cores)
fulltab <- do.call(rbind, results_list)
write.csv(fulltab, "./data/forestdata_fulltab_hex_3000ha_reproj.csv", row.names = FALSE)


## write data to geopackage

fulltab <- read.csv("./data/results/final_polygons/forestdata_fulltab_hex_3000ha_reproj.csv")
#hex_3035 <- hex_3035 %>% select("h3_address")
shp_re_full <- merge(hex_3035, fulltab, by.x="id", by.y="id", all.x=T, sort=F)

fulltab_2017 <- shp_re_full[shp_re_full$year == 2017, ]
fulltab_2018 <- shp_re_full[shp_re_full$year == 2018, ]
fulltab_2019 <- shp_re_full[shp_re_full$year == 2019, ]
fulltab_2020 <- shp_re_full[shp_re_full$year == 2020, ]
fulltab_2021 <- shp_re_full[shp_re_full$year == 2021, ]
fulltab_2022 <- shp_re_full[shp_re_full$year == 2022, ]
fulltab_2023 <- shp_re_full[shp_re_full$year == 2023, ]
st_write(fulltab_2017, "./data/results/final_polygons/fulltab_hex_res6_2017_hex_3000ha_reproj.gpkg", overwrite=T)
st_write(fulltab_2018, "./data/results/final_polygons/fulltab_hex_res6_2018_hex_3000ha_reproj.gpkg", overwrite=T)
st_write(fulltab_2019, "./data/results/final_polygons/fulltab_hex_res6_2019_hex_3000ha_reproj.gpkg", overwrite=T)
st_write(fulltab_2020, "./data/results/final_polygons/fulltab_hex_res6_2020_hex_3000ha_reproj.gpkg", overwrite=T)
st_write(fulltab_2021, "./data/results/final_polygons/fulltab_hex_res6_2021_hex_3000ha_reproj.gpkg", overwrite=T)
st_write(fulltab_2022, "./data/results/final_polygons/fulltab_hex_res6_2022_hex_3000ha_reproj.gpkg", overwrite=T)
st_write(fulltab_2023, "./data/results/final_polygons/fulltab_hex_res6_2023_hex_3000ha_reproj.gpkg", overwrite=T)



