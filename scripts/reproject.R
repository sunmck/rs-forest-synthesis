library(gdalUtilities)

target_crs = 3035

# Reproject the forest structure rasterfiles
file_paths <- list.files("/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_ForestStructure/", pattern = "\\.tif$")

for(file in file_paths){
  gdalwarp(srcfile = paste0("/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_ForestStructure/", file), 
           dstfile = paste0("./data/reproj/final/", substr(file, 1, 9), "_epgs3035.tif"),
           t_srs = paste0("EPSG:",target_crs), 
           r = "near",
           tr = c(10,10),
           overwrite = TRUE,
           co = "COMPRESS=DEFLATE",
           te = c(4008242,2684052,4702432,3561572))
  gc()
}

# Reproject the fcsd rasterfiles
file_paths <- list.files("/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_FractionalCoverStandingDeadwood/", pattern = "\\.tif$")

for(file in file_paths){
  gdalwarp(srcfile = paste0("/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_FractionalCoverStandingDeadwood/", file), 
           dstfile = paste0("./data/reproj/final/", substr(file, 1, 7), "_epgs3035.tif"),
           t_srs = paste0("EPSG:",target_crs), 
           r = "near",
           tr = c(10,10),
           overwrite = TRUE,
           co = "COMPRESS=DEFLATE",
           te = c(4008242,2684052,4702432,3561572))
  gc()
}

# Reproject stocked forest rasterfile
gdalwarp(srcfile = "/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_StockedForestArea/Thuenen_Germany_stockedForest_2018_Int16_epsg32632.tif", 
         dstfile = "./data/reproj/final/Thuenen_Germany_stockedForest_2018_epsg3035.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))


# Reproject fccl
gdalwarp(srcfile = "/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_TreeCanopyCoverLoss/FCCL_fullseries_GER.tif", 
         dstfile = "./data/reproj/final/FCCL_fullseries_GER_epsg3035.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))

# Reproject copernicus tree species
gdalwarp(srcfile = "/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_DominantTreeSpecies/Thuenen_Germany_TreeSpecies_201718_Int16_epsg4326.tif", 
         dstfile = "./data/reproj/final/Thuenen_Germany_TreeSpecies_201718_epsg3035.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))


# Reproject copernicus tree species
gdalwarp(srcfile = "/dss/dsstbyfs02/pn49ci/pn49ci-dss-0017/Germany_CopernicusHRL/ForestType/HRL_FTY2018_Germany_epsg4326.tif", 
         dstfile = "./data/reproj/final/HRL_FTY2018_Germany_epsg3035.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))

# Reproject copernicus tree species
gdalwarp(srcfile = "./data/reproj/final/Thuenen_Germany_stockedForest_TreeSpecies_combined_epsg3035.tif", 
         dstfile = "./data/reproj/final/Thuenen_Germany_stockedForest_TreeSpecies_combined_epsg3035_v2.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))


# Reproject change agents
gdalwarp(srcfile = "./data/reproj/disturbance_agent_aggregated_germany.tif", 
         dstfile = "./data/reproj/final/disturbance_agent_aggregated_germany_epsg3035.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))

# Reproject change agents
gdalwarp(srcfile = "./data/reproj/FCCL_2017_23_epsg3035_binary.tif", 
         dstfile = "./data/reproj/final/FCCL_2017_23_epsg3035_binary.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))

# Reproject size layer
gdalwarp(srcfile = "./data/reproj/FCCL_2020_area_m2_epsg3035_binary.tif", 
         dstfile = "./data/reproj/final/FCCL_2020_area_m2_epsg3035.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572))

# Reproject windthrow event hesse layer
gdalwarp(srcfile = "./frederieke_jan2018_ref_merged_rasterized.tif", 
         dstfile = "./data/reproj/final/frederieke_jan2018_ref_merged_rasterized.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572)
         )

gdalwarp(srcfile = "./data/reproj/frederieke_jan2018_ref_dam3_size_rasterized.tif", 
         dstfile = "./data/reproj/final/frederieke_jan2018_ref_dam3_size_rasterized.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572)
)

# Reproject bavarian forest layers
gdalwarp(srcfile = "./data/reproj/ref_bavarianforest_year.tif", 
         dstfile = "./data/reproj/final/ref_bavarianforest_year.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572)
)

gdalwarp(srcfile = "./data/reproj/ref_bavarianforest_removalcode.tif", 
         dstfile = "./data/reproj/final/ref_bavarianforest_removalcode.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572)
)

gdalwarp(srcfile = "./data/reproj/Bavarian_NP_2020_removed_size.tif", 
         dstfile = "./data/reproj/final/Bavarian_NP_2020_removed_size.tif",
         t_srs = paste0("EPSG:",target_crs), 
         r = "near",
         tr = c(10,10),
         overwrite = TRUE,
         co = "COMPRESS=DEFLATE",
         te = c(4008242,2684052,4702432,3561572)
)

