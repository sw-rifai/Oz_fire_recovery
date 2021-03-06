# Load packages in this order (important)
library(tidyverse);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(arrow)
# source("src/R/functions_time_to_recover.R")

# STAGE 1: Generate reference grid ---------------------------------------------------
tmp1 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000000000.tif", 
                          proxy = F, RasterIO=list(bands=1)) %>% set_names("nir")
tmp2 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000001536.tif", 
                          proxy = F, RasterIO=list(bands=1)) %>% set_names("nir") 
tmp3 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000001536-0000000000.tif", 
                          proxy = F, RasterIO=list(bands=1))  %>% set_names("nir")
tmp4 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000001536-0000001536.tif", 
                          proxy = F, RasterIO=list(bands=1))  %>% set_names("nir")
gc(full=TRUE)
grid <- st_mosaic(tmp1,tmp2,tmp3,tmp4)
rm(tmp1,tmp2,tmp3,tmp4); gc(full=TRUE)


# # Load data ---------------------------------------------------------------
# oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
#   sf::st_simplify(., dTolerance = 0.1) %>% 
#   select(NAME_1)


dat <- arrow::read_parquet(file="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet")
coords_vi <- lazy_dt(dat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_stars <- st_as_stars(coords_vi)

# proc nvis -----------------------------------------------------
nvis <- read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif") %>%
  set_names('vc')
dem <- stars::read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif")
nvis <- st_warp(nvis,dest=dem)
tmp <- st_extract(nvis, coords_vi)
codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt",
                         fwf_widths(c(2,100)), skip = 1) %>%
  set_names(c("vc","vc_name")) %>%
  mutate(vc_name = factor(vc_name))
tmp <- inner_join(tmp,codes,by='vc')
coords_vi$vc <- tmp$vc
coords_vi$vc_name <- tmp$vc_name
nvis <- coords_vi %>%
  select(id,vc,vc_name) %>%
  as.data.table() %>%
  select(-geometry)
# 
# vec_eof <- nvis[vc_name=="Eucalypt Open Forests"]
# vec_ew <- nvis[vc_name=="Eucalypt Woodlands"]
# vec_etof <- nvis[vc_name=="Eucalypt Tall Open Forests"]
# vec_eow <- nvis[vc_name=="Eucalypt Open Woodlands"]
# vec_cfw <- nvis[vc_name=="Callitris Forests and Woodlands"]
# vec_rv <- nvis[vc_name=='Rainforests and Vine Thickets']
# 
# save(nvis, vec_eof,vec_ew, vec_etof,vec_eow,vec_cfw,vec_rv, 
#         file='outputs/pixel_vegClass_groups.rds')

codes2 <- readr::read_fwf("../data_general/NVIS/nvis_majorVegSubClasses.txt",
                         fwf_widths(c(2,100)), skip = 1) %>%
  set_names(c("vc2","vc2_name")) %>% 
  mutate(vc2_name = factor(vc2_name))
nvis2 <- stars::read_stars("../data_general/NVIS/nvis6_0p001.tif", 
                          proxy=F)
st_crs(nvis2) <- st_crs(4326)
names(nvis2) <- 'vc2'
# grid_points <- st_as_sf(st_coordinates(grid), coords = c("x","y"))
# st_crs(grid_points) <- st_crs(grid)
tmp_nvis2 <- st_extract(nvis2, coords_vi)
tmp_nvis2 <- left_join(tmp_nvis2,codes2,by='vc2')
coords_vi$vc2 <- tmp_nvis2$vc2
coords_vi$vc2_name <- tmp_nvis2$vc2_name
nvis2 <- coords_vi %>% 
  select(id,vc2,vc2_name) %>% 
  as.data.table() %>% 
  select(-geometry)


# proc elevation ------------------------------------------------
tmp <- read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif") %>% 
  set_names('elevation')
dem <- st_extract(tmp, coords_vi)
dem$id <- coords_vi$id
dem <- as.data.table(dem) %>% select(-geometry)

# proc slope ------------------------------------------------
tmp <- read_stars("../data_general/Oz_misc_data/slope_wwfhydrosheds_500m_SE_coastal.tif") %>% 
  set_names('slope')
slope <- st_extract(tmp, coords_vi)
slope$id <- coords_vi$id
slope <- as.data.table(slope) %>% select(-geometry)

# proc aspect ------------------------------------------------
tmp <- read_stars("../data_general/Oz_misc_data/aspect_wwfhydrosheds_500m_SE_coastal.tif") %>% 
  set_names('aspect')
aspect <- st_extract(tmp, coords_vi)
aspect$id <- coords_vi$id
aspect <- as.data.table(aspect) %>% select(-geometry)


# proc TPI ------------------------------------------------
tmp <- read_stars("../data_general/Oz_misc_data/mtpi_alos_500m_SE_coastal.tif") %>% 
  set_names('tpi')
tpi <- st_extract(tmp, coords_vi)
tpi$id <- coords_vi$id
tpi <- as.data.table(tpi) %>% select(-geometry)

# proc slga ------------------------------------------------
names_soil <- raster::stack("../data_general/Oz_misc_data/SLGA_500m_SE_coastal.tif") %>% 
  names()
soil <- read_stars("../data_general/Oz_misc_data/SLGA_500m_SE_coastal.tif",
                   proxy=FALSE) %>% 
  st_set_dimensions(.,3,values=names_soil)
# der <- st_extract(soil[,,,str_which(names_soil,'DER_DER_000_999_EV')],coords_vi)
# des <- st_extract(soil[,,,str_which(names_soil,'DES_DES_000_200_EV')],coords_vi)

der <- soil[,,,str_which(names_soil,'DER_DER_000_999_EV')] %>% st_apply(.,c("x","y"),mean)
des <- soil[,,,str_which(names_soil,'DES_DES_000_200_EV')] %>% st_apply(.,c("x","y"),mean)
pto <- st_apply(soil[,,,str_which(names_soil,"PTO_")],c("x","y"),mean)
sand <- st_apply(soil[,,,str_which(names_soil,"SND_")],c("x","y"),mean)
pH <- st_apply(soil[,,,str_which(names_soil,"pHc_")],c("x","y"),mean)
clay <- st_apply(soil[,,,str_which(names_soil,"CLY_")],c("x","y"),mean)
silt <- st_apply(soil[,,,str_which(names_soil,"SLT_")],c("x","y"),mean)
bdw <- st_apply(soil[,,,str_which(names_soil,"BDW_")],c("x","y"),mean)
soc <- st_apply(soil[,,,str_which(names_soil,"SOC_")],c("x","y"),mean)
awc <- st_apply(soil[,,,str_which(names_soil,"AWC_")],c("x","y"),mean)
nto <- st_apply(soil[,,,str_which(names_soil,"NTO_")],c("x","y"),mean)
ece <- st_apply(soil[,,,str_which(names_soil,"ECE_")],c("x","y"),mean)

names(der) <- 'der'
names(des) <- 'des'
names(pto) <- 'pto'
names(sand) <- 'sand'
names(pH) <- 'pH'
names(clay) <- 'clay'
names(silt) <- 'silt'
names(bdw) <- 'bdw'
names(soc) <- 'soc'
names(awc) <- 'awc'
names(nto) <- 'nto'
names(ece) <- 'ece'
d_soil <- st_extract(c(der,des,pto,sand,pH,clay,silt,bdw,soc,awc,nto,ece), coords_vi)
d_soil$id <- coords_vi$id
d_soil <- d_soil %>% st_drop_geometry()
d_soil <- d_soil %>% as.data.table()



nvis <- merge(nvis,nvis2,by='id')
tmp1 <- merge(nvis,dem,by='id')
tmp1_1 <- merge(slope,aspect,by='id')
tmp1 <- merge(tmp1,tmp1_1,by='id'); rm(tmp1_1); gc(full=TRUE)
tmp2 <- merge(tmp1,tpi,by='id')
tmp3 <- merge(tmp2, d_soil, by='id')
write_parquet(tmp3, sink='../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet')



