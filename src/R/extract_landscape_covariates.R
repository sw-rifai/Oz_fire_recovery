# Load packages in this order (important)
library(phenofit);
library(tidyverse);
library(usethis);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(RcppArmadillo)
library(arrow)
source("src/R/functions_time_to_recover.R")

# Load data ---------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)


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

vec_eof <- nvis[vc_name=="Eucalypt Open Forests"]
vec_ew <- nvis[vc_name=="Eucalypt Woodlands"]
vec_etof <- nvis[vc_name=="Eucalypt Tall Open Forests"]
vec_eow <- nvis[vc_name=="Eucalypt Open Woodlands"]
vec_cfw <- nvis[vc_name=="Callitris Forests and Woodlands"]
vec_rv <- nvis[vc_name=='Rainforests and Vine Thickets']

save(nvis, vec_eof,vec_ew, vec_etof,vec_eow,vec_cfw,vec_rv, 
        file='outputs/pixel_vegClass_groups.rds')

# proc elevation ------------------------------------------------
tmp <- read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif") %>% 
  set_names('elevation')
dem <- st_extract(tmp, coords_vi)
dem$id <- coords_vi$id
dem <- as.data.table(dem) %>% select(-geometry)

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

tmp1 <- merge(nvis,dem,by='id')
tmp2 <- merge(tmp1,tpi,by='id')
tmp3 <- merge(tmp2, d_soil, by='id')
write_parquet(tmp3, sink='../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet')



# soil[,,,str_which(names_soil,"PH_")]
# 
# plot(pto, col=viridis::inferno(100),breaks='equal')
# 
# names_soil[str_detect(names_soil,"PTO_")]
# 
# 
# des %>% plot
# tmp <- st_extract(soil, coords_vi)
# as.data.table(tmp)
# 
# st_as_stars(tmp)
# slga$id <- coords_vi$id
# slga <- as.data.table(slga) %>% select(-geometry)
# 
# 
# coords_vi["vc_name"] %>% plot
# coords_vi$vc_name %>% table %>% as_tibble() %>% arrange(desc(n))
# coords_vi %>% ggplot(data=.,aes(vc_name,fill=vc_name))+
#   stat_count()+
#   scale_x_discrete(guide=guide_axis(check.overlap = TRUE))
# 
# 
