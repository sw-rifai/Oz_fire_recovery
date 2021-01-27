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


dat <- arrow::read_parquet(file="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet")
coords_vi <- lazy_dt(dat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)


# proc nvis -----------------------------------------------------
nvis <- read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif") %>% 
  set_names('vc')
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
dem <- read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif")
dem <- st_extract(dem, coords_vi)
coords_vi





coords_vi["vc_name"] %>% plot
coords_vi$vc_name %>% table %>% as_tibble() %>% arrange(desc(n))
coords_vi %>% ggplot(data=.,aes(vc_name,fill=vc_name))+
  stat_count()+
  scale_x_discrete(guide=guide_axis(check.overlap = TRUE))


