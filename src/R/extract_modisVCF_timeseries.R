# Load packages in this order (important)
library(tidyverse);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(arrow)

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



dat <- arrow::read_parquet(file="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet", 
                           col_select = c("x","y","id"))
coords_vi <- lazy_dt(dat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
# coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
# st_crs(coords_vi) <- st_crs(4326)
# coords_stars <- st_as_stars(coords_vi)


# proc trees ------------------------------------------------
tmp <- read_stars("../data_general/proc_data_Oz_fire_recovery/MOD44B_tree_500m_SE_coastal_2000_2019.tif", 
                  proxy = FALSE) %>% 
  set_names('tree_cover') %>% 
  st_set_dimensions(., 3, values=2000:2019, names='year')

# proc grass ------------------------------------------------
tmp2 <- read_stars("../data_general/proc_data_Oz_fire_recovery/MOD44B_grass_500m_SE_coastal_2000_2019.tif", 
                  proxy = FALSE) %>% 
  set_names('grass_cover') %>% 
  st_set_dimensions(., 3, values=2000:2019, names='year')

tmp <- c(tmp,tmp2)
tmp <- as.data.frame(tmp)
setDT(tmp)

tmp <- merge(tmp,coords_vi,by=c("x","y"))
arrow::write_parquet(tmp, sink="../data_general/proc_data_Oz_fire_recovery/MOD44B_tree_grass_500m_SE_coastal_2001_2019.parquet")



# plot(tmp[,,,5],col=viridis::viridis(5),breaks = 'equal')
# dev.new()
# plot(tmp2[,,,5],col=viridis::viridis(100),breaks = 'equal')
