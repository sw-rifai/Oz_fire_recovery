library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
# setDTthreads(threads = 16)

ref_grid <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/grid_500m_SE_coastal.tif")
tmp1 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD11A1_C6_LST_500m_SE_coastal_2002-7_2021-4-0000000000-0000000000.tif", 
                          proxy=F) %>% 
  set_names("lsta") %>% 
  st_warp(., ref_grid, use_gdal = F) %>% 
  st_set_dimensions(., 3, values=seq(ymd("2002-07-01"),ymd("2021-03-31"),by='1 month'), names='date') %>% 
  as.data.table() %>% 
  .[is.na(lsta)==FALSE]; gc(full=TRUE)

tmp2 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD11A1_C6_LST_500m_SE_coastal_2002-7_2021-4-0000000000-0000001792.tif", 
  proxy=F) %>% 
  set_names("lsta") %>% 
  st_warp(., ref_grid, use_gdal = F) %>% 
  st_set_dimensions(., 3, values=seq(ymd("2002-07-01"),ymd("2021-03-31"),by='1 month'), names='date')%>% 
  as.data.table() %>% 
  .[is.na(lsta)==FALSE]; gc(full=TRUE)

tmp3 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD11A1_C6_LST_500m_SE_coastal_2002-7_2021-4-0000001792-0000000000.tif", 
  proxy=F) %>% 
  set_names("lsta") %>% 
  st_warp(., ref_grid, use_gdal = F) %>% 
  st_set_dimensions(., 3, values=seq(ymd("2002-07-01"),ymd("2021-03-31"),by='1 month'), names='date')%>% 
  as.data.table() %>% 
  .[is.na(lsta)==FALSE]; gc(full=TRUE)

tmp4 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD11A1_C6_LST_500m_SE_coastal_2002-7_2021-4-0000001792-0000001792.tif", 
  proxy=F) %>% 
  set_names("lsta") %>% 
  st_warp(., ref_grid, use_gdal = F) %>% 
  st_set_dimensions(., 3, values=seq(ymd("2002-07-01"),ymd("2021-03-31"),by='1 month'), names='date')%>% 
  as.data.table() %>% 
  .[is.na(lsta)==FALSE]; gc(full=TRUE)
gc(full=TRUE)
rm(ref_grid); gc()
gc(full=T)
lst_out <- rbindlist(list(tmp1,tmp2,tmp3,tmp4), use.names = TRUE)
gc(full=TRUE)
arrow::write_parquet(lst_out, 
                     sink="../data_general/proc_data_Oz_fire_recovery/MOD11A1_C6_LST_500m_SE_coastal_2002-7_2021-4_.parquet",
                     compression = 'snappy')

