pacman::p_load(tidyverse, data.table, lubridate,stars, sf,magick)

# Load data ------------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 1000) %>% 
  select(NAME_1)

library(h2o) 
h2o.init()
# top_mod <- h2o.loadModel("../data_general/proc_data_Oz_fire_recovery/xgboost_ttr5-lai-ocat_2021-06-27_/GBM_5_AutoML_20210627_093118")
top_mod <- h2o.loadModel("outputs/GBM_bestModNoSpecies_2022-03-11 14:47:20/GBM_model_R_1646970170483_1")


dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
dpreds <- bind_cols(h2o.predict(top_mod, as.h2o(dpreds)) %>% as.data.table(),
                    dpreds)

dout <- dpreds %>% select(x,y,predict)

dpreds <- dpreds %>% 
  mutate(
    x=round(x*200)/200, 
    y=round(y*200)/200) %>% 
  .[,.(predict = mean(predict,na.rm=T)), by=.(x,y)] 

unique(dpreds$x) %>% sort %>% diff
unique(dpreds$y) %>% sort %>% diff

dfull <- expand_grid(x = seq(min(dpreds$x),max(dpreds$x),by=0.005),
  y= seq(min(dpreds$y),max(dpreds$y),by=0.005)) %>% 
  as.data.table() %>% 
  merge(., dpreds,by =c("x","y"), all.x=T,all.y=T)

sout <- st_as_stars(dpreds,
    dims=c("x","y"), 
    y_decreasing=T,
    ) %>% st_set_crs(4326)

tmp <- raster::rasterFromXYZ(dpreds)
tmp <- tmp %>% st_as_stars() %>% st_set_crs(4326) 
tmp %>% plot(col=viridis::viridis(1000),breaks='equal')
stars::write_stars(tmp, 
  dsn="../data_general/proc_data_Oz_fire_recovery/predicted_TTR_GBM_bestModNoSpecies_2022-03-11.tif")


# 

# tmp <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/GEDI-AUS-v27_500m_SE_coastal.tif")
# 
# gc(full=T)
# sout2 <- st_warp(sout, dest = tmp, use_gdal = T, options=list(no_data_value=NA_real_))
# 
# 
# 
# sout
# # %>% 
# #   st_set_crs(4326)
# st_raster_type(sout)
# st_get_dimension_values(sout,1) %>% diff %>% unique
# st_get_dimension_values(sout,2) %>% diff
# plot(sout,breaks='equal',col=viridis::viridis(100,option = 'B'))
# 
# 
# gc(full=T)
# test <- st_warp(sout, cellsize=0.005, crs=4326, use_gdal = F, no_data_value=NA_real_)
# 
# 
# 
# sout
# plot(sout)
# 
# st_warp(sout, dx=0.004491576, dy=0.004491576, use_gdal = F, crs=4326)
# 
# sout %>% st_raster_type()
# 
# st_warp(sout, dx=0.004491576, dy=0.004491576, use_gdal = F, crs=4326)
# 
# st_get_dimension_values(sout,1) %>% diff %>% unique
# 
# dfull$x %>% unique() %>% sort %>% diff %>% unique
# 
# dpreds %>% 
#   ggplot(data=.,aes(x,y,fill=predict))+
#   geom_tile()+
#   scale_fill_viridis_c(option='B')+
#   coord_sf(crs=4326)
