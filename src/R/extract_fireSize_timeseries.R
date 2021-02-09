library(tidyverse);
library(stars);
library(data.table); 
library(dtplyr); 
library(lubridate) # load AFTER data.table
library(arrow); 


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


# OLD GRID - don't think it worked
# grid <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")
# setDT(grid)
# grid <- grid[date==min(date)] %>% select(x,y,id)
# grid <- st_as_stars(grid)
# st_crs(grid) <- st_crs(4326)
# gc(full=TRUE); 




cc_area1 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000000000-0000000000.tif",
                         proxy=FALSE) %>% 
  slice('band', seq(1,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('ba_m2') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date') %>% 
  st_warp(., dest=grid)
gc(full=TRUE)
cc_label1 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000000000-0000000000.tif",
                              proxy=FALSE) %>% 
  slice('band', seq(2,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('label') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date')%>% 
  st_warp(., dest=grid)
gc(full=TRUE)
cc1 <- c(cc_area1,cc_label1) 
rm(cc_area1,cc_label1); gc(full=TRUE)
cc1 <- as_tibble(cc1)
gc(full=TRUE)
cc1 <- cc1 %>% filter(is.na(ba_m2)==FALSE)
gc(full=TRUE)
cc1 <- cc1 %>% as.data.table() 
gc(full=TRUE)

cc_area2 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000000000-0000001280.tif",
                              proxy=FALSE) %>% 
  slice('band', seq(1,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('ba_m2') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date')%>% 
  st_warp(., dest=grid)
gc(full=TRUE)
cc_label2 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000000000-0000001280.tif",
                               proxy=FALSE) %>% 
  slice('band', seq(2,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('label') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date')%>% 
  st_warp(., dest=grid)
gc(full=TRUE)
cc2 <- c(cc_area2,cc_label2) %>% as_tibble() %>% filter(is.na(ba_m2)==FALSE)
rm(cc_area2,cc_label2); gc(full=TRUE)
cc2 <- as.data.table(cc2)
gc(full=TRUE)

cc_area3 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000001280-0000000000.tif",
                              proxy=FALSE) %>% 
  slice('band', seq(1,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('ba_m2') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date')%>% 
  st_warp(., dest=grid)

cc_label3 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000001280-0000000000.tif",
                               proxy=FALSE) %>% 
  slice('band', seq(2,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('label') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date')%>% 
  st_warp(., dest=grid)

cc3 <- c(cc_area3,cc_label3) %>% as_tibble() %>% filter(is.na(ba_m2)==FALSE)
rm(cc_area3,cc_label3); gc(full=TRUE)
cc3 <- as.data.table(cc3)
gc(full=TRUE)


cc_area4 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000001280-0000001280.tif",
                              proxy=FALSE) %>% 
  slice('band', seq(1,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('ba_m2') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date')%>% 
  st_warp(., dest=grid)

cc_label4 <- stars::read_stars("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000001280-0000001280.tif",
                               proxy=FALSE) %>% 
  slice('band', seq(2,by=2,length.out = dim(.)[3]/2)) %>% 
  set_names('label') %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),ymd("2020-12-01"),by='1 month'), 
                    names = 'date')%>% 
  st_warp(., dest=grid)
cc4 <- c(cc_area4,cc_label4) %>% as_tibble() %>% filter(is.na(ba_m2)==FALSE)
rm(cc_area4,cc_label4); gc(full=TRUE)
cc4 <- as.data.table(cc4)
gc(full=TRUE)

cc <- rbindlist(list(cc1,cc2,cc3,cc4),use.names = TRUE)
arrow::write_parquet(cc, sink='outputs/mcd64_conComp_2001_2020.parquet')


# unique(cc[,.(date,ba_m2,label)]) %>% 
#   ggplot(data=.,aes(ba_m2))+
#   geom_histogram(bins=100)+
#   scale_x_continuous(trans='log10')
# 
# cc %>% lazy_dt() %>% 
#   mutate(hydro_year=year(date+months(3))) %>% 
#   group_by(hydro_year) %>% 
#   summarize(val = sum(ba_m2,na.rm=TRUE)/sum(ba_m2>0,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(hydro_year, val))+
#   geom_line()+
#   geom_point()
# 
# 
# cc1 %>% tibble %>% 
#  mutate(year=year(date)) %>% 
#   filter(year==2019) %>% pull(ba_m2) %>% unique
#   ggplot(data=.,aes(x,y,fill=ba_m2))+
#   geom_tile()+
#   coord_equal()+
#   scale_fill_viridis_c()
