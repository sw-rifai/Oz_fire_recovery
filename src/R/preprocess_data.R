# Load packages in this order (important)
library(phenofit);
library(tidyverse);
library(usethis);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(RcppArmadillo)
source("src/R/functions_time_to_recover.R")

# Load data ---------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

# Import Fire ---------------------------------------------------
tmp1_fire <- read_stars("../data_general/MCD64/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
                        proxy=TRUE) %>% 
  set_names("fire_doy") %>% 
  st_warp(., slice(tmp1,'band',1), use_gdal = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
                    names='date')
gc()
tmp2_fire <- read_stars("../data_general/MCD64/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
                        proxy=TRUE) %>% 
  set_names("fire_doy") %>% 
  st_warp(., slice(tmp2,'band',1), use_gdal = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
                    names='date')
gc()
tmp3_fire <- read_stars("../data_general/MCD64/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
                        proxy=TRUE) %>% 
  set_names("fire_doy") %>% 
  st_warp(., slice(tmp3,'band',1), use_gdal = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
                    names='date')
gc()
tmp4_fire <- read_stars("../data_general/MCD64/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
                        proxy=TRUE) %>% 
  set_names("fire_doy") %>% 
  st_warp(., slice(tmp4,'band',1), use_gdal = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
                    names='date')
gc()
tmp_fire <- rbindlist(list(
  as.data.table(tmp1_fire)[fire_doy>0],
  as.data.table(tmp2_fire)[fire_doy>0],
  as.data.table(tmp3_fire)[fire_doy>0],
  as.data.table(tmp4_fire)[fire_doy>0]))
gc(full=TRUE)
# arrow::write_parquet(tmp_fire, sink="/home/sami/scratch/mcd64_se_coastal_fire.parquet")
# tmp_fire <- arrow::read_parquet(file = "/home/sami/scratch/mcd64_se_coastal_fire.parquet")


# # Import DEM ----------------------------------------------------
# tmp_dem <- stars::read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif",
#                              proxy = F) %>% 
#   st_warp(., tmp1[,,,1], use_gdal = F) %>% 
#   set_names("elevation") %>% 
#   as.data.table()
# gc()

cci1 <- stars::read_stars("../data_general/MCD43/MOD_CCI_500m_SE_coastal2001-01-01_to_2020-12-31-0000000000-0000000000.tif", 
                           proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("cci") %>% 
  as.data.table() %>% 
  .[is.na(cci)==FALSE]; gc(full=TRUE)
cci2 <- stars::read_stars("../data_general/MCD43/MOD_CCI_500m_SE_coastal2001-01-01_to_2020-12-31-0000000000-0000002304.tif", 
                          proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("cci") %>% 
  as.data.table() %>% 
  .[is.na(cci)==FALSE]; gc(full=TRUE)
cci3 <- stars::read_stars("../data_general/MCD43/MOD_CCI_500m_SE_coastal2001-01-01_to_2020-12-31-0000002304-0000000000.tif", 
                          proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("cci") %>% 
  as.data.table() %>% 
  .[is.na(cci)==FALSE]; gc(full=TRUE)
cci4 <- stars::read_stars("../data_general/MCD43/MOD_CCI_500m_SE_coastal2001-01-01_to_2020-12-31-0000002304-0000002304.tif", 
                          proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("cci") %>% 
  as.data.table() %>% 
  .[is.na(cci)==FALSE]; gc(full=TRUE)
cci <- rbindlist(list(cci1,cci2,cci3,cci4), use.names = TRUE)
arrow::write_parquet(cci, sink="/home/sami/scratch/mcd43_se_coastal_cci.parquet")
rm(cci,cci1,cci2,cci3,cci4)

# Import SWIR and Band ------------------------------------------
swir1 <- stars::read_stars("../data_general/MCD43/MCD43A4_swir_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000000000.tif", 
                           proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("swir") %>% 
  as.data.table() %>% 
  .[is.na(swir)==FALSE]; gc(full=TRUE)
swir2 <- stars::read_stars("../data_general/MCD43/MCD43A4_swir_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000002304.tif", 
                           proxy = F)%>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("swir") %>% 
  as.data.table() %>% 
  .[is.na(swir)==FALSE]; gc(full=TRUE)
swir3 <- stars::read_stars("../data_general/MCD43/MCD43A4_swir_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000002304-0000000000.tif", 
                           proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("swir") %>% 
  as.data.table() %>% 
  .[is.na(swir)==FALSE]; gc(full=TRUE)
swir4 <- stars::read_stars("../data_general/MCD43/MCD43A4_swir_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000002304-0000002304.tif", 
                           proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("swir") %>% 
  as.data.table() %>% 
  .[is.na(swir)==FALSE]; gc(full=TRUE)
swir <- rbindlist(list(swir1,swir2,swir3,swir4), use.names = TRUE)
arrow::write_parquet(swir,sink='/home/sami/scratch/mcd43_se_coastal_swir.parquet', 
                     compression = 'snappy')
rm(swir1,swir2,swir3,swir4); gc(full=TRUE)

# Import Band 6------------------------------------------
b6_1 <- stars::read_stars("../data_general/MCD43/MCD43A4_b6_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000000000.tif", 
                           proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("b6") %>% 
  as.data.table() %>% 
  .[is.na(b6)==FALSE]; gc(full=TRUE)
b6_2 <- stars::read_stars("../data_general/MCD43/MCD43A4_b6_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000002304.tif", 
                           proxy = F)%>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("b6") %>% 
  as.data.table() %>% 
  .[is.na(b6)==FALSE]; gc(full=TRUE)
b6_3 <- stars::read_stars("../data_general/MCD43/MCD43A4_b6_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000002304-0000000000.tif", 
                           proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("b6") %>% 
  as.data.table() %>% 
  .[is.na(b6)==FALSE]; gc(full=TRUE)
b6_4 <- stars::read_stars("../data_general/MCD43/MCD43A4_b6_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000002304-0000002304.tif", 
                           proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("b6") %>% 
  as.data.table() %>% 
  .[is.na(b6)==FALSE]; gc(full=TRUE)
b6 <- rbindlist(list(b6_1,b6_2,b6_3,b6_4), use.names = TRUE)
arrow::write_parquet(b6,sink='/home/sami/scratch/mcd43_se_coastal_b6.parquet', 
                     compression = 'snappy')
rm(b6_1,b6_2,b6_3,b6_4); gc(full=TRUE)


# The order of subsequent data imports is to alleviate memory usage
# Import part1 NIR & Red ----------------------------------------------
tmp1 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000000000.tif", 
                          proxy = F) 
tmp2 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000001536.tif", 
                          proxy = F) 
tmp3 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000001536-0000000000.tif", 
                          proxy = F) 
tmp4 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000001536-0000001536.tif", 
                          proxy = F) 
gc()


# proc nir ------------------------------------------------------
tmp1_nir <- tmp1 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp1)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")
tmp2_nir <- tmp2 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp2)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")
tmp3_nir <- tmp3 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp3)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")
tmp4_nir <- tmp4 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp4)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")



# proc red ------------------------------------------------------
tmp1_red <- tmp1 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp1)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")
tmp2_red <- tmp2 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp2)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")
tmp3_red <- tmp3 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp3)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")
tmp4_red <- tmp4 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp4)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")


# cleanup
rm(tmp1,tmp2,tmp3,tmp4)
gc()







# Merge and cast to data.table --------------------------------------
dat1 <- c(tmp1_nir,tmp1_red) %>% as.data.table(); 
rm(tmp1_nir, tmp1_red); gc(full = T)
dat1 <- dat1[is.na(nir)==F]; gc(full = T)

dat2 <- c(tmp2_nir,tmp2_red) %>% as.data.table(); gc(full = T)
rm(tmp2_nir, tmp2_red); gc(full = T)
dat2 <- dat2[is.na(nir)==F]; gc(full = T)

dat3 <- c(tmp3_nir,tmp3_red) %>% as.data.table(); gc(full = T)
rm(tmp3_nir, tmp3_red); gc(full = T)
dat3 <- dat3[is.na(nir)==F]; gc(full = T)

dat4 <- c(tmp4_nir,tmp4_red) %>% as.data.table(); gc(full = T)
rm(tmp4_nir, tmp4_red); gc(full = T)
dat4 <- dat4[is.na(nir)==F]; gc(full = T)

gc(full = T,reset = T)
dat <- data.table::rbindlist(list(dat1,dat2,dat3,dat4)); gc(full = T)
rm(dat1,dat2,dat3,dat4)
gc(full=T,reset = T)

dat <- merge(dat,tmp_fire,by=c("x","y","date"),all=TRUE,allow.cartesian = TRUE)
gc(full=TRUE)
dat <- dat[is.na(nir)==FALSE]
gc()
# dat <- merge(dat,as.data.table(tmp_dem),by=c("x","y"),allow.cartesian = TRUE,all.x = TRUE)
gc(full=TRUE)

# Add id, year, month -------------------------------------------
dat[, `:=`(id = .GRP), keyby = .(x, y)]
dat[, `:=`(year = year(date), month = month(date))]
gc()

# dat[id==95854] %>% 
#   # mutate(ndvi = (nir-red)/(nir+red) ) %>% 
#   ggplot(data=.,aes(date,ndvi))+
#   geom_line()+
#   geom_line(aes(date,sndvi),color='blue')

# Fn: Smooth data with Whittaker filter -----------------------------
smooth_ndvi <- function(din){
  din <- din[order(date)]
  ndvi <- (din$nir-din$red)/(din$nir+din$red)
  x1 <- data.table::nafill(ndvi,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  out <- din
  out$sndvi <- x3
  out$ndvi <- ndvi
  return(out)
}

# Apply smoothing function -------------------------------------------
system.time(dat <- dat[,smooth_ndvi(.SD), by=.(x,y,id)])
gc(full=TRUE)

# Place saver ***
# arrow::write_parquet(dat, sink="/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")
# dat <- arrow::read_parquet(file ="/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")

swir <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_swir.parquet")
gc(full=TRUE)
dat <- merge(dat, swir, by=c("x","y","date"),all.x = TRUE)
dat <- dat %>% lazy_dt() %>% 
  mutate(nbr = (nir-swir)/(swir+nir)) %>% 
  as.data.table()

# Calculate anomalies -----------------------------------------------------
gc(full=TRUE)
dat_norms <- dat[, `:=`(month = month(date))] %>%
  .[, .(ndvi_u = mean(sndvi, na.rm=TRUE),
        ndvi_mmax = max(sndvi,na.rm=TRUE),
        ndvi_sd = sd(sndvi, na.rm=TRUE)),
    keyby = .(x,y,month)]
dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>%
  mutate(ndvi_anom = sndvi-ndvi_u) %>%
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_sd,
         ndvi_fanom = sndvi/ndvi_u) %>%
  mutate(ndvi_fmax = sndvi/ndvi_mmax) %>% 
  as.data.table()
# dat[,`:=`(ndvi_mmax = max(sndvi,na.rm=TRUE)), keyby=.(x,y,month)]
# dat[,`:=`(ndvi_fmax = sndvi/ndvi_mmax), keyby=.(x,y,month)]
gc(full=TRUE)


# Remove bad grid cells ----------------------------------
# some locs have consistently negative ndvi; salt beds?
bad_pix <- dat %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(val = median(sndvi,na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(val <= 0.15) %>% 
  as.data.table()

# dat %>% lazy_dt() %>% 
#   group_by(x,y) %>% 
#   summarize(val = median(sndvi,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(x,y,fill=val))+
#   geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
#   geom_tile()+
#   coord_sf(xlim = c(140,154),
#            ylim = c(-40,-25), expand = FALSE)+
#   scale_fill_gradient2()

dat <- dat[!id %in% bad_pix$id]; gc(full=TRUE)
dat <- arrow::read_parquet(file ="/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")
swir <- arrow::read_parquet(file ="/home/sami/scratch/mcd43_se_coastal_swir.parquet")
cci <- arrow::read_parquet(file ="/home/sami/scratch/mcd43_se_coastal_cci.parquet")
swir <- merge(swir,cci,by=c("x","y","date"),all=TRUE)
rm(cci); gc(full=TRUE)
dat <- merge(dat,swir,by=c("x","y","date"),all=TRUE)
rm(swir); gc(full=TRUE)
dat <- dat %>% mutate(date=as.Date(date))
arrow::write_parquet(dat,
                     sink="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet", 
                     compression = 'snappy')



# STAGE 5: NBR & CCI Anoms -----------------------------------------------------------
# calc normalized burn ratio
dat <- dat %>% lazy_dt() %>% 
  mutate(nbr = (nir-swir)/(nir+swir)) %>%
  as.data.table()


# Fn: Smooth CCI with Whittaker filter 
smooth_cci <- function(din){
  din <- din[order(date)]
  cci <- din$cci
  x1 <- data.table::nafill(cci,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  out <- din
  out$scci <- x3
  # out$cci <- cci
  return(out)
}
vec_no_cci <- dat %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(all_na = all(is.na(cci)==TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
vec_no_cci <- vec_no_cci[is.na(id)==FALSE][all_na==FALSE]$id
# smooth_cci(dat[id==1000]) %>% 
#   ggplot(data=.,aes(date,cci))+
#   geom_line()+
#   geom_line(aes(date,scci),col='red')
system.time(dat <- dat[id%in%vec_no_cci][,smooth_cci(.SD), by=.(x,y,id)])


# Calculate anomalies 
gc(full=TRUE)
dat_norms <- dat[, `:=`(month = month(date))] %>%
  .[, .(ndvi_u = mean(sndvi, na.rm=TRUE),
        nbr_u = mean(nbr,na.rm=TRUE), 
        cci_u = mean(scci,na.rm=TRUE),
        ndvi_sd = sd(sndvi, na.rm=TRUE), 
        nbr_sd = sd(nbr,na.rm=TRUE), 
        cci_sd = sd(cci,na.rm=TRUE)),
    keyby = .(x,y,month)]
dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>%
  mutate(ndvi_anom = sndvi-ndvi_u, 
         nbr_anom = nbr - nbr_u, 
         cci_anom = scci - cci_u) %>%
  as.data.table()
# dat[,`:=`(ndvi_mmax = max(sndvi,na.rm=TRUE)), keyby=.(x,y,month)]
# dat[,`:=`(ndvi_fmax = sndvi/ndvi_mmax), keyby=.(x,y,month)]
gc(full=TRUE)
dat <- dat[is.na(ndvi_u)==FALSE]

arrow::write_parquet(dat,
                     sink="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet", 
                     compression = 'snappy')
