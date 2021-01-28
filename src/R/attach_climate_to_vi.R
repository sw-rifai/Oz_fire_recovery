# library(raster); library(rasterVis)
library(arrow)
library(sf); library(stars)
library(tidyverse); 
library(data.table); library(lubridate);
library(dtplyr)
setDTthreads(threads=8)


#*******************************************************************************
# Get base Coords ---------------------------------------------------------
#*******************************************************************************
base <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet", 
                            col_select = c("x","y","id")) %>% 
  as.data.table() %>% lazy_dt() %>% distinct() %>% as.data.table()
#*******************************************************************************
#* END SECTION
#*******************************************************************************


#*******************************************************************************
# Extract AWAP clim grid cells for southeast coastal Oz ----------------------------------------------
#*******************************************************************************
atmax <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/tmax/AWAP_monthly_tmax_1970_2019.nc")
names(atmax) <- "tmax"
st_crs(atmax) <- st_crs(4326)
atmin <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/tmin/AWAP_monthly_tmin_1970_2019.nc")
names(atmin) <- "tmin"
st_crs(atmin) <- st_crs(4326)
avp9 <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/vph09/AWAP_monthly_vph09_1970_2019.nc")
avp15 <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/vph15/AWAP_monthly_vph15_1970_2019.nc")
names(avp9) <- "vp9"
names(avp15) <- "vp15"
st_crs(avp9) <- st_crs(4326)
st_crs(avp15) <- st_crs(4326)
aprecip <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/rain/AWAP_monthly_rain_1970_2019.nc")
names(aprecip) <- "precip"
st_crs(aprecip) <- st_crs(4326)
apet <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/pet/AWAP_monthly_PriestleyTaylor_PET_1990_2019.nc")
names(apet) <- "pet"
st_crs(apet) <- st_crs(4326)
gc(full=TRUE)


atmax <- filter(atmax,time>=ymd("2000-01-01"))
atmin <- filter(atmin,time>=ymd("2000-01-01"))
atmax <- c(atmax,atmin, try_hard=TRUE)

avp9 <- filter(avp9,time>=ymd("2000-01-01"))
avp15 <- filter(avp15,time>=ymd("2000-01-01"))
avp9 <- c(avp9,avp15)

aprecip <- filter(aprecip,time>=ymd("2000-01-01"))
apet <- filter(apet,time>=ymd("2000-01-01"))
apet <- apet %>% st_set_dimensions(.,3,values=st_get_dimension_values(aprecip,3))
attr(apet,"dimensions")$time <- attr(aprecip,"dimensions")$time # apet does not have interval dates for unknown reasons
aprecip <- c(aprecip,apet)
gc(full=TRUE)

clim <- c(atmax,avp9,aprecip)
rm(atmax,atmin,avp15,avp9,aprecip,apet); gc(full=TRUE)

eoz_box <- st_bbox(c(xmin = min(base$x),
                     ymin = min(base$y),
                     xmax = max(base$x),
                     ymax = max(base$y)), 
                   crs = st_crs(4326))
clim <- st_crop(clim, eoz_box)

clim <- clim %>% as.data.table() %>% 
                units::drop_units()
#*******************************************************************************
#* END SECTION
#*******************************************************************************

#*******************************************************************************
#* Calculate anomalies ---------------------------------------------------
#*******************************************************************************
clim <- clim %>% rename(x=longitude, y=latitude, date=time)
clim <- clim[is.na(tmax)==FALSE]
clim <- clim[is.na(pet)==FALSE]

clim <- clim[, `:=`(month = month(date))] # create month
clim <- clim[, `:=`(year = year(date))]   # create year
clim <- clim[,`:=`("ppet" = precip/pet)]
#' Calculates saturation vapour pressure
#' @return saturation vapour pressure
calc_esat <- function(airtemp){
  #Tair in degrees C
  
  #From Jones (1992), Plants and microclimate: A quantitative approach 
  #to environmental plant physiology, p110
  esat <- 613.75 * exp(17.502 * airtemp / (240.97+airtemp))
  
  return(esat)
}
clim <- clim %>% lazy_dt() %>% 
  mutate(vpd15 = 0.01*(calc_esat(tmax)/10 - vp15), 
         vpd9 = 0.01*(calc_esat(0.5*(tmax+tmin))/10 - vp9)) %>% 
  as.data.table

norms_clim <- lazy_dt(clim) %>% 
  group_by(x,y,month) %>% 
  summarize(tmax_u = mean(tmax,na.rm=TRUE), 
            tmin_u = mean(tmin,na.rm=TRUE), 
            vpd9_u = mean(vpd9,na.rm=TRUE), 
            vpd15_u = mean(vpd15,na.rm=TRUE), 
            precip_u = mean(precip,na.rm=TRUE), 
            pet_u = mean(pet,na.rm=TRUE), 
            ppet_u = mean(ppet,na.rm=TRUE),
            
            tmax_sd = sd(tmax,na.rm=TRUE), 
            tmin_sd = sd(tmin,na.rm=TRUE), 
            vpd9_sd = sd(vpd9,na.rm=TRUE), 
            vpd15_sd = sd(vpd15,na.rm=TRUE), 
            precip_sd = sd(precip,na.rm=TRUE), 
            pet_sd = sd(pet,na.rm=TRUE), 
            ppet_sd = sd(ppet,na.rm=TRUE)) %>% 
  as.data.table()

ma_clim <- lazy_dt(clim) %>% 
  group_by(x,y,year) %>% 
  summarize(tmax_ma = mean(tmax,na.rm=TRUE), 
            tmin_ma = mean(tmin,na.rm=TRUE), 
            vpd9_ma = mean(vpd9,na.rm=TRUE), 
            vpd15_ma = mean(vpd15,na.rm=TRUE), 
            precip_ma = mean(precip,na.rm=TRUE)*12, 
            pet_ma = mean(pet,na.rm=TRUE)*12, 
            ppet_ma = mean(ppet,na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(x,y) %>% 
  summarize(matmax = mean(tmax_ma,na.rm=TRUE), 
            matmin = mean(tmin_ma,na.rm=TRUE), 
            mavpd9 = mean(vpd9_ma,na.rm=TRUE), 
            mavpd15 = mean(vpd15_ma,na.rm=TRUE), 
            map = mean(precip_ma,na.rm=TRUE), 
            mapet = mean(pet_ma,na.rm=TRUE), 
            mappet = mean(ppet_ma,na.rm=TRUE)) %>% 
  as.data.table()

clim <- merge(clim, norms_clim, by=c("x","y","month"))
clim <- merge(clim, ma_clim, by=c("x","y"))

clim <- clim[, `:=`( 
                  precip_anom = precip-precip_u,  # calc raw anomaly 
                  pet_anom = pet-pet_u, 
                  ppet_anom = ppet-ppet_u, 
                  tmax_anom = tmax-tmax_u, 
                  tmin_anom = tmin-tmin_u, 
                  vpd9_anom = vpd9-vpd9_u,
                  vpd15_anom = vpd15 - vpd15_u)]
clim <- clim[, `:=`(
                  precip_anom_sd = precip_anom/precip_sd,  # calc sd anomaly 
                  pet_anom_sd = pet_anom/pet_sd, 
                  ppet_anom_sd = ppet_anom/ppet_sd, 
                  tmax_anom_sd = tmax_anom/tmax_sd, 
                  tmin_anom_sd = tmin_anom/tmin_sd, 
                  vpd15_anom_sd = vpd15_anom/vpd15_sd, 
                  vpd9_anom_sd = vpd9_anom/vpd9_sd)]
clim <- clim %>% lazy_dt() %>% mutate(date=as.Date(date)) %>% as.data.table()
#*******************************************************************************
#* END SECTION
#*******************************************************************************

#*******************************************************************************
# Calculate the multi-year anomalies ------
#*******************************************************************************
# calculate the rolling 12-month sums 
clim <- clim[order(x,y,date)][, precip_12mo := frollsum(precip,n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, pet_12mo := frollsum(pet,n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, ppet_12mo := frollmean(ppet,n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, tmax_anom_12mo := frollapply(tmax_anom,FUN=max,
                                                           n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, tmin_anom_12mo := frollapply(tmin_anom,FUN=max,
                                                           n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_12mo := frollapply(vpd15_anom,FUN=mean,
                                                       n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_anom_12mo := frollapply(vpd15_anom,FUN=mean,
                                                            n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_12mo := frollsum(precip_anom,n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, pet_anom_12mo := frollsum(pet_anom,n = 12,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, ppet_anom_12mo := frollmean(ppet_anom,n = 12,fill = NA,align='right'), by=.(x,y)]


clim <- clim[order(x,y,date)][, post_precip_12mo := frollsum(precip,n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_pet_12mo := frollsum(pet,n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_ppet_12mo := frollmean(ppet,n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_tmax_anom_12mo := frollapply(tmax_anom,FUN=max,
                                                             n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_tmin_anom_12mo := frollapply(tmin_anom,FUN=max,
                                                             n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_vpd15_12mo := frollapply(vpd15_anom,FUN=mean,
                                                         n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_vpd15_anom_12mo := frollapply(vpd15_anom,FUN=mean,
                                                              n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_precip_anom_12mo := frollsum(precip_anom,n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_pet_anom_12mo := frollsum(pet_anom,n = 12,fill = NA,align='left'), by=.(x,y)]
clim <- clim[order(x,y,date)][, post_ppet_anom_12mo := frollmean(ppet_anom,n = 12,fill = NA,align='left'), by=.(x,y)]

arrow::write_parquet(clim, sink="/home/sami/scratch/awap_clim_se_coastal.parquet", 
                     compression='snappy')




#*******************************************************************************
#* Low Memory Merge: Clim & VI
#*******************************************************************************
dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet")
gc(full=TRUE)

coords_vi <- lazy_dt(dat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)

coords_awap <- unique(clim[,.(x,y)])
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
st_crs(coords_awap_sf) <- st_crs(4326)


nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)

# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
coords_keep_awap$id_vi <- coords_vi$id
gc(full=TRUE)
coords_keep_awap %>% select(x,y,idx_awap) %>% distinct()

clim <- merge(clim,
              coords_keep_awap %>% select(x,y,idx_awap) %>% distinct(),by=c("x","y"), 
              all.y=TRUE)

dat <- merge(dat, coords_keep_awap %>% rename(id=id_vi) %>% select(id,idx_awap), by='id')
gc(full=TRUE)
dat <- merge(dat, clim[,.(idx_awap,date,precip_anom_12mo,pet_anom_12mo)], by=c("idx_awap","date"))
#*******************************************************************************
#* END SECTION
#*******************************************************************************










#' coords_awap <- unique(clim[,.(longitude,latitude)])
#' coords_awap <- coords_awap %>% lazy_dt() %>% rename(x=longitude, y=latitude) %>% as.data.table()
#' coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
#' 
#' coords_vi <- lazy_dt(base) %>% select(x,y) %>% distinct() %>% as.data.table()
#' coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
#' st_crs(coords_vi) <- st_crs(4326)
#' nn_coords <- RANN::nn2(
#'   coords_awap_sf %>% st_coordinates(),
#'   coords_vi %>% st_coordinates(), 
#'   k=1
#' )
#' # df of all coords in awap with idx
#' coords_keep_awap <- coords_awap %>% 
#'   mutate(idx_awap = row_number())
#' 
#' # subset df of awap coords to just those identified by nn_coords
#' coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
#' dim(coords_keep_awap)
#' 
#' coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
#'                       y_vi=st_coordinates(coords_vi)[,"Y"], 
#'                       x_clim=coords_keep_awap$x, 
#'                       y_clim=coords_keep_awap$y)
#' coords_dict <- setDT(coords_dict)
#' 
#' # test if awap coords object has equal number of rows as coords_vi
#' assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])
#' clim <- clim %>% rename(x=longitude,y=latitude)
#' # atmin <- atmin %>% as.data.table() %>% units::drop_units() %>%  rename(x=longitude,y=latitude)
#' # avp9 <- avp9 %>% as.data.table() %>% units::drop_units() %>%  rename(x=longitude,y=latitude)
#' # avp15 <- avp15 %>% as.data.table() %>% units::drop_units() %>%  rename(x=longitude,y=latitude)
#' # aprecip <- aprecip %>% as.data.table() %>% units::drop_units() %>%  rename(x=longitude,y=latitude)
#' # apet <- apet %>% as.data.table() %>% units::drop_units() %>%  rename(x=longitude,y=latitude)
#' gc(full=TRUE)
#' coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
#' clim <- merge(clim,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE);gc(full=TRUE)
#' clim <- clim %>% mutate(date=as.Date(time)) %>% select(-time) 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' clim <- rename(x_clim=x, y_clim=y) %>% 
#'   rename(x=x_vi, y=y_vi)
#' 
#' # atmin <- merge(atmin,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE);gc(full=TRUE)
#' # avp9 <- merge(avp9,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE);gc(full=TRUE)
#' # avp15 <- merge(avp15,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE);gc(full=TRUE)
#' # aprecip <- merge(aprecip,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE);gc(full=TRUE)
#' # apet <- merge(apet,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE);gc(full=TRUE)
#' 
#' clim <- clim[is.na(x_vi)==F]; gc(full=TRUE)
#' # atmin <- atmin[is.na(x_vi)==F]; gc(full=TRUE)
#' # avp9 <- avp9[is.na(x_vi)==F]; gc(full=TRUE)
#' # avp15 <- avp15[is.na(x_vi)==F]; gc(full=TRUE)
#' # aprecip <- aprecip[is.na(x_vi)==F]; gc(full=TRUE)
#' # apet <- apet[is.na(x_vi)==F]; gc(full=TRUE)
#' gc(full=TRUE)
#' dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet")
#' gc(full=TRUE)
#' 
#' dat1 <- merge(clim[date<=ymd("2010-01-01")],dat[date<=ymd("2010-01-01")], 
#'               by=c("x","y","date"), all.y=TRUE); gc(full=TRUE)
#' clim <- clim[date>ymd("2010-01-01")]; gc(full=TRUE)
#' dat <- dat[date>ymd("2010-01-01")]; gc(full=TRUE)
#' dat2 <- merge(clim,dat, 
#'               by=c("x","y","date"), all.y=TRUE); gc(full=TRUE)
#' gc(full=TRUE)
#' dat <- rbindlist(list(dat1,dat2),use.names = TRUE);
#' rm(dat1);rm(dat2); gc(full = TRUE)
#' dat <- arrow::write_parquet(dat,
#'                             sink="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_awapClim.parquet",
#'                             compression = 'snappy')
#' 
#' # vis check that vi and clim coords are close
#' coords_dict[sample(.N,1000)] %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
#' coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
#' coords_dict %>% head
#' coords_dict
#' 
#' 
#' # complicated way of doing full join
#' atmax %>% head
#' 
#' round(atmax$x[1:5])==round(atmax$x_vi[1:5])
#' round(atmax$y[1:5])==round(atmax$y_vi[1:5])
#' 
#' # visual check
#' atmax[time%in%c(ymd("2001-01-01",tz='UTC'),
#'                 ymd("2019-12-01",tz='UTC'))] %>%
#'   as_tibble() %>% 
#'   ggplot(data=., aes(x,y,fill=tmax))+
#'   geom_tile()+
#'   coord_equal()+
#'   scale_fill_viridis_c(direction = -1)+
#'   facet_grid(~as.factor(time))
#' gc(full=TRUE)
#' #*******************************************************************************
#' # END SECTION
#' #*******************************************************************************
#' 
#' #*******************************************************************************
#' # Extract AWAP tmin grid cells for east Oz -------------------------------------
#' #*******************************************************************************
#' atmin <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/tmin/AWAP_monthly_tmin_1970_2019.nc")
#' names(atmin) <- "tmin"
#' st_crs(atmin) <- st_crs(4326)
#' 
#' eoz_box <- st_bbox(c(xmin = min(base$x),
#'                      ymin = min(base$y),
#'                      xmax = max(base$x),
#'                      ymax = max(base$y)), 
#'                    crs = st_crs(4326))
#' atmin <- st_crop(atmin, eoz_box)
#' atmin <- atmin %>% as_tibble() %>% as.data.table()
#' atmin <- atmin %>% units::drop_units()
#' coords_awap <- atmin %>% select(longitude,latitude) %>% distinct()
#' coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
#' coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
#' 
#' coords_vi <- base %>% select(x,y) %>% distinct()
#' coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
#' st_crs(coords_vi) <- st_crs(4326)
#' nn_coords <- RANN::nn2(
#'   coords_awap_sf %>% st_coordinates(),
#'   coords_vi %>% st_coordinates(), 
#'   k=1
#' )
#' # df of all coords in awap with idx
#' coords_keep_awap <- coords_awap %>% 
#'   mutate(idx_awap = row_number())
#' 
#' # subset df of awap coords to just those identified by nn_coords
#' coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
#' dim(coords_keep_awap)
#' 
#' coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
#'                       y_vi=st_coordinates(coords_vi)[,"Y"], 
#'                       x_clim=coords_keep_awap$x, 
#'                       y_clim=coords_keep_awap$y)
#' coords_dict <- setDT(coords_dict)
#' 
#' # test if awap coords object has equal number of rows as coords_vi
#' assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])
#' 
#' 
#' 
#' # vis check that vi and clim coords are close
#' # coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
#' # coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
#' # coords_dict %>% head
#' coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
#' coords_dict
#' 
#' atmin <- atmin %>% rename(x=longitude,y=latitude)
#' 
#' # complicated way of doing full join
#' atmin <- merge(atmin,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
#' atmin <- atmin[is.na(x_vi)==F]
#' atmin %>% head
#' 
#' round(atmin$x[1:5])==round(atmin$x_vi[1:5])
#' round(atmin$y[1:5])==round(atmin$y_vi[1:5])
#' 
#' # visual check
#' atmin[time%in%c(ymd("1990-01-01",tz='UTC'),
#'                 ymd("2019-12-01",tz='UTC'))] %>%
#'   as_tibble() %>% 
#'   ggplot(data=., aes(x,y,fill=tmin))+
#'   geom_tile()+
#'   coord_equal()+
#'   scale_fill_viridis_c(direction = -1)+
#'   facet_grid(~as.factor(time))
#' gc()
#' #*******************************************************************************
#' # END SECTION
#' #*******************************************************************************
#' 
#' #*******************************************************************************
#' # Extract AWAP VP3pm (-> VPD) grid cells for east Oz ----------------------------------------------
#' #*******************************************************************************
#' #' Calculates saturation vapour pressure
#' #' @return saturation vapour pressure
#' calc_esat <- function(airtemp){
#'   #Tair in degrees C
#'   
#'   #From Jones (1992), Plants and microclimate: A quantitative approach 
#'   #to environmental plant physiology, p110
#'   esat <- 613.75 * exp(17.502 * airtemp / (240.97+airtemp))
#'   
#'   return(esat)
#' }
#' avp9 <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/vph09/AWAP_monthly_vph09_1970_2019.nc")
#' avp15 <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/vph15/AWAP_monthly_vph15_1970_2019.nc")
#' names(avp9) <- "vp9"
#' names(avp15) <- "vp15"
#' st_crs(avp9) <- st_crs(4326)
#' st_crs(avp15) <- st_crs(4326)
#' 
#' avp <- c(avp9,avp15)
#' 
#' 
#' eoz_box <- st_bbox(c(xmin = min(base$x),
#'                      ymin = min(base$y),
#'                      xmax = max(base$x),
#'                      ymax = max(base$y)), 
#'                    crs = st_crs(4326))
#' avp <- st_crop(avp, eoz_box)
#' avp <- avp %>% as_tibble() %>% as.data.table() %>% units::drop_units()
#' coords_awap <- avp %>% select(longitude,latitude) %>% distinct()
#' coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
#' coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
#' 
#' coords_vi <- unique(base[,.(x,y)])
#' coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
#' st_crs(coords_vi) <- st_crs(4326)
#' nn_coords <- RANN::nn2(
#'   coords_awap_sf %>% st_coordinates(),
#'   coords_vi %>% st_coordinates(), 
#'   k=1
#' )
#' # df of all coords in awap with idx
#' coords_keep_awap <- coords_awap %>% 
#'   mutate(idx_awap = row_number())
#' 
#' # subset df of awap coords to just those identified by nn_coords
#' coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
#' dim(coords_keep_awap)
#' 
#' coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
#'                       y_vi=st_coordinates(coords_vi)[,"Y"], 
#'                       x_clim=coords_keep_awap$x, 
#'                       y_clim=coords_keep_awap$y)
#' coords_dict <- setDT(coords_dict)
#' 
#' # test if awap coords object has equal number of rows as coords_vi
#' assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])
#' 
#' coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
#' 
#' avp <- avp %>% rename(x=longitude,y=latitude)
#' 
#' # complicated way of doing full join
#' avp <- merge(avp,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
#' avp <- avp[is.na(x_vi)==F]
#' 
#' 
#' avp <- merge(avp, atmax, by=c("x","y","x_vi","y_vi","time"))
#' 
#' avp <- avp %>% lazy_dt() %>% 
#'   mutate(vpd15 = 0.01*(calc_esat(tmax)/10 - vp15)) %>% 
#'   as.data.table()
#' gc(full=TRUE)
#' #*******************************************************************************
#' #* END SECTION
#' #*******************************************************************************
#' 
#' 
#' #*******************************************************************************
#' # Extract AWAP precip grid cells for east Oz ----------------------------------------------
#' #*******************************************************************************
#' aprecip <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/rain/AWAP_monthly_rain_1970_2019.nc")
#' names(aprecip) <- "precip"
#' st_crs(aprecip) <- st_crs(4326)
#' 
#' eoz_box <- st_bbox(c(xmin = min(base$x),
#'                      ymin = min(base$y),
#'                      xmax = max(base$x),
#'                      ymax = max(base$y)), 
#'                    crs = st_crs(4326))
#' aprecip <- st_crop(aprecip, eoz_box)
#' aprecip <- aprecip %>% as_tibble() %>% as.data.table()
#' aprecip <- aprecip %>% units::drop_units()
#' coords_awap <- aprecip %>% select(longitude,latitude) %>% distinct()
#' coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
#' coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
#' 
#' coords_vi <- base %>% select(x,y) %>% distinct()
#' coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
#' st_crs(coords_vi) <- st_crs(4326)
#' nn_coords <- RANN::nn2(
#'   coords_awap_sf %>% st_coordinates(),
#'   coords_vi %>% st_coordinates(), 
#'   k=1
#' )
#' # df of all coords in awap with idx
#' coords_keep_awap <- coords_awap %>% 
#'   mutate(idx_awap = row_number())
#' 
#' # subset df of awap coords to just those identified by nn_coords
#' coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
#' dim(coords_keep_awap)
#' 
#' coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
#'                       y_vi=st_coordinates(coords_vi)[,"Y"], 
#'                       x_clim=coords_keep_awap$x, 
#'                       y_clim=coords_keep_awap$y)
#' coords_dict <- setDT(coords_dict)
#' 
#' # test if awap coords object has equal number of rows as coords_vi
#' assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])
#' 
#' 
#' 
#' # vis check that vi and clim coords are close
#' coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
#' coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
#' coords_dict %>% head
#' coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
#' coords_dict
#' 
#' aprecip <- aprecip %>% rename(x=longitude,y=latitude)
#' 
#' # complicated way of doing full join
#' aprecip <- merge(aprecip,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
#' aprecip <- aprecip[is.na(x_vi)==F]
#' aprecip %>% head
#' 
#' round(aprecip$x[1:5])==round(aprecip$x_vi[1:5])
#' round(aprecip$y[1:5])==round(aprecip$y_vi[1:5])
#' 
#' # visual check
#' aprecip[time%in%c(ymd("1990-01-01",tz='UTC'),
#'                   ymd("2019-12-01",tz='UTC'))] %>%
#'   as_tibble() %>% 
#'   ggplot(data=., aes(x,y,fill=precip))+
#'   geom_tile()+
#'   coord_equal()+
#'   scale_fill_viridis_c(direction = -1)+
#'   facet_grid(~as.factor(time))
#' gc()
#' #*******************************************************************************
#' # END SECTION
#' #*******************************************************************************
#' 
#' #*******************************************************************************
#' # Extract AWAP pet grid cells for east Oz ----------------------------------------------
#' #*******************************************************************************
#' apet <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/pet/AWAP_monthly_PriestleyTaylor_PET_1990_2019.nc")
#' names(apet) <- "pet"
#' st_crs(apet) <- st_crs(4326)
#' 
#' eoz_box <- st_bbox(c(xmin = min(base$x),
#'                      ymin = min(base$y),
#'                      xmax = max(base$x),
#'                      ymax = max(base$y)), 
#'                    crs = st_crs(4326))
#' apet <- st_crop(apet, eoz_box)
#' apet <- apet %>% as_tibble() %>% as.data.table()
#' apet <- apet %>% units::drop_units()
#' coords_awap <- apet %>% select(longitude,latitude) %>% distinct()
#' coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
#' coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
#' 
#' coords_vi <- base %>% select(x,y) %>% distinct()
#' coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
#' st_crs(coords_vi) <- st_crs(4326)
#' nn_coords <- RANN::nn2(
#'   coords_awap_sf %>% st_coordinates(),
#'   coords_vi %>% st_coordinates(), 
#'   k=1
#' )
#' # df of all coords in awap with idx
#' coords_keep_awap <- coords_awap %>% 
#'   mutate(idx_awap = row_number())
#' 
#' # subset df of awap coords to just those identified by nn_coords
#' coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
#' dim(coords_keep_awap)
#' 
#' coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
#'                       y_vi=st_coordinates(coords_vi)[,"Y"], 
#'                       x_clim=coords_keep_awap$x, 
#'                       y_clim=coords_keep_awap$y)
#' coords_dict <- setDT(coords_dict)
#' 
#' # test if awap coords object has equal number of rows as coords_vi
#' assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])
#' 
#' 
#' 
#' # vis check that vi and clim coords are close
#' coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
#' coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
#' coords_dict %>% head
#' coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
#' coords_dict
#' 
#' apet <- apet %>% rename(x=longitude,y=latitude)
#' 
#' # complicated way of doing full join
#' apet <- merge(apet,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
#' apet <- apet[is.na(x_vi)==F]
#' apet %>% head
#' 
#' round(apet$x[1:5])==round(apet$x_vi[1:5])
#' round(apet$y[1:5])==round(apet$y_vi[1:5])
#' 
#' # visual check
#' apet[time%in%c(ymd("1990-01-01",tz='UTC'),
#'                ymd("2019-12-01",tz='UTC'))] %>%
#'   as_tibble() %>% 
#'   ggplot(data=., aes(x,y,fill=pet))+
#'   geom_tile()+
#'   coord_equal()+
#'   scale_fill_viridis_c(direction = -1)+
#'   facet_grid(~as.factor(time))
#' 
#' #*******************************************************************************
#' # END SECTION
#' #*******************************************************************************
#' 
#' 
#' 
#' 
#' 
#' #*******************************************************************************
#' # Extract ERA5 pet grid cells for east Oz ----------------------------------------------
#' #*******************************************************************************
#' e5pet <- stars::read_ncdf("../data_general/clim_grid/era5-land/Oz/Oz/Oz_era5-land_pet_1981_2019.nc")
#' names(e5pet) <- "pet"
#' st_crs(e5pet) <- st_crs(4326)
#' 
#' eoz_box <- st_bbox(c(xmin = min(base$x),
#'                      ymin = min(base$y),
#'                      xmax = max(base$x),
#'                      ymax = max(base$y)), 
#'                    crs = st_crs(4326))
#' e5pet <- st_crop(e5pet, eoz_box)
#' e5pet <- e5pet %>% as_tibble() %>% as.data.table()
#' e5pet <- e5pet %>% units::drop_units()
#' coords_awap <- e5pet %>% select(longitude,latitude) %>% distinct()
#' coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
#' coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
#' 
#' coords_vi <- base %>% select(x,y) %>% distinct()
#' coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
#' st_crs(coords_vi) <- st_crs(4326)
#' nn_coords <- RANN::nn2(
#'   coords_awap_sf %>% st_coordinates(),
#'   coords_vi %>% st_coordinates(), 
#'   k=1
#' )
#' # df of all coords in awap with idx
#' coords_keep_awap <- coords_awap %>% 
#'   mutate(idx_awap = row_number())
#' 
#' # subset df of awap coords to just those identified by nn_coords
#' coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
#' dim(coords_keep_awap)
#' 
#' coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
#'                       y_vi=st_coordinates(coords_vi)[,"Y"], 
#'                       x_clim=coords_keep_awap$x, 
#'                       y_clim=coords_keep_awap$y)
#' coords_dict <- setDT(coords_dict)
#' 
#' # test if awap coords object has equal number of rows as coords_vi
#' assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])
#' 
#' 
#' 
#' # vis check that vi and clim coords are close
#' coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
#' coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
#' coords_dict %>% head
#' coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
#' coords_dict
#' 
#' e5pet <- e5pet %>% rename(x=longitude,y=latitude)
#' 
#' # complicated way of doing full join
#' e5pet <- merge(e5pet,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
#' e5pet <- e5pet[is.na(x_vi)==F]
#' e5pet %>% head
#' 
#' round(e5pet$x[1:5])==round(e5pet$x_vi[1:5])
#' round(e5pet$y[1:5])==round(e5pet$y_vi[1:5])
#' 
#' # visual check
#' e5pet[time%in%c(ymd("1990-01-01",tz='UTC'),
#'                 ymd("2019-12-01",tz='UTC'))] %>%
#'   as_tibble() %>% 
#'   ggplot(data=., aes(x,y,fill=pet))+
#'   geom_tile()+
#'   coord_equal()+
#'   scale_fill_viridis_c(direction = -1)+
#'   facet_grid(~as.factor(time))
#' 
#' 
#' # ATTEMPTING TO RECALIBRATE ERA5-LAND PET TO AWAP PET
#' tmp_apet <- apet %>% select(x_vi,y_vi,time,pet) %>% 
#'   mutate(time=as.Date(time))
#' tmp_e5pet <- e5pet %>% select(x_vi,y_vi,time,pet) %>% 
#'   mutate(time=as.Date(time))
#' tmp_pet <- tmp_apet[tmp_e5pet,on=.(x_vi,y_vi,time)]
#' 
#' system.time(
#'   fit_pet <- tmp_pet %>% 
#'     filter(is.na(pet)==F) %>% 
#'     filter(is.na(i.pet)==F) %>% 
#'     sample_frac(0.35) %>% 
#'     group_by(x_vi,y_vi) %>% 
#'     summarize(beta0 = coef(lm(pet~i.pet+I(i.pet**2)))[1], 
#'               beta1 = coef(lm(pet~i.pet+I(i.pet**2)))[2],
#'               beta2 = coef(lm(pet~i.pet+I(i.pet**2)))[3]) %>% 
#'     ungroup()
#' )
#' 
#' setDT(fit_pet)
#' tmp_pet <- fit_pet[tmp_pet,on=.(x_vi,y_vi)]
#' tmp_pet <- tmp_pet[, `:=`(pet_pred = beta0 + beta1*i.pet + beta2*i.pet**2)]
#' 
#' # Visual Checks ***************
#' # tmp_pet %>% 
#' #   filter(is.na(pet_pred)==F) %>% 
#' #   sample_frac(0.005) %>% 
#' #   ggplot(data=., aes(pet_pred,pet))+
#' #   ggpointdensity::geom_pointdensity(alpha=0.1)+
#' #   scale_color_viridis_c()+
#' #   geom_smooth(method='lm')+
#' #   geom_abline(aes(intercept=0,slope=1),color='red')
#' 
#' # tmp_pet %>% 
#' #   filter(time==ymd("1995-11-01")) %>% 
#' #   ggplot(data=., aes(x_vi,y_vi,fill=pet_pred-pet))+
#' #   geom_tile()+
#' #   coord_equal()+
#' #   scale_fill_gradient2()
#' gc()
#' #*******************************************************************************
#' # END SECTION
#' #*******************************************************************************
#' 
#' #*******************************************************************************
#' # Join climate data to maintain date structure ---------
#' #*******************************************************************************
#' tmp_pet[time>=ymd("1981-01-01") & time<=ymd("1989-12-31")] %>% 
#'   .[,.(x_vi,y_vi,time,pet_pred)] %>% 
#'   .[,.(x_vi, y_vi,time,pet=pet_pred)] %>% 
#'   .[is.na(pet)==F]
#' 
#' jpet <- rbindlist(list(
#'   tmp_pet[time>=ymd("1981-01-01") & time<=ymd("1989-12-31")] %>% 
#'     .[,.(x_vi,y_vi,time,pet_pred)] %>% 
#'     .[,.(x_vi, y_vi,time,pet=pet_pred)],
#'   apet[,.(x_vi,y_vi,time,pet)][,`:=`(time=as.Date(time))]
#' ))
#' 
#' jpet <- jpet[is.na(pet)==F]
#' gc(reset = T, full = T)
#' 
#' # atmax <- atmax[,.(x_vi,y_vi,time,tmax)][,`:=`(time=as.Date(time))]
#' avp <- avp[,.(x_vi,y_vi,time,tmax,vp9,vp15,vpd15)][,`:=`(time=as.Date(time))]
#' attmax <- attmax[,`:=`(time=as.Date(time))]
#' atmin <- atmin[,`:=`(time=as.Date(time))]
#' aprecip <- aprecip[,.(x_vi,y_vi,time,precip)][,`:=`(time=as.Date(time))]
#' tmp_clim <- merge(avp,jpet,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)
#' tmp_clim <- tmp_clim[is.na(pet)==F][order(time,x_vi,y_vi)]
#' tmp_clim <- merge(tmp_clim,aprecip,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)
#' tmp_clim <- merge(tmp_clim,attmax,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)
#' tmp_clim <- merge(tmp_clim,atmin,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)
#' 
#' gc(reset = T, full=T)
#' #*******************************************************************************
#' # END SECTION
#' #*******************************************************************************
#' 
#' 
#' #*******************************************************************************
#' # base and NDVI (forest & woodlands) ------------------------------------------------
#' #*******************************************************************************
#' base <- arrow::read_parquet("../data_general/MCD43/MCD64_AVHRR_NDVI_hybrid_2020-05-18.parquet") %>% 
#'   as.data.table()
#' # base <- base %>% lazy_dt() %>% 
#' #   mutate(x = round(x, 2) ,
#' #          y = round(y, 2)) %>% 
#' #   as.data.table()
#' # names(base) <- "nirv"
#' # vec_dates <- seq(ymd("1982-01-01"),ymd("2019-12-01"),by="1 month")
#' # base <- st_set_dimensions(base, 3, values=vec_dates, names='date')
#' # base <- stars::read_stars("../data_general/base/base51_majorVegClass_0p05.tif")
#' # base2 <- st_warp(src=base, dest=base[,,,1], use_gdal = T)
#' # names(base2) <- "veg_class"
#' # base <- base2 %>% as_tibble() %>% as.data.table()
#' # codes <- readr::read_fwf("../data_general/base/base51_majorVegClass_codes.txt", 
#' #                          fwf_widths(c(2,100)), skip = 1) %>% 
#' #   set_names(c("veg_class","veg_class_descrip")) %>% 
#' #   mutate(vc = as.factor(veg_class_descrip))
#' # base <- inner_join(base, codes, by='veg_class')
#' # base <- base %>% filter(veg_class <= 15) # !!! only forests and woodlands !!!
#' # base <- base %>% select(-veg_class_descrip)
#' # base <- base %>% as_tibble() %>% as.data.table()
#' # base_vc <- base <- base[base,on=.(x,y)]
#' # base <- base_vc; rm(base_vc)
#' 
#' #*******************************************************************************
#' # END SECTION
#' #*******************************************************************************
#' 
#' 
#' #*******************************************************************************
#' # JOIN ALL THE PIECES -----------------------------------------------------
#' #*******************************************************************************
#' tmp_clim <- tmp_clim %>% rename(date=time)
#' base <- base %>% rename(x_vi=x,y_vi=y)
#' 
#' tmp_clim <- merge(tmp_clim, 
#'                   base,
#'                   by=c("x_vi","y_vi","date"), 
#'                   all=TRUE,allow.cartesian=TRUE)
#' tmp_clim %>% head
#' tmp_clim <- tmp_clim[is.na(pet)==F][order(date,x_vi,y_vi)]
#' 
#' arrow::write_parquet(tmp_clim, 
#'                      sink=paste0("../data_general/Oz_misc_data/ARD_ndvi_aclim_",Sys.Date(),".parquet"),
#'                      compression='gzip',compression_level = 9)
#' #*******************************************************************************
#' #*
#' #*******************************************************************************
#' rm(fit_pet, tmp_pet, tmp_apet, apet);
#' rm(aprecip, atmax,e5pet,jpet);
#' rm(tmp_e5pet)
#' gc(reset = T, full = T)
