pacman::p_load(tidyverse, data.table, lubridate, stars)

# Eucalypts (including the genera Angophora, Coimbra, and Eucalyptus) comprise X% of southeast Australia's woody ecosystems
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")

dat[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=malai))+
  geom_tile()+
  coord_equal()

bb <- st_as_stars(dat[date==min(date)][,.(x,y)], coords=c('x','y'), crs=st_crs(4326))
st_crs(bb) <- st_crs(4326)

nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif", 
                          proxy=F)
names(nvis) <- 'vc_code'
nvis_codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("vc_code","veg_class_descrip")) %>% 
  mutate(vc_name = as.factor(veg_class_descrip)) %>% 
  select(-veg_class_descrip)
nvis <- st_crop(nvis, bb)
plot(nvis)
lc <- merge(as.data.table(nvis), nvis_codes, by='vc_code')

lc[vc_code %in% c(2,3,4,5,11)] # Euc woody
lc[vc_code %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,23,30,31,32)] # all woody

lc[vc_code %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,23,30,31,32)] %>% 
  ggplot(data=., aes(x,y,fill=vc_name))+
  geom_tile()+
  geom_tile(data=lc[vc_code %in% c(2,3,4,5,11)], 
    aes(x,y),fill='red')+
  coord_equal()

val1 <- 100*nrow(lc[vc_code %in% c(2,3,4,5,11)])/nrow(lc[vc_code %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,23,30,31,32)])
paste0("Eucalypts (including the genera Angophora, Coimbra, and Eucalyptus) comprise ", 
format(val1,digits=2),"% of southeast Australia's woody ecosystems")

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

# Coastal and montane forests are typically wetter with a mean annual precipitation spanning X-Y mm yr-1, while the more inland Eucalypt woodlands only receive X-Y mm yr-1.
pacman::p_load(tidyverse, stars, data.table, lubridate, arrow, mgcv, broom, patchwork)

# Load vulnerable species 
vs <- fread("../data_general/ALA/Gallagher_cleaned/threatened_status.csv")
fn_simp <- function(x){
  x <- str_replace(x, " x ", " ")
  x <- str_remove(x, " sp\\.")
  x <- str_remove(x, " subsp\\.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}
vs <- vs[,species := fn_simp(Taxon),by=seq_len(nrow(vs))]
# tmp_sp_obs <- ala_mq[,species:=str_replace(species," ",".")]


# load fits ----
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-16 09:31:56.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>=0][r<0.024][r2>0.333][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
 ]
# estimate TTR from the logistic function
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
fits[,fire_year := year(date_fire1-months(3))]

# predicted dominant species ---- 
# sout <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-06-25 09:34:41.tiff")
# sout_rat <- fread("../data_general/proc_data_Oz_fire_recovery/predicted_species-distribution-ala-mq_top-40-species_LUT.csv")
# sout_rat

out <- read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-14 16:45:26.parquet")
sout <- set_names(st_as_stars(out[,.(x,y,predict)], dims = c("x","y"),crs=4326),c("species"))
st_crs(sout) <- st_crs(4326)
vv <- st_extract(sout,
  st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))

# Attach climate
clim_ma <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y",
                                           "map","mapet","mappet",
                                           "matmax","matmin","mavpd15")) %>% 
  distinct() %>% 
  as.data.table()
clim_ma <- st_as_stars(clim_ma, dims = c('x','y'))
st_crs(clim_ma) <- st_crs(4326)
vv <- st_extract(clim_ma,
                 st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))

vv <- st_extract(nvis,
                 st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))

lc[vc_code %in% c(2,3,4,5,11)] # Euc woody
lc[vc_code %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,23,30,31,32)] # all woody

lc[vc_code %in% c(2,3,4,5,11)] %>% 
  ggplot(data=., aes(x,y,fill=vc_name))+
  geom_tile()+
  coord_equal()

fits[vc_code %in% c(5)]$map %>% quantile(., c(0.01,0.99))
fits[vc_code %in% c(3,2)]$map %>% quantile(., c(0.01,0.99))


fits[vc_code %in% c(5)] %>% 
  ggplot(data=., aes(x,y,color=map))+
  geom_point()+
  coord_equal()+
  scale_color_viridis_b(option='H', n.breaks=15,direction = -1)+
  facet_wrap(~vc_code)

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

# Only a small proportion of pixel locations in the data dataset burned twice (X%) or even three times (X%) prior to the 2019/20 megafires.

# Isolate slow recovering pixels --------------------------------
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
               col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
nfires <- dat[date<=ymd("2019-08-01")][,.(nobs=.N, fires=sum(is.na(fire_doy)==F)),by=.(id,x,y)]

nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif", 
                          proxy=F)
names(nvis) <- 'vc_code'
nvis_codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("vc_code","veg_class_descrip")) %>% 
  mutate(vc_name = as.factor(veg_class_descrip)) %>% 
  select(-veg_class_descrip)

vv <- st_extract(nvis,
                 st_as_sf(nfires[,.(x,y)],coords=c("x","y"),crs=4326))
nfires <- bind_cols(nfires, st_drop_geometry(vv))
nfires[vc_code %in% c(2,3,4,5,11)]

lai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
lai_coords <- unique(lai[date==ymd('2005-01-01')][,.(x,y,lai)])

nfires <- merge(lai_coords,nfires,by=c("x","y"))
nfires <- nfires[vc_code %in% c(2,3,4,5,11)]

100*nrow(nfires[fires==2])/nrow(nfires[fires>=1])
100*nrow(nfires[fires>=3])/nrow(nfires[fires>=1])


#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
# More than X% of the study region had burned at least once between the 2001-2019 fire years
100*nrow(nfires[fires>0])/nrow(nfires)

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
top_mod <- h2o.loadModel("../data_general/proc_data_Oz_fire_recovery/xgboost_ttr5-lai-ocat_2021-06-27_/GBM_5_AutoML_20210627_093118")
top_mod


#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
# Annual burn area was notably correlated with both annual precipitation anomalies (rho = ;Fig. 2a) and/or low seasonal VPD (rho = ;Fig. 2b). High VPD also correlated with greater LAI loss from fire (rho = ; Fig. 2b).

tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
               col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
nfires <- dat[date<=ymd("2021-01-01")][,.(nobs=.N, fires=sum(is.na(fire_doy)==F)),by=.(id,x,y)]

nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif", 
                          proxy=F)
names(nvis) <- 'vc_code'
nvis_codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("vc_code","veg_class_descrip")) %>% 
  mutate(vc_name = as.factor(veg_class_descrip)) %>% 
  select(-veg_class_descrip)

vv <- st_extract(nvis,
                 st_as_sf(nfires[,.(x,y)],coords=c("x","y"),crs=4326))
nfires <- bind_cols(nfires, st_drop_geometry(vv))
nfires <- nfires[vc_code %in% c(2,3,4,5,11)]


clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]

# Attach AWAP pixel id to VI ***
coords_vi <- unique(nfires[,.(x,y,id)])
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
coords_awap <- coords_awap %>% mutate(idx_awap = row_number()) %>% as.data.table()
gc(full=TRUE)
coords_vi <- coords_vi %>% st_drop_geometry() %>% as.data.table()
coords_vi$idx_awap <- coords_awap[nn_coords$nn.idx,]$idx_awap
gc(full=TRUE)

# # merges
gc(full=TRUE)
clim <- merge(clim,coords_awap,by=c('x','y'))
gc(full=TRUE)
nfires <- merge(nfires, coords_vi, by='id')
gc(full=TRUE)


clim
clim <- clim[idx_awap%in%unique(nfires$idx_awap)][,fire_year:=lubridate::year(date-months(3))]
s_ba <- dat[date<=ymd("2021-01-01")][id%in%nfires$id][,fire_year:=lubridate::year(date-months(3))][,.(ba_ha = sum(is.na(fire_doy)==F)*25),by=.(fire_year)]
s_precip <- clim[,.(p_anom = mean(precip_anom_12mo)), by=.(fire_year)]
s_vpd15 <- clim[month%in%c(9,10,11,12,1,2)][,.(vpd15_anom = mean(vpd15_anom_3mo)), by=.(fire_year)]

merge(s_ba,s_precip)[fire_year>=2001&fire_year<=2019][,.(ba_ha,p_anom)] %>% cor
merge(s_ba,s_vpd15)[fire_year>=2001&fire_year<=2019][,.(ba_ha,vpd15_anom)] %>% cor

# Get cor(LAI_anom, VPD) ***********************************************
lai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
  col_select = c("x","y","id","date","slai_anom","month"))
locs <- unique(lai[,.(x,y,id)])

tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
               col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
ba <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
ba <- ba[is.na(fire_doy)==F]
ba[,`:=`(fire_month = month(date), 
         fire_year = year(date-months(3)))]
ba[,`:=`(date_fire1 = date)]

ssdat <- lai[id %in% nfires[fires==1]$id][id %in% ba$id]
ssdat <- merge(ssdat, 
  ba[,.(x,y,id,date_fire1)], 
  by=c("x","y","id"), 
  allow.cartesian = T)
ssdat <- ssdat[date >= date_fire1]
ssdat <- ssdat[date<= (date_fire1+months(6))]
ssdat[,`:=`(fire_year = year(date_fire1-months(3)), 
         fire_month = month(date_fire1))]
ssdat <- ssdat[fire_month %in% c(9,10,11,12,1,2)]
la2 <- ssdat[,.(min_slai_anom = min(slai_anom,na.rm=T)), 
  by=.(id,fire_year)]

s_lai_anom <- la2[,.(lai_anom = mean(min_slai_anom,na.rm=T)),by=.(fire_year)]
merge(s_lai_anom,s_vpd15)[fire_year>=2001&fire_year<=2019][,.(lai_anom,vpd15_anom)] %>% cor
#*******************************************************************************
