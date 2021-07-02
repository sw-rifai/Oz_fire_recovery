library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(arrow)
library(lubridate)

# WRANGLING ORDER:
# 1) Climate data
# 2) LAI data
# 3) Landscape data

# Required variables for predictions ------------------------------------------
# min_nbr_anom+
#   fire_month +
#   malai+
#   map+mapet+mavpd15+
#   matmax+matmin+
#   elevation+sand+pH+silt+
#   der+slope+aspect+
#   vpd15_anom_3mo+
#   post_vpd15_anom_12mo+
#   pre_fire_slai_anom_3mo + 
#   pre_fire_slai_anom_12mo + 
#   post_precip_anom_12mo+
#   precip_anom_12mo,

bs_start_date <- ymd("2019-09-01")
bs_end_date <- ymd("2020-02-01")

# Load initial pre & post_fire climate -------------------------
post_clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                                 col_select = c("x","y","date",
                                                "map","mapet","mappet","mavpd15",
                                                "matmin","matmax",
                                                "precip",
                                                "precip_u",
                                                "precip_anom",
                                                "precip_anom_12mo",
                                                "post_precip_anom_12mo",
                                                "vpd15",
                                                "vpd15_u",
                                                "vpd15_anom",
                                                "vpd15_anom_12mo",
                                                "tmax_anom_12mo",
                                                "tmax",
                                                "tmax_u",
                                                "tmin",
                                                "tmin_u",
                                                "post_vpd15_anom_12mo", 
                                                "post_tmax_anom_12mo",
                                                "post_pet_anom_12mo")) %>% 
  as.data.table()
post_clim <- post_clim[,`:=`(x = round(x,2), 
                   y = round(y,2))]

post_clim <- post_clim[,month:=month(date)]
post_clim[date==ymd("2020-12-01")] %>% 
  .[,.(x,y,precip_anom_12mo)] %>% 
  st_as_stars(., c("x","y")) %>% 
  plot


post_clim2 <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/rescaled-era5land_202101_202103.parquet")
post_clim2 <- post_clim2[,`:=`(precip_anom_12mo = pred_precip_12mo-map, 
                 vpd15_anom_12mo = pred_vpd15_12mo - mavpd15)]

post_clim[is.na(vpd15)==F]$date %>% max
# END SECTION ******************************************************************


# Prepare VPD -------------------------------------------------------------
p_vpd15 <- post_clim[date<ymd("2021-01-01")][,.(x,y,date,vpd15,vpd15_u,vpd15_anom,vpd15_anom_12mo,mavpd15)]

p2_vpd15 <- post_clim2[date>max(p_vpd15[is.na(vpd15)==F]$date)][,.(x,y,date,pred_vpd15,pred_vpd15_12mo)] %>% 
  rename(vpd15 = pred_vpd15, 
         vpd15_12mo = pred_vpd15_12mo) %>% 
  as.data.table()
p2_vpd15 <- p2_vpd15[,month:=month(date)]

p2_vpd15 <- merge(p2_vpd15,
      unique(post_clim[,.(x,y,month,mavpd15,vpd15_u)]), 
      by=c("x","y","month"))
p2_vpd15 <- p2_vpd15[,vpd15_anom := vpd15 - vpd15_u]

tmp <- rbindlist(list(p_vpd15, p2_vpd15), 
                 fill=TRUE)
tmp <- tmp[,vpd15_anom := vpd15 - vpd15_u]
# unique(tmp$date) %>% sort %>% plot
tmp <- tmp[order(x,y,date)][
  , vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                 n = 3,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][
  , post_vpd15_anom_12mo := frollapply(vpd15_anom,FUN=mean,
                                 n = 12,fill = NA,align='left'), by=.(x,y)]
p_vpd15 <- tmp
rm(tmp)

p_vpd15 <- p_vpd15[,.(x,y,date,mavpd15, vpd15_anom_3mo,vpd15_anom_12mo,post_vpd15_anom_12mo)]
# END SECTION ******************************************************************


# Prepare Precip -------------------------------------------------------------
p_precip <- post_clim[date<ymd("2021-01-01")][
  ,.(x,y,date,precip,precip_u,precip_anom,precip_anom_12mo,map)]

p2_precip <- post_clim2[date>max(p_precip[is.na(precip)==F]$date)][
  ,.(x,y,date,pred_precip,pred_precip_12mo)] %>% 
  rename(precip = pred_precip, 
         precip_12mo = pred_precip_12mo) %>% 
  as.data.table()
p2_precip <- p2_precip[,month:=month(date)]

p2_precip <- merge(p2_precip,
                  unique(post_clim[,.(x,y,month,map,precip_u)]), 
                  by=c("x","y","month"))
p2_precip <- p2_precip[,precip_anom := precip - precip_u]

tmp <- rbindlist(list(p_precip, p2_precip), 
                 fill=TRUE)
tmp <- tmp[,precip_anom := precip - precip_u]
# unique(tmp$date) %>% sort %>% plot
# tmp <- tmp[order(x,y,date)][
#   , precip_anom_3mo := frollapply(precip_anom,FUN=mean,
#                                  n = 3,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][
  , precip_anom_12mo := frollapply(precip_anom,FUN=mean,
                                        n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][
  , post_precip_anom_12mo := frollapply(precip_anom,FUN=mean,
                                       n = 12,fill = NA,align='left'), by=.(x,y)]
p_precip <- tmp
rm(tmp)

p_precip <- p_precip[,.(x,y,date,map, 
                        precip_anom_12mo,post_precip_anom_12mo)]
# END SECTION ******************************************************************


# Merge clim vars ---------------------------------------------------------
#   map+mapet+mavpd15+
#   matmax+matmin+
#   vpd15_anom_3mo+
#   post_vpd15_anom_12mo+
#   post_precip_anom_12mo+
#   precip_anom_12mo,

#   bs_start_date <- ymd("2019-08-01")

post_clim <- merge(p_vpd15, p_precip, by=c("x","y","date"))

ma_clim <- unique(arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                                 col_select = c("x","y",
                                                "matmax","matmin","mapet","mappet",
                                 )))
ma_clim <- ma_clim[,`:=`(x=round(x,2),
              y=round(y,2))]
post_clim <- merge(post_clim, ma_clim,by=c("x","y"))
post_clim <- post_clim[date>=bs_start_date]
# END SECTION ******************************************************************


# Prepare pre/post LAI data ----------------------------------------------------
# min_nbr_anom
#   fire_month 
#   malai
#   pre_fire_slai_anom_3mo  
#   pre_fire_slai_anom_12mo 

dlai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
id_bs <- dlai %>% 
  lazy_dt() %>%
  filter(date >= bs_start_date) %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y,id) %>%
  summarize(nburns = n()) %>%
  as.data.table()

firedate_bs <- dlai[id%in%id_bs$id] %>% lazy_dt() %>%
  filter(date>=bs_start_date) %>% 
  filter(fire_doy > 0) %>%
  group_by(id) %>%
  summarize(date_fire1 = first(date)) %>%
  ungroup() %>%
  select(id,date_fire1) %>%
  as.data.table()

dlai2 <- merge(dlai[date>=bs_start_date][id%in%firedate_bs$id], 
      firedate_bs, by='id')

dlai2 <- dlai2 %>% lazy_dt() %>%
  mutate(days_since_fire = as.double(date - date_fire1)) %>%
  as.data.table()

d_min_nbr <- dlai2 %>% 
  lazy_dt() %>% 
  # mutate(nbr_anom = nbr - nbr_u) %>% 
  filter(days_since_fire <= 100) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

dlai3 <- merge(dlai2[date==date_fire1][,.(x,y,id,date_fire1,slai_anom_3mo,slai_anom_12mo,malai)], 
      d_min_nbr, 
      by='id')

dlai3 <- dlai3[,fire_month := lubridate::month(date_fire1)] %>% 
              .[fire_month %in% c(8,9,10,11,12,1,2,3)]
dlai3
# id_bs %>% ggplot(data=.,aes(x,y,fill=nburns))+geom_tile()

# END SECTION ******************************************************************



# PART : PREP landscape data ----------------------------------------------------------------
# Data load 
# landscape covars 
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
d_soil <- merge(d_soil, dlai3[,.(x,y,id)], by='id')

# species 
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-SEOZcoastal-grid.parquet")
dom <- dom[,.(x,y,id,dom_sp)]

# Do the small merges up-front to catch which ds are missing data
dat <- merge(d_soil, dom %>% select(-x,-y) %>% as.data.table(), by='id',allow.cartesian = TRUE)
dat <- merge(dat, dlai3, by=c("x","y","id"))

#  Attach clim pixel id to VI ------------------------------------
coords_vi <- lazy_dt(dat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_clim <- unique(post_clim[,.(x,y)])
coords_clim_sf <- st_as_sf(coords_clim, coords = c('x','y'))
st_crs(coords_clim_sf) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_clim_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
coords_clim <- coords_clim %>% mutate(idx_clim = row_number()) %>% as.data.table()
gc(full=TRUE)
coords_vi <- coords_vi %>% st_drop_geometry() %>% as.data.table()
coords_vi$idx_clim <- coords_clim[nn_coords$nn.idx,]$idx_clim
gc(full=TRUE)

# merges
gc(full=TRUE)
post_clim <- merge(post_clim,
                   coords_clim,by=c('x','y'))
gc(full=TRUE)
dat <- merge(dat, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
post_clim <- post_clim[idx_clim %in% unique(dat$idx_clim)]


# Create 'dat' for writing to disk --------------------------------------------------------
dat <- merge(dat, post_clim %>% select(-x,-y) %>% as.data.table, 
             by=c('idx_clim'),allow.cartesian = TRUE)
dat <- dat[date==date_fire1]
dat <- dat[,`:=`(post_precip_anom_frac = post_precip_anom_12mo/map,
                 precip_anom_frac = precip_anom_12mo/map,
                 post_vpd15_anom_frac = post_vpd15_anom_12mo/mavpd15)]
dat[,fire_month:=month(date_fire1)]
dat[,month := month(date_fire1)]


# Filter dat to just Eucs in the Sep-Feb burning -------------------
dat <- dat[date>=bs_start_date & date<=bs_end_date]
dat <- dat[vc %in% c(2,3,5)][month %in% c(9,10,11,12,1,2)]
nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
# sp_fac <-
#   unique(dat[dom_sp %in% nobs[rank <= 30]$dom_sp][,.(dom_sp,hab_hnd)]) %>%
#   .[order(hab_hnd)] %>%
#   .[,sp_fac := forcats::fct_inorder(dom_sp)]
# Final mutations
# dat[,fire_month_f := lubridate::month(date_fire1,label = TRUE,abbr = TRUE)][
#   ,fire_month_f := factor(fire_month_f,
#                           levels=c("Sep","Oct","Nov","Dec","Jan","Feb"),
#                           ordered = TRUE)][,dom_sp_f := factor(dom_sp)]
# dat[,vc_name_f := factor(vc_name)]


# Write Black Summer pred data to disk -----------------------------------------
arrow::write_parquet(dat, 
                     sink=paste0("../data_general/proc_data_Oz_fire_recovery/BlackSummer_pred-data_",Sys.Date(),".parquet"), 
                     compression = 'snappy')


