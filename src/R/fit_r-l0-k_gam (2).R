library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(fasterize); 
library(furrr)
library(arrow)
library(mgcv)

oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

# fits --- 
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-05-19 18:18:29.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>0][r<0.024][r2>0.5][month%in%c(9,10,11,12,1,2)][
  ,delta:=(L0/K)
][delta<=0.75]

# dominant species ---- 
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
dom <- dom[,.(x,y,id,dom_sp)]

# species bioclim ranges ---
hab <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_species-median-clim-topo.parquet")
names(hab) <- c("species",paste0("hab_",names(hab)[-1]))

# mean annual climate -------------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y","map","mapet","mappet",
                                           "matmax","matmin","mavpd15")) %>% 
  distinct() %>% 
  as.data.table()
rclim <- st_as_stars(clim)
st_crs(rclim) <- st_crs(4326)



# merge into dat  ---------------------
dat <- merge(fits,dom,by=c("id","x","y"))
dat <- merge(dat, setnames(hab, "species","dom_sp"), by=c("dom_sp"))

# attach climate to dat ---------------
coords <- st_as_sf(dat[,.(x,y)],coords=c("x","y"))
st_crs(coords) <- st_crs(4326)
tmp <- st_extract(rclim, coords) %>% st_drop_geometry()
dat <- bind_cols(dat,tmp) %>% as.data.table()
rm(tmp)

# attach landscape covars ------------
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
dat <- merge(dat,d_soil,by='id')

# attach height above nearest drainage ------
rhnd <- read_stars("../data_general/Oz_misc_data/mert_hnd_500m_SE_coastal.tif") %>% 
  set_names("hnd")
tmp <- st_extract(rhnd, coords) %>% st_drop_geometry()
dat <- bind_cols(dat,tmp) %>% as.data.table()
rm(tmp)


# mean annual climate -------------------------
post_clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                                 col_select = c("x","y","date",
                                                "precip_anom_12mo", "post_precip_anom_12mo", 
                                                "vpd15_anom","post_vpd15_anom_12mo", 
                                                "post_tmax_anom_12mo",
                                                "post_pet_anom_12mo")) %>% 
  as.data.table()
post_clim <- post_clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                             n = 3,fill = NA,align='right'), by=.(x,y)]

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
post_clim <- merge(post_clim,coords_clim,by=c('x','y'))
gc(full=TRUE)
dat <- merge(dat, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
post_clim <- post_clim[idx_clim %in% unique(dat$idx_clim)]
dat <- merge(dat, post_clim, by=c('idx_clim','date'),allow.cartesian = TRUE)

dat <- dat[,`:=`(post_precip_anom_frac = post_precip_anom_12mo/map,
                 precip_anom_frac = precip_anom_12mo/map,
          post_vpd15_anom_frac = post_vpd15_anom_12mo/mavpd15)]


# Filter dat to just Eucs in the Oct-Feb burning -------------------
dat <- dat[vc %in% c(2,3,5)][month %in% c(10,11,12,1,2)]
nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
sp_fac <- 
  unique(dat[dom_sp %in% nobs[rank <= 30]$dom_sp][,.(dom_sp,hab_hnd)]) %>% 
  .[order(hab_hnd)] %>% 
  .[,sp_fac := forcats::fct_inorder(dom_sp)]


b1 <- bam(r~s(delta,k=5,bs='cs')+
            s( I(lai_yr_sd/malai))+
            I(post_precip_anom_12mo/map)
            # s(precip_anom_12mo,k=5,bs='cs')
          , 
          data=dat, 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(b1)
plot(b1, scheme=2, pages=1,rug=TRUE, 
     all.terms=TRUE)

b1 <- bam(r~s(delta,k=5)+
            s(post_vpd15_anom_frac,k=5)+
            # s(I(post_pet_anom_12mo/mapet),k=5)+
            s(post_precip_anom_frac,k=5)+
            s(precip_anom_frac,k=5)+
          # s(precip_anom_12mo,k=5,bs='cs')
            s(log1p(elevation),k=5)+
            s(malai,k=5)+
            s(lai_yr_sd,k=5)
          , 
          data=dat[,post_precip_anom_frac := post_precip_anom_12mo/map], 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(b1)
plot(b1, scheme=2, pages=1,rug=T, 
     all.terms=TRUE)

expand_grid(delta = seq(0.05,0.75,length.out=50), 
            precip_anom_frac = 0,
            post_precip_anom_frac = c(-0.4,0,0.4), 
            post_vpd15_anom_frac = c(-0.1,0,0.1), 
            malai = c(1,3,6), 
            lai_yr_sd = 0.5, 
            elevation = 0) %>% 
  mutate(pred = predict(b1, newdata=., type='response')) %>% 
  ggplot(data=.,aes(delta, pred,color=malai,group=malai))+
  geom_line()+
  facet_grid(post_precip_anom_frac ~ post_vpd15_anom_frac, labeller = label_both)



r1 <- bam(r~s(delta,k=5)+
            s(post_vpd15_anom_frac,k=5)+
            # s(I(post_pet_anom_12mo/mapet),k=5)+
            s(post_precip_anom_frac,k=5)+
            s(precip_anom_frac,k=5)+
            # s(precip_anom_12mo,k=5,bs='cs')
            s(sqrt(elevation),k=5)+
            s(malai,k=5)+
            s(lai_yr_sd,k=5)
          , 
          data=dat[,post_precip_anom_frac := post_precip_anom_12mo/map], 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)

delta1 <- bam(delta~
                  te(vpd15_anom_3mo,precip_anom_12mo,malai,lai_yr_sd,k=5),
          data=dat[sample(.N, 10000)], 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(delta1)
plot(delta1, scheme=2, all.terms=TRUE, pages=1)


k1 <- bam(malai ~ 
             s(log1p(elevation),
              aspect, 
              bs=c("tp","cc"))+
            s(map,mapet,matmax,k=5)+
            # s(map,k=5)+
            # s(mavpd15,k=5)+
            # s(mapet,k=5)+
            # s(aspect,k=5,bs='cc')+
            # te(mavpd15,mappet,k=5)+
            # te(map,mapet,k=5)+
            # s(map,mavpd15,k=5)+
            # s(awc,k=5)+
            s(log1p(der),k=5),
            # s(pH,k=5),
            # te(map,mappet,mavpd15,k=5), 
            data=dat[elevation>=0][sample(.N,50000)], 
            # family=Gamma(link='log'),
            select=TRUE, 
            discrete=TRUE)
summary(k1)
plot(k1, scheme=2, all.terms=TRUE, pages=1)

dat[sample(.N,1000)] %>% ggplot(data=.,aes(malai,K))+geom_point()+geom_smooth(method='lm')

# Scratch ----------------------------------------------------------------------
r1 <- ranger::ranger(r~delta+
                     post_precip_anom_frac+ 
                     post_vpd15_anom_frac+
                     malai+
                     lai_yr_sd+
                     elevation+
                     precip_anom_12mo, 
                     data=dat[is.na(elevation)==F][is.na(aspect)==F],
                     importance='impurity_corrected')
vip::vip(r1)                     
r1


lime::plot_features(r1)
