library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(arrow)
library(mgcv)
library(mgcViz)


dttr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-05-24 16:20:29.parquet")

# landscape covars ------------
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")

# species ----------------------
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
dom <- dom[,.(x,y,id,dom_sp)]
hab <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_species-median-clim-topo.parquet")



# pre & post_fire climate -------------------------
post_clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                                 col_select = c("x","y","date",
                                                "map","mapet","mappet","mavpd15",
                                                "matmin","matmax",
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
post_clim <- post_clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                         n = 3,fill = NA,align='right'), by=.(x,y)]


#  Attach clim pixel id to VI ------------------------------------
coords_vi <- lazy_dt(dttr) %>% select(x,y,id) %>% distinct() %>% as.data.table()
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
dat <- merge(dttr, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
post_clim <- post_clim[idx_clim %in% unique(dat$idx_clim)]


# Create 'dat'
dat <- merge(dat, post_clim %>% select(-x,-y) %>% as.data.table, 
             by=c('idx_clim','date'),allow.cartesian = TRUE)

dat <- dat[,`:=`(post_precip_anom_frac = post_precip_anom_12mo/map,
                 precip_anom_frac = precip_anom_12mo/map,
                 post_vpd15_anom_frac = post_vpd15_anom_12mo/mavpd15)]
dat[,fire_month:=month(date_fire1)]

dat <- merge(dat, d_soil, by='id')
dat <- merge(dat, dom, by='id')


# Filter dat to just Eucs in the Oct-Feb burning -------------------
dat <- dat[vc %in% c(2,3,5)][month %in% c(9,10,11,12,1,2)]
dom50obs <- dom[,.(nobs=.N), by=dom_sp][nobs>=50]
dat <- dat[is.na(ttr5_lai)==F][is.na(der)==F][elevation>0][is.na(dom_sp)==F][dom_sp%in%dom50obs$dom_sp]
dat <- dat[,dom_sp_f := factor(dom_sp)]
# nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
# sp_fac <-
#   unique(dat[dom_sp %in% nobs[rank <= 30]$dom_sp][,.(dom_sp,hab_hnd)]) %>%
#   .[order(hab_hnd)] %>%
#   .[,sp_fac := forcats::fct_inorder(dom_sp)]
# Final mutations
dat[,fire_month_f := lubridate::month(date_fire1,label = TRUE,abbr = TRUE)][
  ,fire_month_f := factor(fire_month_f,
                          levels=c("Sep","Oct","Nov","Dec","Jan","Feb"),
                          ordered = TRUE)]
dat[,vc_name_f := factor(vc_name)]
# dat[,r365 := r*365]



t1 <- bam(ttr5_lai ~ 
            s(x.x,y.x)+
            s(dom_sp_f,bs='re')+
            s(min_nbr_anom,k=3,bs='cs')+
            fire_month_f + 
            s(malai, k=5,bs='cs')+
            # s(map, k=5,bs='cs')+
            # s(mapet, k=5,bs='cs')+
            s(mavpd15,k=5,bs='cs')+
            s(I(matmax-matmin),k=5,bs='cs')+
            log(des)+
            s(pH,k=5,bs='cs')+
            vpd15_anom_3mo+
            precip_anom_frac + 
            post_precip_anom_frac + 
            vc_name_f, 
          family=Gamma(link='identity'),
          data=dat[sample(.N,75000)],
          discrete=TRUE, 
          select=TRUE)
summary(t1)
print(plot(getViz(t1),allTerms=T),pages=1)

dat[sample(.N,1000)] %>% 
  ggplot(data=.,aes(nbr_anom, slai_anom))+
  geom_point()+
  geom_smooth()



t2 <- bam(ttr5_lai ~ 
            te(log1p(elevation+10), 
               sand,
               silt)+
            s(min_nbr_anom,k=3,bs='cs')+
            fire_month_f + 
            s(malai, k=5,bs='cs')+
            # s(map, k=5,bs='cs')+
            # s(mapet, k=5,bs='cs')+
            s(mavpd15,k=5,bs='cs')+
            s(I(matmax-matmin),k=5,bs='cs')+
            # log(des)+
            s(pH,k=5,bs='cs')+
            vpd15_anom_3mo+
            # precip_anom_frac + 
            post_precip_anom_frac + 
            vc_name_f
          , 
          family=Gamma(link='identity'),
          data=dat[sample(.N,50000)],
          discrete=TRUE, 
          select=TRUE)
summary(t2)
rmse_vec(dat$ttr5_lai,
         predict(t2,newdata=dat,type='response'))
rsq_trad_vec(dat$ttr5_lai,
         predict(t2,newdata=dat,type='response'))


t3 <- bam(ttr5_lai ~ 
            s(dom_sp_f,bs='re')+
            s(min_nbr_anom,k=3,bs='cs')+
            fire_month_f + 
            s(malai, k=5,bs='cs')+
            s(mavpd15,k=5,bs='cs')+
            s(I(matmax-matmin),k=5,bs='cs')+
            log(des)+
            s(pH,k=5,bs='cs')+
            vpd15_anom_3mo+
            precip_anom_frac + 
            post_precip_anom_frac + 
            s(post_vpd15_anom_frac,k=5)+
            vc_name_f, 
          family=Gamma(link='identity'),
          data=dat[sample(.N,75000)],
          discrete=TRUE, 
          select=TRUE)
summary(t3)
print(plot(getViz(t3),allTerms=T),pages=1)


t4 <- bam(ttr5_lai ~ 
            dom_sp_f+
            te(mavpd15,
               vpd15_anom_3mo,
               post_vpd15_anom_12mo,
               k=5)+
            te(map, 
               post_precip_anom_12mo,
               precip_anom_12mo,
              k=5)+
            # te(map,
            #    precip_anom_12mo, 
            #    k=5)+
            s(min_nbr_anom,k=3,bs='cs')+
            fire_month_f + 
            # te(malai,
            #    lai_yr_sd,
            #    k=5)+
            s(matmin, 
              I(matmax-matmin),k=5)+
            # s(log(der))+
            # s(log(elevation))+
            # s(pH,k=5,bs='cs')+
            # s(vpd15_anom_3mo)+
            vc_name_f, 
          family=Gamma(link='identity'),
          data=dat[sample(.N,75000)],
          discrete=TRUE, 
          select=TRUE)
summary(t4)
print(plot(getViz(t4),allTerms=T),pages=1)
plot(t4,select = 1)
rmse_vec(dat$ttr5_lai,
         predict(t4,newdata=dat,type='response'))
rsq_trad_vec(dat$ttr5_lai,
             predict(t4,newdata=dat,type='response'))


dat[sample(.N, 1000)] %>% ggplot(data=.,aes(mavpd15,post_precip_anom_12mo))+geom_point()+
  geom_smooth

apply(dat, 2, function(x) sum(is.na(x)))

predict(t4,newdata=dat,type='response',na.action = na.pass) %>% 
  is.na %>% 
  table

f_dom_sp <- levels(dat$dom_sp_f)
vec_c <- coef(t4)
vec_c[str_detect(names(vec_c),pattern = "dom_sp_f")]


t5 <- bam(ttr5_lai ~ 
            dom_sp_f+
            te(mavpd15,vpd15_anom_3mo)+
            te(mavpd15, map, post_precip_anom_12mo, precip_anom_12mo),
          # family=Gamma(link='identity'),
          data=dat[sample(.N,75000)],
          discrete=TRUE, 
          select=TRUE)
summary(t5)

z1 <- coef(t4)[str_detect(names(coef(t4)),pattern = "dom_sp_f")]
z2 <- coef(t5)[str_detect(names(coef(t5)),pattern = "dom_sp_f")]
tibble(z1=z1[1:170],
       z2=z2[1:170]) %>% 
  ggplot(data=.,aes(z1,z2))+
  geom_point()+
  geom_smooth(method='lm')



tmp <- tibble(dom_sp_f = f_dom_sp)
tmp$xp <- c(0,coef(t5)[str_detect(names(coef(t5)),pattern = ":post")])
tmp <- tmp %>% mutate(species=as.character(dom_sp_f))
tmp <- left_join(tmp, hab, by='species')
tmp %>%
  filter(between(xp,-20,20)) %>% 
  ggplot(data=.,aes(map, xp))+
  geom_point()+
  geom_smooth(method='lm')
