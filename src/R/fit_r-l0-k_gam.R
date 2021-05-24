library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(furrr)
library(arrow)
library(mgcv)

oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

# fits --- 
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-05-19 18:18:29.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>0][r<0.024][r2>0.5][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
][ldk<=0.75]

dttr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")

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
# post_clim <- post_clim[order(x,y,date)][, precip_anom_3mo := frollapply(precip_anom,FUN=mean,
#                                                                        n = 3,fill = NA,align='right'), by=.(x,y)]
# # 24 month
# post_clim <- post_clim[order(x,y,date)][, precip_anom_24mo := frollsum(precip_anom,n = 24,fill = NA,align='right'), by=.(x,y)]
# # 36 month
# post_clim <- post_clim[order(x,y,date)][, precip_anom_36mo := frollsum(precip_anom,n = 36,fill = NA,align='right'), by=.(x,y)]

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
dat <- merge(dat, post_clim %>% select(-x,-y) %>% as.data.table, 
             by=c('idx_clim','date'),allow.cartesian = TRUE)

dat <- dat[,`:=`(post_precip_anom_frac = post_precip_anom_12mo/map,
                 precip_anom_frac = precip_anom_12mo/map,
          post_vpd15_anom_frac = post_vpd15_anom_12mo/mavpd15)]
dat[,fire_month:=month(date_fire1)]

# Filter dat to just Eucs in the Oct-Feb burning -------------------
dat <- dat[vc %in% c(2,3,5)][month %in% c(9,10,11,12,1,2)]
nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
sp_fac <- 
  unique(dat[dom_sp %in% nobs[rank <= 30]$dom_sp][,.(dom_sp,hab_hnd)]) %>% 
  .[order(hab_hnd)] %>% 
  .[,sp_fac := forcats::fct_inorder(dom_sp)]
# Final mutations
dat[,delta:=K-L0]
dat[,fire_month_f := lubridate::month(date_fire1,label = TRUE,abbr = TRUE)][
  ,fire_month_f := factor(fire_month_f,
                          levels=c("Sep","Oct","Nov","Dec","Jan","Feb"),
                          ordered = TRUE)][,dom_sp_f := factor(dom_sp)]
dat[,vc_name_f := factor(vc_name)]
dat[,r365 := r*365]


# Predict delta: How much leaf area was lost due to fire
# delta = K-L0
delta1 <- bam(delta~
                vc_name+
                fire_month_f+
                malai+
                lai_yr_sd+
                te(mavpd15,vpd15_anom_3mo,k=5)+
                te(map,precip_anom_frac,k=5),
              family=Gamma(link='log'),
              data=dat[sample(.N, 10000)], 
              select=TRUE, 
              discrete=TRUE)
summary(delta1)
plot(delta1, scheme=2, all.terms=F, pages=1)

# Predict l0: How much leaf area is left after the fire
l01 <- bam(L0 ~
                vc_name+
                fire_month_f+
                malai+
                lai_yr_sd+
                te(mavpd15,vpd15_anom_3mo,k=5)+
                te(map,precip_anom_frac,k=5),
              family=Gamma(link='log'),
              data=dat[sample(.N, 10000)], 
              select=TRUE, 
              discrete=TRUE)
summary(l01)
plot(l01, scheme=2, all.terms=F, pages=1)


# Predict recovery LAI (k: carrying capacity)
k1 <- bam(K ~ 
            lai_yr_sd+
            malai,
          data=dat[sample(.N,50000)], 
          select=TRUE, 
          discrete=TRUE)
summary(k1)
print(plot(getViz(k1),allTerms=TRUE),pages=1)
plot(k1, scheme=2, all.terms=TRUE, pages=1)

# Predict r: Akin to the growth rate in the linear segment of the logistic function
r1 <- bam(r ~
            fire_month_f+ 
            log(ldk)+
            lai_yr_sd+
            I(matmax-matmin) +
            malai +
            # I(post_tmax_anom_12mo/matmax)+
            # I(post_vpd15_anom_12mo/mavpd15) +
            # post_tmin_anom_12mo+
            post_vpd15_anom_frac + 
            post_precip_anom_frac
          , 
          data=dat[sample(.N, 50000)], 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(r1)
print(plot(getViz(r1),allTerms=TRUE),pages=1)
gratia::observed_fitted_plot(r1)+
  geom_abline(col='red')
# plot(r1, scheme=2, pages=1,rug=TRUE, 
#      all.terms=TRUE)
visreg::visreg(r1,xvar = 'mavpd15')

# Predict r: Akin to the growth rate in the linear segment of the logistic function
r2 <- bam(exp(r)~fire_month_f+ 
          te(ldk,malai,
             lai_yr_sd,
             post_precip_anom_frac,k=5), 
          data=dat[sample(.N, 50000)], 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(r2)
plot(r2, scheme=2, pages=1,rug=TRUE, 
     all.terms=TRUE)
visreg::visreg(r2,xvar = 'mavpd15')



# beta regression approach
ldk1 <- bam(ldk~
              vc_name_f+
              # s(tmax_u,tmin_u, k=5)+
              
              fire_month_f+
              s(malai,k=3,bs='cs')+
              s(lai_yr_sd, k=3,bs='cs')+
              s(precip_anom_frac,k=3,bs='cs')+
              s(vpd15_anom_3mo,k=3,bs='cs'),
            data=dat[sample(.N, 50000)], 
            family=betar(),
            select=TRUE, 
            discrete=TRUE)
summary(ldk1)
plot(ldk1, scheme=2,all.terms=T,pages=1)
visreg::visreg(ldk1,scale='response',xvar='vpd15_anom_3mo',by='precip_anom_frac')
dat[sample(.N,10000)] %>% 
  .[,pred := predict(ldk1,newdata=.,type='response')] %>% 
  ggplot(data=.,aes(pred,ldk))+geom_point()+geom_abline(col='red')
predict(ldk1,type='response') %>% hist

# Predictions -------------------------------------------------------------
pdat <- dat
pdat$pred_delta <- predict(delta1, newdata=pdat, type='response')
pdat$pred_k <- predict(k1, newdata=pdat, type='response')
pdat$pred_l0 <- predict(l01, newdata=pdat, type='response')
pdat$pred_ldk <- predict(ldk1, newdata=pdat, type='response')
pdat$pred_r <- log(predict(r1, 
                           newdata=pdat, 
                           type='response'))
tmp <- pdat %>%
  lazy_dt() %>% 
  select(-ldk) %>% 
  rename(ldk = pred_ldk) %>% 
  as.data.table()
pdat$pred2_r <- log(predict(r1, 
                           newdata=tmp, 
                           type='response'))



# Plot predictions --------------------------------------------------------

# SM Fig: GOF L0-K-r ------------------------------------------------------
fig_k <- pdat %>% 
  ggplot(data=.,aes(pred_k,K))+
  geom_density_2d_filled()+
  scale_fill_viridis_d(option='A')+
  geom_smooth(method='lm', 
              col='lightblue')+
  geom_abline(col='#CF0000',lty=1)+
  scale_x_continuous(expand=c(0,0), limits=c(0,7))+
  scale_y_continuous(expand=c(0,0), limits=c(0,7))+
  labs(x=expression(paste("Predicted K ",bgroup("(",frac(m**2,m**2),")"),"")), 
       y=expression(paste("K ",bgroup("(",frac(m**2,m**2),")"),"")))+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'none')

fig_r <- pdat %>% 
  ggplot(data=.,aes(pred_r,r))+
  geom_density_2d_filled()+
  scale_fill_viridis_d(option='A')+
  geom_smooth(method='lm', 
              col='lightblue')+
  geom_abline(col='#CF0000',lty=1)+
  scale_x_continuous(expand=c(0,0), limits=c(0,0.01))+
  scale_y_continuous(expand=c(0,0), limits=c(0,0.01))+
  labs(x=expression(paste("Predicted r ",bgroup("(",frac(m**2,day),")"),"")), 
       y=expression(paste("r ",bgroup("(",frac(m**2,day),")"),"")))+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'none')

fig_l0 <- pdat %>% 
  ggplot(data=.,aes(pred_l0,L0))+
  geom_density_2d_filled()+
  scale_fill_viridis_d(option='A')+
  geom_smooth(method='lm', 
              col='lightblue')+
  geom_abline(col='#CF0000',lty=1)+
  scale_x_continuous(expand=c(0,0), limits=c(0,1.5))+
  scale_y_continuous(expand=c(0,0), limits=c(0,1.5))+
  labs(x=expression(paste("Predicted ",L[0]," ", bgroup("(",m**2,")"),"")), 
       y=expression(paste(L[0],bgroup("(",m**2,")"),"")))+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'none')
fig_l0

library(patchwork)
ggsave(fig_l0+fig_k+fig_r+plot_layout(ncol=3)+plot_annotation(tag_levels = 'a', 
                                                              tag_prefix = '(',
                                                              tag_suffix = ')'), 
       filename="figures/figure_pred-obs-gof_log-growth-l0-k-r.png",
       width=30,
       height=10,
       units='cm',
       dpi=350)
# END Figure *******************************************************************

# Plot r mod covars -----------------------------------------------
library(mgcViz)
p1 <- plot(pterm(getViz(r1),1))+
  l_fitPoints()+
  l_ciBar()+
  scale_x_discrete(limits=c("Sep","Oct","Nov","Dec","Jan","Feb"))+
  labs(x='Month of Fire', 
       y='Estimate'); p1
p2 <- plot(pterm(getViz(r1),2))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(log~bgroup("(",frac(L[0],K),")"))), 
       y='Estimate')
p3 <- plot(pterm(getViz(r1),3))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(VPD[MA]~bgroup("(",kPa,")"))), 
       y='Estimate')
p4 <- plot(pterm(getViz(r1),4))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(frac("Post Fire Precip. Anom."['12 mo'],
                               "Mean Annual Precip.")~("mm"))), 
       y='Estimate')
p_r <- mgcViz::gridPrint(p1,p2,p3,p4, ncol=2, 
                  top=grid::textGrob(expression(paste(bolditalic(r)~~"Model Covariates"))))

ggsave(p_r, 
       filename="figures/figure_mgcViz-covars_log-growth-r.png",
       width=20,
       height=18,
       units='cm',
       dpi=350)
# END Figure *******************************************************************


# Plot k mod covars -----------------------------------------------
library(mgcViz)
p1 <- plot(pterm(getViz(k1),1))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(sigma~LAI[yr])), 
       y='Estimate'); p1
p2 <- plot(pterm(getViz(k1),2))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(LAI[MA])), 
       y='Estimate')
p_k <- mgcViz::gridPrint(p2,p1,ncol=2, 
                         top=grid::textGrob(expression(paste(bolditalic(K)~~"Model Covariates"))))

ggsave(p_k, 
       filename="figures/figure_mgcViz-covars_log-growth-k.png",
       width=20,
       height=10,
       units='cm',
       dpi=350)
# END Figure *******************************************************************


# Plot L0/K mod covars -----------------------------------------------
p1 <- plot(pterm(getViz(ldk1),1))+
  l_fitPoints()+
  l_ciBar()+
  labs(x='NVIS Major Veg. Group', 
       y='Estimate')+
  theme(axis.text.x=element_text(angle=45,hjust=1)); p1

p2 <- plot(pterm(getViz(ldk1),2))+
  l_fitPoints()+
  l_ciBar()+
  scale_x_discrete(limits=c("Sep","Oct","Nov","Dec","Jan","Feb"))+
  labs(x='Month of Fire', 
       y='Estimate'); p2

p3 <- plot(sm(getViz(ldk1),1))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(LAI[MA])), 
       y='Estimate'); p3

p4 <- plot(sm(getViz(ldk1),2))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(sigma~LAI[yr])), 
       y='Estimate'); p4

p5 <- plot(sm(getViz(ldk1),3))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(bgroup("(",frac("Precip. Anom."["12 mo"], 
                                          "Precip. "[MA]),
                                 ")"))), 
       y='Estimate'); p5

p6 <- plot(sm(getViz(ldk1),4))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste("VPD Anom."["3 mo"])), 
       y='Estimate'); p6

p_ldk <- mgcViz::gridPrint(p1,p2,p5,p3,p4,p6,ncol=3, 
             top=grid::textGrob(expression(paste(bolditalic(frac(L[0],K))~~"Model Covariates"))))

ggsave(p_ldk, 
       filename="figures/figure_mgcViz-covars_log-growth-l0dk.png",
       width=30,
       height=20,
       units='cm',
       dpi=350)

# END figure *******************************************************************







dat[sample(.N,10000)] %>%
  # ggplot(data=.,aes((malai-delta)/malai, ldk))+
  ggplot(data=.,aes( (malai-delta)/malai ,ldk, color=post_precip_anom_12mo))+
  geom_point()+
  geom_smooth()+
  scale_color_gradient2()

expand_grid(ldk = seq(0.05,0.75,length.out=50), 
            precip_anom_frac = 0,
            post_precip_anom_frac = c(-0.4,0,0.4), 
            post_vpd15_anom_frac = c(-0.1,0,0.1), 
            malai = c(1,3,6), 
            lai_yr_sd = 0.5, 
            elevation = 0) %>% 
  mutate(pred = predict(r1, newdata=., type='response')) %>% 
  ggplot(data=.,aes(ldk, pred,color=malai,group=malai))+
  geom_line()+
  facet_grid(post_precip_anom_frac ~ post_vpd15_anom_frac, labeller = label_both)


# Scratch ----------------------------------------------------------------------
r1 <- ranger::ranger(ldk ~ 
                       aspect+
                       elevation+
                       hnd+
                       map+matmax+matmin+mavpd15+
                       precip_anom_12mo+
                       vpd15_anom_3mo,
                     data=dat[is.na(elevation)==F][is.na(aspect)==F],
                     importance='impurity_corrected')
vip::vip(r1)                     
r1


