library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(furrr)
library(arrow)
library(mgcv)
library(mgcViz)

# oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
#   sf::st_simplify(., dTolerance = 0.1) %>% 
#   select(NAME_1)

# fits --- 
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-22 09:37:07.parquet")
fits <- fits[isConv==TRUE][r2>0.3][L0<K][L0>0][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
]

dttr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS_2021-06-05 13:01:38.parquet")

# dominant species ---- 
out <- read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-14 16:45:26.parquet")
sout <- set_names(st_as_stars(out[,.(x,y,predict)], dims = c("x","y"),crs=4326),c("species"))
st_crs(sout) <- st_crs(4326)
vv <- st_extract(sout,
  st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))

# # species bioclim ranges ---
# hab <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_species-median-clim-topo.parquet")
# names(hab) <- c("species",paste0("hab_",names(hab)[-1]))

# mean annual climate -------------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y","map","mapet","mappet",
                                           "matmax","matmin","mavpd15")) %>% 
  distinct() %>% 
  as.data.table()
rclim <- st_as_stars(clim)
st_crs(rclim) <- st_crs(4326)



# merge into dat  ---------------------
# dat <- merge(fits,dom,by=c("id","x","y"))
dat <- fits
# dat <- merge(dat, setnames(hab, "species","dom_sp"), by=c("dom_sp"))

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


# pre & post_fire climate -------------------------
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
dat[,`:=`(vpd15_anom_frac = vpd15_anom_12mo/mavpd15)]

# Seasonal LAI vars -------------------------------------
smalai <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_mean-annual-lai.tif", 
                            proxy=F)
m1 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-1.tif", 
                        proxy=F)
m2 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-2.tif", 
                        proxy=F)
m3 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-3.tif", 
                        proxy=F)
m4 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-4.tif", 
                        proxy=F)
m5 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-5.tif", 
                        proxy=F)
m6 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-6.tif", 
                        proxy=F)
m7 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-7.tif", 
                        proxy=F)
m8 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-8.tif", 
                        proxy=F)
m9 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-9.tif", 
                        proxy=F)
m10 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-10.tif", 
                         proxy=F)
m11 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-11.tif", 
                         proxy=F)
m12 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-12.tif", 
                         proxy=F)
m_djf <- (m12+m1+m2)/3
m_mam <- (m3+m4+m5)/3
m_jja <- (m6+m7+m8)/3
m_son <- (m9+m10+m11)/3
slai <- c(smalai, m_djf, m_mam, m_jja, m_son)
names(slai) <- c("malai",
                 "lai-djf","lai-mam","lai-jja", 
                 "lai-son")
rm(smalai, m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12, 
   m_djf, m_mam, m_jja, m_son)
vv <- st_extract(slai,st_as_sf(dat[,.(x,y)],coords=c("x","y"),crs=st_crs(4326)))
dat <- bind_cols(dat,as.data.table(st_drop_geometry(vv)) %>% select(-'malai') %>% as.data.table())
tmp <- dat[,.(lai_range_median = range(c(`lai-djf`,`lai-mam`,`lai-jja`,`lai-son`))),by=.(x,y)]
dat <- merge(dat, tmp, by=c("x","y"))


# Filter dat to just Eucs in the Oct-Feb burning -------------------
dat <- dat[vc %in% c(2,3,5)][month %in% c(9,10,11,12,1,2)]
# nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
# sp_fac <- 
#   unique(dat[dom_sp %in% nobs[rank <= 30]$dom_sp][,.(dom_sp,hab_hnd)]) %>% 
#   .[order(hab_hnd)] %>% 
#   .[,sp_fac := forcats::fct_inorder(dom_sp)]
# Final mutations
dat[,delta:=K-L0]
dat[,fire_month_f := lubridate::month(date_fire1,label = TRUE,abbr = TRUE)][
  ,fire_month_f := factor(fire_month_f,
                          levels=c("Sep","Oct","Nov","Dec","Jan","Feb"),
                          ordered = TRUE)]#[,dom_sp_f := factor(dom_sp)]
dat[,vc_name_f := factor(vc_name)]
dat[,r365 := r*365]


# # Predict delta: How much leaf area was lost due to fire
# # delta = K-L0
# delta1 <- bam(delta~
#                 vc_name+
#                 fire_month_f+
#                 malai+
#                 lai_yr_sd+
#                 te(mavpd15,vpd15_anom_3mo,k=5)+
#                 te(map,precip_anom_frac,k=5),
#               family=Gamma(link='log'),
#               data=dat[sample(.N, 10000)], 
#               select=TRUE, 
#               discrete=TRUE)
# summary(delta1)
# plot(delta1, scheme=2, all.terms=F, pages=1)
# 
# # Predict l0: How much leaf area is left after the fire
# l01 <- bam(L0 ~
#                 vc_name+
#                 fire_month_f+
#                 malai+
#                 lai_yr_sd+
#                 te(mavpd15,vpd15_anom_3mo,k=5)+
#                 te(map,precip_anom_frac,k=5),
#               family=Gamma(link='log'),
#               data=dat[sample(.N, 10000)], 
#               select=TRUE, 
#               discrete=TRUE)
# summary(l01)
# plot(l01, scheme=2, all.terms=F, pages=1)
# 
# 
# # Predict recovery LAI (k: carrying capacity)
# k1 <- bam(K ~ 
#             s(malai),
#           data=dat[sample(.N,50000)], 
#           select=TRUE, 
#           discrete=TRUE)
# summary(k1)
# print(plot(getViz(k1),allTerms=TRUE),pages=1)
# plot(k1, scheme=2, all.terms=TRUE, pages=1)
# 
# Predict r: Akin to the growth rate in the linear segment of the logistic function
r1 <- bam(r ~
            fire_month_f+ 
            s(log(ldk))+
            # lai_yr_sd+
            # I(matmax-matmin) +
            s(malai,k=5,bs='cs') +
            # s(mapet, k=5, bs='cs')+
            s(mappet,k=5,bs='cs')+
            # I(post_tmax_anom_12mo/matmax)+
            # I(post_vpd15_anom_12mo/mavpd15) +
            # post_tmin_anom_12mo+
            s(post_vpd15_anom_frac,k=5,bs='cs') + 
            s(post_precip_anom_frac,k=5,bs='cs')
          , 
          data=dat, 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(r1)
print(plot(getViz(r1),allTerms=TRUE),pages=1)
# gratia::observed_fitted_plot(r1)+
#   geom_abline(col='red')
# plot(r1, scheme=2, pages=1,rug=TRUE,
#      all.terms=F)
# visreg::visreg(r1,xvar = 'mavpd15')

# Predict r: Akin to the growth rate in the linear segment of the logistic function
r2 <- bam(exp(r)~fire_month_f+ 
          te(ldk,malai,
             post_precip_anom_frac,k=5), 
          data=dat, 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(r2)
# plot(r2, scheme=2, pages=1,rug=TRUE, 
#      all.terms=TRUE)
# visreg::visreg(r2,xvar = 'mavpd15')

r3 <- bam(r ~
          s(mappet, k=5, bs='cs')+
          s(post_precip_anom_frac, k=5, bs='cs')+
          s(ldk, k=5, bs='cs')+
          fire_month_f +
          s(species, bs='re')
          , 
          data=dat, #[between(post_vpd15_anom_frac,-0.05,0.075)], 
          family=Gamma(link='log'),
          select=TRUE, 
          discrete=TRUE)
summary(r3)
# print(plot(getViz(r3),allTerms=TRUE),pages=1)


# Predictions -------------------------------------------------------------
pdat <- dat
# pdat$pred_delta <- predict(delta1, newdata=pdat, type='response')
# pdat$pred_k <- predict(k1, newdata=pdat, type='response')
# pdat$pred_l0 <- predict(l01, newdata=pdat, type='response')
# pdat$pred_ldk <- predict(ldk1, newdata=pdat, type='response')
pdat$pred_r <- log(predict(r1,
                           newdata=pdat,
                           type='response'))
# tmp <- pdat %>%
#   lazy_dt() %>% 
#   select(-ldk) %>% 
#   rename(ldk = pred_ldk) %>% 
#   as.data.table()
# pdat$pred2_r <- log(predict(r1, 
#                            newdata=tmp, 
#                            type='response'))



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
  ggplot(data=.,aes(exp(pred_r),r))+
  geom_density_2d_filled()+
  scale_fill_viridis_d(option='A')+
  geom_smooth(method='lm', 
              col='lightblue')+
  geom_abline(col='#CF0000',lty=1)+
  scale_x_continuous(expand=c(0,0), limits=c(0,0.01))+
  scale_y_continuous(expand=c(0,0), limits=c(0,0.01))+
  labs(x=expression(paste("Predicted r ",bgroup("(",frac(m**2,day),")"),"")), 
       y=expression(paste("r ",bgroup("(",frac(m**2,day),")"),"")))+
  # coord_cartesian(#xlim = c(0,0.02),
  #                 #ylim=c(0, 0.02)
  #                 )+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'none')
fig_r
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
p2 <- plot(sm(getViz(r1),1))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(log~bgroup("(",frac(L[0],K),")"))), 
       y='Estimate');p2
p3 <- plot(sm(getViz(r1),2))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(LAI[MA]~"(m²/m²)")), 
       y='Estimate');p3
p4 <- plot(sm(getViz(r1),3))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(P:PET[MA])), 
       y='Estimate');p4
p5 <- plot(sm(getViz(r1),4))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(frac("Post Fire VPD. Anom."['12 mo'],
                               "Mean Annual VPD."))), 
       y='Estimate');p5
p6 <- plot(sm(getViz(r1),5))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(frac("Post Fire Precip. Anom."['12 mo'],
                               "Mean Annual Precip."))), 
       y='Estimate');p6
p_r <- mgcViz::gridPrint(p1,p2,p3,p4,p5,p6, ncol=3, 
                  top=grid::textGrob(expression(paste(bolditalic(r)~~"Model Covariates"))))

ggsave(p_r, 
       filename="figures/fig-S7_r-gam-covars.png",
       width=25,
       height=15,
       units='cm',
       dpi=350)
# END Figure *******************************************************************




# SM Fig 8 L0/K ----------------------------------------------------------------
# beta regression approach
ldk1 <- bam(ldk~
              # s(species, bs='re')+
              # s(pre_fire_slai_anom_12mo,k=4,bs='cs')+
              s(I(pre_fire_slai_anom_3mo/malai),bs='cs',k=5)+
              s(precip_anom_frac, bs='cs',k=5)+
              # vc_name_f+
              # s(tmax_u,tmin_u, k=5)+
              fire_month_f+
              # s(lai_range_median,k=5,bs='cs')+
              # scale(malai)+
              # scale(mappet)+
              s(malai,k=5,bs='cs')+
              s(mappet,k=5,bs='cs'),
              # s(lai_yr_sd, k=3,bs='cs')+
              # s(precip_anom_frac,k=3,bs='cs'),
              # s(vpd15_anom_12mo,k=3,bs='cs'),
              # s(vpd15_anom_3mo,k=3,bs='cs'),
              # species+
              # te(vpd15_anom_3mo, mavpd15,k=5)+
              # s(mavpd15,bs='cs',k=5)+
              # scale(vpd15_anom_frac)+
              # scale(precip_anom_frac),
              # s(precip_anom_frac,bs='cs',k=5),
            data=dat,#[sample(.N, 50000)],
            family=betar(),
            select=TRUE,
            discrete=TRUE)
summary(ldk1)
plot(ldk1, scheme=2,all.terms=F,pages=1)
getViz(ldk1) %>% pterm(3) %>% plot %>% print(pages=1)
getViz(ldk1) %>% plot(allTerms=T) %>% print(pages=1)

p1 <- plot(pterm(getViz(ldk1),1))+
  l_fitPoints()+
  l_ciBar()+
  scale_x_discrete(limits=c("Sep","Oct","Nov","Dec","Jan","Feb"))+
  labs(x='Month of Fire', 
       y='Estimate'); p1
p2 <- plot(sm(getViz(ldk1),1))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(bgroup("",
    frac("LAI Anom."["3 mo"],
    "LAI"["MA"]),""))),
       y='Estimate')+
  coord_cartesian(xlim=c(-0.7,0.7));p2
p3 <- plot(sm(getViz(ldk1),2))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(bgroup("",
    frac("Precip. Anom."["12 mo"],
    "Precip."["MA"]),""))),
       y='Estimate')+
    coord_cartesian(xlim=c(-0.6,0.3), expand = F);p3
p4 <- plot(sm(getViz(ldk1),3))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste(LAI[MA])), 
       y='Estimate')+
  coord_cartesian(xlim=c(0.5,6));p4
p5 <- plot(sm(getViz(ldk1),4))+
  l_fitLine()+
  l_ciLine()+
  l_rug()+
  labs(x=expression(paste("P:PET"["MA"])), 
       y='Estimate')+
  coord_cartesian(xlim=c(0.5,3));p5
p_ldk <- mgcViz::gridPrint(p1,p2,p3,p4,p5,p6, ncol=3, 
                  top=grid::textGrob(expression(paste(
                    bolditalic(frac(L[0],K))~~"Model Covariates"))))
ggsave(p_ldk, 
       filename="figures/fig-S8_ldk-gam-covars.png",
       width=25,
       height=15,
       units='cm',
       dpi=350)


# How much longer is TTR by month? 
dat[,.(ttr_mw = median(ttr5_lai),
       ttr_lf = median(pred_ttr)), by=.(month)][
         ,.(ratio_mw = max(ttr_mw)/min(ttr_mw), 
           ratio_lf = max(ttr_lf)/min(ttr_lf))]

dat[,.(ttr_mw = median(ttr5_lai),
       ttr_lf = median(pred_ttr)), by=.(fire_month_f)][
         ,.(ratio_mw = max(ttr_mw)/min(ttr_mw), 
           ratio_lf = max(ttr_lf)/min(ttr_lf))]


# How much longer is TTR by fire year? 
dat[,.(ttr_mw = median(ttr5_lai),
       ttr_lf = median(pred_ttr)), by=.(month)][
         ,.(ratio_mw = max(ttr_mw)/min(ttr_mw), 
           ratio_lf = max(ttr_lf)/min(ttr_lf))]

dat[,`:=`(fire_year = year(date-months(3)))][,.(ttr_mw = median(ttr5_lai),
       ttr_lf = median(pred_ttr)), by=.(fire_year)][fire_year>2000][order(fire_year)][
         ,.(ratio_mw = max(ttr_mw)/min(ttr_mw), 
           ratio_lf = max(ttr_lf)/min(ttr_lf))]
