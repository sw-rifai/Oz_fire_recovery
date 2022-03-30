library(tidyverse); 
# library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(furrr)
library(arrow)
library(mgcv)
library(mgcViz)
library(gratia)
library(patchwork)

# oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
#   sf::st_simplify(., dTolerance = 0.1) %>% 
#   select(NAME_1)

# fits --- 
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-16 09:31:56.parquet")
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
post_clim <- post_clim[order(x,y,date)][, 
vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                     n = 3,fill = NA,align='right'), by=.(x,y)]
post_clim <- post_clim[order(x,y,date)][, 
precip_anom_3mo := frollapply(precip_anom,FUN=mean,
                                               n = 3,fill = NA,align='right'), by=.(x,y)]
# # 24 month
# post_clim <- post_clim[order(x,y,date)][, precip_anom_24mo := frollsum(precip_anom,n = 24,fill = NA,align='right'), by=.(x,y)]
# # 36 month
# post_clim <- post_clim[order(x,y,date)][, precip_anom_36mo := frollsum(precip_anom,n = 36,fill = NA,align='right'), by=.(x,y)]

#  Attach clim pixel id to VI ------------------------------------
coords_vi <- dat %>% select(x,y,id) %>% distinct() %>% as.data.table()
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
dat[,`:=`(bs = 1 - (min_slai_anom+slai_u)/slai_u)]
dat[,`:=`(bs = ifelse(bs>1,1,bs))]

pre_fire_slai_range <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/pre-fire_slai_range.parquet")

dat[,`:=`(frac_pf_slai_anom_3mo = pre_fire_slai_anom_3mo/malai, 
  frac_pf_slai_anom_12mo = pre_fire_slai_anom_12mo/malai)]
dat[,`:=`(frac_p12_anom = precip_anom_12mo/map,
         frac_vpd3_anom = vpd15_anom_3mo/mavpd15)]

# burn severity model -----------------------------------------------------
dat_bs <- merge(dat[between(pre_fire_slai_anom_12mo/malai,-1,1)][
                bs>0
              ][between(frac_pf_slai_anom_3mo,quantile(dat$frac_pf_slai_anom_3mo,c(0.01,0.99),na.rm=T)[1],quantile(dat$frac_pf_slai_anom_3mo,c(0.01,0.99),na.rm=T)[2])], 
pre_fire_slai_range,by=c("x","y","id"))
dat_bs[,`:=`(bs = ifelse(bs>1,1,bs))]


tmp1 <- dat_bs[sample(.N, 50000)]
dat_bs_eval <- data.table::fsetdiff(dat_bs,tmp1)
dat_bs <- tmp1

# dat_bs[sample(.N,1000)] %>% ggplot(data=., aes(slai_range,lai_yr_sd))+
#   geom_point()

bs1 <- bam(
  bs ~ # fire_month_f + 
   # I(1-(min_slai_anom+slai_u)/slai_u) ~
                # fire_month_f+
                # s(malai)+
                # s(I(vpd15_anom_12mo/mavpd15),k=5)+
                s(aspect,bs='cc')+
                s(hnd,k=5)+
    # s(malai,k=5)+
                # s(species,bs='re')+
                # scale(pH)+scale(elevation)+
                # s(I(vpd15_anom_3mo/mavpd15),k=5)+
                s(I(precip_anom_3mo/precip_u),k=5)+
                s(frac_pf_slai_anom_3mo,k=5)+
                s(frac_pf_slai_anom_12mo,k=5)+
                # s(pre_fire_slai_anom_3mo,k=5)+
                # s(pre_fire_slai_anom_12mo,k=5)+
                # s(matmin,k=5)+
                s(mapet,k=5),
                # s(pre_fire_slai_anom_3mo,slai_u)+
                # s(precip_anom_12mo,k=5),
                # te(mavpd15,vpd15_anom_3mo,k=5)+
                # te(map,precip_anom_12mo,k=5),
              family=gaussian(link='identity'),
              # family = betar(),
              data=dat_bs, 
              select=TRUE, 
              discrete=TRUE)
summary(bs1)
plot(bs1, scheme=2, all.terms=T, pages=1,hcolors=scico::scico(10,palette='roma',direction = -1),scale=0,rug=T)
dat_bs_eval %>% 
  mutate(pred = predict(bs1,newdata=.,type='response')) %>%
  summarize(r2 = yardstick::rsq_trad_vec(bs,pred),
            mae = yardstick::mae_vec(bs,pred))


# r365 model --------------------------------------------------------------
dat_r <- dat[between(pre_fire_slai_anom_12mo/malai,-1,1)][
    r<quantile(dat$r,0.99)
  ][r2>0.1]

tmp1 <- dat_r[sample(.N, 50000)]
dat_r_eval <- data.table::fsetdiff(dat_r,tmp1)
dat_r <- tmp1

r1 <- bam(r365~
    s(bs,k=5)+
    # s(I((min_slai_anom+slai_u)/slai_u),k=5) +
    s(frac_pf_slai_anom_3mo,k=5)+
    s(mapet,k=5)+
    fire_month_f+
    # s(I(post_pet_anom_12mo/mapet),k=5)+
    s(post_precip_anom_frac,k=5)+ # seems real
    s(log(mappet),k=5),
  family=Gamma(link='identity'),
  # method='fREML',
  data=dat_r, 
  select=T,
  discrete=TRUE
  )
summary(r1)
plot(r1, scheme=2, all.terms=T, pages=1,hcolors=scico::scico(10,palette='roma',direction = -1),rug=T)

dat_r_eval %>% 
  mutate(pred = predict(r1,newdata=.,type='response')) %>% 
  summarize(r2 = yardstick::rsq_trad_vec(r365,pred),
            mae = yardstick::mae_vec(r365,pred))


# TTR5 model --------------------------------------------------------------
dat_ttr <- dat[between(pre_fire_slai_anom_12mo/malai,-1,1)][
    between(frac_pf_slai_anom_3mo,-1,1)
  ]
tmp1 <- dat_ttr[sample(.N, 50000)]
dat_ttr_eval <- data.table::fsetdiff(dat_ttr,tmp1)
dat_ttr <- tmp1

t1 <- bam(
  ttr5_lai ~
    s( bs, k=5)+
    # s(I((min_slai_anom+slai_u)/slai_u),k=5)+
    fire_month_f+
    # s(malai,k=5,m=c(1,-1),bs='gp')+
    # s(sand,k=5)+
    # s(pH,k=5)+
    s(post_precip_anom_frac,k=5)+
    # s(frac_pf_slai_anom_3mo,k=5)+
    s(frac_pf_slai_anom_12mo,k=5)+
    # s(mapet,k=5)+
    # s(log(mappet),k=5)+
    s(mavpd15,k=5)+
    # scale(I(post_vpd15_anom_12mo/mavpd15))+
    # s(malai) +
    # s(pre_fire_slai_anom_3mo,slai_u)+
    # te(mavpd15,vpd15_anom_3mo,k=5)+
    # te(map,precip_anom_12mo,k=5)
    s(frac_p12_anom,k=5
      )
  ,
  family=Gamma(link='identity'),
  data=dat_ttr, 
  select=T, 
  discrete=TRUE)
summary(t1)
plot(t1, scheme=2, all.terms=T, pages=1,hcolors=scico::scico(10,palette='roma',direction = -1),rug=T)

dat_ttr_eval %>% 
  mutate(pred = predict(t1,newdata=.,type='response')) %>% 
  summarize(r2 = yardstick::rsq_trad_vec(ttr5_lai,pred),
            mae = yardstick::mae_vec(ttr5_lai,pred))


# Plotting ----------------------------------------------------------------
## Plot burn severity -----
library(gratia)
summary(bs1)
p_bs_1 <-  draw(bs1,select="frac_pf_slai_anom_3mo",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
  labs(
    y="Effect",
    title=NULL,
       x = expression(paste("(Pre-fire ",LAI["3 mo. anom"],"):",
         LAI["normal"])))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_bs_1

p_bs_2 <- draw(bs1,select="frac_pf_slai_anom_12mo",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         x = expression(paste("(Pre-fire"~LAI["12 mo. anom"],"):",
           LAI["normal"])))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_bs_2

p_bs_3 <- draw(bs1,select="precip_anom_3mo",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         x = expression(paste("(Pre-fire"~P["3 mo. anom"],"):",
           P["normal"])))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_bs_3

p_bs_4 <- draw(bs1,select="mapet",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  # geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y='Effect',
    # y = expression(paste("Effect on ",
    # LAI["post-fire"]:LAI["norm."])),
    title=NULL,
         x = expression(paste("Mean Annual PET"~(mm~yr**-1))))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_bs_4

p_bs_5 <- draw(bs1,select="aspect",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  # geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         x = expression(paste(Aspect~(degrees))))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_bs_5

p_bs_6 <- draw(bs1,select="hnd",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey')+
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         x = expression(paste("Height Above Nearest Drainage (m)")))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_bs_6

# ((p_bs_1|p_bs_2|p_bs_3)+plot_layout(ncol=3))/
# ((p_bs_4|p_bs_5|p_bs_6)+plot_layout(ncol=3))

# ggsave(filename = "figures/fig_S-X_univar-gam-smooth_burn-severity.png",
#   width=25,
#   height=15,
#   units='cm',
#   dpi=350)
# END section ********************************************


## Plot growth rate 365 --------------------------------
summary(r1)
p_r1_1 <- draw(r1,select="mapet",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  # geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y='Effect',
      # y = expression(paste("Effect on ",
      #   bolditalic("r ")~(LAI~yr**-1)
      #   )),
    title=NULL,
         x = expression(paste("Mean Annual PET"~(mm~yr**-1))))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_r1_1



p_r1_2 <-  draw(r1,select="frac_pf_slai_anom_3mo",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
  # smooth_estimates(bs1, smooth="frac_pf_slai_anom_3mo", partial_match = T) %>%   # add_confint() %>% 
  # ggplot(aes(y = est, x = frac_pf_slai_anom_3mo)) +
  #   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
  #               alpha = 0.2, fill = "black") +
  #   geom_line(colour = "black", size = 1.5) +
  labs(y=NULL,
    title=NULL,
       # title = expression("Partial effect of" ~ "Aspect"),
       x = expression(paste("(Pre-fire"~LAI["3 mo. anom"],"):",
         LAI["normal"])))+
  scale_x_continuous(expand=c(0,0), 
    # limits=c(-0.8,0.75)
    )+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_r1_2

p_r1_3 <- draw(r1,select="post_precip_anom_frac",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         x = expression(paste("(Post-fire"~P["12 mo. anom"],"):",
           P["Mean Annual"])))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_r1_3

p_r1_4 <- draw(r1,select="s(log(mappet))",partial_match = F) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y="Effect",
      # y = expression(paste("Effect on ",
      #   bolditalic("r ")~(LAI~yr**-1)
      #   )),
    title=NULL,
         x = expression(paste(Log~bgroup("(","Mean Annual P:PET",")"))))+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_r1_4

p_r1_5 <- draw(r1,select='bs',partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey')+
  geom_vline(xintercept=0,lty=2,col='grey')+
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         x =expression(paste("Burn Severity ",
      "(",1 - LAI["post-fire"],":",LAI["normal"],")"))) +
  scale_x_continuous(
    # expand=c(0,0),
    # limits=c(0,1.01)
    )+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_r1_5


p_r1_6 <- evaluate_parametric_term(r1,'fire_month_f') %>% 
  group_by(value) %>% 
  summarize(partial = median(partial),
            lower=median(partial-1.96*se),
            upper=median(partial+1.96*se)) %>% 
  ungroup() %>% 
  ggplot(data=., aes(value,partial))+
  geom_pointrange(aes(ymin=lower,ymax=upper))+
    labs(
      y=NULL,
      x="Month of Fire",
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         # x = expression(paste(Aspect~(degrees)))
      )+
  # scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_r1_6


# ((p_r1_1|p_r1_2|p_r1_3)+plot_layout(ncol=3))/
# ((p_r1_4|p_r1_5|p_r1_6)+plot_layout(ncol=3))

# ggsave(filename = "figures/fig_S-X_univar-gam-smooth_lai-recovery-rate.png",
#   width=25,
#   height=18,
#   units='cm',
#   dpi=350)
# END section ********************************************

## Plot growth rate 365 --------------------------------
summary(t1) 
# (1) bs, (2) frac_p12_anom, (3) frac_pf_slai_anom_12mo
# (4) mavpd15, (5) post_precip_anom_frac, (6) fire_month_f
p_t1_1 <- draw(t1,select="bs",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  # geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
         x =expression(paste("Burn Severity ",
      "(",1 - LAI["post-fire"],":",LAI["normal"],")")), 
      y='Effect',
      # y = expression(paste("Effect on ",
      #   "TTR "~(days)
      #   )),
    title=NULL,
      )+
  scale_x_continuous(expand=c(0,0), limits=c(0,0.95))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_t1_1



p_t1_2 <-  draw(t1,select="frac_p12_anom",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
  # smooth_estimates(bs1, smooth="frac_pf_slai_anom_3mo", partial_match = T) %>%   # add_confint() %>% 
  # ggplot(aes(y = est, x = frac_pf_slai_anom_3mo)) +
  #   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
  #               alpha = 0.2, fill = "black") +
  #   geom_line(colour = "black", size = 1.5) +
    labs(
      y=NULL,
      x=expression(paste("(Pre-fire ",P["12 mo. anom."],"):",P["Mean Annual"])),
      # y = expression(paste("Effect on ",
      #   "TTR "~(days)
      #   )),
    title=NULL,
      )+
  scale_x_continuous(expand=c(0,0), 
    # limits=c(-0.8,0.75)
    )+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_t1_2

p_t1_3 <- draw(t1,select="frac_pf_slai_anom_12mo",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         x = expression(paste("(Pre-fire"~LAI["12 mo. anom"],"):",
         LAI["normal"]))
      )+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_t1_3

p_t1_4 <- draw(t1,select="mavpd15",partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey') + 
  geom_vline(xintercept=0,lty=2,col='grey') + 
    labs(
      y='Effect',
      # y = expression(paste("Effect on ",
      #   "TTR"~(days)
      #   )),
    title=NULL,
      x="Mean Annual VPD (kPa)"
      )+
  scale_x_continuous(
    limits=quantile(dat$mavpd15,c(0,0.999),na.rm=T),
    expand=c(0.01,0)
    )+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_t1_4

p_t1_5 <- draw(t1,select='post_precip_anom_frac',partial_match = T) +
  geom_hline(yintercept=0,lty=2,col='grey')+
  geom_vline(xintercept=0,lty=2,col='grey')+
    labs(
      y=NULL,
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
        x = expression(paste("(Post-fire"~P["12 mo. anom"],"):",
           P["Mean Annual"]))
      )+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_t1_5


p_t1_6 <- evaluate_parametric_term(t1,'fire_month_f') %>% 
  group_by(value) %>% 
  summarize(partial = median(partial),
            lower=median(partial-1.96*se),
            upper=median(partial+1.96*se)) %>% 
  ungroup() %>% 
  ggplot(data=., aes(value,partial))+
  geom_pointrange(aes(ymin=lower,ymax=upper))+
    labs(
      y=NULL,
      x="Month of Fire",
      # y = expression(paste("Partial effect on ",
      # frac(LAI["post-fire"],LAI["normal"]))),
    title=NULL,
         # x = expression(paste(Aspect~(degrees)))
      )+
  # scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_t1_6


# ((p_t1_1|p_t1_2|p_t1_3)+plot_layout(ncol=3))/
# ((p_t1_4|p_t1_5|p_t1_6)+plot_layout(ncol=3))

# ggsave(filename = "figures/fig_S-X_univar-gam-smooth_ttr5-lai.png",
#   width=25,
#   height=18,
#   units='cm',
#   dpi=350)




# Joint figure ------------------------------------------------------------

height_dim <- 12
width_dim <- 30
p_top <- (((p_bs_1|p_bs_2|p_bs_3)+plot_layout(ncol=3))/
((p_bs_4|p_bs_5|p_bs_6)+plot_layout(ncol=3))  +
  plot_annotation(
  tag_prefix = '(', tag_levels = list(letters[1:6]),tag_suffix = ')',
 title = 
      expression(paste('Partial effects of the drivers of Burn Severity ',
        (1-~LAI["post-fire"]:LAI["norm"])))
    )) 
ggsave(p_top,
  filename = "figures/tmp1.png",
    width=width_dim,
  height=height_dim,
  units='cm',
  dpi=300)

p_mid <- (((p_r1_1|p_r1_2|p_r1_3)+plot_layout(ncol=3))/
((p_r1_4|p_r1_5|p_r1_6)+plot_layout(ncol=3))+
  plot_annotation(
  tag_prefix = '(', tag_levels = list(letters[7:12]),tag_suffix = ')',
 title = expression(paste('Partial effects of the drivers of Post-fire LAI Recovery Rate ',
   (LAI~yr**-1)))))
ggsave(p_mid,
  filename = "figures/tmp2.png",
    width=width_dim,
  height=height_dim,
  units='cm',
  dpi=300)

p_bot <- ((p_t1_1|p_t1_2|p_t1_3)+plot_layout(ncol=3))/
((p_t1_4|p_t1_5|p_t1_6)+plot_layout(ncol=3)) + 
  plot_annotation(
  tag_prefix = '(', tag_levels = list(letters[13:18]),tag_suffix = ')',
title = 'Partial effects of the drivers of Time to Recover (days)')
ggsave(p_bot,
  filename = "figures/tmp3.png",
    width=width_dim,
  height=height_dim,
  units='cm',
  dpi=300)
 
library(magick)
tmp1 <- image_read("figures/tmp1.png")
tmp2 <- image_read("figures/tmp2.png")
tmp3 <- image_read("figures/tmp3.png")

image_out <- image_append(c(tmp1,tmp2,tmp3),stack=T)
magick::image_write(image_out,path = "figures/Fig_SX_uniResponseGAMS_BS_r_TT_v2.png",format = 'png')

