library(tidyverse)
library(data.table)
library(lubridate)
library(stars)
library(patchwork)
library(magic)

# Figures -------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 1000) %>% 
  select(NAME_1)

dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")

r_burned <- dat[is.na(fire_doy)==F][,.(x,y,date)] %>% 
  unique %>% 
  .[,`:=`(date_dec = decimal_date(date))] %>% 
  .[,.(val = median(date_dec), 
       nobs = .N), 
    by=.(x,y)]

# pals::parula(n=3)
# pals::pal.bands(pals::parula(10))
parula10 <- pals::parula(10)
col_kilmore <- parula10[3]
col_act <- parula10[6]
col_moogem <- parula10[10]

# Prep Kilmore East fire --------------------------------------------------
kilmore <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/kilmore_east_LT05_L1TP_092086_20090216_20161029_01_T1.tif")
p_kilmore <- RStoolbox::ggRGB(as(kilmore,Class = 'Raster'),r = 5, g=4, b=3, 
                              stretch='none', 
                              limits=c(100,4500))
p_kilmore <- p_kilmore+
  labs(x=NULL,y=NULL)+
  geom_label(aes(
    label='Kilmore E.',
    fontface='bold',
    x=145.07+0.05,
    y=-37.485-0.025    
  ), 
  label.size=NA,
  color = col_kilmore,
  # color='#c316d9',
  fill='grey10',
  alpha=0.95,
  size=8,
  )+
  theme_linedraw()+
  scale_x_continuous(expand=c(0,0),limits=c(145,145.45))+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        plot.background = element_rect(fill='grey10'))
p_kilmore


# Prep ACT 2003 fire ------------------------------------------------------
act <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/act_LE07_L1TP_090085_20030415_20170125_01_T1.tif")
p_act <- RStoolbox::ggRGB(as(act,Class = 'Raster'),r = 5, g=4, b=3)
# p_act
p_act <- p_act+
  labs(x=NULL,y=NULL)+
  geom_label(
           aes(label='ACT',
           fontface='bold',
           x=148.725+0.025,
           y=-35.586 + 0.025),
           label.size=NA,
           alpha=0.75,
           fill='black',
           color=col_act,
           # color='#12e1fc',
           size=8)+
  # annotate(geom='text',
  #          label='ACT',
  #          fontface='bold',
  #          color='#12e1fc',
  #          size=8,
  #          x=148.725+0.025,
  #          y=-35.586)+
  theme_linedraw()+
  scale_x_continuous(expand=c(0,0),
                     # limits=c(145,145.45)
                     )+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(-35.6,-35.3))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        plot.background = element_rect(fill='grey10'))
p_act

# Prep moogem 2002 fire ------------------------------------------------------
moogem <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/moogem_LE07_L1TP_089081_20021030_20170127_01_T1.tif")
p_moogem <- RStoolbox::ggRGB(as(moogem,Class = 'Raster'),
                             r = 5, g=4, b=3,
                             # stretch = 'hist',
                             # quantiles=c(0.1,0.9), 
                             limits=c(100,4500))
# p_moogem
p_moogem <- p_moogem+
  labs(x=NULL,y=NULL)+
  geom_label(
         aes(
             label='Moogem',
           fontface='bold',
           x=152.215+0.02,
           y=-29.57), 
         label.size=NA,
         fill='grey10',
         alpha=0.75,
         color=col_moogem,
         size=8)+
  theme_linedraw()+
  scale_x_continuous(expand=c(0,0),
                     # limits=c(145,145.45)
  )+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        plot.background = element_rect(fill='grey10'))
p_moogem


p_main <- r_burned %>% 
  ggplot(data=.,aes(x,y,fill=factor(nobs)))+
  geom_sf(data=oz_poly, inherit.aes = F, 
          fill='gray90',color='gray30')+
  geom_raster(data=unique(dat[date==ymd("2015-01-01")][,.(x,y)]), 
    aes(x,y),fill='grey50')+
  geom_tile()+
  coord_sf(xlim = c(144.5,153.75),
           ylim = c(-39.1,-27.25), expand = FALSE)+
  ggspatial::annotation_scale(location='bl')+
  scale_fill_manual(limits=factor(c(0,1,2,3)),
                    values=c("grey",viridis::rocket(3,end=0.85)), 
                    labels=c('0','1','2','≥3'))+
  # scale_fill_gradientn(colors=c("grey50",viridis::inferno(n = 3,end=0.9)),
  #                      # end=0.9,
  #                      limits=c(0,3),
  #                      breaks=c(0,1,2,3),
  #                      labels=c('0','1','2','3+'),
  #                      oob=scales::squish)+
  labs(x=NULL, y=NULL, fill='burns')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill='lightblue'),
        axis.text = element_text(size=17),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        legend.position = c(1,0), 
        legend.justification = c(0.99,0.01)); p_main

p_main <- p_main+
  geom_rect(aes(xmin=st_bbox(kilmore)$xmin, 
                   xmax=st_bbox(kilmore)$xmax,
                   ymin=st_bbox(kilmore)$ymin,
                   ymax=st_bbox(kilmore)$ymax),
               fill=NA,lwd=3.5,
               col='black')+
  geom_rect(aes(xmin=st_bbox(act)$xmin, 
                xmax=st_bbox(act)$xmax,
                ymin=st_bbox(act)$ymin,
                ymax=st_bbox(act)$ymax),
            fill=NA,lwd=3.5,
            col='black')+
  geom_rect(aes(xmin=st_bbox(moogem)$xmin, 
                xmax=st_bbox(moogem)$xmax,
                ymin=st_bbox(moogem)$ymin,
                ymax=st_bbox(moogem)$ymax),
            fill=NA,lwd=3.5,
            col='black')+
  geom_rect(aes(xmin=st_bbox(kilmore)$xmin, 
                     xmax=st_bbox(kilmore)$xmax,
                     ymin=st_bbox(kilmore)$ymin,
                     ymax=st_bbox(kilmore)$ymax),
                 fill=NA,lwd=2,
                 # col='#c316d9'
                 col=col_kilmore,
    )+
  geom_rect(aes(xmin=st_bbox(act)$xmin, 
                xmax=st_bbox(act)$xmax,
                ymin=st_bbox(act)$ymin,
                ymax=st_bbox(act)$ymax),
            fill=NA,lwd=2,
                     col=col_act,
            # col='#12e1fc'
    )+
  geom_rect(aes(xmin=st_bbox(moogem)$xmin, 
                xmax=st_bbox(moogem)$xmax,
                ymin=st_bbox(moogem)$ymin,
                ymax=st_bbox(moogem)$ymax),
            fill=NA,lwd=2,
            col=col_moogem
    )

p_ss <- p_moogem/p_act/p_kilmore
# p_ss
p_map <- p_main+inset_element(p_ss,
                              clip=TRUE,
                              align_to = 'panel',
                              left = 0,
                              bottom = 0.35,
                              right=0.425,
                              top=1 #left bottom right top
                              ) 

ggsave(p_map, 
       filename = "figures/p_map.png",
       height=35, 
       width=23,
       units='cm',
       dpi = 350)

# END map plot *******************************************************

# cleanup --- 
rm(dat); gc(full=T)

# Plot time series --------------------------------------------------------
# Data import ---------------------------------------------------
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
sdat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS_2021-06-05 13:01:38.parquet")
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-05-19 18:18:29.parquet")
d_soil <- arrow::read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")

gc(full=T)
# calculate the rolling metrics ------------------------------------------------
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_3mo := frollapply(precip_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]


#! Only fitting locations where the recovery was at least one year
sdat <- sdat[is.na(ttr5_lai)==FALSE][date_fire1<ymd('2015-01-01')][ttr5_lai>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_fire1,ttr5_lai)], 
               by=c("x","y","id"))
gc(full=T)

mdat <- ssdat[ssdat[, .I[date <= date_fire1 + days(2030)], by = .(x, 
    y, id)]$V1][, `:=`(post_days = as.double(date - date_fire1))]

# mdat <- ssdat %>% 
#   lazy_dt() %>%
#   mutate(recovery_date = date_fire1+days(ttr5_lai)) %>% 
#   group_by(x,y,id) %>% 
#   filter(date >= date_fire1) %>% 
#   filter(date <= date_fire1+days(2030)) %>% 
#   ungroup() %>%
#   mutate(post_days = as.double(date - date_fire1)) %>% 
#   as.data.table() 

rm(dat); gc(full=TRUE)
fits <- merge(fits,d_soil,by='id')


# STAGE 2: Attach AWAP pixel id to VI ------------------------------------
# coords_vi <- lazy_dt(fits) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- unique(fits[,.(x,y,id)])
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
fits <- merge(fits, coords_vi, by='id')
gc(full=TRUE)
#*************************************************************************

# ACT quantile plot -------------------------------------------------------
vec_act <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
     between(x,148.7,148.93)][
    between(y,-35.6,-35.315)][,rank:=frank(ttr5_lai,ties.method = 'first')]

clim_act <- clim[, .(val = mean(precip_anom_3mo, na.rm = TRUE) * 3), keyby = .(date)][, 
    `:=`(post_days = as.double(date - median(vec_act$date_fire1)))]

# clim_act <- clim[idx_awap %in% (vec_act$idx_awap %>% unique)] %>% 
#   lazy_dt() %>%
#   group_by(date) %>% 
#   summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
#   ungroup() %>% 
#   mutate(post_days = as.double(date-median(vec_act$date_fire1))) %>% show_query()
#   as.data.table()


q_act <- mdat[id %in% vec_act$id][
      ,post_days := round(post_days/20)*20
    ][, 
      .(q_mid = median(slai_anom, na.rm=T), 
        q_hi = quantile(slai_anom,0.975, na.rm=T), 
        q_lo = quantile(slai_anom,0.025), 
        m_date = median(date, na.rm=T),
        nobs=.N), 
      by='post_days']

vec_hi_low <- floor(quantile(vec_act$rank,c(0.025,0.95)))
traj_act <- rbindlist(
  list(mdat[id==vec_act[rank==vec_hi_low[1]]$id][,recov:='fast'], 
       mdat[id==vec_act[rank==vec_hi_low[2]]$id][,recov:='slow']))

p_act <- q_act[nobs>10][post_days<=2000] %>% 
  ggplot(data=., aes(post_days, q_mid))+
  geom_rect(data=clim_act[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=T),
                fill=val))+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_ribbon(aes(ymin=q_lo, 
                  ymax=q_hi), 
              alpha=0.5, 
              fill='grey50')+
  geom_line(data=traj_act[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov), 
            lwd=3,alpha=0.25)+
  geom_line(data=traj_act[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov))+
  geom_line(lwd=1)+
  # geom_line(data=traj_act[recov=='fast'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov),
  #           color='#542788',
  #           alpha=0.25,
  #           lwd=2,
  #           lty=1)+
  # geom_line(data=traj_act[recov=='fast'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov),
  #           color='#542788',
  #           lty=1)+
  # geom_line(data=traj_act[recov=='slow'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov), 
  #           color='#b35806',
  #           lwd=2,
  #           alpha=0.25,
  #           lty=1)+
  # geom_line(data=traj_act[recov=='slow'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov), 
  #           color='#b35806',
  #           lty=1)+
  scale_color_manual(values = c("fast"='#542788', 
                                "slow"='#b35806'), 
                     breaks = c("fast","slow"), 
                     labels = c("Fast","Slow")
                        )+
  scale_fill_gradient2(limits=c(-150,150),
                       breaks=c(-150,0,150),
                       labels=c("≤ -150 ",  "0"," ≥ +150"),
                       oob=scales::squish)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2000), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(#x='Days after fire', 
       x=NULL,
       y='LAI Anomaly (m²/m²)', 
       color='Recovery', 
       title='ACT', 
       fill="Precip. Anom. 3-mo (mm)     ")+
  # facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
  #            ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'none'); p_act
# END ACT ******************************************************************


# Kilmore E quantile plot -----------------------------------------------------------
vec_ke <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
  between(x,145.0,145.335)][
    between(y,-37.45,-37.27)][,rank:=frank(ttr5_lai,ties.method = 'first')]

clim_ke <- clim[, .(val = mean(precip_anom_3mo, na.rm = TRUE) * 3), keyby = .(date)][, 
    `:=`(post_days = as.double(date - median(vec_ke$date_fire1)))][]

# clim_ke <- clim[idx_awap %in% (vec_ke$idx_awap %>% unique)] %>% 
#   lazy_dt() %>%
#   group_by(date) %>% 
#   summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
#   ungroup() %>% 
#   mutate(post_days = as.double(date-median(vec_ke$date_fire1))) %>% show_query()
#   as.data.table()

q_ke <- mdat[id %in% fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
  between(x,145.0,145.335)][
    between(y,-37.45,-37.27)]$id][
      ,post_days := round(post_days/20)*20
    ][, 
      .(q_mid = median(slai_anom, na.rm=T), 
        q_hi = quantile(slai_anom,0.975, na.rm=T), 
        q_lo = quantile(slai_anom,0.025), 
        m_date = median(date, na.rm=T),
        nobs=.N), 
      by='post_days']

vec_hi_low <- floor(quantile(vec_ke$rank,c(0.025,0.95)))
traj_ke <- rbindlist(
  list(mdat[id==vec_ke[rank==vec_hi_low[1]]$id][,recov:='fast'], 
       mdat[id==vec_ke[rank==vec_hi_low[2]]$id][,recov:='slow']))

p_ke <- q_ke[nobs>10][post_days<=2000] %>% 
  ggplot(data=., aes(post_days, q_mid))+
  geom_rect(data=clim_ke[post_days>= 0][post_days<=2000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=T),
                fill=val))+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_ribbon(aes(ymin=q_lo, 
                  ymax=q_hi), 
              alpha=0.5, 
              fill='grey50')+
  geom_line(data=traj_ke[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov), 
            lwd=3,alpha=0.25)+
  geom_line(data=traj_ke[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov))+
  geom_line(lwd=1)+
  scale_color_manual(values = c("fast"='#542788',
                                "slow"='#b35806'),
                     breaks = c("fast","slow"),
                     labels = c("Fast","Slow")
  )+
  scale_fill_gradient2(limits=c(-150,150),
                       breaks=c(-150,0,150),
                       labels=c("≤ -150 ",  "0"," ≥ +150"),
                       oob=scales::squish)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2000), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='LAI Anomaly (m²/m²)', 
       color='Recovery', 
       title='Kilmore East', 
       fill="Precip. Anom. 3-mo (mm)     ")+
  guides(color=guide_legend(title.position = 'top'),
         fill=guide_colorbar(title.position = 'top'))+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
    legend.position = 'bottom'); p_ke
# END **************************************************************************

# Moogem quantile plot -----------------------------------------------------------
vec_mo <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
    between(x,152.1842,152.382)][
    between(y,-29.581,-29.439)][,rank:=frank(ttr5_lai,ties.method = 'first')]

clim_mo <- clim[, .(val = mean(precip_anom_3mo, na.rm = TRUE) * 3), keyby = .(date)][, 
    `:=`(post_days = as.double(date - median(..vec_mo$date_fire1)))]
# clim_mo <- clim[idx_awap %in% (vec_mo$idx_awap %>% unique)] %>% 
#   lazy_dt() %>%
#   group_by(date) %>% 
#   summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
#   ungroup() %>% 
#   mutate(post_days = as.double(date-median(vec_mo$date_fire1))) %>% show_query()
#   as.data.table()

q_mo <- mdat[id %in% vec_mo$id][
      ,post_days := round(post_days/20)*20
    ][, 
      .(q_mid = median(slai_anom, na.rm=T), 
        q_hi = quantile(slai_anom,0.975, na.rm=T), 
        q_lo = quantile(slai_anom,0.025), 
        m_date = median(date, na.rm=T),
        nobs=.N), 
      by='post_days']

vec_hi_low <- floor(quantile(vec_mo$rank,c(0.025,0.95)))
traj_mo <- rbindlist(
  list(mdat[id==vec_mo[rank==vec_hi_low[1]]$id][,recov:='fast'], 
       mdat[id==vec_mo[rank==vec_hi_low[2]]$id][,recov:='slow']))

(p_mo <- q_mo[nobs>10][post_days<=2000] %>% 
  ggplot(data=., aes(post_days, q_mid))+
  geom_rect(data=clim_mo[post_days>= 0][post_days<=2000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=T),
                fill=val))+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_ribbon(aes(ymin=q_lo, 
                  ymax=q_hi), 
              alpha=0.5, 
              fill='grey50')+
    geom_line(data=traj_mo[post_days<ttr5_lai], 
              aes(post_days, slai_anom, color=recov), 
              lwd=3,alpha=0.25)+
    geom_line(data=traj_mo[post_days<ttr5_lai], 
              aes(post_days, slai_anom, color=recov))+
    geom_line(lwd=1)+
    scale_color_manual(values = c("fast"='#542788',
                                  "slow"='#b35806'),
                       breaks = c("fast","slow"),
                       labels = c("Fast","Slow")
    )+
  scale_fill_gradient2(limits=c(-150,150),
                       breaks=c(-150,0,150),
                       labels=c("≤ -150 ",  "0"," ≥ +150"),
                       oob=scales::squish)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2000), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x=NULL,
    # x='Days after fire', 
       y='LAI Anomaly (m²/m²)', 
       color='Recovery', 
       title='Moogem', 
       fill="Precip. Anom. 3-mo (mm)     ")+
  theme_linedraw()+
  theme(
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'none'))
# END **************************************************************************

library(patchwork)
p_ts <- p_mo/p_act/p_ke+plot_annotation(
  tag_levels = list(c("(b)","(c)","(d)"))
   )+plot_layout(
  guides='keep')&theme(
  # legend.position = 'bottom',
  legend.key.width = unit(1.45,'cm'),
  text = element_text(size=18),
  legend.title = element_text(size=15),
  plot.margin = margin(t = 1,r = 10,b = 1,l = 0,unit = 'pt'))
p_ts



ggsave(p_ts, 
       filename = "figures/p_ts.png",
       height=35, 
       width=20.5,
       units='cm',
       dpi = 350)
library(magick)
f_1 <- image_read("figures/p_map.png")
f_1 <- image_annotate(f_1, text="(a)", location = "+10+10",size = 150)
f_2 <- image_read("figures/p_ts.png")
f_out <- image_append(c(f_1,f_2), stack=F)
image_write(f_out, path="figures/figure_combo_map-time-series_moogem-act-kilmoreE_v3.png")

