library(tidyverse)
library(data.table)
library(lubridate)
library(stars)
library(patchwork)

# Figures -------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.025) %>% 
  select(NAME_1)

dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")

r_burned <- dat[is.na(fire_doy)==F][,.(x,y,date)] %>% 
  unique %>% 
  .[,`:=`(date_dec = decimal_date(date))] %>% 
  .[,.(val = median(date_dec), 
       nobs = .N), 
    by=.(x,y)]


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
  color='#c316d9',
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
p_act
p_act <- p_act+
  labs(x=NULL,y=NULL)+
  geom_label(
           aes(label='ACT',
           fontface='bold',
           x=148.725+0.025,
           y=-35.586 + 0.025),
           label.size=NA,
           alpha=0.25,
           color='#12e1fc',
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
p_moogem
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
         color='yellow',
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
  ggplot(data=.,aes(x,y,fill=nobs))+
  geom_sf(data=oz_poly, inherit.aes = F, 
          fill='gray90',color='gray30')+
  geom_tile()+
  coord_sf(xlim = c(144.5,153.75),
           ylim = c(-39.1,-27.25), expand = FALSE)+
  scale_fill_viridis_c(option='B',
                       end=0.9,
                       limits=c(1,3),
                       breaks=c(1,2,3),
                       labels=c('1','2','3+'),
                       oob=scales::squish)+
  labs(x=NULL, y=NULL, fill='burns')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill='lightblue'),
        axis.text = element_text(size=17),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        legend.position = c(1,0), 
        legend.justification = c(0.99,0.01)); p_main

p_main <- p_main+geom_rect(aes(xmin=st_bbox(kilmore)$xmin, 
                               xmax=st_bbox(kilmore)$xmax,
                               ymin=st_bbox(kilmore)$ymin,
                               ymax=st_bbox(kilmore)$ymax),
                           fill=NA,lwd=2.1,
                           col='black')+
  geom_rect(aes(xmin=st_bbox(act)$xmin, 
                xmax=st_bbox(act)$xmax,
                ymin=st_bbox(act)$ymin,
                ymax=st_bbox(act)$ymax),
            fill=NA,lwd=2.1,
            col='black')+
  geom_rect(aes(xmin=st_bbox(moogem)$xmin, 
                xmax=st_bbox(moogem)$xmax,
                ymin=st_bbox(moogem)$ymin,
                ymax=st_bbox(moogem)$ymax),
            fill=NA,lwd=2.1,
            col='black')+
  geom_rect(aes(xmin=st_bbox(kilmore)$xmin, 
                     xmax=st_bbox(kilmore)$xmax,
                     ymin=st_bbox(kilmore)$ymin,
                     ymax=st_bbox(kilmore)$ymax),
                 fill=NA,lwd=1,
                 col='#c316d9')+
  geom_rect(aes(xmin=st_bbox(act)$xmin, 
                xmax=st_bbox(act)$xmax,
                ymin=st_bbox(act)$ymin,
                ymax=st_bbox(act)$ymax),
            fill=NA,lwd=1,
            col='#12e1fc')+
  geom_rect(aes(xmin=st_bbox(moogem)$xmin, 
                xmax=st_bbox(moogem)$xmax,
                ymin=st_bbox(moogem)$ymin,
                ymax=st_bbox(moogem)$ymax),
            fill=NA,lwd=1,
            col='yellow')

p_ss <- p_moogem/p_act/p_kilmore
# p_ss
p_out <- p_main+inset_element(p_ss,
                              clip=TRUE,
                              align_to = 'panel',
                              left = 0,
                              bottom = 0.35,
                              right=0.425,
                              top=1 #left bottom right top
                              ) 
# p_out
ggsave(p_out, filename='map-test.png',dpi=350,height=35,width=30,units='cm')
ggsave(p_out, 
        filename = "figures/figure_map-n-burns_insets-moogem-act-kilmore.png", 
       height=35, 
       width=25,
       units='cm',
    dpi = 350)


# SCRAP -------------------------------------------------------------------
# fsize <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
# 
# fsize[rev(order(ba_m2))] %>% 
#   mutate(ba_km2 = ba_m2/(1000**2)) %>% 
#   filter(ba_km2 >= 500000) %>% 
#   filter(date==min(date))
# 
# fsize[rev(order(ba_m2))] %>% 
#   mutate(ba_km2 = ba_m2/(1000**2)) %>% 
#   filter(ba_km2 >= 100000) %>% 
#   filter(date<=ymd("2017-01-01")) %>% 
#   filter(x>150)
# 
# 
# fsize[rev(order(ba_m2))] %>% 
#     mutate(ba_km2 = ba_m2/(1000**2)) %>% 
#     filter(ba_km2 >= 100000) %>% 
#   ggplot(data=.,aes(x,y,fill=decimal_date(date)))+
#   geom_tile()+
#   scale_fill_viridis_c()+
#   coord_equal()
# 
# 
# fsize[near(label,1.452429e+14,tol=1e10)] %>% 
#   ggplot(data=.,aes(x,y,fill=ba_m2))+
#   geom_tile()+
#   scale_fill_viridis_c()+
#   coord_equal()
# 
# 
# dat[date==ymd("2009-03-01")] %>% 
#   ggplot(data=.,aes(x,y,fill=nbr_anom))+
#   geom_tile()+
#   coord_sf(xlim = c(145,153.75),
#            ylim = c(-39.1,-27.5), expand = FALSE)+
#   geom_point(x=145.195,y=-37.37,col='red')
#   
# dat[near(x,145.195,tol=0.01)&near(y,-37.3735,tol=0.01)] %>% 
#   .[date<=ymd("2013-01-01")] %>% 
#   ggplot(data=.,aes(date, nbr_anom,group=id))+
#   geom_line()
# 
# dat[near(x,145.195,tol=0.1)&near(y,-37.3735,tol=0.1)] %>% 
#   .[date==ymd("2009-10-01")] %>% 
#   ggplot(data=.,aes(x,y,fill=nbr_anom))+
#   geom_tile()+
#   coord_equal()+
#   scico::scale_fill_scico(direction = -1)
# 
# 
# 
# plot(kilmore,rgb=c(5,4,3))
# 
# 
# as(kilmore,Class = 'Raster')
# plot(kilmore[,,,1])
# 
# calc_nbr <- function(x) (x[4]-x[7])/(x[4]+x[7]) 
# kilmore_nbr <- st_apply(kilmore, MARGIN = c("x","y"),FUN=calc_nbr) 
# plot(kilmore_nbr,col=viridis::inferno(10,direction = 1),breaks = 'equal')
