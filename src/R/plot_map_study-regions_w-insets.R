library(tidyverse)
library(data.table)
library(lubridate)
library(stars)
library(patchwork)

# Figures -------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet", 
                            col_select = c("x","y","id","date","sndvi","ndvi_anom","fire_doy",
                                           "nbr_anom",
                                           "ndvi_u","ndvi_sd","ndvi_anom_sd"))
gc(full=T)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet",
                            col_select = c("x","y","id","date","sndvi","ndvi_anom","fire_doy",
                                           "nbr_anom",
                                           "ndvi_u","ndvi_sd","ndvi_anom_sd"))
gc(full=T)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE); gc(full=TRUE)
rm(tmp1,tmp2)
gc(full=T)

fsize <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")


r_burned <- dat[is.na(fire_doy)==F][,.(x,y,date)] %>% 
  unique %>% 
  .[,`:=`(date_dec = decimal_date(date))] %>% 
  .[,.(val = median(date_dec), 
       nobs = .N), 
    by=.(x,y)]


# Prep Kilmore East fire --------------------------------------------------
kilmore <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/kilmore_east_LT05_L1TP_092086_20090216_20161029_01_T1.tif")
p_kilmore <- RStoolbox::ggRGB(as(kilmore,Class = 'Raster'),r = 5, g=4, b=3, stretch='hist')
p_kilmore <- p_kilmore+
  labs(x=NULL,y=NULL)+
  annotate(geom='text',
           label='Kilmore E.',
           fontface='bold',
           color='black',
           size=9.7,
           x=145.07,
           y=-37.485)+
  # annotate(geom='text',
  #          label='Kilmore East',
  #          fontface='bold',
  #          color='yellow',
  #          size=9.5,
  #          x=145.06,
  #          y=-37.485)+
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
p_act <- p_act+
  labs(x=NULL,y=NULL)+
  annotate(geom='text',
           label='ACT',
           fontface='bold',
           color='yellow',
           size=10,
           x=148.725,
           y=-35.6)+
  theme_linedraw()+
  scale_x_continuous(expand=c(0,0),
                     # limits=c(145,145.45)
                     )+
  scale_y_continuous(expand=c(0,0))+
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
  annotate(geom='text',
           label='Moogem',
           fontface='bold',
           color='yellow',
           size=10,
           x=152.2,
           y=-29.57)+
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
  coord_sf(xlim = c(145,153.75),
           ylim = c(-39.1,-27.5), expand = FALSE)+
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
        legend.position = c(1,0), 
        legend.justification = c(1,0)); p_main

p_main <- p_main+geom_rect(aes(xmin=st_bbox(kilmore)$xmin, 
                               xmax=st_bbox(kilmore)$xmax,
                               ymin=st_bbox(kilmore)$ymin,
                               ymax=st_bbox(kilmore)$ymax),
                           fill=NA,lwd=2,
                           col='lightblue')+
  geom_rect(aes(xmin=st_bbox(act)$xmin, 
                xmax=st_bbox(act)$xmax,
                ymin=st_bbox(act)$ymin,
                ymax=st_bbox(act)$ymax),
            fill=NA,lwd=2,
            col='lightblue')+
  geom_rect(aes(xmin=st_bbox(moogem)$xmin, 
                xmax=st_bbox(moogem)$xmax,
                ymin=st_bbox(moogem)$ymin,
                ymax=st_bbox(moogem)$ymax),
            fill=NA,lwd=2,
            col='lightblue')+
  geom_rect(aes(xmin=st_bbox(kilmore)$xmin, 
                     xmax=st_bbox(kilmore)$xmax,
                     ymin=st_bbox(kilmore)$ymin,
                     ymax=st_bbox(kilmore)$ymax),
                 fill=NA,
                 col='blue')+
  geom_rect(aes(xmin=st_bbox(act)$xmin, 
                xmax=st_bbox(act)$xmax,
                ymin=st_bbox(act)$ymin,
                ymax=st_bbox(act)$ymax),
            fill=NA,
            col='blue')+
  geom_rect(aes(xmin=st_bbox(moogem)$xmin, 
                xmax=st_bbox(moogem)$xmax,
                ymin=st_bbox(moogem)$ymin,
                ymax=st_bbox(moogem)$ymax),
            fill=NA,
            col='blue')

p_main

p_ss <- p_moogem/p_act/p_kilmore

p_out <- p_main+inset_element(p_ss,0,0.5,0.5,1) #left bottom right top
p_out

ggsave(p_out, 
        filename = "figures/figure_map-n-burns_insets-moogem-act-kilmore.png", 
       height=20, 
       width=15,
       units='cm',
    dpi = 350)


# SCRAP -------------------------------------------------------------------
fsize[rev(order(ba_m2))] %>% 
  mutate(ba_km2 = ba_m2/(1000**2)) %>% 
  filter(ba_km2 >= 500000) %>% 
  filter(date==min(date))

fsize[rev(order(ba_m2))] %>% 
  mutate(ba_km2 = ba_m2/(1000**2)) %>% 
  filter(ba_km2 >= 100000) %>% 
  filter(date<=ymd("2017-01-01")) %>% 
  filter(x>150)


fsize[rev(order(ba_m2))] %>% 
    mutate(ba_km2 = ba_m2/(1000**2)) %>% 
    filter(ba_km2 >= 100000) %>% 
  ggplot(data=.,aes(x,y,fill=decimal_date(date)))+
  geom_tile()+
  scale_fill_viridis_c()+
  coord_equal()


fsize[near(label,1.452429e+14,tol=1e10)] %>% 
  ggplot(data=.,aes(x,y,fill=ba_m2))+
  geom_tile()+
  scale_fill_viridis_c()+
  coord_equal()


dat[date==ymd("2009-03-01")] %>% 
  ggplot(data=.,aes(x,y,fill=nbr_anom))+
  geom_tile()+
  coord_sf(xlim = c(145,153.75),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  geom_point(x=145.195,y=-37.37,col='red')
  
dat[near(x,145.195,tol=0.01)&near(y,-37.3735,tol=0.01)] %>% 
  .[date<=ymd("2013-01-01")] %>% 
  ggplot(data=.,aes(date, nbr_anom,group=id))+
  geom_line()

dat[near(x,145.195,tol=0.1)&near(y,-37.3735,tol=0.1)] %>% 
  .[date==ymd("2009-10-01")] %>% 
  ggplot(data=.,aes(x,y,fill=nbr_anom))+
  geom_tile()+
  coord_equal()+
  scico::scale_fill_scico(direction = -1)



plot(kilmore,rgb=c(5,4,3))


as(kilmore,Class = 'Raster')
plot(kilmore[,,,1])

calc_nbr <- function(x) (x[4]-x[7])/(x[4]+x[7]) 
kilmore_nbr <- st_apply(kilmore, MARGIN = c("x","y"),FUN=calc_nbr) 
plot(kilmore_nbr,col=viridis::inferno(10,direction = 1),breaks = 'equal')
