library(data.table); 
library(tidyverse);
library(lubridate) # load AFTER data.table
library(arrow)
library(stars)
library(patchwork)

# NDVI -------------------------------------------------------------------------
ndvi <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet")

ndvi <- expand_grid(ndvi, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(post_days>ttr5,na.rm=TRUE)/n()) %>% 
  ungroup()

p_ndvi <- ndvi %>% 
  filter(between(hydro_year,2001,2015)) %>% 
  mutate(valid = ifelse(post_days <= 365*(2021.3-hydro_year), T,F)) %>% 
  mutate(val = ifelse(valid==TRUE, val,NA)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
               filter(post_days==min(post_days,na.rm=TRUE)) %>% 
               ungroup(), 
             inherit.aes = FALSE, 
             aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish, 
                       na.value='grey50')+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Fire Year',
       title='Normalized Difference Vegetation Index',
       fill='Fraction Recovered    ')+
  # guides(fill=ggplot2::guide_colorbar(title.position = 'left',
  #                                     title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'right', 
        # legend.key.width = unit(1,'cm'), 
        # legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_ndvi
# END NDVI *********************************************************************



# CCI --------------------------------------------------------------------------
cci <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_cci_ttrDef5_preBS2021-04-19 13:06:11.parquet")
cci <- expand_grid(cci, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(post_days>ttr5_cci,na.rm=TRUE)/n()) %>% 
  ungroup()
p_cci <- cci %>% 
  filter(between(hydro_year,2002,2015)) %>% 
  mutate(valid = ifelse(post_days <= 365*(2021.3-hydro_year), T,F)) %>% 
  mutate(val = ifelse(valid==TRUE, val,NA)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
               filter(post_days==min(post_days,na.rm=TRUE)) %>% 
               ungroup(), 
             inherit.aes = FALSE, 
             aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish, 
                       na.value='grey50')+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2002,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Fire Year',
       title='Chlorphyll Carotenoid Index', 
       fill='Fraction Recovered    ')+
  # guides(fill=ggplot2::guide_colorbar(title.position = 'left',
  #                                     title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'right', 
        # legend.key.width = unit(1,'cm'), 
        # legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_cci


# kNDVI ---------------------------------------------------
kndvi <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_kndvi_ttrDef5_preBS2021-04-21 16:20:09.parquet")

kndvi <- expand_grid(kndvi, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(post_days>ttr5_kn,na.rm=TRUE)/n()) %>% 
  ungroup()

p_kndvi <- kndvi %>% 
  filter(between(hydro_year,2001,2015)) %>% 
  mutate(valid = ifelse(post_days <= 365*(2021.3-hydro_year), T,F)) %>% 
  mutate(val = ifelse(valid==TRUE, val,NA)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
               filter(post_days==min(post_days,na.rm=TRUE)) %>% 
               ungroup(), 
             inherit.aes = FALSE, 
             aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish, 
                       na.value='grey50')+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Fire Year',
       title='Kernel Normalized Difference Vegetation Index',
       fill='Fraction Recovered    ')+
  # guides(fill=ggplot2::guide_colorbar(title.position = 'left',
  #                                     title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'right', 
        # legend.key.width = unit(1.5,'cm'), 
        # legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_kndvi


# LAI ---------------------------------------------------------------------
lai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-05-24 16:20:29.parquet")
lai <- expand_grid(lai, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  mutate(ttr5_lai = ifelse(is.na(ttr5_lai)==T, 5000, ttr5_lai)) %>% 
  mutate(recovered = ifelse(post_days >= ttr5_lai, 1, 0)) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(recovered,na.rm=TRUE)/n()) %>% 
  ungroup()

p_lai <- lai %>% 
  filter(between(hydro_year,2001,2015)) %>% 
  mutate(valid = ifelse(post_days <= 365*(2021.3-hydro_year), T,F)) %>% 
  mutate(val = ifelse(valid==TRUE, val,NA)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
               filter(post_days==min(post_days,na.rm=TRUE)) %>% 
               ungroup(), 
             inherit.aes = FALSE, 
             aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish, 
                       na.value='grey50')+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Fire Year',
       title='Leaf Area Index',
       fill='Fraction Recovered    ')+
  # guides(fill=ggplot2::guide_colorbar(title.position = 'left',
  #                                     title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'right', 
        # legend.key.width = unit(1,'cm'), 
        # legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_lai



# deltaT ------------------------------------------------------------------
deltaT <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_delta_t_ttrDef5_preBS2021-04-25 16:28:57.parquet")
deltaT <- expand_grid(deltaT, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(post_days>ttr5_delta_t,na.rm=TRUE)/n()) %>% 
  ungroup()

p_deltaT <- deltaT %>% 
  filter(between(hydro_year,2001,2015)) %>% 
  mutate(valid = ifelse(post_days <= 365*(2021.3-hydro_year), T,F)) %>% 
  mutate(val = ifelse(valid==TRUE, val,NA)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
               filter(post_days==min(post_days,na.rm=TRUE)) %>% 
               ungroup(), 
             inherit.aes = FALSE, 
             aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish, 
                       na.value='grey50')+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Fire Year',
       title='(Land Skin Temp. - Air Temp.)',
       fill='Fraction Recovered    ')+
  # guides(fill=ggplot2::guide_colorbar(title.position = 'left',
  #                                     title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'right', 
        # legend.key.width = unit(1,'cm'), 
        # legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_deltaT



# tree cover --------------------------------------------------------------
treecover <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_tree_cover_ttrDef5_preBS2021-04-20 16:21:42.parquet")
treecover <- expand_grid(treecover, 
                  post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = fire_year1 ) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum( (post_days/365) >= ttr5_tree_cover,na.rm=TRUE)/n()) %>% 
  ungroup()

p_treecover <- treecover %>% 
  filter(between(hydro_year,2001,2015)) %>% 
  mutate(valid = ifelse(post_days <= 365*(2021.3-hydro_year), T,F)) %>% 
  mutate(val = ifelse(valid==TRUE, val,NA)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
               filter(post_days==min(post_days,na.rm=TRUE)) %>% 
               ungroup(), 
             inherit.aes = FALSE, 
             aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish, 
                       na.value='grey50')+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Fire Year',
       title = "Tree Cover",
       fill='Fraction Recovered    ')+
  # guides(fill=ggplot2::guide_colorbar(title.position = 'left',
  #                                     title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'right', 
        # legend.key.width = unit(1,'cm'), 
        # legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_treecover



# Precip Zonal Anomaly ----------------------------------------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-26 15:23:33.parquet")

# Attach AWAP pixel id to VI ------------------------------------
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


# merges
gc(full=TRUE)
clim <- merge(clim,coords_awap,by=c('x','y'))
gc(full=TRUE)
fits <- merge(fits, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
clim <- clim[idx_awap %in% unique(fits$idx_awap)]
# END Attach climate ***********************************************************

clim2 <- expand_grid(date_fire = seq(ymd("2001-07-01"),length.out = 16, by='1 year'),
                     clim)
dim(clim2)
clim2 <- clim2 %>% 
  mutate(fire_year = year(date_fire)) %>% 
  mutate(post_days = as.double(date-date_fire)) %>% 
  as.data.table()


# PLOT PRECIP ---------------------------------------------------------------------
d1 <- clim2 %>% 
  mutate(val = precip_anom_12mo/map) %>% 
  group_by(fire_year,post_days) %>% 
  summarize(val2 = mean(val,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

d1 %>% dim
d1[fire_year==2002][post_days>=0]

p_precip <- d1 %>% 
  filter(between(fire_year,2001.9,2015.9)) %>% 
  filter(post_days >= 0) %>% 
  # mutate(valid = ifelse(post_days <= 365*(2021.3-fire_year), T,F)) %>%
  # mutate(val2 = ifelse(valid==TRUE, val2,NA)) %>%
  as_tibble() %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+31,
                ymin=fire_year-0.51,
                ymax=fire_year+0.51,
                fill=val2*100))+
  scale_fill_distiller(type='div',palette='RdYlBu',direction = 1,na.value = 'grey')+
  scale_x_continuous(limits=c(365,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2001.25,2015.75), 
                     breaks=seq(2002,by=2,length.out=8))+
  labs(x='Days after July-1 of the Fire Year', 
       y='Fire Year',
       fill = " 12-month Anomaly (%)   ",
       title='Regional Precipitation Anomaly'
       # fill=' 12-month VPD Anomaly (%)   '
  )+
  # guides(fill=ggplot2::guide_colorbar(title.position = 'left',
  #                                     title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(0.8,'cm'),
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey')); p_precip


# VPD Zonal Anomaly ---------------------------------------------------------------------
d2 <- clim2 %>% 
  mutate(val = vpd15_anom_12mo/mavpd15) %>% 
  group_by(fire_year,post_days) %>% 
  summarize(val2 = mean(val,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

p_vpd <- d2 %>% 
  filter(between(fire_year,2001.9,2015.9)) %>% 
  filter(post_days >= 0) %>% 
  # mutate(valid = ifelse(post_days <= 365*(2021.3-fire_year), T,F)) %>%
  # mutate(val2 = ifelse(valid==TRUE, val2,NA)) %>%
  as_tibble() %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+31,
                ymin=fire_year-0.51,
                ymax=fire_year+0.51,
                fill=val2*100))+
  scale_fill_distiller(type='div',palette='PuOr',direction = -1,na.value = 'grey')+
  scale_x_continuous(limits=c(365,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2001.25,2015.75), 
                     breaks=seq(2002,by=2,length.out=8))+
  labs(x='Days after July-1 of the Fire Year', 
       y='Fire Year',
       fill = " 12-month Anomaly (%)   ",
       title='Regional Vapor Pressure Deficit Anomaly'
       # fill=' 12-month VPD Anomaly (%)   '
       )+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(0.8,'cm'),
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey')); p_vpd



# Assemble the plots ------------------------------------------------------
p_top <- p_precip|p_vpd+plot_layout(guides='keep')


p_bot <- (p_ndvi|p_kndvi)/
  (p_lai|p_cci)/
  (p_deltaT|p_treecover)+plot_layout(guides='collect') & 
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1,'cm'),
        legend.key.height = unit(0.2,'cm'))

p_out <- p_top/p_bot+plot_layout(heights=c(1,3.5))+
  plot_annotation(tag_levels = 'a',
    tag_prefix = '(',
    tag_suffix = ')',
    theme=theme(plot.margin = margin(
      t=0,
      r = 0,
      b = 0,
      l=0)))

scale_factor <- 0.333
ggsave(p_out, 
       filename = "figures/plot_tile-fraction_precip-vpd-ndvi-kndvi-lai-cci-deltaT-treecover_labelled.png",
       width=800*scale_factor, 
       height=900*scale_factor, 
       units='mm',
       dpi=350)
