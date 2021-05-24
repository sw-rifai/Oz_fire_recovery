library(tidyverse); 
library(dtplyr); 
library(data.table); 
library(lubridate); 

clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-26 15:23:33.parquet")

# Attach AWAP pixel id to VI ------------------------------------
coords_vi <- lazy_dt(fits) %>% select(x,y,id) %>% distinct() %>% as.data.table()
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

clim2 <- clim2 %>% 
  lazy_dt() %>% 
  mutate(fire_year = year(date_fire)) %>% 
  mutate(post_days = as.double(date-date_fire)) %>% 
  # filter(post_days >= 0) %>%
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

d1 %>% 
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
       y='Bushfire Year',
       fill=' 12-month Precip. Anomaly (%)   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey'))

ggsave(filename = 'figures/figure_tile_mean-frac-precip-12mo-anom.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)

# PLOT VPD ---------------------------------------------------------------------
d2 <- clim2 %>% 
  mutate(val = vpd15_anom_12mo/mavpd15) %>% 
  group_by(fire_year,post_days) %>% 
  summarize(val2 = mean(val,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

d2 %>% 
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
       y='Year of Bushfire',
       fill=' 12-month VPD Anomaly (%)   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey'))

ggsave(filename = 'figures/figure_tile_mean-frac-vpd15-12mo-anom.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)




d3 <- clim2 %>% 
  mutate(val = ppet_anom_12mo/mappet) %>% 
  group_by(fire_year,post_days) %>% 
  summarize(val2 = mean(val,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

d3 %>% 
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
  # scale_fill_distiller(type='div',palette='PiYG',direction = 1,na.value = 'grey')+
  scico::scale_fill_scico(palette='broc',direction=-1, limits=c(-50,50), oob=scales::squish)+
  scale_x_continuous(limits=c(365,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2001.25,2015.75), 
                     breaks=seq(2002,by=2,length.out=8))+
  labs(x='Days after July-1 of the Fire Year', 
       y='Year of Bushfire',
       fill=' 12-month P:PET Anomaly (%)   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey'))

ggsave(filename = 'figures/figure_tile_mean-frac-ppet-12mo-anom.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)
