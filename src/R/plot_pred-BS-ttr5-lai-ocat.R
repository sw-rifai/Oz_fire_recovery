library(raster)
library(tidyverse);
library(stars);
library(dtplyr); 
library(data.table); 
library(lubridate); 

dpreds <- arrow::read_parquet("outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet")

# Plot Proportion Recovered by 2021-03 ----------------------------------------------------
bs_start_date <- ymd("2019-09-01")

d0 <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                            col_select = c("x","y","id","malai","lai_yr_sd")) %>% 
  unique()

dlai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
dlai <- dlai[date>=bs_start_date]
dlai <- dlai[order(x,y,date)][
  , lai_12mo := frollapply(lai, FUN=mean,
                                 n = 12,fill = NA,align='right'), by=.(x,y)]

tmp <- merge(d0, dlai[date==max(date)], by=c("x","y"))

dat <- merge(tmp, as.data.table(dpreds %>% select(-malai)), by=c("x","y","id"))
dat[,lai_anom_12mo := lai_12mo-malai]

dat %>% 
  .[,recov := if_else(lai_anom_12mo > -0.25*lai_yr_sd, T,F)] %>% 
  ggplot()+
  geom_bar(aes(y=.pred_class, 
               x=..prop..), 
           stat='count')
  
dat %>% 
  .[,recov := if_else(lai_anom_12mo > -0.25*lai_yr_sd, T,F)] %>% 
  ggplot(data=.,aes(y=.pred_class,
                    # fill=..prop..,
                    fill=factor(recov)))+
  geom_bar(position='fill')+
  scale_y_discrete(limits=rev, 
                   labels=c("≥ 5","4","3","2", "≤ 1"))+
  labs(x='Proportion',
       y="Time to Recover (years)", 
       fill="Recovered\nby 2021-03")+
  scale_fill_viridis_d(end=0.8)+
  facet_wrap(~vc_name)+
  theme_linedraw()
ggsave(filename = "figures/prop-plot_pred-black-summer_recoveredBy202103.png",
       width=20,
       height=10, 
       units='cm',
       dpi=350)
# END SECTION ******************************************************************

# Plot Black Summer time to recover predictions --------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.025) %>% 
  select(NAME_1)


dpreds %>% 
  select(x,y,.pred_class) %>% 
  st_as_stars(., dims=c("x","y")) %>% 
  plot

d_metro <- tibble(city=c("Sydney","Canberra","Mallacoota","Port Macquarie"), 
                  x=c(151.21, 149.13, 149.749, 152.9), 
                  y=c(-33.87, -35.28, -37.549, -31.433)) %>% 
  st_as_sf(., coords=c("x","y"),crs=st_crs(4326))
plot(d_metro)

f1 <- dpreds %>% 
  ggplot(data=.,aes(x,y,fill=.pred_class))+
  geom_sf(data=oz_poly, 
          fill='grey70',
          color='grey30',
          inherit.aes = F)+
  geom_tile()+
  geom_sf_label(data=d_metro, 
                inherit.aes = F, 
                col='black',
                alpha=0.75,
                label.size = NA,
                aes(label=city), 
                # Syd/Can/Malla/PMacq
                nudge_x=c(0,-0.1,0.75,0.5), 
                nudge_y=c(0,0.15,-0.1,-0.1))+
  coord_sf(crs = st_crs(4326),
           xlim = c(146.5,154),
           ylim = c(-38,-28)
  )+
  scale_x_continuous(breaks=seq(148,154,by=2))+
  scale_fill_viridis_d(option='B',
                       begin = 0.1,
                       end=0.9, 
                       breaks=c(1:5), 
                       labels=c("≤1","2","3","4","≥5"))+
  labs(x=NULL, 
       y=NULL, 
       fill='Years')+
  theme_minimal()+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill='lightblue'))
ggsave(f1, 
       filename = "figures/map_pred-black-summer-recovery-classes.png",
       units='cm',
       width=15,
       height=25,
       dpi=350)
# END SECTION ******************************************************************

dpreds <- arrow::read_parquet("outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet")
r <- dpreds %>% 
  mutate(pred_class = as.numeric(.pred_class)) %>% 
  select(x,y,pred_class) %>% 
  st_as_sf(., coords=c('x','y'),crs=st_crs(4326))
  # st_as_stars(., dims=c('x','y'), CRS=st_crs(4326)) %>% 
  # select(pred_class)
# st_crs(r) <- st_crs(4326)
# plot(r)

rotation = function(a){
  r = a * pi / 180 #degrees to radians
  matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
} 
rotation(-90)[1,1:2]

rw


(grd = st_as_stars(st_bbox(r), 
                   nx = 1000, 
                   ny = 1000, 
                   # xlim = c(0, 1.0), ylim = c(0, 1),
                   values = NA_real_))
grd
rw <- st_rasterize(r, template = grd)
attr(attr(rw, "dimensions"), "raster")$affine = c(-0.25,-0.25)# rotation(45)[1,1:2]
plot(rw,col=viridis::inferno(5))
atan(c(-0.25,-0.25))*180/pi

curve(atan2(c(x,1))*180/pi,0,10)


library(grid)
p <- plot(rw)
rotation <- 45
grid.draw(p, vp=viewport(angle=rotation,  width = unit(.75, "npc"), height = unit(.75, "npc")))

g <- ggplot()+geom_stars(data=rw)+coord_sf()+scale_fill_viridis_c()
g <- g + annotation_custom(
  grob = textGrob(label = "", rot = -75,
                  x = unit(1.1, "npc"), y = unit(1.1, "npc")))
print(g, vp = viewport(width = unit(0.5, "npc"),
                       height = unit(0.5, "npc"), angle = -75))

plot(rw,breaks='equal',col=viridis::inferno(5))

plot(as(rw,'Raster'))

rw <- st_warp(r, crs=st_crs(4326), cellsize=0.01,use_gdal = F)

st_get_dimension_values(r,1) %>% diff

as(r,Class = "Raster")
raster::raster(r)

raster::rasterFromXYZ(dpreds %>% select(x,y,.pred_class))



dat %>% 
  .[,recov := if_else(lai_anom_12mo > -0.25*lai_yr_sd, 1,0)] %>% 
  # .[,.(val = sum(recov)/.N), 
  #   by=.(.pred_class, vc_name)] %>% 
  ggplot(., aes(x=.pred_class, fill=factor(recov)))+
  geom_bar()+
  facet_wrap(~vc_name)

dat %>% 
  .[,recov := if_else(lai_anom_12mo > -0.25*lai_yr_sd, 1,0)] %>% 
  .[,.(val = sum(recov)/.N), 
    by=.(.pred_class, vc_name)] %>% 
  ggplot(data=.,aes(y=.pred_class,
                    group=paste(.pred_class,vc_name), 
                    fill=vc_name, 
                    x=val))+
  geom_bar()






dlai3 <- dlai3[,fire_month := lubridate::month(date_fire1)] %>% 
  .[fire_month %in% c(8,9,10,11,12,1,2,3)]
dlai3
# id_bs %>% ggplot(data=.,aes(x,y,fill=nburns))+geom_tile()

dlai4 <- dlai2[date==max(date)]
# END SECTION ******************************************************************

tmp <- merge(as.data.table(dpreds)[,.(x,y,id,.pred_class)], 
      dlai4, 
      by='id')

tmp %>% 
  ggplot(data=.,aes(y=.pred_class, x=slai_anom_12mo/malai, group=.pred_class))+
  geom_boxplot(outlier.colour = NA)

dat %>% 
  .[,recov := if_else(slai_anom_12mo > -0.25*lai_yr_sd, 1,0)] %>% 
  .[,.(val = sum(recov)/.N), 
    by=.(.pred_class)]

