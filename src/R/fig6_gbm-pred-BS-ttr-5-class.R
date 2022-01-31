pacman::p_load(tidyverse, data.table, lubridate,stars, sf,magick)

# Load data ------------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 1000) %>% 
  select(NAME_1)

library(h2o) 
h2o.init()
top_mod <- h2o.loadModel("../data_general/proc_data_Oz_fire_recovery/xgboost_ttr5-lai-ocat_2021-06-27_/GBM_5_AutoML_20210627_093118")

dpreds <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_ttr5-lai-5class_gbm_2021-06-27 13:48:15_.parquet")
# out_raster <- raster::rasterFromXYZ(dpreds[,.(x,y,predict)][,`:=`(predict=as.integer(predict))] %>% 
#   rename(z=predict) %>% 
#   as_tibble())
# raster::crs(out_raster) <- raster::crs("EPSG:4326")
# raster::writeRaster(out_raster,filename = "outputs/predicted_ttr5-lai-5class_gbm_2021-06-27.tif",overwrite=T,format="GTiff", 
#   options=c("COMPRESS=ZSTD"))
# raster::raster("outputs/predicted_ttr5-lai-5class_gbm_2021-06-27.tif") %>% 
#   raster::plot()

d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")
d_rf <- d_rf %>% filter(date_fire1 <= lubridate::ymd("2017-03-01")) %>% 
  mutate(fire_year = lubridate::year(date_fire1-months(3))) %>% 
  mutate(fire_year = factor(fire_year))
d_rf <- d_rf %>% mutate(ttr_ocat = factor(ttr_ocat,levels = 1:5, labels = 1:5, ordered = F))
# END LOAD ******************************************************************

vec_labels = c(post_precip_anom_12mo='Post fire 12-mo Precip. anom.',
               vpd15_anom_3mo='Pre fire 3-mo VPD anom.',
               min_nbr_anom='NBR anomaly', 
               elevation="Elevation",
               precip_anom_12mo="Pre fire 12-month Precip. anom.",
               post_vpd15_anom_12mo="Post fire 12-month VPD anom.",
               matmin="Mean Annual Tmin",
               mapet="Mean Annual PET",
               pre_fire_slai_anom_12mo="Pre fire 12-month LAI anom.",
               matmax='Mean Annual Tmax')

p_i2 <- h2o.varimp(top_mod) %>% 
  as_tibble() %>% 
  arrange(-scaled_importance) %>% 
  first(10) %>% 
  mutate(ord = row_number()) %>% 
  mutate(f = fct_reorder(variable,ord)) %>% 
  ggplot(data=.,aes(y=f,
                    x=scaled_importance))+
  geom_col(fill='grey70')+
  scale_y_discrete(limits=rev, 
                   labels=c("post_precip_anom_12mo"='Post fire 12-mo Precip. anom.',
                            "vpd15_anom_3mo"='Pre fire 3-mo VPD anom.',
                            "min_nbr_anom"='NBR anomaly', 
                            "elevation"="Elevation",
                            "precip_anom_12mo"="Pre fire 12-month Precip. anom.",
                            "post_vpd15_anom_12mo"="Post fire 12-month VPD anom.",
                            "matmin"="Mean Annual Tmin",
                            "mapet"="Mean Annual PET",
                       "pre_fire_slai_anom_12mo"="Pre fire 12-month LAI anom.",
                       "matmax"='Mean Annual Tmax'))+
  labs(y=NULL,
       x='Scaled Variable Importance')+
  scale_x_continuous(expand=c(0,0),limits = c(-0.06,1.01))+
  theme_linedraw()+
  theme(text = element_text(size=16),
    panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.length.y = unit(0, 'mm'),
        axis.text.y.left = element_text(hjust=0,
                                        size=14,
                                        colour = 'black',
                                        face='bold',
                            margin = margin(t = 0,
                                            r = -250,
                                            b = 0,
                                            l = 0)), 
        plot.margin = margin(l=30,r=10))
p_i2
ggsave(p_i2,filename = "figures/vip10_ttr-class_BS-pred.png",
       width=12, 
       height=14,
       units='cm',
       dpi=350)






p_i1 <- d_rf %>% 
  # mutate(fire_year = year(date_fire1-months(3))) %>% 
  filter(fire_year %in% c(2002,2006,2008,2013)) %>% 
  mutate(fire_year_f = factor(fire_year)) %>% 
  mutate(ttr_ocat = as.numeric(ttr_ocat)) %>% 
  ggplot(data=.,aes(ttr_ocat,color=fire_year_f))+
  geom_freqpoly(bins=5,aes(y=after_stat(density)))+
  geom_freqpoly(data=dpreds %>% 
                  mutate(ttr_ocat = as.numeric(predict)) %>% 
                  mutate(fire_year_f = 'Pred. 2019'),
                bins=5, aes(ttr_ocat, y=after_stat(density)), 
                lwd=1)+
  scale_x_continuous(breaks=c(1,2,3,4,5),
                     labels=c('≤ 1','2','3','4','≥ 5'),
                     limits=c(1,5),
                     expand=c(0,0.1))+
  scale_color_viridis_d(option='H')+
  labs(x='Time to Recover (year)', 
       color=NULL)+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.99,0.99), 
        legend.justification = c(0.99,0.99), 
        legend.background = element_rect(fill=alpha('white',0)), 
        legend.margin = margin(0,0,0,0));
ggsave(p_i1,filename = "figures/freqpoly_ttr-class_BS-pred.png",
       width=10, 
       height=6.5,
       units='cm',
       dpi=350)

d_metro <- tibble(city=c("Sydney","Canberra","Mallacoota","P. Macquarie"), 
                  x=c(151.21, 149.13, 149.749, 152.9), 
                  y=c(-33.87, -35.28, -37.549, -31.433)) %>% 
  st_as_sf(., coords=c("x","y"),crs=st_crs(4326))

f1 <- dpreds %>% 
  ggplot(data=.,aes(x,y,fill=predict))+
  geom_sf(data=oz_poly, 
          fill='grey70',
          color='grey30',
          inherit.aes = F)+
  geom_raster()+
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
           xlim = c(146,154),
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
  theme_linedraw()+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill='lightblue'))
f1

ggsave(f1, 
       filename = "figures/map_gbm-pred-black-summer-5-classes.png",
       units='cm',
       width=17,
       height=25,
       dpi=350)

pbase <- magick::image_read("figures/map_gbm-pred-black-summer-5-classes.png")
p1 <- magick::image_read("figures/freqpoly_ttr-class_BS-pred.png")
p2 <- magick::image_read( "figures/vip10_ttr-class_BS-pred.png")

pout <- image_composite(pbase,image_scale(p1,"x900"),offset="+127+70")
pout2 <- image_composite(pout, image_scale(p2,"x1200"),offset="+125+970")
# image_write(image = pout2, 
#             path = "test.png")
image_write(image = pout2, 
            path = "figures/map_gbm-pred-black-summer-5-classes-2insets.png")

# SCRATCH -----------------------------------------------------------------


d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")
d_rf <- d_rf %>% filter(date_fire1 <= lubridate::ymd("2017-03-01")) %>% 
  mutate(fire_year = lubridate::year(date_fire1-months(3))) %>% 
  mutate(fire_year = factor(fire_year))
d_rf <- d_rf %>% mutate(ttr_ocat = factor(ttr_ocat,levels = 1:5, labels = 1:5, ordered = F))
d_rf <- d_rf %>% select(ttr_ocat, 
                        fire_year,
                        
                        malai, 
                        fire_month, 
                        min_nbr_anom, 
                        pre_fire_slai_anom_3mo, 
                        pre_fire_slai_anom_12mo, 
                        
                        elevation,
                        slope, 
                        aspect,
                        der,
                        des,
                        sand,silt,clay,
                        pH,
                        bdw,
                        awc,
                        nto,ece,pto,
                        
                        map,mapet,mavpd15,mappet,matmin,matmax,
                        precip_anom_12mo,vpd15_anom_3mo,
                        post_precip_anom_12mo,
                        post_vpd15_anom_12mo)

# drop the NA's, necessary for xgboost
d_rf <- d_rf %>% drop_na()

colnames(d_rf)
response <- "ttr_ocat"
predictors <- setdiff(colnames(d_rf), response)
predictors <- predictors[str_detect(predictors,'fire_year')==F]

d_rf %>% select(predictors)

vip::vi_permute(top_mod,10,geom='boxplot')


h2o.varimp_heatmap(top_mod)

# Mean absolute error
mae <- function(actual, predicted) {
  mean(abs(actual - predicted))
}

vip:::get_training_data(top_mod)
# Permutation-based VIP with user-defined MAE metric
set.seed(1101)  # for reproducibility
out_vip <- vip::vip(top_mod, 
         method = "permute", 
         train=as.h2o(d_rf),
         target = "ttr_ocat", 
         metric = mae,
         nsim=3,
         sample_frac=0.75,
         replace=T,
    smaller_is_better = TRUE,
    pred_wrapper = function(object, newdata) predict(object, newdata)
) + ggtitle("PPR")
