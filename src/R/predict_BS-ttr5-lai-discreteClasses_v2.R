#*******************************************************************************
#* PART 2: FITTING And Predicting THE RF 
# 
library(janitor)
library(dplyr)
library(tidyverse)
library(rsample)
library(recipes)
library(parsnip)
library(tune)
library(dials)
library(workflows)
library(yardstick)
library(themis)

# Load, filter, and split data -------------------------
d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")
d_rf <- d_rf %>% filter(date_fire1 <= lubridate::ymd("2017-03-01"))
d_rf <- d_rf %>% mutate(ttr_ocat = factor(ttr_ocat,levels = 1:5, labels = 1:5, ordered = TRUE))

d_split <- initial_split(d_rf, 
                         strata = ttr_ocat, 
                         breaks=5, 
                         pool=0.1, 
                         prop = 0.1)
d_split2 <- initial_split(d_rf, 
                          strata = ttr_ocat, 
                          breaks=5, 
                          pool=0.1, 
                          prop = 0.5)

d_train <- training(d_split)
d_test <- testing(d_split)

# Specify RF mod specs ------------------
rf_spec <- rand_forest(mode='classification', 
                       mtry=tune(), 
                       trees=tune(),
                       min_n=tune()) %>% 
  set_engine("ranger",
             respect.unordered.factors='order',
             # regularization.factor = tune("regularization"), 
             num.threads=8)


# Specify CV fold data ------------------------
d_folds <- vfold_cv(d_train, v=3)

# Specify 'Recipe' --------------------
rf_rec <- recipe(ttr_ocat ~ min_nbr_anom+
                   fire_month +
                   malai+
                   map+mapet+mavpd15+
                   matmax+matmin+
                   elevation+sand+pH+silt+
                   der+slope+aspect+
                   vpd15_anom_3mo+
                   pre_fire_slai_anom_3mo + 
                   pre_fire_slai_anom_12mo + 
                   post_precip_anom_12mo+
                   post_vpd15_anom_12mo+
                   precip_anom_12mo, 
                 data=d_train) %>% 
  step_smote(ttr_ocat) %>% 
  prep()

# Specify 'Workflow' ------------------
rf_wf <- workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rf_rec)



doParallel::registerDoParallel()
final_rf <- finalize_model(rf_spec, parameters=tibble(mtry=3, trees=998, min_n=3))
final_rf %>% translate()


# Fit Final Model -------------------------
final_wf <- workflow() %>%
  add_model(final_rf) %>% 
  add_recipe(rf_rec)

pdat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/BlackSummer_pred-data_2021-06-14.parquet")
pdat <- pdat %>% rename(pre_fire_slai_anom_3mo = slai_anom_3mo, 
                         pre_fire_slai_anom_12mo = slai_anom_12mo) %>% 
  as.data.table()
pdat <- pdat[is.na(sand)==F & is.na(pH)==F & is.na(silt)==F & is.na(der)==F &
       is.na(slope)==F & is.na(aspect)==F]

dpreds <- final_wf %>%
  # fit the model on all the training data
  fit(
    data    = training(d_split)
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = pdat) %>%
  bind_cols(pdat)

arrow::write_parquet(dpreds,sink=paste0("outputs/pred-black-summer-ttr5-lai-ocat_RF_",Sys.Date(),".parquet"))

# Plot Black Summer time to recover predictions --------------------------------
library(tidyverse); 
library(sf);
library(stars);
library(lubridate)
library(patchwork)
d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")
d_rf <- d_rf %>% mutate(fire_year = year(date_fire1-months(3)))
dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet')
dpreds <- dpreds %>% 
  mutate(fire_month_f = factor(month(date_fire1,abbr = TRUE,label=T), 
                               levels=c("Sep","Oct","Nov","Dec","Jan","Feb"), 
                               ordered = TRUE))
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  select(NAME_1) %>% 
  filter(NAME_1 %in% c("Australian Capital Territory", 
                       "New South Wales",
                       "Queensland",
                       "Victoria")) %>% 
  sf::st_simplify(., dTolerance = 1000)  
st_crs(oz_poly)

d_rf %>% 
  mutate(fire_year = year(date_fire1-months(3))) %>% 
  group_by(fire_year) %>% summarize(nobs=n()) %>% arrange(desc(nobs))

d_rf %>% 
  mutate(fire_year = year(date_fire1-months(3))) %>% 
  mutate(fire_year_f = factor(fire_year)) %>% 
  ggplot(data=.,aes(ttr_ocat,fill=fire_year_f))+
  geom_histogram(binwidth = 1, 
                 aes(y=after_stat(density)))

p_fp <- d_rf %>% 
  mutate(fire_year = year(date_fire1-months(3))) %>% 
  filter(fire_year %in% c(2002,2006,2008,2013)) %>% 
  mutate(fire_year_f = factor(fire_year)) %>% 
  ggplot(data=.,aes(ttr_ocat,color=fire_year_f))+
  geom_freqpoly(bins=5,aes(y=after_stat(density)))

p_i <- p_fp+
  geom_freqpoly(data=dpreds %>% 
                  mutate(ttr_ocat = as.numeric(.pred_class)) %>% 
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
        legend.background = element_rect(fill=NA)); p_i


d_metro <- tibble(city=c("Sydney","Canberra","Mallacoota","P. Macquarie"), 
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
  coord_sf(#crs = st_crs(4326),
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
        panel.background = element_rect(fill='lightblue')); f1
ggsave(f1, 
       filename = "figures/map_pred-black-summer-recovery-classes.png",
       units='cm',
       width=15,
       height=25,
       dpi=350)

p_out <- f1+inset_element(p_i, left = 0,bottom = 0.70,right = 0.6,top=1)
ggsave(p_out, 
       filename = "figures/map_pred-black-summer-recovery-classes-inset.png",
       units='cm',
       width=18,
       height=25,
       dpi=350)


dpreds %>% 
  ggplot(data=.,aes(x,y,fill=.pred_class))+
  geom_tile()+
  coord_sf(xlim = c(147.5,151),
           ylim = c(-38,-35))+
  scale_fill_viridis_d(option='H',end=0.95)

r1 <- st_as_stars(dpreds %>% select(x,y,.pred_class), 
                  dims=c('x','y'))
st_crs(r1) <- st_crs(4326)
plot(r1)
st_crs("+proj=laea")

target_crs <- st_crs("+proj=laea 
             +lat_0=-30 
             +lon_0=140 
             +x_0=432100 
             +y_0=321000 
             +ellps=GRS80 
             +towgs84=0,0,0,0,0,0,0 
             +units=m 
             +no_defs")


s1 <- sf_transform_xy(dpreds %>% select(x,y,.pred_class), 
                target_crs=target_crs, 
                source_crs = st_crs(4326))
s2 <- sf::st_as_sf(s1, coords=c('x','y'))
st_crs(s2) <- target_crs
plot(s2, 
     pch=15, cex=0.15)

ggplot()+
  geom_sf(aes(color=.pred_class), 
          data=s2)



dpreds %>% 
  mutate(pred = as.numeric(.pred_class)) %>% 
  group_by(pred) %>% 
  summarize(across(.fns = mean)) %>% 
  select(pred, 
         malai,
         elevation,
         slope,
         aspect,
         pH,
         min_nbr_anom,
         post_precip_anom_12mo,
         post_vpd15_anom_12mo,
         pre_fire_slai_anom_3mo,
         pre_fire_slai_anom_12mo) %>% 
  ggplot(data=.,aes(pred,post_vpd15_anom_12mo, color=factor(pred)))+
  geom_point(position='jitter')

dpreds %>% 
  filter(between(x,150,152)&
           between(y,-34,-32)) %>% 
  # filter(.pred_class==5) %>% 
  ggplot(data=.,aes(x,y,fill=.pred_class))+
  geom_tile()+
  coord_equal()

dpreds %>% 
  # filter(between(x,150,152)&
  #          between(y,-34,-32)) %>% 
  ggplot(data=., aes(x,y,fill=post_precip_anom_12mo))+
  geom_tile()+
  coord_sf()+
  scale_fill_gradient2(mid='grey70', limits=c(-50,50),oob=scales::squish)

dpreds %>% 
  # filter(between(x,150,152)&
  #          between(y,-34,-32)) %>% 
  ggplot(data=., aes(x,y,fill=post_vpd15_anom_frac))+
  geom_tile()+
  coord_sf()+
  scale_fill_gradient2(mid='grey70', limits=c(-0.1,0.1),oob=scales::squish)


# min_nbr_anom+
#   fire_month +
#   malai+
#   map+mapet+mavpd15+
#   matmax+matmin+
#   elevation+sand+pH+silt+
#   der+slope+aspect+
#   vpd15_anom_3mo+
#   pre_fire_slai_anom_3mo + 
#   pre_fire_slai_anom_12mo + 
#   post_precip_anom_12mo+
#   post_vpd15_anom_12mo+
#   precip_anom_12mo
dpreds %>% 
  # filter(between(x,150,152)&
  #          between(y,-34,-32)) %>% 
  ggplot(data=., aes(x,y,fill=pre_fire_slai_anom_12mo/malai))+
  geom_tile()+
  coord_sf()+
  scale_fill_gradient2(limits=c(-1,1))

dpreds %>% 
  ggplot(data=., aes(x,y,fill=matmax-matmin))+
  geom_tile()+
  coord_sf()+
  scale_fill_viridis_c()

dpreds %>% 
  ggplot(data=., aes(x=matmax-matmin,y=mavpd15))+
  geom_point()+
  geom_smooth()

library(mgcv)
b1 <- bam(ttr5_lai ~ te(min_nbr_anom,
                        mavpd15,
                        post_precip_anom_12mo,
                        by=fire_month_f)
            , 
          data=d_rf %>% filter(fire_y),
          discrete=T)
summary(b1)
plot(b1,scheme = 2,pages = 1)

names(d_rf)

train <- d_rf %>% filter(fire_year != 2013) %>% filter(fire_year != 2002)
b2 <- mgcv::bam(log(ttr5_lai) ~
                  min_nbr_anom*malai,
                  # te(min_nbr_anom,malai, k=5, bs='cs')+
                  # te(pre_fire_slai_anom_3mo, pre_fire_slai_anom_12mo,k=5)+
                  # te(map,post_precip_anom_12mo, k=5,bs='cs')+
                  # fire_month_f,
                  # te(log(elevation+10), 
                  #    log(der),k=5)+
                  # s(pH,k=5)+
                  # s( I(matmax-matmin), k=5)+
                  # te(min_nbr_anom,malai, k=5)+
                  # s(mavpd15,k=5)+
                  # te(map, precip_anom_12mo, post_precip_anom_12mo,k=4),
                data= train, 
                family=Gamma(link='log'),
                # family=Tweedie(p=1.67, 
                #                link = 'log'),
                select=TRUE, 
                discrete=TRUE)
summary(b2)
plot(b2, scheme = 2, pages = 1)

test <- d_rf %>% 
  filter(is.na(ttr5_lai)==F) %>% 
  filter(fire_year==2013 | fire_year==2002) %>% 
  mutate(pred = exp(predict(b2,newdata=., newdata.guaranteed = TRUE, type='response'))) 
table(is.na(test$ttr5_lai))
table(is.na(test$pred))
yardstick::rsq_trad_vec(truth = test$ttr5_lai, 
                        estimate = test$pred, na_rm = T)
yardstick::rmse_vec(truth = test$ttr5_lai, 
                        estimate = test$pred)
yardstick::rmse_vec(truth = test$ttr5_lai, 
                    estimate = rep(mean(train$ttr5_lai,na.rm=T),dim(test)[1]))
test$pred %>% hist
test$ttr5_lai %>% hist
test$pred %>% summary
test$ttr5_lai %>% summary

test %>% ggplot(data=.,aes(pred, ttr5_lai))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')
test %>% select(pred,ttr5_lai) %>% drop_na() %>% cor %>% `^`(2)

dpreds <- dpreds %>% 
  mutate(pred_ttr5_lai = exp(predict(b2, newdata=., type='response')))
dpreds$pred_ttr5_lai %>% hist(breaks=100)

dpreds %>% 
  mutate(b2_p = cut(pred_ttr5_lai, 
                    c(0,367,365*2,365*3,365*4,Inf), 
                    labels=c(1,2,3,4,5))) %>% 
  mutate(b2_p = factor(b2_p,ordered = F)) %>% 
  ggplot(data=.,aes(x,y,fill=b2_p))+
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
  coord_sf(#crs = st_crs(4326),
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
