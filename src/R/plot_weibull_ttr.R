library(phenofit);
library(data.table); 
library(dtplyr, warn.conflicts = FALSE); 
library(tidyverse);
library(usethis);
library(stars);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
source("src/R/functions_time_to_recover.R")

# Isolate slow recovering pixels --------------------------------
load("outputs/pixel_vegClass_groups.rds")
sdat <- arrow::read_parquet("outputs/weibull_fits_pre2005_fires_2021-01-18.parquet")
tmp2 <- expand_grid(merge(sdat,nvis, by='id') %>% 
                      filter(vc!=25) %>%
                      filter(vc %in% c(2,3,4,5,11)) %>% 
                      filter(is.na(vc)==FALSE) %>% 
                      sample_n(100) %>% 
                      as.data.table(), 
                    pred_days=seq(1,2000,length.out=2000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table()

vec_slowgrow <- tmp2[near(pred_days,100)][p_diff<0.0001]$id
vec_slowgrow <- tmp2 %>% lazy_dt() %>% 
  filter(pred_days <= 100) %>% 
  group_by(id) %>%
  summarize(low_grow_days = sum(p_diff <= 0.00005)) %>% 
  ungroup() %>% 
  as.data.table()

vec_inflection <- tmp2[pred_days==1] %>% 
  mutate(inflection = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))^(1.0/pwr)) %>% 
  select(id,inflection)

tmp2 %>% 
  merge(., vec_inflection, by='id') %>% 
  mutate(slowgrow = if_else(id%in%vec_slowgrow[low_grow_days>=100]$id,
                            'delayed','instant')) %>% 
  ggplot(data=.,aes(pred_days, pred, group=factor(id)))+
  geom_vline(aes(xintercept=inflection,group=factor(id),color=inflection,alpha=inflection),lwd=0.25,lty=1)+
  geom_hline(aes(yintercept=-0.5*Drop,group=factor(id),color=inflection,alpha=inflection),lwd=0.25,lty=1)+
  geom_line(lwd=0.1)+
  # geom_vline(aes(xintercept=100),col='red')+
  scale_x_continuous(expand=c(0,0))+
  labs(x='Days post fire',
       y='NDVI Anom.',
       title='Weibull Fit - Fires 2003/2004')+
  scale_color_viridis_c('Linear TTR (days)',end=0.99,option='B',direction=-1)+
  facet_wrap(~vc_name)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        # legend.position = c(1,0), 
        # legend.justification = c(1,0)
  )
ggsave(filename = 'figures/ObsTTR_delayed_inflection_WeibullTTR_2003-04_fires.png', 
       width=20, height=12, units='cm')


# TO DO: Use inflection (~ -0.5*Drop), number of days with ~0 growth, and TTR to 
#        construct classes of vegetation recovery
# Maybe with a PCA, or a k-means
# END ****************************************************************


# Largscale 1burn recovery shapes ------------------------------------
w1 <- read_parquet(file = 'outputs/weibull_fits_1burn_2001-2014fires_2021-01-27.parquet')

d_delayed <- expand_grid(merge(w1,nvis, by='id') %>% 
              filter(vc!=25) %>%
              filter(vc %in% c(2,3,4,5,11)) %>% 
              filter(is.na(vc)==FALSE), 
            pred_days=seq(1,1000,length.out=1000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table() %>% 
  lazy_dt() %>% 
  filter(pred_days <= 100) %>% 
  group_by(id) %>%
  summarize(low_grow_days = sum(p_diff <= 0.00005)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)

d_inflection <- w1 %>% lazy_dt() %>% 
  mutate(inflection = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))^(1.0/pwr)) %>% 
  select(id,inflection) %>% 
  filter(inflection < 2000) %>% 
  as.data.table()

merge(w1,nvis, by='id') %>%
  lazy_dt() %>% 
  filter(vc!=25) %>%
  filter(vc %in% c(2,3,5,11)) %>% 
  filter(is.na(vc)==FALSE) %>% 
  filter(id %in% sample(d_inflection[inflection>100]$id,1000)) %>% 
  as.data.table() %>% 
  expand_grid(., pred_days=seq(1,1000,length.out=20) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(pred_days, pred, group=factor(id)))+
  geom_line(lwd=0.1)+
  geom_smooth(aes(pred_days,pred,color=vc_name),inherit.aes = FALSE)+
  facet_wrap(~vc_name)

p_1 <- merge(w1,nvis, by='id') %>% merge(., d_inflection,by='id') %>%
  ggplot(data=.,aes(date_first_fire, inflection))+
  ggpointdensity::geom_pointdensity(adjust=5,size=0.5)+
  geom_smooth(col='#CF0000')+
  labs(x='Date of Fire', y='Days until inflection point')+
  scale_color_viridis_c(option='A')+
  theme_linedraw()
ggsave(p_1, filename='figures/timeseries_daysUntilInflection_1burn_2001-2014.png', 
       width=14,height=7,units='cm')

nvis %>% select(vc,vc_name) %>% distinct()
colorspace::hcl_palettes(plot = T)
library(mgcv)
(p_2 <- merge(w1,nvis, by='id') %>% merge(., d_inflection,by='id') %>%
  filter(vc %in% c(2,3,5,11)) %>% 
  mutate(inflection = if_else(inflection<1,1,inflection)) %>% 
  ggplot(data=.,aes(date_first_fire, inflection,color=vc_name))+
  geom_point(alpha=0.1)+
  geom_smooth(method='loess',se=FALSE)+
  # geom_smooth(method='bam', 
  #             formula=y~s(x,bs='cs'),
  #             method.args=list(method='fREML', family=Gamma(link='log')))+
  labs(x='Date of Fire', y='Days until inflection point')+
  scico::scale_color_scico_d(palette = 'batlow')+
  theme_linedraw())

(p_2 <- merge(w1,nvis, by='id') %>% merge(., d_inflection,by='id') %>%
    filter(vc %in% c(2,3,5,11)) %>% 
    mutate(inflection = if_else(inflection<1,1,inflection)) %>% 
    ggplot(data=.,aes(date_first_fire, inflection,color=vc_name))+
    stat_summary()+
        # geom_point(alpha=0.1)+
    # geom_smooth(method='bam', 
    #             formula=y~s(x,bs='cs'),
    #             method.args=list(method='fREML', family=Gamma(link='log')))+
    labs(x='Date of Fire', y='Days until inflection point')+
    scico::scale_color_scico_d('NVIS class', 
                               palette = 'batlow',end = 0.8)+
    theme_linedraw()+
    theme(legend.position = 'bottom'))
ggsave(p_2, filename='figures/timeseries_pointRange_daysUntilInflection_1burn_2001-2014.png', 
       width=20,height=10,units='cm')



# Connect time to inflection point with climate --------------------------
load("outputs/pixel_vegClass_groups.rds")
w1 <- read_parquet(file = 'outputs/weibull_fits_1burn_2001-2014fires_2021-01-27.parquet')
d_delayed <- expand_grid(merge(w1,nvis, by='id') %>% 
                           filter(vc!=25) %>%
                           filter(vc %in% c(2,3,4,5,11)) %>% 
                           filter(is.na(vc)==FALSE), 
                         pred_days=seq(1,1000,length.out=1000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table() %>% 
  lazy_dt() %>% 
  filter(pred_days <= 100) %>% 
  group_by(id) %>%
  summarize(low_grow_days = sum(p_diff <= 0.00005)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)

d_inflection <- w1 %>% lazy_dt() %>% 
  mutate(inflection = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))^(1.0/pwr)) %>% 
  select(id,inflection) %>% 
  filter(inflection < 2000) %>% 
  as.data.table()

w1 <- merge(w1,nvis, by='id') %>% merge(., d_inflection,by='id') %>% 
  filter(vc %in% c(2,3,5,11)) 
  
coords_vi <- st_as_sf(w1, coords = c("x","y"))

clim <- arrow::read_parquet(file = "/home/sami/scratch/awap_clim_se_coastal.parquet")
# col_select=c("x","y","date","precip_anom_12mo"))
coords_awap <- st_as_sf(unique(clim[,.(x,y)]), coords=c("x","y"))
nn_coords <- RANN::nn2(
  coords_awap %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)

# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] #%>% as.data.table()
coords_keep_awap$x <- st_coordinates(coords_keep_awap)[,'X']
coords_keep_awap$y <- st_coordinates(coords_keep_awap)[,'Y']
coords_keep_awap$id_vi <- coords_vi$id

gc(full=TRUE)
# coords_keep_awap %>% select(x,y,idx_awap) %>% distinct()
coords_keep_awap <- coords_keep_awap %>% select(x,y,idx_awap,id_vi) %>% as.data.table() %>% 
  select(-geometry) %>% distinct() %>% 
  rename(id=id_vi)

clim <- merge(clim,
              coords_keep_awap %>% select(x,y,idx_awap) %>% distinct() %>% as.data.table(),
              by=c("x","y"), 
              all.y=TRUE)

dat <- merge(w1, coords_keep_awap %>% select(id,idx_awap), by='id')
gc(full=TRUE)
dat <- dat %>% rename(date = date_first_fire)
dat <- merge(dat, clim, by=c("idx_awap","date"))
# dat <- merge(dat, clim[,.(idx_awap,date,
#                           precip_anom_12mo,pet_anom_12mo)], by=c("idx_awap","date"))

dat %>% 
  ggplot(data=.,aes(mappet, inflection,color=vc_name))+
  geom_point(alpha=0.1)+
  geom_smooth(method='lm')

lm(ttr~ppet_anom_12mo+post_ppet_anom_12mo, data=dat) %>% summary


# mgcv time to inflection ---------------------------------------
library(mgcv)
fg0 <- bam(inflection~s(delta_vi_12mo,bs='cs',k=5),
           data=dat, select=TRUE)

fg1 <- bam(inflection~s(precip_anom_12mo,bs='cs',k=5)+
             s(post_precip_anom_12mo,bs='cs',k=5)+
             s(vpd15_anom_12mo,bs='cs',k=5)+
             s(post_vpd15_anom_12mo,bs='cs',k=5), 
           data=dat, select=TRUE)
summary(fg1)
plot(fg1)
gratia::draw(fg1)





library(tidymodels)
set.seed(123)
xgdat <- dat %>% select(inflection, 
               starts_with('post'), 
               contains('anom'),
               contains("ma"),
               contains("_u"),
               vc_name, 
               delta_vi_12mo) %>% 
  select(-post_fire_vi_36mo, -post_fire_vi_12mo) %>% 
  sample_n(10000)

vb_split <- initial_split(xgdat, strata = inflection)
vb_train <- training(vb_split)
vb_test <- testing(vb_split)

xgb_spec <- boost_tree(
  trees = 1000, 
  tree_depth = tune(), min_n = tune(), 
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), mtry = tune(),         ## randomness
  learn_rate = tune(),                         ## step size
) %>% 
  set_engine("xgboost") %>% 
  set_mode("regression")

xgb_spec

xgb_grid <- grid_latin_hypercube(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), vb_train),
  learn_rate(),
  size = 10
)

xgb_grid

xgb_wf <- workflow() %>%
  add_formula(inflection ~ .) %>%
  add_model(xgb_spec)

xgb_wf

set.seed(123)
vb_folds <- vfold_cv(vb_train, strata = inflection)

vb_folds

doParallel::registerDoParallel()

set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = vb_folds,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)

xgb_res


collect_metrics(xgb_res)

xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")


show_best(xgb_res, "rmse")


best_rmse <- select_best(xgb_res, "rmse")
best_rmse

final_xgb <- finalize_workflow(
  xgb_wf,
  best_rmse
)

final_xgb


library(vip)

final_xgb %>%
  fit(data = vb_train) %>%
  pull_workflow_fit() %>%
  vip(geom = "point")







tmp2 %>% 
  merge(., vec_inflection, by='id') %>% 
  mutate(slowgrow = if_else(id%in%vec_slowgrow[low_grow_days>=100]$id,
                            'delayed','instant')) %>% 
  ggplot(data=.,aes(pred_days, pred, group=factor(id)))+
  geom_vline(aes(xintercept=inflection,group=factor(id),color=inflection,alpha=inflection),lwd=0.25,lty=1)+
  geom_hline(aes(yintercept=-0.5*Drop,group=factor(id),color=inflection,alpha=inflection),lwd=0.25,lty=1)+
  geom_line(lwd=0.1)+
  # geom_vline(aes(xintercept=100),col='red')+
  scale_x_continuous(expand=c(0,0))+
  labs(x='Days post fire',
       y='NDVI Anom.',
       title='Weibull Fit - Fires 2003/2004')+
  scale_color_viridis_c('Linear TTR (days)',end=0.99,option='B',direction=-1)+
  facet_wrap(~vc_name)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        # legend.position = c(1,0), 
        # legend.justification = c(1,0)
  )




sdat <- arrow::read_parquet("outputs/linear_ttr_2003-2004_2021-01-19.parquet")
ss <- sdat[date_first_fire < ymd('2005-01-01')][is.na(ttr)==F][fire_count==1]
ssdat <- dat[id%in%ss$id]

ssdat <- merge(ssdat, 
               ss[,.(x,y,id,date_first_fire,recovery_date,ttr)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>% 
  group_by(x,y,id) %>% 
  filter(date > date_first_fire) %>% 
  filter(date <= recovery_date) %>% 
  ungroup() %>% 
  as.data.table() 

mdat <- mdat %>% lazy_dt() %>% 
  mutate(post_days = as.double(date - date_first_fire)) %>% 
  as.data.table()

n_w <- arrow::read_parquet("outputs/weibull_fits_pre2005_fires_2021-01-18.parquet")


expand_grid(merge(n_w[isConv==TRUE][sample(.N,100)], mdat, by=c("x","y","id")), 
            pred_days=seq(1,2000,length.out=100)) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
  # mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  # mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
  # pull(delta_growth) %>% plot
  # mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
  # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
  # filter(pred_ttr < 2500) %>%
  ggplot(data=.,aes(pred_days, pred, group=factor(id), color=ttr))+
  geom_line(lwd=0.1)+
  scale_x_continuous(expand=c(0,0))+
  labs(x='Days post fire',
       y='NDVI Anom.',
       title='Weibull Fit - Fires 2003/2004')+
  scale_color_viridis_c('Linear TTR (days)')+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))




expand_grid(merge(n_w[isConv==TRUE][sample(.N,100)], mdat, by=c("x","y","id")), 
            pred_days=seq(1,2000,length.out=100)) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
  # mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  # mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
  # pull(delta_growth) %>% plot
  # mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
  # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
  # filter(pred_ttr < 2500) %>%
  ggplot(data=.,aes(pred_days, pred, group=factor(id), color=ttr))+
  geom_line(lwd=0.1)+
  scale_x_continuous(expand=c(0,0))+
  labs(x='Days post fire',
       y='NDVI Anom.',
       title='Weibull Fit - Fires 2003/2004')+
  scale_color_viridis_c('Linear TTR (days)')+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))
ggsave(filename = 'figures/timeseries100_ndviAnom_WeibullTTR_2003-04_fires.png')




expand_grid(merge(n_w[isConv==TRUE][Drop<0.6][sample(.N,100)], mdat, by=c("x","y","id")), 
            pred_days=seq(1,2000,length.out=100)) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
  # mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  # mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
  # pull(delta_growth) %>% plot
  # mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
  # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
  # filter(pred_ttr < 2500) %>%
  ggplot(data=.,aes(pred_days, pred, group=factor(id)))+
  geom_line(lwd=0.1)
# ggplot(data=.,aes(pred_days, pred, group=as_factor(id)))+
# geom_hline(aes(yintercept=0),col='grey',lwd=1.5)+
# geom_line(lwd=1,col='red')+
# geom_line(aes(pred_days, delta_growth),col='blue')+
# geom_point(aes(post_days, ndvi_anom),alpha=0.025)+
# geom_vline(aes(xintercept=pred_ttr),lty=3)+
# scale_color_viridis_c(end=0.9)+
# labs(x='Days post fire', 
#      y='NDVI Anom.', 
#      title='Weibull Fit - Fires 2003/2004')+
# facet_wrap(~id,ncol=4)+
# theme_linedraw()





expand_grid(merge(n_w[isConv==TRUE][Drop<0.6][sample(.N,1)], mdat, by=c("x","y","id")), 
            pred_days=seq(1,2000,length.out=100)) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
  mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
  pull(delta_growth) %>% plot
# mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
# mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
filter(pred_ttr < 2500) %>%
  ggplot(data=.,aes(pred_days, delta_growth, color=factor(id)))+
  geom_line()
# ggplot(data=.,aes(pred_days, pred, group=as_factor(id)))+
# geom_hline(aes(yintercept=0),col='grey',lwd=1.5)+
# geom_line(lwd=1,col='red')+
# geom_line(aes(pred_days, delta_growth),col='blue')+
# geom_point(aes(post_days, ndvi_anom),alpha=0.025)+
# geom_vline(aes(xintercept=pred_ttr),lty=3)+
# scale_color_viridis_c(end=0.9)+
# labs(x='Days post fire', 
#      y='NDVI Anom.', 
#      title='Weibull Fit - Fires 2003/2004')+
# facet_wrap(~id,ncol=4)+
# theme_linedraw()

expand_grid(n_w[isConv==TRUE][Drop<0.6][sample(.N,10)], post_days=1:1000) %>% 
  mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
  ggplot(data=.,aes(post_days, delta_growth,group=id))+
  geom_line()



expand_grid(merge(n_w[isConv==TRUE][Drop<0.6][sample(.N,14)], mdat, by=c("x","y","id")), 
            pred_days=seq(1,2000,length.out=20)) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
  mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>%
  # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
  filter(pred_ttr < 2500) %>%
  ggplot(data=.,aes(pred_days, pred, group=as_factor(id)))+
  geom_hline(aes(yintercept=0),col='grey',lwd=1.5)+
  geom_line(lwd=1,col='red')+
  geom_point(aes(post_days, ndvi_anom),alpha=0.025)+
  geom_vline(aes(xintercept=pred_ttr),lty=3)+
  scale_color_viridis_c(end=0.9)+
  labs(x='Days post fire', 
       y='NDVI Anom.', 
       title='Weibull Fit - Fires 2003/2004')+
  facet_wrap(~id,ncol=4)+
  theme_linedraw()
ggsave(filename = 'figures/timeseries_ndviAnom_WeibullTTR_2003-04_fires.png')

junk <- n_w[sample(.N,10)]
junk
expand_grid(junk, post_days=1:1500) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
  ggplot(data=.,aes(post_days, pred, group=id,color=pred_ttr50/Drop))+
  geom_line()+
  scale_color_viridis_c()

expand_grid(merge(n_w[sample(.N,10)], mdat, by=c("x","y","id")), 
            pred_days=seq(1,1500,length.out=10))



merge(sdat, n_w, by=c("x","y","id")) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  ggplot(data=.,aes(factor(month_of_fire), delta_vi_12mo))+
  geom_boxplot()+
  labs(x='Month of Fire', y=expression(Annual~NDVI[post-fire]~-~NDVI[pre-fire]), 
       title='Fires between 2003-2004')+
  theme_linedraw()
ggsave(filename = 'figures/monthOfFire_deltaNDVI_2003-04_fires.png')


merge(sdat, n_w, by=c("x","y","id")) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  ggplot(data=.,aes(factor(month_of_fire), ttr))+
  geom_boxplot()+
  labs(x='Month of Fire', y='Time to Recover (days)', 
       title='Fires between 2003-2004')+
  theme_linedraw()
ggsave(filename = 'figures/monthOfFire_TTR_2003-04_fires.png')

merge(sdat, n_w, by=c("x","y","id")) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  ggplot(data=.,aes(delta_vi_12mo, ttr,color=month_of_fire))+
  geom_point(alpha=0.05,size=1)+
  geom_smooth(color='#fc0000', se=F, 
              data=. %>% filter(between(delta_vi_12mo,-0.4,0.2)))+
  scale_color_viridis_c('Month of Fire', option='B', end=0.9)+
  labs(x=expression(Annual~NDVI[post-fire]~-~NDVI[pre-fire]), 
       y='Time to Recover (days)', 
       subtitle='Fires between 2003-2004')+
  theme_linedraw()
ggsave(filename = 'figures/fireDeltaNDVI12mo_TTR_2003-04_fires.png')


merge(sdat, n_w, by=c("x","y","id")) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
  mutate(pred_ttr = w_ttr(Asym, Drop, lrc, pwr)) %>% 
  filter(isConv==TRUE) %>%
  # filter(pred_ttr50 < 2000) %>%  #pull(pred_ttr50) %>% hist
  filter(pred_ttr < 10000) %>%  #pull(pred_ttr50) %>% hist
  filter(delta_vi_12mo < 0) %>% 
  ggplot(data=.,aes(pred_ttr, ttr,color=month_of_fire))+
  geom_hline(aes(yintercept=365*9),col='grey',lwd=2)+
  geom_point(alpha=0.05,size=1)+
  geom_smooth(color='#fc0000', se=F, method='lm' 
              # data=. %>% filter(between(delta_vi_12mo,-0.4,0.2))
  )+
  geom_abline(color='blue',lwd=1)+
  scale_color_viridis_c('Month of Fire', option='B', end=0.9)+
  labs(x='Time to Recovery (days): Weibull Function Fit',
       y='Observed Time to Recover (days)',
       subtitle='Fires between 2003-2004')+
  # facet_wrap(~cut_interval(delta_vi_12mo, 4))+
  theme_linedraw()
ggsave(filename = 'figures/ObsTTR_WeibullTTR_2003-04_fires.png')


n_w %>% 
  mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
  mutate(pred_ttr = w_ttr(Asym, Drop, lrc, pwr)) %>% 
  filter(isConv==TRUE) %>%
  filter(pred_ttr50 < 2000) %>%  #pull(pred_ttr50) %>% hist
  filter(pred_ttr < 5000) %>%  #pull(pred_ttr50) %>% hist
  ggplot(data=.,aes(pred_ttr50, pred_ttr))+
  geom_point()



n_w %>% 
  mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
  mutate(pred_ttr50 = w_ttr(Asym, Drop, lrc, pwr)) %>%   
  filter(isConv==TRUE) %>%
  filter(pred_ttr50 < 1000) %>%  #pull(pred_ttr50) %>% hist
  ggplot(data=.,aes(x,y,fill=pred_ttr50))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c()



expand_grid(n_w[sample(.N,10)], post_days=1:1500) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>% #pull(pred_ttr50) %>% summary
  ggplot(data=.,aes(post_days, pred, group=id))+
  geom_line()


merge(mdat, n_w, by=c("x","y","id")) %>% 
  lazy_dt() %>% 
  group_by(id) %>% 
  filter(ndvi_anom==min(ndvi_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(-Drop,ndvi_anom))+geom_point()+geom_abline()



tmp <- mdat[id%in%sample(vec_ids,100)][ttr>100][,fn_w(.SD), by=.(x,y,id)]
tmp %>% summary

expand_grid(tmp, post_days=1:1500) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(post_days, pred, group=id,color=lrc))+
  geom_line(lwd=0.1)+
  scale_color_viridis_c(end=0.9)


with(tmp[1,], w_ttr(0.2,Drop,lrc,pwr))

(exp(26.877)*log(10.0*x/(10.0*0.00001 - 1.0)))^(1.0/4.86)
curve((exp(26.877)*log(10.0*x/(10.0*0.00001 - 1.0)))^(1.0/4.86), 0,1)


tmp %>% filter(Drop<1) %>% 
  mutate(pred_ttr = (exp(-lrc)*log(10.0*Drop/(10.0*0.1)))^(1.0/pwr)) %>% #pull(pred_ttr) %>% hist
  ggplot(data=.,aes(Drop, pred_ttr, color=lrc))+
  geom_point()+
  scale_color_viridis_c(end=0.9)


Asym*exp(-b2*b3^x)
curve(SSgompertz(x,Asym = -0.1, b2 = -1,b3=0.005),0,1000)
Asym+(R0-Asym)*exp(-exp(lrc)*input)
curve(SSasymp(x,Asym = 0, R0 = -1,lrc = 700), 0,1000)


Asym-Drop*exp(-exp(lrc)*x^pwr)
curve(SSweibull(x, Asym = 0, Drop = 1, lrc = 0.01, pwr = 10), 0,1000)
