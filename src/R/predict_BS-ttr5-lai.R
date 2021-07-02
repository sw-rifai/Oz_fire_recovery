# To do!
# Get the min nbr anom and date of fire for the BS epoch
# Merge with the LAI anom etc
# Merge with the d_soil, etc
# Check all the predictors are there for the GAM and RF models to predict TTR

library(dtplyr)
library(data.table)
library(tidyverse)
library(lubridate)

# Load and merge data sets --------------------------------
bs_start <- ymd("2019-09-01")
# fires and NBR
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
                            col_select = c("x","y","date","id","fire_doy","nbr_anom"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id","fire_doy","nbr_anom"))
gc(full=TRUE)
tmp3 <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
tmp3 <- tmp3[date>=bs_start]

dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                           col_select = c("x","y","id","date","slai","slai_anom_12mo","malai","slai_anom_3mo"))
dat[,`:=`(slai_12mo = slai_anom_12mo+malai)]

dat <- dat[date>=bs_start]
dat <- merge(dat,tmp3 %>% select(-x,-y) %>% as.data.table(),by=c("id","date"))

# landscape covars 
d_soil <- arrow::read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
dat <- merge(dat,d_soil, by='id', allow.cartesian = TRUE)

# species ----------------------
dom <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
dom <- dom[,.(x,y,id,dom_sp)]
dat <- merge(dat, dom, by=c("id","x","y"))

# filtering for fires only in 2019!
# NEEDS UPDATE: Precip and other climate data not yet available for full 2020
dat_bs <- dat[date<=ymd("2020-05-01")] # not subsetting to end of 2019 yet so I can get the min_nbr_anom(up to 3mo post fire)

dat_bs <- dat_bs[id %in% dat_bs[date<=ymd("2019-12-31")][fire_doy>0][,.(nobs = .N), by='id'][nobs==1]$id]


# Calculate min_nbr_anom
# dat2: Black Summer fires -----------------------------------------------------
gc(full=TRUE)

# get the date of the fire
firedate <- dat_bs %>% lazy_dt() %>%
  filter(date >= bs_start & date <= ymd("2019-12-31")) %>% 
  filter(fire_doy > 0) %>%
  group_by(id) %>%
  mutate(date_fire1 = date) %>%
  ungroup() %>%
  select(id,date_fire1) %>%
  as.data.table()
dat_bs <- merge(dat_bs,firedate,by='id',allow.cartesian=TRUE)

gc(full=TRUE)
dat_bs <- dat_bs %>% lazy_dt() %>%
  mutate(days_since_fire = as.double(date - date_fire1)) %>%
  as.data.table()
gc(full=TRUE)
d_min_nbr <- dat_bs %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
d_L0 <- dat_bs %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 365) %>% 
  group_by(id) %>% 
  filter(slai == min(slai,na.rm=TRUE)) %>% 
  summarize(L0_rs = slai, 
            K_rs = malai) %>% 
  # summarize(L0_rs = min(slai,na.rm=TRUE), 
  #           K_rs = max(slai_12mo)) %>% 
  ungroup() %>% 
  mutate(ldk = L0_rs/K_rs) %>% 
  as.data.table()
d_min_nbr <- merge(d_min_nbr, d_L0, by='id')
gc(full=TRUE)
dat_bs <- merge(dat_bs,d_min_nbr,by='id')
gc(full=TRUE)

# filter down to the date of the fire
dat <- dat_bs[date==date_fire1]

# pre & post_fire climate -------------------------
post_clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                                 col_select = c("x","y","date",
                                                "map","mapet","mappet","mavpd15",
                                                "matmin","matmax",
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

# next step breaks tune_grid???
post_clim <- post_clim[date>=ymd("2018-01-01")][order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                                       n = 3,fill = NA,align='right'), by=.(x,y)]


#  Attach clim pixel id to VI ------------------------------------
coords_vi <- dat %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- sf::st_as_sf(coords_vi, coords = c("x","y"))
sf::st_crs(coords_vi) <- sf::st_crs(4326)
coords_clim <- unique(post_clim[,.(x,y)])
coords_clim_sf <- sf::st_as_sf(coords_clim, coords = c('x','y'))
sf::st_crs(coords_clim_sf) <- sf::st_crs(4326)
nn_coords <- RANN::nn2(
  coords_clim_sf %>% sf::st_coordinates(),
  coords_vi %>% sf::st_coordinates(), 
  k=1
)
coords_clim <- coords_clim %>% mutate(idx_clim = row_number()) %>% as.data.table()
gc(full=TRUE)
coords_vi <- coords_vi %>% sf::st_drop_geometry() %>% as.data.table()
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


# Create 'dat'
dat <- merge(dat, post_clim %>% select(-x,-y) %>% as.data.table, 
             by=c('idx_clim','date'),allow.cartesian = TRUE)

dat <- dat[,`:=`(post_precip_anom_frac = post_precip_anom_12mo/map,
                 precip_anom_frac = precip_anom_12mo/map,
                 post_vpd15_anom_frac = post_vpd15_anom_12mo/mavpd15)]
dat[,fire_month:=month(date_fire1)]



# Final data subset ----------------------------------------------------------
dat <- dat[vc %in% c(2,3,5)][fire_month %in% c(9,10,11,12,1,2)]
# nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
# sp_fac <-
#   unique(dat[dom_sp %in% nobs[rank <= 30]$dom_sp][,.(dom_sp,hab_hnd)]) %>%
#   .[order(hab_hnd)] %>%
#   .[,sp_fac := forcats::fct_inorder(dom_sp)]
# Final mutations
dat[,fire_month_f := lubridate::month(date_fire1,label = TRUE,abbr = TRUE)][
  ,fire_month_f := factor(fire_month_f,
                          levels=c("Sep","Oct","Nov","Dec","Jan","Feb"),
                          ordered = TRUE)][,dom_sp_f := factor(dom_sp)]
dat[,vc_name_f := factor(vc_name)]




# Write data for preds to disk -------------------------------------------------
arrow::write_parquet(dat, sink=paste0("../data_general/proc_data_Oz_fire_recovery/dataForBSpreds_",Sys.Date(),".parquet"))



################################################################################
# GAM preds ---------------------------------------------------------------
# MGCV comparison -----------
library(mgcv)
d_rf <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai.parquet") %>% 
  as_tibble() %>% 
  filter(is.na(ttr5_lai)==F) %>% 
  filter(is.na(sand)==F & is.na(pH)==F & is.na(elevation)==F&
           is.na(slope)==F & is.na(aspect)==F & is.na(pre_fire_slai_anom_12mo)==F) %>% 
  mutate(target = ttr5_lai-364)

d_train <- d_rf %>% 
  filter(!year %in% c(2003, 2007, 2012, 2015)) %>% 
  # filter((year%%2)==1) %>%
  # slice_sample(prop=0.5) %>%
  as_tibble()
d_test <- d_rf %>% 
  filter(year %in% c(2003, 2007, 2012, 2015)) %>% 
  # filter((year%%2)==0) %>% 
  # slice_sample(prop=0.5) %>% 
  as_tibble()

b1 <- mgcv::bam(target ~
                  te(log(elevation+10), 
                     log(der),k=5)+
                  s(pH,k=5)+
                  s( I(matmax-matmin), k=5)+
                  te(min_nbr_anom,malai, k=5)+
                  te(map, mavpd15,k=5)+
                  te(map, precip_anom_12mo, post_precip_anom_12mo,k=5),
                  # fire_month_f +
                  # vc_name_f + 
                  # s(map,k=5,bs='cs')+
                  # s(mavpd15,k=5,bs='cs')+
                  # s(I(matmax-matmin),k=5,bs='cs')+
                  # s(min_nbr_anom, k=5, bs='cs')+
                  # s(malai,k=5,bs='cs')+
                  # s(elevation,k=5,bs='cs')+
                  # s(pH, k=5,bs='cs')+
                  # s(sand, k=5, bs='cs')+
                  # s(silt, k=5, bs='cs')+
                  # s(slai_anom_12mo, k=5, bs='cs')+
                  # s(precip_anom_12mo, k=5,bs='cs')+
                  # s(post_vpd15_anom_12mo,k=5,bs='cs')+
                  # s(post_precip_anom_12mo,k=5,bs='cs'),
                   # fire_month_f+
                   # te(min_nbr_anom, elevation,malai)+ 
                   # te(min_nbr_anom, pre_fire_slai_anom_12mo,post_precip_anom_12mo),
                 data=d_rf, 
                family=Tweedie(p=1.67, 
                          link = 'log'),
                # family=Gamma(link='identity'),
                # family=nb(link = 'identity'),
                # family=negbin(theta=1, link='log'),
                 # family=Tweedie(p=1.15),
                 # weights = d_train$ttr5_lai,
                 # weights = as.numeric(scale(d_train$ttr5_lai, center = F)),
                 # family=Gamma(link='log'),
                 select=TRUE, 
                 discrete=TRUE)
summary(b1)
b1$
plot(b1, scheme=2, pages=1)
qq.gam(b1)
rsq_vec(d_test$ttr5_lai, 
             (predict(b1, newdata=d_test,type='response')+364))
rmse_vec(d_test$ttr5_lai, 
         (predict(b1, newdata=d_test,type='response')+364))
rmse_vec(truth=d_test$ttr5_lai,
          estimate = rep(mean(d_train$ttr5_lai,na.rm=T),dim(d_test)[1]))

d_rf %>% 
  sample_n(10000) %>% 
  mutate(pred = (predict(b1, newdata=., type='response'))+364) %>% #pull(pred) %>% summary
  ggplot(data=.,aes(ttr5_lai,ttr5_lai-pred))+
  ggpointdensity::geom_pointdensity()+
  geom_hline(col='red', aes(yintercept=0))+
  geom_smooth()+
  scale_color_viridis_c(option='H')+
  labs(x='Observations', 
       y='Obs. - Pred.')



# Plot GAM BS predictions -----------------------------------------------------
bs <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/dataForBSpreds_2021-06-03.parquet")
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.01) %>% 
  select(NAME_1)

bs <- bs %>% 
  as_tibble() %>% 
  # rename(x=x.x, y=y.x) %>% 
  mutate(pred = predict(b1, newdata=., type='response')+364) 

bs %>% 
  ggplot(data=.,aes(x,y,fill=pred))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c(option='H', 
                       limits=c(365,3000), 
                       oob=scales::squish)+
  labs(x=NULL,y=NULL,
       fill="days")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))


p_vic <- bs %>% 
  ggplot(data=.,aes(x,y,fill=pred))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(144.75,149.5),
           ylim = c(-38.5,-35), expand = FALSE)+
  # scico::scale_fill_scico(palette='roma', 
  #                         direction = -1,
  #                         limits=c(365,2500), 
  #                         oob=scales::squish)+
  scale_fill_viridis_c(option='H',
                       limits=c(365,2500),
                       oob=scales::squish)+
  labs(x=NULL,y=NULL,
       fill="days", 
       title='Victoria & ACT')+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(0,1),
        legend.direction = 'vertical',
        legend.justification = c(0,1)); p_vic
p_nsw <- bs %>% 
  ggplot(data=.,aes(x,y,fill=pred))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(149.75,151.5),
           ylim = c(-35.5,-32.5), expand = FALSE)+
  # scico::scale_fill_scico(palette='roma', 
  #                         direction = -1,
  #                         limits=c(365,2500), 
  #                         oob=scales::squish)+
  scale_fill_viridis_c(option='H',
                       limits=c(365,2500),
                       oob=scales::squish)+
  labs(x=NULL,y=NULL,
       fill="days", 
       title="Coastal Central NSW")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = 'none', #c(1,0), 
        legend.justification = c(1,0)); p_nsw

library(patchwork)
p_vic + p_nsw + plot_annotation(tag_levels = 'a', 
                                tag_prefix = '(', 
                                tag_suffix = ')', 
                                title = 'Black Summer Fires- Time to Recover Forecast')+
  plot_layout()
ggsave(filename = "figures/predict_BS-TTR-GAM.png", 
       width=23, 
       height=15, 
       units='cm', 
       dpi=350)


################################################################################
# PART 3: RANDOM FOREST  -----
# restart R
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

# load pre-processed data
d_rf <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai.parquet") %>% 
  as_tibble() %>% 
  filter(is.na(sand)==F & is.na(pH)==F & is.na(elevation)==F&
           is.na(slope)==F & is.na(aspect)==F & is.na(pre_fire_slai_anom_12mo)==F)

d_train <- d_rf %>% filter((year%%2)==0) %>% 
  slice_sample(prop=0.5) %>% as_tibble()
d_test <- d_rf %>% filter((year%%2)==1) %>% 
  slice_sample(prop=1) %>% as_tibble()


# Set up hyp-param grid ---------------------------------------------------
set.seed(123)
d_split <- initial_split(d_train, 
                         strata = 'ttr5_lai', 
                         breaks=5, 
                         pool=0.1, 
                         prop = 0.8)
d_train <- training(d_split)
d_test <- d_test # testing(d_split)

# Specify RF mod specs ------------------
rf_spec <- rand_forest(mode='regression', 
                       mtry=tune(), 
                       trees=tune(),
                       min_n=tune()) %>% 
  set_engine("ranger",
             regularization.factor = tune("regularization"), 
             num.threads=5)

# Specify CV fold data ------------------------
d_folds <- vfold_cv(d_train, v=3)

# Specify 'Recipe' --------------------
rf_rec <- recipe(ttr5_lai ~ min_nbr_anom+
                   fire_month_f + 
                   malai+
                   map+mapet+mavpd15+
                   matmax+matmin+
                   elevation+sand+pH+silt+
                   der+
                   slope+aspect+vpd15_anom_3mo+
                   pre_fire_slai_anom_12mo + 
                   post_precip_anom_12mo+
                   post_vpd15_anom_12mo+
                   precip_anom_12mo, data=d_train) %>% 
  # step_log(ttr5_lai) %>% 
  # update_role(year, new_role = 'ID') %>% 
  prep(training = d_train)

# Specify 'Workflow' ------------------
rf_wf <- workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rf_rec)

gc(full=T)

# set up the grid of tuning parameters
rf_grid <- grid_latin_hypercube(finalize(mtry(range=c(2,8)), d_train),
                                trees(range=c(200,2500)),
                                min_n(range=c(1,10)),
                                regularization = regularization_factor(range = c(0.5,1)),
                                size=20)
rf_grid

doParallel::registerDoParallel()
# doParallel::stopImplicitCluster()
set.seed(234)
rf_res <- tune_grid(
  rf_wf,
  resamples = d_folds,
  grid = rf_grid,
  control = control_grid(save_pred = TRUE)
)

rf_res
collect_metrics(rf_res)


# Plot the model fits to the tuning parameters
rf_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")

rf_res %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry, trees, min_n, regularization) %>%
  pivot_longer(mtry:regularization,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "R2")

best_rmse <- select_best(rf_res, metric = 'rmse')
best_rmse

final_rf <- finalize_model(rf_spec, parameters = best_rmse)
final_rf %>% translate()

# Examine variable importance ------------------
library(vip)

test <- workflow() %>% 
  add_recipe(., recipe = rf_rec) %>% 
  add_model(final_rf %>% 
              set_engine("ranger", 
                         importance = "permutation", 
                         num.threads=10)) %>% 
  fit(data=d_train)

test$fit$fit %>% vip::vip(geom='point')


# Fit Final Model -------------------------
final_wf <- workflow() %>%
  add_model(final_rf) %>% 
  add_recipe(rf_rec)

# final_res <- fit(final_wf, training(d_split))

final_res <- final_wf %>%
  last_fit(d_split)

final_res %>%
  collect_metrics()

final_res %>% 
  collect_predictions() %>% 
  sample_n(1000) %>% 
  ggplot(data=.,aes(ttr5_lai, .pred))+
  # ggplot(data=.,aes(exp(ttr5_lai), exp(.pred)-exp(ttr5_lai)))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')

final_res %>% 
  collect_predictions() %>% 
  sample_n(1000) %>% 
  ggplot(data=.,aes(ttr5_lai, ttr5_lai-.pred))+
  # ggplot(data=.,aes(exp(ttr5_lai), exp(.pred)-exp(ttr5_lai)))+
  geom_point()+
  geom_hline(aes(yintercept=0),col='red')+
  geom_smooth()


test_prediction <- final_wf %>%
  # fit the model on all the training data
  fit(
    data    = d_train
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = d_test) %>%
  bind_cols(d_test)

rsq_trad_vec(test_prediction$ttr5_lai, 
             test_prediction$.pred)

# SCRATCH ------------------------------------------

r1 <- ranger::ranger(log(ttr5_lai) ~ 
                       min_nbr_anom+
                       fire_month_f+
                       malai+
                       map+mapet+mavpd15+
                       matmax+matmin+
                       elevation+sand+pH+silt+
                       der+
                       slope+aspect+vpd15_anom_3mo+
                       pre_fire_slai_anom_12mo + 
                       post_precip_anom_12mo+
                       post_vpd15_anom_12mo+
                       precip_anom_12mo, 
                     data=d_train, 
                     # splitrule = 'maxstat',
                     # minprop = 0.5,
                     # alpha = 0.5,
                     # replace = F,
                     mtry=5,
                     min.node.size = 2,
                     regularization.factor = 0.75,
                     # # treetype='regression',
                     # max.depth = 10,
                     num.trees = 500,
                     # regularization.usedepth = TRUE,
                     # case.weights = s,
                     # case.weights = abs(fn_dn(d_train$ttr5_lai) - 8e-04),
                     num.threads = 10)
r1
rsq_trad_vec(d_test$ttr5_lai,
             exp(predict(r1,data=d_test)$predictions))


d_test %>% 
  sample_n(10000) %>% 
  mutate(pred = predict(r1, data=.)$predictions) %>% 
  ggplot(data=.,aes(ttr5_lai, ttr5_lai - exp(pred)))+
  ggpointdensity::geom_pointdensity()+
  geom_hline(col='red',aes(yintercept=0))+
  geom_smooth()+
  scale_color_viridis_c(option='H')

d_test %>% 
  sample_n(10000) %>% 
  mutate(pred = predict(r1, data=.)$predictions) %>% 
  ggplot(data=.,aes(exp(pred), ttr5_lai - exp(pred)))+
  ggpointdensity::geom_pointdensity()+
  geom_hline(col='red',aes(yintercept=0))+
  geom_smooth()+
  scale_color_viridis_c(option='H')
