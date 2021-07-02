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
                           col_select = c("x","y","id","date","slai","slai_anom_12mo","malai"))
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
  mutate(target = ttr5_lai-364) %>% 
  mutate(ttr_ocat = case_when(ttr5_lai <= 366 ~ 1,
                           between(ttr5_lai, 367, 365*2)~2,
                           between(ttr5_lai, 731, 1096)~3,
                           between(ttr5_lai, 1097, 1460)~4,
                           between(ttr5_lai, 1461, Inf)~5,
  ))

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

b1 <- mgcv::gam(ttr_ocat ~
                  # te(log(elevation+10),
                  #    log(der),k=5)+
                  s(pH,k=5)+
                  # s( I(matmax-matmin), k=5)+
                  te(min_nbr_anom,malai, k=5)+
                  # te(map, mavpd15,k=5)+
                  te(mavpd15, post_vpd15_anom_12mo, k=5)+
                  te(map, precip_anom_12mo, post_precip_anom_12mo,k=5)+
                  fire_month_f,
                
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
                data=d_rf %>% mutate(ttr_ocat = as.numeric(ttr_ocat)), 
                family=ocat(R=5),
                # family=Tweedie(p=1.67, 
                #                link = 'log'),
                # family=Gamma(link='identity'),
                # family=nb(link = 'identity'),
                # family=negbin(theta=1, link='log'),
                # family=Tweedie(p=1.15),
                # weights = d_train$ttr5_lai,
                # weights = as.numeric(scale(d_train$ttr5_lai, center = F)),
                # family=Gamma(link='log'),
                select=TRUE)
summary(b1)
plot(b1, scheme=2, pages=1)
m_pred <- predict(b1,newdata = d_rf,type="response",se=F)
colMeans(m_pred)
predict(b1) %>% dim
dim(d_rf)

yardstick::average_precision_vec(truth=factor(d_rf$ttr_ocat), 
                                 m_pred)


library(probably)
tibble(m_pred)
as.data.frame(m_pred) %>% 
  as_tibble() %>% 
  set_names(paste0("c",1:5)) %>% 
  bind_cols(class= factor(d_rf$ttr_ocat,levels = c(1:5),labels=paste0("c",1:5))) %>% 
  mutate(.class_pred = make_class_pred(
  c1, c2, c3, c4, c5, 
  levels=(paste0("c",1:5)), 
  min_prob=0.3
  )) %>% 
  select(class, .class_pred) %>% 
  count(truth=class, .class_pred)



species_probs %>%
  mutate(
    .class_pred = make_class_pred(
      .pred_bobcat, .pred_coyote, .pred_gray_fox,
      levels = levels(Species),
      min_prob = .4
    )
  )


make_class_pred(m_pred, levels=c(1,2,3,4,5),min_prob = 0.2) %>% class

yardstick::average_precision_vec(truth=factor(d_rf$ttr_ocat), 
                             m_pred)
caret::confusionMatrix(factor(d_rf$ttr_ocat), 
                       m_pred)
caret::multiClassSummary(factor(d_rf$ttr_ocat),m_pred)



caret::confusionMatrix(iris$Species, sample(iris$Species))
newPrior <- c(.05, .8, .15)
names(newPrior) <- levels(iris$Species)
caret::confusionMatrix(iris$Species, sample(iris$Species))

qq.gam(b1)
rsq_vec(d_test$ttr5_lai, 
        (predict(b1, newdata=d_test,type='response')+364))
rmse_vec(d_test$ttr5_lai, 
         (predict(b1, newdata=d_test,type='response')+364))


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
library(themis)

# load pre-processed data
d_rf <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai.parquet") %>% 
  as_tibble() %>% 
  filter(is.na(ttr5_lai)==F) %>% 
  filter(is.na(sand)==F & is.na(pH)==F & is.na(elevation)==F&
           is.na(slope)==F & is.na(aspect)==F & is.na(pre_fire_slai_anom_12mo)==F) %>% 
  mutate(target = ttr5_lai-364) %>% 
  mutate(ttr_ocat = case_when(ttr5_lai <= 366 ~ 1,
                              between(ttr5_lai, 367, 365*2)~2,
                              between(ttr5_lai, 731, 1096)~3,
                              between(ttr5_lai, 1097, 1460)~4,
                              between(ttr5_lai, 1461, Inf)~5,
  )) %>% 
  mutate(ttr_ocat = factor(ttr_ocat, levels = c(1,2,3,4,5), ordered = TRUE))

d_train <- d_rf %>% filter((year%%2)==0) %>% 
  slice_sample(prop=0.5) %>% as_tibble()
d_test <- d_rf %>% filter((year%%2)==1) %>% 
  slice_sample(prop=1) %>% as_tibble()


# Set up hyp-param grid ---------------------------------------------------
set.seed(123)
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
  # step_log(ttr5_lai) %>% 
  # update_role(year, new_role = 'ID') %>% 
  # step_dummy(all_nominal(), -ttr_ocat) %>% 
  step_smote(ttr_ocat) %>% 
  prep()#%>% 
  # prep(training = d_train) 

# Specify 'Workflow' ------------------
rf_wf <- workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rf_rec)

gc(full=T)
doParallel::registerDoParallel()
# set.seed(74403)
# ranger_rs <-
#   fit_resamples(rf_wf,
#                 resamples = d_folds,
#                 control = control_resamples(save_pred = TRUE)
#   )


# set up the grid of tuning parameters
rf_grid <- grid_latin_hypercube(finalize(mtry(range=c(5,6)), d_train),
                                trees(range=c(1000,1500)),
                                min_n(range=c(5,6)),
                                # regularization = regularization_factor(range = c(0.5,1)),
                                size=5)
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

rf_res$.notes[[1]]
collect_metrics(rf_res)


# Plot the model fits to the tuning parameters
rf_res %>%
  collect_metrics() %>%
  filter(.metric == "accuracy") %>%
  select(mean, mtry, trees, min_n) %>%
  pivot_longer(mtry:min_n,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "Accuracy")

rf_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry, trees, min_n) %>%
  pivot_longer(mtry:min_n,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "roc_auc")

best_accuracy <- select_best(rf_res, metric = 'accuracy')
best_accuracy

final_rf <- finalize_model(rf_spec, parameters = best_accuracy)
final_rf %>% translate()

# Examine variable importance ------------------
library(vip)

test <- workflow() %>% 
  add_recipe(., recipe = rf_rec) %>% 
  add_model(final_rf %>% 
              set_engine("ranger", 
                         respect.unordered.factors='order',
                         importance = "permutation", 
                         num.threads=10)) %>% 
  fit(data=training(d_split2))

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
  roc_curve(ttr_ocat, .pred_1:.pred_5) %>% 
  ggplot(data=., aes(1-specificity, sensitivity))+
  geom_abline(lty=3,col='gray80')+
  geom_path()+
  facet_wrap(~.level)

conf_mat_resampled(final_res, tidy = FALSE) %>%
  autoplot(type='heatmap')


final_res %>% 
  collect_predictions() %>% 
  pull(ttr_ocat) %>% as.numeric %>% hist
final_res %>% 
  collect_predictions() %>% 
  pull(.pred_class) %>% as.numeric %>% hist

final_res %>% 
  collect_predictions() %>% 
  roc_curve(ttr_ocat, .pred_class) %>% 
  autoplot()

final_res %>% 
  collect_predictions() %>% 
  sample_n(1000) %>% 
  ggplot(data=.,aes(ttr_ocat, .pred_class))+
  geom_smooth()
  # ggplot(data=.,aes(exp(ttr5_lai), exp(.pred)-exp(ttr5_lai)))+
  geom_point(position = 'jitter')+
  # geom_abline(col='red')+
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

r1 <- ranger::ranger(ttr_ocat ~ 
                       min_nbr_anom+
                       fire_month_f+
                       vc_name_f + 
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
                     respect.unordered.factors = TRUE,
                     # splitrule = 'maxstat',
                     # minprop = 0.5,
                     # alpha = 0.5,
                     # replace = F,
                     mtry=6,
                     min.node.size = 6,
                     # regularization.factor = 0.75,
                     # # treetype='regression',
                     # max.depth = 10,
                     num.trees = 1400,
                     # regularization.usedepth = TRUE,
                     # case.weights = s,
                     # case.weights = abs(fn_dn(d_train$ttr5_lai) - 8e-04),
                     num.threads = 10)
r1
rsq_trad_vec(d_test$ttr5_lai,
             exp(predict(r1,data=d_test)$predictions))

predict(r1, data=d_test) %>% str

d_test %>% 
  mutate(pred = predict(r1, data=.)$predictions) %>% 
  mutate(ttr_ocat = as.numeric(ttr_ocat)) %>% 
  roc_curve(ttr_ocat, pred)
roc_curve(d_test[1:100,]$ttr_ocat, 
          predict(r1, data = d_test[1:100,])$predictions
          )
hpc_cv$obs %>% class
hpc_cv$VF %>% class
hpc_cv %>%
  # head() %>% 
  filter(Resample == "Fold01") %>%
  roc_curve(obs, VF:L) %>%
  autoplot()

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


### SCRATCH ###################################################################
library(usemodels)
use_ranger(ttr_ocat~., data=d_train) %>% 
  translate()


model <- gam(mpg ~ s(wt) + cyl, data = mtcars)
performance::model_performance(model)

model <- gam(vs ~ s(wt) + (mpg), data = mtcars, family = "binomial")
performance::model_performance(model)
summary(model)
