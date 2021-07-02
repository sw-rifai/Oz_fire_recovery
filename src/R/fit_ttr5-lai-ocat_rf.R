library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(arrow)
library(lubridate)

# PART 1: PREP  ----------------------------------------------------------------
# Data load --------------------------
dttr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS_2021-06-05 13:01:38.parquet")

# landscape covars ------------
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")

# species ----------------------
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-SEOZcoastal-grid.parquet")
dom <- dom[,.(x,y,id,dom_sp)]

# Do the small merges up-front to catch which ds are missing data
dat <- merge(dttr, d_soil %>% select(-x,-y) %>% as.data.table(), by='id',allow.cartesian = TRUE)
dat <- merge(dat, dom %>% select(-x,-y) %>% as.data.table(), by='id',allow.cartesian = TRUE)


# pre & post_fire climate -------------------------
post_clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                                 col_select = c("x","y","date",
                                                "map","mapet","mappet","mavpd15",
                                                "matmin","matmax",
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
post_clim <- post_clim[order(x,y,date)][
  , vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                 n = 3,fill = NA,align='right'), by=.(x,y)]


#  Attach clim pixel id to VI ------------------------------------
coords_vi <- lazy_dt(dttr) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_clim <- unique(post_clim[,.(x,y)])
coords_clim_sf <- st_as_sf(coords_clim, coords = c('x','y'))
st_crs(coords_clim_sf) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_clim_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
coords_clim <- coords_clim %>% mutate(idx_clim = row_number()) %>% as.data.table()
gc(full=TRUE)
coords_vi <- coords_vi %>% st_drop_geometry() %>% as.data.table()
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



# Filter dat to just Eucs in the Oct-Feb burning -------------------
dat <- dat[vc %in% c(2,3,5)][month %in% c(9,10,11,12,1,2)]
nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
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

# Convert to tibbles ---------------------------------------------------------
d_fit <- dat %>% as_tibble() %>% 
  filter(is.na(sand)==F & is.na(pH)==F & is.na(elevation)==F&
           is.na(slope)==F & is.na(aspect)==F & is.na(pre_fire_slai_anom_12mo)==F) %>% 
  mutate(target = ttr5_lai-364) %>% 
  mutate(ttr_ocat = case_when(ttr5_lai <= 366 ~ 1,
                              between(ttr5_lai, 367, 365*2)~2,
                              between(ttr5_lai, 731, 1096)~3,
                              between(ttr5_lai, 1097, 1460)~4,
                              between(ttr5_lai, 1461, Inf)~5,
                              is.na(ttr5_lai)==T ~ 5
  ))



arrow::write_parquet(d_fit, sink="../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")


#*******************************************************************************
#* PART 2: FITTING THE RF 
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


# Specify grid of tuning parameters ---------------------------------
rf_grid <- grid_latin_hypercube(finalize(mtry(range=c(3,8)), d_train),
                                trees(range=c(500,1500)),
                                min_n(range=c(3,9)),
                                # regularization = regularization_factor(range = c(0.5,1)),
                                size=30)
rf_grid


# Start tuning in parallel -----------------------------------------
doParallel::registerDoParallel()
# doParallel::stopImplicitCluster()
set.seed(333)
rf_res <- tune_grid(
  rf_wf,
  resamples = d_folds,
  grid = rf_grid,
  control = control_grid(save_pred = TRUE)
)


# Plot the model fits to the tuning parameters ----------------------
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

#        mtry trees min_n .config              
# <int> <int> <int> <chr>                
#   1     3   998     3 Preprocessor1_Model26

best_roc <- select_best(rf_res, metric = 'roc_auc')
best_roc

#        mtry trees min_n .config              
# <int> <int> <int> <chr>                
#   1     4   745     4 Preprocessor1_Model06

# final_rf <- finalize_model(rf_spec, parameters=tibble(mtry=3, trees=998, min_n=3))
final_rf <- finalize_model(rf_spec, parameters = best_roc)
final_rf %>% translate()

# Examine variable importance ------------------
test <- workflow() %>% 
  add_recipe(., recipe = rf_rec) %>% 
  add_model(final_rf %>% 
              set_engine("ranger", 
                         respect.unordered.factors='order',
                         importance = "permutation", 
                         num.threads=10)) %>% 
  fit(data=training(d_split))

test$fit$fit %>% 
  vip::vip(geom='point')+
  theme_linedraw()
ggsave(filename = paste0("figures/vip_ttr5-lai-ocat_randomForest_",Sys.Date(),".png"), 
       width=15, 
       height=10, 
       units='cm', 
       dpi=350)

# Fit Final Model -------------------------
final_wf <- workflow() %>%
  add_model(final_rf) %>% 
  add_recipe(rf_rec)

final_res <- final_wf %>%
  last_fit(d_split2)

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
  autoplot(type='heatmap')+
  scale_fill_viridis_c(alpha=0.3,direction = -1, 
                       limits=c(0, 25000))+
  scale_x_discrete(breaks=1:5, 
                   labels=c("≤1","2","3","4","≥5"))+
  scale_y_discrete(breaks=1:5, 
                   labels=c("≤1","2","3","4","≥5"))+
  labs(x='Observed Time to Recover (years)', 
       y='Predicted Time to Recover (years)')
ggsave(filename = paste0("figures/confusion-matrix_ttr5-lai-ocat_randomForest_",Sys.Date(),".png"), 
       width=15, 
       height=10, 
       units='cm', 
       dpi=350)

test_prediction <- final_wf %>%
  # fit the model on all the training data
  fit(
    data    = training(d_split2)
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = d_rf) %>%
  bind_cols(d_rf)

test_prediction %>% 
  select(.pred_class, ttr_ocat) %>% 
  mutate(p = as.numeric(.pred_class), 
         o = as.numeric(ttr_ocat)) %>% 
  ggplot(data=.,aes(p,o))+
  geom_point(position='jitter',size=0.1,alpha=0.1)+
  geom_smooth(method='lm')


sp_count <- test_prediction %>% 
  rename(dom_sp = dom_sp.x) %>% 
  group_by(dom_sp) %>% 
  summarize(nobs = n()) %>% 
  ungroup() %>% 
  mutate(rank = rank(-nobs))

test_prediction %>% filter(
  dom_sp.x %in% (sp_count %>% 
    filter(rank <= 40) %>% pull(dom_sp))
  ) %>% 
  mutate(dom_sp = factor(dom_sp.x)) %>% 
  mutate(ttr_cat = fct_rev(ttr_ocat)) %>% 
  # select(dom_sp.x, ttr_ocat, .pred_class) %>% 
  # mutate(.pred_class = factor(.pred_class, levels=1:5,labels=1:5,ordered = TRUE)) %>% 
  # pivot_longer(cols=c("ttr_ocat",".pred_class")) %>% 
  ggplot(data=.,aes(y=dom_sp, 
                    fill=ttr_ocat))+
  geom_bar(position=position_fill(reverse=T))+
  scale_y_discrete(limits = rev)+
  scale_x_continuous(expand=c(0,0), 
                     limits=c(0,1))+
  scale_fill_viridis_d(option='B', end=0.95, direction = 1)+
  labs(y=NULL,
       x='proportion',
       fill='Years')+
  # guides(fill = guide_colorsteps())
  theme(legend.position = 'bottom')

test_prediction %>% 
  filter(dom_sp.x == "Eucalyptus mackintii") %>% 
  pull(ttr_ocat) %>% 
  table

test_prediction %>% 
  filter(dom_sp.x == "Eucalyptus regnans") %>% 
  pull(ttr_ocat) %>% 
  table


# Predict the Black Summer Recovery Time ----------------------------------
library(data.table)
library(lubridate)
d_pred <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/dataForBSpreds_2021-06-06.parquet") %>% 
  rename(pre_fire_slai_anom_3mo = slai_anom_3mo, 
         pre_fire_slai_anom_12mo = slai_anom_12mo) %>% 
  as.data.table()
d_pred <- d_pred[is.na(elevation)==F & 
                 is.na(sand)==F &
                 is.na(pH)==F &
                 is.na(silt)==F & 
                 is.na(der)==F &
                 is.na(slope)==F & 
                 is.na(aspect)==F]

d_pred2 <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/rescaled-era5land_202101_202103.parquet")
unique(d_pred[,.(x,y,elevation,sand,pH,silt,der,slope,aspect)])




d_pred <- d_pred %>% as_tibble()

pred_mod <- final_wf %>%
  # fit the model on all the training data
  fit(
    data    = training(d_split2)
  )
write_rds(pred_mod,file = paste0('outputs/rf-mod-ttr5-ocat_',Sys.Date(),'.rds'),compress = 'gz')


d_pred <- pred_mod %>% 
  # use the training model fit to predict the test data
  predict(new_data = d_pred) %>%
  bind_cols(d_pred)

d_pred$.pred_class %>% as.numeric() %>% hist(6)

