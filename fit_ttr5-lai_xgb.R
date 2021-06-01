library(data.table)
# library(dtplyr); 
# library(sf); library(stars)
library(arrow)
library(tidyverse); 


# Data load and prep --------------------------
dttr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-05-24 16:20:29.parquet")

# landscape covars ------------
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")

# species ----------------------
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
dom <- dom[,.(x,y,id,dom_sp)]


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

# next step breaks tune_grid???
post_clim <- post_clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                                       n = 3,fill = NA,align='right'), by=.(x,y)]


#  Attach clim pixel id to VI ------------------------------------
coords_vi <- dttr %>% select(x,y,id) %>% distinct() %>% as.data.table()
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
dat <- merge(dttr, coords_vi, by='id')
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

dat <- merge(dat, d_soil, by='id')
dat <- merge(dat, dom, by='id')


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

dat <- dat[
  ,.(ttr5_lai, 
     year,
     # fire_month_f,
     # vc_name_f,
     min_nbr_anom, 
     malai, 
     map,
     mapet,
     mavpd15,
     matmax,
     matmin,
     des,
     der,
     pH,
     silt,
     sand,
     clay,
     elevation,
     slope,
     aspect,
     pre_fire_slai_anom_12mo,
     vpd15_anom_3mo,
     precip_anom_12mo,
     post_precip_anom_12mo,
     post_vpd15_anom_12mo,
     post_tmax_anom_12mo)][
       is.na(ttr5_lai)==F &
         is.na(elevation)==F &
         is.na(pH)==F &
         is.na(des)==F &
         is.na(slope)==F &
         is.na(aspect)==F &
         is.na(pre_fire_slai_anom_12mo)==F] %>%
  as_tibble()

arrow::write_parquet(dat, 
                     sink="../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai.parquet")

################################################################################
# Set up hyp-param grid ---------------------------------------------------
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai.parquet")
# data cleaning
library(janitor)

# data prep
library(dplyr)
library(tidyverse)

# tidymodels
library(rsample)
library(recipes)
library(parsnip)
library(tune)
library(dials)
library(workflows)
library(yardstick)

# speed up computation with parrallel processing (optional)
# library(doParallel)
# all_cores <- parallel::detectCores(logical = FALSE)
# registerDoParallel(cores = all_cores)

# tidymodels_prefer()
set.seed(123)

d_split <- initial_split(dat, strata = ttr5_lai)
d_train <- training(d_split)
d_test <- testing(d_split)

# Specify CV fold data ------------------------
d_folds <- vfold_cv(d_train, strata = year, v=10)

# Specify RF mod specs ------------------
# xgb_spec <- boost_tree(mtry=4,
#                        trees = 100,
#                        min_n = tune(),
#                        # tree_depth = tune(),
#                        # learn_rate = tune(),
#                        # loss_reduction = tune(),
#                        # sample_size = tune()
#                        ) %>% 
#   set_mode(mode='regression') %>% 
#   set_engine("xgboost", nthread=12)

xgb_spec <- boost_tree(
  trees = tune(), 
  tree_depth = tune(), 
  min_n = tune(), 
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), 
  mtry = tune(),         ## randomness
  learn_rate = tune(),                         ## step size
  ) %>% 
  set_engine("xgboost") %>% 
  set_mode("regression")

# Specify 'Recipe' --------------------
xgb_rec <- recipe(ttr5_lai ~ ., data=d_train) #%>% 
  # update_role(year, new_role = 'ID')

# Specify 'Workflow' ------------------
xgb_wf <- workflow() %>%
  add_formula(ttr5_lai ~ .) %>% 
  # add_recipe(xgb_rec) %>% 
  add_model(xgb_spec) 
  
# test_fit <- xgb_wf %>% fit(data=d_train)
#   
# gc(full=T)
# 
# test_fit <- workflow() %>% 
#   add_recipe(xgb_rec) %>% 
#   add_model(boost_tree() %>% 
#               set_mode('regression') %>% 
#               update(mtry=3,
#                      min_n=3,
#                      tree_depth=4,
#                      trees=100) %>% 
#               set_engine("xgboost",nthread=12)) %>% 
#   fit(d_train)
# test_fit
# apply(d_train,2,class)
# apply(d_train,2,function(x) sum(is.infinite(x)))


# set up the grid of tuning parameters
xgb_grid <- grid_max_entropy(   finalize(mtry(), d_train),
                                trees(range=c(800,800)),
                                min_n(range=c(3,5)),
                                tree_depth(),
                                learn_rate(range = c(0.005,0.05)),
                                loss_reduction(range = c(1e-3,3e-3)),
                                sample_size = sample_prop(c(0.6,0.9)),
                                # finalize(sample_size(range=c(0,1)), d_train),
                                size=20)
xgb_grid

doParallel::registerDoParallel()
# doParallel::stopImplicitCluster()
set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = d_folds,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)
xgb_res
xgb_res$.notes[[1]]$.notes
xgb_res$.notes[[2]]$.notes
collect_metrics(xgb_res)


# Plot the model fits to the tuning parameters
xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry, tree_depth, min_n,learn_rate,loss_reduction,  sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")

xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry, tree_depth, min_n,learn_rate,loss_reduction,  sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "R2")

best_rmse <- select_best(xgb_res, metric = 'rmse')
best_rmse

final_fit <- finalize_model(xgb_spec, parameters = best_rmse)
final_fit

# Examine variable importance ------------------
library(vip)

test <- workflow() %>%
  add_formula(ttr5_lai ~ .) %>% 
  # add_recipe(xgb_rec) %>% 
  # add_model(xgb_spec) %>% 
  add_model(final_fit %>% 
              set_engine("xgboost", importance = "permutation")) %>% 
  fit(data=d_train)

test$fit$fit %>% vip(geom='point')


# Fit Final Model -------------------------
final_wf <- workflow() %>%
  add_model(final_fit) %>% 
  add_recipe(xgb_rec)


test_prediction <- final_wf %>%
  # fit the model on all the training data
  fit(
    data    = d_train
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = d_test) %>%
  bind_cols(testing(d_split))

# measure the accuracy of our model using `yardstick`
xgboost_score <- 
  test_prediction %>%
  yardstick::metrics(ttr5_lai, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))

knitr::kable(xgboost_score)

final_res <- last_fit(final_wf, d_split)

final_res %>%
  collect_metrics()

final_res %>% 
  collect_predictions() %>% 
  ggplot(data=.,aes(.pred, ttr5_lai))+
  geom_point()+
  geom_abline(col='red')

# MGCV comparison -----------
library(mgcv)
b1 <- bam(ttr5_lai~
            te(malai, min_nbr_anom,mapet, k=5, bs='cs')+
            te(map, post_precip_anom_12mo,post_vpd15_anom_12mo, k=5, bs='cs'),
          # te(matmax,matmin, k=5, bs='cs')+
          # s(mapet, k=5, bs='cs')+
          # s(mavpd15, k=5, bs='cs')+
          # s(post_vpd15_anom_12mo, k=5, bs='cs')+
          # s(elevation, k=5, bs='cs'),
          # te(map, 
          #    min_nbr_anom,
          #    post_precip_anom_12mo)+
          # s(matmax,matmin)+
          # s(mavpd15, post_vpd15_anom_12mo)+
          # te(elevation,mapet),
          data=d_train, 
          select=TRUE, 
          discrete=TRUE)
summary(b1)
plot(b1, scheme=2, pages=1)

yardstick::rsq_vec(d_test$ttr5_lai, 
                   estimate=predict(b1, newdata=d_test, type='response'))
