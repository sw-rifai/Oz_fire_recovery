library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(arrow)
library(xgboost)
library(tidymodels); 

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
post_clim <- post_clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
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


# Set up hyp-param grid ---------------------------------------------------
set.seed(123)
d_split <- initial_split(dat[
  ,.(ttr5_lai, 
     # year,
     fire_month_f, 
     vc_name_f, 
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
     post_tmax_anom_12mo)
 ][sample(.N,50000)], strata = fire_month_f)
d_train <- training(d_split)
d_test <- testing(d_split)

xgb_spec <- boost_tree(
  trees = tune(),
  tree_depth = tune(), 
  min_n = tune(), 
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), 
  mtry = tune(),         ## randomness
  learn_rate = tune()                         ## step size
   ) %>% 
  set_engine("xgboost") %>% 
  set_mode("regression")

xgb_spec

xgb_grid <- grid_latin_hypercube(
  trees(),
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), d_train),
  learn_rate(),
  size = 20
)

xgb_grid

xgb_wf <- workflow() %>%
  add_formula(ttr5_lai ~ .) %>%
  add_model(xgb_spec) #%>%
  # step_rm('year')
  # update_role(year, new_role = 'ID')
  
xgb_wf


set.seed(123)
d_folds <- vfold_cv(d_train, strata = fire_month_f)

d_folds
gc()


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

xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rsq")


show_best(xgb_res, "rmse")


best_rmse <- select_best(xgb_res, "rmse")
best_rmse

best_rsq <- select_best(xgb_res, "rsq")
best_rsq

final_xgb <- finalize_workflow(
  xgb_wf,
  best_rmse
)

final_xgb


library(vip)

final_xgb %>%
  fit(data = d_train) %>%
  pull_workflow_fit() %>%
  vip(geom = "point")


final_res <- last_fit(final_xgb, d_split)
collect_metrics(final_res)

# END First tuning stage

# 2nd Tuning --------------------------------------------------------------
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

xgb_spec2 <- boost_tree(
  trees = 1200, 
  tree_depth = c(3,4,12), 
  min_n = tune(), 
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = c(0.25,0.5), 
  mtry = tune(),         ## randomness
  learn_rate = 0.01,                         ## step size
) %>% 
  set_engine("xgboost") %>% 
  set_mode("regression")

xgb_spec2 <- xgb_spec %>% 
  update(trees=1200, 
         tree_depth=c(3,4,12),
         sample_size=c(0.25,0.5),
         learn_rate=0.1)

# xgb_grid2 <- grid_latin_hypercube(
#   trees(),
#   tree_depth(),
#   min_n(),
#   loss_reduction(),
#   sample_size = sample_prop(),
#   finalize(mtry(), d_train),
#   learn_rate(),
#   size = 20
# )

xgb_wf2 <- workflow() %>%
  add_formula(ttr5_lai ~ .) %>%
  add_model(xgb_spec2)

xgb_wf2


set.seed(123)
d_folds <- vfold_cv(d_train, strata = fire_month_f)

d_folds

doParallel::registerDoParallel()
# doParallel::stopImplicitCluster()
set.seed(234)
xgb_res2 <- tune_grid(
  xgb_wf2,
  resamples = d_folds,
  grid = xgb_grid2,
  control = control_grid(save_pred = TRUE)
)

xgb_res2

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

xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "rsq") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rsq")


show_best(xgb_res, "rmse")


final_xgb %>% class
pull_workflow_fit(final_xgb)
fit(final_xgb, 
    data=d_split)


o %>% class
xgb.importance(pull_workflow_fit(o))
