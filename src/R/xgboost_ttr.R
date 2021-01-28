library(tidyverse);
library(usethis);
library(stars);
library(data.table); 
library(dtplyr); 
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
source("src/R/functions_time_to_recover.R")
library(tidymodels)
set.seed(123)

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

p5 <- final_xgb %>%
  fit(data = vb_train) %>%
  pull_workflow_fit()
p5 %>%
  vip(geom = "col", 
      num_features = 20)+
  theme_linedraw()
ggsave(filename = 'figures/xgb_vip_timeToInflectionWeibull_1burn_2001-2015.png')
