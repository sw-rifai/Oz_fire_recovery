library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(arrow)

# PART 1: PREP  ----------------------------------------------------------------
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


arrow::write_parquet(dat, "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai.parquet")

################################################################################
# PART 2: RANDOM FOREST  -----
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
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai.parquet")

# Set up hyp-param grid ---------------------------------------------------
set.seed(123)
d_split <- initial_split(dat, 
 strata = ttr5_lai, prop = 0.25,breaks = 30)
d_train <- training(d_split)
d_test <- testing(d_split)

# Specify RF mod specs ------------------
rf_spec <- rand_forest(mode='regression', 
                      mtry=tune(), 
                      trees=tune(),
                      min_n=tune()) %>% 
  set_engine("ranger",
             regularization.factor = tune("regularization"), 
             num.threads=4)

# Specify CV fold data ------------------------
d_folds <- vfold_cv(d_train, strata = year, v=6)

# Specify 'Recipe' --------------------
rf_rec <- recipe(ttr5_lai ~ ., data=d_train) %>% 
  # step_log(ttr5_lai) %>% 
  update_role(year, new_role = 'ID') %>% 
  prep(training = d_train)

# Specify 'Workflow' ------------------
rf_wf <- workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rf_rec)

gc(full=T)

# set up the grid of tuning parameters
rf_grid <- grid_latin_hypercube(finalize(mtry(range=c(2,8)), d_train),
                                trees(range=c(200,1600)),
                                min_n(range=c(2,7)),
                                regularization = regularization_factor(range = c(0.5,2)),
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
final_rf

# Examine variable importance ------------------
library(vip)

test <- workflow() %>% 
  add_recipe(., recipe = rf_rec) %>% 
  add_model(final_rf %>% 
  set_engine("ranger", 
             importance = "permutation", 
             num.threads=10)) %>% 
  fit(data=d_train)

test$fit$fit %>% vip(geom='point')


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
  ggplot(data=.,aes(ttr5_lai, .pred-ttr5_lai))+
  # ggplot(data=.,aes(exp(ttr5_lai), exp(.pred)-exp(ttr5_lai)))+
  geom_point()


# MGCV comparison -----------
library(mgcv)
tmp %>% ggplot(data=.,aes(min_nbr_anom,ttr5_lai))+
  geom_smooth()
dat <- dat %>%  
  mutate(y_trans = log(ttr5_lai-363)) %>% 
  mutate(w =  (abs(min_nbr_anom)))
dat$w %>% summary
dat$w %>% hist
tmp <- dat %>% sample_n(10000)
b1 <- bam(y_trans~
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
            data=tmp,
          weights = tmp$w,
          family=tw(),
          select=TRUE, 
          discrete=TRUE)
summary(b1)
plot(b1, scheme=2, pages=1)

yardstick::rsq_vec(d_test$ttr5_lai, 
                   estimate=predict(b1, newdata=d_test, type='response'))
gratia::simulate(b1)
simulate(b1) %>% as.numeric() %>% exp() %>% hist
dat$ttr5_lai %>% quantile()
simulate(b1) %>% as.numeric() %>% `^`(2) %>% `+`(363) %>% hist
gam.check(b1)
qq.gam(b1)

5 %>% `^`(2)

dat %>% sample_n(10000) %>% 
  mutate(y_trans = (ttr5_lai-363)) %>% pull(y_trans) %>% summary
dat %>% 
  sample_n(10000) %>% 
  mutate(pred = exp(predict(b1,type='response',newdata=.)) + 363) %>% 
  filter(pred < 20000) %>% 
  ggplot(data=.,aes(ttr5_lai, ttr5_lai-pred))+
  geom_point()+
  geom_smooth()+
  geom_hline(aes(yintercept=0),col='red')

# SCRATCH --------------------------------------------------------------------
test_results <- 
  d_test %>%
  select(ttr5_lai) %>%
  mutate(ttr5_lai = log(ttr5_lai)) %>%
  bind_cols(
    predict(rf_xy_fit, new_data = ames_test[, preds])
  )
test_normalized <- bake(rf_rec, new_data = d_test, all_predictors())
test_results <- 
  test_normalized %>%
  rename(`random forest` = .pred) %>%
  bind_cols(
    predict(final_res, new_data = test_normalized)
  )

d_train$ttr5_lai %>% hist

fn <- ecdf(d_train$ttr5_lai)
curve(fn(x),365,5000)
ws <- fn(d_train$ttr5_lai)
l1 <- lm(log(ttr5_lai)~scale(min_nbr_anom), 
         data=d_train, 
        weights = ws)
summary(l1)
d_train %>% 
  sample_n(100) %>% 
  mutate(pred = exp(as.double(predict(l1, newdata=., type='response')))) %>% 
  mutate(res = pred-ttr5_lai) %>% 
  ggplot(data=.,aes(ttr5_lai, res))+
  geom_point()


d_train %>% 
  sample_n(1000) %>% 
  mutate(pred = exp(as.double(predict(l1, newdata=., type='response')))) %>% 
  mutate(res = pred-ttr5_lai) %>% 
  lm(res~ttr5_lai, data=.) -> l2

d_train %>% 
  sample_n(1000) %>% 
  mutate(pred = exp(as.double(predict(l1, newdata=., type='response')))) %>% 
  mutate(res = pred-ttr5_lai) %>% 
  mutate(p_res = predict(l2)) %>% 
  mutate(pred2 = pred+p_res) %>% 
  ggplot(data=.,aes(ttr5_lai, pred2-ttr5_lai))+
  geom_point()+
  geom_hline(aes(yintercept=0), col='red')+
  geom_smooth(method='lm')+
  geom_smooth(method='lm',aes(ttr5_lai,pred-ttr5_lai),col='yellow')


fn_dn <- approxfun(density(d_train$ttr5_lai))
s <- case_when(d_train$ttr5_lai <= 366 ~ 0.75, 
               between(d_train$ttr5_lai,366.1,750) ~ 1, 
               between(d_train$ttr5_lai,750.1,1000) ~ 2, 
               between(d_train$ttr5_lai,1001,1201) ~ 3, 
               between(d_train$ttr5_lai,1201,1500) ~ 6, 
               between(d_train$ttr5_lai,1500,2500) ~ 18, 
               d_train$ttr5_lai > 2500 ~ 36)
sum(is.na(s))
hist(s)
curve( exp((x-365)/750), 365,4000)
s <- exp((d_train$ttr5_lai-365)/750)/100
r1 <- ranger::ranger(log(ttr5_lai) ~ 
                   min_nbr_anom+malai+
                   map+mapet+mavpd15+
                     matmax+matmin+
                   elevation+sand+pH+silt+
                   der+
                   pre_fire_slai_anom_12mo + 
                   post_precip_anom_12mo+
                   post_vpd15_anom_12mo+
                   precip_anom_12mo, 
                 data=dat %>% filter((year%%2)==0) %>% 
                   slice_sample(prop=0.5), 
                 # splitrule = 'maxstat',
                 # minprop = 0.5,
                 # alpha = 0.5,
                 replace = TRUE,
                 mtry=6,
                 regularization.factor = 1,
                 # treetype='regression',
                 num.trees = 1000,
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
  ggplot(data=.,aes(ttr5_lai, exp(pred)-ttr5_lai ))+
  ggpointdensity::geom_pointdensity()+
  geom_hline(col='red',aes(yintercept=0))+
  geom_smooth()+
  scale_color_viridis_c(option='H')

d_test %>% 
  mutate(pred = exp(predict(r1, data=.)$predictions)) %>% 
  select(pred, ttr5_lai, year) %>% 
  pivot_longer(cols=c("pred","ttr5_lai")) %>% 
  ggplot(data=.,aes(x=year, y=value, 
                    group=paste(name,year), 
                    fill=name))+
  geom_boxplot(outlier.colour = NA)+
  coord_cartesian(ylim=c(365,2500), expand=F)+
  scale_x_continuous(breaks=seq(2001,2015,by=2))+
  scale_fill_manual(values=c("pred"='lightblue',
                             "ttr5_lai"="#CF0000"), 
                    labels=c("pred"='Pred.',
                             "ttr5_lai"='Obs.'))+
  labs(x=NULL,
       y='LAI: Time to recover (days)', 
       fill=NULL, 
       caption = "Untuned Random Forest, trained on even years only")+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave(filename = "figures/boxplot_pred-randForest-untuned.png")

d_test %>% 
  sample_n(1000) %>% 
  mutate(pred = predict(r1, data=.)$predictions) %>% 
  ggplot(data=.,aes(exp(pred), (ttr5_lai)))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')

d_test %>% 
  sample_n(1000) %>% 
  mutate(pred = predict(r1, data=.)$predictions) %>% 
  ggplot(data=.,aes(ttr5_lai, exp(pred)))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')


d_test %>% 
  sample_n(1000) %>% 
  mutate(pred = predict(r1, data=.)$predictions) %>% 
  mutate(res = ttr5_lai-exp(pred)) %>% 
  ggplot(data=.,aes(exp(pred), res))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_hline(aes(yintercept=0),col='red')

d_test %>% 
  sample_n(1000) %>% 
  mutate(pred = predict(r1, data=.)$predictions) %>% 
  mutate(res = exp(pred) - ttr5_lai) %>% 
  ggplot(data=.,aes(ttr5_lai, res))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_hline(aes(yintercept=0),col='red')+
  labs(x='Obs.', 
       y='Residuals')

  


jj <- d_test %>% 
  sample_n(1000) 
jj$pred <- predict(r1, data=jj)$predictions
jj$res = exp(jj$pred) - jj$ttr5_lai 

mpe_vec(truth=jj$ttr5_lai, 
        estimate=exp(jj$pred))
mpe_vec(truth=jj$ttr5_lai, 
        estimate=exp(predict(l1,newdata=jj)))
mpe_vec(truth=jj$ttr5_lai, 
        estimate=runif(365,5000,n = dim(jj)[1]))


jj %>% ggplot(data=.,aes(ttr5_lai, res))+
  geom_point()+
  geom_smooth(method='lm')


jj %>% mutate(pred = predict(r1, data=.)$predictions) %>% 
  ggplot(data=.,aes(ttr5_lai, exp(pred),color=res))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')+
  geom_vline(aes(xintercept=1000),col='orange')+
  geom_hline(aes(yintercept=1000),col='orange')+
  scale_color_gradient2(mid = 'black',high='red',low='blue',
                        limits=c(-1000,1000),oob=scales::squish)+
  coord_equal(xlim=c(365,3000),
              ylim=c(365,3000))

jj %>% 
  ggplot(data=.,aes(ttr5_lai, res))+
  geom_point()+
  geom_hline(aes(yintercept=0), col='red')+
  geom_smooth(method='lm')

jj %>% 
  ggplot(data=.,aes(exp(pred), ttr5_lai))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')



library(xgboost)
x <- Matrix::sparse.model.matrix(~, -1, d_train %>% 
  select(min_nbr_anom,malai,
  map,mapet,
  pre_fire_slai_anom_12mo , 
  post_precip_anom_12mo,
  post_vpd15_anom_12mo) %>% as.matrix)
y <- d_train$ttr5_lai %>% as.matrix()



xd_train <- xgb.DMatrix(data = x, 
                        # label = y,
                        missing = NA)

# the tweedie_variance_power parameter determines the shape of
# distribution
# - closer to 1 is more poisson like and the mass
#   is more concentrated near zero
# - closer to 2 is more gamma like and the mass spreads to the
#   the right with less concentration near zero

params <- list(
  objective = 'reg:tweedie',
  eval_metric = 'rmse',
  tweedie_variance_power = 1.4,
  max_depth = 6,
  eta = 1)

bst <- xgb.train(
  data = d_train,
  params = params,
  maximize = FALSE,
  watchlist = list(train = d_train),
  nrounds = 20)

var_imp <- xgb.importance(attr(x, 'Dimnames')[[2]], model = bst)

preds <- predict(bst, d_train)

rmse <- sqrt(sum(mean((y - preds) ^ 2)))


plot(y=preds-y, 
     x=y)
abline(h=0,col='red')
