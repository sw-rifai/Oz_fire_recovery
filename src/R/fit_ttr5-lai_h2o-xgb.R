pacman::p_load(tidyverse, 
               data.table,
               h2o,
               lubridate)

# Load, filter, and split data -------------------------
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


# H2O xgboost model fitting ----------------------------------------------------
h2o.init(nthreads = 22)
Sys.setenv(sys.ai.h2o.auc.maxClasses=81)
Sys.getenv("sys.ai.h2o.auc.maxClasses")

# Set predictors and response; set response as a factor
hdat <- as.h2o(d_rf)

colnames(d_rf)
response <- "ttr_ocat"
predictors <- setdiff(colnames(d_rf), response)
predictors <- predictors[str_detect(predictors,'fire_year')==F]

# Split the dataset into train and valid

splits <- h2o.splitFrame(data =  hdat, ratios = .75, seed = 1234)
train <- splits[[1]]
valid <- splits[[2]]




# auto ml attempt -------------------------------------------------------------
a1 <- h2o.automl(x = predictors, 
                 y = response,
                 training_frame = train, 
                 validation_frame = valid,
                 include_algos = c('XGBoost','GBM'),
                 fold_column = 'fire_year',
                 max_runtime_secs = 1800,
                 seed = 1234)
lb <- h2o.get_leaderboard(a1,extra_columns = 'ALL')
print(lb,n = nrow(lb))
head(lb)
top_mod <- h2o.get_best_model(a1,criterion = 'mean_per_class_error')
top_mod
h2o.varimp(top_mod)
vip::vip(top_mod, 25)
h2o.saveModel(top_mod,
              path=paste0("../data_general/proc_data_Oz_fire_recovery/",
                          "xgboost_ttr5-lai-ocat_",Sys.Date(),"_"))

top_mod@parameters$ntrees     # 69
top_mod@parameters$max_depth  # 20
top_mod@parameters$col_sample_rate #0.8
top_mod@parameters$reg_lambda # 10
top_mod@parameters$distribution
top_mod@parameters
top_mod@allparameters
h2o.performance(top_mod)@metrics
h2o.performance(top_mod,newdata = valid)


dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
dpreds <- bind_cols(h2o.predict(top_mod, as.h2o(dpreds)) %>% as.data.table(),
                    dpreds)

dpreds %>% 
  ggplot(data=.,aes(x,y,fill=predict))+
  geom_raster()+
  coord_equal()+
  scale_fill_viridis_d(option='B')

# grid search GBM -------------------------------------------------------------
grid1  <- h2o.grid("gbm", 
                 x = predictors, 
                 y = response, 
                 training_frame = train,
                 # validation_frame = valid,
                 # fold_column = 'fire_year',
                 hyper_params = list(
                   ntrees = floor(seq.int(20,100,length.out=10)),
                   max_depth = c(5)
                   ))
summary(grid1)

# GBM hyperparameters
gbm_params1 <- list(learn_rate = c(0.01, 0.1),
                    max_depth = c(3, 5, 9),
                    sample_rate = c(0.8, 1.0),
                    col_sample_rate = c(0.2, 0.5, 1.0))

# Train and validate a cartesian grid of GBMs
gbm_grid1 <- h2o.grid("gbm", x = predictors, y = response,
                      grid_id = "gbm_grid1",
                      training_frame = train,
                      validation_frame = valid,
                      ntrees = 100,
                      seed = 1,
                      hyper_params = gbm_params1)

# Get the grid results, sorted by validation AUC
gbm_gridperf1 <- h2o.getGrid(grid_id = "gbm_grid1",
                             sort_by = "mean_per_class_error",
                             decreasing = F)
print(gbm_gridperf1)

# Grab the top GBM model, chosen by validation AUC
best_gbm1 <- h2o.getModel(gbm_gridperf1@model_ids[[1]])

# Now let's evaluate the model performance on a test set
# so we get an honest estimate of top model performance
best_gbm_perf1 <- h2o.performance(model = best_gbm1,
                                  newdata = valid)
h2o.mean_per_class_error(best_gbm_perf1)

# Look at the hyperparameters for the best model
print(best_gbm1@model[["model_summary"]])
print(best_gbm1@model[["validation_metrics"]])
print(best_gbm1@model[["variable_importances"]])
best_gbm1@allparameters[["learn_rate"]]
best_gbm1@allparameters[["max_depth"]]
best_gbm1@allparameters[["sample_rate"]]
best_gbm1@allparameters[["col_sample_rate"]]


gbm_params2 <- list(learn_rate = c(0.09,0.1,0.2),
                    max_depth = c(9,10,11,12,13),
                    sample_rate = c(1.0),
                    col_sample_rate = c(1.0))
gbm_grid2 <- h2o.grid("gbm", x = predictors, y = response,
                      grid_id = "gbm_grid2",
                      training_frame = train,
                      validation_frame = valid,
                      ntrees = 100,
                      seed = 1,
                      hyper_params = gbm_params2)
gbm_gridperf2 <- h2o.getGrid(grid_id = "gbm_grid2",
                             sort_by = "mean_per_class_error",
                             decreasing = F)
print(gbm_gridperf2)
best_gbm2 <- h2o.getModel(gbm_gridperf2@model_ids[[1]])
best_gbm_perf2 <- h2o.performance(model = best_gbm2,
                                  newdata = valid)
h2o.mean_per_class_error(best_gbm_perf2)
print(best_gbm2@model[["model_summary"]])
print(best_gbm2@model[["validation_metrics"]])
print(best_gbm2@model[["variable_importances"]])
best_gbm2@allparameters[["learn_rate"]]
best_gbm2@allparameters[["max_depth"]]
best_gbm2@allparameters[["sample_rate"]]
best_gbm2@allparameters[["col_sample_rate"]]



gbm_params3 <- list(learn_rate = c(0.2,0.3,0.4,0.5),
                    max_depth = c(10,15,20),
                    sample_rate = c(1.0),
                    col_sample_rate = c(1.0))
gbm_grid3 <- h2o.grid("gbm", x = predictors, y = response,
                      grid_id = "gbm_grid3",
                      training_frame = train,
                      validation_frame = valid,
                      ntrees = 100,
                      seed = 1,
                      hyper_params = gbm_params3)
gbm_gridperf3 <- h2o.getGrid(grid_id = 'gbm_grid3',
                             sort_by = "mean_per_class_error",
                             decreasing = F)
print(gbm_gridperf3)
best_gbm3 <- h2o.getModel(gbm_gridperf3@model_ids[[1]])
best_gbm_perf3 <- h2o.performance(model = best_gbm3,
                                  newdata = valid)
h2o.mean_per_class_error(best_gbm_perf3)
print(best_gbm3@model[["model_summary"]])
print(best_gbm3@model[["validation_metrics"]])
print(best_gbm3@model[["variable_importances"]])
best_gbm3@allparameters[["learn_rate"]]
best_gbm3@allparameters[["max_depth"]]
best_gbm3@allparameters[["sample_rate"]]
best_gbm3@allparameters[["col_sample_rate"]]


gbm_params4 <- list(learn_rate = c(0.2),
                    ntrees=floor(seq(20,300,by=50)),
                    max_depth = c(20),
                    sample_rate = c(1.0),
                    col_sample_rate = c(1.0))
gbm_grid4 <- h2o.grid("gbm", x = predictors, y = response,
                      grid_id = "gbm_grid4",
                      training_frame = train,
                      validation_frame = valid,
                      seed = 1,
                      hyper_params = gbm_params4)
gbm_gridperf4 <- h2o.getGrid(grid_id = "gbm_grid4",
                             sort_by = "mean_per_class_error",
                             decreasing = F)
print(gbm_gridperf4)
best_gbm4 <- h2o.getModel(gbm_gridperf4@model_ids[[1]])
best_gbm_perf4 <- h2o.performance(model = best_gbm4,
                                  newdata = valid)
h2o.mean_per_class_error(best_gbm_perf4)
print(best_gbm4@model[["model_summary"]])
print(best_gbm4@model[["validation_metrics"]])
print(best_gbm4@model[["variable_importances"]])
best_gbm4@allparameters[["learn_rate"]]
best_gbm4@allparameters[["max_depth"]]
best_gbm4@allparameters[["sample_rate"]]
best_gbm4@allparameters[["col_sample_rate"]]

top_mod <- best_gbm4
h2o.saveModel(top_mod,
              path=paste0("../data_general/proc_data_Oz_fire_recovery/",
                          "gbm_ttr5-lai-ocat_",Sys.Date(),"_"))

dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
dpreds <- bind_cols(h2o.predict(top_mod, as.h2o(dpreds)) %>% as.data.table(),
         dpreds)
arrow::write_parquet(dpreds,
                   sink=paste0("../data_general/proc_data_Oz_fire_recovery/predicted_ttr5-lai-5class_gbm_",Sys.time(),"_.parquet"),
                   compression='snappy')





# SCRATCH ---------------------------------------------------------------------
# dpreds$predict %>% as.numeric() %>% hist
# dpreds %>%
#   ggplot(data=.,aes(x,y,fill=predict))+
#   # geom_sf(data=oz_poly,
#   #         fill='grey70',
#   #         color='grey30',
#   #         inherit.aes = F)+
#   geom_raster()+
#   coord_equal()+
#   scale_fill_viridis_d(option='B',end=0.95)+
#   coord_sf(xlim = c(146,153.25),
#            ylim = c(-38,-27.75))
# 
# dpreds %>% 
#   mutate(y_group = cut(y,c(-40,-35,-31,-27), 
#                        labels=c("south","mid","north"))) %>% 
#   ggplot(data=.,aes(y=predict,
#                     fill=y_group,
#                     x=elevation))+
#   geom_violin(outlier.colour = NA)+
#   facet_wrap(cut_number(min_nbr_anom,3)~y_group,scales='free')
# 
# 
# # Set predictors and response; set response as a factor
# colnames(d_rf)
# response <- "ttr_ocat"
# predictors <- setdiff(colnames(d_rf), response)
# predictors <- predictors[str_detect(predictors,'fire_year')==F]
# 
# hdat2 <- as.h2o(d_rf %>% filter(fire_year!=2013))
# splits2 <- h2o.splitFrame(data =  hdat2, ratios = .75, seed = 1234)
# train2 <- splits2[[1]]
# valid2 <- splits2[[2]]
# 
# rtop <- h2o.gbm(x = predictors, 
#                  y = response,
#                  training_frame = train2, 
#                  validation_frame = valid2,
#                  fold_column = 'fire_year',
#                 ntrees = 100,#81
#                 max_depth=15,#7
#                 min_rows=5,
#                 nbins = 20,
#                 nbins_top_level = 1024,
#                 nbins_cats = 1024,
#                 stopping_metric = 'deviance',
#                 seed=10999,
#                 learn_rate=0.1,
#                 learn_rate_annealing = 1,
#                 distribution='multinomial',
#                 huber_alpha = 0.9,
#                 sample_rate=1,
#                 col_sample_rate = 0.4,
#                 col_sample_rate_per_tree = 0.4)
# 
# h2o.performance(rtop,newdata = valid2)
# 
# 
# vec_tru <- d_rf %>% filter(fire_year==2013) %>% pull(ttr_ocat)
# vec_pred <- h2o.predict(rtop,newdata=as.h2o(d_rf %>% 
#                                   filter(fire_year==2013)))
# vec_pred <- vec_pred %>% as.data.table
# 
# 
# data.table(pred=vec_pred$predict, 
#            obs=vec_tru) %>% 
#   ggplot(data=.,aes(as.numeric(pred),as.numeric(obs)))+
#   geom_point(position='jitter',
#              size=0.1)+
#   scale_color_viridis_c()+
#   # geom_abline(col='red')+
#   geom_smooth(method='lm')
# 
# 
# 
# 
# 
# 
