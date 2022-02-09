pacman::p_load(tidyverse, 
               data.table,
               h2o,
               lubridate, 
  stars) 

# Load, filter, and split data -------------------------
d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")
d_rf <- d_rf %>% filter(date_fire1 <= lubridate::ymd("2017-03-01")) %>% 
  mutate(fire_year = lubridate::year(date_fire1-months(3))) %>% 
  mutate(fire_year = factor(fire_year))
d_rf <- d_rf %>% mutate(ttr_ocat = factor(ttr_ocat,levels = 1:5, labels = 1:5, ordered = F)) 
d_rf <- d_rf %>% 
  rename(x=x.x,y=y.x) %>% 
  # mutate(x = round(x.x*10)/10,
  # y=round(y.x*10)/10) %>% 
  as.data.table()

out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-14 16:45:26.parquet")
sout <- set_names(st_as_stars(out[,.(x,y,predict)], dims = c("x","y"),crs=4326),c("species"))
st_crs(sout) <- st_crs(4326)
vv <- st_extract(sout,
  st_as_sf(d_rf[,.(x,y)],coords=c("x","y"),crs=4326))
d_rf <- bind_cols(d_rf, st_drop_geometry(vv))

d_rf <- d_rf %>% select(
  ttr5_lai, 
  fire_year,
  species,
    x,y,
    malai, 
    fire_month, 
    min_nbr_anom, 
    pre_fire_slai_anom_3mo, 
    pre_fire_slai_anom_12mo, 
    
    elevation, 
    # x, y,
    # slope, 
    # aspect,
    # der,
    # des,
    sand,silt,clay,
    # pH,
    # bdw,
    # awc,
    # nto,ece,pto,
    
    map,mapet,mavpd15,mappet,matmin,matmax,
    precip_anom_12mo,vpd15_anom_3mo,
    post_precip_anom_12mo,
    post_vpd15_anom_12mo)

# drop the NA's, necessary for xgboost
d_rf <- d_rf %>% drop_na()

tmp_cv <- d_rf %>% 
  group_by(fire_year) %>% 
  summarize(nobs=n()) %>% 
  ungroup() %>% 
  arrange(nobs) %>% 
  mutate(id = 1:nrow(.)) %>% 
  mutate(cv_group=id%%4) %>% 
  mutate(cv_group = factor(cv_group)) %>% 
  select(fire_year,cv_group)
d_rf <- inner_join(d_rf,tmp_cv,by='fire_year') %>% 
  select(-fire_year)

set.seed(1)
d_out <- d_rf %>% sample_n(10000)
d_rf <- anti_join(d_rf,d_out)

# Prep pred data ----------------------------------------------------------
dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
vv <- st_extract(sout,
  st_as_sf(dpreds[,.(x,y)],coords=c("x","y"),crs=4326))
dpreds <- bind_cols(dpreds, st_drop_geometry(vv))





# Setup H2O cluster and data ----------------------------------------------------
h2o.init(nthreads = 20, max_mem_size = "32G")

# Set predictors and response; set response as a factor
hdat <- as.h2o(d_rf)

colnames(d_rf)
response <- "ttr5_lai"
predictors <- setdiff(colnames(d_rf %>% select(-x,-y)), response)
predictors <- predictors[str_detect(predictors,'cv_group')==F]

# Split the dataset into train and valid

splits <- h2o.splitFrame(data =  hdat, ratios = .75, seed = 1234)
train <- splits[[1]]
valid <- splits[[2]]

# Fit single models --------------------------------------------------
## ntrees = 70, max_depth = 50, min_rows = 10
s0 <- h2o.gbm(x = predictors[!str_detect(predictors,"species")],
        y = response,
        ntrees = 70,
        max_depth = 50,
        min_rows = 10,
        fold_column = 'cv_group',
        seed = 1111,
        keep_cross_validation_predictions = F,
        training_frame = train)
s0.1 <- h2o.gbm(x = predictors[!str_detect(predictors,"species")],
        y = response,
        ntrees = 70,
        max_depth = 50,
        min_rows = 10,
        # fold_column = 'cv_group',
        seed = 1111,
        keep_cross_validation_predictions = F,
        training_frame = train, 
        validation_frame = valid)
h2o.performance(s0,newdata = train)
h2o.performance(s0,newdata = valid)

h2o.performance(s0,newdata = valid)
h2o.performance(s0.1,newdata = valid)

h2o.performance(s0,newdata = as.h2o(d_out))
h2o.performance(s0.1,newdata = as.h2o(d_out))




# With species
# ntrees = 100, max_depth=30, min_rows=5
s1 <- h2o.gbm(x = c("species","min_nbr_anom","elevation",
       "pre_fire_slai_anom_12mo",
       "post_precip_anom_12mo",
       "vpd15_anom_3mo",
       "post_vpd15_anom_12mo", 
       "malai","matmin","matmax",
       "mavpd15",
  "map","mapet", # keep
     # "mappet",
       "sand",
  "silt",
  "clay",
       "pre_fire_slai_anom_3mo", 
       "precip_anom_12mo", 
       "fire_month"),
        y = response,
        ntrees = 100,
        max_depth = 30,
        min_rows=5,
        fold_column = 'cv_group',
        seed = 1111,
        keep_cross_validation_predictions = F,
        training_frame = train)
s1.1 <- h2o.gbm(x = c("species","min_nbr_anom","elevation",
       "pre_fire_slai_anom_12mo",
       "post_precip_anom_12mo",
       "vpd15_anom_3mo",
       "post_vpd15_anom_12mo", 
       "malai","matmin","matmax",
       "mavpd15",
  "map","mapet", # keep
     # "mappet",
       "sand",
  "silt",
  "clay",
       "pre_fire_slai_anom_3mo", 
       "precip_anom_12mo", 
       "fire_month"),
        y = response,
        ntrees = 100,
        max_depth = 30,
        min_rows=5,
        # fold_column = 'cv_group',
        seed = 1111,
        keep_cross_validation_predictions = F,
        training_frame = train, 
        validation_frame = valid)
h2o.performance(s1,newdata = as.h2o(d_out))#@metrics$RMSE # beat 342
h2o.performance(s1.1,newdata = as.h2o(d_out))#@metrics$r2


h2o.performance(s0,newdata = as.h2o(d_out))#@metrics$RMSE # beat 342
h2o.performance(s1,newdata = as.h2o(d_out))#@metrics$r2
h2o.performance(s0,newdata = valid)#@metrics$RMSE # beat 342
h2o.performance(s1,newdata = valid)#@metrics$r2


h2o.saveModel(s0,path=paste0("outputs/GBM_bestModNoSpecies_",Sys.time()))
h2o.saveModel(s1,path=paste0("outputs/GBM_bestModWithSpecies_",Sys.time()))


# Plot Variable Importance Comparison -------------------------------------
h2o.saveModel(s0,path=paste0("outputs/GBM_bestModNoSpecies_",Sys.time()))
h2o.saveModel(s1,path=paste0("outputs/GBM_bestModWithSpecies_",Sys.time()))


mets_s0 <- h2o.performance(s0,newdata = valid)
mets_s1 <- h2o.performance(s1,newdata = valid)

p_s0 <- h2o.varimp(s0) %>% 
  as_tibble() %>% 
  arrange(-scaled_importance) %>% 
  first(10) %>% 
  mutate(ord = row_number()) %>% 
  mutate(f = fct_reorder(variable,ord)) %>% 
  ggplot(data=.,aes(y=f,
                    x=scaled_importance))+
  geom_col(fill='grey70')+
  scale_y_discrete(limits=rev, 
                   labels=c("post_precip_anom_12mo"='Post fire 12-mo Precip. anom.',
                            "vpd15_anom_3mo"='Pre fire 3-mo VPD anom.',
                            "min_nbr_anom"='NBR anomaly', 
                            "elevation"="Elevation",
                            "malai" = "Mean Annual LAI",
                            "precip_anom_12mo"="Pre fire 12-month Precip. anom.",
                            "post_vpd15_anom_12mo"="Post fire 12-month VPD anom.",
                            "matmin"="Mean Annual Tmin",
                            "mapet"="Mean Annual PET",
                       "sand"="Sand",
                       "pre_fire_slai_anom_3mo"="Pre fire 3-month LAI anom.",
                       "pre_fire_slai_anom_12mo"="Pre fire 12-month LAI anom.",
                       "matmax"='Mean Annual Tmax'))+
  geom_rect(aes(xmin=0.5,xmax=1,ymin=0,ymax=5), 
    fill='grey90')+
  geom_text(aes(0.7, 4.2),
    label=paste0("R²: ",format(mets_s0@metrics$r2,digits=2)), 
    size=9,
    color='black')+
  geom_text(aes(0.73, 2.8),
    label=paste0("RMSE: ",format(mets_s0@metrics$RMSE,digits=3)), 
    size=9,
    color='black')+
  geom_text(aes(0.735, 1.5),
    label=paste0("MAE: ",format(mets_s0@metrics$mae,digits=3)), 
    size=9,
    color='black')+
  labs(y=NULL,
       x='Scaled Variable Importance',
        title='Best Model Without Species')+
  scale_x_continuous(expand=c(0,0),limits = c(-0.06,1.05))+
  theme_linedraw()+
  theme(text = element_text(size=16),
    # panel.grid.minor = element_blank(), 
    #     panel.grid.major = element_blank(),
    #     axis.ticks.length.y = unit(0, 'mm'),
    #     axis.text.y.left = element_text(hjust=0,
    #                                     size=14,
    #                                     colour = 'black',
    #                                     face='bold',
    #                         margin = margin(t = 0,
    #                                         r = -250,
    #                                         b = 0,
    #                                         l = 0)), 
    #     plot.margin = margin(l=30,r=10)
    ); p_s0

p_s1 <- h2o.varimp(s1) %>% 
  as_tibble() %>% 
  arrange(-scaled_importance) %>% 
  first(10) %>% 
  mutate(ord = row_number()) %>% 
  mutate(f = fct_reorder(variable,ord)) %>% 
  ggplot(data=.,aes(y=f,
                    x=scaled_importance))+
  geom_col(fill='grey70')+
  scale_y_discrete(limits=rev, 
                   labels=c("post_precip_anom_12mo"='Post fire 12-mo Precip. anom.',
                            "vpd15_anom_3mo"='Pre fire 3-mo VPD anom.',
                            "min_nbr_anom"='NBR anomaly', 
                            "elevation"="Elevation",
                            "malai" = "Mean Annual LAI",
                            "precip_anom_12mo"="Pre fire 12-month Precip. anom.",
                            "post_vpd15_anom_12mo"="Post fire 12-month VPD anom.",
                            "matmin"="Mean Annual Tmin",
                            "mapet"="Mean Annual PET",
                       "sand"="Sand",
                       "species"="Species",
                       "clay"="Clay",
                       "pre_fire_slai_anom_3mo"="Pre fire 3-month LAI anom.",
                       "pre_fire_slai_anom_12mo"="Pre fire 12-month LAI anom.",
                       "matmax"='Mean Annual Tmax'))+
  geom_rect(aes(xmin=0.5,xmax=1,ymin=0,ymax=5), 
    fill='grey90')+
  geom_text(aes(0.7, 4.2),
    label=paste0("R²: ",format(mets_s1@metrics$r2,digits=2)), 
    size=9,
    color='black')+
  geom_text(aes(0.73, 2.8),
    label=paste0("RMSE: ",format(mets_s1@metrics$RMSE,digits=3)), 
    size=9,
    color='black')+
  geom_text(aes(0.735, 1.5),
    label=paste0("MAE: ",format(mets_s1@metrics$mae,digits=3)), 
    size=9,
    color='black')+
  labs(y=NULL,
       x='Scaled Variable Importance',
        title='Best Model Without Species')+
  labs(y=NULL,
       x='Scaled Variable Importance',
      title='Best Model Including Species')+
  scale_x_continuous(expand=c(0,0),limits = c(-0.06,1.05))+
  theme_linedraw()+
  theme(text = element_text(size=16),
    # panel.grid.minor = element_blank(), 
    #     panel.grid.major = element_blank(),
    #     axis.ticks.length.y = unit(0, 'mm'),
    #     axis.text.y.left = element_text(hjust=0,
    #                                     size=14,
    #                                     colour = 'black',
    #                                     face='bold',
    #                         margin = margin(t = 0,
    #                                         r = -250,
    #                                         b = 0,
    #                                         l = 0)), 
    #     plot.margin = margin(l=30,r=10)
    )
p_s1

# dev.new(width=12*2,height=15,units='cm')
ggsave(p_s0|p_s1, 
  filename = paste0("figures/fig-SX_compare_var-imp_",Sys.Date(),".png"),
  width=40,
  height=20,
  units='cm',
  dpi=300)
# END section ****************************************************************

# eval block ***
dp1 <- bind_cols(h2o.predict(s0, as.h2o(d_out)) %>% as.data.table(),
                    d_out)
dp4 <- bind_cols(h2o.predict(s1.1, as.h2o(d_out)) %>% as.data.table(),
                    d_out)
merge(
dp1 %>% select(x,y,predict) %>% rename(p1 = predict) %>% mutate(idx=1:nrow(dp1)),
dp4 %>% select(x,y,predict) %>% rename(p4 = predict) %>% mutate(idx=1:nrow(dp1)),
  by=c("x","y")) %>% 
  mutate(diff = p1-p4) %>% #pull(diff) %>% hist(breaks=100)
  ggplot(data=.,aes(x,y,
    # fill=diff
    color=diff
    ))+
  geom_point(size=0.25)+
  # geom_tile()+
  # scico::scale_fill_scico(palette='roma',midpoint=0)+
  scico::scale_color_scico(palette='roma',midpoint=0, 
    limits=c(-300,300),oob=scales::squish)+
  coord_sf()+
  theme_dark()

# eval block on BS preds ***
dp1 <- bind_cols(h2o.predict(s0, as.h2o(dpreds)) %>% as.data.table(),
                    dpreds)
dp2 <- bind_cols(h2o.predict(s1.1, as.h2o(dpreds)) %>% as.data.table(),
                    dpreds)
merge(
dp1 %>% select(x,y,predict) %>% rename(p1 = predict) %>% mutate(idx=1:nrow(dp1)),
dp2 %>% select(x,y,predict) %>% rename(p2 = predict) %>% mutate(idx=1:nrow(dp1)),
  by=c("x","y")) %>% 
  mutate(diff = p1-p2) %>% #pull(diff) %>% hist(breaks=100)
  ggplot(data=.,aes(x,y,
    fill=diff
    # color=diff
    ))+
  # geom_point(size=0.25)+
  geom_tile()+
  # scico::scale_fill_scico(palette='roma',midpoint=0)+
  scico::scale_fill_scico(palette='roma',midpoint=0, 
    limits=c(-500,500),oob=scales::squish)+
  coord_sf()+
  theme_dark()




# s1 <- h2o.gbm(x = predictors[!str_detect(predictors,"species")],
#         y = response,
#         ntrees = 55,
#         max_depth = 7,
#         fold_column = 'cv_group',
#         seed = 1111,
#         keep_cross_validation_predictions = TRUE,
#         training_frame = train)
# # 315 min
# yardstick::rmse_vec(
#   truth=d_out$ttr5_lai,
#   estimate=as.data.table(h2o.predict(s1,newdata = as.h2o(d_out)))$predict)
# 
# h2o.performance(s1,newdata = valid) # 
# h2o.performance(s1,newdata = valid)@metrics$r2
# h2o.performance(s1,newdata = train)@metrics$RMSE
# h2o.performance(s1,newdata = valid)@metrics$RMSE
# h2o.performance(s1,newdata = as.h2o(d_out))@metrics$RMSE
# s1@allparameters
# vip::vip(s1,10)
# 
# # eval block ***
# dp1 <- bind_cols(h2o.predict(s0, as.h2o(d_out)) %>% as.data.table(),
#                     d_out)
# dp4 <- bind_cols(h2o.predict(s4, as.h2o(d_out)) %>% as.data.table(),
#                     d_out)
# merge(
# dp1 %>% select(x,y,predict) %>% rename(p1 = predict) %>% mutate(idx=1:nrow(dp1)),
# dp4 %>% select(x,y,predict) %>% rename(p4 = predict) %>% mutate(idx=1:nrow(dp1)),
#   by=c("idx")) %>% 
#   mutate(diff = p1-p4) %>% 
#   ggplot(data=.,aes(x.x,y.y,fill=diff,color=diff))+
#   geom_point()+
#   scico::scale_fill_scico(palette='roma',midpoint=0)+
#   scico::scale_color_scico(palette='roma',midpoint=0)+
#   coord_sf()+
#   theme_dark()
# 
# 
# 
# s3 <- h2o.gbm(x = predictors,
#         y = response,
#         ntrees = 555,
#         max_depth = 7,
#         fold_column = 'cv_group',
#         seed = 1111,
#         keep_cross_validation_predictions = TRUE,
#         training_frame = train)
# h2o.performance(s3,newdata = valid) 
# h2o.performance(s3,newdata = valid)@metrics$r2
# h2o.performance(s3,newdata = valid)@metrics$RMSE
# h2o.performance(s3,newdata = as.h2o(d_out))@metrics$RMSE
# h2o.performance(s3,newdata = as.h2o(d_out))@metrics$r2
# s1@allparameters
# vip::vip(s3,15)+vip::vip(s1,15)
# 
# # parameter set for models with species: ------------------
# # ntrees = 100, max_depth=30, min_rows=5
# s4 <- h2o.gbm(x = c("species","min_nbr_anom","elevation",
#        "pre_fire_slai_anom_12mo",
#        "post_precip_anom_12mo",
#        "vpd15_anom_3mo",
#        "post_vpd15_anom_12mo", 
#        "malai","matmin","matmax",
#        "mavpd15",
#   "map","mapet", # keep
#      # "mappet",
#        "sand",
#   "silt",
#   "clay",
#        "pre_fire_slai_anom_3mo", 
#        "precip_anom_12mo", 
#        "fire_month"),
#         y = response,
#         ntrees = 100,
#         max_depth = 30,
#         min_rows=5,
#         fold_column = 'cv_group',
#         seed = 1111,
#         keep_cross_validation_predictions = TRUE,
#         training_frame = train)
# h2o.performance(s4,newdata = as.h2o(d_out))#@metrics$RMSE # beat 342
# h2o.performance(s4,newdata = as.h2o(d_out))#@metrics$r2
# vip::vip(s4,10)
# predictors
# 
# h2o.varimp_plot(s4,num_of_features = 25)
# 
# dp1 <- bind_cols(h2o.predict(s1, as.h2o(dpreds)) %>% as.data.table(),
#                     dpreds)
# dp4 <- bind_cols(h2o.predict(s4, as.h2o(dpreds)) %>% as.data.table(),
#                     dpreds)
# 
# merge(
# dp1 %>% select(x,y,predict) %>% rename(p1 = predict),
# dp4 %>% select(x,y,predict) %>% rename(p4 = predict),
#   by=c("x","y")) %>% 
#   mutate(diff = p1-p4) %>% 
#   ggplot(data=.,aes(x,y,fill=diff))+
#   geom_tile()+
#   scico::scale_fill_scico(palette='roma',midpoint=0)+
#   coord_sf()+
#   theme_dark()
# 
# tmp <- bind_cols(as.data.table(h2o.predict(s1,newdata = as.h2o(d_out))),d_out)
# 
# tmp %>% 
#   group_by(species) %>% 
#   summarize(r2 = cor(predict,ttr5_lai)**2, 
#             nobs = n()) %>%
#   ungroup() %>% 
#   merge(., tmp, by='species') %>% 
#   ggplot(data=.,aes(x,y,color=r2))+
#   geom_point()+
#   scale_color_viridis_c(option='H',direction = -1)+
#   coord_sf()
# 
# 
# 
# # S2 ----------------------------------------------------------------------
# s2 <- h2o.gbm(x = predictors,
#   y = response,
#   ntrees = 55,  #55
#   max_depth = 7, # 7
#   tweedie_power = 1.1,
#   min_rows=10,
#   balance_classes = T,
#   # monotone_constraints = list("min_nbr_anom"=-1, 
#   #   "elevation"=1,
#   #   "mappet"=1,
#   #   "pre_fire_slai_anom_12mo"=1,
#   #   "sand"=1),
#   distribution = 'tweedie',
#   fold_column = 'cv_group',
#   seed = 1111,
#   keep_cross_validation_predictions = F,
#   training_frame = train)
# vip::vip(s2)
# s2@allparameters$max_depth
# h2o.performance(s2,newdata = valid) # 347
# h2o.performance(s2,newdata = valid)@metrics$r2
# h2o.performance(s2,newdata = valid)@metrics$RMSE
# h2o.predict(s2,newdata=valid) %>%
#   as.data.table() %>% 
#   bind_cols(., as.data.table(valid)) %>% 
#   group_by(cv_group) %>% 
#   summarize(
#     r2 = yardstick::rmse_vec(truth=ttr5_lai,
#     estimate=predict), 
#     r2 = yardstick::rsq_trad_vec(truth=ttr5_lai,
#     estimate=predict))
# 
# 
# h2o.predict(s2,newdata=valid) %>% 
#   as.data.table() %>% 
#   pull(predict) %>% hist()
# dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
# dpreds <- bind_cols(h2o.predict(s2, as.h2o(dpreds)) %>% as.data.table(),
#                     dpreds)
# dpreds %>%
#   ggplot(data=.,aes(x,y,fill=predict))+
#   geom_raster()+
#   coord_equal()+
#   scale_fill_viridis_b(
#     option='B',
#     limits=c(365,2000),
#     oob=scales::squish,
#     n.breaks=6)+
#   # scico::scale_fill_scico(palette = 'batlowK',
#   #   direction = 1,
#   #   limits=c(365,2000),oob=scales::squish)+
#   theme_void()+
#   theme(panel.background = element_rect(fill='grey'))
# library(mgcv)
# dpreds %>% 
#   sample_n(10000) %>% 
#   select(contains(predictors),predict) %>% 
#   pivot_longer(cols = -predict) %>% 
#   ggplot(data=.,aes(value,predict))+
#   # geom_point()
#   geom_smooth(method='bam',se=F, 
#     formula=y~s(x,bs='cs'),
#     method.args=list(discrete=T))+
#   geom_rug()+
#   facet_wrap(~name,scales='free_x')
#   # scale_color_viridis_d(option='H')+
#   # facet_wrap(~cut(min_nbr_anom,c(-Inf,-0.75,-0.5,-0.25,Inf)))
# predictors
# dpreds[sample(.N,10000)] %>% 
#   ggplot(data=.,aes(min_nbr_anom,predict, 
#     color=cut(malai,c(0,1,2,3,4,5,Inf))))+
#   # geom_point()+
#   geom_smooth(method='gam',se=F)+
#   scale_color_viridis_d(option='H')+
#   facet_wrap(~cut(min_nbr_anom,c(-Inf,-0.75,-0.5,-0.25,Inf)))
# 
# 
# # Grid search madness ... -------------------------------------------------
# ### GBM grid search with species -----------------------------------------
# 
# grid0 <-
#   h2o.grid(algorithm = 'gbm',
#     # x = predictors,
#   x = c("species","min_nbr_anom","elevation",
#        "pre_fire_slai_anom_12mo",
#        "post_precip_anom_12mo",
#        "vpd15_anom_3mo",
#        "post_vpd15_anom_12mo", 
#        "malai","matmin","matmax",
#        "mavpd15",
#   "map","mapet", # keep
#      # "mappet",
#        "sand",
#   "silt",
#   "clay",
#        "pre_fire_slai_anom_3mo", 
#        "precip_anom_12mo", 
#        "fire_month"),
#   y = response,
#   # grid_id = 'grid_1',
#   hyper_params=list(
#   # ntrees = c(25,35,45,55,100,150,200,250,500,1200,2400,5000,10000),
#   # ntrees = c(350, 450, 500),
#   ntrees = c(100)
#   # max_depth = c(7,10,15,20,25,30,35),
#   # max_depth = floor(seq(20,25,length.out=3)),
#   # min_rows=c(5)
#   ),
#   max_depth = 30,
#   min_rows = 5,
#   tweedie_power = 1.1,
#   balance_classes = T,
#   # monotone_constraints = list("min_nbr_anom"=-1,
#   #   "elevation"=1,
#   #   "mappet"=1,
#   #   "pre_fire_slai_anom_12mo"=1,
#   #   "sand"=1),
#   distribution = 'tweedie',
#   # fold_column = 'cv_group',
#   seed = 1111,
#   keep_cross_validation_predictions = T,
#   training_frame = train,
#   validation_frame = valid)
# grid0
# 
# gm1 <- h2o.getModel(grid0@model_ids[[1]])
# 
# bind_cols(as.data.table(train), as.data.table(
#   h2o.getFrame(gm1@model[["cross_validation_holdout_predictions_frame_id"]][["name"]]))) %>% 
#   group_by(cv_group) %>% 
#   summarize(rmse = yardstick::rmse_vec(ttr5_lai, predict))
# 
# 
# gbm_gridperf0 <- h2o.getGrid(grid0@grid_id,
#                              sort_by = "mae",
#                              decreasing = F)
# gbm_gridperf0@summary_table %>% 
#   as.data.table() %>% 
#   ggplot(data=.,aes(ntrees,mae,color=factor(ntrees)))+
#   geom_point()
# # 
# # h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[1]]),newdata=valid)
# i <- 3
# h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[i]]),newdata=as.h2o(d_out))
# h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[i]]),newdata=valid)
# vip::vip(h2o.getModel(gbm_gridperf0@model_ids[[1]]))
# 
# 
# ## Grid search, no species -----------------------------
# grid1 <-
#   h2o.grid(algorithm = 'gbm',
#     x = predictors[!str_detect(predictors,'species')],
#   y = response,
#   hyper_params=list(
#   # ntrees = c(25,35,45,55,100,150,200,250,500,1200,2400,5000,10000),
#   ntrees = floor(seq.int(70)),
#   # ntrees = 250,
#   max_depth = c(50),
#   # max_depth = seq(10,100,by=5),
#   tweedie_power = 1.1,
#   min_rows=c(10)),
#   balance_classes = T,
#   # monotone_constraints = list("min_nbr_anom"=-1,
#   #   "elevation"=1,
#   #   "mappet"=1,
#   #   "pre_fire_slai_anom_12mo"=1,
#   #   "sand"=1),
#   distribution = 'tweedie',
#   # fold_column = 'cv_group',
#   seed = 1111,
#   keep_cross_validation_predictions = T,
#   training_frame = train,
#   validation_frame = valid)
# gbm_gridperf1 <- h2o.getGrid(grid1@grid_id,
#                              sort_by = "mae",
#                              decreasing = F)
# gbm_gridperf1@summary_table %>% 
#   as.data.table() %>% 
#   ggplot(data=.,aes(ntrees,mae,color=factor(max_depth)))+
#   geom_point()
# # 
# # h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[1]]),newdata=valid)
# i <- 1
# h2o.performance(h2o.getModel(gbm_gridperf1@model_ids[[i]]),newdata=as.h2o(d_out))
# h2o.performance(h2o.getModel(gbm_gridperf1@model_ids[[i]]),newdata=valid)
# h2o.performance(h2o.getModel(gbm_gridperf1@model_ids[[5]]),newdata=as.h2o(d_out))
# vip::vip(h2o.getModel(gbm_gridperf1@model_ids[[1]]))






# # RMSE 315 (no species)
# 
# vip::vip(h2o.getModel(gbm_gridperf0@model_ids[[1]])) + 
# vip::vip(h2o.getModel(gbm_gridperf0@model_ids[[2]]))
# 
# h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[1]]),
#   newdata=as.h2o(d_out))@metrics
# 
# h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[1]]),newdata=valid)
# 
# 
# h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[2]]),newdata=valid)
# h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[3]]),newdata=valid)
# h2o.performance(h2o.getModel(gbm_gridperf0@model_ids[[4]]),newdata=valid)
# 
# 
# grid1  <- h2o.grid("gbm", 
#                  x = predictors, 
#                  y = response, 
#             grid_id = 'grid1',
#                  training_frame = train,
#                  validation_frame = valid,
#                  fold_column = 'cv_group',
#                  hyper_params = list(
#                    ntrees = floor(seq.int(50,60,length.out=3)),
#                    max_depth = c(7)
#                    ))
# summary(grid1)
# gbm_gridperf1 <- h2o.getGrid(grid_id = grid1@grid_id,
#                              sort_by = "rmse",
#                              decreasing = F)
# print(gbm_gridperf1)
# best_gbm1 <- h2o.getModel(gbm_gridperf1@model_ids[[1]])
# best_gbm1
# 
# best_gbm1@allparameters
# h2o.performance(best_gbm1,newdata = valid)@metrics$r2
# # h2o.performance(top_mod,newdata = valid)@metrics$r2
# h2o.performance(best_gbm1,newdata = valid)@metrics$RMSE
# # h2o.performance(top_mod,newdata = valid)@metrics$RMSE
# vip::vip(best_gbm1,20)
# 
# dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
# dpreds <- bind_cols(h2o.predict(best_gbm1, as.h2o(dpreds)) %>% as.data.table(),
#                     dpreds)
# dpreds %>% 
#   ggplot(data=.,aes(x,y,fill=predict))+
#   geom_raster()+
#   coord_equal()+
#   scale_fill_viridis_b(
#     option='B',
#     limits=c(365,2000),
#     oob=scales::squish,
#     n.breaks=6)+
#   # scico::scale_fill_scico(palette = 'batlowK',
#   #   direction = 1,
#   #   limits=c(365,2000),oob=scales::squish)+
#   theme_void()+
#   theme(panel.background = element_rect(fill='grey'))
# 
# dpreds %>% 
#   sample_n(10000) %>% 
#   ggplot(data=.,aes(pre_fire_slai_anom_3mo/malai,predict))+
#   geom_smooth()+
#   geom_rug()
# 
# # END SECTION *****************************************************************
# 
# 
# 
# 
# 
# grid2  <- h2o.grid("gbm", 
#                  x = predictors, 
#                  y = response, 
#                 grid_id = 'grid2',
#                  training_frame = train,
#                  validation_frame = valid,
#                  fold_column = 'cv_group',
#                  hyper_params = list(
#                    distribution='tweedie',
#                    tweedie_power = c(1.05,1.1,1.5),
#                    ntrees = c(55),
#                    max_depth = c(7)
#                    ))
# summary(grid2)
# 
# gbm_gridperf2 <- h2o.getGrid(grid_id = grid2@grid_id,
#                              sort_by = "rmse",
#                              decreasing = F)
# print(gbm_gridperf2)
# gbm_gridperf2@summary_table %>% as.data.frame
# best_gbm2 <- h2o.getModel(gbm_gridperf2@model_ids[[1]])
# h2o.performance(best_gbm2,newdata = valid)@metrics$r2
# h2o.performance(best_gbm2,newdata = valid)@metrics$RMSE
# # rm(grid2)
# 
# dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
# dpreds <- bind_cols(h2o.predict(best_gbm2, as.h2o(dpreds)) %>% as.data.table(),
#                     dpreds)
# 
# 
# dpreds %>%
#   ggplot(data=.,aes(x,y,fill=predict))+
#   geom_raster()+
#   coord_equal()+
#   scale_fill_viridis_b(
#     option='B',
#     limits=c(365,2000),
#     oob=scales::squish,
#     n.breaks=6)+
#   # scico::scale_fill_scico(palette = 'batlowK',
#   #   direction = 1,
#   #   limits=c(365,2000),oob=scales::squish)+
#   theme_void()+
#   theme(panel.background = element_rect(fill='grey'))
# 
# gbm_gridperf2@model_ids
# h2o.varimp_heatmap(gbm_gridperf2@model_ids[1:3])
# 
# 
# vip::vip(best_gbm2, 30)
# dpreds %>%
#   # filter(pre_fire_slai_anom_3mo < 0.5) %>% 
#   ggplot(data=.,aes(pre_fire_slai_anom_3mo,predict,
#     color=cut_interval(malai,4)))+
#   # geom_point()
#   # geom_smooth()+
#   geom_smooth(method=MASS::rlm)+
#   geom_rug()+
#   facet_wrap(~cut_interval(pre_fire_slai_anom_12mo,4))+
#   theme_linedraw()
# 
# vip::vi_model(best_gbm2)
# pc <- h2o.predict_contributions(best_gbm2,train,top_n = 5)
# pc <- as.data.table(pc)
# pc$top_feature_1 %>% table %>% sort
# pc$top_feature_2 %>% table %>% sort
# pc$top_feature_3 %>% table %>% sort
# pc$top_feature_4 %>% table %>% sort
# pc$top_feature_5 %>% table %>% sort
# 
# best_gbm2
# x1 <- h2o.feature_interaction(best_gbm2, max_interaction_depth = 5, 
#   max_tree_depth = 5)
# x1[2] %>% 
#   as.data.table() %>% 
#   .[str_detect(interaction,"pre_fire_slai_anom_3mo")] %>% 
#   .[order(-gain)] %>% 
#   head()
# 
# x1[1:2] %>% 
#   map_dfr()
# 
# 
# 
# # # auto ml attempt -------------------------------------------------------------
# # a1 <- h2o.automl(x = predictors, 
# #                  y = response,
# #                  training_frame = train, 
# #                  validation_frame = valid,
# #                    include_algos = c('GBM','StackedEnsemble'),
# #                  # include_algos = c('XGBoost','GBM'),
# #                  fold_column = 'fire_year',
# #                  max_runtime_secs = 60*60,
# #                  seed = 1234)
# # lb <- h2o.get_leaderboard(a1,extra_columns = 'ALL')
# # print(lb,n = nrow(lb))
# # head(lb)
# # top_mod <- h2o.get_best_model(a1,criterion = 'RMSE')
# # h2o.saveModel(top_mod,
# #               path=paste0("../data_general/proc_data_Oz_fire_recovery/",
# #                           "xgboost_ttr5-lai_",Sys.time(),"_"))
# # 
# # h2o.varimp(top_mod)
# # vip::vip(top_mod, 25)
# # 
# # top_mod@parameters$ntrees     # 69
# # top_mod@parameters$max_depth  # 20
# # top_mod@parameters$col_sample_rate #0.8
# # top_mod@parameters$reg_lambda # 10
# # top_mod@parameters$distribution
# # top_mod@parameters
# # top_mod@allparameters
# # h2o.performance(top_mod)@metrics
# # h2o.performance(top_mod,newdata = valid)@metrics$r2
# # h2o.performance(top_mod,newdata = valid)@metrics$RMSE
# # 
# # 
# # top_mod <- h2o.loadModel(list.files(list.files("../data_general/proc_data_Oz_fire_recovery/", 
# #   pattern = "xgboost_ttr5-lai_",full.names = T) %>% 
# #   sort() %>% 
# #   last(),full.names = T))
# # 
# # dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
# # dpreds <- bind_cols(h2o.predict(top_mod, as.h2o(dpreds)) %>% as.data.table(),
# #                     dpreds)
# # 
# # 
# # dpreds %>% 
# #   ggplot(data=.,aes(x,y,fill=predict))+
# #   geom_raster()+
# #   coord_equal()+
# #   scale_fill_viridis_b(
# #     option='B',
# #     limits=c(365,2000),
# #     oob=scales::squish,
# #     n.breaks=6)+
# #   # scico::scale_fill_scico(palette = 'batlowK',
# #   #   direction = 1,
# #   #   limits=c(365,2000),oob=scales::squish)+
# #   theme_void()+
# #   theme(panel.background = element_rect(fill='grey'))
# # 
# # vip::vip(top_mod, 10)
# # dpreds %>% 
# #   ggplot(data=.,aes(sand,predict,
# #     color=cut_interval(malai,4)))+
# #   geom_smooth(method=MASS::rlm)
# # 
# 
# 
# # GBM hyperparameters
# gbm_params1 <- list(learn_rate = c(0.1),
#                     ntrees = 65,
#                     max_depth = c(7),
#   min_rows = c(10,20,40),
#                     sample_rate = c(1),
#                     col_sample_rate = c(1))
# 
# # Train and validate a cartesian grid of GBMs
# gbm_grid1 <- h2o.grid("gbm", x = predictors, y = response,
#                       grid_id = "gbm_grid1",
#                       training_frame = train,
#                       validation_frame = valid,
#                    fold_column = 'fire_year',
#                       seed = 1,
#                       hyper_params = gbm_params1)
# h2o.getGrid(grid_id = "gbm_grid1",
#                              sort_by = "rmse",
#                              decreasing = F)
# 
# best_gbm1 <- h2o.getModel(h2o.getGrid(grid_id = "gbm_grid1",
#                              sort_by = "rmse",
#                              decreasing = F)@model_ids[[1]])
# best_gbm1@allparameters
# h2o.performance(best_gbm1,newdata = valid)@metrics$r2
# h2o.performance(top_mod,newdata = valid)@metrics$r2
# h2o.performance(best_gbm1,newdata = valid)@metrics$RMSE
# h2o.performance(top_mod,newdata = valid)@metrics$RMSE
# vip::vip(best_gbm1)
# 
# vec_names <- h2o.getModel(h2o.getGrid(grid_id = "gbm_grid1",
#                              sort_by = "rmse",
#                              decreasing = F)@model_ids %>% unlist())
# 
# vec_names <- unlist(h2o.getGrid(grid_id = "gbm_grid1",
#                              sort_by = "rmse",
#                              decreasing = F)@model_ids)[1:5]
# ens <- h2o.stackedEnsemble(x = predictors, y = response, 
#   training_frame = train, 
#   validation_frame = valid, 
#   base_models = vec_names)
# 
# # Get the grid results, sorted by validation AUC
# gbm_gridperf1 <- h2o.getGrid(grid_id = "gbm_grid1",
#                              sort_by = "mean_per_class_error",
#                              decreasing = F)
# print(gbm_gridperf1)
# 
# # Grab the top GBM model, chosen by validation AUC
# best_gbm1 <- h2o.getModel(gbm_gridperf1@model_ids[[1]])
# 
# # Now let's evaluate the model performance on a test set
# # so we get an honest estimate of top model performance
# best_gbm_perf1 <- h2o.performance(model = best_gbm1,
#                                   newdata = valid)
# h2o.mean_per_class_error(best_gbm_perf1)
# 
# # Look at the hyperparameters for the best model
# print(best_gbm1@model[["model_summary"]])
# print(best_gbm1@model[["validation_metrics"]])
# print(best_gbm1@model[["variable_importances"]])
# best_gbm1@allparameters[["learn_rate"]]
# best_gbm1@allparameters[["max_depth"]]
# best_gbm1@allparameters[["sample_rate"]]
# best_gbm1@allparameters[["col_sample_rate"]]
# 
# 
# gbm_params2 <- list(learn_rate = c(0.09,0.1,0.2),
#                     max_depth = c(9,10,11,12,13),
#                     sample_rate = c(1.0),
#                     col_sample_rate = c(1.0))
# gbm_grid2 <- h2o.grid("gbm", x = predictors, y = response,
#                       grid_id = "gbm_grid2",
#                       training_frame = train,
#                       validation_frame = valid,
#                       ntrees = 100,
#                       seed = 1,
#                       hyper_params = gbm_params2)
# gbm_gridperf2 <- h2o.getGrid(grid_id = "gbm_grid2",
#                              sort_by = "mean_per_class_error",
#                              decreasing = F)
# print(gbm_gridperf2)
# best_gbm2 <- h2o.getModel(gbm_gridperf2@model_ids[[1]])
# best_gbm_perf2 <- h2o.performance(model = best_gbm2,
#                                   newdata = valid)
# h2o.mean_per_class_error(best_gbm_perf2)
# print(best_gbm2@model[["model_summary"]])
# print(best_gbm2@model[["validation_metrics"]])
# print(best_gbm2@model[["variable_importances"]])
# best_gbm2@allparameters[["learn_rate"]]
# best_gbm2@allparameters[["max_depth"]]
# best_gbm2@allparameters[["sample_rate"]]
# best_gbm2@allparameters[["col_sample_rate"]]
# 
# 
# 
# gbm_params3 <- list(learn_rate = c(0.2,0.3,0.4,0.5),
#                     max_depth = c(10,15,20),
#                     sample_rate = c(1.0),
#                     col_sample_rate = c(1.0))
# gbm_grid3 <- h2o.grid("gbm", x = predictors, y = response,
#                       grid_id = "gbm_grid3",
#                       training_frame = train,
#                       validation_frame = valid,
#                       ntrees = 100,
#                       seed = 1,
#                       hyper_params = gbm_params3)
# gbm_gridperf3 <- h2o.getGrid(grid_id = 'gbm_grid3',
#                              sort_by = "mean_per_class_error",
#                              decreasing = F)
# print(gbm_gridperf3)
# best_gbm3 <- h2o.getModel(gbm_gridperf3@model_ids[[1]])
# best_gbm_perf3 <- h2o.performance(model = best_gbm3,
#                                   newdata = valid)
# h2o.mean_per_class_error(best_gbm_perf3)
# print(best_gbm3@model[["model_summary"]])
# print(best_gbm3@model[["validation_metrics"]])
# print(best_gbm3@model[["variable_importances"]])
# best_gbm3@allparameters[["learn_rate"]]
# best_gbm3@allparameters[["max_depth"]]
# best_gbm3@allparameters[["sample_rate"]]
# best_gbm3@allparameters[["col_sample_rate"]]
# 
# 
# gbm_params4 <- list(learn_rate = c(0.2),
#                     ntrees=floor(seq(20,300,by=50)),
#                     max_depth = c(20),
#                     sample_rate = c(1.0),
#                     col_sample_rate = c(1.0))
# gbm_grid4 <- h2o.grid("gbm", x = predictors, y = response,
#                       grid_id = "gbm_grid4",
#                       training_frame = train,
#                       validation_frame = valid,
#                       seed = 1,
#                       hyper_params = gbm_params4)
# gbm_gridperf4 <- h2o.getGrid(grid_id = "gbm_grid4",
#                              sort_by = "mean_per_class_error",
#                              decreasing = F)
# print(gbm_gridperf4)
# best_gbm4 <- h2o.getModel(gbm_gridperf4@model_ids[[1]])
# best_gbm_perf4 <- h2o.performance(model = best_gbm4,
#                                   newdata = valid)
# h2o.mean_per_class_error(best_gbm_perf4)
# print(best_gbm4@model[["model_summary"]])
# print(best_gbm4@model[["validation_metrics"]])
# print(best_gbm4@model[["variable_importances"]])
# best_gbm4@allparameters[["learn_rate"]]
# best_gbm4@allparameters[["max_depth"]]
# best_gbm4@allparameters[["sample_rate"]]
# best_gbm4@allparameters[["col_sample_rate"]]
# 
# top_mod <- best_gbm4
# h2o.saveModel(top_mod,
#               path=paste0("../data_general/proc_data_Oz_fire_recovery/",
#                           "gbm_ttr5-lai-ocat_",Sys.Date(),"_"))
# 
# dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
# dpreds <- bind_cols(h2o.predict(top_mod, as.h2o(dpreds)) %>% as.data.table(),
#          dpreds)
# arrow::write_parquet(dpreds,
#                    sink=paste0("../data_general/proc_data_Oz_fire_recovery/predicted_ttr5-lai-5class_gbm_",Sys.time(),"_.parquet"),
#                    compression='snappy')
# 
# 
# 
# 
# 
# # SCRATCH ---------------------------------------------------------------------
# # dpreds$predict %>% as.numeric() %>% hist
# # dpreds %>%
# #   ggplot(data=.,aes(x,y,fill=predict))+
# #   # geom_sf(data=oz_poly,
# #   #         fill='grey70',
# #   #         color='grey30',
# #   #         inherit.aes = F)+
# #   geom_raster()+
# #   coord_equal()+
# #   scale_fill_viridis_d(option='B',end=0.95)+
# #   coord_sf(xlim = c(146,153.25),
# #            ylim = c(-38,-27.75))
# # 
# # dpreds %>% 
# #   mutate(y_group = cut(y,c(-40,-35,-31,-27), 
# #                        labels=c("south","mid","north"))) %>% 
# #   ggplot(data=.,aes(y=predict,
# #                     fill=y_group,
# #                     x=elevation))+
# #   geom_violin(outlier.colour = NA)+
# #   facet_wrap(cut_number(min_nbr_anom,3)~y_group,scales='free')
# # 
# # 
# # # Set predictors and response; set response as a factor
# # colnames(d_rf)
# # response <- "ttr_ocat"
# # predictors <- setdiff(colnames(d_rf), response)
# # predictors <- predictors[str_detect(predictors,'fire_year')==F]
# # 
# # hdat2 <- as.h2o(d_rf %>% filter(fire_year!=2013))
# # splits2 <- h2o.splitFrame(data =  hdat2, ratios = .75, seed = 1234)
# # train2 <- splits2[[1]]
# # valid2 <- splits2[[2]]
# # 
# # rtop <- h2o.gbm(x = predictors, 
# #                  y = response,
# #                  training_frame = train2, 
# #                  validation_frame = valid2,
# #                  fold_column = 'fire_year',
# #                 ntrees = 100,#81
# #                 max_depth=15,#7
# #                 min_rows=5,
# #                 nbins = 20,
# #                 nbins_top_level = 1024,
# #                 nbins_cats = 1024,
# #                 stopping_metric = 'deviance',
# #                 seed=10999,
# #                 learn_rate=0.1,
# #                 learn_rate_annealing = 1,
# #                 distribution='multinomial',
# #                 huber_alpha = 0.9,
# #                 sample_rate=1,
# #                 col_sample_rate = 0.4,
# #                 col_sample_rate_per_tree = 0.4)
# # 
# # h2o.performance(rtop,newdata = valid2)
# # 
# # 
# # vec_tru <- d_rf %>% filter(fire_year==2013) %>% pull(ttr_ocat)
# # vec_pred <- h2o.predict(rtop,newdata=as.h2o(d_rf %>% 
# #                                   filter(fire_year==2013)))
# # vec_pred <- vec_pred %>% as.data.table
# # 
# # 
# # data.table(pred=vec_pred$predict, 
# #            obs=vec_tru) %>% 
# #   ggplot(data=.,aes(as.numeric(pred),as.numeric(obs)))+
# #   geom_point(position='jitter',
# #              size=0.1)+
# #   scale_color_viridis_c()+
# #   # geom_abline(col='red')+
# #   geom_smooth(method='lm')
# # 
# # 
# # 
# # 
# # 
# # 
