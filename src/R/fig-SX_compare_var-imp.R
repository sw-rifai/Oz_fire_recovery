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

# Plot Variable Importance Comparison -------------------------------------
s0 <- h2o.loadModel("outputs/GBM_bestModNoSpecies_2022-02-08 13:43:04/")
s2 <- h2o.loadModel("outputs/GBM_bestModWithSpecies_2022-02-08 13:43:04/")


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
