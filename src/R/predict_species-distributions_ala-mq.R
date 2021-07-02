library(tidyverse)
library(stars); 
library(dtplyr)
library(data.table) 
library(lubridate)
library(arrow)
library(furrr)

# Oz poly -------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 1000) %>% 
  select(NAME_1)
oz_poly <- oz_poly %>% filter(NAME_1 %in% c("Queensland","New South Wales","Australian Capital Territory","Victoria"))

# nvis ---------------------
nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif", 
                          proxy=F)
names(nvis) <- 'vc_code'
nvis_codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("vc_code","veg_class_descrip")) %>% 
  mutate(vc_name = as.factor(veg_class_descrip)) %>% 
  select(-veg_class_descrip)


# ALA MQ cleaned --------------
ala_mq <- fread("../data_general/ALA/Gallagher_cleaned/euc_occurrence_clean.csv") %>% 
  rename( 
    x=longitude,
    y=latitude) %>% 
  as.data.table() %>% 
  .[x>=140 & y <= -23 & y>=-40]

vv <- st_extract(nvis, 
           st_as_sf(ala_mq[,.(x,y)], coords=c("x","y"), crs=4326))
ala_mq <- bind_cols(ala_mq, st_drop_geometry(vv))
ala_mq <- merge(ala_mq[is.na(vc_code)==F], nvis_codes, by='vc_code')
ala_mq[vc_code %in% c(2,3,4,5,11)]
unique(ala_mq[,.(vc_name,vc_code)])

## Filter ALA by NVIS Euc classes 
# ala_mq[,.(vc_name,vc_code)][,.(nobs=.N),by='vc_code'][order(-nobs)][nobs>=400]$vc_code -> junk
# ala_mq %>% .[vc_code %in% c(2,3,4,5,11)] %>%
# # ala_mq[vc_code%in%junk] %>% 
#   ggplot(data=.,aes(x,y,color=vc_name))+
#   geom_point(size=0.25)+coord_sf()+scale_color_viridis_d(option='H')
ala_mq <- ala_mq %>% .[vc_code %in% c(2,3,4,5,11)]

## Simplify species names
fn <- function(x){
  x <- str_remove(x, " sp\\.")
  x <- str_remove(x, " subsp\\.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}
ala_mq <- ala_mq[,species := fn(taxon),by=seq_len(nrow(ala_mq))]


# Topography vars ---------------
rtpi <- read_stars("../data_general/Oz_misc_data/mtpi_alos_500m_SE_coastal.tif") %>% 
  set_names('tpi')
rdem <- read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif") %>% 
  set_names('elevation')
rslope <- read_stars("../data_general/Oz_misc_data/slope_wwfhydrosheds_500m_SE_coastal.tif") %>%   set_names('slope')
raspect <- read_stars("../data_general/Oz_misc_data/aspect_wwfhydrosheds_500m_SE_coastal.tif") %>% 
  set_names('aspect')
rhnd <- read_stars("../data_general/Oz_misc_data/mert_hnd_500m_SE_coastal.tif") %>% 
  set_names("hnd")

# Subset NVIS -------------------
nvis <- st_crop(nvis, rtpi)
nvis %in% c(2,3,4,5)
alt <- nvis
alt <- alt %>% mutate(vc_code = ifelse(vc_code%in%c(2,3,4,5), vc_code,NA))

# Load Climate -----------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y","month",
                                           "map","mapet","mappet",
                                           "matmax","matmin","mavpd15", 
                                           "precip_u",
                                           "tmax_u",
                                           "tmin_u",
                                           "vpd15_u",
                                           "vpd9_u")) %>% 
  distinct() %>% 
  as.data.table()

clim_ma <- clim[,.(x,y,map,mapet,mappet,matmax,matmin,mavpd15)] %>% unique()


lut_seasons <- data.table(month=1:12)[,season := case_when(month%in%c(12,1,2)~"DJF",
                                            month%in%c(3,4,5)~"MAM",
                                            month%in%c(6,7,8)~"JJA",
                                            month%in%c(9,10,11)~"SON")]

clim <- merge(clim, lut_seasons,by='month')
clim <- clim[,.(precip_u = sum(precip_u), 
        tmax_u = mean(tmax_u), 
        tmin_u = mean(tmin_u), 
        vpd15_u = mean(vpd15_u), 
        vpd9_u = mean(vpd9_u)),.(x,y,season)]


# clim_precip <- dcast(clim[,.(x,y,season,precip_u)], x+y~paste0("precip_month",clim[,seq_len(.N),by=.(x,y)]$V1), 
#       value.var = 'precip_u')

clim_precip <- dcast(clim[,.(x,y,season,precip_u)], x+y~paste0("precip_season_",season), 
      value.var = 'precip_u')
clim_tmax <- dcast(clim[,.(x,y,season,tmax_u)], x+y~paste0("tmax_season_",season),                      value.var = 'tmax_u')
clim_tmin <- dcast(clim[,.(x,y,season,tmin_u)], x+y~paste0("tmin_season_",season),                    value.var = 'tmin_u')
clim_vpd15 <- dcast(clim[,.(x,y,season,vpd15_u)], x+y~paste0("vpd15_u_season_",season),value.var = 'vpd15_u')
clim_vpd9 <- dcast(clim[,.(x,y,season,vpd9_u)], x+y~paste0("vpd9_u_season_",season),value.var = 'vpd9_u')

## Cast clim wide --------------------------------
clim_wide <- merge(clim_ma, clim_precip, by=c('x','y'))
clim_wide <- merge(clim_wide, clim_tmax, by=c('x','y'))
clim_wide <- merge(clim_wide, clim_tmin, by=c('x','y'))
clim_wide <- merge(clim_wide, clim_vpd15, by=c('x','y'))
clim_wide <- merge(clim_wide, clim_vpd9, by=c('x','y'))



# LAI vars -------------------------------------
smalai <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_mean-annual-lai.tif", 
                            proxy=F)
m1 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-1.tif", 
                        proxy=F)
m2 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-2.tif", 
                        proxy=F)
m3 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-3.tif", 
                        proxy=F)
m4 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-4.tif", 
                        proxy=F)
m5 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-5.tif", 
                        proxy=F)
m6 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-6.tif", 
                        proxy=F)
m7 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-7.tif", 
                        proxy=F)
m8 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-8.tif", 
                        proxy=F)
m9 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-9.tif", 
                        proxy=F)
m10 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-10.tif", 
                         proxy=F)
m11 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-11.tif", 
                         proxy=F)
m12 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_median-lai-month-12.tif", 
                         proxy=F)

m_djf <- (m12+m1+m2)/3
m_mam <- (m3+m4+m5)/3
m_jja <- (m6+m7+m8)/3
m_son <- (m9+m10+m11)/3

slai <- c(smalai, m_djf, m_mam, m_jja, m_son)
names(slai) <- c("malai",
                 "lai-djf","lai-mam","lai-jja", 
                 "lai-son")
rm(smalai, m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12, 
   m_djf, m_mam, m_jja, m_son)

# slai <- c(smalai, m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)
# names(slai) <- c("malai",
#                  "lai-m1","lai-m2","lai-m3", 
#                  "lai-m4","lai-m5","lai-m6",
#                  "lai-m7","lai-m8","lai-m11",
#                  "lai-m9","lai-m10","lai-m12")
# rm(smalai, m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)
gc(full=T)


# Soil variables -----------------------------------
s_soil <- read_stars("../data_general/Oz_misc_data/SLGA_500m_SE_coastal.tif", 
                     proxy = F) 
st_crs(s_soil) <- st_crs(4326)
s_soil <- split(s_soil,f = 'band')


# prepares stars extractions ----------------------------------------
dat <- st_as_sf(ala_mq, coords = c("x","y"))
st_crs(dat) <- st_crs(4326)

tmp <- st_as_stars(clim_wide,dims=c('x','y'))
st_crs(tmp) <- st_crs(4326)
vv <- st_extract(tmp, 
                 at=dat)
vv
dat <- bind_cols(dat,st_drop_geometry(vv))

s_soil

vv <- st_extract(s_soil, 
                 dat)

dat <- bind_cols(dat, st_drop_geometry(vv))

vv <- st_drop_geometry(st_extract(c(rtpi,rdem,rslope,raspect), at=dat))
dat <- bind_cols(dat,vv)

vv <- st_drop_geometry(st_extract(c(rhnd), at=dat))
dat <- bind_cols(dat,vv)

dat <- bind_cols(dat, 
          as.data.table(st_coordinates(dat)) %>% set_names(c("x","y")))



fdat <- as.data.table(dat)
fdat <- fdat[is.na(hnd)==F]  
fdat <- fdat[is.na(slope)==F]  
fdat <- fdat[is.na(aspect)==F]  
fdat <- fdat[is.na(map)==F]
fdat <- fdat[is.na(AWC_AWC_000_005_EV)==F]  
fdat <- fdat[is.na(BDW_BDW_100_200_EV)==F]  
fdat <- fdat[is.na(species)==F][species!="Eucalyptus x"]
apply(fdat, 2, function(x) sum(is.na(x)))

# Dropping LAI given the difference between ALA sample age and LAI record
# vv <- st_extract(slai, fdat %>% select(x,y) %>% 
#                    as.data.table() %>% 
#                    st_as_sf(., coords=c('x','y'),crs=4326))
# fdat <- bind_cols(fdat,st_drop_geometry(vv))
fdat <- fdat[year>=1950]

# fdat[species == "Eucalyptus regnans"] %>% dim
# fdat %>% select(x,y) %>% as_tibble() %>% 
#   ggplot(data=.,aes(x,y))+
#   geom_point()+coord_sf()

# Filter to species with 80+ nob ------------------------------------------
vec_nobs_80 <- fdat[,.(nobs=.N),by='species'][,rank:=frank(-nobs,ties.method = 'first')][order(rank)][nobs>=80]$species
fdat2 <- fdat[species%in%vec_nobs_80]
fdat2[,species:=as.factor(species)]
fdat2 <- fdat2[,-c("geometry","year","taxon","vc_name","vc_code")]
dim(fdat2)

# H2O xgboost model fitting ----------------------------------------------------
library(h2o)
h2o.init()
Sys.setenv(sys.ai.h2o.auc.maxClasses=81)
Sys.getenv("sys.ai.h2o.auc.maxClasses")

# Set predictors and response; set response as a factor
hdat <- as.h2o(fdat2)
hdat['species'] %>% as.data.table() %>% unique

colnames(fdat2)
response <- "species"
predictors <- setdiff(colnames(fdat2), response)


# Split the dataset into train and valid

splits <- h2o.splitFrame(data =  hdat, ratios = .75, seed = 1234)
train <- splits[[1]]
valid <- splits[[2]]


# auto ml attempt
a1 <- h2o.automl(x = predictors, 
                 y = response,
                 training_frame = train, 
                 validation_frame = valid,
                 balance_classes = TRUE,
                 include_algos = 'XGBoost',
                 max_models = 30,
                 sort_metric = "mean_per_class_error",
                 # max_runtime_secs = 3600,
                 seed = 1234)
lb <- h2o.get_leaderboard(a1,extra_columns = 'ALL')
print(lb,n = nrow(lb))
head(lb)
top_mod <- h2o.get_best_model(a1,criterion = 'mean_per_class_error')
h2o.varimp(top_mod)
vip::vip(top_mod, 25)
h2o.saveModel(top_mod,
              path=paste0("../data_general/proc_data_Oz_fire_recovery/",
                          "xgboost_euc-sdm-gte80nobs_",Sys.Date(),"_"))

top_mod@parameters$ntrees     # 69
top_mod@parameters$max_depth  # 20
top_mod@parameters$col_sample_rate #0.8
top_mod@parameters$reg_lambda # 10
top_mod@parameters$distribution
top_mod@parameters
top_mod@allparameters
h2o.performance(top_mod)@metrics$mean_per_class_error #0.153 # 0.2(80 forest)
h2o.performance(top_mod)@metrics$logloss # 0.614
h2o.performance(top_mod)@metrics$MSE
h2o.performance(top_mod)@metrics$RMSE
# MSE 0.21
# RMSE 0.459



# 
# hyperparameters_xgboost <- list(ntrees = seq(20, 200, 20),
#                                 max_depth = seq(4,20,1),
#                                 learn_rate = seq(0.05, 0.3, 0.05),
#                                 sample_rate = c(0.5,1,0.25),
#                                 col_sample_rate = seq(0.5, 1, 0.1))
# hyperparameters_xgboost
# search_criteria <- list(strategy = "RandomDiscrete", max_models = 30, seed = 1)
# xg_grid1 <- h2o.grid('xgboost', 
#                      y=response, 
#                      x=predictors, 
#                      training_frame = train, 
#                      validation_frame = valid,
#                      seed=3,
#                      hyper_params = hyperparameters_xgboost,
#                      search_criteria = search_criteria)
# print(xg_grid1)
# xg_grid1@grid_id
# xg_grid1 %>% str
# xg_grid1@summary_table$logloss
# g_perf <- h2o.getGrid(xg_grid1@grid_id, 
#             sort_by = "logloss", 
#             decreasing = F)
# print(g_perf)
# 
# best_xg <- h2o.getModel(g_perf@model_ids[[1]])
# best_xg@model$cross_validation_predictions
# best_xg@model$variable_importances
# best_xg@model$domains
# best_xg@model$native_parameters
# best_xg@parameters
# 
# 
# vip::vip(best_xg, 30)
# 
# best_xg <- h2o.getModel(xg_grid1@model_ids[[1]])
# best_xg@model$cross_validation_predictions
# xg_grid1@model_ids
# h2o.performance(best_xg)@metrics$logloss
# h2o.performance(best_xg)@metrics$mean_per_class_error
# h2o.performance(best_xg) %>% str
# 
# xg_old@allparameters
# # Train the XGB model
# # xg1 <- h2o.xgboost(x = predictors, 
# #                         y = response,
# #                         training_frame = train, 
# #                         validation_frame = valid,
# #                         distribution='multinomial',
# #                         booster = "dart", 
# #                    normalize_type = "tree",
# #                    seed = 1234)
# # xg_old <- xg1
# 
xg1 <- h2o.xgboost(x = predictors,
                   y = response,
                   training_frame = train,
                   validation_frame = valid,
                   seed=1340,
                   stopping_tolerance = 0.008571778,
                   fold_assignment = 'Modulo',
                   categorical_encoding = "OneHotInternal",
                   nfolds = 5,
                   ntrees = 500,
                   max_depth = 69,
                   min_rows = 10,
                   min_child_weight = 10,
                   col_sample_rate_per_tree = 0.8, 
                   colsample_bytree = 0.8,
                   score_tree_interval = 5,
                   reg_lambda = 10,
                   reg_alpha = 0.01,
                   dmatrix_type = 'dense',
                   learn_rate = 0.05,
                   eta = 0.6,
                   sample_rate=0.6,
                   distribution='multinomial',
                   tree_method='exact',
                   backend='cpu',
                   sample_type = 'uniform',
                   normalize_type = 'tree',
                   max_bins = 256,
                   calibrate_model = F,
                   tweedie_power = 1.5,
                   keep_cross_validation_predictions = T,
                   grow_policy = 'depthwise',
                   # booster = "dart",
                   # normalize_type = "tree",
                   # auc_type="MACRO_OVR",
                   stopping_metric = "mean_per_class_error")
# 
# # Error: n-incorrect/n-total
# # Rate: n-incorrect/n-total
# # hit ratio: correct/n-total
p1 <- h2o.performance(xg1,valid)
methods(class=class(p1))
p1_cm <- p1@metrics$cm$table %>% as.data.table()
p1@metrics$mean_per_class_error # 0.5
p1@metrics$MSE # 0.5
p1@metrics$logloss # 2.28
p1@metrics$r2 # 0.999
p1@metrics$RMSE # 0.709
p1@metrics
p1@metrics$multinomial_aucpr_table
h2o.gainsLift(p1)
ff <- h2o.feature_frequencies(xg1,newdata = valid)
vip::vip(xg1,30)
# 
# 
# h2o.multinomial_auc_table(xg1)
# z <- p1_cm[,-c("Error","Rate")] %>% as.matrix()
# z
# image(z[,ncol(z):1], axes=FALSE)
# 
# 
# p1@metrics$frame
# p1@metrics$mean_per_class_error
# p1@metrics$cm
# p1@metrics$hit_ratio_table
# 
# 
# 
# fdat2[,.(nobs=.N),by=species][,rank:=frank(-nobs)][order(rank)][rank<=5]
# pr <- h2o.predict(xg1,newdata = valid)
# pr <- pr %>% as.data.table()
# pr0 <- valid %>% as.data.table()
# pr1 <- bind_cols(pr0,pr)
# pr1[species=='Eucalyptus melliodora'] %>% 
#   ggplot(data=.,aes(x,y))+
#   geom_sf(data=oz_poly, inherit.aes = F)+
#   geom_point(size=2)+
#   geom_point(data=pr1[predict=='Eucalyptus melliodora'], 
#              aes(x,y),col='red')+
#   coord_sf(ylim=c(-40,-27), 
#            xlim=c(145,154))
# 
# 
gc()
pc <- expand_grid(x=seq(145,154,by=0.01), 
            y=seq(-40,-27,by=0.01)) %>% 
  st_as_sf(., coords=c("x","y"), crs=4326)
tmp <- st_as_stars(clim_wide,dims=c('x','y'))
st_crs(tmp) <- st_crs(4326)
vv <- st_extract(tmp, 
                 at=pc)
pc <- bind_cols(pc,st_drop_geometry(vv))
vv <- st_extract(s_soil, 
                 pc)
pc <- bind_cols(pc, st_drop_geometry(vv))
vv <- st_drop_geometry(st_extract(c(rtpi,rdem,rslope,raspect), at=pc))
pc <- bind_cols(pc,vv)
vv <- st_drop_geometry(st_extract(c(rhnd), at=pc))
pc <- bind_cols(pc,vv)
vv <- st_extract(slai, pc)
pc <- bind_cols(pc, st_drop_geometry(vv))
pc <- bind_cols(st_drop_geometry(pc), 
          st_coordinates(vv) %>% as.data.table() %>% set_names(c("x","y")))
pc <- na.omit(pc)


pc1 <- h2o.predict(top_mod,newdata = as.h2o(pc))
pc1 <- pc1 %>% as.data.table()
out <- bind_cols(pc,pc1) %>% as.data.table()
# sout <- out[,.(x,y,predict)] %>% 
#   # .[,predict := as.character(predict)] %>% 
#   st_as_stars(., dims=c('x','y'))
# st_crs(sout) <- st_crs(4326)
# # plot(sout)
# mask <- st_warp(alt, sout)
# mask <- mask %>% mutate(val = ifelse(is.na(vc_code)==T,F,T)) %>% select(val)
# mask$val %>% dim
# sout$predict %>% dim
# sout$predict[mask$val==F] <- NA
# plot(sout)
# sout <- c(mask,sout) %>% mutate(predict = ifelse(val==T,predict,NA)) %>%
#   select(predict)


arrow::write_parquet(out, 
                     sink=paste0("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_",Sys.time(),".parquet"))

write_csv(
  data.table(species=forcats::fct_unique(fdat2$species)) %>% 
    .[,species_fac_n:=as.numeric(species)], 
  file = "../data_general/proc_data_Oz_fire_recovery/predicted_species-distribution-ala-mq_nobs80-species_LUT.csv")
stars::write_stars(sout, 
                   dsn = paste0("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_",Sys.time(),".tiff"))



pc1 <- h2o.predict(top_mod,newdata = as.h2o(pc))
pc1 <- pc1 %>% as.data.table()
out <- bind_cols(as.data.table(pc)[,.(x,y)],pc1) %>% as.data.table()
sout <- out %>% 
  st_as_stars(., dims=c('x','y'))
st_crs(sout) <- st_crs(4326)
mask <- st_warp(alt, sout)
mask <- mask %>% mutate(val = ifelse(is.na(vc_code)==T,F,T)) %>% select(val)
sout[names(sout)][mask$val==F]
sout[mask$val==F] <- NA
plot(sout)
sout <- c(mask,sout) %>% mutate(predict = ifelse(val==T,predict,NA)) %>%
  select(predict)
stars::write_stars(sout, 
                   layer=names(sout),
                   dsn = paste0("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_",Sys.time(),".tiff"))

junk <- read_stars("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-06-25 09:34:41.tiff")


pc1[,.(nobs = .N), by=predict][,rank:=frank(-nobs)][order(rank)][rank<=20]
vec_top40 <- pc1[,.(nobs = .N), by=predict][,rank:=frank(-nobs)][order(rank)][rank<=100]$predict
out %>% #[predict %in% vec_top20] %>% 
  ggplot(data=.,aes(x,y,fill=predict))+
  geom_sf(data=oz_poly, 
          inherit.aes = F)+
  geom_raster()+
  scale_fill_viridis_d(option='H')+
  coord_sf(ylim=c(-38.75,-28), 
           xlim=c(145.5,153.5))+
  labs(x=NULL, 
       y=NULL)+
  guides(fill = guide_legend(ncol=2))+
  theme(legend.position = 'right', 
        panel.background = element_rect(fill='lightblue'))
ggsave(filename=paste0('figures/map_top40species_xgboost-prediction_',Sys.Date(),".png"), 
       width=35, 
       height=30, 
       units='cm', 
       dpi=350)

out[,.(x,y,predict)] %>% st_as_stars(., dims=c('x','y')) %>% plot
plot(sout)

unique(out$predict)
levels(out$predict)
tibble(species = unique(out$predict)) %>% 
  mutate(sp = as.numeric(species)) %>% 
  filter(sp==77)

out[,.(nobs=.N),by=predict][order(-nobs)]$nobs %>% plot
ala_mq[,.(nobs=.N),by=species][order(-nobs)]$nobs %>% plot

out %>% #[predict %in% vec_top20] %>% 
  filter(predict=="Eucalyptus blakelyi") %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=predict))+
  geom_sf(data=oz_poly, 
          inherit.aes = F)+
  geom_raster()+
  geom_point(data=ala_mq[species=="Eucalyptus blakelyi"], 
             aes(x,y),color='red',inherit.aes = F)+
  scale_fill_viridis_d(option='H')+
  coord_sf(ylim=c(-38.75,-28), 
           xlim=c(145.5,153.5))+
  labs(x=NULL, 
       y=NULL)+
  guides(fill = guide_legend(ncol=2))+
  theme(legend.position = 'right', 
        panel.background = element_rect(fill='lightblue'))

# h2o.shutdown()

