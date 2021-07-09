pacman::p_load(tidyverse, data.table, lubridate, patchwork)

dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                           col_select = c("x","y","id","date","slai","slai_anom_12mo","malai"))
dat[,`:=`(slai_12mo = slai_anom_12mo+malai)]
sdat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS_2021-06-05 13:01:38.parquet")
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-1mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-06-19 17:33:16.parquet")
fits <- fits[isConv==TRUE][r2>0][,ldk:=L0/K]
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
fits[,fire_year := lubridate::year(date_fire1 - months(3))]

# nvis ---------------------
nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif", 
                          proxy=F)
names(nvis) <- 'vc_code'
nvis_codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("vc_code","veg_class_descrip")) %>% 
  mutate(vc_name = as.factor(veg_class_descrip)) %>% 
  select(-veg_class_descrip)



# ML SDM predictions ---------------------
out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-07 14:54:12.parquet")
out %>% names

s_sdm <- st_as_stars(out[,.(x,y,predict)], dims=c('x','y'))
st_crs(s_sdm) <- st_crs(4326)

s_fits <- st_as_sf(fits[,.(x,y,pred_ttr,fire_year,r,L0,K)], coords=c('x','y'), 
  crs=4326)

vv <- st_extract(s_sdm, s_fits)
s_fits <- bind_cols(s_fits, st_drop_geometry(vv))

vv <- st_extract(nvis, s_fits)
s_fits <- bind_cols(s_fits, st_drop_geometry(vv))

s_fits <- st_drop_geometry(s_fits) %>% as.data.table()
s_fits <- merge(s_fits, nvis_codes, by='vc_code')
s_fits <- s_fits[vc_code %in% c(2,3,5,11)]

s_fits$vc_name %>% unique

dp <- s_fits %>% 
  as_tibble() %>% 
  mutate(ttr_ocat = case_when(pred_ttr <= 366 ~ 1,
                              between(pred_ttr, 367, 365*2)~2,
                              between(pred_ttr, 731, 1096)~3,
                              between(pred_ttr, 1097, 1460)~4,
                              between(pred_ttr, 1461, Inf)~5,
                              is.na(pred_ttr)==T ~ 5
  )) %>% 
  filter(is.na(predict)==F) %>% 
  filter(is.na(ttr_ocat)==F) %>% 
  mutate(ttr_ocat = factor(ttr_ocat)) %>% 
  mutate(predict = fct_drop(predict), 
         ttr_ocat = fct_drop(ttr_ocat))


obs_count <- dp %>% 
  group_by(predict, fire_year, vc_name) %>% 
  summarize(nobs = n()) %>% 
  group_by(predict, vc_name) %>% 
  summarize(tot_nobs = sum(nobs), 
            n_fy = length(unique(nobs))) %>% 
  ungroup()


merge(dp, obs_count, by=c("predict","vc_name")) %>% 
  filter(L0/K <1) %>% 
  ggplot(data=.,aes(x=(L0/K), 
    y=vc_name))+
  ggridges::stat_density_ridges(quantile_lines=T, 
    bandwidth = 0.01)
  
scale_x_continuous(expand=c(0,0), 
                     limits=c(0,1))+
  # scale_fill_viridis_d(option='B', end=0.95, direction = 1, 
  #                      labels=c("≤ 1","2","3","4","≥ 5"))+
  labs(y=NULL,
       x='L0/K',
       fill='Years to Recover')+
  # guides(fill = guide_colorsteps())
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(0,0.5,0,0),'cm'))

merge(dp, obs_count, by=c("predict","vc_name")) %>% 
  # filter(n_fy >= 5 & tot_nobs >= 50) %>% 
  # filter(is.na(predict)==F) %>% 
  # mutate(predict = fct_drop(predict)) %>% 
  ggplot(data=.,aes(fill=factor(ttr_ocat), 
    y=vc_name))+
  geom_bar(position=position_fill(reverse=T))+
    scale_y_discrete(limits = rev)+
  scale_x_continuous(expand=c(0,0), 
                     limits=c(0,1))+
  scale_fill_viridis_d(option='B', end=0.95, direction = 1, 
                       labels=c("≤ 1","2","3","4","≥ 5"))+
  labs(y=NULL,
       x='proportion',
       fill='Years to Recover')+
  # guides(fill = guide_colorsteps())
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(0,0.5,0,0),'cm'))



merge(dp, obs_count, by=c("predict","vc_name")) %>% 
  filter(n_fy >= 5 & tot_nobs >= 50) %>% 
  # filter(is.na(predict)==F) %>% 
  # mutate(predict = fct_drop(predict)) %>% 
  ggplot(data=.,aes(fill=factor(ttr_ocat), 
    y=predict))+
  geom_bar(position=position_fill(reverse=T))+
    scale_y_discrete(limits = rev)+
  scale_x_continuous(expand=c(0,0), 
                     limits=c(0,1))+
  scale_fill_viridis_d(option='B', end=0.95, direction = 1, 
                       labels=c("≤ 1","2","3","4","≥ 5"))+
  labs(y=NULL,
       x='proportion',
       fill='Years to Recover')+
  # guides(fill = guide_colorsteps())
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(0,0.5,0,0),'cm'))+
  facet_wrap(~vc_name, scales='free', ncol=5)



  
  
d_rf %>% 
  mutate(ttr_cat = fct_rev(ttr_ocat)) %>% 
  ggplot(data=.,aes(y=vc_name_f, 
                    fill=ttr_ocat))+
  geom_bar(position=position_fill(reverse=T))+
  scale_y_discrete(limits = rev)+
  scale_x_continuous(expand=c(0,0), 
                     limits=c(0,1))+
  scale_fill_viridis_d(option='B', end=0.95, direction = 1, 
                       labels=c("≤ 1","2","3","4","≥ 5"))+
  labs(y=NULL,
       x='proportion',
       fill='Years to Recover')+
  # guides(fill = guide_colorsteps())
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(0,0.5,0,0),'cm'))

