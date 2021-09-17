pacman::p_load(tidyverse, data.table, lubridate)

delta_t <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_delta_t_ttrDef5_preBS2021-04-25 16:28:57.parquet")

# check if this is the right one
ndvi <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet")

lai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS_2021-06-05 13:01:38.parquet")

cci <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_cci_ttrDef5_preBS2021-04-19 13:06:11.parquet")

kndvi <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_kndvi_ttrDef5_preBS2021-04-21 16:20:09.parquet")



delta_t <- delta_t[,.(id,ttr5_delta_t)] %>% rename(delta_t = ttr5_delta_t)
ndvi <- ndvi[,.(id,ttr5)] %>% rename(ndvi=ttr5)
kndvi <- kndvi[,.(id,ttr5_kn)] %>% rename(kndvi=ttr5_kn)
lai <- lai[,.(id,ttr5_lai)] %>% rename(lai=ttr5_lai)
cci <- cci[,.(id,ttr5_cci)] %>% rename(cci=ttr5_cci)

merge(delta_t, ndvi, by='id') %>% 
  merge(., kndvi, by='id') %>% 
  merge(., lai, by='id') %>% 
  merge(., cci,  by='id') %>% 
  drop_na() %>% 
  select(-id) %>% 
  rename(Tc_Ta = delta_t) %>% 
  select(cci, kndvi,lai,ndvi,Tc_Ta) %>% 
  cor() %>% 
  kableExtra::kable(.,format='pipe',digits=2)
