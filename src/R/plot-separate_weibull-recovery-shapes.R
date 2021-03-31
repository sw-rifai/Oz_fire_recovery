library(phenofit);
library(tidyverse);
library(usethis);
library(stars);
library(data.table); 
# library(dtplyr); 
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
source("src/R/functions_time_to_recover.R")

# Prep data -------------------------------------------------------------------
load("outputs/pixel_vegClass_groups.rds")
sdat <- read_parquet("outputs/weibull_fits_pre2005_fires_2021-01-18.parquet")
tmp2 <- expand_grid(merge(sdat,nvis, by='id') %>% 
                      filter(vc!=25) %>%
                      filter(vc %in% c(2,3,4,5,11)) %>% 
                      filter(is.na(vc)==FALSE) %>% 
                      sample_n(100), 
                    pred_days=seq(1,2000,length.out=2000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table()


pfit1 <- sdat %>% 
  summarize(across(contains(c("Asym","Drop",'lrc','pwr')), ~scale(.x))) %>% 
  prcomp(~Asym+Drop+lrc+pwr, data=.) 
pfit1 %>% summary
sdat2 <- bind_cols(sdat,data.frame(predict(pfit1)) %>% as_tibble())

tmp3 <- inner_join(tmp2, sdat2 %>% select(id,PC1,PC2,PC3,PC4))
  
tmp3 %>% ggplot(data=.,aes(pred_days, pred, color=PC1, group=id))+
  geom_line()+
  scale_color_viridis_c(option='B')+
  facet_wrap(~vc)

tmp3 %>% ggplot(data=.,aes(pred_days, p_diff, color=PC2, group=id))+
  geom_line()+
  scale_color_viridis_c(option='B')+
  facet_wrap(~cut_interval(PC1,4), scales = 'free')


tmp3 %>% select(Asym,Drop,pwr,lrc,PC1,PC2,PC3,PC4) %>% 
  sample_n(1000) %>% 
  plot


tmp3 %>% ggplot(data=.,aes(pred_days, pred, color=pwr, group=id))+
  geom_line()+
  scale_color_viridis_c(option='B')+
  facet_wrap(~cut_number(lrc,4), scales = 'fixed')

tmp3 %>% 
  filter(lrc > -10) %>% 
  ggplot(data=.,aes(pred_days, pred, color=pwr, group=id))+
  geom_line()

tmp3 %>% ggplot(data=.,aes(pred_days, pred, color=lrc, group=id))+
  geom_line()+
  scale_color_viridis_c(option='B')+
  facet_wrap(~cut_number(pwr,4), scales = 'fixed')


# Asym-Drop*exp(-exp(lrc)*x^pwr) 
expand_grid(data.frame(t(colMeans(sdat))) %>% select(-pwr),
            pwr=seq(-0.1,0.1,length.out = 5),
            pred_days=seq(1,2000,length.out=2000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table() %>% 
  ggplot(data=.,aes(pred_days,pred,color=factor(pwr)))+
  geom_line()+
  scale_color_viridis_d()+
  facet_wrap(~pwr,scales='free')

expand_grid(data.frame(t(colMeans(sdat))) %>% select(-lrc),
            lrc=seq(-50,10,length.out = 10),
            pred_days=seq(1,2000,length.out=2000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table() %>% 
  ggplot(data=.,aes(pred_days,pred,color=factor(lrc)))+
  geom_line()+
  scale_color_viridis_d()+
  facet_wrap(~lrc,scales='free')

# Should pwr>0 ?: yes!
sdat$pwr %>% hist




expand_grid(merge(sdat[pwr>=0],nvis, by='id') %>% 
              filter(vc!=25) %>%
              filter(vc %in% c(2,3,4,5,11)) %>% 
              filter(is.na(vc)==FALSE) %>% 
              sample_n(100), 
            pred_days=seq(1,2000,length.out=2000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table() %>% 
  ggplot(data=.,aes(pred_days,pred,color=pwr,group=id))+
  geom_line()+
  scale_color_viridis_c()
  # facet_wrap(~lrc,scales='free')

