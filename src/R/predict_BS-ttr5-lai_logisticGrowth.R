library(tidyverse); 
library(mgcv); 
library(data.table); 
library(lubridate)
library(yardstick)

# logistic growth function fits
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-1mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-06-19 17:33:16.parquet")
fits <- fits[isConv==T]
fits <- fits[L0<K][,fire_month := month(date_fire1)][fire_month %in% c(9,10,11,12,1,2)] %>% 
  .[,fire_year := year(date_fire1 - months(3))] %>%
  .[r<0.05]

d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet") %>% 
  as.data.table()
fits <- merge(fits[,.(id,K,L0,r,r2,rmse,fire_year)],
      d_rf, by='id')
fits <- fits %>% rename(dom_sp = dom_sp.x, 
                        x = x.x, 
                        y = y.x)
fits[,dom_sp_fac := factor(dom_sp)]
fits[,fy := factor(fire_year)]


# # species 
# dom <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-SEOZcoastal-grid.parquet")
# dom <- dom[,.(x,y,id,dom_sp)]
# fits <- merge(fits,dom,by=c("x","y","id"))

# core data for BS preds
dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet')
dpreds <- dpreds %>% 
  mutate(fire_month_f = factor(month(date_fire1,abbr = TRUE,label=T), 
                               levels=c("Sep","Oct","Nov","Dec","Jan","Feb"), 
                               ordered = TRUE), 
         dom_sp_fac = factor(dom_sp))

# Extract the min LAI following the fire
dlai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
o1 <- dlai[year==2003][,.(l_amp = diff(range(slai_u))), by=.(id)]
dlai <- dlai[date>=ymd("2019-09-01")][date<=ymd("2020-05-01")][id %in% dpreds$id][
  ,.(min_lai_anom = min(slai_anom_3mo,na.rm=TRUE), 
     lai_yr_sd = unique(lai_yr_sd)),by=.(id)]

dpreds <- merge(dlai, dpreds,by='id')
dpreds[,est_L0 := malai+min_lai_anom]
dpreds[,est_L0 := ifelse(est_L0<0,0.1,est_L0)]
dpreds <- dpreds %>% filter(dom_sp_fac %in% unique(fits$dom_sp_fac)) # problem: losing the rare species to predict!!!

# hist( (exp(fits$r*100))**0.05)

b1 <- bam(r ~ s(I(L0/K),k=5,bs='cs') +
              s(precip_anom_12mo, map, k=5)+
              s(vpd15_anom_12mo, mavpd15, k=5)+
              s(dom_sp_fac, bs='re') + 
              s(fy, bs='re')+
              fire_month_f, 
              family=Gamma(link='log'),  
              select=TRUE,
              data=fits)
summary(b1)
plot(b1,scheme=2,pages = 1)

dpreds %>% dim

dpreds2 <- dpreds %>% 
  select(est_L0, malai, dom_sp_fac,id,x,y,lai_yr_sd,fire_month_f, 
         precip_anom_12mo, map, 
         vpd15_anom_12mo, mavpd15, 
         dom_sp_fac) %>% 
  rename(L0 = est_L0, 
         K = malai) %>% 
  mutate(pred_r = predict(b1, 
                          newdata=., 
                          newdata.guaranteed = TRUE, 
                          exclude=c("s(fy)"),
                          type='response'))

dpreds2 <- expand_grid(dpreds2 %>% select(L0,K,pred_r,id,x,y,lai_yr_sd), 
            post_days = seq(30,3000, length.out=30)) %>% 
  mutate(pred_L = K/(1+((K-L0)/L0)*exp(-pred_r*post_days)))


dpreds3 <- dpreds2 %>% 
  # filter(id %in% sample(dpreds2$id,10000)) %>% 
  mutate(recov = ifelse(pred_L > (K-0.25*lai_yr_sd), T, F)) %>% 
  filter(recov ==T) %>% 
  group_by(id,x,y) %>% 
  summarize(pred_ttr5_lai = min(post_days)) %>% 
  ungroup()


dpreds3 %>% 
  ggplot(data=.,aes(x,y,fill=pred_ttr5_lai/365))+
  geom_raster()+
  coord_sf()+
  scale_fill_viridis_c(option='B', 
                       end=0.95,
                       limits=c(1,5), 
                       oob=scales::squish)


# Model evaluation -------------------------------------------------------------
deval <- fits %>%
  # select(L0, malai, dom_sp_fac,id,x,y,lai_yr_sd,ttr5_lai) %>% 
  # rename(# L0 = est_L0,
  #        K = malai) %>%
  mutate(pred_r = predict(b1, newdata=., newdata.guaranteed = TRUE, type='response'))

deval <- expand_grid(deval %>% select(L0,K,pred_r,id,x,y,fire_year, lai_yr_sd,ttr5_lai), 
            post_days = seq(365,6000, by=30)) %>% 
  mutate(pred_L = K/(1+((K-L0)/L0)*exp(-pred_r*post_days))) %>% 
  mutate(recov = ifelse(pred_L > (K-0.25*lai_yr_sd), T, F)) %>% 
  filter(recov ==T) %>% 
  group_by(id,x,y) %>% 
  summarize(pred_ttr5_lai = min(post_days), 
            ttr5_lai = unique(ttr5_lai), 
            fire_year = unique(fire_year)) %>% 
  ungroup()

rsq_trad_vec(deval$ttr5_lai, 
             deval$pred_ttr5_lai)
rmse_vec(deval$ttr5_lai, 
             deval$pred_ttr5_lai)
rmse_vec(deval$ttr5_lai, 
         rep(mean(deval$ttr5_lai,na.rm=T),dim(deval)[1]))

deval %>% 
  ggplot(data=.,aes(pred_ttr5_lai,ttr5_lai))+
  geom_bin2d()+
  geom_smooth(method='lm')+
  geom_abline(col='red')+
  scale_fill_viridis_c()+
  facet_wrap(~fire_year)



fits %>% ggplot(data=.,aes(-slai_anom,(malai-L0)))+geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')

# SCRATCH ----------------------------------------------------------------------
fits %>% ggplot(data=., aes(precip_anom_12mo/map, ttr5_lai))+
  geom_smooth()

# fits[L0<K][,fire_month := month(date_fire1)][fire_month %in% c(9,10,11,12,1,2)] %>% 
#   .[,fire_year := year(date_fire1 - months(3))] %>% 
#   ggplot(aes(L0/K,r,color=fire_month,group=fire_month))+
#   geom_smooth()+
#   facet_wrap(~fire_year)
fits %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_raster()+
  coord_equal()+
  scale_fill_viridis_c(option='B', 
                       limits=c(0,0.01), 
                       oob=scales::squish)


b1 <- bam(log(malai/r) ~ 
            te(malai,L0,x,y)+
            # te(x,y, k=15)+
            fire_month_f, 
          family=Gamma(link='log'),
          select=F,
          discrete = TRUE,
          data=fits)
summary(b1)
yhat <- predict(b1,type='response')

cor(log(fits$malai/fits$r), yhat)**2
rsq_trad_vec(log(fits$malai/fits$r), estimate = yhat)
rsq_trad_vec(log(fits$malai/fits$r), estimate = yhat)

rhat <- fits %>% 
  mutate(pred = predict(b1,type='response')) %>% 
  mutate(pred = exp(pred)) %>% 
  mutate(pred = pred/malai) %>% 
  mutate(pred = 1/pred) %>% 
  pull(pred)
rsq_trad_vec(fits$r, estimate = rhat)

plot(b1, scheme=2,pages=1)

fits %>% 
  mutate(pred = predict(b1,type='response')) %>% 
  mutate(pred = exp(pred)) %>% 
  mutate(pred = pred/malai) %>% 
  mutate(pred = 1/pred) %>% 
  pull(pred) %>% hist(breaks=100)
fits$r %>% hist(breaks=100)

fits %>% 
  ggplot(data=.,aes(log(malai/r)))+
  geom_histogram(bins=100)


deval <- fits %>%
  # select(L0, malai, dom_sp_fac,id,x,y,lai_yr_sd,ttr5_lai) %>% 
  # rename(# L0 = est_L0,
  #        K = malai) %>%
  mutate(pred_r = exp(predict(b1, newdata=., newdata.guaranteed = TRUE, type='response'))) %>% 
  mutate(pred_r = pred_r/malai) %>% 
  mutate(pred_r = 1/pred_r) 

deval <- expand_grid(deval %>% select(L0,K,r,pred_r,id,x,y,fire_year, lai_yr_sd,ttr5_lai), 
                     post_days = seq(365,6000, by=30)) %>% 
  mutate(pred_L = K/(1+((K-L0)/L0)*exp(-r*post_days))) %>% 
  mutate(recov = ifelse(pred_L > (K-0.25*lai_yr_sd), T, F)) %>% 
  filter(recov ==T) %>% 
  group_by(id,x,y) %>% 
  summarize(pred_ttr5_lai = min(post_days), 
            ttr5_lai = unique(ttr5_lai), 
            fire_year = unique(fire_year)) %>% 
  ungroup()

rsq_trad_vec(deval$ttr5_lai, 
             deval$pred_ttr5_lai)
rmse_vec(deval$ttr5_lai, 
         deval$pred_ttr5_lai)
rmse_vec(deval$ttr5_lai, 
         rep(mean(deval$ttr5_lai,na.rm=T),dim(deval)[1]))

deval %>% 
  ggplot(data=.,aes(pred_ttr5_lai,ttr5_lai))+
  geom_bin2d()+
  geom_smooth(method='lm')+
  geom_abline(col='red')+
  scale_fill_viridis_c()+
  facet_wrap(~fire_year)






deval <- expand_grid(deval %>% select(L0,K,r,pred_r,id,x,y,fire_year, lai_yr_sd,ttr5_lai), 
                     post_days = seq(365,6000, by=30)) %>% 
  mutate(pred_L = K/(1+((K-L0)/L0)*exp(-r*post_days))) %>% 
  mutate(recov = ifelse(pred_L > (K-0.25*lai_yr_sd), T, F)) %>% 
  filter(recov ==T) %>% 
  group_by(id,x,y) %>% 
  summarize(pred_ttr5_lai = min(post_days), 
            ttr5_lai = unique(ttr5_lai), 
            fire_year = unique(fire_year)) %>% 
  ungroup()

deval %>% 
  ggplot(data=.,aes(min_nbr_anom,r))+
  geom_smooth(formula=y~s(x,k=5,bs='cs'))

deval %>% 
  ggplot(data=.,aes(L0/K,ttr5_lai))+
  geom_smooth(formula=y~s(x,k=5,bs='cs'))

deval <- merge(deval,o1,by='id')

bam(r~s(min_nbr_anom,k=5,bs='cs'),data=deval) %>% summary
bam(r~s(I(L0/malai),k=5,bs='cs') + 
      s(I(l_amp/malai), k=5,bs='cs')+
    +s(lai_yr_sd,k=5,bs='cs'),data=deval) %>% summary

fits %>% ggplot(data=.,aes(L0,K,color=r))+
  geom_point(size=0.1)+
  scale_color_viridis_c(option='H')

fits %>% ggplot(data=.,aes(L0,K))+
  geom_density2d_filled()


fits[r2>0.3] %>% ggplot(data=.,aes(L0,K,color=r))+
  geom_point(size=0.1)+
  scale_color_viridis_c(option='H')+
  facet_wrap(~cut_interval(post_precip_anom_12mo,4))

fits[r2>0.3] %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  scale_fill_viridis_b(option='H')+
  # coord_equal()+
  facet_grid(~cut_width(y,4),scales = 'free')


fits[r2>0.3] %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  geom_point(data=. %>% .[r>0.015], aes(x,y),col='red',size=0.25)+
  scale_fill_viridis_b(option='H')


fits[,`:=`(x1=round(x/0.125)*125, 
           y1=round(y/0.125)*125)][,.(val = mean(r,na.rm=T)), 
                        by=.(x1,y1)] %>% 
  ggplot(data=.,aes(x1,y1,fill=val))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='H')

fits[,.(val = mean(post_precip_anom_12mo)), by=cut_number(r,10)]

fits %>% #[post_precip_anom_12mo %between% c(-500,500)] %>% 
  ggplot(data=.,aes(post_precip_anom_12mo/map,r))+
  geom_smooth(formula=y~s(x,k=6,bs='cs'))+
  geom_rug()


# Estimate species level TTR from r and malai ----------------------------------
fits[,fy := factor(fire_year)]
b3 <- bam(log(malai/r) ~ 
            dom_sp_fac + 
            s(fy, bs='re') + 
            s(fire_month_f,bs='re'), 
          data=fits)
plot(b3)
betas <- coef(b3) %>% as.data.table(keep.rownames = T) %>% 
  set_names(c("effect","estimate"))
betas <- betas[,effect:=str_remove(betas$effect, "dom_sp_fac")]
betas2 <- betas[str_detect(betas$effect, c("Angophora","Corymbia","Eucalyptus")),]
g2 <- dlai[date==ymd("2003-01-01")][,.(x,y,malai,id)]
g3 <- merge(betas2 %>% rename(dom_sp=effect), dom, by='dom_sp')
g3 <- merge(g3,g2,by=c('x','y','id'))

g3 %>% 
  .[,p1 := (exp(estimate+betas[effect=="(Intercept)"]$estimate))] %>% 
  .[,p2 := p1/malai] %>% 
  .[,p3 := 1/p2] %>% 
  ggplot(data=.,aes(x,y,fill=p1/365))+
  geom_tile()+
  coord_sf()+
  scale_fill_viridis_c(option='B')

g3 %>% 
  .[,p1 := (exp(estimate+betas[effect=="(Intercept)"]$estimate))] %>% 
  .[,p2 := p1/malai] %>% 
  ggplot(data=.,aes(y=dom_sp,
                    x=p2,
                    group=dom_sp))+
  geom_boxplot(outlier.colour = NA)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim = c(0,2000), 
                  expand=F)+
  labs(y=NULL, 
       x='Time to Recover (days)')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave(filename = "figures/gam_ttr-estimate_by-species_from_r-malai.png", 
       width=15, 
       height=25,
       units='cm')

# Estimate species level TTR from r and malai ----------------------------------
library(mgcv)
b3 <- bam(ttr5_lai ~ 
            dom_sp_fac + 
            s(fy, bs='re') + 
            s(fire_month_f,bs='re'), 
          data=fits)
plot(b3)
betas <- coef(b3) %>% as.data.table(keep.rownames = T) %>% 
  set_names(c("effect","estimate"))
betas <- betas[,effect:=str_remove(betas$effect, "dom_sp_fac")]
betas2 <- betas[str_detect(betas$effect, c("Angophora","Corymbia","Eucalyptus")),]
g2 <- dlai[date==ymd("2003-01-01")][,.(x,y,malai,id)]
g3 <- merge(betas2 %>% rename(dom_sp=effect), dom, by='dom_sp')
g3 <- merge(g3,g2,by=c('x','y','id'))

g3 %>% 
  .[,p1 := (estimate+betas[effect=="(Intercept)"]$estimate)] %>% 
  # .[,p2 := p1/malai] %>% 
  # .[,p3 := 1/p2] %>% 
  ggplot(data=.,aes(x,y,fill=p1/365))+
  geom_tile()+
  coord_sf()+
  scale_fill_viridis_c(option='B')

g3 %>% 
  .[,p1 := ((estimate+betas[effect=="(Intercept)"]$estimate))] %>% 
  .[,p2 := p1/malai] %>% 
  ggplot(data=.,aes(y=dom_sp,
                    x=p1,
                    group=dom_sp))+
  geom_boxplot(outlier.colour = NA)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim = c(0,2000), 
                  expand=F)+
  labs(y=NULL, 
       x='Time to Recover (days)')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave(filename = "figures/gam_ttr-estimate_by-species_from_r-malai.png", 
       width=15, 
       height=25,
       units='cm')

fits %>% 
  ggplot(data=.,aes(malai/r, ttr5_lai))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')+
  scale_x_continuous(limits=c(0,5000))


j <- fits[id==960299]
# j <- fits[sample(.N,10)]
j1 <- expand_grid( j %>% 
               select(L0,K,r,id,x,y,fire_year, lai_yr_sd,ttr5_lai,malai), 
             post_days = seq(1,6000, by=30)) %>% 
  mutate(pred_L = K/(1+((K-L0)/L0)*exp(-r*post_days))) %>% 
  mutate(recov = ifelse(pred_L > (K-0.25*lai_yr_sd), T, F)) %>% 
  filter(id==960299)  
expand_grid( j %>% 
              select(L0,K,r,id,x,y,fire_year, lai_yr_sd,ttr5_lai,malai), 
            post_days = seq(1,6000, by=30)) %>% 
  mutate(pred_L = K/(1+((K-L0)/L0)*exp(-r*post_days))) %>% 
  mutate(recov = ifelse(pred_L > (K-0.25*lai_yr_sd), T, F)) %>% 
  ggplot(data=.,aes(post_days, pred_L, group=id,color=id))+
  geom_line()+
  geom_point(aes(ttr5_lai,malai,color=id,group=id))+
  scale_color_viridis_c(option='H')+
  facet_wrap(~id)
j2 <- dlai[id%in%j$id]

dlai[id==960299][date<=ymd("2009-01-01")][date>=ymd("2002-06-01")] %>% 
  ggplot(data=.,aes(date, slai_anom))+
  geom_line()+
  geom_vline(data=j[id==960299], 
             aes(xintercept=date_fire1), 
             col='red')+
  geom_vline(data=j[id==960299], 
             aes(xintercept=date_fire1+days(ttr5_lai)), 
             col='blue')+
  geom_hline(data=j[id==960299], 
             aes(yintercept=malai), 
             col='grey')+
  geom_line(data=j1, 
            aes(x=ymd("2002-10-01")+post_days, 
                y=pred_L), 
            col='orange')

j[id==960299]$date_fire1




mdat[id==960299] %>% 
  ggplot(data=.,aes(post_days,slai_1mo))+
  geom_line()

fits[r2>0.5][ttr5_lai > 366][sample(.N,1000)] %>% 
  ggplot(data=.,aes((malai-0.25*lai_yr_sd)/r, ttr5_lai))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')

fits[r2>0.5][ttr5_lai > 366][sample(.N,1000)] %>% 
  ggplot(data=.,aes(r, ttr5_lai))+
  geom_point()+
  geom_smooth()

fits[r2>0.5][sample(.N,10000)] %>% 
  ggplot(data=.,aes(K,malai))+geom_point()+
  geom_smooth()+
  geom_abline(col='red')+
  coord_equal()

dlai$date %>% min
fits$date_fire1 %>% summary
fits[r2>0.5][sample(.N,100000)] %>% 
  ggplot(data=.,aes(L0/K,r))+
  ggpointdensity::geom_pointdensity()+
  geom_smooth()+
  scale_color_viridis_c(option='H')
  # geom_abline(col='red')+

fits$r2 %>% hist

jj <- mdat[id%in%fits[r>=0.03]$id]


fn_pred <- function(d){
  with(d,
  out <- K/((1+((K-L0)/L0)*exp(-r*post_days)))
  )
  return(out)
}

jj[jj[, .I[date >= date_fire1], by = .(id)]$V1] %>% 
  .[id %in% sample(fits[r>=0.03]$id, 10)] %>% 
  .[post_days <= (ttr5_lai+365)] %>% 
  merge(., fits[,.(id,r,L0,K)]) %>% 
  mutate(pred_L =  K/((1+((K-L0)/L0)*exp(-r*post_days)))) %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(post_days, slai_1mo))+
  geom_line()+
  geom_line(aes(post_days, pred_L), col='red')+
  geom_vline(aes(xintercept=ttr5_lai),lty=1,col='navy')+
  geom_hline(aes(yintercept=K),col='blue',lty=3)+
  geom_hline(aes(yintercept=malai),col='red',lty=3)+
  facet_wrap(~id, scales='free')

fits[r>=0.03]$ttr5_lai %>% hist

df <- data.frame(x=1:10, 
                 y=exp(-0.4*(1:10)+0.3*rnorm(10)))
fit <- nls(y~exp(-a*x), df, list(a=-1))
ggplot(data=df, aes(x,y))+
  geom_point()+
  geom_function(fun=~predict(fit, data.frame(x=.x)))

fn_pred(expand_grid(fits[1,],post_days=seq(1:100)))


fits %>% 
  as_tibble() %>% 
  sample_n(1000) %>% 
  mutate(ttr_ = log(1/(1+(K-L0)/K))/-r) %>% 
  ggplot(data=.,aes(ttr_,ttr5_lai))+
  geom_point()+
  geom_smooth()+
  geom_abline(col='red')

junk <- fits[sample(.N,1000)] %>% as_tibble()
expand_grid(post_days=1:2000, 
            junk) %>% 
  mutate(pred_L = K/((1+((K-L0)/L0)*exp(-r*post_days))) ) %>% 
  mutate(ttr_ = -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r) %>% 
  filter(pred_L >= 0.95*K) %>% 
  group_by(id) %>% 
  summarize(est_ttr = min(post_days), 
            ttr_ = unique(ttr_)) %>% 
  ungroup() %>% 
  ggplot(data=., aes(ttr_, est_ttr))+
  geom_point()+
  geom_abline()

fits[r2>0.5][L0<K][K>malai][sample(.N,1000)] %>% as_tibble() %>% 
  mutate(ttr_ = -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r) %>% 
  ggplot(data=., aes(ttr_+180, ttr5_lai))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')


fits %>% mutate(val = K>malai) %>% as_tibble() %>% 
  filter(val ==F) %>% 
  select(L0, K, r, malai) %>% 
  pull(r) %>% hist



fits[r2>0.55][r<0.05][L0<K][K>malai] %>%
  .[fire_year>=2001 & fire_year<=2014] %>% 
  as_tibble() %>% 
 ggplot(data=.,aes(L0/K,r,color=factor(fire_year)))+
  # geom_point()+
  geom_smooth(se=F)

fits[r2>0.75][L0<K][K>malai][sample(.N,10000)] %>% as_tibble() %>% 
  ggplot(data=.,aes(L0/malai,ttr5_lai))+
  # geom_point()+
  geom_smooth()

fits[,.(nobs=.N),by=fire_year][,rank:=frank(-nobs)][rank<=10]
fits[r2>0.75][L0<K][K>malai][,`:=`(fire_year = year(date_fire1-months(3)))] %>% 
  .[fire_year>=2001 & fire_year<=2014] %>% 
  .[fire_year %in% fits[,.(nobs=.N),by=fire_year][,rank:=frank(-nobs)][rank<=8]$fire_year] %>% 
  ggplot(data=.,aes(L0/malai,ttr5_lai,color=factor(fire_year)))+
  geom_smooth(se=F)+
  scale_color_viridis_d(option='H')

fits %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  geom_point(data=. %>% filter(r>0.03), 
             aes(x,y),col='red', 
             size=0.1)+
  scale_fill_viridis_b(option='B')+
  coord_equal()

