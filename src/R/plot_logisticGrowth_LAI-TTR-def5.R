library(tidyverse);
library(stars);
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)
library(data.table); 
library(lubridate) # load AFTER data.table
library(mgcv)

# Data import ---------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                           col_select = c("x","y","id","date","slai","slai_anom_12mo","malai"))
dat[,`:=`(slai_12mo = slai_anom_12mo+malai)]
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-1mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-06-19 17:33:16.parquet")
fits <- fits[isConv==TRUE][r2>0]
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")

# calculate the rolling metrics ------------------------------------------------
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_24mo := frollapply(precip_anom,FUN=mean,
                                                              n = 24,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_36mo := frollapply(precip_anom,FUN=mean,
                                                              n = 36,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_3mo := frollapply(precip_anom,FUN=mean,
                                                              n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_3mo := frollapply(precip,FUN=sum,
                                                n = 3,fill = NA,align='right'), by=.(x,y)]


#! Only fitting locations where the recovery was at least one year
sdat <- sdat[is.na(ttr5_lai)==FALSE][date_fire1<ymd('2015-01-01')][ttr5_lai>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_fire1,ttr5_lai)], 
               by=c("x","y","id"))

# mdat <- ssdat %>% 
#   lazy_dt() %>%
#   mutate(recovery_date = date_fire1+days(ttr5_lai)) %>% 
#   group_by(x,y,id) %>% 
#   filter(date >= date_fire1) %>% 
#   filter(date <= recovery_date + years(1)) %>% 
#   ungroup() %>%
#   mutate(post_days = as.double(date - date_fire1)) %>% 
#   as.data.table() 

rm(dat); gc(full=TRUE)
fits <- merge(fits,d_soil,by='id')


#  Attach AWAP pixel id to VI ------------------------------------
coords_vi <- lazy_dt(fits) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_awap <- unique(clim[,.(x,y)])
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
st_crs(coords_awap_sf) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
coords_awap <- coords_awap %>% mutate(idx_awap = row_number()) %>% as.data.table()
gc(full=TRUE)
coords_vi <- coords_vi %>% st_drop_geometry() %>% as.data.table()
coords_vi$idx_awap <- coords_awap[nn_coords$nn.idx,]$idx_awap
gc(full=TRUE)


# merges
gc(full=TRUE)
clim <- merge(clim,coords_awap,by=c('x','y'))
gc(full=TRUE)
fits <- merge(fits, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
clim <- clim[idx_awap %in% unique(fits$idx_awap)]


# clim[,.(idx_awap,date,precip,precip_anom,map)]
jj <- merge(fits[,date_recovery:=(date_fire1+days(ttr5_lai))], 
            clim[,.(idx_awap,date,precip,
                    precip_anom,precip_anom_3mo, precip_anom_12mo,precip_anom_24mo,precip_anom_36mo,
                    vpd15_anom_3mo, vpd15_anom_12mo,
                    ppet_anom_12mo,
                    pet,map,mapet,mappet, pet_anom_12mo,
                    matmax,matmin,mavpd9,
                    mavpd15,vpd15,vpd15_anom)], 
            by='idx_awap', allow.cartesian = TRUE)
jj <- jj %>% rename(date=date.y, junk=date.x) %>% as.data.table()
jj <- jj[jj[, .I[date >= date_fire1 & date <= date_recovery], 
                by = .(idx_awap)]$V1]

f_cwd <- function(cwd_et,precip,et, min_threshold){
  # No reset during the wettest month of the year
  for(i in seq(2,length(precip))){
    
    cwd_et[i] <-  min(0, cwd_et[i-1] + (precip[i]) - max(et[i],1, na.rm=T), na.rm=T)
    cwd_et[i] <- ifelse(cwd_et[i] < min_threshold, min_threshold, cwd_et[i])
    
  }
  cwd_et
}

f_cwd_r7 <- function(cwd_et,precip,et, dates){
  #reset on month 7
  for(i in seq(2,length(precip))){
    
    cwd_et[i] <-  min(0, cwd_et[i-1] + (precip[i]) - max(et[i],1, na.rm=T), na.rm=T)
    cwd_et[i] <- ifelse(month(dates[i])==7, 0, cwd_et[i])
    
  }
  cwd_et
}

jj <- jj[,`:=`(cwd = 0)] #
jj <- jj[order(id,date)]
jj <- jj[(date>=date_fire1 & date<=date_recovery)]
jj <- jj[,`:=`(cwd = f_cwd_r7(cwd_et=cwd, 
                             precip=precip, 
                             et=pet,
                             date=date)), 
           by=.(id)]



# PLOTS -------------------------------------------------------------------
# Composite plot: ttr5~L0/K by year ----------------
pan_b <- fits[year %in% fits[,.(nobs = .N),by='year'][,rank := frank(-nobs)][order(rank)][rank<=10]$year] %>% 
  # sample_n(10000) %>% 
  as_tibble() %>% 
  mutate(fire_year = factor(year)) %>% 
  ggplot(data=.,aes(ldk, ttr5_lai,color=fire_year,group=fire_year))+
  # geom_point()+
  geom_smooth(method='gam',
              formula=y~s(x,k=5,bs='cs'),
              method.args=list(select=TRUE))+
  scale_color_viridis_d(option='H')+
  labs(x=expression(paste(L[0]/K)), 
       y='Time to Recover (days)',
       color='Fire year')+
  coord_cartesian(xlim=c(0,0.75), 
                  expand=F)+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom'); pan_b

pan_c <- expand_grid(tibble(K=c(3,3,3),
                   r=c(0.01, 0.002, 0.0075),
                   L0=c(0.3,1.5,2)), 
            post_days = floor(seq(1,3000,length.out=100))) %>% 
  mutate(reduction = round(100*(1 - L0/K),digits=2)) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days)) ) %>% 
  mutate(recovered = ifelse(between(pred, 0.975*K, 0.98*K),pred,NA)) %>% 
  ggplot(data=.,aes(post_days,pred,
                    group=paste(r,L0,K), 
                    color=factor(r)))+
  geom_line()+
  # geom_point(aes(post_days,recovered))+
  labs(x='Days after fire',
       y='LAI', 
       color=expression(paste(bolditalic(r))))+
  scale_color_brewer(palette='Set2')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom'); pan_c

pan_a <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-26 15:23:33.parquet") %>% 
  as.data.table() %>% 
  .[isConv==TRUE] %>% 
  .[r<0.05] %>% 
  .[r2>0] %>% 
  .[L0<K] %>% 
  # .[r2>0.75] %>% 
  ggplot(data=.,aes(L0/K,r,color=ttr5_lai))+
  geom_point(alpha=0.5,size=0.5)+
  geom_smooth(fullrange=T,
              aes(L0/K,r),
              color='#cf0000',
              weight=3,
              level=0.999)+
  scale_color_viridis_c(option='F',
                        limits=c(365,3000),
                        direction = 1,
                        end = 0.85,
                        oob=scales::squish
  )+
  scale_x_continuous(limits=c(0,1),expand=c(0,0))+
  scale_y_continuous(limits=c(0,0.02),expand=c(0,0))+
  labs(y=expression(paste("r  (",m**2~day**-1,")")), 
       x=expression(paste(L[0]/K,"  (",m**2/m**2,")")), 
       color="TTR (days)")+
  theme_linedraw() 

library(patchwork)
pan_a/(pan_b|pan_c)+plot_annotation(tag_levels = 'a',
                                    tag_prefix = '(',
                                    tag_suffix = ')')
ggsave(pan_a/(pan_b|pan_c)+plot_annotation(tag_levels = 'a',
                                           tag_prefix = '(',
                                           tag_suffix = ')'), 
       filename = 'figures/multipanel_logfitParams_r-ldk_ttr-ldk_logFigDiagram.png', 
       width=20, 
       height=20,
       units='cm',
       dpi=300)





jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(11,12,1,2,3)][date==(date_fire1+years(1))] %>% 
  ggplot(data=.,aes(map,r,color=vc_name))+
  geom_point(alpha=0.1)+
  # scale_color_viridis_c(option='A')+
  geom_smooth(method='lm')+
  facet_wrap(~cut_number(L0/K,n=10),scales='free')

# jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(10,11,12,1,2,3)][date==(date_fire1+years(0))] %>% 
#   lm(K~malai,data=.) %>% summary
#   ggplot(data=.,aes(malai,K))+
#   geom_smooth(method='lm')+
#   geom_abline()

b_K <- jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(10,11,12,1,2,3)][date==(date_fire1+years(0))] %>% 
  bam(malai~
        # te(pH,k=5)+
        te(matmax,matmin,mapet, k=5)+
        te(map, mappet,mavpd15,k=5),
        # s(map, k=5)+
        # s(mavpd15,k=5)+
        # s(mappet,k=5)+
        # s(pH,k=5)+
        # s(der,k=5)+
        # s(elevation,k=5),
        family=Gamma(link='log'),
      data=., 
      select=T, 
      discrete = T)
summary(b_K)
plot(b_K,scheme = 2)

b_L0 <- jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(10,11,12,1,2,3)][date==(date_fire1)] %>% 
  .[,delta:=(K-L0)] %>% 
  bam(delta ~
        # s(mappet,k=5)+
        # s(I(precip_anom_3mo/map),k=5)+
        # s(vpd15_anom_3mo,k=5)+
        # s(I(slai_anom_12mo/malai), k=5)+
        # s(precip_anom_12mo, k=5)+
        # s(pet_anom_12mo,k=5)+
        # s(I(precip_anom_24mo/map), k=5)+
        # s(precip_anom_36mo, k=5)+
        # s(I(ppet_anom_12mo/mappet),k=5)+
        # s(I(vpd15_anom_12mo/mavpd15),k=5)+
        te(vpd15_anom_3mo, precip_anom_12mo,
           malai,lai_yr_sd, k=5),
        # s(soc,k=5)+
        # s(pH,k=5)+
        # s(tpi,k=5)+
        # s(mappet,k=5),
        family=Gamma(link='log'),
      data=., 
      select=T, 
      discrete = T)
summary(b_L0)
plot(b_L0,scheme = 2)

b_r <- jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(10,11,12,1,2,3)][date==(date_fire1+years(0))] %>% 
  bam(r~s(I(L0/K), k=5)+
        s(lai_yr_sd,k=5)+
        # s(malai,k=5)+
        # s(mavpd15,k=5)+
        # s(elevation,by=aspect,k=5)+
        # s(awc,k=5)+
        s(silt,k=5)+
        s(pH,k=5)+
        s(mappet,k=5),
        # s(map,k=5),
      family=scat(),
      data=., 
      select=T, 
      discrete = T)
summary(b_r)
plot(b_r,scheme = 2)

b_r2 <- jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(11,12,1,2,3)][date==(date_fire1+days(ttr5_lai))] %>% 
  bam(r~te(I(map+precip_anom_12mo),I(L0/K)), 
      family=scat(),
      data=., 
      select=T, 
      discrete = T)
summary(b_r2)
plot(b_r2,scheme=2)
mgcViz::getViz(b_r2) %>% plot


b_r3 <- jj[r<0.1][L0<K][vc%in%c(2,3,5)][month(date_fire1)%in%c(11,12,1,2,3)][date==(date_fire1+days(ttr5_lai))] %>% 
  bam(r~te(map,I(L0/K), k=5), 
      family=scat(),
      data=., 
      select=T, 
      discrete = T)
library(mgcViz)
plotSlice(sm(getViz(b_r3),1),fix=list('map'=c(500,1000,1500)))+
  l_fitRaster()
plot(mgcViz::getViz(b_r3))+l_fitRaster()+l_fitContour()+
  l_rug()+
  scale_fill_gradient2(limits=c(-0.01,0.01),oob=scales::squish)

b0_map <- jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(11,12,1,2,3)][date==(date_fire1+years(1))] %>% 
  bam(r~te(map), 
      family=scat(),
      data=., 
      select=T, 
      discrete = T)
summary(b0_map)
b0_delta <- jj[r<0.1][L0<K][(L0/K)<0.5][vc%in%c(2,3,5)][month(date_fire1)%in%c(11,12,1,2,3)][date==(date_fire1+years(1))] %>% 
  bam(r~te(I(L0/K)), 
      family=scat(),
      data=., 
      select=T, 
      discrete = T)
summary(b0_delta)

jj %>% group_by(id) %>% summarize(ttr=unique(ttr5_lai), 
                                  val = median(cwd,na.rm=TRUE)) %>% as.data.table()

expand_grid(expand_grid(K=c(2,4,6), 
                        r=c(0.0015,0.0025,0.005),
                        reduction=c(0.95,0.5,0.25)) %>% 
              mutate(L0=K*(1-reduction)) %>% 
              filter(L0<0.95*K),
            post_days=floor(seq(1,5000,length.out=100))) %>% 
  # left_join(., sdat, by='id') %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days)) ) %>% 
  mutate(recovered = ifelse(between(pred, 0.975*K, 0.98*K),pred,NA)) %>% 
  ggplot(data=.,aes(post_days,pred,
                    group=paste(r,L0,K), 
                    color=factor(reduction)))+
  geom_line()+
  geom_point(aes(post_days,recovered))+
  labs(x='days after fire',
       y='LAI', 
       color='Fire LAI reduction (%)')+
  scale_color_brewer(palette='Set2')+
  facet_grid(r~K,labeller = label_both)+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom')
ggsave(filename = "figures/figure_logistic_growth_function_diagram.png", 
       width=25*0.65,
       height=25*0.55, 
       units='cm',
       dpi=300)

merge(fits, d_soil,by='id')[r<0.1][L0>0][L0<K][r2>0.5][sample(.N,5000)] %>% 
  .[,'delta_L':=(L0/K)] %>% 
  .[,.(elevation,der,pto,sand,silt,pH,K,delta_L,r)] %>% 
  GGally::ggpairs(mapping=list(alpha=0.01))+
  theme_linedraw()

fits[between(fire_doy,300,366) || between(fire_doy,0,90)] %>% 
  .[,fire_year := year(date_fire1-months(3))] %>% 
  ggplot(data=.,aes(ttr5_lai))+
  stat_bin()+
  geom_vline(aes(xintercept=1000),color='yellow')+
  geom_vline(aes(xintercept=2000),color='red')+
  facet_wrap(~fire_year,scales='free_y')

jj[date>=date_fire1][date<=date_recovery][r<0.1][vc %in% c(2,3,5)][,
                       .(r=mean(r), 
                         pu=mean(precip), 
                         ptot = sum(precip,na.rm=T),
                         ttr=mean(ttr5_lai)), 
                       by=.(x,y,idx_awap,vc_name)] %>% 
  # [sample(.N,20000)] %>% 
  ggplot(data=.,aes(pu,r,color=vc_name))+
  geom_point(size=1,alpha=0.05)+
  geom_smooth(method='lm')+
  scico::scale_color_scico_d(end=0.9)+
  labs(x='Mean monthly precip during recovery period', 
       y=expression(paste(r~"(",m**2/day,")")), 
       color="Veg. Class")+
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = "figures/figure_lai-logGrowthFit_r_meanP-during-recovery.png", 
       width=25*1,
       height=25*0.55, 
       units='cm',
       dpi=300)

jj[date>=date_fire1][date<=date_recovery][r<0.1][L0<K][vc %in% c(2,3,5)][,.(r=mean(r),
         delta = mean(K-L0),
         val =mean(precip*12/map,na.rm=T),                     
         p_anom = mean(12*precip_anom/map,na.rm=T),                     
         ttr=mean(ttr5_lai)), 
      by=.(x,y,idx_awap,vc_name)] %>%
  .[delta > 0] %>% 
  ggplot(data=.,aes(p_anom,r,color=vc_name))+
  geom_point(size=0.2,alpha=0.2)+
  geom_smooth(method='lm')+
  scico::scale_color_scico_d(end=0.9)+
  labs(x='Average annual precip anom relative to MAP during recovery period', 
     y=expression(paste(r~"(",m**2/day,")")), 
     color="Veg. Class")+
  scale_y_continuous(limits=c(0,0.03))+
  scale_x_continuous(limits=c(-0.5,0.75))+
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = "figures/figure_lai-logGrowthFit_r_meanPanom-rel-MAP-during-recovery.png", 
       width=25*1,
       height=25*0.5, 
       units='cm',
       dpi=300)


# vec_sel <- sample(unique(fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0]$id),25)
vec_sel <- fits[isConv==TRUE][is.na(r2)==FALSE][(L0/K)<0.25][L0>0][r2>0.5][vc %in% c(2,3,5)][
  ,.SD[sample(.N,5)],by=vc_name]$id

# plot GOF of logistic growth model per grid cell ------------------------------
fits[id %in% vec_sel] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,2000,length.out=300))) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days))) %>% 
  filter(post_days <= ttr5_lai+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  mutate(r2 = format(r2,2)) %>%
  ggplot(data=.,aes(post_days,pred,group=id,color=vc_name))+
  geom_point(data=merge(mdat[id%in%vec_sel][post_days < (ttr5_lai+365)],
                        fits[,.(id,start_day)],by='id',all.x = T,all.y=F) %>% 
               .[,`:=`(post_days2 = post_days-start_day)] %>% .[post_days2 >= 0],
             inherit.aes = F,
             aes(post_days2, slai_12mo),
             alpha=0.5,size=0.5)+
  scale_color_viridis_d(option='D',end=0.8,direction = 1)+
  geom_hline(aes(yintercept=0),col='grey')+
  geom_line()+
  geom_vline(aes(xintercept=ttr5_lai),col='black',lty=3)+
  geom_hline(aes(yintercept=malai),lty=3)+
  geom_text(data=. %>% filter(post_days==min(post_days)),
            aes(2000,
                1,
                label=paste('r^2=',format(r2,digits=2))),
            color='black')+
  labs(x='Days after fire', 
       y='LAI', 
       color='NVIS Vegetation Group')+
  facet_wrap(~id)+
  # facet_grid(rows=vars(vc_name), 
  #            # cols = vars(id), 
  #            margins=F)+
  # facet_wrap(~paste0("x:",format(x,digits=4),
  #                    "  y:",format(y,digits=4)),scales = 'free_x',nrow = 3)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')

ggsave(filename = 'figures/figure_15timeseriesExample_logisticGrowth_LAI12mo-TTR-Def5.png',
       width=15*2.1,
       height=10*2.1,
       units='cm',
       dpi=350)










jj[date>=date_fire1][date<=date_recovery][r<0.1][vc %in% c(2,3,5)][,.(r=mean(r), 
      pu=mean(precip), 
      ptot = sum(precip,na.rm=T),
      ttr=mean(ttr5_lai)), 
   by=.(x,y,idx_awap)] %>% 
  ggplot(data=.,aes(ptot,r))+
  geom_point(size=0.1,alpha=0.1)+
  geom_smooth(method='glm', 
              method.args=list(family=Gamma(link='inverse')))





jj[r<0.1][vc %in% c(2,3,5)][,.(r=mean(r),
                               L0=mean(L0),
                               K=mean(K),
                               delta = mean(K-L0),
                               lai = mean(malai),
                               v_anom =mean(vpd15_anom/mavpd15,na.rm=T),                     
                               p_anom = mean(12*precip_anom/map,na.rm=T),
                               p_tot = sum(precip,na.rm=TRUE),
                               ttr=mean(ttr5_lai)), 
                            by=.(x,y,idx_awap,vc_name)][sample(.N,30000)] %>%
  .[delta > 0] %>% 
  # gam(ttr~scale(p_anom)*vc_name+scale(val)*vc_name, data=., family=Gamma(link='log')) %>% summary
  ggplot(data=.,aes(p_anom,r,color=L0/lai))+
  geom_point(alpha=0.5)+
  geom_smooth(method='lm')+
  scale_color_viridis_c(option='B',limits=c(0,1),oob=scales::squish,direction = 1,end=0.9) 
# facet_wrap(~cut_interval(lai,n=4))




fits[vc %in% c(2,3,5)][r<0.1] %>% 
  ggplot(data=.,aes(date_fire1,r,color=vc_name))+
  # geom_point(alpha=0.1)+
  geom_smooth(formula=y~s(x,bs='cs',k=10))








# Kilmore FIG -------------------------------------------------------
fits[,`:=`(fire_year = year(date_fire1-months(3)))]
fits$fire_year %>% table %>% sort


vec_sel3 <- fits[isConv==TRUE][is.na(r2)==FALSE][r2>0][vc %in% c(2,3,5)]
fits[isConv==TRUE][is.na(r2)==FALSE][L0>0][vc %in% c(2,3,5)][
  between(r2,0.2,1)]

vec_sel3[between(x,145.0,145.335)][
  between(y,-37.45,-37.27)][K>L0] %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  coord_sf()+
  scale_fill_viridis_c(option='F')

vec_sel4 <- fits[isConv==TRUE][is.na(r2)==FALSE][L0>0][vc %in% c(2,3,5)][
  between(r2,0.2,1)][between(x,145.0,145.335)][
    between(y,-37.45,-37.27)][sample(.N,100)]

median(vec_sel4$date_fire1)

clim_sel4 <- clim[idx_awap %in% (vec_sel4$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_12mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel4$date_fire1))) %>% 
  as.data.table()

fits[id %in% vec_sel4$id] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  # left_join(., sdat, by='id') %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days)) ) %>% 
  # filter(post_days <= ttr5_lai+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  # mutate(r2 = format(r2,2)) %>% 
  ggplot(data=.,aes(post_days,pred,
                    # group=id,
                    # color=id,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel4[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=0,
                ymax=7,
                fill=val))+
  # scale_color_viridis_d(option='D',end=0.9,direction = -1)+
  # scale_color_viridis_c(option='D',end=0.9,direction = -1)+
  # scale_fill_viridis_c()+
  scale_fill_gradient2()+
  # geom_smooth()+
  geom_hline(aes(yintercept=0),col='grey90',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_line(lwd=0.1)+
  # geom_smooth(inherit.aes = F,
  #             aes(post_days,pred))+
  # geom_quantile(method = "rqss", lambda = 1000, 
  #               quantiles=c(0.05,0.5,0.95))+
  # geom_vline(aes(xintercept=ttr5_lai),col='black',lty=3)+
  scale_x_continuous(limits=c(0,3000), expand=c(0,0))+
  scale_y_continuous(limits=c(0,7), expand=c(0,0))+
  labs(x='Days after fire', 
       y='LAI ', 
       color='NVIS Vegetation Group', 
       fill="3 month Precip. Anomaly (mm)", 
       title='Kilmore East')+
  # facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
  #            ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')
ggsave(filename = 'figures/figure_killmoreEAst_logistic_lai12mo-TTR-Def5.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)


vec_sel4 <- fits[isConv==TRUE][is.na(r2)==FALSE][L0>0][vc %in% c(2,3,5)][
  between(r2,0.2,1)][between(x,145.0,145.335)][
    between(y,-37.45,-37.27)][sample(.N,100)]

median(vec_sel4$date_fire1)

clim_sel4 <- clim[idx_awap %in% (vec_sel4$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_12mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel4$date_fire1))) %>% 
  as.data.table()

p_kilmore <- mdat[id %in% vec_sel4$id][date_fire1==median(vec_sel4$date_fire1)] %>% 
  ggplot(data=.,aes(post_days,slai,
                    # group=id,
                    color=ttr5_lai,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel4[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_sel4$id]$kn_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_sel4$id]$kn_anom_3mo,na.rm=T),
                fill=val))+
  # scale_color_viridis_d(option='D',end=0.9,direction = -1)+
  # scale_color_viridis_c(option='D',end=0.9,direction = -1)+
  scale_color_viridis_c(option='B',
                        end=0.9,
                        direction = 1, 
                        limits=c(365,1500), 
                        oob=scales::squish)+
  scale_fill_gradient2(#limits=c(-1000,1000),
                       # breaks=c(-200,-100,0,100,200),
                       # labels=c("< -2000 ", "-1000", "0","1000"," > +2000"),
                       oob=scales::squish)+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2050), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_sel4$id]$slai,na.rm=TRUE),
             max(mdat[id %in% vec_sel4$id]$slai,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='LAI', 
       color='TTR (days)', 
       title='Kilmore East', 
       fill="P Anom. 12-mo (mm)")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'right'); p_kilmore

ggsave(filename = 'figures/figure_kilmoreEAst_timeSeries_lai_2000days.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)





# Plot lai anomaly date_fire1->recovery_date -----------------------------------
vec_sel <- fits[isConv==TRUE][is.na(r2)==FALSE][L0>0][L0<K][r<0.1][r2>0.25][vc %in% c(2,3,5)][
  ,.SD[sample(.N,1000)],by=vc_name]$id
mdat[id %in% vec_sel] %>% 
  as_tibble() %>% 
  mutate(val = slai_anom_12mo) %>% 
  filter(post_days <= ttr5_lai+365) %>% 
  # filter(post_days>=365) %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  mutate(months_since_july = if_else(fire_month<7, fire_month+6, fire_month-7)) %>%
  # pull(months_since_july) %>% table
  mutate(season = case_when(fire_month %in% c(9,10,11)~"SON",
                            fire_month %in% c(12,1,2) ~ "DJF", 
                            fire_month %in% c(3,4,5)~"MAM",
                            fire_month %in% c(6,7,8)~"JJA")) %>% 
  # filter(season %in% c("SON","DJF")) %>% 
  inner_join(., fits %>% select(id,vc_name) %>% as_tibble(),by='id') %>% 
  ggplot(data=.,aes(post_days,val,group=id,color=months_since_july))+
  # scale_color_viridis_d(option='D',end=0.8,direction = 1)+
  # scale_color_viridis_c(option='B')+
  scico::scale_color_scico(palette = 'batlow',end=0.99,direction = -1)+
  geom_line(lwd=0.2)+
  geom_hline(aes(yintercept=0),col='black')+
  # geom_vline(aes(xintercept=ttr5_lai),col='black',lty=3)+
  # geom_hline(aes(yintercept=malai),lty=3)+
  labs(x='Days after fire', 
       y=expression(paste(LAI["12 mo"]~ Anomaly)), 
       color='Months after July')+
  scale_x_continuous(limits=c(365,3500),
                     breaks=c(365,1000,2000,3000),
                     expand=c(0,0))+
  facet_wrap(~vc_name,ncol = 1)+
  # facet_grid(cut_number(delta,3)~cut_number(r,3),labeller = label_both)+
  # facet_grid(rows=vars(vc_name), 
  #            # cols = vars(id), 
  #            margins=F)+
  # facet_wrap(~paste0("x:",format(x,digits=4),
  #                    "  y:",format(y,digits=4)),scales = 'free_x',nrow = 3)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')
ggsave(filename = 'figures/figure_timeseries_lai-anom_by_season-vc.png',
       width=30,
       height=20,
       units='cm',
       dpi=350)



# NeXT figure --------------------------------------------------------------
# vec_sel <- sample(unique(fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0]$id),25)
vec_sel <- merge(fits,d_soil,by='id')[isConv==TRUE][is.na(r2)==FALSE][L0>0][L0<K][r<0.1][r2>0.25][vc %in% c(2,3,5)][
  ,.SD[sample(.N,10)],by=vc_name]$id


tmp1 <- merge(mdat[id%in%vec_sel][,.(post_days,slai_12mo,id)],fits,by='id')[,delta:=L0/K]

tmp1[,.(r,K,L0,id,ttr5_lai,date_fire1)] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,2000,length.out=50))) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days))) %>% 
  mutate(delta = L0/K) %>% 
  filter(post_days <= ttr5_lai+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  mutate(season = case_when(fire_month %in% c(9,10,11)~"SON",
                            fire_month %in% c(12,1,2) ~ "DJF", 
                            fire_month %in% c(3,4,5)~"MAM",
                            fire_month %in% c(6,7,8)~"JJA")) %>% 
  filter(season %in% c("DJF")) %>% 
  ggplot(data=.,aes(post_days,pred,group=id))+
  # scale_color_viridis_d(option='D',end=0.8,direction = 1)+
  # scale_color_viridis_c(option='B')+
  scico::scale_color_scico_d(palette = 'batlow',end=0.8)+
  geom_hline(aes(yintercept=0),col='grey')+
  geom_line(lwd=0.2,col='red')+
  geom_point(data=tmp1, 
             aes(post_days, slai_12mo, group=id), inherit.aes = F, 
             alpha=0.1)+
  # geom_vline(aes(xintercept=ttr5_lai),col='black',lty=3)+
  # geom_hline(aes(yintercept=malai),lty=3)+
  labs(x='Days after fire', 
       y='LAI', 
       color='')+
  facet_grid(cut_number(delta,3)~cut_number(r,3),labeller = label_both)+
  # facet_grid(rows=vars(vc_name), 
  #            # cols = vars(id), 
  #            margins=F)+
  # facet_wrap(~paste0("x:",format(x,digits=4),
  #                    "  y:",format(y,digits=4)),scales = 'free_x',nrow = 3)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')

fits[id%in%vec_sel][r<0.005] %>% ggplot(data=.,aes(abs(fire_doy-180),K-L0,color=r))+
  geom_point()+
  scale_color_viridis_c()


fits[r<0.005][L0>0][L0<K][r<0.1][r2>0.25][vc %in% c(2,3,5)] %>% as_tibble() %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  mutate(season = case_when(fire_month %in% c(9,10,11)~"SON",
                          fire_month %in% c(12,1,2) ~ "DJF", 
                          fire_month %in% c(3,4,5)~"MAM",
                          fire_month %in% c(6,7,8)~"JJA")) %>% 
  # pull(season) %>% table
  filter(season %in% c("SON","DJF")) %>% 
  select(season,L0,K,r,vc_name) %>% 
  pivot_longer(cols=c(L0,K,r)) %>% 
  ggplot(data=., aes(vc_name,value,group=paste(season,vc_name),fill=paste(season)))+
  geom_boxplot(outlier.colour = NA)+
  facet_wrap(~name,ncol = 1,scales='free')






fits %>% 
  as_tibble() %>% 
  mutate(val = slai_anom) %>% 
  # filter(post_days <= ttr5_lai+365) %>% 
  # filter(post_days>=365) %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  mutate(months_since_july = if_else(fire_month<7,fire_month+7,fire_month)) %>%
  mutate(season = case_when(fire_month %in% c(9,10,11)~"SON",
                            fire_month %in% c(12,1,2) ~ "DJF", 
                            fire_month %in% c(3,4,5)~"MAM",
                            fire_month %in% c(6,7,8)~"JJA")) %>% 
  # filter(season %in% c("SON","DJF")) %>% 
  inner_join(., fits %>% select(id,vc_name) %>% as_tibble(),by='id') %>% 
  ggplot(data=.,aes(y=ttr5_lai,
                    x=fire_month, 
                    group=fire_month))+
  geom_boxplot(outlier.colour = NA)+
  coord_cartesian(ylim=c(365,2500))


s_r <- stars::st_as_stars(fits[r<0.05][,.(x,y,r)])

ggplot()+
  geom_sf(data=filter(oz_poly,NAME_1%in%c('New South Wales','Victoria')),
          inherit.aes = F,
          fill=NA,
          color='black')+
  geom_stars(data=s_r)+
  scale_fill_viridis_c(limits=c(0, quantile(fits[r<0.05]$r,0.975,na.rm=TRUE)), 
                       na.value = "#00000000")+
  theme_linedraw()
  

fits[r<0.05] %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  # geom_sf(data=filter(oz_poly,NAME_1%in%c('New South Wales','Victoria')), 
  #         inherit.aes = F,
  #         fill=NA,
  #         color='black')+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0, quantile(fits[r<0.05]$r,0.975,na.rm=TRUE)))+
  coord_sf(crs="+proj=omerc +lon_0=148 + lat_0=-34 +k=1")
  # coord_sf(crs = "+proj=laea +lat_0=-36 +lon_0=148 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ", 
  # coord_sf(crs = "+proj=aeqd +lon_0=155", 
  #          # xlim=c(148-4,148.8+5),
  #          # ylim=c(-40, -29),
  #          expand = F)
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         # pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)

p$coordinates$limits
p$coordinates$clip  
p$coordinates  
ggplot()+
  geom_sf(data=filter(oz_poly,NAME_1%in%c('New South Wales','Victoria')))+
  coord_sf(crs="+proj=aeqd +lon_0=155")

filter(oz_poly,NAME_1%in%c('New South Wales','Victoria'))

mapproj::mapproject()


library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$name)), " countries)"))

ggplot(data = world) +
  geom_sf() +
  coord_sf(crs = "+proj=laea +lat_0=-30 +lon_0=145 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")+
  annotation_north_arrow(location = "bl",
                          which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)
  

library(sf)
ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97))
