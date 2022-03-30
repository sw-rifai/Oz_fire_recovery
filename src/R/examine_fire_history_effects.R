library(tidyverse)
library(data.table)
library(lubridate)
library(stars)
library(mgcv)
library(patchwork)
library(magic)
library(mgcViz)

# Rasterize NSW 20th century fire history -------------------------------------
fp <- list.files("../data_general/fire/NSW_fire_history/",'.gpkg',full.names = T)
fires <- st_read(fp[1])
(grd = st_as_stars(st_bbox(fires), nx = 1000, ny = 750, 
  # xlim = c(0, 1.0), ylim = c(0, 1), 
  values = 0))

# fire year
fy = st_rasterize(fires %>% 
    filter(FireYear < 200101) %>% 
    arrange(FireYear) %>% 
    select(FireYear) %>% 
    mutate(FireYear = as.numeric(substr(FireYear,1,4))), 
  grd, 
  options = c("MERGE_ALG=REPLACE", "ALL_TOUCHED=FALSE")) %>% 
  set_names('last_fire_year_lte2000')

# fire frequency
ff = st_rasterize(fires %>% 
    filter(FireYear < 200101) %>% 
    mutate(fire_bin = 1) %>%
    select(fire_bin)
  , 
  grd, 
  options = c("MERGE_ALG=ADD", "ALL_TOUCHED=FALSE")) %>% 
  set_names('fire_count_pre2000')
# fh %>% as.data.table() %>%
#   filter(FireYear>0) %>% 
#   mutate(year=substr(FireYear,1,4) %>% as.numeric) %>% 
#   pull(year) %>% summary()
#   ggplot(data=.,aes(x,y,fill=year))+
#   geom_tile()+
#   scale_fill_viridis_c(option='H')+
#   coord_sf()

fh_rec = st_rasterize(fires %>% 
    filter(FireYear < 200101) %>% 
    filter(FireYear > 199412) %>% 
    # select(FireYear) %>% 
    # mutate(fire_bin = sum(is.na(FireYear)==F)) %>% 
    # select(fire_bin)
    select(Branch) %>% 
    mutate(fire_bin = ifelse(Branch>=1,1,0)) %>% 
    select(fire_bin)
  , 
  grd, 
  options = c("MERGE_ALG=ADD", "ALL_TOUCHED=FALSE")) %>% 
  set_names('burn_90s')
# fh_rec %>% as_tibble() %>% pull(burn_90s) %>% table()

# Rasterize VID fire history ----------------------------------------------
fires_vic <- st_read("../data_general/fire/VIC_fire_history//VIC_fire_history.gpkg")
names(fires_vic) <- tolower(names(fires_vic))

(grd_vic = st_as_stars(st_bbox(fires_vic), nx = 1000, ny = 750, 
  # xlim = c(0, 1.0), ylim = c(0, 1), 
  values = 0))


tmp <- fires_vic %>% 
    rowwise() %>% 
    mutate(FireYear = 
        last(as.numeric(str_split(fireyears,pattern = ' ',simplify = T)))) %>% 
  ungroup()
tmp

# fire year
fy_vic <- st_rasterize(
    tmp %>% 
    filter(FireYear < 2001) %>% 
    arrange(FireYear) %>% 
    select(FireYear), 
  grd, 
  options = c("MERGE_ALG=REPLACE", "ALL_TOUCHED=FALSE")) %>% 
  set_names('last_fire_year_lte2000')

# ggplot()+geom_stars(data=fy_vic)+scale_fill_viridis_c(limits=c(1927,2001))

# fire frequency
fn <- function(x){
  vec <- as.numeric(str_split(x,pattern = ' ',simplify = T))
  out <- length(vec[vec<2001])
  return(out)}
tmp <- tmp %>% 
    filter(FireYear < 2001) %>% 
  rowwise() %>% 
    mutate(fire_count_pre2000 = fn(allyears)) %>% 
  ungroup()
ff_vic <- st_rasterize(tmp %>% select(fire_count_pre2000), 
  grd, 
  options = c("MERGE_ALG=ADD", "ALL_TOUCHED=FALSE")) %>% 
  set_names('fire_count_pre2000')
# fh %>% as.data.table() %>%
#   filter(FireYear>0) %>% 
#   mutate(year=substr(FireYear,1,4) %>% as.numeric) %>% 
#   pull(year) %>% summary()
#   ggplot(data=.,aes(x,y,fill=year))+
#   geom_tile()+
#   scale_fill_viridis_c(option='H')+
#   coord_sf()

fh_rec_vic <- st_rasterize(tmp %>% 
    filter(FireYear < 2001) %>% 
    filter(FireYear > 1994) %>% 
    # select(FireYear) %>% 
    # mutate(fire_bin = sum(is.na(FireYear)==F)) %>% 
    # select(fire_bin)
    select(firecount) %>% 
    mutate(fire_bin = ifelse(firecount>=1,1,0)) %>% 
    select(fire_bin)
  , 
  grd, 
  options = c("MERGE_ALG=ADD", "ALL_TOUCHED=FALSE")) %>% 
  set_names('burn_90s')
plot(fh_rec_vic,breaks='equal')
# ggplot()+geom_stars(data=fh_rec_vic)+scale_fill_viridis_c(limits=c(0,9))
# fh_rec_vic %>% as_tibble() %>% pull(burn_90s) %>% unique



# mosaic rasters ----------------------------------------------------------
## Fire frequency
ff %>% as_tibble()
ff <- ff %>% 
  mutate(fire_count_pre2000 = ifelse(fire_count_pre2000==0,NA,fire_count_pre2000))
ff_vic <- ff_vic %>% 
  mutate(fire_count_pre2000 = ifelse(fire_count_pre2000==0,NA,fire_count_pre2000))
# st_mosaic(ff, ff_vic) %>% plot(breaks='equal',col=viridis::viridis(9))
ff_all <- st_mosaic(ff, ff_vic)

## Fire year 
fy <- fy %>% mutate(last_fire_year_lte2000 = ifelse(last_fire_year_lte2000==0,NA,last_fire_year_lte2000))
fy_vic <- fy_vic %>% mutate(last_fire_year_lte2000 = ifelse(last_fire_year_lte2000==0,NA,last_fire_year_lte2000))
# st_mosaic(fy, fy_vic) %>% plot(breaks='equal',col=viridis::viridis(9))
fy_all <- st_mosaic(fy, fy_vic)


fh_rec <- fh_rec %>% mutate(burn_90s = ifelse(burn_90s==0,NA,1))
fh_rec_vic <- fh_rec_vic %>% mutate(burn_90s = ifelse(burn_90s==0,NA,1))
# st_mosaic(fh_rec,fh_rec_vic) %>% plot(breaks='even')
fh_rec_all <- st_mosaic(fh_rec,fh_rec_vic)

# load TTR fits ---------------------------------------------------------------
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-22 09:37:07.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>=0][r<0.024][r2>0.333][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
 ]
# estimate TTR from the logistic function
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]

xy_fits <- fits[,.(x,y)] %>% st_as_sf(coords=c('x','y')) %>% st_set_crs(4326)
fr <- c(ff_all,fy_all,fh_rec_all)
fr_4326 <- st_warp(fr, crs = st_crs(4326))
vv <- st_extract(fr_4326,xy_fits)
fits <- bind_cols(vv,fits)

fits <- fits %>% as.data.table() %>% 
  .[,fire_month := factor(month)]

# fh_rec_4326 <- st_warp(fh_rec, crs = st_crs(4326))
# vv <- st_extract(fh_rec_4326,xy_fits)
# fits <- bind_cols(vv %>% st_drop_geometry(),fits)


# Using NSW -------------------------------------------------------------------
# pre_fire_freq NA values are from outside NSW
# More frequently burned forests have   
tmp <- fits %>% 
  # filter(is.na(pre_fire_freq)==F) %>% 
  mutate(burn_90s = ifelse(is.na(burn_90s)==T,0,burn_90s)) %>% 
  mutate(burn_90s = factor(burn_90s,ordered = F)) %>% 
  mutate(fmonth = factor(month)) %>% 
  mutate(last_fire_year_lte2000 = ifelse(is.na(last_fire_year_lte2000)==T,1934,last_fire_year_lte2000)) %>% 
  mutate(tsd_pre2001 = 2001-last_fire_year_lte2000) %>%
  mutate(tsd_pre2001 = ifelse(tsd_pre2001==2001,100,tsd_pre2001))
  # filter(tsd_pre2001 < 50)
  # pull(tsd_pre2001)
  # filter(pre_fire_freq < 5) %>% 
  # filter(tsd_pre2001 < 100) 
tmp[,`:=`(bs = 1 - (min_slai_anom+slai_u)/slai_u)][,`:=`(bs = ifelse(bs>1,1,bs))]
tmp[,`:=`(frac_pf_slai_anom_3mo = pre_fire_slai_anom_3mo/malai, 
  frac_pf_slai_anom_12mo = pre_fire_slai_anom_12mo/malai)]
# tmp[,`:=`(frac_p12_anom = precip_anom_12mo/map,
#          frac_vpd3_anom = vpd15_anom_3mo/mavpd15)]

# Time since disturbance effect --------------------------------------
b0 <- 
     bam(
        bs ~
       te(x,y)+
      s(tsd_pre2001,k=5),
      data=tmp, 
    family=Gamma(link='identity'),
    select=T,
    discrete = T)
summary(b0)
b0 %>% plot(pages=1,scale=0,rug=T,scheme=2)

b1 <- 
     bam(
      I(r*365) ~
     te(x,y)+
      s(tsd_pre2001,bs='cs',k=5), 
      data=tmp, 
    family=Gamma(link='identity'),
    select=T,
    discrete = T) 
b1 %>% summary
b1 %>% plot(pages=1,scale=0,rug=T,scheme=2)
library(gratia); library(patchwork)
p_1 <- draw(b1,select = 'tsd_pre2001',partial_match = T,overall_uncertainty = F)+
  geom_hline(aes(yintercept=0),col='grey')+
  labs(x="Time Since Last Fire (years)", 
    title=NULL)+
  scale_x_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid = element_blank())
p_2 <- tmp %>%
  ggplot(data=.,aes(tsd_pre2001))+
  geom_histogram(bins=30,fill='red')+
  # geom_vline(aes(xintercept=66),lty=2)+
  geom_rect(aes(xmin=64,xmax=71,ymin=0,ymax=50000), 
    fill='grey',
    alpha=0.01,
    color=NA)+
  scale_x_continuous(limits=c(-1,71),
    expand=c(0,0), 
    breaks = c(0,20,40,60,67), 
    labels=c("0","20","40","60","≥67"))+
  scale_y_continuous(expand = c(0,0))+
    labs(x="Time Since Last Fire (years)", 
    title=NULL)+
  theme_linedraw()+
  theme(panel.grid = element_blank())
p_out <- (p_1+p_2)+plot_layout(ncol=2)+
  plot_annotation(tag_levels = 'a',tag_prefix = '(',tag_suffix = ')')
ggsave(p_out, 
  filename = 'figures/fig_fire_history_effects.png',
  width=16,
  height=8,
  units='cm',
  dpi=350)


# b2 <- bam(I(r*365)~
#     s(tsd_pre2001,k=5) +
#     s(bs,k=5)+
#     s(frac_pf_slai_anom_3mo,k=5)+
#     # s(mapet,k=5)+
#     fire_month+
#     # s(post_precip_anom_frac,k=5)+ # seems real
#     # s(log(mappet),k=5),
#     te(x,y),
#   family=Gamma(link='identity'),
#   # method='fREML',
#   data=tmp, 
#   select=T,
#   discrete=TRUE
#   )
# summary(b2)
# # dev.new()
# plot(b2,pages=1,scale=0,rug=T,scheme=2)


# b2_l <- bam(I(r*365)~
#     tsd_pre2001 + 
#     # s(tsd_pre2001,k=5) + 
#     s(bs,k=5)+
#     s(frac_pf_slai_anom_3mo,k=5)+
#     # s(mapet,k=5)+
#     fire_month+
#     # s(post_precip_anom_frac,k=5)+ # seems real
#     # s(log(mappet),k=5),
#     te(x,y),
#   # family=Gamma(link='identity'),
#   # method='fREML',
#   data=tmp, 
#   select=T,
#   discrete=TRUE
#   )
# summary(b2_l)
# plot(b2_l,all.terms = T,pages = 1)


# Effect of fire frequency in 1900s -----
f1 <- 
   bam(
      I(r*365) ~
       te(x,y)+
       # fire_count_pre2000,
      s(fire_count_pre2000,k=5,bs='cs'),
      data=tmp, 
    family=Gamma(link='identity'),
    select=T,
    discrete = T)
summary(f1)
plot(f1,pages=1,scale=0,rug=T,scheme=2)

f2 <- 
   bam(bs~
       te(x,y)+
       # fire_count_pre2000,
      s(fire_count_pre2000,k=5,bs='cs'),
      data=tmp, 
    family=Gamma(link='identity'),
    select=T,
    discrete = T)
summary(f2)
plot(f2,pages=1,scale=0,rug=T,scheme=2)



# Effect of fire in the 1990s --------------
b4 <- 
   bam(
     # bs ~ 
      I(r*365) ~
     te(x,y)+
       burn_90s
      # s(burn_90s,k=5,bs='cs')
     , 
      data=tmp, 
    family=Gamma(link='identity'),
    select=T,
    discrete = T)
summary(b4)
plot(b4,pages=1,scale=0,rug=T,scheme=2,all.terms = T)

b5 <- 
   bam(
     bs ~
     te(x,y)+
       burn_90s
      # s(burn_90s,k=5,bs='cs')
     , 
      data=tmp, 
    family=Gamma(link='identity'),
    select=T,
    discrete = T)
summary(b5)
plot(b5,pages=1,scale=0,rug=T,scheme=2,all.terms = T)



tmp %>% 
  # filter(tsd_pre2001 < 100) %>%
  # filter(last_fire_year_lte2000 > 1900) %>% 
  filter(tsd_pre2001 <= 50) %>% 
  ggplot(data=.,aes(tsd_pre2001, slai_anom_3mo))+
  # geom_boxplot(outlier.colour = NA)
  geom_point(alpha=0.25)+
  geom_smooth(method='lm')
  # geom_smooth(method='bam',
  #   formula=y~s(x,bs='cs'),
  #   method.args=list(select=T,
  #     discrete=T))
junk <- tmp %>% 
  filter(tsd_pre2001 <= 50) %>% 
  lm(ttr5_lai ~ scale(log(ldk))+
      scale(slai_anom_3mo)+
      scale(lai_ma)+
      scale(tsd_pre2001),data=.)  
junk %>% summary

fires %>%
  mutate(pb = str_detect(Label,"Prescribed")) %>% 
  select(pb, Shape_Area) %>% 
  ggplot(data=.,aes(Shape_Area,fill=pb))+
  geom_histogram(bins=100)+
  scale_y_log10()+
  facet_wrap(~pb)


fits %>% 
  filter(is.na(pre_fire_freq)==F) %>% 
  filter(pre_fire_freq < 5) %>% 
  # mutate(pre_fire_freq = if_else(is.na(pre_fire_freq),0,pre_fire_freq)) %>% 
  ggplot(data=.,aes(pre_fire_freq, r))+
  # geom_point(alpha=0.01, position = 'jitter', 
  #   color='grey25')+
  # geom_violin(aes(x=pre_fire_freq,y=lai_ma,group=pre_fire_freq))+
  geom_smooth(method='lm')+
  # geom_smooth(method='bam', 
  #   formula=y~s(x,bs='cs'),
  #   method.args=list(
  #     select=T, 
  #     discrete=T),
  #   color='#cf0000')+
  labs(x="N. Fires 1900-2000") + 
  coord_cartesian(expand = F) +
  facet_wrap(~cut(lai_ma,c(0:5,Inf)))+
  theme_linedraw()+
  theme(panel.grid = element_blank())


r %>% 
  # select(fire_bin) %>% 
  plot(nbreaks=11,
  breaks='equal',
  col=viridis::turbo(10)
  )
ggplot()+
  geom_stars(data=r)+
  scale_fill_viridis_c()
r %>% as.data.table() %>% 
  .[fire_bin>0] %>% 
  .[order(fire_bin)]



fp
st_read(fp)
fires

fires %>% 
  # head() %>% 
  select(StartDate) %>% 
  plot()

fires %>% 
  mutate(year=substr(FireYear,1,4) %>% as.numeric) %>% 
  select(year,AreaHa) %>% 
  sf::st_drop_geometry() %>% 
  group_by(year) %>% 
  summarize(val = sum(AreaHa,na.rm=T)) %>% 
  ggplot(data=.,aes(year,val))+
  geom_point()


fits[,`:=`(fire_month = factor(month,ordered = F))]
fits <- fits %>% mutate(fmonth=month(date_fire1,label = T))
dat_train <- fits[sample(.N, floor(nrow(fits)*0.5))]
dat_test <- fsetdiff(
  dat_train %>% select(-geometry),
  fits %>% select(-geometry), all=T) 
b2 <- bam(
      I(r*1000) ~
          fmonth +
      # s(ldk, k=5)+
      # s(nbr_anom,k=5)+
      # s(log(ldk))+
      s( I(1-(min_slai/lai_ma)))+
      # (fire_count_pre2000)+
      # factor(fire_count_pre2000, ordered = T) +
      # factor(burn_90s,ordered = T) +
      # s(log(tsd_pre2001), bs='ad')+
      # s(log10(tsd_pre2001),k=4)+
      # s(lai_ma, k=5) +
      # s(ldk)+
      # scale(min_nbr_anom)+
      s(slai_anom_3mo,bs='cs',k=5) + 
      # s(slai_anom_3mo,bs='cs')+
      s(lai_ma,bs='cs') + 
      # te(lai_ma, slai_anom_3mo, k=5, bs='cs')+
      te(x,y,fx = F, k=c(30,30)),
  # te(x,y, k=c(10,10)),
  # data=tmp,
    # data=fits[lai_ma > 3],
    data=dat_train[last_fire_year_lte2000 < 1995],
    select=T,
    discrete = T) 
plot(b2,scheme=2, pages=1,
  hcolors=scico::scico(10,palette='roma'),all=T,rug=F)
summary(b2)

plot(sm(getViz(b2), select=1 ))+l_fitLine()+l_rug()
  labs(y='Predicted TTR (days)', 
    x="Pre-fire 3-mo LAI Anomaly")
ggsave(filename = 'figures/lai-anom-3mo_ttr_effect.png', 
  width=12,
  height=10,
  units='cm',
  dpi=350)

yardstick::rsq_trad_vec(truth=exp(dat_test$r), 
  estimate=(predict(b2,newdata=dat_test,type='response')))




# tmp %>% 
#   mutate(val = min_slai/lai_ma) %>% 
#   select(val,ldk) %>% cor
