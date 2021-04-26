library(data.table); 
library(tidyverse);
library(stars);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(dtplyr)


# Import data ------------------------------------------------------------------
dat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD44B_tree_grass_500m_SE_coastal_2001_2019.parquet")
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_tree_cover_ttrDef5_preBS2021-04-20 16:21:42.parquet")
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
                            col_select = c("x","y","id","date","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","id","date","fire_doy"))
gc(full=TRUE)
base <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)



# Calc anoms --------------------------------------------------------------
norms <- dat %>% lazy_dt() %>% 
  group_by(id,year) %>% 
  summarize(tree_cover_yr = mean(tree_cover, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  summarize(matree_cover = mean(tree_cover_yr,na.rm=TRUE), 
            tree_cover_yr_sd = sd(tree_cover_yr,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

dat <- merge(dat,norms,by=c("id"))
dat[,`:=`(tree_cover_anom = tree_cover - matree_cover)]

# Get fire doy -----------------------------------------------------------------
base <- base[is.na(fire_doy)==FALSE][,`:=`(year=year(date))]
id_train <- base[date < ymd("2019-08-01")
                 ][,.(nburns=.N),keyby=.(x,y,id)
                   ][nburns==1]


# dat1: Pre Black Summer fires ---------------------------------------
base_fire <- base %>% 
  group_by(id) %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  ungroup() %>% 
  mutate(fire_year1 = year) %>% 
  as.data.table()
base_fire <- base_fire[,.(id,fire_doy,fire_year1)]
dat1 <- dat[id%in%id_train$id]
dat1 <- merge(dat1,base_fire,on=c("id"),allow.cartesian = TRUE)
dat1[,`:=`(years_since_fire = year-fire_year1)]

d_min_tree_cover_anom <- dat1[is.na(tree_cover_anom)==FALSE][fire_year1 < 2019
                              ][year >= (fire_year1-1) & year<=(fire_year1+1)
                  ][,.(min_tree_cover_anom = min(tree_cover_anom,na.rm=TRUE)), 
                  by=.(id)]

gc(full=TRUE)
dat1 <- merge(dat1,d_min_tree_cover_anom,by='id')
gc(full=TRUE)

d_year1 <- dat1 %>% group_by(id) %>% 
  filter(tree_cover_anom==min_tree_cover_anom) %>% 
  mutate(year0 = year) %>% 
  select(id,year0) %>% 
  as.data.table()

dat1 <- merge(dat1,d_year1,by='id',allow.cartesian = TRUE) %>% 
  mutate(years_since_fire = year-year0) %>% 
  as.data.table()

dat1 <- merge(dat1, sdat[,.(id,ttr5_tree_cover)], by='id',allow.cartesian = TRUE)
dat1 <- dat1[is.na(ttr5_tree_cover)==FALSE]
mdat <- dat1[years_since_fire>=0]
mdat <- mdat[years_since_fire <= (ttr5_tree_cover+3)]


# Fit Weibull -------------------------------------------------------------
fn_w4 <- function(din){
  # set.seed(333)
  try(fit <- nls_multstart(tree_cover ~ Asym - Drop*exp(-exp(lrc)*years_since_fire^(pwr)), 
                           data=din,
                           # iter=1,
                           iter=1000,
                           supp_errors = 'Y',
                           start_lower = c(Asym=0,Drop=0, lrc=-10, pwr=-10),
                           start_upper = c(Asym=0.1, Drop=0.5, lrc=-5, pwr=10), 
                           lower= c(Asym=0,Drop=-0.7*100, lrc=-1000, pwr=-100), 
                           upper = c(Asym=100, Drop=100, lrc=2000, pwr=200))
      ,silent = TRUE)
  if(exists('fit')==FALSE){
    out <- data.table(Asym=NA_real_,Drop=NA_real_,lrc=NA_real_,pwr=NA_real_,isConv=FALSE)
  }
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(Asym=NA_real_,Drop=NA_real_,lrc=NA_real_,pwr=NA_real_,isConv=FALSE)
  }
  ,silent=TRUE)
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <- fit %>% coef(.) %>% t() %>% as.data.table()
    out$isConv <- fit$convInfo$isConv
  },silent=TRUE)
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out$r2 <- yardstick::rsq_trad_vec(truth = din$tree_cover, 
                                      estimate = predict(fit))
    out$rmse <- yardstick::rmse_vec(truth = din$tree_cover,
                                      estimate = predict(fit))
  },silent=F)
  out$nobs_til_recovery <- nrow(din)
  return(out)
}
vec_ids <- mdat[min_tree_cover_anom <= -10]$id %>% unique

fn_w4(mdat[id==673])
fn_w4(mdat[id==163400])
fn_w4(mdat[id==947745])

mdat[id==947745] %>% 
  ggplot(data=.,aes(years_since_fire,tree_cover))+geom_point()+geom_line()+
  geom_hline(aes(yintercept=matree_cover),col='red')

mdat$min_tree_cover_anom %>% hist

# furrr approach -----------------------------------------------
library(furrr)
plan(multisession, workers=8) # how many workers causes the parellelization to fail?
system.time(out <- mdat[id %in% sample(vec_ids,1000)] %>% 
              split(.$id) %>%
              future_map(~fn_w4(.x), .progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
setDT(out)
out[,`:=`(id=as.integer(id))]
arrow::write_parquet(out,
                     sink=paste0("outputs/weibull4Param_MODIS-VCF-TreeCover_fits_1burn_2001-2014fires_",Sys.time(),".parquet"))
# END ****************************************************************

# mdat[id %in% sample(unique(mdat$id),10)] %>% 
#   split(.$id)

out[is.na(r2)==F][isConv==T]$rmse %>% hist


fn_w3(dat1[id==70][years_since_fire>=0])

ggplot(data=dat1[id==70],aes(years_since_fire, tree_cover_anom,group=id))+
  geom_line()

dat1[id %in% sample(unique(dat1$id), 10)][] %>% 
 ggplot(data=.,aes(years_since_fire, tree_cover_anom,group=id))+
  geom_line()

dat1[id %in% sample(unique(dat1$id), 5)][years_since_fire>=0] %>% 
  ggplot(data=.,aes(years_since_fire, tree_cover_anom,group=id))+
  geom_line()+
  geom_vline(aes(xintercept=ttr5_tree_cover),col='blue')+
  facet_wrap(~id)


out[is.na(r2)==F][r2 > 0.2][isConv==T]$Asym %>% summary
out[is.na(r2)==F][r2 > 0.2][isConv==T]$Drop %>% summary

out[is.na(r2)==F][isConv==T][Drop < 0] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
  .[sample(.N,5)] %>% 
  expand_grid(.,
              post_years=seq(0,10,length.out = 100)) %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_years^(pwr))) %>% 
  # left_join(., mdat, by=c('id','post_days')) %>% 
  as_tibble() %>% 
  inner_join(., 
             distinct(sdat[,.(id,matree_cover)]) %>% as_tibble(), 
             by='id') %>% 
  # mutate(recovered = ifelse(pred > 0, 1,0)) %>%
  # filter(recovered == 0) %>% 
  # mutate(form = ifelse(lrc <= -10, 'slow','fast')) %>% 
  # group_by(form) %>% 
  # summarize(val = mean(ttr5))
 ggplot(data=.,aes(post_years,pred,group=id,color=factor(id)))+
  # ggpointdensity::geom_pointdensity()+
  geom_line(lwd=0.5)+
  geom_hline(aes(yintercept=matree_cover,group=id,color=factor(id)),lty=3)+
  # geom_vline(inherit.aes=F, 
  #            aes(xintercept=mean(ttr5),
  #                group=form))+
  # geom_abline()+
  # geom_smooth(inherit.aes = F, 
  #             aes(post_days,pred))+
  scico::scale_color_scico_d(end=0.7,
                           begin=0,
                           palette = 'batlow')+
  # scico::scale_color_scico(end=0.7,
  #                            begin=0,
  #                            palette = 'batlow')+
  scale_x_continuous(limits=c(0,10), 
                     expand=c(0,0))+
  labs(x='post_days', 
       y='Tree Cover anom',
       title=' Bushfires')+
  # facet_wrap(~form,ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank())

vec_ids <- unique(mdat$id)
mdat[id %in% vec_ids[1:3]] %>% 
  ggplot(data=.,aes(years_since_fire,tree_cover_anom,group=id,color=factor(id)))+
  geom_line()

fn_w4(mdat[id==673])
