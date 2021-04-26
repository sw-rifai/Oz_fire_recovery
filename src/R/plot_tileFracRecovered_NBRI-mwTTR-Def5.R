library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
# setDTthreads(threads = 16)

# Isolate slow recovering pixels --------------------------------
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
                            col_select = c("x","y","date","id",
                                           "ndvi","sndvi","ndvi_u","ndvi_sd", 
                                           "nbr","nbr_u","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id",
                                           "ndvi","sndvi","ndvi_u","ndvi_sd", "nbr","nbr_u","fire_doy"))
gc(full=TRUE)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)

dat <- dat[,.(x,y,date,id,ndvi,sndvi,ndvi_u,ndvi_sd, nbr,nbr_u,fire_doy)]
gc(full=TRUE)
norms <- dat %>% lazy_dt() %>% 
  mutate(year=year(date)) %>% 
  group_by(id,year) %>% 
  summarize(nbr_yr = mean(nbr,na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  summarize(manbr = mean(nbr_yr,na.rm=TRUE), 
            nbr_yr_sd = sd(nbr_yr,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

gc(full=TRUE)
dat <- merge(dat,norms,by='id')
rm(norms)
gc(full=TRUE)
dat[,`:=`(nbr_anom = nbr-nbr_u)]
dat <- dat[order(x,y,date)][,nbr_anom_3mo := frollmean(nbr_anom,
                                                        n = 3,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(x,y,date)][,nbr_anom_12mo := frollmean(nbr_anom,
                                                         n = 12,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(x,y,date)][,nbr_anom_24mo := frollmean(nbr_anom,
                                                         n = 24,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(x,y,date)][,nbr_anom_36mo := frollmean(nbr_anom,
                                                         n = 36,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(id,date)]
gc(full=TRUE)


id_train <- dat %>% 
  lazy_dt() %>%
  filter(date < ymd("2019-08-01")) %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y,id) %>%
  summarize(nburns = n()) %>%
  as.data.table() %>% 
  .[nburns==1]
gc(full=TRUE)

id_test <- dat %>% 
  lazy_dt() %>%
  filter(date >= ymd("2019-08-01")) %>% 
  filter(!(id %in% id_train$id)) %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y,id) %>%
  summarize(nburns = n()) %>%
  as.data.table() %>% 
  .[nburns==1]
gc(full=TRUE)


# dat1: Pre Black Summer fires ---------------------------------------
dat1 <- dat[id%in%id_train$id][date <= ymd("2019-08-01")]
gc(full=TRUE)
firedate_train <- dat1 %>% lazy_dt() %>%
  filter(fire_doy > 0) %>%
  group_by(id) %>%
  mutate(date_fire1 = date) %>%
  ungroup() %>%
  select(id,date_fire1) %>%
  as.data.table()
gc(full=TRUE)
dat1 <- left_join(dat1,firedate_train,by='id') %>% as.data.table()
gc(full=TRUE)
dat1 <- dat1 %>% lazy_dt() %>%
  mutate(days_since_fire = as.double(date - date_fire1)) %>%
  as.data.table()
gc(full=TRUE)
d_min_nbr <- dat1 %>% 
  lazy_dt() %>% 
  mutate(nbr_anom = nbr - nbr_u) %>% 
  filter(days_since_fire <= 100) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)
dat1 <- merge(dat1,d_min_nbr,by='id')
gc(full=TRUE)



gc(full=TRUE,reset = TRUE)

fn_min <- function(x){ 
  suppressWarnings({
    out <- min(x,na.rm=TRUE)
    if(is.infinite(out)==TRUE){out <- NA_real_}
  })
  return(out)}

fn_ttr5 <- function(din){
  ttr <- din[days_since_fire>=365][nbr_anom_12mo>= -0.25*nbr_yr_sd]$days_since_fire
  ttr <- fn_min(ttr)
  din$ttr5_nbr <- ttr
  din$pre_fire_nbr_anom_3mo <- din[date == floor_date(date_fire1-1,'month')]$nbr_anom_3mo[1]
  din$pre_fire_nbr_anom_12mo <- din[date == floor_date(date_fire1-1,'month')]$nbr_anom_12mo[1]
  din$pre_fire_nbr_anom_24mo <- din[date == floor_date(date_fire1-1,'month')]$nbr_anom_24mo[1]
  din$pre_fire_nbr_anom_36mo <- din[date == floor_date(date_fire1-1,'month')]$nbr_anom_36mo[1]
  din <- din[is.na(fire_doy)==FALSE]
  return(din)
}
grpn <- uniqueN(id_train$id)
system.time(
  out <- dat1[,
              {cat("progress",.GRP/grpn*100,"%\n"); 
                fn_ttr5(.SD)}, 
              by=.(id,x,y)]
)
arrow::write_parquet(out[date==date_fire1][,.(x,y,id,date_fire1,fire_doy,
                                              min_nbr_anom,
                                              manbr,
                                              nbr_anom_3mo,
                                              pre_fire_nbr_anom_3mo,
                                              pre_fire_nbr_anom_12mo,
                                              pre_fire_nbr_anom_24mo,
                                              pre_fire_nbr_anom_36mo, 
                                              ttr5_nbr)], 
                     sink = paste0("../data_general/proc_data_Oz_fire_recovery/fit_nbr_ttrDef5_preBS",Sys.time(),".parquet"))

# cleanup 
rm(out); gc(full=TRUE)




# FASTER Method -----------------------------------------------------------
out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_nbr_ttrDef5_preBS2021-04-19 11:39:00.parquet")
d3 <- expand_grid(out, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(post_days>ttr5_nbr,na.rm=TRUE)/n()) %>% 
  ungroup()


d3 %>% 
  filter(between(hydro_year,2001,2015)) %>% 
  mutate(valid = ifelse(post_days <= 365*(2021.3-hydro_year), T,F)) %>% 
  mutate(val = ifelse(valid==TRUE, val,NA)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
               filter(post_days==min(post_days,na.rm=TRUE)) %>% 
               ungroup(), 
             inherit.aes = FALSE, 
             aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish, 
                       na.value='grey50')+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Year of Bushfire',
       fill='Norm. Burn Ratio Fraction Recovered   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank())

ggsave(filename = 'figures/figure_cumulativeRecovered_byYear_TTR-Def5-NBR.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)
