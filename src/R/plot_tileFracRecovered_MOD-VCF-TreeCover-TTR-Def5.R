library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
# setDTthreads(threads = 16)

# Import data ------------------------------------------------------------------
dat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD44B_tree_grass_500m_SE_coastal_2001_2019.parquet")
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
id_train <- base[date < ymd("2019-08-01")] %>% 
  .[,.(nburns=.N),keyby=.(x,y,id)][nburns==1]

# dat1: Pre Black Summer fires ---------------------------------------
base_fire <- base %>% 
  group_by(id) %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  ungroup() %>% 
  mutate(fire_year1 = year) %>% 
  as.data.table()
base_fire <- base_fire[,.(id,fire_doy,fire_year1)]
dat1 <- dat[id%in%id_train$id][year < 2019]
dat1 <- merge(dat1,base_fire,on=c("id"),allow.cartesian = TRUE)
dat1[,`:=`(years_since_fire = year-fire_year1)]

d_min_tree_cover_anom <- dat1 %>% 
  lazy_dt() %>% 
  filter(is.na(tree_cover_anom)==FALSE) %>% 
  group_by(id) %>% 
  filter(between(year,fire_year1-1,fire_year1+1)) %>% 
  summarize(min_tree_cover_anom = min(tree_cover_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)
dat1 <- merge(dat1,d_min_tree_cover_anom,by='id')
gc(full=TRUE)



gc(full=TRUE,reset = TRUE)

fn_min <- function(x){ 
  suppressWarnings({
    out <- min(x,na.rm=TRUE)
    if(is.infinite(out)==TRUE){out <- NA_real_}
  })
  return(out)}

fn_ttr5 <- function(din){
  ttr <- din[years_since_fire>=1][tree_cover_anom >= -0.25*tree_cover_yr_sd]$years_since_fire
  ttr <- fn_min(ttr)
  din$ttr5_tree_cover <- ttr
  # din$pre_fire_delta_t_anom_3mo <- din[date == floor_date(date_fire1-1,'month')]$delta_t_anom_3mo[1]
  # din$pre_fire_delta_t_anom_12mo <- din[date == floor_date(date_fire1-1,'month')]$delta_t_anom_12mo[1]
  # din$pre_fire_delta_t_anom_24mo <- din[date == floor_date(date_fire1-1,'month')]$delta_t_anom_24mo[1]
  # din$pre_fire_delta_t_anom_36mo <- din[date == floor_date(date_fire1-1,'month')]$delta_t_anom_36mo[1]
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
arrow::write_parquet(out[year==fire_year1][,.(x,y,id,fire_year1,fire_doy,
                                              min_tree_cover_anom,
                                              matree_cover,
                                              tree_cover_anom,
                                              # pre_fire_tree_cover_anom_3mo,
                                              # pre_fire_tree_cover_anom_12mo,
                                              # pre_fire_tree_cover_anom_24mo,
                                              # pre_fire_tree_cover_anom_36mo, 
                                              ttr5_tree_cover)], 
                     sink = paste0("../data_general/proc_data_Oz_fire_recovery/fit_tree_cover_ttrDef5_preBS",Sys.time(),".parquet"))

# cleanup 
rm(out); gc(full=TRUE)


# FASTER Method -----------------------------------------------------------

out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_tree_cover_ttrDef5_preBS2021-04-19 11:39:00.parquet")
d3 <- expand_grid(out[year==fire_year1], 
                  post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = fire_year1 ) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum( (post_days/365) >= ttr5_tree_cover,na.rm=TRUE)/n()) %>% 
  ungroup()

d3 %>% write_parquet(.,
                     sink="../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_tree_cover_ttrDef5_preBS.parquet")

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
       fill='Tree Cover Fraction Recovered   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank())

ggsave(filename = 'figures/figure_cumulativeRecovered_byYear_TTR-Def5-tree_cover.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)


