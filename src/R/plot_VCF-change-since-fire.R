library(tidyverse)
library(dtplyr)
library(data.table)
library(lubridate)
library(stars); library(sf)
library(mgcv)

hydro_year_offset <- 3 # (months)

vcf <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD44B_tree_grass_500m_SE_coastal_2001_2019.parquet")
dttr <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet")
grid <- read_stars("../data_general/proc_data_Oz_fire_recovery/grid_500m_SE_coastal.tif")
d_soil <- arrow::read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")

"
| vc|vc_name                                                |
|--:|:------------------------------------------------------|
|  1|Rainforests and Vine Thickets                          |
|  2|Eucalypt Tall Open Forests                             |
|  3|Eucalypt Open Forests                                  |
|  4|Eucalypt Low Open Forests                              |
|  5|Eucalypt Woodlands                                     |
|  6|Acacia Forests and Woodlands                           |
|  7|Callitris Forests and Woodlands                        |
|  8|Casuarina Forests and Woodlands                        |
|  9|Melaleuca Forests and Woodlands                        |
| 10|Other Forests and Woodlands                            |
| 11|Eucalypt Open Woodlands                                |
| 14|Mallee Woodlands and Shrublands                        |
| 15|Low Closed Forests and Tall Closed Shrublands          |
| 16|Acacia Shrublands                                      |
| 17|Other Shrublands                                       |
| 18|Heathlands                                             |
| 19|Tussock Grasslands                                     |
| 20|Hummock Grasslands                                     |
| 21|Other Grasslands, Herblands, Sedgelands and Rushlands  |
| 22|Chenopod Shrublands, Samphire Shrublands and Forblands |
| 23|Mangroves                                              |
| 24|Inland Aquatic - freshwater, salt lakes, lagoons       |
| 25|Cleared, Non-Native Vegetation, Buildings              |
| 26|Unclassified Native Vegetation                         |
| 27|Naturally Bare - sand, rock, claypan, mudflat          |
| 28|Sea and Estuaries                                      |
| 29|Regrowth, Modified Native Vegetation                   |
| 31|Other Open Woodlands                                   |
| 99|Unknown/No Data                                        |
"

# Unnecessary???
# # Prep MCD64 --------------------------------------------------------------------
# ba <- read_stars("../data_general/proc_data_Oz_fire_recovery/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
#                  proxy=TRUE) %>% 
#   set_names("fire_doy") %>% 
#   st_warp(., grid, use_gdal = F) %>% 
#   st_set_dimensions(., 3, 
#                     values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
#                     names='date')
# 
# ba <- ba %>% as.data.frame()
# setDT(ba)
# gc(full=TRUE)
# ba <- ba[fire_doy > 0]
# ba[,`:=`(year=year(date), 
#          month=month(date))]
# ba[,`:=`(hydro_year = year(date-months(3)))]
# # Done *************************************************************************

dttr <- merge(dttr,d_soil, by=c('id'))
dttr[,`:=`(hydro_year_fire1 = year(date_fire1-months(hydro_year_offset)))]
vcf <- merge(vcf, dttr, all.x = TRUE, all.y = FALSE, 
      by.x = c("x","y","id"), 
      by.y = c("x", "y","id"), 
      allow.cartesian = TRUE)
vcf <- vcf[is.na(hydro_year_fire1)==FALSE]
vcf[vc%in%1:14]


vec_ids <- vcf[vc%in%c(1,2,3,4,5)][between(min_nbr_anom,-1,-0.5)][hydro_year_fire1==2002][sample(.N,10000)]$id
         
dttr[id%in%vec_ids][,.(val=decimal_date(date_fire1+days(ttr5))),by=.(id,vc_name)] %>% 
  .[.(mval = mean(val)),by=.(vc_name)]
lazy_dt(dttr) %>% filter(id%in%vec_ids) %>% 
  mutate(val = decimal_date(date_fire1+days(ttr5))) %>% 
  group_by(vc_name) %>% 
  summarize(ttr5 = mean(val)) %>% 
  ungroup() %>% 
  show_query()
dttr[id %in% vec_ids][, `:=`(val = decimal_date(date_fire1 + 
    days(ttr5)))][, .(ttr5 = mean(val)), keyby = .(vc_name)]


vcf[year<2019][id %in% vec_ids] %>% 
  as_tibble() %>% 
  select(id,year,vc_name, tree_cover,grass_cover) %>% 
  mutate(period=ifelse(year<2002,'pre','post')) %>% 
  pivot_longer(cols=c(tree_cover,grass_cover), 
               names_to = "veg", 
               values_to='cover') %>% 
  ggplot(data=.,aes(year, cover,group=paste(id,veg),color=veg))+
  # geom_line()+
  geom_smooth(inherit.aes = FALSE, 
              data=. %>% filter(year%in%2000:2002),
              aes(year,cover,color=veg), 
              method='lm')+
  geom_smooth(inherit.aes = FALSE, 
              data=. %>% filter(period=='post'),
              aes(year,cover,color=veg), 
              method='gam', 
              formula=y~s(x))+
  geom_vline(aes(xintercept=2002))+
  geom_vline(data=dttr[id %in% vec_ids][, `:=`(val = decimal_date(date_fire1 + 
    days(ttr5)))][, .(ttr5 = quantile(val, 0.25,na.rm=TRUE)), keyby = .(vc_name)], 
    aes(xintercept=ttr5),col='blue',lty=3)+
  geom_vline(data=dttr[id %in% vec_ids][, `:=`(val = decimal_date(date_fire1 + 
            days(ttr5)))][, .(ttr5 = median(val,na.rm=TRUE)), keyby = .(vc_name)], 
             aes(xintercept=ttr5),col='blue')+
  geom_vline(data=dttr[id %in% vec_ids][, `:=`(val = decimal_date(date_fire1 + 
            days(ttr5)))][, .(ttr5 = quantile(val,0.75,na.rm=TRUE)), keyby = .(vc_name)], 
             aes(xintercept=ttr5),col='blue',lty=3)+
  facet_wrap(~vc_name,ncol = 1,scales = 'free')

vcf[vc%in%2:7][near(year,hydro_year_fire1, tol=1.1)][id==1115] # window of fire
vcf[vc%in%2:7][near(year,hydro_year_fire1+(ttr5/365), tol=1.1)][id==1115] # window of recovery

d1 <- vcf[vc%in%2:7][near(year,hydro_year_fire1, tol=1.1)][,
              .(tree_min_fire = min(tree_cover),
                tree_max_fire = max(tree_cover),
                grass_max_fire = max(grass_cover),
                grass_min_fire = min(grass_cover),
                tree_loss = max(tree_cover)-min(tree_cover), 
                grass_gain = max(grass_cover)-min(grass_cover)),
              by=.(x,y,id,hydro_year_fire1)]

d2 <- vcf[vc%in%2:7][near(year,hydro_year_fire1+(ttr5/365),tol=1.1)] %>% 
                           .[,.(tree_max_ttr = max(tree_cover), 
                           grass_max_ttr = max(grass_cover),
                           tree_min_ttr = min(tree_cover), 
                           grass_min_ttr = min(grass_cover),
                           min_nbr_anom = unique(min_nbr_anom)), 
                        by=.(x,y,id,ttr5,vc,vc_name)]
d3 <- vcf[vc%in%2:7][hydro_year_fire1 %in% 2001:2015] %>% 
         .[,.(hydro_year = median(hydro_year_fire1)),by=.(id)] %>% 
         .[,.(nobs = .N), by=.(hydro_year)]

merge(d1,d2,by=c('x','y','id')) %>% 
  .[between(hydro_year_fire1,2001,2015)] %>% 
  .[vc %in% c(2,3,5)] %>% 
  # .[sample(.N,100)] %>% 
  as_tibble() %>% 
  group_by(id) %>% 
  summarize(`Tree Cover` = tree_max_ttr-tree_max_fire,
            `Non-tree Veg. Cover` = grass_min_ttr-grass_min_fire,
            hydro_year = unique(hydro_year_fire1)) %>% 
  ungroup() %>% 
  pivot_longer(cols=c(`Tree Cover`, `Non-tree Veg. Cover`), 
               names_to="veg",
               values_to="cover") %>% 
  inner_join(., as_tibble(d3), by='hydro_year') %>% 
  ggplot(data=.,aes(x=cover,
                    y=factor(hydro_year),
                    fill=nobs*0.25 #stat(x)
                    ))+
  geom_vline(aes(xintercept=0),col='grey60',lwd=1)+
  ggridges::geom_density_ridges_gradient(
    rel_min_height=0.01,
    scale=1,
    quantile_lines = TRUE, 
    quantiles=c(0.5),
    # quantiles = c(0.05,0.25,0.5,0.75,0.95), 
    alpha = 0.7)+
  labs(x='(% at Recovery)-(% at pre-Fire)', 
       y='Hydraulic Year', 
       fill=expression(paste(km**2)))+
  scale_x_continuous(expand=c(0,0))+
  scico::scale_fill_scico(end=0.9)+
  coord_cartesian(xlim=c(-25,25))+
  facet_wrap(~veg,ncol=2)+
  theme_linedraw()+
  theme(text = element_text(family = 'Palatino'),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(linetype = 3, color = 'black', size = 0.1),
        # panel.border = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black', 
                                  face='bold'))
ggsave(filename="figures/figure_change-vcf-pre-post-fire.png", 
       width=20,
       height=15, 
       units='cm', 
       dpi=350)


merge(d1,d2,by=c('x','y','id')) %>% 
  .[hydro_year_fire1 %in% 2003] %>% 
  .[vc %in% c(2,3,5)] %>% 
  .[sample(.N,100)] %>% 
  ggplot(data=.)+
  geom_segment(aes(x=hydro_year_fire1,
                   y=tree_min_at_fire, 
                   xend=hydro_year_fire1+(ttr5/365), 
                   yend=tree_at_ttr, 
                   color=min_nbr_anom), 
               lwd=0.2)+
  labs(y='Tree Cover (%)')+
  scale_x_continuous(expand=c(0,0))+
  scico::scale_color_scico(direction=-1)+
  facet_wrap(~vc_name,ncol=1,scales = 'free_y')+
  theme_linedraw()

d1$hydro_year_fire1 %>% table %>% sort

merge(d1,d2,by=c('x','y','id')) %>% 
  .[hydro_year_fire1 <= 2015] %>% 
  .[hydro_year_fire1 >= 2001] %>% 
  # .[hydro_year_fire1 %in% 2003] %>% 
  .[vc %in% c(2,3,5)] %>% 
  # .[sample(.N,10000)] %>% 
  ggplot(data=., aes(min_nbr_anom, tree_loss,color=vc_name))+
  geom_smooth()+
  facet_wrap(~hydro_year_fire1)


merge(d1,d2,by=c('x','y','id')) %>% 
  .[hydro_year_fire1 <= 2015] %>% 
  .[hydro_year_fire1 >= 2001] %>% 
  # .[hydro_year_fire1 %in% 2003] %>% 
  .[vc %in% c(2,3,5)] %>% 
  bam(tree_loss~min_nbr_anom+s(hydro_year_fire1,bs='re'),
      data=.) %>% 
  summary

# approximately: tree_loss ~ 0.71 + -47*min_nbr_anom


# SCRAP -------------------------------------------------------------------


lazy_dt(vcf) %>% group_by(id) %>% filter(year>=hydro_year_fire1) %>% ungroup() %>% show_query()
tmp <- lazy_dt(vcf) %>% 
  filter(is.na(hydro_year_fire1)==FALSE) %>% 
  group_by(id) %>% 
  filter(year>=hydro_year_fire1) %>% 
  ungroup() %>% 
  mutate(post_years = )
  as.data.table()


d_soil %>% select(vc,vc_name) %>% distinct() %>% as_tibble() %>% 
  arrange(vc) %>% 
  knitr::kable()



fn <- function(nobs, alpha=-1, beta=1){
  x <- rnorm(nobs)
  y <- rnorm(nobs, alpha+x*beta, sd=1)
  out <- tibble(x=x,
                y=y, 
                group = factor(sample.int(1e6,1)))
  return(out)
}
bind_rows(fn(100), 
          fn(100), 
          fn(100), 
          fn(100), 
          fn(25,alpha=5,beta=-5)) %>% #ggplot(data=.,aes(x,y,color=group))+geom_point()
  lm(y~x, data=.) %>% 
  summary

bind_rows(fn(100), 
          fn(100), 
          fn(100), 
          fn(100), 
          fn(25,alpha=5,beta=-5)) %>% #ggplot(data=.,aes(x,y,color=group))+geom_point()
  lme4::lmer(y~x+(1|group), data=.) %>% 
  summary









