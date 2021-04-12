library(data.table); 
library(tidyverse);
library(stars);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(dtplyr)

# Data import ---------------------------------------------------
dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet", 
                           col_select = c("x","y","id","date","ndvi_anom"))
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet")

#! Only fitting locations where the recovery was at least one year
sdat <- sdat[is.na(ttr5)==FALSE][date_fire1<ymd('2015-01-01')][ttr5>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_fire1,ttr5)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>%
  mutate(recovery_date = date_fire1+days(ttr5)) %>% 
  group_by(x,y,id) %>% 
  filter(date > date_fire1) %>% 
  filter(date <= recovery_date + years(3)) %>% 
  ungroup() %>%
  mutate(post_days = as.double(date - date_fire1)) %>% 
  as.data.table() 



out <- arrow::read_parquet("outputs/weibull4Param_fits_1burn_2001-2014fires_2021-04-08 12:48:07.parquet")


sdat <- sdat[,`:=`(id=as.integer(id))]
out <- out[,`:=`(id=as.integer(id))]
dd <- merge(out, sdat, by=c("id")) %>% as.data.table



sum(test$lrc < -10,na.rm=TRUE)/dim(test)[1] # fraction with long recovery
sum(test$r2 >0.2,na.rm=TRUE)/dim(test)[1] # fraction with reasonable fit


vec_post_days <- sort(unique(mdat$post_days))
test[isConv==TRUE][Drop>0][r2>0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
  .[sample(.N,1000)] %>% 
  expand_grid(.,
              post_days=c(0,vec_post_days)) %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
  # left_join(., mdat, by=c('id','post_days')) %>% 
  as_tibble() %>% 
  mutate(recovered = ifelse(pred > 0, 1,0)) %>%
  filter(recovered == 0) %>% 
  mutate(form = ifelse(lrc <= -10, 'slow','fast')) %>% 
  group_by(form) %>% 
  summarize(val = mean(ttr5))
  ggplot(data=.,aes(post_days,pred,group=id,color=ttr5))+
  # ggpointdensity::geom_pointdensity()+
  geom_line(lwd=0.05)+
  geom_hline(aes(yintercept=0),col='#CF0000')+
  # geom_vline(inherit.aes=F, 
  #            aes(xintercept=mean(ttr5),
  #                group=form))+
  # geom_abline()+
  # geom_smooth(inherit.aes = F, 
  #             aes(post_days,pred))+
  scico::scale_color_scico(end=0.7,
                           begin=0,
                           palette = 'batlow')+
  scale_x_continuous(limits=c(0,2500), 
                     expand=c(0,0))+
  labs(x='post_days', 
       y='NDVI anom',
       title=' Bushfires')+
  facet_wrap(~form,ncol = 1)+
  theme_linedraw()+
  theme()

scico::scico_palette_show()
test[ttr5>2000][min_nbr_anom < -0.5]


fn <- ecdf(sdat$ttr5)
tibble(post_days = c(0,vec_post_days)) %>% 
  mutate(frac = fn(post_days)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,ymin=0,ymax=1,fill=frac))+
  scale_fill_viridis_c()

vec_post_days <- mdat$post_days %>% unique %>% sort
d2 <- expand_grid(sdat,post_days=vec_post_days) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(post_days>ttr5)/n()) %>% 
  ungroup()


d2 %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=val))+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish)+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(365,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2014.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Year of Bushfire',
       fill='Fraction Recovered   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1.5,'cm'), 
        legend.key.height = unit(0.2,'cm'))
ggsave(filename = 'figures/figure_cumulativeRecovered_byYear_TTR-Def5.png', 
       width=15, 
       height=8, 
       units='cm',
       dpi=350)
  





