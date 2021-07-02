library(tidyverse); 
library(data.table)
library(arrow); 
library(lubridate)
start_year <- 2002
end_year <- 2016

# proc lai ----------------
d_lai <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-05-24 16:20:29.parquet") %>% 
  .[,fire_year := year(date_fire1-months(3))]
d_lai <- merge(d_lai[fire_year>= start_year & fire_year<=end_year], 
      d_lai[fire_year>= start_year & fire_year<=end_year][,.(nobs=.N),by='fire_year'],
                by='fire_year')
d_lai <- d_lai[,.(val = .N, 
    nobs = unique(nobs)), 
  by=.(fire_year, ttr5_lai)] %>% 
  .[,val2 := val/nobs] %>% 
  .[order(ttr5_lai), csum := cumsum(val2), by=list(fire_year)] %>% 
  .[,days_past := as.double(ymd("2019-08-01")-ymd(paste(fire_year,1,1)))] %>% 
  .[,csum := ifelse(ttr5_lai > days_past, NA_real_, csum)] %>% 
  .[,vi:='LAI'] %>% 
  .[,ttr := ttr5_lai]


# proc ndvi ------------------
d_ndvi <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet") %>% 
  .[,fire_year := year(date_fire1-months(3))] %>% 
  rename(ttr5_ndvi = ttr5) %>% 
  as.data.table()
d_ndvi <- merge(d_ndvi[fire_year>= start_year & fire_year<=end_year], 
               d_ndvi[fire_year>= start_year & fire_year<=end_year][,.(nobs=.N),by='fire_year'],
               by='fire_year')
d_ndvi <- d_ndvi[,.(val = .N, 
                  nobs = unique(nobs)), 
               by=.(fire_year, ttr5_ndvi)] %>% 
  .[,val2 := val/nobs] %>% 
  .[order(ttr5_ndvi), csum := cumsum(val2), by=list(fire_year)] %>% 
  .[,days_past := as.double(ymd("2019-08-01")-ymd(paste(fire_year,1,1)))] %>% 
  .[,csum := ifelse(ttr5_ndvi > days_past, NA_real_, csum)] %>% 
  .[,vi:='NDVI'] %>% 
  .[,ttr := ttr5_ndvi]


# proc kndvi ------------------
d_kndvi <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_kndvi_ttrDef5_preBS2021-04-21 16:20:09.parquet") %>% 
  .[,fire_year := year(date_fire1-months(3))] %>% 
  rename(ttr5_kndvi = ttr5_kn) %>% 
  as.data.table()
d_kndvi <- merge(d_kndvi[fire_year>= start_year & fire_year<=end_year], 
                d_kndvi[fire_year>= start_year & fire_year<=end_year][,.(nobs=.N),by='fire_year'],
                by='fire_year')
d_kndvi <- d_kndvi[,.(val = .N, 
                    nobs = unique(nobs)), 
                 by=.(fire_year, ttr5_kndvi)] %>% 
  .[,val2 := val/nobs] %>% 
  .[order(ttr5_kndvi), csum := cumsum(val2), by=list(fire_year)] %>% 
  .[,days_past := as.double(ymd("2019-08-01")-ymd(paste(fire_year,1,1)))] %>% 
  .[,csum := ifelse(ttr5_kndvi > days_past, NA_real_, csum)] %>% 
  .[,vi:='kNDVI'] %>% 
  .[,ttr := ttr5_kndvi]


# proc CCI ------------------
d_cci <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_cci_ttrDef5_preBS2021-04-19 13:06:11.parquet") %>% 
  .[,fire_year := year(date_fire1-months(3))] %>% 
  as.data.table()
d_cci <- merge(d_cci[fire_year>= start_year & fire_year<=end_year], 
                 d_cci[fire_year>= start_year & fire_year<=end_year][,.(nobs=.N),by='fire_year'],
                 by='fire_year')
d_cci <- d_cci[,.(val = .N, 
                      nobs = unique(nobs)), 
                   by=.(fire_year, ttr5_cci)] %>% 
  .[,val2 := val/nobs] %>% 
  .[order(ttr5_cci), csum := cumsum(val2), by=list(fire_year)] %>% 
  .[,days_past := as.double(ymd("2019-08-01")-ymd(paste(fire_year,1,1)))] %>% 
  .[,csum := ifelse(ttr5_cci > days_past, NA_real_, csum)] %>% 
  .[,vi:='CCI'] %>% 
  .[,ttr := ttr5_cci]



# proc deltaT ------------------
d_deltat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_delta_t_ttrDef5_preBS2021-04-25 16:28:57.parquet") %>% 
  .[,fire_year := year(date_fire1-months(3))] %>% 
  rename(ttr5_deltat = ttr5_delta_t) %>% 
  as.data.table()
d_deltat <- merge(d_deltat[fire_year>= start_year & fire_year<=end_year], 
               d_deltat[fire_year>= start_year & fire_year<=end_year][,.(nobs=.N),by='fire_year'],
               by='fire_year')
d_deltat <- d_deltat[,.(val = .N, 
                  nobs = unique(nobs)), 
               by=.(fire_year, ttr5_deltat)] %>% 
  .[,val2 := val/nobs] %>% 
  .[order(ttr5_deltat), csum := cumsum(val2), by=list(fire_year)] %>% 
  .[,days_past := as.double(ymd("2019-08-01")-ymd(paste(fire_year,1,1)))] %>% 
  .[,csum := ifelse(ttr5_deltat > days_past, NA_real_, csum)] %>% 
  .[,vi:='deltaT'] %>% 
  .[,ttr := ttr5_deltat]


# proc tree cover ------------------
d_treecover <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_tree_cover_ttrDef5_preBS2021-04-20 16:21:42.parquet") %>% 
  .[,fire_year := fire_year1] %>% 
  rename(ttr5_treecover = ttr5_tree_cover) %>% 
  as.data.table() %>% 
  .[,ttr5_treecover := ttr5_treecover*365]
d_treecover <- merge(d_treecover[fire_year>= start_year & fire_year<=end_year], 
                  d_treecover[fire_year>= start_year & fire_year<=end_year][,.(nobs=.N),by='fire_year'],
                  by='fire_year')
d_treecover <- d_treecover[,.(val = .N, 
                        nobs = unique(nobs)), 
                     by=.(fire_year, ttr5_treecover)] %>% 
  .[,val2 := val/nobs] %>% 
  .[order(ttr5_treecover), csum := cumsum(val2), by=list(fire_year)] %>% 
  .[,days_past := as.double(ymd("2019-08-01")-ymd(paste(fire_year,1,1)))] %>% 
  .[,csum := ifelse(ttr5_treecover > days_past, NA_real_, csum)] %>% 
  .[,vi:='treecover'] %>% 
  .[,ttr := ttr5_treecover]



p_out <- rbindlist(list(d_ndvi,d_kndvi,d_cci,d_lai,d_deltat, d_treecover), 
          fill=T) %>% 
  filter(fire_year > start_year) %>%
  ggplot(data=., aes(ttr, csum, color=vi))+
  geom_line(lwd=1)+
  scale_x_continuous(breaks=c(400,800,1200,1600))+
  coord_cartesian(xlim=c(365,2000), 
                  expand=F)+
  scale_color_viridis_d(option='H',direction = -1)+
  labs(x="Time to Recover (days)", 
       y="Fraction Recovered",
       color='Metric')+
  facet_wrap(~fire_year,ncol=3)+
  guides(color=guide_legend(nrow=1))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom'); p_out
  
ggsave(p_out, 
       filename = "figures/plot_vi-cumulative-recovery-trajectory.png",
       width=15,
       height=20, 
       units='cm',
       dpi=350)
# 
# p_out$theme$text
# 
# d_treecover %>% ggplot(data=.,aes(ttr, csum,color=factor(fire_year)))+
#   geom_line()+
#   scale_color_viridis_d(option='H')
# 
# 
# 
# tmp <- merge(out[fire_year>2001 & fire_year<=2016], out[,.(nobs = .N), by='fire_year']) 
# vec_days <- unique(tmp$ttr5_lai) %>% na.omit() %>% sort
# 
# merge(data.table(post_days = vec_days), 
#       tmp, by='post_days', allow.cartesian = TRUE)
# 
# tmp[][,.(val = .N, 
#        nobs = unique(nobs)), 
#     by=.(fire_year, ttr5_lai)] %>% 
#   .[,val2 := val/nobs] %>% 
#   .[order(ttr5_lai), csum := cumsum(val2), by=list(fire_year)] %>% 
#   .[,days_past := as.double(ymd("2019-08-01")-ymd(paste(fire_year,1,1)))] %>% 
#   .[,csum := ifelse(ttr5_lai > days_past, NA_real_, csum)] %>% 
#   ggplot(data=.,aes(ttr5_lai, csum,color=factor(fire_year)))+
#   geom_line()+
#   scale_color_viridis_d(option='H')
# 
# 
# 
# 
# 
# d_ndvi <- read_parquet("../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_ndvi_ttrDef5_preBS.parquet") %>% 
#   mutate(vi = 'ndvi')
# d_cci <- read_parquet("../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_cci_ttrDef5_preBS.parquet") %>% 
#   mutate(vi = 'cci') 
# d_deltat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_delta_t_ttrDef5_preBS.parquet") %>% 
#   mutate(vi = 'deltaT')
# d_kndvi <- read_parquet("../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_kndvi_ttrDef5_preBS.parquet") %>% 
#   mutate(vi = 'kndvi')
# d_treecover <- read_parquet("../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_tree_cover_ttrDef5_preBS.parquet") %>% 
#   mutate(vi = 'tree_cover')
# d_lai <- read_parquet("../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_lai_ttrDef5_preBS.parquet") %>% 
#   mutate(vi = 'lai')
# 
# bind_rows(d_ndvi, d_cci, d_deltat, d_kndvi, d_treecover, d_lai) %>% 
#   mutate(days_since_end = (ymd("2019-08-01")-ymd(paste(hydro_year,1,1)))) %>% 
#   mutate(val = 
#          ifelse(post_days  > days_since_end, NA_real_, val)) %>% 
#   filter(hydro_year>2001 & hydro_year<=2017) %>%
#   ggplot(data=., aes(post_days,val,color=vi))+
#   geom_line(lwd=2)+
#   geom_hline(aes(yintercept=0.95),lty=3)+
#   scale_x_continuous(limits=c(375,3000), 
#                      expand=c(0,0))+
#   scale_y_continuous(limits=c(0,1), 
#                      expand=c(0,0))+
#   scale_color_viridis_d(option='H')+
#   # scico::scale_color_scico_d(end=0.95)+
#   facet_wrap(~hydro_year,ncol = 3)+
#   theme_minimal()
# 
# 
# out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-05-24 16:20:29.parquet")
# out[,hydro_year := year(date_fire1-months(3))]
# out[,date_recovered :=date_fire1+days(ttr5_lai)]
# out[date_recovered<=ymd("2021-03-01")] %>% 
#   .[date_recovered>=ymd("2016-01-01")] %>% 
#   ggplot(data=.,aes(date_recovered))+geom_bar()
# 
# d3 <- expand_grid(out, post_days=seq.int(30,3000,by=30)) %>% 
#   mutate(hydro_year = year(date_fire1 - months(3))) %>% 
#   mutate(ttr5_lai = ifelse(is.na(ttr5_lai)==T, 9999, ttr5_lai)) %>%
#   mutate(recovered = ifelse(post_days >= ttr5_lai, 1, 0)) %>% 
#   # mutate(recovered = ifelse(is.na(ttr5_lai)==T, 0, 1)) %>% 
#   group_by(post_days,hydro_year) %>% 
#   summarize(val = sum(recovered,na.rm=TRUE), 
#             nobs = n()) %>% 
#   mutate(val = val/nobs) %>%
#   mutate(post_date = ymd(paste(hydro_year,1,1))+days(post_days)) %>% 
#   ungroup()
# 
# 
# d3 <- expand_grid(out[hydro_year>2001 & hydro_year<=2015], 
#             post_days=seq.int(30,3000,by=30))%>% 
#   as_tibble() %>% 
#   mutate(recovered = ifelse(is.na(ttr5_lai)==F, 1,0)) %>% 
#   arrange(date_recovered) %>% 
#   group_by(hydro_year) %>% 
#   mutate(nobs = n()) %>% 
#   mutate(csum = cumsum(recovered)) %>% 
#   mutate(val = csum/nobs)
# 
# d3 <- expand_grid(out[hydro_year>2001 & hydro_year<=2015], 
#                   post_days=seq.int(30,3000,by=30))%>% 
#   as_tibble() %>% 
#   mutate(recovered = ifelse(is.na(ttr5_lai)==F, 1,0)) %>% 
#   arrange(post_days) %>% 
#   group_by(hydro_year,post_days) %>% 
#   summarize(nobs = n(), 
#             csum = sum(recovered)) %>% 
#   ungroup()
# 
# d3 %>% 
#   ggplot(data=.,aes(post_days, csum, color=factor(hydro_year)))+
#   geom_point()
# 
# 
# d3 %>% group_by(hydro_year, post_days) %>% 
#   summarize(val = mean(csum/nobs))
# 
# d3 %>%   
#   ggplot(data=.,aes(date_recovered, csum/nobs, group=year,color=year))+
#   geom_point()+
#   geom_hline(aes(yintercept=1))
# 
# 
# out[hydro_year==2015] %>% ggplot(aes(date_recovered))+geom_bar()
# out[hydro_year==2015] %>% 
#   lazy_dt() %>% 
#   group_by(date_recovered) %>% 
#   summarize(val = n()/1264) %>% 
#   ungroup() %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(date_recovered,val))+
#   geom_line()
# 
# out[hydro_year==2015] %>% as_tibble() %>% 
#   mutate(recovered = ifelse(is.na(ttr5_lai)==F, 1,0)) %>% 
#   arrange(date_recovered) %>% 
#   group_by(hydro_year) %>% 
#   mutate(nobs = n()) %>% 
#   mutate(csum = cumsum(recovered)) %>% 
#   ggplot(data=.,aes(date_recovered, csum/nobs))+
#   geom_point()+
#   geom_hline(aes(yintercept=0.876))
# 
# 
# out[hydro_year>2001 & hydro_year<=2015] %>% 
#   as_tibble() %>% 
#   mutate(recovered = ifelse(is.na(ttr5_lai)==F, 1,0)) %>% 
#   arrange(date_recovered) %>% 
#   group_by(hydro_year) %>% 
#   mutate(nobs = n()) %>% 
#   mutate(csum = cumsum(recovered)) %>% 
#   ggplot(data=.,aes(date_recovered, csum/nobs, group=year,color=year))+
#   geom_point()+
#   geom_hline(aes(yintercept=1))
# 
# out[hydro_year>2001 & hydro_year<=2015] %>% 
#   as_tibble() %>% 
#   mutate(recovered = ifelse(is.na(ttr5_lai)==F, 1,0)) %>% 
#   arrange(date_recovered) %>% 
#   group_by(hydro_year) %>% 
#   mutate(nobs = n()) %>% 
#   mutate(csum = cumsum(recovered)) %>% 
#   arrange(hydro_year,date_recovered) %>% 
#   ggplot(data=.,aes(date_recovered, csum/nobs, 
#                     group=factor(year),color=factor(year)))+
#   geom_line()+
#   geom_hline(aes(yintercept=1))
# 
# 
# out[hydro_year==2015]$ttr5_lai %>% is.na %>% table
# 
# 1 - 156/(1108+156)
# 
# d3 %>% filter(hydro_year==2015) %>% 
#   ggplot(aes(post_days,val))+
#   geom_line()+
#   geom_hline(aes(yintercept=1))
# 
# expand_grid(out[hydro_year==2015], post_days=seq.int(30,3000,by=30)) %>% 
#   mutate(ttr5_lai = ifelse(is.na(ttr5_lai)==T, 9999, ttr5_lai)) %>%
#   mutate(recovered = ifelse(post_days >= ttr5_lai, 1, 0)) %>% 
#   group_by(post_days) %>% 
#   summarize(val = sum(recovered)) %>% 
#   ungroup() %>% 
#   # filter(post_days > 800) %>% 
#   View
# 
# out[hydro_year==2017]$ttr5_lai %>% is.na %>% table
# out[hydro_year==2017]$date_fire1+days(810)
# 
# d3 %>% filter(hydro_year==2017)
# d3 %>% filter(hydro_year==2015) %>% ggplot(data=.,aes(post_days,val))+geom_line()
#   
# d3 %>% filter(hydro_year==2017) %>% 
#   filter(date_recovered >= ymd("2019-01-01") & is.na(date_recovered)==F) %>% 
#   pull(recovered)
# out[,hydro_year := year(date_fire1 - months(3))]
# out[hydro_year==2017]
# 
# 
# d3 %>% filter(hydro_year==2007)
# d3 %>% filter(hydro_year==20115) %>% 
#   mutate(post_date = ymd("2016-09-01")+days(post_days)) %>% 
#   ggplot(data=.,aes(post_date,val))+
#   geom_line()
# d3 %>% filter(hydro_year==2017) %>% filter(post_date>ymd("2018-01-01")) %>% View
# out[between(date_fire1, ymd("2016-09-01"),ymd("2017-03-01"))] %>% 
#   ggplot(data=.,aes(date_recovered))+
#   geom_bar()
# d3 %>% filter(hydro_year==2017) %>% 
#   mutate(post_date = ymd("2016-09-01")+days(post_days)) %>% 
#   filter(post_date >= ymd("2019-01-01"))
# 
# d3 %>% mutate(vi='lai') %>% 
#   mutate(days_since_end = (ymd("2021-03-01")-ymd(paste(hydro_year,1,1)))) %>% 
#   mutate(val =
#            ifelse(post_days  > days_since_end, NA_real_, val)) %>%
#   filter(hydro_year>2001 & hydro_year<=2017) %>%
#   ggplot(data=., aes(post_days,val,color=vi))+
#   geom_line(lwd=2)+
#   geom_hline(aes(yintercept=0.95),lty=3)+
#   scale_x_continuous(limits=c(375,3000), 
#                      expand=c(0,0))+
#   scale_y_continuous(limits=c(0,1), 
#                      expand=c(0,0))+
#   scale_color_viridis_d(option='H')+
#   # scico::scale_color_scico_d(end=0.95)+
#   facet_wrap(~hydro_year,ncol = 3)+
#   theme_minimal()
# 
# d3 %>% filter(hydro_year==2017) %>% ggplot(data=.,aes(post_days,val))+geom_line()
# out[between(date_fire1, ymd("2016-09-01"),ymd("2017-03-01"))] %>% 
#   ggplot(data=.,aes(date_recovered))+
#   geom_bar()
# 
# 
# 
# d_ndvi %>% 
#   # filter(hydro_year==2001) %>% 
#   ggplot(data=.,aes(post_days, y=val,group=hydro_year,color=hydro_year))+
#   geom_vline(aes(xintercept=365),lty=3)+
#   geom_line()+
#   geom_hline(aes(yintercept=0.95),lty=3)+
#   scale_color_viridis_c(option='B', 
#                         end=0.9, 
#                         direction = -1)+
#   scale_x_continuous(limits=c(370,3000), 
#                      expand=c(0,0))
# 
# d_lai %>% 
#   # filter(hydro_year==2001) %>% 
#   ggplot(data=.,aes(post_days, y=val,group=hydro_year,color=hydro_year))+
#   geom_vline(aes(xintercept=365),lty=3)+
#   geom_line()+
#   geom_hline(aes(yintercept=0.95),lty=3)+
#   scale_color_viridis_c(option='B', 
#                         end=0.9, 
#                         direction = -1)+
#   scale_x_continuous(limits=c(370,3000), 
#                      expand=c(0,0))
# 
#   
#   
#   
#   
# d3 %>% 
#   filter(between(hydro_year,2001,2015)) %>% 
#   ggplot()+
#   geom_rect(aes(xmin=post_days,xmax=post_days+30,
#                 ymin=hydro_year-0.5,
#                 ymax=hydro_year+0.5,
#                 fill=val))+
#   geom_point(data=. %>% filter(val>=0.95) %>% group_by(hydro_year) %>% 
#                filter(post_days==min(post_days,na.rm=TRUE)) %>% 
#                ungroup(), 
#              inherit.aes = FALSE, 
#              aes(x=post_days,y=hydro_year),shape=20, col='#5555FF',size=3)+
#   scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
#                        oob=scales::squish)+
#   # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
#   scale_x_continuous(limits=c(400,2600), 
#                      expand=c(0,0))+
#   scale_y_continuous(expand=c(0,0), 
#                      limits=c(2000.5,2015.5), 
#                      breaks=seq(2001,by=2,length.out=8))+
#   labs(x='Days post fire', 
#        y='Year of Bushfire',
#        fill='NDVI Fraction Recovered   ')+
#   guides(fill=ggplot2::guide_colorbar(title.position = 'left',
#                                       title.hjust = 1000))+
#   theme_linedraw()+
#   theme(legend.position = 'bottom', 
#         legend.key.width = unit(1.5,'cm'), 
#         legend.key.height = unit(0.2,'cm'),
#         panel.grid = element_blank())
