library(tidyverse)


dry <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDrought_preBS2021-04-07 17:48:39.parquet") %>% 
  as_tibble()
hot <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef4_preBS2021-04-07 14:05:44.parquet") %>% 
  as_tibble()

d_fig <- bind_rows(dry %>% select(date_driest,ttr2) %>% mutate(mech='drought') %>% 
            rename(date = date_driest, ttr=ttr2) %>% filter(ttr>=365), 
          hot %>% select(date_fire1,ttr4) %>% mutate(mech='fire') %>% 
            rename(date = date_fire1, ttr = ttr4)
) %>% 
  mutate(hydro_year = year(date-months(3))) 

# d_fill <- d_fig %>% group_by(hydro_year, mech) %>% summarize(frac = sum(is.na(ttr)==F)/n()) %>% 
#   ungroup()
d_fill <- d_fig %>% group_by(hydro_year, mech) %>% summarize(frac = sum(is.na(ttr)==F)) %>% 
  ungroup()


inner_join(d_fig, d_fill) %>% 
  filter(hydro_year %in% 2002:2015) %>% 
  ggplot(data=.,aes(x=hydro_year,y=ttr,group=paste(hydro_year,mech),
                    color=mech,
                    fill=frac))+
  geom_boxplot(outlier.colour = NA)+
  geom_segment(data=tibble(fire_year=2000:2016) %>%
                 mutate(obs_days = 365*(2021-fire_year)),
               inherit.aes = F,
               aes(x=fire_year,
                   xend=fire_year+1,
                   y=obs_days,
                   yend=obs_days), lty=3)+
  annotate(geom = 'text',x = 2015.5, y=2750, label='observation window limit',
           angle=-52.5, size=5)+
  geom_hline(aes(yintercept=365),lty=3)+
  scale_color_manual(values=c("fire"='#CF0000',
                              "drought"="#AAAAFF"))+
  scale_fill_viridis_c(option='B',begin = 0.1,end=0.9)+
  # scale_fill_gradientn(colors=c("white","grey60"), 
  #                      # limits=c(0,1)
  #                      )+
  scale_x_continuous(limits=c(2001.5,2015.5),
                     breaks=seq(2002,2014,by=2), 
                     expand = c(0,0))+
  scale_y_continuous(limits=c(0,3000))+
  labs(x='Hydraulic Year', 
       y='Time to recover (days)', 
       fill='N')+
  theme_linedraw()
ggsave(filename = "figures/draft_fig_barplot_drought-2sd-ttr_fire-ttrDef4.png")
