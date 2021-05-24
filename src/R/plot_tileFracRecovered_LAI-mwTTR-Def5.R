library(tidyverse); 

# FASTER Method -----------------------------------------------------------
out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")
d3 <- expand_grid(out, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(val = sum(post_days>ttr5_lai,na.rm=TRUE)/n()) %>% 
  ungroup()

d3 %>% write_parquet(.,
                     sink="../data_general/proc_data_Oz_fire_recovery/cumulativeRecovery_lai_ttrDef5_preBS.parquet")

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
       fill='LAI Fraction Recovered   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank())

ggsave(filename = 'figures/figure_cumulativeRecovered_byYear_TTR-Def5-lai.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)
