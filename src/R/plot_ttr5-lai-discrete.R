library(tidyverse); 

d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")
d_rf <- d_rf %>% filter(date_fire1 <= lubridate::ymd("2017-03-01"))
d_rf <- d_rf %>% mutate(ttr_ocat = factor(ttr_ocat,levels = 1:5, labels = 1:5, ordered = TRUE))

d_rf %>% filter(
  dom_sp.x %in% (sp_count %>% 
                   filter(rank <= 40) %>% pull(dom_sp))
  ) %>% 
  mutate(dom_sp = factor(dom_sp.x)) %>% 
  mutate(ttr_cat = fct_rev(ttr_ocat)) %>% 
  # select(dom_sp.x, ttr_ocat, .pred_class) %>% 
  # mutate(.pred_class = factor(.pred_class, levels=1:5,labels=1:5,ordered = TRUE)) %>% 
  # pivot_longer(cols=c("ttr_ocat",".pred_class")) %>% 
  ggplot(data=.,aes(y=dom_sp, 
                    fill=ttr_ocat))+
  geom_bar(position=position_fill(reverse=T))+
  scale_y_discrete(limits = rev)+
  scale_x_continuous(expand=c(0,0), 
                     limits=c(0,1))+
  scale_fill_viridis_d(option='B', end=0.95, direction = 1, 
                       labels=c("≤ 1","2","3","4","≥ 5"))+
  labs(y=NULL,
       x='proportion',
       fill='Years to Recover')+
  # guides(fill = guide_colorsteps())
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(0,0.5,0,0),'cm'))

ggsave(filename = "figures/prop-plot_dom-sp-ttr5-lai_5year-discrete_2001-2017.png", 
       dpi=350, 
       height=20, 
       width=15, 
       units='cm')


d_rf %>% 
  mutate(ttr_cat = fct_rev(ttr_ocat)) %>% 
  ggplot(data=.,aes(y=vc_name_f, 
                    fill=ttr_ocat))+
  geom_bar(position=position_fill(reverse=T))+
  scale_y_discrete(limits = rev)+
  scale_x_continuous(expand=c(0,0), 
                     limits=c(0,1))+
  scale_fill_viridis_d(option='B', end=0.95, direction = 1, 
                       labels=c("≤ 1","2","3","4","≥ 5"))+
  labs(y=NULL,
       x='proportion',
       fill='Years to Recover')+
  # guides(fill = guide_colorsteps())
  theme(legend.position = 'bottom', 
        plot.margin = unit(c(0,0.5,0,0),'cm'))

