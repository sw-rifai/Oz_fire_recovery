pacman::p_load(tidyverse, data.table, patchwork)

fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-22 09:37:07.parquet")
ba <- arrow::read_parquet('outputs/mcd64_conComp_2001_2020.parquet')
ba <- ba %>% rename(date_fire1 = date)
ba <- merge(fits, ba, by=c("x","y","date_fire1"))
ba_rank <- ba[,.(ba_m2 = median(ba_m2)), by=label][,size_rank := frank(-ba_m2,ties.method = 'first')][order(size_rank)]
ba_rank
vec_sel <- ba_rank[size_rank<=15]$label

pan_lf <- merge(ba, ba_rank, by='label') %>% 
  filter(size_rank <= 15) %>% as.data.table() %>%  
  ggplot(data=.,aes(x,y,fill=pred_ttr))+
  geom_raster()+
  scale_fill_viridis_c(option='B', limits=c(365,2000), oob=scales::squish, 
        breaks=c(365, 365*2, 365*3, 365*4, 365*5), 
    labels=c("≤ 1", " 2", " 3", " 4", "≥ 5"))+
  labs(x=NULL, y=NULL, 
    title="Time-to-recover: Logistic Function Estimate",
    fill='years   ')+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)+
  facet_wrap(~size_rank, scales = 'free',ncol=5,labeller = label_both)+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='grey90')); pan_lf

pan_mw <- merge(ba, ba_rank, by='label') %>% 
  filter(size_rank <= 15) %>% as.data.table() %>% 
  ggplot(data=.,aes(x,y,fill=ttr5_lai))+
  geom_raster()+
  scale_fill_viridis_c(option='B', limits=c(365,2000), oob=scales::squish, 
        breaks=c(365, 365*2, 365*3, 365*4, 365*5), 
    labels=c("≤ 1", " 2", " 3", " 4", "≥ 5"))+
  labs(x=NULL, y=NULL, 
    title="Time-to-recover: 12-month Moving Window Estimate",
    fill='years   ')+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)+
  facet_wrap(~size_rank, scales = 'free',ncol=5,labeller = label_both)+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='grey90'))

pan_mw/pan_lf+plot_annotation(tag_levels = 'a',tag_prefix = '(',tag_suffix = ')')+
  plot_layout(guides = 'collect')&theme(
    legend.position = 'bottom', 
    legend.key.height = unit(0.2,'cm'),
    legend.key.width = unit(0.6,'cm'))
scale_factor <- 3
ggsave(filename = 'figures/map_ttr-comparison_15-large-fires.png', 
  height=10*scale_factor,
  width=9*scale_factor,
  units='cm',
  dpi=500)
