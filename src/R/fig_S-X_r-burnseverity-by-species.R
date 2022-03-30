tmp <- merge(fits,rank_species,by='species')[,`:=`(bs=ifelse(bs<0,0,bs))][r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.5][,hi_lo := ifelse(bs <= 0.5,'lo','hi')][] 

tmp <- tmp[between(rank,1,80)][,`:=`(rank_group=cut_interval(rank,length=10))] %>% 
  .[,`:=`(species_short = str_replace(species,"Eucalyptus","E."))] %>% 
  .[,`:=`(species_short = str_replace(species_short,"Angophora","A."))] %>% 
  .[,`:=`(species_short = str_replace(species_short,"Corymbia","C."))]

dfs <- split(tmp, f=tmp$rank_group)
gg_l <- lapply(dfs, function(x){
  ggplot(data=x, aes(bs,r,color=species_short))+
    geom_smooth(method='gam',
    formula = y~s(x,bs='cs',k=5),
    # formula = y~s(x,bs='ad',k=10),
    se=F)+
  scale_color_viridis_d(option='H')+
  labs(color=NULL, 
    x=expression(paste(Burn~Severity~(1-LAI["post-fire"]/LAI['norm'])))
    )+
  coord_cartesian(xlim=c(0,1.05),
    ylim=c(0,0.018),
    expand = F)+
  facet_wrap(~rank_group,labeller = label_both)+
  theme_linedraw()+
  theme(
    plot.margin = margin(0,0,0,0),
    legend.position = 'right',
    # legend.position = c(0.01,0.99),
    # legend.justification = c(0.01,0.99),
    legend.text = element_text(size=7),
    legend.key.height = unit(0.75,'line'))})
p_out <- wrap_plots(gg_l, ncol=2)
ggsave(p_out,filename = 'figures/fig_S-X_r-burnseverity-by-species.png',
  height = 25,
  width = 25,
  units='cm',
  dpi=350)

