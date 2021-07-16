pacman::p_load(tidyverse, data.table, lubridate, patchwork)

dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                           col_select = c("x","y","id","date","slai","slai_anom_12mo","malai"))
dat[,`:=`(slai_12mo = slai_anom_12mo+malai)]
sdat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS_2021-06-05 13:01:38.parquet")

# preprocess logistic function fits ------
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-1mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-08 21:41:09.parquet")
fits <- fits[isConv==TRUE][r2>0][,ldk:=L0/K]
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
fits[,fire_year := lubridate::year(date_fire1 - months(3))]



# PLOTS -------------------------------------------------------------------
# Composite plot: ttr5~L0/K by year ----------------
v_top8 <- fits[,.(nobs=.N),by='fire_year'][order(-nobs)]$fire_year[1:8]

pan_b <- fits[fire_year %in% v_top8] %>%  
  as_tibble() %>% 
  mutate(fire_year=factor(fire_year)) %>% 
  ggplot(data=.,aes(ldk, ttr5_lai,color=fire_year,group=fire_year))+
  # geom_point()+
  geom_smooth(method='gam',
              formula=y~s(x,k=5,bs='cs'),
              method.args=list(select=TRUE))+
  scale_color_viridis_d(option='H')+
  # scico::scale_color_scico_d()+
  # scale_color_manual(values=pals::glasbey(8))+
  labs(       x=expression(paste("Fraction remaining leaf area:"~italic(L[0]/K),"  (",m**2/m**2,")")), 
       y='Time to Recover (days)',
       color='Fire year')+
  coord_cartesian(xlim=c(0,1),
                  ylim=c(0,1500),
                  expand=F)+
  theme_linedraw()+
    theme(panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill=NA),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = c(0.25,0), 
        legend.justification = c(0.25,0), 
    legend.direction = 'horizontal'); pan_b

pan_c <- expand_grid(tibble(K=c(3,3,3),
                   r=c(0.01, 0.0025, 0.0075),
                   L0=c(0.3,1.5,0.8*3)), 
            post_days = floor(seq(1,2000,length.out=100))) %>% 
  mutate(reduction = round(100*(1 - L0/K),digits=2)) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days)) ) %>% 
  mutate(recovered = ifelse(between(pred, 0.975*K, 0.98*K),pred,NA)) %>% 
  ggplot(data=.,aes(post_days,pred,
                    group=paste(r,L0,K), 
                    color=factor(r)))+
  geom_line(lwd=1)+
  # geom_point(aes(post_days,recovered))+
  coord_cartesian(expand=F, xlim=c(0,2030),
    ylim=c(0,3.1))+
  labs(x='Days after fire',
       y='LAI (m²/m²)', 
       color=expression(paste(italic(r))))+
  # scale_color_brewer(palette='Set1')+
  scale_color_manual(values=c("#cf0000", "#1505fc", "#05ecfc"))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill=NA),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = c(0.99,0), 
        legend.justification = c(0.99,0), 
    legend.direction = 'horizontal'); pan_c

(pan_a <- fits %>% 
  as.data.table() %>% 
  .[fire_year>=2001 & fire_year<=2014] %>% 
  .[isConv==TRUE] %>% 
  .[r<0.05] %>% 
  .[r2>0] %>% 
  .[L0<K] %>% 
  # .[r2>0.75] %>% 
  ggplot(data=.,aes(L0/K,r,color=pred_ttr))+
  geom_point(alpha=0.25,size=0.25)+
  geom_smooth(fullrange=T,
              aes(L0/K,r),
              formula=y~s(x,bs='cs',k=5),
              color="white",
              weight=5,
              level=0.999)+
  geom_smooth(fullrange=T,
              aes(L0/K,r),
              formula=y~s(x,bs='cs',k=5),
              color="#cf0000",
              weight=2.5,
              level=0.999)+
  scale_color_viridis_c(option='F',
                        limits=c(365,2000),
                        direction = 1,
                        end = 0.85,
                        breaks=c(365,730,1095,1460,1825),
                        labels=c("≤ 1","2","3","4","≥5"),
                        oob=scales::squish
  )+
  scale_x_continuous(limits=c(0,1),expand=c(0,0))+
  scale_y_continuous(limits=c(0,0.05),expand=c(0,0))+
  labs(y=expression(paste("Growth rate: ",italic(r),"  (m² day¯¹)")), 
       x=expression(paste("Fraction remaining leaf area:"~italic(L[0]/K),"  (",m**2/m**2,")")), 
       color="TTR (years)   ")+
  theme_linedraw()+
    theme(panel.grid = element_blank(), 
      legend.direction = 'horizontal',
      legend.position = c(0.5,0.99),
      legend.justification = c(0.5,0.99),
      legend.background = element_rect(fill=alpha(colour = 'white',alpha=0.85)),
      legend.key.height = unit(2,'mm'), 
      legend.key.width = unit(13,'mm')))

((pan_b|pan_c)+plot_layout(widths = c(1.3,1)))/pan_a+plot_layout(heights = c(1.35,1))+
  plot_annotation(tag_levels = 'a',
                                    tag_prefix = '(',
                                    tag_suffix = ')')
scale_factor <- 0.333
ggsave(
       filename = 'figures/plot_composite_r-L0K_ttr-L0K_logFunDiagram.png', 
       width=815*scale_factor, # 815
       height=815*scale_factor, # 915
       units='mm',
       dpi=300)


# (fits %>% 
#   as.data.table() %>% 
#   .[fire_year>=2001 & fire_year<=2014] %>% 
#   .[isConv==TRUE] %>% 
#   .[r<0.05] %>% 
#   .[r2>0] %>% 
#   .[L0<K] %>% 
#   # .[r2>0.75] %>% 
#   ggplot(data=.,aes(L0/K,r,color=cut_interval(r2,4)))+
#   geom_point(alpha=0.25,size=0.25)+
#   geom_smooth(fullrange=F,
#               # aes(L0/K,r),
#               # formula=y~s(x,bs='cs',k=5),
#               # color="#cf0000",
#               # weight=2.5,
#               # level=0.999
#     )+
#   scale_color_viridis_d(option='H')+
#   #     scale_color_viridis_c(option='F',
#   #                       limits=c(365,2000),
#   #                       direction = 1,
#   #                       end = 0.85,
#   #                       breaks=c(365,730,1095,1460,1825),
#   #                       labels=c("≤ 1","2","3","4","≥5"),
#   #                       oob=scales::squish
#   # )+
#   scale_x_continuous(limits=c(0,1),expand=c(0,0))+
#   scale_y_continuous(limits=c(0,0.05),expand=c(0,0))+
#   labs(y=expression(paste("Growth rate: ",italic(r),"  (m² day¯¹)")), 
#        x=expression(paste("Fraction remaining leaf area:"~italic(L[0]/K),"  (",m**2/m**2,")")), 
#        color="TTR (years)   ")+
#   theme_linedraw()+
#     theme(panel.grid = element_blank(), 
#       legend.direction = 'horizontal',
#       legend.position = c(0.5,0.99),
#       legend.justification = c(0.5,0.99),
#       legend.background = element_rect(fill=alpha(colour = 'white',alpha=0.85)),
#       legend.key.height = unit(2,'mm'), 
#       legend.key.width = unit(13,'mm')))
# 
# fits$r2 %>% hist
# fits$L0 %>% hist
# fits$K %>% hist
# fits %>% 
#   mutate(val = L0/K) %>% 
#   pull(val) %>% hist
