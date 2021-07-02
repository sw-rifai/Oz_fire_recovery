pacman::p_load(tidyverse, 
               data.table)
out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-06-25 12:26:19.parquet")

vec_names <- out[,102:210] %>% names
vec_names

mout <- melt(id.vars=c("x","y"),
     out[,99:210][,-c("predict")], 
     variable.name = 'species')

idx <- seq(1,length(vec_names),by=9)
i <- 1
for(i in 1:length(idx)){
  print(i)
 p_out <-  ggplot(data=mout[species%in%vec_names[idx[i]:idx[i+1]]],
         aes(x,y,fill=value))+
    geom_sf(data=oz_poly, 
            inherit.aes = F, 
            fill='grey20',
            color='black')+
    geom_raster()+
    geom_sf(data=oz_poly, 
            inherit.aes = F, 
            fill=NA,
            color='grey80')+
    scale_fill_gradientn(colors=c('black',pals::turbo(12)))+
    # scale_fill_viridis_c(option='H',
    #                      # n.breaks=6, 
    #                      limits=c(0,1))+
    coord_sf(ylim=c(-38.75,-28), 
             xlim=c(145.5,153.5))+
    scale_x_continuous(breaks=c(146,149,152))+
    labs(x=NULL, 
         y=NULL, 
         fill='prob.' 
         # title=str_replace(eval(vec_names[i]),"\\."," ")
    )+
    facet_wrap(~species,nrow=2)+
    # guides(fill = guide_colorbar(ncol=1))+
    theme(
      #  legend.position = c(1,0),
      # legend.justification = c(1,0),
      strip.background = element_rect(fill='grey90'),
      panel.background = element_rect(fill='lightblue'), 
      panel.grid = element_blank())
im_name <- paste0(c(str_replace(c("sdm_eucs_xgboost_",vec_names[idx[i]],"_",vec_names[idx[i+1]]),pattern="\\.",replacement = "_"),".png"),collapse = "")
ggsave(p_out, 
       filename=paste0('figures/euc_sdm/',im_name), 
         width=16*2,
         height=9*2,
         units='cm',
         dpi=350)
}



out[,c('x','y',vec_names[i]), with=F] %>% 
  ggplot(data=.,aes(x,y,fill= get(vec_names[i]) ))+
  geom_sf(data=oz_poly, 
          inherit.aes = F, 
          fill='grey20',
          color='black')+
  geom_raster()+
  geom_sf(data=oz_poly, 
          inherit.aes = F, 
          fill=NA,
          color='grey80')+
  scale_fill_viridis_b(option='F',
                       n.breaks=6, 
                       limits=c(0,1))+
  coord_sf(ylim=c(-38.75,-28), 
           xlim=c(145.5,153.5))+
  labs(x=NULL, 
       y=NULL, 
       fill='prob.', 
       title=str_replace(eval(vec_names[i]),"\\."," "))+
  # guides(fill = guide_colorbar(ncol=1))+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        panel.background = element_rect(fill='lightblue'), 
        panel.grid = element_blank())

