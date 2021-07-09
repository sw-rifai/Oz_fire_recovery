pacman::p_load(tidyverse, 
               data.table)

dir.create(paste0('figures/euc_sdm_v0.3_',Sys.Date()))

# ML predictions ---------------------
out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-07 14:54:12.parquet")


# nvis ---------------------
nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif", 
                          proxy=F)
names(nvis) <- 'vc_code'
nvis_codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("vc_code","veg_class_descrip")) %>% 
  mutate(vc_name = as.factor(veg_class_descrip)) %>% 
  select(-veg_class_descrip)


# ALA MQ cleaned --------------
ala_mq <- fread("../data_general/ALA/Gallagher_cleaned/euc_occurrence_clean.csv") %>% 
  rename( 
    x=longitude,
    y=latitude) %>% 
  as.data.table() %>% 
  .[x>=140 & y <= -23 & y>=-40]

vv <- st_extract(nvis, 
           st_as_sf(ala_mq[,.(x,y)], coords=c("x","y"), crs=4326))
ala_mq <- bind_cols(ala_mq, st_drop_geometry(vv))
ala_mq <- merge(ala_mq[is.na(vc_code)==F], nvis_codes, by='vc_code')
ala_mq[vc_code %in% c(2,3,4,5,11)]
unique(ala_mq[,.(vc_name,vc_code)])

tmp_sp_obs <- ala_mq[,species:=str_replace(taxon," ",".")]


vec_names <- names(out)[str_which(names(out),"Angophora")[1]:dim(out)[2]]


mout <- melt(id.vars=c("x","y"),
  out[,c("x","y",eval(vec_names)),with=F],
  variable.name = 'species')

# mout <- melt(id.vars=c("x","y"),
#      out[,99:210][,-c("predict")], 
#      variable.name = 'species')

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
   geom_point(data=tmp_sp_obs[species%in%vec_names[idx[i]:idx[i+1]]], 
     inherit.aes = F, 
     aes(x,y), 
     col='pink',
     shape=21,
     fill=NA,
     alpha=0.333,
     size=2)+
    scale_fill_stepsn(
      colors=c("black",viridis::turbo(100)),
      breaks=c(0,0.01,0.1,0.3,0.5,0.7,1))+
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
       filename=file.path(paste0('figures/euc_sdm_v0.3_',Sys.Date()),im_name), 
         width=16*2,
         height=9*2,
         units='cm',
         dpi=350)
}


# SCRATCH **********************************************************************
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
# 
# 
# 
# 
# ggplot(data=mout[species%in%vec_names[idx[i]:2]],
#          aes(x,y,fill=value))+
#     geom_sf(data=oz_poly, 
#             inherit.aes = F, 
#             fill='grey20',
#             color='black')+
#     geom_raster()+
#     geom_sf(data=oz_poly, 
#             inherit.aes = F, 
#             fill=NA,
#             color='grey80')+
#    geom_point(data=tmp_sp_obs[species%in%vec_names[idx[i]:2]], 
#      inherit.aes = F, 
#      aes(x,y), 
#      col='red',
#      fill='white',
#      alpha=0.25,
#      size=2)+
#     scale_fill_stepsn(
#       colors=c("black",viridis::turbo(100)),
#       breaks=c(0,0.01,0.1,0.3,0.5,0.7,1),
#     )+
#     # scale_fill_viridis_b(option='H',
#     #                      # n.breaks=6,
#     #                      breaks=c(-Inf,0.01,0.1,0.3,0.5,0.7,Inf),
#     #                      limits=c(0,1))+
#     coord_sf(ylim=c(-38.75,-28), 
#              xlim=c(145.5,153.5))+
#     scale_x_continuous(breaks=c(146,149,152))+
#     labs(x=NULL, 
#          y=NULL, 
#          fill='prob.' 
#          # title=str_replace(eval(vec_names[i]),"\\."," ")
#     )+
#     facet_wrap(~species,nrow=1)+
#     # guides(fill = guide_colorbar(ncol=1))+
#     theme(
#       #  legend.position = c(1,0),
#       # legend.justification = c(1,0),
#       strip.background = element_rect(fill='grey90'),
#       panel.background = element_rect(fill='lightblue'), 
#       panel.grid = element_blank())
