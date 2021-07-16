pacman::p_load(tidyverse, 
  stars,
 data.table)

dir.create(paste0('figures/euc_sdm_v0.3_',Sys.Date()))

# ML predictions ---------------------
out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-14 16:45:26.parquet")

# Oz poly -------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 1000) %>% 
  select(NAME_1)
oz_poly <- oz_poly %>% filter(NAME_1 %in% c("Queensland","New South Wales","Australian Capital Territory","Victoria"))

# nvis ---------------------
nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif", 
                          proxy=F)
names(nvis) <- 'vc_code'
nvis_codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("vc_code","veg_class_descrip")) %>% 
  mutate(vc_name = as.factor(veg_class_descrip)) %>% 
  select(-veg_class_descrip)

nvis_mask <- nvis %>% 
  # st_crop(., st_bbox(c(xmin=145.5 , ymin= -38.75, xmax=153.5, ymax=-28),
  # crs=4326)) %>% 
  as_tibble() %>% 
  mutate(vc_code = ifelse(vc_code%in%c(2,3,4,5,11),NA,0))

# ALA MQ cleaned --------------
ala_mq <- fread("../data_general/ALA/Gallagher_cleaned/euc_occurrence_clean_v2.csv") %>% 
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

## Simplify species names -----
ala_mq$taxon %>% unique

fn_simp <- function(x){
  x <- str_replace(x, " x ", " ")
  x <- str_remove(x, " sp\\.")
  x <- str_remove(x, " subsp\\.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}

str_replace("Eucalyptus x macmahonii", " x ", " ")
fn_simp("Eucalyptus x macmahonii")


ala_mq <- ala_mq[,species := fn_simp(taxon),by=seq_len(nrow(ala_mq))]
tmp_sp_obs <- ala_mq[,species:=str_replace(species," ",".")]


vec_names <- names(out)[str_which(names(out),"Angophora")[1]:dim(out)[2]]


mout <- melt(id.vars=c("x","y"),
  out[,c("x","y",eval(vec_names)),with=F],
  variable.name = 'species')

# mout <- melt(id.vars=c("x","y"),
#      out[,99:210][,-c("predict")], 
#      variable.name = 'species')


vec_split <- split(vec_names, cut(seq_along(vec_names), 13, labels = FALSE))


# idx <- seq(1,length(vec_names),by=9)
i <- 1
for(i in 1:length(vec_split)){
  print(i)
 p_out <-  ggplot(data=mout[species%in%vec_split[[i]]],
         aes(x,y,fill=value))+
    geom_sf(data=oz_poly, 
            inherit.aes = F, 
            fill='grey20',
            color='black')+
    geom_raster()+
    geom_raster(
     inherit.aes = F,
      data=nvis_mask %>% filter(is.na(vc_code)==F), 
      aes(x,y))+
    geom_sf(data=oz_poly, 
            inherit.aes = F, 
            fill=NA,
            color='grey80')+
    geom_sf(data=oz_poly, 
            inherit.aes = F, 
            fill=NA,
            color='grey80')+
   geom_point(data=tmp_sp_obs[species%in%vec_split[[i]]], 
     inherit.aes = F, 
     aes(x,y), 
     col='pink',
     shape=21,
     fill=NA,
     alpha=0.1,
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
    facet_wrap(~species,nrow=3)+
    # guides(fill = guide_colorbar(ncol=1))+
    theme(
       legend.position = 'right',
      # legend.justification = c(1,0),
      strip.background = element_rect(fill='grey90'),
      panel.background = element_rect(fill='lightblue'), 
      panel.grid = element_blank())
      
  im_name <- paste0(c(str_replace(c("sdm_eucs_xgboost_",vec_split[[i]][1],"_",last(vec_split[[i]])),pattern="\\.",replacement = "_"),".png"),collapse = "")
  ggsave(p_out, 
       filename=file.path(paste0('figures/euc_sdm_v0.3_',Sys.Date()),im_name), 
         width=8*2,
         height=11*2,
         units='cm',
         dpi=350)
}


# SCRATCH **********************************************************************
# out[,c('x','y',vec_names[i]), with=F] %>%
#   ggplot(data=.,aes(x,y,fill= get(vec_names[i]) ))+
#   geom_sf(data=oz_poly,
#           inherit.aes = F,
#           fill='grey20',
#           color='black')+
#   geom_raster()+
#   geom_sf(data=oz_poly,
#           inherit.aes = F,
#           fill=NA,
#           color='grey80')+
#   scale_fill_viridis_b(option='F',
#                        n.breaks=6,
#                        limits=c(0,1))+
#   coord_sf(ylim=c(-38.75,-28),
#            xlim=c(145.5,153.5))+
#   labs(x=NULL,
#        y=NULL,
#        fill='prob.',
#        title=str_replace(eval(vec_names[i]),"\\."," "))+
#   # guides(fill = guide_colorbar(ncol=1))+
#   theme(legend.position = c(1,0),
#         legend.justification = c(1,0),
#         panel.background = element_rect(fill='lightblue'),
#         panel.grid = element_blank())
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
