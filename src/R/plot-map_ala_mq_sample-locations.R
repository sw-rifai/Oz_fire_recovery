pacman::p_load(tidyverse, 
  stars,
 data.table)

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

vec_nobs_80 <- ala_mq[,.(nobs=.N),by='species'][,rank:=frank(-nobs,ties.method = 'first')][order(rank)][nobs>=80]$species
fdat2 <- ala_mq[species%in%vec_nobs_80]
fdat2[,species:=as.factor(species)]


ggplot(data=fdat2,aes(x,y))+
  # geom_stars(data=nvis,aes(fill=factor(vc_code)))+
  # scale_fill_viridis_d(option='D',direction = -1,alpha=1)+
  geom_sf(data=oz_poly,
    inherit.aes = F,
    color='grey70',
    aes(fill=NAME_1))+
  scale_fill_viridis_d(option='H',begin=0.1, end=0.7,direction = -1)+
  coord_sf(ylim=c(-38.75,-28), 
             xlim=c(145.5,153.5))+
  scale_x_continuous(breaks=c(146,149,152))+
  geom_point(col='black',pch=20,size=0.05,alpha=0.25)+
  labs(x=NULL,y=NULL,fill=NULL)+
  theme_linedraw()+
  theme(legend.position = c(1,0),
    legend.justification = c(1,0),
    panel.grid = element_blank(),
    panel.background = element_rect(fill='lightblue'))

ggsave(filename = 'figures/plot-map_ala-mq_sample-locations.png',
  device=grDevices::png,
  width=15,
  height=20,
  units='cm',
  dpi=350)
