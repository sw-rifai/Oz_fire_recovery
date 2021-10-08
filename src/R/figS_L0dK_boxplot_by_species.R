pacman::p_load(tidyverse, stars, data.table, lubridate, arrow, patchwork)

# Load vulnerable species 
vs <- fread("../data_general/ALA/Gallagher_cleaned/threatened_status.csv")
fn_simp <- function(x){
  x <- str_replace(x, " x ", " ")
  x <- str_remove(x, " sp\\.")
  x <- str_remove(x, " subsp\\.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}
vs <- vs[,species := fn_simp(Taxon),by=seq_len(nrow(vs))]
# tmp_sp_obs <- ala_mq[,species:=str_replace(species," ",".")]


# load fits ----
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-16 09:31:56.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>=0][r<0.024][r2>0.333][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
 ]
# estimate TTR from the logistic function
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]


# predicted dominant species ---- 
# sout <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-06-25 09:34:41.tiff")
# sout_rat <- fread("../data_general/proc_data_Oz_fire_recovery/predicted_species-distribution-ala-mq_top-40-species_LUT.csv")
# sout_rat

out <- read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-14 16:45:26.parquet")
sout <- set_names(st_as_stars(out[,.(x,y,predict)], dims = c("x","y"),crs=4326),c("species"))
st_crs(sout) <- st_crs(4326)
vv <- st_extract(sout,
  st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))

# Plot boxplot TTR by most probably species ------------------------------------
rank_nobs <- fits[is.na(species)==F][,.(nobs=.N),by=species][,nob_rank:=frank(-nobs)]
rank_ldk <- fits[is.na(species)==F][,.(val = median(ldk,na.rm=T)), by=species][order(-val)][,fct_order:=.I][]


rank_nobs$species %in% vs[`listed EPBC/NSW/VIC/ACT/QLD`==1]$species %>% table
vs[`listed EPBC/NSW/VIC/ACT/QLD`==1]
vs_listed <- vs[,`:=`(listed=`listed EPBC/NSW/VIC/ACT/QLD`)][,.(species,listed)]
vs_listed <- vs_listed[species%in%rank_ldk$species] %>% 
  group_by(species) %>% 
  summarize(threat = max(listed,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(color=if_else(threat==1, '#cf0000', 'black')) %>% 
  as.data.table()
vs_listed <- merge(vs_listed, rank_ldk, by='species')[order(fct_order)]

fits <- merge(vs_listed, fits, by='species',all.x=F,all.y=T,allow.cartesian = T)
fits <- fits[is.na(species)==F]



# Species level cor
merge(clim_species, fits[,.(val=median(r)),by=species],by='species') %>% 
  select(-species) %>% 
  cor
# Pixel level cor
fits %>% select(mapet, matmax,matmin,mavpd15,mappet,elevation,ldk) %>% cor

str_replace("Eucalyptus cunninghamii",pattern = "Eucalyptus",replacement = "E.")

rank_ldk <- rank_ldk %>% 
  .[,`:=`(species = str_replace(species,pattern = "Corymbia",replacement = "C."))] %>% 
  .[,`:=`(species = str_replace(species,pattern = "Eucalyptus",replacement = "E."))] %>% 
  .[,`:=`(species = str_replace(species,pattern = "Angophora",replacement = "A."))] 

p_right <- fits[is.na(species)==F][ldk>0] %>% 
  merge(., clim_species, by='species') %>% 
  .[,`:=`(species = str_replace(species,pattern = "Corymbia",replacement = "C."))] %>% 
  .[,`:=`(species = str_replace(species,pattern = "Eucalyptus",replacement = "E."))] %>% 
  .[,`:=`(species = str_replace(species,pattern = "Angophora",replacement = "A."))] %>%
  .[species %in% rank_ldk$species[1:55]] %>% 
  ggplot(data=., aes(y=species, 
                     x=ldk, 
                     fill=mapet_sp_u))+
  geom_boxplot(outlier.colour = NA,color='grey60')+
  geom_vline(aes(xintercept=quantile(fits[ldk>0]$ldk,0.1,na.rm=T)), 
             col="black",lty=2,lwd=0.5)+
  geom_vline(aes(xintercept=median(fits[ldk>0]$ldk,na.rm=T)), 
             col="black",lwd=1)+
  geom_vline(aes(xintercept=quantile(fits[ldk>0]$ldk,0.9,na.rm=T)), 
             col="black",lty=2,lwd=0.5)+
  scale_y_discrete(#breaks=rank_ldk$species[1:55], 
                   limits=rev(rank_ldk$species[1:55]))+
  # scale_y_discrete(breaks=rank_ldk$species[1:55], 
  #                  limits=rank_ldk$species[1:55])+
  scale_fill_viridis_c(option='A',
                       direction=1,
                       begin = 0, 
                       limits=c(900,1600),
                       oob=scales::squish)+
  scale_x_continuous(breaks=c(0,500,1000,1500,2500))+
  coord_cartesian(xlim=c(0,2550),
                  expand=F)+
  labs(x="Time to Recover (days)",
       y=NULL,
       fill="PET\n(mm)")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'none',
        axis.text.y = element_text(color = rev(vs_listed$color[1:55]), 
          size=13),
            axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=13),
        # legend.position = c(0.99,0.1), 
        # legend.justification = c(0.99,0.1),
        legend.background = element_rect(fill="#FFFFFFBB")); p_right

p_left <- fits[is.na(species)==F][ldk>0] %>% 
  merge(., clim_species, by='species') %>% 
  .[,`:=`(species = str_replace(species,pattern = "Corymbia",replacement = "C."))] %>% 
  .[,`:=`(species = str_replace(species,pattern = "Eucalyptus",replacement = "E."))] %>% 
  .[,`:=`(species = str_replace(species,pattern = "Angophora",replacement = "A."))] %>%
  .[species %in% rank_ldk$species[56:109]] %>%
  ggplot(data=., aes(y=species, 
                     x=ldk, 
                     fill=mapet_sp_u))+
  geom_boxplot(outlier.colour = NA,color='grey60')+
  geom_vline(aes(xintercept=quantile(fits[ldk>0]$ldk,0.1,na.rm=T)), 
             col="black",lty=2,lwd=0.5)+
  geom_vline(aes(xintercept=median(fits[ldk>0]$ldk,na.rm=T)), 
             col="black",lwd=1)+
  geom_vline(aes(xintercept=quantile(fits[ldk>0]$ldk,0.9,na.rm=T)), 
             col="black",lty=2,lwd=0.5)+
  scale_y_discrete(#breaks=rank_ldk$species[1:55], 
    limits=rev(rank_ldk$species[56:106]))+
  # scale_y_discrete(breaks=rank_ldk$species[1:55], 
  #                  limits=rank_ldk$species[1:55])+
  scale_fill_viridis_c(option='A',
                       direction=1,
                       begin = 0, 
                       limits=c(900,1600),
                       oob=scales::squish)+
  scale_x_continuous(breaks=c(0,500,1000,1500,2000))+
  coord_cartesian(xlim=c(0,2650),
                  expand=F)+
  labs(x="Time to Recover (days)",
       y=NULL,
       fill="PET\n(mm yr¯¹)")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
            axis.text.y = element_text(color = rev(vs_listed$color[56:106]), 
              size=13),
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=13),
        legend.position = 'none',#c(0.99,0.01), 
        legend.justification = c(0.99,0.01),
        legend.key.height = unit(1,'cm'),
        legend.background = element_rect(fill="#FFFFFFBB")); p_left

(p_left+patchwork::inset_element(p_legend,left = 0.6,bottom = 0,right = 1,top = 0.333))+
  p_right+
  plot_layout(guides='keep')

scale_factor <- 3
ggsave(
       filename="figures/boxplot_L0dK_by_most-probable-species.png",
       width=8.1*scale_factor,
       height=9.1*scale_factor,
       units='cm',
       dpi=350)


