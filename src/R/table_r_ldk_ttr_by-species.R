pacman::p_load(tidyverse, stars, data.table, lubridate, arrow, kableExtra)

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
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-22 09:37:07.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>=0][r<0.024][r2>0.333][month%in%c(9,10,11,12,1,2)][
  ,TTR_LF:=(L0/K)
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

# Attach climate
clim_ma <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y",
                                           "map","mapet","mappet",
                                           "matmax","matmin","mavpd15")) %>% 
  distinct() %>% 
  as.data.table()
clim_ma <- st_as_stars(clim_ma, dims = c('x','y'))
st_crs(clim_ma) <- st_crs(4326)
vv <- st_extract(clim_ma,
                 st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))

# Attach elevation
s_dem <- stars::read_stars("../data_general/Oz_misc_data/SRTM_elevation_500m_EastOz_.tif")
names(s_dem) <- 'elevation'
# locs <- st_as_stars(fits, dims = c('x','y'))
# st_crs(locs) <- st_crs(4326)
vv <- st_extract(s_dem,
                 st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))


# Plot boxplot TTR by most probably species ------------------------------------
rank_nobs <- fits[is.na(species)==F][,.(nobs=.N),by=species][,nob_rank:=frank(-nobs)]
rank_ttr <- fits[is.na(species)==F][,.(val = median(pred_ttr,na.rm=T)), by=species][order(-val)][,fct_order:=.I][]


rank_nobs$species %in% vs[`listed EPBC/NSW/VIC/ACT/QLD`==1]$species %>% table
vs[`listed EPBC/NSW/VIC/ACT/QLD`==1]
vs_listed <- vs[,`:=`(listed=`listed EPBC/NSW/VIC/ACT/QLD`)][,.(species,listed)]
vs_listed <- vs_listed[species%in%rank_ttr$species] %>% 
  group_by(species) %>% 
  summarize(threat = max(listed,na.rm=T)) %>% 
  ungroup() %>% 
  mutate(color=if_else(threat==1, '#cf0000', 'black')) %>% 
  as.data.table()
vs_listed <- merge(vs_listed, rank_ttr, by='species')[order(fct_order)]

fits <- merge(vs_listed, fits, by='species',all.x=F,all.y=T,allow.cartesian = T)
fits <- fits[is.na(species)==F]


clim_species <- fits[,.(map_sp_u = mean(map), 
                        matmax_sp_u = mean(matmax), 
                        matmin_sp_u = mean(matmin), 
                        mavpd15_sp_u = mean(mavpd15), 
                        mappet_sp_u = mean(mappet),
                        mapet_sp_u = mean(mapet),
                        matrange_sp_u = mean(range(matmax,matmin)), 
                        elevation_u = mean(elevation,na.rm=T)), by=species]
(p_legend <- fits %>% 
  ggplot(data=.,aes(y=mapet, 
    fill=cut_interval(mapet, 100)))+
  geom_histogram(show.legend = F)+
  scale_fill_viridis_d(option='A')+
  coord_cartesian(ylim=c(900,1600), 
    expand=F)+
  scale_x_continuous(breaks=c(0,10000))+
  labs(y=expression(paste("Potential Evapotranspiration  (",mm~yr**-1,")" )),
    x=NULL)+
  theme_minimal()+
  theme(panel.grid = element_blank()
  ,axis.text.x = element_blank()))

# Species level cor
merge(clim_species, fits[,.(val=median(r)),by=species],by='species') %>% 
  select(-species) %>% 
  cor
# Pixel level cor
fits %>% select(mapet, matmax,matmin,mavpd15,mappet,elevation,pred_ttr) %>% cor


library(kableExtra)
fits[is.na(species)==F][pred_ttr>0] %>% 
  .[,.(nobs = .N, 
     r_lo = quantile(r, 0.05), 
     r_med = quantile(r,0.5),
     r_hi = quantile(r,0.95), 
     K_lo = quantile(K, 0.05), 
     K_med = quantile(K,0.5),
     K_hi = quantile(K,0.95), 
     ldk_lo = quantile(ldk, 0.05), 
     ldk_med = quantile(ldk,0.5),
     ldk_hi = quantile(ldk,0.95), 
     pred_ttr_lo = quantile(pred_ttr, 0.05), 
     pred_ttr_med = quantile(pred_ttr,0.5),
     pred_ttr_hi = quantile(pred_ttr,0.95), 
     ttr5_lai_lo = quantile(ttr5_lai, 0.05), 
     ttr5_lai_med = quantile(ttr5_lai,0.5),
     ttr5_lai_hi = quantile(ttr5_lai,0.95)
    ),by=species] %>% 
  .[,species:=as.character(species)] %>% 
  .[order(species)] %>% 
  .[,r := paste0(format(r_med,digits=2)," (",format(r_lo,digits=2),"-",format(r_hi,digits=2),")" )] %>% 
  .[,K := paste0(format(K_med,digits=2)," (",format(K_lo,digits=2),"-",format(K_hi,digits=2),")" )] %>% 
  .[,ldk := paste0(format(ldk_med,digits=2)," (",format(ldk_lo,digits=2),"-",format(ldk_hi,digits=2),")" )] %>% 
  .[,TTR_LF := paste0(format(pred_ttr_med,digits=2)," (",format(pred_ttr_lo,digits=2),"-",format(pred_ttr_hi,digits=2),")" )] %>% 
  .[,TTR_MW := paste0(format(ttr5_lai_med,digits=2)," (",format(ttr5_lai_lo,digits=2),"-",format(ttr5_lai_hi,digits=2),")" )] %>% 
  select(species,N=nobs, r,K,`L0/K`= ldk, TTR_LF, TTR_MW) %>% 
  knitr::kable(., 'html') %>% 
  kable_styling("striped") %>%
  kableExtra::save_kable(., file = "doc/table_species-level_LF_TTR.pdf") 




