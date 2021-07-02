library(tidyverse)
library(stars); 
# library(dtplyr)
library(data.table) 
library(lubridate)
library(arrow)
# library(furrr)
library(patchwork)


oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 1000) %>% 
  select(NAME_1)
oz_poly <- oz_poly %>% filter(NAME_1 %in% c("Queensland","New South Wales","Australian Capital Territory","Victoria"))

# load fits ----
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-1mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-06-19 17:33:16.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>=0][r<0.024][r2>0.333][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
 ]
# estimate TTR from the logistic function
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]


# predicted dominant species ---- 
# sout <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-06-25 09:34:41.tiff")
sout_rat <- fread("../data_general/proc_data_Oz_fire_recovery/predicted_species-distribution-ala-mq_top-40-species_LUT.csv")
sout_rat

out <- read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-06-25 12:26:19.parquet")
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



# dominant species ---- 
# dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
# dom <- dom[,.(x,y,id,dom_sp)]
# dom_n <- dom[,.(nobs = .N), by='dom_sp'][,rank := frank(-nobs)][order(rank)][rank<=40]

# dominant species clim habitats
# hab <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_species-median-clim-topo.parquet")
# hab[,dom_sp:=species]
# 
# dom_n <- merge(dom_n,hab,by='dom_sp')

# d13 ---------
d13 <- fread("../data_general/AusTraits/wcornwell-leaf13C.csv")


# AusTraits ----
at <- fread("../data_general/AusTraits/austraits-2.1.0/data/traits.csv")

# functions for species cleanup ----
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

simp_name <- function(x){
  x <- str_remove(x, " sp\\.")
  x <- str_remove(x, " subsp\\.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}

# simplify taxon tames from AusTraits
at_s <- at[,species:=simp_name(taxon_name), by=seq_len(nrow(at))]
at_m <- at_s[trait_name %in% c(
  "leaf_area","leaf_delta13C","wood_density",
  "specific_leaf_area","photosynthetic_rate_per_area_maximum",
  "photosynthetic_rate_per_area_saturated"
)][,.(val = median(as.numeric(value))), by=.(species, trait_name)]

# Get fire traits from AusTraits using Nicolle 2006
nic <- at[dataset_id=="Nicolle_2006"][trait_name %in% c('fire_response',"regen_strategy")]
nic <- nic[,species := simp_name(taxon_name),by=seq_len(nrow(nic))][
  ,.(value = get_mode(value)), by=.(species,trait_name)
]
missing_sp_fr <- dom_n$dom_sp[!dom_n$dom_sp%in%unique(nic[trait_name=='fire_response']$species)]
missing_sp_rs <- dom_n$dom_sp[!dom_n$dom_sp%in%unique(nic[trait_name=='regen_strategy']$species)]

tmp_fr <- at[trait_name == "fire_response"][value %in% c("fire_killed","resprouts")][taxon_name%in%missing_sp_fr][
  ,species := taxon_name][
  ,.(value = get_mode(value)), by=.(species)
]

tmp_rs <- at[trait_name == "regen_strategy"][taxon_name%in%missing_sp_rs][
  ,species := taxon_name][
    ,.(value = get_mode(value)), by=.(species)
  ]


d_fr <- rbindlist(list(nic[trait_name=='fire_response'][,.(species,value)], tmp_fr))
d_rs <- rbindlist(list(nic[trait_name=='regen_strategy'][,.(species,value)], tmp_rs))


dom_n <- merge(dom_n[,species := dom_sp], 
      d_fr %>% rename(fire_response=value) %>% as.data.table(), 
      by='species')

dom_n <- merge(dom_n[,species := dom_sp], 
               d_rs %>% rename(regen_strategy=value) %>% as.data.table(),
               by='species') %>% 
  as.data.table()


dat <- merge(fits,dom,by=c("id","x","y"))
dat <- merge(dat, dom_n, by=c("dom_sp"))


# Plotting --------------------------------------------------------------------
pals::pal.bands(pals::parula)
colorspace::hcl_palettes(plot=T)
fits[r2>0.333][isConv==T] %>% 
  ggplot(data=.,aes((L0/K),r,color=cut_width(mappet,width = 0.5)))+
  geom_density_2d_filled(color=NA,
                         aes(contour_var="ndensity"), 
                         bins=15)+
  geom_smooth(se=F,
              formula=y~s(x,bs='cs',k=3))+
  scale_color_manual(values =rev(pals::parula(8)))+
  colorspace::scale_fill_discrete_sequential(palette='Grays')+
  labs(x=expression(paste(L[0]/K~~(m**2/m**2))), 
       y="r (m²/day)", 
       color="P:PET ratio")+
  guides(fill = guide_none())+
  coord_cartesian(expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave(filename='figures/gam-smooth_r-l0K_by-PPET-range.png', 
       width=15, 
       height=10,
       units='cm', 
       dpi=350)

# Plot boxplot TTR by most probably species ------------------------------------
rank_nobs <- fits[is.na(species)==F][,.(nobs=.N),by=species][,nob_rank:=frank(-nobs)]
rank_ttr <- fits[is.na(species)==F][,.(val = median(pred_ttr,na.rm=T)), by=species][order(-val)][,fct_order:=.I][]
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
  labs(y=expression(paste("PET (",mm~yr**-1,")" )),
    x=NULL)+
  theme_minimal()+
  theme(panel.grid = element_blank()))

# Species level cor
merge(clim_species, fits[,.(val=median(r)),by=species],by='species') %>% 
  select(-species) %>% 
  cor
# Pixel level cor
fits %>% select(mapet, matmax,matmin,mavpd15,mappet,elevation,pred_ttr) %>% cor
p_right <- fits[is.na(species)==F][pred_ttr>0] %>% 
  merge(., clim_species, by='species') %>% 
  .[species %in% rank_ttr$species[1:55]] %>% 
  ggplot(data=., aes(y=species, 
                     x=pred_ttr, 
                     fill=mapet_sp_u))+
  geom_boxplot(outlier.colour = NA,color='grey60')+
  geom_vline(aes(xintercept=quantile(fits[pred_ttr>0]$pred_ttr,0.1,na.rm=T)), 
             col="#CF0000",lty=1,lwd=0.5)+
  geom_vline(aes(xintercept=median(fits[pred_ttr>0]$pred_ttr,na.rm=T)), 
             col="#CF0000",lwd=1)+
  geom_vline(aes(xintercept=quantile(fits[pred_ttr>0]$pred_ttr,0.9,na.rm=T)), 
             col="#CF0000",lty=1,lwd=0.5)+
  scale_y_discrete(#breaks=rank_ttr$species[1:55], 
                   limits=rev(rank_ttr$species[1:55]))+
  # scale_y_discrete(breaks=rank_ttr$species[1:55], 
  #                  limits=rank_ttr$species[1:55])+
  scale_fill_viridis_c(option='A',
                       direction=1,
                       begin = 0, 
                       limits=c(900,1600),
                       oob=scales::squish)+
  scale_x_continuous(breaks=c(0,500,1000,1500,2500))+
  coord_cartesian(xlim=c(0,2500),
                  expand=F)+
  labs(x="Time to Recover (days)",
       y=NULL,
       fill="PET\n(mm)")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'none',
        # legend.position = c(0.99,0.1), 
        # legend.justification = c(0.99,0.1),
        legend.background = element_rect(fill="#FFFFFFBB"))

p_left <- fits[is.na(species)==F][pred_ttr>0] %>% 
  merge(., clim_species, by='species') %>% 
  .[species %in% rank_ttr$species[56:109]] %>% 
  ggplot(data=., aes(y=species, 
                     x=pred_ttr, 
                     fill=mapet_sp_u))+
  geom_boxplot(outlier.colour = NA,color='grey60')+
  geom_vline(aes(xintercept=quantile(fits[pred_ttr>0]$pred_ttr,0.1,na.rm=T)), 
             col="#CF0000",lty=1,lwd=0.5)+
  geom_vline(aes(xintercept=median(fits[pred_ttr>0]$pred_ttr,na.rm=T)), 
             col="#CF0000",lwd=1)+
  geom_vline(aes(xintercept=quantile(fits[pred_ttr>0]$pred_ttr,0.9,na.rm=T)), 
             col="#CF0000",lty=1,lwd=0.5)+
  scale_y_discrete(#breaks=rank_ttr$species[1:55], 
    limits=rev(rank_ttr$species[56:109]))+
  # scale_y_discrete(breaks=rank_ttr$species[1:55], 
  #                  limits=rank_ttr$species[1:55])+
  scale_fill_viridis_c(option='A',
                       direction=1,
                       begin = 0, 
                       limits=c(900,1600),
                       oob=scales::squish)+
  scale_x_continuous(breaks=c(0,500,1000,1500,2000))+
  coord_cartesian(xlim=c(0,2500),
                  expand=F)+
  labs(x="Time to Recover (days)",
       y=NULL,
       fill="PET\n(mm yr¯¹)")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'none',#c(0.99,0.01), 
        legend.justification = c(0.99,0.01),
        legend.key.height = unit(1,'cm'),
        legend.background = element_rect(fill="#FFFFFFBB"))

(p_left+patchwork::inset_element(p_legend,left = 0.6,bottom = 0,right = 1,top = 0.333))+
  p_right+
  plot_layout(guides='keep')

scale_factor <- 3.5
ggsave(filename="figures/boxplot_pred-ttr-logfit_by_most-probable-species.png",
       width=8.1*scale_factor,
       height=9.1*scale_factor,
       units='cm',
       dpi=350)

fits$species %>% table
fits[str_detect(species,'camaldulensis')]
fits[is.na(species)==F] %>% 
  mutate(prop_remain = cut_interval((L0/K), 3)) %>% 
  as_tibble() %>% 
  ggplot(data=., aes(y=species, 
                     x=r,
                     # group=paste(species,prop_remain), 
                     # color=prop_remain,
                     fill=prop_remain))+
  geom_boxplot(outlier.colour = NA)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(0,0.025))+
  facet_wrap(~prop_remain)


# Plot 1: ----
library(colorspace)
colorspace::hcl_palettes(plot=T)
ggplot(data=out, aes(x,y,fill=predict))+
  geom_sf(data=oz_poly, inherit.aes = F)+
  geom_raster()+
  scale_fill_viridis_d(option='H')+
  # scale_fill_discrete_qualitative(palette='Dynamic')+
  coord_sf(ylim=c(-40,-28), 
           xlim=c(145.5,153.5))+
  guides(fill = guide_legend(ncol=2))+
  theme(legend.position = 'right')


p1 <- dat %>% ggplot(data=.,aes(fire_response,r))+
  geom_boxplot(outlier.colour = NA)+
    labs(x='Fire Response')+
    coord_cartesian(ylim=c(0,0.015),expand=c(0,0))

p2 <- dat %>% ggplot(data=.,aes(regen_strategy,r))+
  geom_boxplot(outlier.colour = NA)+
  labs(x='Regeneration Strategy')+
  coord_cartesian(ylim=c(0,0.015),expand=c(0,0))
p1+p2+plot_layout(widths = c(1,2))
ggsave(p1+p2+plot_layout(widths = c(1,2)), 
       filename = 'figures/boxplot_r_by_fire-response_regen-strategy.png',
       width=25,
       height=10,
       dpi=350,
       units='cm')
# END Plot **************************************************************************


# Plot : ----
library(forcats)
dat[species %in% dom_n$dom_sp][,.(n_fires = length(unique(date_fire1))),by=species][
  order(n_fires)
]

fits %>% 
  ggplot(data=., 
         aes(y=species, 
             x=(L0/K)/r))+
  geom_boxplot(outlier.colour = NA)+
  coord_cartesian(xlim=c(0,550))


fits[is.na(species)==F] %>% 
  ggplot(data=.,aes(x=r, 
                    y=species,
                    # y=fct_reorder(species, ttr5_lai, .fun=median), 
                    fill=malai))+
  geom_boxplot(outlier.colour = NA)+
  labs(y=NULL,
       x='Time to recover LAI (days)',
       fill=expression(paste(MAP~bgroup("(",frac(mm,yr),")"))))+
  geom_rect(aes(xmin=0, xmax=365, 
                ymin=0,ymax=40),
            lty=1,fill='grey')+
  coord_cartesian(xlim=c(0,3200), 
                  expand = c(0,0))+
  scale_fill_viridis_c(option='E',direction=-1)+
  theme_linedraw()


dat[species %in% dom_n$dom_sp] %>% 
  merge(., dom_n[,.(species, map, mavpd15)], by='species') %>% 
  ggplot(data=.,aes(x=ttr5_lai, 
                    y=fct_reorder(species, ttr5_lai, .fun=median), 
                    fill=map))+
  geom_boxplot(outlier.colour = NA)+
  labs(y=NULL,
       x='Time to recover LAI (days)',
       fill=expression(paste(MAP~bgroup("(",frac(mm,yr),")"))))+
  geom_rect(aes(xmin=0, xmax=365, 
                ymin=0,ymax=40),
            lty=1,fill='grey')+
  coord_cartesian(xlim=c(0,3200), 
                  expand = c(0,0))+
  scale_fill_viridis_c(option='E',direction=-1)+
  theme_linedraw()
ggsave(filename = "figures/boxplot_ala-mq_dom_nsp_ttr5-lai.png", 
       width=20,
       height=30,
       dpi=350,
       units='cm')

# END Plot **************************************************************************

# Plot : ----
p_l <- merge(dat, at_m[trait_name=='wood_density'], by='species') %>%
  ggplot(data=.,aes(x=val, 
                    y=r))+
  geom_point(alpha=0.1)+
  geom_smooth(method='lm')+
  labs(x="Wood Density",
       y=expression(paste(bolditalic(r))))+
  theme_linedraw()
p_r <- merge(dat, at_m[trait_name=='specific_leaf_area'], by='species') %>%
  ggplot(data=.,aes(x=val, 
                    y=r))+
  geom_point(alpha=0.1)+
  geom_smooth(method='lm')+
  labs(x="Specific Leaf Area",
       y=expression(paste(bolditalic(r))))+
  theme_linedraw()
p_l|p_r
ggsave(p_l|p_r,
       filename = "figures/scatterplot_r_wood-density_SLA.png", 
       width=20,
       height=10,
       dpi=350,
       units='cm')

# END Plot **************************************************************************

# Plot : ----
at_m[trait_name=='leaf_delta13C']
at_m[trait_name=='ci_over_ca']
at[trait_name=='ci_over_ca']$species %>% table

merge(dat, at_m[trait_name=='leaf_delta13C'], by='species') %>%
  ggplot(data=.,aes(x=val, 
                    y=r))+
  geom_point(alpha=0.1)+
  geom_smooth(method='lm')+
  labs(x="leaf_delta13C",
       y=expression(paste(bolditalic(r))))+
  theme_linedraw()

# END Plot **************************************************************************

# Plot : ----

# END Plot **************************************************************************

# Plot : ----

# END Plot **************************************************************************



# fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]

clim_species <- fits[,.(map_sp_u = mean(map), 
                        matmax_sp_u = mean(matmax), 
                        matmin_sp_u = mean(matmin), 
                        mavpd15_sp_u = mean(mavpd15), 
                        mappet_sp_u = mean(mappet),
                        matrange_sp_u = mean(range(matmax,matmin))), by=species]


fits[is.na(species)==F][isConv==T] %>% 
  merge(., clim_species, by='species') %>% 
  merge(., fits[,.(nobs=.N),by='species'][,rank:=frank(-nobs)][], 
        by='species') %>% 
  ggplot(data=.,aes((L0/K),r,group=species,color=mappet_sp_u))+
  geom_smooth(se=F, 
              formula=y~s(x,bs='cs',k=3))+
  scale_color_viridis_c(option='H',direction = -1,end=0.9)


fits[is.na(species)==F][isConv==T] %>% 
  ggplot(data=.,aes((L0/K),r,color=cut_width(map,width = 200)))+
  geom_smooth(se=T, 
              formula=y~s(x,bs='cs',k=3))+
  scale_color_viridis_d(option='H',direction = -1,end=0.9)


fits[is.na(species)==F][isConv==T] %>% 
  ggplot(data=.,aes((L0/K),r,color=cut_number(map,n=5)))+
  geom_smooth(se=T, 
              formula=y~s(x,bs='cs',k=3))+
  scale_color_viridis_d(option='H',direction = -1,end=0.9)


