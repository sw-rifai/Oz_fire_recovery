pacman::p_load(tidyverse, stars, data.table, lubridate, arrow, mgcv, broom, patchwork)

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
fits[,fire_year := year(date_fire1-months(3))]

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



d1 <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=30 & ldk_range >= 0.5][,hi_lo := ifelse(L0/K<=0.5,'lo','hi')] %>% 
  filter(hi_lo == 'lo') %>% 
  nest(data=-species) %>% 
  mutate(
    fit = map(data, ~lm(r~I(L0/K), data=.x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)
d2 <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.7][,hi_lo := ifelse(L0/K<=0.5,'lo','hi')] %>% 
  filter(hi_lo == 'hi') %>% 
  nest(data=-species) %>% 
  mutate(
    fit = map(data, ~lm(r~I(L0/K), data=.x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)

vec_inc <- d2 %>% filter(term=='I(L0/K)') %>% 
  filter(estimate > 0) %>% 
  filter(p.value < 0.05) %>% 
  pull(species)

vec_species <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.7] %>% pull(species) %>% unique
vec_not_inc <- vec_species[!vec_species %in% vec_inc]

# plot species that increase r with less burn severity 
p_inc <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.5][,hi_lo := ifelse(L0/K<=0.5,'lo','hi')] %>% 
  .[species %in% vec_inc] %>% 
  ggplot(data=.,aes(L0/K, r))+
  ggpointdensity::geom_pointdensity(adjust=0.25)+
  scico::scale_color_scico(palette='grayC',begin=0.25)+
  # scale_color_viridis_c(option='F',end=0.9)+
  # geom_point(alpha=0.25,col='grey70')+
  geom_smooth(method='gam',
    col='#0569FF', 
    lwd=1.5,
    formula=y~s(x,bs='cs',k=5), 
    method.args=list(family=Gamma(link='log')))+
  geom_rug()+
  labs(x=expression(paste("Post Fire Remaining Leaf Fraction: ", frac(L[0],K)~"   (m²/m²)")), 
       y=expression(paste("Leaf Growth Rate: ", bolditalic(r)~"(m²/day)")))+
  facet_wrap(~species, scales = 'free_y',ncol=3)+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
    legend.position = 'none', 
    strip.background = element_blank(), 
    strip.text = element_text(color='black',face = 'bold'))
scale_fraction <- 0.75
ggsave(p_inc, filename = 'figures/r_increase-with-L0dK.png', 
  height=60*scale_fraction, 
  width=45*scale_fraction, 
  units='cm', dpi=350)

# plot species that don't increase r with less burn severity 
p_dec <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.5] %>% 
  .[species %in% vec_not_inc] %>% 
  ggplot(data=.,aes(L0/K, r))+
  ggpointdensity::geom_pointdensity(adjust=0.25)+
  scico::scale_color_scico(palette='grayC',begin=0.25)+
  # geom_point(alpha=0.25,col='grey70')+
  geom_smooth(method='gam',
    color='#cf0000', 
    lwd=1.5,
    formula=y~s(x,bs='cs',k=5), 
    method.args=list(family=Gamma(link='log')))+
  geom_rug()+
  labs(x=expression(paste("Remaining Leaf Fraction", frac(L[0],K)~"(m²/m²)")), 
       y=expression(paste("Leaf Growth Rate: ", bolditalic(r)~"(m²/day)")))+
  facet_wrap(~species, scales = 'free_y',ncol=4)+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
    legend.position = 'none', 
    strip.background = element_blank(), 
    strip.text = element_text(color='black'))
scale_fraction <- 0.75
ggsave(p_dec, filename = 'figures/r_no-increase-with-L0dK.png', 
  height=60*scale_fraction, 
  width=45*scale_fraction, 
  units='cm', dpi=350)


# 
# # Plot boxplot TTR by most probably species ------------------------------------
# rank_nobs <- fits[is.na(species)==F][,.(nobs=.N),by=species][,nob_rank:=frank(-nobs)]
# rank_ttr <- fits[is.na(species)==F][,.(val = median(pred_ttr,na.rm=T)), by=species][order(-val)][,fct_order:=.I][]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# rank_species <- fits[r2>0.2][,.(nobs=.N, 
#   ldk_range = diff(range(ldk)), 
#   n_fy = length(unique(fire_year))),by=species][,rank:=frank(-nobs,ties.method = 'first')][]
# merge(fits,rank_species,by='species')[r2>0.2][rank <= 12] %>% 
#   ggplot(data=.,aes(ldk,r))+
#   geom_smooth(se=F)+
#   facet_wrap(~species)
# 
# 
# b1 <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.7 & n_fy>=3] %>% 
#   bam(r~s(I(L0/K), species, bs='fs',k=4),
#       # s(pre_fire_slai_anom_12mo, bs='cs')+
#       # s(species, bs='re'), 
#     family=Gamma(link='log'),
#     data=., 
#     select=T, discrete=T)
# plot(b1, pages=5)
# summary(b1)
# v1 <- mgcViz::getViz(b1)
# plot(sm(v1, 1)) + l_fitLine(alpha = 0.6) + labs(title = "Smooth factor interactions")
# plot(v1)+l_fitLine()
# 
# 
# p_test <- merge(fits,rank_species,by='species')[r2>0.2][is.na(species)==F][nobs>=100 & ldk_range >= 0.5 & n_fy>=3] %>% 
#   ggplot(data=.,aes(ldk,r))+
#   geom_smooth(se=F)+
#   geom_smooth(aes(min_slai/malai,r), 
#     se=F, col=c("#cf0000"))+
#   geom_rug()+
#   facet_wrap(~species)
# 
# p_test <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.7 & n_fy>=3] %>% 
#   ggplot(data=.,aes(L0/K, r))+
#   geom_point()+
#   geom_smooth(
#     method='gam',
#     formula=y~s(x,bs='cs',k=5),
#     se=F, col=c("#cf0000"))+
#   geom_rug()+
#   theme_linedraw()+
#   facet_wrap(~species, scales = 'free_y')
# ggsave(p_test, filename = 'test.png', 
#   width=60, height=40, units='cm', dpi=400)
# 
# p_test <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=30 & ldk_range >= 0.5][,hi_lo := ifelse(L0/K<=0.5,'lo','hi')] %>% 
#   ggplot(data=.,aes(L0/K, r,group=hi_lo,color=hi_lo))+
#   geom_point(alpha=0.25,col='grey70')+
#   geom_smooth(method=MASS::rlm)+
#   geom_rug()+
#   scale_color_viridis_d(option='H',end=0.9)+
#   facet_wrap(~species, scales = 'free_y')+
#   theme_linedraw()+
#   theme(panel.grid = element_blank())
# ggsave(p_test, filename = 'test.png', 
#   width=60, height=40, units='cm', dpi=400)
# 
# merge(fits,rank_species,by='species')[r2>0.2][rank <= 30] %>%
#   ggplot(aes(y=species,x=malai,group=species))+
#   geom_boxplot(outlier.colour = NA)
# 
# 
# 
# merge(fits,rank_species,by='species')[r2>0.2][is.na(species)==F][nobs>=100 & ldk_range >= 0.5 & n_fy>=3] %>% 
#   ggplot(data=.,aes(r, malai/ttr5_lai))+
#   geom_smooth()+
#   geom_abline(col='red')
# 
# 
# d1 <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=30 & ldk_range >= 0.5][,hi_lo := ifelse(L0/K<=0.5,'lo','hi')] %>% 
#   filter(hi_lo == 'lo') %>% 
#   nest(data=-species) %>% 
#   mutate(
#     fit = map(data, ~lm(r~I(L0/K), data=.x)),
#     tidied = map(fit, tidy)
#   ) %>% 
#   unnest(tidied)
# d2 <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.5][,hi_lo := ifelse(L0/K<=0.5,'lo','hi')] %>% 
#   filter(hi_lo == 'hi') %>% 
#   nest(data=-species) %>% 
#   mutate(
#     fit = map(data, ~lm(r~I(L0/K), data=.x)),
#     tidied = map(fit, tidy)
#   ) %>% 
#   unnest(tidied)
# 
# d1 %>% filter(term=='I(L0/K)') %>% 
#   pull(estimate) %>% 
#   hist
# vec_inc <- d2 %>% filter(term=='I(L0/K)') %>% 
#   filter(estimate <= 0) %>% 
#   filter(p.value <= 0.05) %>% 
#   pull(species)
# 
# 
# 
# p_test <- merge(fits,rank_species,by='species')[r2>0.3][is.na(species)==F][nobs>=50 & ldk_range >= 0.5][,hi_lo := ifelse(L0/K<=0.5,'lo','hi')] %>% 
#   .[species %in% vec_inc] %>% 
#   ggplot(data=.,aes(L0/K, r,group=hi_lo,color=hi_lo))+
#   geom_point(alpha=0.25,col='grey70')+
#   geom_smooth(method='lm')+
#   geom_rug()+
#   scale_color_viridis_d(option='H',end=0.9)+
#   facet_wrap(~species, scales = 'free_y')+
#   theme_linedraw()+
#   theme(panel.grid = element_blank())
# ggsave(p_test, filename = 'test.png', 
#   width=60, height=40, units='cm', dpi=400)
