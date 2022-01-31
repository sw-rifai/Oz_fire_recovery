pacman::p_load(tidyverse, stars, data.table, lubridate, arrow, patchwork,mgcv)

# load fits ----
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-22 09:37:07.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>=0][r<0.024][r2>0.333][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
 ]
# estimate TTR from the logistic function
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]

out <- read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-14 16:45:26.parquet")
sout <- set_names(st_as_stars(out[,.(x,y,predict)], dims = c("x","y"),crs=4326),c("species"))
st_crs(sout) <- st_crs(4326)
vv <- st_extract(sout,
  st_as_sf(fits[,.(x,y)],coords=c("x","y"),crs=4326))
fits <- bind_cols(fits, st_drop_geometry(vv))


# vulnerable species 
fn_simp <- function(x){
  x <- str_replace(x, " x ", " ")
  x <- str_remove(x, " sp\\.")
  x <- str_remove(x, " subsp\\.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}


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

d_fire_response <- at_s[trait_name=='fire_response'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(fire_response = get_mode(value)),by=species]

merge(fits,d_fire_response,by='species') %>% 
  .[,.(n = .N),by=fire_response]

merge(fits,d_fire_response,by='species')$species %>% unique %>% length
merge(fits,d_fire_response,by='species')[,.(val = sum(fire_response=='fire_killed')),by=species][val > 0]

# How much of the burned region was dominated by 'fire_killed' species
100*nrow(merge(fits,d_fire_response,by='species')[fire_response=='fire_killed'])/
nrow(merge(fits,d_fire_response,by='species'))

# How much larger is r for resprouting than fire_killed species? [for main text]
m_fire_response <- merge(fits,d_fire_response,by='species') %>% 
  gam(r~fire_response,
    data=.,
    family=Gamma(link='log'), 
    method='REML')
tibble(fire_response = c("fire_killed","resprouts")) %>% 
  predict(m_fire_response, newdata=., type='response')

# How large is the difference between regeneration strategies? [for main text]
m_regen <- at_s[trait_name=='regen_strategy'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(regen = get_mode(value)),by=species] %>% 
  merge(., fits, by='species') %>%  
  gam(r~regen,
    data=.,
    family=Gamma(link='log'), 
    method='REML')
tibble(regen = at_s[trait_name=='regen_strategy'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(regen = get_mode(value)),by=species] %>% 
  merge(., fits, by='species') %>% pull(regen) %>% unique
  ) %>% 
  mutate(r_pred = predict(m_regen, newdata=., type='response')) %>% 
  arrange(r_pred)


p_1 <- merge(fits,d_fire_response,by='species') %>% 
  ggplot(data=.,aes(fire_response,r))+
  geom_boxplot(outlier.colour = NA)+
  labs(x="Fire Response",
    y=expression(paste(italic(r)~"(m²/day)")))+
  coord_cartesian(ylim=c(0,0.0125))+
  scale_x_discrete(guide=guide_axis(n.dodge = 2))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_1


p_2 <- at_s[trait_name=='regen_strategy'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(regen = get_mode(value)),by=species] %>% 
  merge(., fits, by='species') %>% 
  ggplot(data=.,aes(regen,r))+
  geom_boxplot(outlier.colour = NA)+
  labs(x="Regeneration Strategy",
    y=expression(paste(italic(r)~"(m²/day)")))+
  scale_x_discrete(guide=guide_axis(n.dodge = 2))+
  coord_cartesian(ylim=c(0,0.0125))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p_2


library(mgcViz)
# Wood density -----------------------------------------------------------------
m_wd_spat <- at_s[trait_name=='wood_density'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(trait = get_mode(value) %>% as.numeric()),by=species] %>% 
  merge(., fits, by='species') %>% 
  gam(r~trait+te(x,y),data=.,method='REML')
m_wd_simp <- at_s[trait_name=='wood_density'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(trait = get_mode(value) %>% as.numeric()),by=species] %>% 
  merge(., fits, by='species') %>% 
  gam(r~trait,data=.,method='REML')

summary(m_wd_spat)
summary(m_wd_simp)
methods(class=class(m_wd_spat))

tibble(trait = seq(0.4,1,length.out=10)) %>% 
  mutate(r_pred = predict(m_wd_simp, newdata=.)) %>% 
  summarize(p_diff = max(r_pred)/min(r_pred))

predict(m_wd_spat,
  terms='terms',
  exclude="te(x,y)",
  # type='iterms',
  newdata=tibble(trait = seq(0.4,1,length.out=10),
    x=runif(10),
    y=runif(10)
    )) %>% summary

p_wd <- tibble(trait = seq(0.4,1,length.out=10), 
  x=145,
  y=-30) %>% 
  mutate(pred_spat = predict(m_wd_spat,
  terms='terms',
  exclude="te(x,y)",
  # type='iterms',
  newdata=.
    )) %>% 
  mutate(pred_simp = predict(m_wd_simp,type='response',newdata=.)) %>% 
  select(trait, pred_spat,pred_simp) %>%
  pivot_longer(cols = c("pred_spat","pred_simp")) %>% 
  mutate(name = case_when(name=="pred_spat"~"linear model w/ spatial effect",
    name=="pred_simp"~"linear model")) %>% 
  ggplot(data=.,aes(trait,value,color=name))+
  geom_line()+
  labs(x=expression(paste(Wood~Density~"(g/cm³)")),
    y=expression(paste(italic(r)~"(m²/day)")), 
    color='model')+
    theme_linedraw()+
  scale_color_viridis_d(option='G',end=0.5)+
  coord_cartesian(ylim=c(0.004,0.007))+
  theme(panel.grid = element_blank())+
  theme(legend.position = c(0.99,0.01),
    legend.justification = c(0.99,0.01)); p_wd

# SLA -------------------------------------------------------------------------
at_s[trait_name=='specific_leaf_area'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(trait = get_mode(value) %>% as.numeric()),by=species] %>% 
  merge(., fits, by='species') %>% 
  bam(r~trait+te(x,y),data=.) %>% 
  getViz() %>%
  pterm(1) %>% 
  plot(allTerms=T) %>% 
  print(pages=1)
m_sla_spat <- at_s[trait_name=='specific_leaf_area'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(trait = get_mode(value) %>% as.numeric()),by=species] %>% 
  merge(., fits, by='species') %>% 
  gam(r~trait+te(x,y),data=.) 
m_sla_simp <- at_s[trait_name=='specific_leaf_area'][species%in%unique(fits$species)][
  ,.(species,trait_name,value)
][,.(trait = get_mode(value) %>% as.numeric()),by=species] %>% 
  merge(., fits, by='species') %>%
  lm(r~trait,data=.) 

summary(m_sla_simp)
tibble(trait = seq(2,17,length.out=10)) %>% 
  mutate(r_pred = predict(m_sla_simp, newdata=.)) %>% 
  summarize(p_diff = max(r_pred)/min(r_pred))


predict(m_sla_spat) %>% summary
predict(m_sla_spat,
  terms='terms',
  exclude="te(x,y)",
  # type='iterms',
  newdata=tibble(trait = seq(2,17,length.out=10),
    x=runif(10),
    y=runif(10)
    )) %>% summary
p_sla <- tibble(trait = seq(2,17,length.out=10)) %>% 
  mutate(pred_spat = predict(m_sla_spat,
  terms='terms',
  exclude="te(x,y)",
  # type='iterms',
  newdata=tibble(trait = seq(2,17,length.out=10),
    x=runif(10),
    y=runif(10)
    ))) %>% 
  mutate(pred_simp = predict(m_sla_simp,type='response',newdata=.)) %>% 
  select(trait, pred_spat,pred_simp) %>% 
  pivot_longer(cols = c("pred_spat","pred_simp")) %>% 
  mutate(name = case_when(name=="pred_spat"~"linear model w/ spatial effect",
    name=="pred_simp"~"linear model")) %>% 
  ggplot(data=.,aes(trait,value,color=name))+
  geom_line()+
  labs(x=expression(paste(Specific~Leaf~Area~"(m²/kg)")),
    y=expression(paste(italic(r)~"(m²/day)")), 
    color='model')+
    theme_linedraw()+
  scale_color_viridis_d(option='G',end=0.5)+
  coord_cartesian(ylim=c(0.004,0.007))+
  theme(panel.grid = element_blank())+
  theme(legend.position = c(0.99,0.01),
    legend.justification = c(0.99,0.01)); p_sla



# Join panels -------------------------------------------------------------
p_out <- (p_wd|p_sla)/(p_1|p_2)+plot_annotation(tag_levels='a',
  tag_prefix = '(',
  tag_suffix = ')')
p_out
ggsave(p_out,
  filename = "figures/fig-s9_r_by_trait.png",
  width=25,
  height=15,
  units='cm',
  dpi = 350)



# p_3 <- at_s[trait_name=='wood_density'][species%in%unique(fits$species)][
#   ,.(species,trait_name,value)
# ][,.(trait = get_mode(value) %>% as.numeric()),by=species] %>% 
#   merge(., fits, by='species') %>% #lm(r~trait, data=.) %>% summary
#   ggplot(data=.,aes(trait,r))+
#   labs(x=expression(paste(Wood~Density~"(g/cm³)")),
#     y=expression(paste(italic(r)~"(m²/m²)")))+
#   geom_point(alpha=0.25)+
#   geom_smooth(method='lm')+
#   # coord_cartesian(ylim=c(0,0.0125))+
#   theme_linedraw()+
#   theme(panel.grid = element_blank()); p_3
# 
# 
# merge(fits,d_fire_response,by='species') %>% 
#   gam(r~fire_response,data=.,method='REML') %>% 
#   summary()
#     # data=. %>% .[,`:=`(fire_response:=factor(fire_response))]) %>% 
#   summary()
# 
# 
# at_s[trait_name=='fire_response'][species%in%unique(fits$species)][
#   ,.(species,trait_name,value)
# ] %>% unique %>%
#   .[,.(val = length(unique(value))),by=.(species)] %>% pull(val)
# 
# at_s[trait_name=='fire_response'][species%in%unique(fits$species)][
#   species=='Eucalyptus delegatensis'
# ][,.(val = get_mode(value)),by=species]
# 
# library(mgcv)
# b1 <- bam(r~s(L0,K)+species,
#   family=Gamma(link='log'),
#   data=fits)
# summary(b1)
# 
# 
# tmp_1 <- merge(unique(fits[,.(species)])[,`:=`(L0=0)], 
#   fits[,.(K=median(K)),by=species], 
#   by='species') %>% 
#   .[,`:=`(pred_r = predict(b1,type='response',newdata = .))]
# tmp_1[species=='Eucalyptus regnans']
# tmp_1 %>% 
#   ggplot(data=.,aes(K,pred_r))+
#   geom_point()+
#   geom_label(aes(label=species),size=3)
# 
# fits[species=='Eucalyptus regnans']$r %>% median
# merge(fits %>% 
#   filter(is.na(species)==F) %>% 
#   filter(is.na(r)==F) %>% 
#   group_by(species) %>% 
#   summarize(r_med = median(r,na.rm=T)) %>% 
#   ungroup() %>% 
#   arrange(r_med) %>% 
#   mutate(rank = row_number()) %>% 
#   mutate(r_species = factor(species,
#     levels = species,
#     ordered = T)), 
#   fits, by='species') %>% 
#   ggplot(data=.,aes(x=r,y=r_species))+
#   geom_boxplot(outlier.colour = NA)+
#   theme(text = element_text(size=6))
