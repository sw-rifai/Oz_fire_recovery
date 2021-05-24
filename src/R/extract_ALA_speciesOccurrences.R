library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(fasterize); 
library(furrr)
library(arrow)

oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

# fits --- 
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-26 15:23:33.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>0][r<0.1]
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
fits <- merge(fits,d_soil,by='id')
rm(d_soil)

# AusTraits --- 
at <- fread("../data_general/AusTraits/austraits-2.1.0/data/traits.csv")
at[dataset_id=="Nicolle_2006"]$trait_name %>% unique
at[dataset_id=="Nicolle_2006"][trait_name=='fire_response']$value %>% table


# species spatial dats --- 
grid <- raster::raster("../data_general/proc_data_Oz_fire_recovery/grid_500m_SE_coastal.tif")
ala <- data.table::fread("../data_general/ALA/euc_records/euc-records-2021-05-04.csv")
ala <- ala %>% mutate(x=decimalLongitude, y=decimalLatitude) %>% as.data.table()
ala <- ala[between(x, min(fits$x), max(fits$x))]
ala <- ala[between(y, min(fits$y), max(fits$y))]
ala[,`:=`(speciesf = factor(species))]
ala[,`:=`(species_id = as.numeric(speciesf))]
ala <- ala[,.(x,y,species,speciesf,species_id)]
ala <- ala[str_length(species)>11] # filter out missing species records
d_nobs <- ala[,.(nobs = .N), by=species]
# d_nobs$nobs %>% hist(breaks=100)
d_nobs <- d_nobs[str_length(species)>0]
d_10 <- d_nobs[order(desc(nobs))][,rank := rank(-nobs)][rank <= 10] # top 10
d_20 <- d_nobs[order(desc(nobs))][,rank := rank(-nobs)][rank <= 20] # top 20 




fits[,`:=`(nsp_1 = NA_real_,
           nsp_2 = NA_real_,
           nsp_3 = NA_real_,
           nsp_4 = NA_real_,
           nsp_5 = NA_real_)]
find5 <- function(din){
  dl <- 0.1
  px <- unique(din$x)
  py <- unique(din$y)
  # vs <- NA_real_
  vs <- ala[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species_id %>% as.double()
  if(length(vs)==0){vs <- -1:-5}
  vs <- as.data.table(table(vs))[order(desc(N))]
  setnames(vs,c("spid","N"))
  vs[,rank:=frank(-N,ties.method = 'first')]
  tmp <- vs[rank<=5]
  if(nrow(tmp)<5){
    dl <- dl*2
    px <- unique(din$x)
    py <- unique(din$y)
    # vs <- NA_real_
    vs <- ala[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species_id %>% as.double()
    if(length(vs)==0){vs <- -1:-5}
    vs <- as.data.table(table(vs))[order(desc(N))]
    setnames(vs,c("spid","N"))
    vs[,rank:=frank(-N,ties.method = 'first')]
    tmp <- vs[rank<=5]
  }
  if(nrow(tmp)<5){
    dl <- dl*2
    px <- unique(din$x)
    py <- unique(din$y)
    # vs <- NA_real_
    vs <- ala[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species_id %>% as.double()
    if(length(vs)==0){vs <- -1:-5}
    vs <- as.data.table(table(vs))[order(desc(N))]
    setnames(vs,c("spid","N"))
    vs[,rank:=frank(-N,ties.method = 'first')]
    tmp <- vs[rank<=5]
  }
  din$nsp_1 <- tmp[rank==1,]$spid
  din$nsp_2 <- tmp[rank==2,]$spid
  din$nsp_3 <- tmp[rank==3,]$spid
  din$nsp_4 <- tmp[rank==4,]$spid
  din$nsp_5 <- tmp[rank==5,]$spid
  # append(0,vs)
  # if(length(vs)==0){vs <- vector(mode='numeric',length = 10)}
  # din$nsp[[1]] <- c(0,vs)
  return(din)
}

vec_ids <- unique(fits$id)
vec1 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[1]]
vec2 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[2]]
vec3 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[3]]
vec4 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[4]]

gc(full=TRUE)
plan(multisession, workers=20)
system.time(out1 <- fits[id%in%vec1] %>% 
              split(.$id) %>%
              future_map(~find5(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=TRUE)
system.time(out2 <- fits[id%in%vec2] %>% 
              split(.$id) %>%
              future_map(~find5(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=TRUE)
system.time(out3 <- fits[id%in%vec3] %>% 
              split(.$id) %>%
              future_map(~find5(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=TRUE)
system.time(out4 <- fits[id%in%vec4] %>% 
              split(.$id) %>%
              future_map(~find5(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=TRUE)
setDT(out1);setDT(out2);setDT(out3);setDT(out4)
out <- rbindlist(list(out1,out2,out3,out4),use.names = TRUE)
out[,`:=`(id=as.integer(id))]
out[,`:=`(nsp_1=as.integer(nsp_1))]
out[,`:=`(nsp_1=as.integer(nsp_2))]
out[,`:=`(nsp_1=as.integer(nsp_3))]
out[,`:=`(nsp_1=as.integer(nsp_4))]
out[,`:=`(nsp_1=as.integer(nsp_5))]

write_parquet(out[,.(x,y,id,nsp_1,nsp_2,nsp_3,nsp_4,nsp_5)], 
              sink = "../data_general/ala_5most-frequent-sp_refGrid.parquet")
# out <- read_parquet("../data_general/ala_5most-frequent-sp_refGrid.parquet")
# grpn <- uniqueN(fits$id)
# system.time(
#   out <- fits[,
#                       find5(.SD),
#               # {cat("progress",.GRP/grpn*100,"%\n"); 
#                 # find5(.SD)}, 
#               by=.(id)]
# )
# beat 0.859

get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
find1 <- function(din){
  dl <- 0.1
  px <- unique(din$x)
  py <- unique(din$y)
  # vs <- NA_real_
  vs <- ala[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species_id %>% as.integer()
  # vs <- as.data.table(table(vs))[order(desc(N))]
  # setnames(vs,c("spid","N"))
  # vs[,rank:=frank(-N,ties.method = 'first')]
  # tmp <- vs[rank==1]
  # if(nrow(tmp)<1){
  #   dl <- dl*2
  #   px <- unique(din$x)
  #   py <- unique(din$y)
  #   # vs <- NA_real_
  #   vs <- ala[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species_id %>% as.double()
  #   if(length(vs)==0){vs <- -1:-5}
  #   vs <- as.data.table(table(vs))[order(desc(N))]
  #   setnames(vs,c("spid","N"))
  #   vs[,rank:=frank(-N,ties.method = 'first')]
  #   tmp <- vs[rank<=5]
  # }
  din$dom_sp <- get_mode(vs)
  return(din)
}
system.time(
  out1 <- fits[,
                      find1(.SD),
              # {cat("progress",.GRP/grpn*100,"%\n");
                # find5(.SD)},
              by=.(id)]
)

write_parquet(out1[,.(x,y,id,dom_sp,r,K,L0,malai)], 
              sink = "../data_general/ala_1most-frequent-sp_refGrid.parquet")
# out1 <- read_parquet("../data_general/ala_1most-frequent-sp_refGrid.parquet")

count_sp <- function(din){
  dl <- 0.1
  px <- unique(din$x)
  py <- unique(din$y)
  # vs <- NA_real_
  vs <- ala[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species_id %>% as.integer()
  # vs <- as.data.table(table(vs))[order(desc(N))]
  # setnames(vs,c("spid","N"))
  # vs[,rank:=frank(-N,ties.method = 'first')]
  # tmp <- vs[rank==1]
  # if(nrow(tmp)<1){
  #   dl <- dl*2
  #   px <- unique(din$x)
  #   py <- unique(din$y)
  #   # vs <- NA_real_
  #   vs <- ala[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species_id %>% as.double()
  #   if(length(vs)==0){vs <- -1:-5}
  #   vs <- as.data.table(table(vs))[order(desc(N))]
  #   setnames(vs,c("spid","N"))
  #   vs[,rank:=frank(-N,ties.method = 'first')]
  #   tmp <- vs[rank<=5]
  # }
  din$num_sp <- length(unique(vs))
  return(din)
}
grpn <- uniqueN(fits$id)
system.time(
  out3 <- fits[,
               {cat("progress",.GRP/grpn*100,"%\n");
               count_sp(.SD)},
               by=.(id)]
)


# Dat merges for plotting --------------------------
out1 <- merge(out1, fits[,.(id,date_fire1)], by='id')
out1[,fire_month := month(date_fire1)]

sp_fac <- merge(out1[,.(dom_sp, r,L0,K)],
                unique(ala[,.(species_id,species)]), 
                by.x='dom_sp', by.y='species_id',
                allow.cartesian = T)[,.(species)] %>% 
  table %>% 
  as.data.table() %>% 
  setnames(., c("species","N")) %>% 
  .[order(-N)] %>% 
  .[,nr := .I] %>% 
  .[,sp_fac := forcats::fct_inorder(species)]

# What fraction of pixels covered by top 20 sp? --- 
sum(sp_fac[nr<=20]$N)/sum(sp_fac$N) 



# Start plotting -----------------------------------

# Plot 1: boxplot r ~ fire_response --------------------------------------------------------------
d_fire_response <- at[trait_name%in%c('fire_response','serotiny','wood_density')] %>% 
  .[,species:=taxon_name] %>% .[,.(species,trait_name,value,dataset_id)]
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
d_fire_response <- d_fire_response[,.(v1 = get_mode(value), 
                                      n_response = length(unique(value))),
                                   by=.(species,trait_name)]
at[trait_name=='fire_response'][,species:=taxon_name][,.(species,value,dataset_id)] %>% 
  .[species=="Eucalyptus regnans"]
d_fire_response[n_response==3]
at[taxon_name=="Banksia serrata"][trait_name=='fire_response']
d_fire_response[species=="Banksia serrata"]

merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), 
      by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[nr <= 30] %>% 
  merge(., d_fire_response[trait_name=='fire_response'],
        by='species',all.x=T,all.y=F,allow.cartesian=T) %>% 
  .[is.na(v1)==F] %>% 
  # .[,nobs:=.N,by=species] %>% 
  # .[,rank:=frankv(-nobs,ties.method = 'dense',na.last='keep')] %>% 
  # .[rank <= 20] %>% 
  # .[,.(rank)] %>% unique
  # .[order(rank)] %>% # 148 ranks, 179 speices
  .[(L0/K) %between% c(0,0.5)] %>% 
  ggplot(data=.,aes(sp_fac ,r,fill=v1))+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_viridis_d(option='H',direction = -1, begin = 0.2)+
  coord_flip(ylim=c(0,0.02))+
  scale_x_discrete(limits=rev)+
  labs(fill='fire response',
       x=NULL,
       y='r')+
  # coord_cartesian(ylim=c(0,0.01))+
  # scale_y_continuous(limits=c(0,0.02))+
  # facet_wrap(~cut_number(dom_sp,2),ncol=2,scales='free_y')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()) -> p1;p1
ggsave(p1, 
       filename = 'figures/boxplot_r_by_fire-response.png', 
       width=20, height=16,units='cm',dpi=350)


# Plot 1: boxplot r ~ regen_strategy --------------------------------------------------------------
d_fire_response <- at[trait_name%in%c('fire_response','serotiny','wood_density', 
                                      'regen_strategy')] %>% 
  .[,species:=taxon_name] %>% .[,.(species,trait_name,value,dataset_id)]
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
d_fire_response <- d_fire_response[,.(v1 = get_mode(value), 
                                      n_response = length(unique(value))),
                                   by=.(species,trait_name)]

merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), 
      by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[nr <= 30] %>% 
  merge(., d_fire_response[trait_name=='regen_strategy'],
        by='species',all.x=T,all.y=F,allow.cartesian=T) %>% 
  .[is.na(v1)==F] %>% 
  # .[,nobs:=.N,by=species] %>% 
  # .[,rank:=frankv(-nobs,ties.method = 'dense',na.last='keep')] %>% 
  # .[rank <= 20] %>% 
  # .[,.(rank)] %>% unique
  # .[order(rank)] %>% # 148 ranks, 179 speices
  .[(L0/K) %between% c(0,0.5)] %>% 
  ggplot(data=.,aes(sp_fac ,r,fill=v1))+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_viridis_d(option='H',direction = -1)+
  coord_flip(ylim=c(0,0.02))+
  scale_x_discrete(limits=rev)+
  labs(fill='regen_strategy',
       x=NULL,
       y='r')+
  # coord_cartesian(ylim=c(0,0.01))+
  # scale_y_continuous(limits=c(0,0.02))+
  # facet_wrap(~cut_number(dom_sp,2),ncol=2,scales='free_y')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()) -> p2;p2
ggsave(p2, 
       filename = 'figures/boxplot_r_by_regen_strategy.png', 
       width=20, height=16,units='cm',dpi=350)

# Plot 2 --------------------------------------------------
d_fire_response <- at[trait_name%in%c('fire_response','serotiny','wood_density', 
                                      'specific_leaf_area')] %>% 
  .[,species:=taxon_name] %>% .[,.(species,trait_name,value,dataset_id)] %>% 
  . [,.(v1 = get_mode(value), 
              n_response = length(unique(value))),
           by=.(species,trait_name)]

merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), 
      by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[nr <= 30] %>% 
  merge(., d_fire_response[trait_name=='specific_leaf_area'],
        by='species',all.x=T,all.y=F,allow.cartesian=T) %>% 
  .[is.na(v1)==F] %>% 
  .[,v1 := as.numeric(v1)] %>% 
  # .[,nobs:=.N,by=species] %>% 
  # .[,rank:=frankv(-nobs,ties.method = 'dense',na.last='keep')] %>% 
  # .[rank <= 20] %>% 
  # .[,.(rank)] %>% unique
  # .[order(rank)] %>% # 148 ranks, 179 speices
  .[(L0/K) %between% c(0,0.5)] %>% 
  ggplot(data=.,aes(sp_fac ,r,fill=v1))+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_viridis_c(option='G',direction = -1)+
  coord_flip(ylim=c(0,0.02))+
  scale_x_discrete(limits=rev)+
  labs(fill='wood density',
       x=NULL,
       y='r')+
  # coord_cartesian(ylim=c(0,0.01))+
  # scale_y_continuous(limits=c(0,0.02))+
  # facet_wrap(~cut_number(dom_sp,2),ncol=2,scales='free_y')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()) -> p2;p2

merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), 
      by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[nr <= 30] %>% 
  merge(., d_fire_response[trait_name=='specific_leaf_area'],
        by='species',all.x=T,all.y=F,allow.cartesian=T) %>% 
  .[is.na(v1)==F] %>% 
  .[,v1 := as.numeric(v1)] %>% 
  # .[,nobs:=.N,by=species] %>% 
  # .[,rank:=frankv(-nobs,ties.method = 'dense',na.last='keep')] %>% 
  # .[rank <= 20] %>% 
  # .[,.(rank)] %>% unique
  # .[order(rank)] %>% # 148 ranks, 179 speices
  .[(L0/K) %between% c(0,0.5)] %>% 
  ggplot(data=.,aes(v1,r))+
  geom_point()+
  geom_smooth(method='lm')

# Plot 4 -------------------------------------------
p4 <- out3 %>% ggplot(data=.,aes(x,y,fill=num_sp))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='H',
                       limits=c(0,90), 
                       oob=scales::squish)+
  labs(title='# species records within 0.1 degree radius')+
  theme_linedraw()+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1)); p4
ggsave(p4, filename = 'figures/map_malai_num-sp_lte0p1degree.png', 
       width=10, height=12,units='cm')

# Plot 5 -------------------------------------------
out3[num_sp<90] %>% lm(r~num_sp,data=.) %>% summary
p5 <- out3[num_sp<90] %>% ggplot(data=.,aes(num_sp,malai))+
  ggpointdensity::geom_pointdensity()+
  geom_smooth(se=F,lwd=2)+
  geom_rug()+
  scale_color_viridis_c(option='F')+
  labs(x='N Species',y='Mean Annual LAI')+
  theme(legend.position = 'none'); p5
ggsave(p5, filename = 'figures/scatterplot_malai_num-sp.png', 
       width=10, height=8,units='cm')

# Plot 6 -------------------------------------------
p6 <- merge(out1[,.(x,y,dom_sp, r,L0,K,fire_month)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[nr <= 20] %>% 
  ggplot(data=.,aes(x,y,fill=species))+
  geom_tile()+
  coord_sf(xlim=c(138,154))+
  scale_fill_manual(values=pals::turbo(20))+
  theme_linedraw()+
  guides(fill = guide_legend(ncol=2))+
  labs(x=NULL, 
       y=NULL,
       title="Dominant species record within 0.1 degree radius",
       fill='Species Record')+
  theme(legend.position = c(0,1), 
        legend.justification = c(0,1)); p6
ggsave(p6, 
       filename = 'figures/map_dom-species-record_lte0p1degree.png', 
       width=20, height=16,units='cm',dpi=350)


# Plot 7 -------------------------------------------
# overly complicated
merge(out1[,.(dom_sp, r,L0,K,fire_month)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[,delta:=L0/K] %>%
  .[delta <= 0.5] %>% 
  bam(r~s(delta,k=3,m=3)+
        s(delta, by=sp_fac, k=3) +
        s(sp_fac,bs='re'), 
      data=., 
      select=T,
      discrete = T) %>% 
  mgcViz::getViz() %>% 
  plot(allTerms=T) %>% 
  print(pages=1)

d_fit <- merge(out1[,.(dom_sp, r,L0,K,fire_month)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[,delta:=L0/K] %>% 
  lm(r~delta*sp_fac,data=.)
coef(d_fit) %>% as_tibble(., rownames=names(coef(d_fit))) %>% head
d_pos <- broom::tidy(d_fit) %>% 
  filter(str_detect(term,"delta:")) %>% 
  filter(estimate > 0) %>% 
  mutate(sp_pos = str_remove(term, "delta:sp_fac"))
d_neg <- broom::tidy(d_fit) %>% 
  filter(str_detect(term,"delta:")) %>% 
  filter(estimate < 0) %>% 
  mutate(sp_pos = str_remove(term, "delta:sp_fac"))
spx <- bind_rows(
  d_pos %>% select(sp_pos,estimate),
  d_neg %>% select(sp_pos,estimate)) %>% 
  mutate(xx = ifelse(estimate>0,'pos','neg')) %>% 
  rename(species=sp_pos) %>% 
  as.data.table()

merge(out1[,.(dom_sp, r,L0,K,fire_month)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  merge(., spx, by='species') %>% 
  # merge(., sp_fac, by='species') %>% 
  .[,delta:=L0/K] %>% 
  .[N >= 200] %>% 
  ggplot(data=.,aes(delta,r,color=species))+
  geom_point(inherit.aes = F,
             aes(delta,r),
             color='grey30',
             alpha=0.1)+
  geom_smooth(method='lm',se=F)+
  geom_smooth(inherit.aes = F,
              aes(delta,r),
              color='yellow',
              lwd=2,
              se=F)+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))+
  scale_x_continuous(expand=c(0.05,0.05))+
  geom_rug()+
  scico::scale_colour_scico_d(palette = 'romaO')+
  labs(x=expression(paste(L[0]/K)), 
       y='r', 
       caption = 'Only pixel locations with a dominant species containing > 200 obs.')+
  theme(legend.position = 'none',
        panel.background = element_rect(fill='grey10'), 
        panel.grid = element_blank()) -> p7; p7
ggsave(p7, 
       filename = 'figures/trend_r_L0K_by-dom-species.png', 
       width=20, height=16,units='cm',dpi=350)

# Plot 8 --------------------------------------------------------------
merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
  unique(ala[,.(species_id,species)]), 
  by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[nr <= 30] %>% 
  # .[,nobs:=.N,by=species] %>% 
  # .[,rank:=frankv(-nobs,ties.method = 'dense',na.last='keep')] %>% 
  # .[rank <= 20] %>% 
   # .[,.(rank)] %>% unique
  # .[order(rank)] %>% # 148 ranks, 179 speices
  .[(L0/K) %between% c(0,0.5)] %>% 
  ggplot(data=.,aes(sp_fac ,r,fill=N))+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_viridis_c(option='G',direction = -1)+
  coord_flip(ylim=c(0,0.02))+
  scale_x_discrete(limits=rev)+
  labs(fill='N pixels',
       x=NULL,
       y='r')+
  # coord_cartesian(ylim=c(0,0.01))+
  # scale_y_continuous(limits=c(0,0.02))+
  # facet_wrap(~cut_number(dom_sp,2),ncol=2,scales='free_y')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()) -> p8;p8
ggsave(p8, 
       filename = 'figures/boxplot_r_by-20dom-species.png', 
       width=20, height=16,units='cm',dpi=350)



# Scratch -----------------------------------------------------------------


merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), 
      by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  .[,.(med_r = median(r,na.rm=T)),by='species'] %>% 
  .[,r_rank := frank(-med_r,ties.method = 'dense')][]

r_fac <- merge(out1[,.(dom_sp, r,L0,K)],
                unique(ala[,.(species_id,species)]), 
                by.x='dom_sp', by.y='species_id',
                allow.cartesian = T)[,.(species,r)] %>% 
  .[,.(med_r = median(r,na.rm=T), 
       nobs = .N),by='species'] %>% 
  .[nobs >= 100] %>% 
  .[order(-med_r)] %>% 
  .[,nr := .I] %>%  
  .[,sp_fac := forcats::fct_inorder(species)]

merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), 
      by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., r_fac, by='species') %>% 
  .[nr <= 30] %>% 
  # .[,nobs:=.N,by=species] %>% 
  # .[,rank:=frankv(-nobs,ties.method = 'dense',na.last='keep')] %>% 
  # .[rank <= 20] %>% 
  # .[,.(rank)] %>% unique
  # .[order(rank)] %>% # 148 ranks, 179 speices
  .[(L0/K) %between% c(0,0.5)] %>%
  ggplot(data=.,aes(sp_fac ,r,fill=med_r))+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_viridis_c(option='G',direction = -1)+
  coord_flip(ylim=c(0,0.025))+
  scale_x_discrete(limits=rev)+
  labs(fill='median_r',
       x=NULL,
       y='r')+
  # coord_cartesian(ylim=c(0,0.01))+
  # scale_y_continuous(limits=c(0,0.02))+
  # facet_wrap(~cut_number(dom_sp,2),ncol=2,scales='free_y')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()) -> p9;p9


merge(out1[,.(fire_month,dom_sp, r,L0,K)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), 
      by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species') %>% 
  .[nr <= 30] %>% 
  # .[,nobs:=.N,by=species] %>% 
  # .[,rank:=frankv(-nobs,ties.method = 'dense',na.last='keep')] %>% 
  # .[rank <= 20] %>% 
  # .[,.(rank)] %>% unique
  # .[order(rank)] %>% # 148 ranks, 179 speices
  .[(L0/K) %between% c(0,0.2)]



tmp1 <- merge(out[,.(nsp_1, r,L0,K)],
 unique(ala[,.(species_id,species)]), by.x='nsp_1', by.y='species_id',allow.cartesian = T) 
tmp1[,nobs:=.N,by=species][,rank:=frank(-nobs,ties.method = 'dense')][order(rank)] %>% 
  .[(L0/K) %between% c(0,0.2)] %>% 
  .[,r_med := median(r,na.rm=T),by='species'] %>% 
  # .[,r_rank := frank(r_med,method='dense')] %>% 
  .[rank <= 20] %>%
  ggplot(data=.,aes(species,r,fill=nobs))+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_viridis_c()+
  coord_flip(ylim=c(0,0.02))+
  # coord_cartesian(ylim=c(0,0.01))+
  # scale_y_continuous(limits=c(0,0.02))+
  facet_wrap(~cut_interval(rank,2),ncol=2,scales='free_y')+
  theme()

tmp1[,nobs:=.N,by=species][,rank:=frank(-nobs,ties.method = 'dense')][order(rank)] %>% 
  # .[(L0/K) %between% c(0,0.2)] %>% 
  .[,r_med := median(r,na.rm=T),by='species'] %>% 
  ggplot(data=.,aes(nobs,r))+
  # geom_point(alpha=0.1)+
  geom_smooth(method='lm')+
  facet_wrap(~cut_number((L0/K),4),ncol=1,scales='free_y')


merge(as.data.table(sort(table(out$nsp_1)))[,rank:=frank(-N,ties.method = 'first')] %>% 
  rename(nsp_1=V1) %>% as.data.table(),
  out,by='nsp_1') %>% as_tibble() %>% 
  merge(., unique(ala[,.(species_id,species)]) %>% rename(nsp_1=species_id), by='nsp_1') %>% 
  filter(rank <= 20) %>% 
  ggplot(data=.,aes(x,y,fill=factor(species)))+
  geom_tile()+
  scale_fill_manual(values = pals::glasbey(20))+
  coord_equal()


ncol(find5(fits[1,]))
ncol(find5(fits[2,]))

dcast(vs[rank<=5],N~rank+spid)[]
tmp <- as_tibble(vs[rank<=5])
pivot_wider(tmp %>% select(-N),names_from=rank,values_from=spid,names_prefix = "nsp_")


find5(fits[7,]) %>% str
find5(fits[8,]) %>% str
find5(fits[9,]) %>% str
fits$vs <- list()
fits[1:10,]$id
for(i in 1:100){
  print(find5(fits[i,]))
  # print(paste(i,"+ ", ncol(fits)))
  # print(i)
  # if(ncol(find5(fits[i,]))<59){break}
}
fits$id
fits[9,]$id
din <- fits[9,]
din <- fits[7,]

sala <- sf::st_as_sf(ala[species%in%d_20$species][,.(x,y,species)], 
                     coords=c("x","y")) 
st_crs(sala) <- st_crs(4326)
st_coordinates(sala)


getmode <- function(v, ...) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
ala[,`:=`(speciesf = factor(species))]
ala[,`:=`(species_id = as.numeric(speciesf))]
sm <- raster::rasterize(ala[,.(x,y)], 
                        grid, 
                        ala$species_id,
                        getmode)

ala[species==d_20$species[1]][,.(x,y)] %>% as.matrix()
s1 <- raster::rasterize(ala[species==d_20$species[1]][,.(x,y)],
                        grid,
                        fun=function(x,...)length(x))
s1 <- raster::rasterize(ala[species==d_20$species[1]][,.(x,y,species)],
                        grid,
                        fun=function(x,...)length(x))
st_as_stars(s1) %>% 
  as_tibble() %>% 
  filter(is.na(layer)==F) %>% 
  pull(layer) %>% sum
ggplot()+
  geom_stars(data=st_as_stars(s1))+
  scale_fill_viridis_c()+
  coord_sf()

st_as_stars(s1) %>% 
  as_tibble() %>% 
  filter(is.na(layer)==F) %>% 
  pull(layer) %>% sum
plot(grid)

getmode <- function(v, ...) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}   
s1 <- terra::rasterize(x=st_coordinates(subset(sala, species==d_20$species[1])),
                 y=grid, 
                 fun='sum')
s1 <- terra::rasterize(x=ala[species==d_20$species[1]][,.(x,y)],
                       y=grid, 
                       fun=mean)
st_as_stars(s1) %>% 
  as_tibble() %>% 
  filter(is.na(layer)==F) %>% 
  pull(layer) %>% sum
d_20$species[1]
ala[species==d_20$species[1]] %>% nrow


bad_grid <- as(grid,Class = 'Raster')
s1 <- fasterize::fasterize(subset(sala, species==d_20$species[1]),grid,fun='sum')

s1 <- st_rasterize(subset(sala, species==d_20$species[1]), 
             template = grid, FUN=sum)
plot(s1,breaks='equal')
filter(sala,species==d_20$species[1])
sala %>% filter(species==d_20$species[1])
subset(sala, species==d_20$species[1])


ala[species%in%d_20$species][,.(x,y,species)] %>% 
  ggplot(data=.,aes(x,y,color=species))+
  geom_sf(data=oz_poly,color='black',fill='grey',inherit.aes = F)+
  geom_point(data=fits, aes(x,y),color='black')+
  geom_point(size=0.01)+
  scale_color_manual(values = pals::glasbey(20))+
  coord_sf(xlim=c(143.5,154),
           ylim=c(-39,-27.9))+
  labs(x=NULL,
       y=NULL,
       color='Species')+
  guides(color=guide_legend(override.aes = list(size=2)))+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank())


s_20 <- ala[species%in%d_20$species][,.(x,y,species)]
xy <- fits[,.(id,x,y)]


coords_vi <- lazy_dt(xy) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_covar <- unique(s_20[,.(x,y)])
coords_covar_sf <- st_as_sf(coords_covar, coords = c('x','y'))
st_crs(coords_covar_sf) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_covar_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=3
)
coords_covar <- coords_covar %>% mutate(idx_covar = row_number()) %>% as.data.table()
gc(full=TRUE)
coords_vi <- coords_vi %>% st_drop_geometry() %>% as.data.table()

coords_vi$idx1 <- coords_covar[nn_coords$nn.idx[,1],]$idx_covar
coords_vi$idx2 <- coords_covar[nn_coords$nn.idx[,2],]$idx_covar
coords_vi$idx3 <- coords_covar[nn_coords$nn.idx[,3],]$idx_covar
gc(full=TRUE)

# merges
gc(full=TRUE)
xy <- merge(xy,coords_vi,by='id')
s_20 <- merge(s_20,coords_covar,by=c("x","y"))

s_20
coords_covar

tmp1 <- merge(s_20[,.(species,idx_covar)],
      coords_vi,by.x='idx_covar',by.y='idx1') %>% 
  rename(sp1=species) %>% 
  select(id,sp1) %>% 
  as.data.table()
tmp2 <- merge(s_20[,.(species,idx_covar)],
              coords_vi,by.x='idx_covar',by.y='idx2') %>% 
  rename(sp2=species) %>% 
  select(id,sp2) %>% 
  as.data.table()
tmp3 <- merge(s_20[,.(species,idx_covar)],
              coords_vi,by.x='idx_covar',by.y='idx3') %>% 
  rename(sp3=species) %>% 
  select(id,sp3) %>% 
  as.data.table()

out <- merge(tmp1,tmp2,by='id',allow.cartesian = T)
out <- merge(out,tmp3,by='id',allow.cartesian = T)

out[,.(n_sp := length(unique(c(sp1,sp2,sp3)))),
    by='id']
out$id %>% length
out$id %>% unique %>% length()



# apply(my.df, 1, min)
apply(out,1,FUN = length)

out[n_sp==13]

out[id==130647]

clim <- merge(clim,coords_covar,by=c('x','y'))
gc(full=TRUE)
fits <- merge(fits, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
clim <- clim[idx_awap %in% unique(fits$idx_awap)]
