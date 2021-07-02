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
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-05-19 18:18:29.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>0][r<0.024][r2>0.5][month%in%c(9,10,11,12,1,2)][
  ,delta:=(L0/K)
  ][delta<=0.75]

# dominant species ---- 
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
dom <- dom[,.(x,y,id,dom_sp)]
dom20 <- dom[,.(nobs = .N), by='dom_sp'][,rank := frank(-nobs)][order(rank)]

# AusTraits --- 
at <- fread("../data_general/AusTraits/austraits-2.1.0/data/traits.csv")
at[dataset_id=="Nicolle_2006"]$trait_name %>% unique
at[dataset_id=="Nicolle_2006"][trait_name=='fire_response']$value %>% table
dom20[rank<=20]$dom_sp %in% unique(at[dataset_id=="Nicolle_2006"][trait_name=='fire_response']$taxon_name)

at[dataset_id=="Nicolle_2006"][trait_name=='fire_response']


# species bioclim ranges ---
hab <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_species-median-clim-topo.parquet")
names(hab) <- c("species",paste0("hab_",names(hab)[-1]))

# mean annual climate -------------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y","map","mapet","mappet",
                                           "matmax","matmin","mavpd15")) %>% 
  distinct() %>% 
  as.data.table()
rclim <- st_as_stars(clim)
st_crs(rclim) <- st_crs(4326)



# merge into dat  ---------------------
dat <- merge(fits,dom,by=c("id","x","y"))
dat <- merge(dat, setnames(hab, "species","dom_sp"), by=c("dom_sp"))

# attach climate to dat ---------------
coords <- st_as_sf(dat[,.(x,y)],coords=c("x","y"))
st_crs(coords) <- st_crs(4326)
tmp <- st_extract(rclim, coords) %>% st_drop_geometry()
dat <- bind_cols(dat,tmp) %>% as.data.table()
rm(tmp)

# attach landscape covars ------------
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
dat <- merge(dat,d_soil,by='id')

# attach height above nearest drainage ------
rhnd <- read_stars("../data_general/Oz_misc_data/mert_hnd_500m_SE_coastal.tif") %>% 
  set_names("hnd")
tmp <- st_extract(rhnd, coords) %>% st_drop_geometry()
dat <- bind_cols(dat,tmp) %>% as.data.table()
rm(tmp)


# mean annual climate -------------------------
post_clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y","date","post_precip_anom_12mo")) %>% 
  as.data.table()


# Filter dat to just Eucs -------------------
dat <- dat[vc %in% c(2,3,5)]
nobs <- dat[,.(nobs = .N), by=dom_sp][,rank:=frank(-nobs)]
sp_fac <- 
  unique(dat[dom_sp %in% nobs[rank <= 30]$dom_sp][,.(dom_sp,hab_hnd)]) %>% 
  .[order(hab_hnd)] %>% 
  .[,sp_fac := forcats::fct_inorder(dom_sp)]

dat[dom_sp %in% nobs[rank <= 30]$dom_sp] %>% 
  merge(., sp_fac,
        by='dom_sp',
        allow.cartesian=T) %>% 
  ggplot(data=.,aes(r,sp_fac,fill=hab_hnd.x))+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_viridis_c(option='H',direction = -1)+
  coord_cartesian(xlim=c(0,0.015))


dat[dom_sp %in% nobs[rank <= 20]$dom_sp] %>% 
  ggplot(data=.,aes(delta,r,color=dom_sp))+
  geom_point(alpha=0.1)+
  geom_smooth(color='black')+
  scale_color_viridis_d(option='H')+
  facet_wrap(~dom_sp)


dat %>% bam(r~hab_mavpd15+hab_mappet+hab_hnd,data=.) %>% getViz() %>% plot(allTerms=T)
dat %>% bam(r~hab_mavpd15+hab_mappet+hab_hnd,data=.) %>% getViz() %>% plot(allTerms=T)
dat %>% bam(r~scale(hab_mappet)+
              # scale(hab_map)+
              scale(hab_hnd),data=.) %>% summary

dat %>% bam(r~scale(hab_mappet),data=.) %>% summary
dat %>% bam(r~scale(mappet),data=.) %>% summary

dat$hnd %>% hist

b1 <- bam(r ~ s(delta, k=5) + 
            s(lai_yr_sd,k=5) + 
            s(silt,k=5)+
            s(pH,k=5)+
            s(mappet,k=5), 
          data=dat, 
          family=Gamma(link='log'),
          select=TRUE,
          discrete = TRUE)
summary(b1)
plot(b1)

b2 <- bam(r ~ s(delta, k=5) +
            s(malai,lai_yr_sd,k=5), 
            # te(aspect, mapet,bs=c('cc', 'cs'))+
            # s(elevation,k=5)+
            # s(pH,k=5)+
            # s(sand,k=5)+
            # s(map,k=5)+
            # s(mavpd15,k=5)+
            # s(matmax,k=5)+
            # te(map,mapet,elevation, k=5), 
          data=dat, 
          family=Gamma(link='log'),
          select=TRUE,
          discrete = TRUE)
summary(b2)
plot(b2,scheme=2)

dat[sample(.N,10000)] %>% 
  ggplot(data=.,aes(delta,r,color=lai_yr_sd))+
  geom_point(alpha=0.5)+
  geom_smooth(method='lm')+
  scale_color_viridis_c(option='H',direction=-1)+
  facet_wrap(~vc_name)


dat[,.(vc,vc_name)][order(vc)] %>% unique

