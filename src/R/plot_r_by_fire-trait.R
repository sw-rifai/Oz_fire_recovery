library(tidyverse)
library(stars); 
library(dtplyr)
library(data.table) 
library(lubridate)
library(arrow)
library(furrr)

oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

# load fits ----
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-05-19 18:18:29.parquet")
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>0][r<0.024][r2>0.5][month%in%c(9,10,11,12,1,2)][
  ,delta:=(L0/K)
][delta<=0.75]

# dominant species ---- 
dom <- read_parquet("../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-fit-1burn-locs.parquet")
dom <- dom[,.(x,y,id,dom_sp)]
dom50 <- dom[,.(nobs = .N), by='dom_sp'][,rank := frank(-nobs)][order(rank)][nobs>=25]

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

# Get fire traits from AusTraits using Nicolle 2006
nic <- at[dataset_id=="Nicolle_2006"][trait_name %in% c('fire_response',"regen_strategy")]
nic <- nic[,species := simp_name(taxon_name),by=seq_len(nrow(nic))][
  ,.(value = get_mode(value)), by=.(species,trait_name)
]
missing_sp_fr <- dom50$dom_sp[!dom50$dom_sp%in%unique(nic[trait_name=='fire_response']$species)]
missing_sp_rs <- dom50$dom_sp[!dom50$dom_sp%in%unique(nic[trait_name=='regen_strategy']$species)]

tmp_fr <- at[trait_name == "fire_response"][value %in% c("fire_killed","resprouts")][taxon_name%in%missing_sp][
  ,species := taxon_name][
  ,.(value = get_mode(value)), by=.(species)
]

tmp_rs <- at[trait_name == "regen_strategy"][taxon_name%in%missing_sp][
  ,species := taxon_name][
    ,.(value = get_mode(value)), by=.(species)
  ]


d_fr <- rbindlist(list(nic[trait_name=='fire_response'][,.(species,value)], tmp_fr))
d_rs <- rbindlist(list(nic[trait_name=='regen_strategy'][,.(species,value)], tmp_rs))


dom50 <- merge(dom50[,species := dom_sp], 
      d_fr %>% rename(fire_response=value) %>% as.data.table(), 
      by='species')

dom50 <- merge(dom50[,species := dom_sp], 
               d_rs %>% rename(regen_strategy=value) %>% as.data.table(),
               by='species') %>% 
  as.data.table()


dat <- merge(fits,dom,by=c("id","x","y"))
dat <- merge(dat, dom50, by=c("dom_sp"))

library(patchwork)
p1 <- dat %>% ggplot(data=.,aes(fire_response,r))+
  geom_boxplot(outlier.colour = NA)+
    labs(x='Fire Response')+
    coord_cartesian(ylim=c(0,0.015),expand=c(0,0))

p2 <- dat %>% ggplot(data=.,aes(regen_strategy,r))+
  geom_boxplot(outlier.colour = NA)+
  labs(x='Regeneration Strategy')+
  coord_cartesian(ylim=c(0,0.015),expand=c(0,0))

ggsave(p1+p2+plot_layout(widths = c(1,2)), 
       filename = 'figures/boxplot_r_by_fire-response_regen-strategy.png',
       width=25,
       height=10,
       dpi=350,
       units='cm')

