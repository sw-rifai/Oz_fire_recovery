library(tidyverse)
library(stars); 
library(dtplyr)
library(data.table) 
library(lubridate)
library(arrow)
library(furrr)

# Load reference grid 
grid <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                           col_select = c("x","y","id","date"))
grid <- unique(grid[,.(x,y,id)])

# Load ALA MQ cleaned species occurence records
ala_mq <- fread("../data_general/ALA/Gallagher_cleaned/euc_occurrence_clean.csv") %>% 
  rename( 
    x=longitude,
    y=latitude) %>% 
  as.data.table() %>% 
  .[x>=140 & y <= -23 & y>=-40]

# 
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)


## Simplify species names 
fn_simp <- function(x){
  x <- str_remove(x, " sp\\.")
  x <- str_remove(x, " subsp\\.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}
ala_mq <- ala_mq[,species := fn_simp(taxon),by=seq_len(nrow(ala_mq))]



## Find the most frequent species record within ±0.15° ------
find1_mq <- function(din){
  dl <- 0.15
  px <- unique(din$x)
  py <- unique(din$y)
  vs <- ala_mq[x %between% c(px-dl,px+dl)][y %between% c(py-dl,py+dl)]$species
  din$dom_sp <- get_mode(vs)
  return(din)
}
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

vec_ids <- unique(grid$id)
vec1 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[1]]
vec2 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[2]]
vec3 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[3]]
vec4 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[4]]

grpn <- nrow(grid[id%in%vec1])
system.time(
  out1_mq <- grid[id%in%vec1][,
                              # {cat("progress",.GRP/grpn*100,"%\n");
                              find1_mq(.SD),
                              # },
                              by=.(id)]
)

grpn <- nrow(grid[id%in%vec2])
system.time(
  out2_mq <- grid[id%in%vec2][,
                              # {cat("progress",.GRP/grpn*100,"%\n");
                              find1_mq(.SD),
                              # },
                              by=.(id)]
)

grpn <- nrow(grid[id%in%vec3])
system.time(
  out3_mq <- grid[id%in%vec3][,
                              # {cat("progress",.GRP/grpn*100,"%\n");
                              find1_mq(.SD),
                              # },
                              by=.(id)]
)


grpn <- nrow(grid[id%in%vec4])
system.time(
  out4_mq <- grid[id%in%vec4][,
                              # {cat("progress",.GRP/grpn*100,"%\n");
                              find1_mq(.SD),
                              # },
                              by=.(id)]
)

out_mq <- rbindlist(list(out1_mq,out2_mq,out3_mq,out4_mq))
rm(out1_mq, out2_mq,out3_mq,out4_mq)

write_parquet(out_mq, sink="../data_general/proc_data_Oz_fire_recovery/ala-mq_1dom-sp_weibull-SEOZcoastal-grid.parquet")
