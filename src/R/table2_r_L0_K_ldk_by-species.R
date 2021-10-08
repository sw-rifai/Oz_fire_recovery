pacman::p_load(tidyverse, stars, data.table, lubridate, arrow, patchwork, gt)

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

fn_quantile_lo <- function(x) quantile(x,c(0.05),na.rm=T)
fn_quantile_med <- function(x) quantile(x,c(0.5),na.rm=T)
fn_quantile_hi <- function(x) quantile(x,c(0.95),na.rm=T)
t2 <- fits %>% 
  filter(is.na(species)==F) %>% 
  select(species,r,L0,K,ldk) %>% 
  group_by(species) %>% 
  summarize_all(.funs = c(fn_quantile_lo,fn_quantile_med,fn_quantile_hi)) %>% 
  ungroup()

t2 <- t2 %>% select(species,
  starts_with('r'),
  starts_with("L0"),
  starts_with("K"),
  starts_with("ldk"))

vec_names <- str_replace(names(t2),"_fn1","_p05")
vec_names <- str_replace(vec_names,"_fn2","_p50")
vec_names <- str_replace(vec_names,"_fn3","_p95")

names(t2) <- vec_names

str_replace(vec_names,"_"," ")
t_out <- t2 %>% 
  gt() %>% 
  fmt_number(columns = 2:4, decimals=3) %>% 
  fmt_number(columns = 5:13, decimals=2)

t_out %>% 
  gtsave("doc/Table_S2.rtf")

t2 %>% 
  write_csv(., "doc/Table_S2.csv")

