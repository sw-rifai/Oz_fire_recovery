pacman::p_load(tidyverse, stars, data.table, lubridate, patchwork)
setDTthreads(threads=24)
# Extract total obs fire pixels ------------------------------------------------
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
               col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
ba <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)

setkeyv(ba,cols=c("x","y"))

b_gt0 <- ba[date>=ymd("2000-09-01") & 
    date<=ymd("2020-02-01")][,.(val = sum(is.na(fire_doy)==F)),by=.(x,y)]

# what fraction burned at least once
100*nrow(b_gt0[val>0])/nrow(b_gt0)

# burned twice
100*nrow(b_gt0[val==2])/nrow(b_gt0)

# burned 3x
100*nrow(b_gt0[val==3])/nrow(b_gt0)

# burned 4x - very odd
100*nrow(b_gt0[val==4])/nrow(b_gt0)
