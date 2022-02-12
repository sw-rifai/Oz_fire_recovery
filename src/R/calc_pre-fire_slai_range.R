pacman::p_load(data.table,lubridate,magrittr)
tmp <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
tmp


tmp1 <- tmp[is.na(fire_doy)==F][,.(fire_date = min(date,na.rm=T)),by=.(id)]
gc()
setkeyv(tmp,cols=c("id"))

tmp <- merge(tmp,tmp1,by='id')

tmp_pre <- tmp[tmp[, .I[date < fire_date], by = .(id)]$V1]

tmp_sd <- tmp_pre[,.(lai_range = abs(diff(range(rnorm(10))))),by=.(x,y,id,year)]
tmp_range <- tmp_sd[,.(slai_range = median(lai_range)),by=.(id,x,y)]

arrow::write_parquet(tmp_range, 
  sink = "../data_general/proc_data_Oz_fire_recovery/pre-fire_slai_range.parquet", 
  compression='snappy')

