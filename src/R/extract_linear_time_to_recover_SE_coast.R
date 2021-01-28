# Description: Prototype script to calculate the time to recover from
# fire in SE Australian Eucalyptus dominant forests.
# Author: Sami Rifai
# Date (init): 2020-01-06


# Load packages in this order (important)
library(phenofit);
library(tidyverse);
library(usethis);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(RcppArmadillo)
source("src/R/functions_time_to_recover.R")

# Load data ---------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
 sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)


dat <- arrow::read_parquet(file ="/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")

# Calculate anomalies -----------------------------------------------------
gc(full=TRUE)
dat_norms <- dat[, `:=`(month = month(date))] %>%
  .[, .(ndvi_u = mean(sndvi, na.rm=TRUE),
        ndvi_mmax = max(sndvi,na.rm=TRUE),
        ndvi_sd = sd(sndvi, na.rm=TRUE)),
    keyby = .(x,y,month)]
dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>%
  mutate(ndvi_anom = sndvi-ndvi_u) %>%
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_sd,
         ndvi_fanom = sndvi/ndvi_u) %>%
  mutate(ndvi_fmax = sndvi/ndvi_mmax) %>% 
  as.data.table()
# dat[,`:=`(ndvi_mmax = max(sndvi,na.rm=TRUE)), keyby=.(x,y,month)]
# dat[,`:=`(ndvi_fmax = sndvi/ndvi_mmax), keyby=.(x,y,month)]
gc(full=TRUE)


# Remove bad grid cells ----------------------------------
# some locs have consistently negative ndvi; salt beds?
bad_pix <- dat %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(val = median(sndvi,na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(val <= 0.15) %>% 
  as.data.table()
  
# dat %>% lazy_dt() %>% 
#   group_by(x,y) %>% 
#   summarize(val = median(sndvi,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(x,y,fill=val))+
#   geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
#   geom_tile()+
#   coord_sf(xlim = c(140,154),
#            ylim = c(-40,-25), expand = FALSE)+
#   scale_fill_gradient2()

dat <- dat[!id %in% bad_pix$id]; gc(full=TRUE)
ldat <- dat %>% lazy_dt()

# Apply Linear TTR Fns --------------------------------------------------
nburns <- ldat %>% group_by(id) %>% 
  summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>%
  as.data.table()
gc(full=TRUE)
burn_dates <- dat[id%in%unique(nburns$id)][is.na(fire_doy)==FALSE] %>% 
  .[, .(date_first_burn = min(date, na.rm = TRUE), 
        date_last_burn = max(date,na.rm = TRUE)), keyby = .(id)]
gc(full=TRUE)
burn_dates <- merge(burn_dates, nburns, by='id')
arrow::write_parquet(burn_dates, sink = "outputs/first_last_burn_dates.parquet")



# Process 2003-2004 fires --------------------------------------------
vec_2003 <- dat %>% lazy_dt() %>% 
  filter(between(year,2003,2004)) %>%
  # filter(is.na(fire_doy)==FALSE) %>% 
  group_by(id) %>% 
  summarize(
    count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  filter(count==1) 

grpn <- uniqueN(vec_2003$id)
out03 <- dat[id%in%vec_2003$id][,{cat("progress",.GRP/grpn*100,"%\n"); time_to_recover_vi_v6(.SD)}, by=.(id,x,y)]
arrow::write_parquet(out03, sink = paste0("outputs/linear_ttr_2003-2004_",Sys.Date(),".parquet"))
#*********************************************************************

# Process 2006-2007 fires --------------------------------------------
vec_2006 <- dat %>% lazy_dt() %>% 
  filter(between(year,2006,2007)) %>%
  # filter(is.na(fire_doy)==FALSE) %>% 
  group_by(id) %>% 
  summarize(
    count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  filter(count==1) 

# grpn <- uniqueN(dat$id)
grpn <- uniqueN(dat$id)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out06 <- dat[id%in%vec_2006$id][,{setTxtProgressBar(pb, .GRP); time_to_recover_vi_v6(.SD)}, by=.(x,y,id)]
# dt[, {setTxtProgressBar(pb, .GRP); Sys.sleep(0.5); sum(a)}, b]
close(pb)
arrow::write_parquet(out06, sink = paste0("outputs/linear_ttr_2003-2004_",Sys.Date(),".parquet"))

# Process 2012-2015 fires --------------------------------------------
vec_2012 <- dat %>% lazy_dt() %>% 
  filter(between(year,2012,2015)) %>%
  # filter(is.na(fire_doy)==FALSE) %>% 
  group_by(id) %>% 
  summarize(
    count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  filter(count==1) 

grpn <- uniqueN(vec_2012$id)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out12 <- dat[id%in%vec_2012$id][,{setTxtProgressBar(pb, .GRP); time_to_recover_vi_v6(.SD)}, by=.(x,y,id)]
close(pb)
arrow::write_parquet(out12, sink = paste0("outputs/linear_ttr_2015-2015_",Sys.Date(),".parquet"))


# Multi fires ---------------------------------------------------
vec_tmp <- nburns[nburns>0]
grpn <- uniqueN(vec_tmp$id)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out_tmp <- dat[id%in%vec_tmp$id][,{setTxtProgressBar(pb, .GRP); time_to_recover_vi_v7(.SD)}, by=.(x,y,id)]
close(pb)
arrow::write_parquet(out_tmp, sink = paste0("outputs/linear_ttr_multiBurns_2001-2020_",Sys.Date(),".parquet"))



# Plot post burn linear trend by month --------------------------
mb <- arrow::read_parquet("outputs/linear_ttr_multiBurns_2001-2020_2021-01-20.parquet")
mb %>% lazy_dt() %>% 
  mutate(month_fire = month(date_first_fire)) %>% 
  as.data.table() %>% 
  drop_na(month_fire) %>% 
  ggplot(data=.,aes(factor(month_fire),ttr_beta_12mo))+
  geom_boxplot()+
  labs(x='Month', y=expression(paste(Delta~NDVI~'12 mo'**-1)))




# Plot TTR across epochs --------------------------------------------
out03$epoch <- '03-04'
out06$epoch <- '06-07'
out12$epoch <- '12-15'
d_ttr <- rbindlist(list(out03,out06,out12), use.names = TRUE)

d_ttr %>% 
  filter(between(delta_vi_12mo,-0.35,0)) %>% 
  ggplot(data=.,aes(delta_vi_12mo,ttr,color=epoch))+
  geom_smooth()+
  scale_color_viridis_d()+
  labs(x='NDVI change from fire', y='Time to recover (days)')+
  theme_linedraw()
ggsave(filename = 'figures/compare_epochs_ttr.png')

ggplot()+
  geom_point(data=out03, aes(delta_vi_12mo,ttr),color="#CF0000",alpha=0.01)+
  geom_point(data=out06, aes(delta_vi_12mo,ttr),color="black",alpha=0.01)+
  geom_point(data=out12, aes(delta_vi_12mo,ttr),color="blue",alpha=0.01)





# Plot timeseries of burn counts --------------------------------
gc()
clim <- arrow::read_parquet("/home/sami/scratch/ARD_ndvi_aclim_anoms.parquet", 
                           col_select = c("x","y","date",
                                          'pe_anom_3mo',
                                          'pe_anom_12mo',
                                          'precip_anom_12mo',
                                          'vpd15_anom_12mo')) %>% 
      setDT()
gc(full=TRUE)
clim <- clim %>% lazy_dt() %>% 
  mutate(x=round(x),
         y=round(y)) %>% 
  as.data.table()
vec_xy <- ldat %>% select(x,y) %>% 
  distinct() %>% 
  mutate(x=round(x),
         y=round(y)) %>% 
  distinct() %>% 
  as.data.table()
clim <- clim[x%in%vec_xy$x][y%in%vec_xy$y][date>=ymd("2001-01-01")]
clim <- merge(vec_xy,clim,by=c("x","y"),all = FALSE)
clim_ts <- clim %>% 
  lazy_dt() %>% 
  filter(is.infinite(pe_anom_12mo)==FALSE) %>% 
  group_by(date) %>% 
  summarize(ppet_anom_3mo = mean(pe_anom_3mo,na.rm=TRUE), 
            ppet_anom_12mo = mean(pe_anom_12mo,na.rm=TRUE), 
            precip_anom_12mo = mean(precip_anom_12mo,na.rm=TRUE), 
            vpd15_anom_12mo = mean(vpd15_anom_12mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
ldat %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  group_by(date) %>% 
  summarize(val = sum(fire_doy>0)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  ggplot(data=.,aes(date, val*25))+
  geom_rect(data=clim_ts, inherit.aes = FALSE, 
            aes(xmin=date,xmax=date+months(1), 
            ymin=0,ymax=80000*30,fill=precip_anom_12mo))+
  geom_line(col='black',lwd=1)+
  labs(y='Burn Area (ha)')+
  scale_fill_gradient2(expression(paste('Precip anom'['12mo']~' (mm)')))+
  theme_linedraw()+
  theme(legend.position = 'bottom')
ggsave(filename = "figures/timeSeries_burnArea_wPrecipAnom_2001-2020.png", 
       width=17,height=12,units='cm')

# Time series of 12-mo recovery rate

mb %>% lazy_dt() %>% 
  mutate(year_fire = year(date_first_fire)) %>% 
  as.data.table() %>% 
  drop_na(year_fire) %>% 
  filter(between(year_fire, 2001,2019)) %>% 
  ggplot(data=.,aes(factor(year_fire),ttr_beta_12mo))+
  geom_rect(data=clim_ts %>% mutate(year_fire=factor(year(date))) %>% 
              group_by(year_fire) %>% summarize(val = mean(precip_anom_12mo,na.rm=TRUE)) %>% 
              ungroup(), inherit.aes = FALSE, 
            aes(xmin=year_fire,xmax=year_fire+1, 
                ymin=-0.1,ymax=0.1,fill=val))+
  # geom_boxplot()+
  labs(x='Year', y=expression(paste(Delta~NDVI~'12 mo'**-1)))


ggplot()+
  geom_rect(data=clim_ts %>% mutate(year_fire=year(date)) %>% 
              group_by(year_fire) %>% summarize(val = mean(precip_anom_12mo,na.rm=TRUE)) %>% 
              ungroup(), inherit.aes = FALSE, 
            aes(xmin=year_fire,xmax=year_fire+1, 
                ymin=-0.1,ymax=0.100,fill=val))+
  geom_boxplot(data=mb %>% lazy_dt() %>% 
                mutate(year_fire = year(date_first_fire)) %>% 
                as.data.table() %>% 
                drop_na(year_fire) %>% 
                filter(between(year_fire, 2001,2019)), 
              aes(year_fire+0.5, ttr_beta_12mo, group=year_fire)) +
  scale_fill_gradient2("Precip Anom. (mm yr-1)")+
  labs(x='year',y='NDVI anom trend post fire (12 mo)')
ggsave(filename = "figures/timeSeries_ndviTrend_wPrecipAnom_2001-2020.png",
       width=16,height=8,units='cm')

mb %>% 
  filter(between(date_first_fire, ymd("2003-01-01"),ymd("2004-12-31"))) %>% 
  ggplot(data=.,aes(delta_vi_12mo, ttr))+
  ggpointdensity::geom_pointdensity()+
  scale_color_viridis_c()



# Plot time between burns ---------------------------------------
tmp <- dat[id%in%burn_dates$id] %>% 
  lazy_dt() %>%
  filter(is.na(fire_doy)==FALSE) %>% 
  group_by(id) %>% 
  arrange(date) %>% 
  mutate(tbb = decimal_date(date)-decimal_date(lag(date,1))) %>%
  ungroup() %>% 
  as.data.table()

tmp %>% 
  filter(is.na(tbb)==F) %>% 
  group_by(x,y) %>% 
  filter(tbb == min(tbb,na.rm=TRUE)) %>% 
  ungroup() %>% 
  ggplot(data=.,aes(x,y,fill=tbb))+
    geom_sf(data=oz_poly, inherit.aes = F, fill='gray70',color='gray30')+
    geom_tile()+
    coord_sf(xlim = c(145,154),
             ylim = c(-40,-27), expand = FALSE)+
  scale_fill_viridis_c('years',option='B', direction=-1)+
  labs(x=NULL,y=NULL,title='Time between burns (2001-2020)')+
  theme(panel.background = element_rect(fill='lightblue'))
ggsave(filename = 'figures/timeBetweenBurns_2001-2020.png', 
       units='cm', width=10, height=15)

tmp %>% 
  group_by(x,y) %>% 
  summarize(ff = decimal_date(min(date,na.rm=TRUE))) %>% 
  ungroup() %>% 
  ggplot(data=.,aes(x,y,fill=ff))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray70',color='gray30')+
  geom_tile()+
  coord_sf(xlim = c(145,154),
           ylim = c(-40,-27), expand = FALSE)+
  scale_fill_viridis_c('Year',option='B', direction=-1)+
  labs(x=NULL,y=NULL,title='First Fire (2001-2020)')+
  theme(panel.background = element_rect(fill='lightblue'))
ggsave(filename = 'figures/firstFireDate_2001-2020.png', 
       units='cm', width=10, height=15)

tmp %>% 
  group_by(x,y) %>% 
  summarize(ff = decimal_date(max(date,na.rm=TRUE))) %>% 
  ungroup() %>% 
  ggplot(data=.,aes(x,y,fill=ff))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray70',color='gray30')+
  geom_tile()+
  coord_sf(xlim = c(145,154),
           ylim = c(-40,-27), expand = FALSE)+
  scale_fill_viridis_c('Year',option='B', direction=-1)+
  labs(x=NULL,y=NULL,title='Last Fire (2001-2020)')+
  theme(panel.background = element_rect(fill='lightblue'))
ggsave(filename = 'figures/LastFireDate_2001-2020.png', 
       units='cm', width=10, height=15)




library(data.table)
dt = data.table(a=1:4, b=c("a","b"))
grpn = uniqueN(dt$b)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
dt[, {setTxtProgressBar(pb, .GRP); Sys.sleep(0.5); sum(a)}, b]
close(pb)

# system.time(sdat1 <- dat[,time_to_recover(.SD), by=.(x,y,id)]) # 1197
# system.time(sdat2 <- dat[,time_to_recover_fmax(.SD), by=.(x,y,id)]) # 3010
# system.time(sdat3 <- dat[,time_to_recover_fmax_v2(.SD), by=.(x,y,id)]) # 
# 
# arrow::write_parquet(sdat1, sink=paste0('outputs/time_to_recover_',Sys.Date()))
# arrow::write_parquet(sdat2, sink=paste0('outputs/time_to_recover_fmax_',Sys.Date()))
# arrow::write_parquet(sdat3, sink=paste0('outputs/time_to_recover_fmax_v2_',Sys.Date()))
# sdat1 <- arrow::read_parquet("outputs/time_to_recover_2021-01-10")
# sdat2 <- arrow::read_parquet("outputs/time_to_recover_fmax_2021-01-10")
# sdat3 <- arrow::read_parquet("outputs/time_to_recover_fmax_v2_2021-01-11")





data.table::setDTthreads(20)
vec_ids <- sort(unique(dat$id))
system.time(sdat1 <- dat[id%in%seq(1,floor(length(vec_ids)/4))][,time_to_recover(.SD), by=.(x,y)]) # 290
gc(full=TRUE)
# 



library(foreach); library(iterators)
mc <- parallel::makeForkCluster(n.cores=20)
print(mc)
doParallel::registerDoParallel(cl=mc)
foreach::getDoParRegistered()


junk <- dat[id %in% sample.int(1e6, 100000)]
vec_idx <- floor(seq.int(1,1e6,length.out = 21))

system.time(
out <- foreach::foreach(i=2:21, 
                        # .combine=data.table::rbindlist,
                        .packages=c("data.table")) %dopar%
  {
    junk[between(id,vec_idx[i],vec_idx[i+1])][,time_to_recover_fmax_v4(.SD), by=.(x,y,id)]
  }) # 65
# rbindlist(out)
system.time(junk[,time_to_recover_fmax_v5(.SD), by=.(x,y,id)]) # ~120 ; 

doParallel::stopImplicitCluster()
parallel::stopCluster(cl=mc)
gc(full=TRUE)

# 
# split(dat[id%in%c(1:5)],)
# 
# system.time(dat[id%in%c(1:100)][,time_to_recover(.SD),by='id'])
# system.time(dat[id%in%c(1:100)][,time_to_recover(.SD),by=.(x,y)])
# system.time(dat[id%in%c(1:100)][,time_to_recover(.SD),by=.(x,y,id)])
# 
# (1+length(vec_ids)/4):(2*length(vec_ids)/4)
# (1+2*length(vec_ids)/4):(3*length(vec_ids)/4)
# (1+3*length(vec_ids)/4):(4*length(vec_ids)/4)
# 
# library(sparklyr)
# spark_install()
# sc <- spark_connect(master='local')
# dat_s <- copy_to(dest=sc, df=dat, name='dat', overwrite = TRUE)
# 
# 
# 

#*******************************************************************************
# Extract AWAP tmax grid cells for east Oz coastal ROI -------------------------
#*******************************************************************************
tmp <- arrow::read_parquet("/home/sami/scratch/ARD_ndvi_aclim_anoms.parquet")
setDT(tmp)
gc(full=TRUE)

xmin <- min(unique(dat$x),na.rm=TRUE)
xmax <- max(unique(dat$x),na.rm=TRUE)
ymin <- min(unique(dat$y),na.rm=TRUE)
ymax <- max(unique(dat$y),na.rm=TRUE)

atmax <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/tmax/AWAP_monthly_tmax_1970_2019.nc")
names(atmax) <- "tmax"
st_crs(atmax) <- st_crs(4326)

eoz_box <- st_bbox(c(xmin = xmin,
                     ymin = ymin,
                     xmax = xmax,
                     ymax = ymax), 
                   crs = st_crs(4326))
atmax <- st_crop(atmax, eoz_box)
atmax <- atmax %>% as.data.table()
atmax <- atmax %>% units::drop_units()
coords_awap <- atmax %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))


coords_vi <- dat %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])



# vis check that vi and clim coords are close
coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
coords_dict %>% head

coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

atmax <- atmax %>% rename(x=longitude,y=latitude)
gc(full=TRUE)
atmax <- atmax[time > ymd("1995-01-01")]

atmax <- atmax[,.(x_vi,y_vi,time,tmax)]
atmax <- atmax %>% rename(x=x_vi, y=y_vi)
atmax <- atmax %>% rename(date=time)
atmax <- atmax %>% lazy_dt() %>% 
  mutate(date = as.Date(date)) %>% 
  as.data.table()

gc(full=TRUE)






# complicated way of doing full join
dat <- merge(atmax,
             dat,
             by=c("x","y","date"),all=TRUE, allow.cartesian = TRUE)

atmax <- merge(atmax,dat,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
atmax <- atmax[is.na(x_vi)==F]
atmax %>% head

round(atmax$x[1:5])==round(atmax$x_vi[1:5])
round(atmax$y[1:5])==round(atmax$y_vi[1:5])

# visual check
atmax[time%in%c(ymd("1990-01-01",tz='UTC'),
                ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=tmax))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))
gc()

unique(atmax$x_vi) %in% unique(dat$x) %>% table
#*******************************************************************************
# END SECTION
#*******************************************************************************
