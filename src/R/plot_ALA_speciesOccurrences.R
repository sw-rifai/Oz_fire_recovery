library(tidyverse); 
library(dtplyr); 
library(data.table)
library(sf); library(stars)
library(fasterize)

oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)
grid <- raster::raster("../data_general/proc_data_Oz_fire_recovery/grid_500m_SE_coastal.tif")
ala <- data.table::fread("../data_general/ALA/euc_records/euc-records-2021-05-04.csv")
ala <- ala %>% mutate(x=decimalLongitude, y=decimalLatitude) %>% as.data.table()
ala <- ala[between(x, min(fits$x), max(fits$x))]
ala <- ala[between(y, min(fits$y), max(fits$y))]
d_nobs <- ala[,.(nobs = .N), by=species]
d_nobs$nobs %>% hist


d_nobs <- d_nobs[str_length(species)>0]
d_10 <- d_nobs[order(desc(nobs))][,rank := rank(-nobs)][rank <= 10]
d_20 <- d_nobs[order(desc(nobs))][,rank := rank(-nobs)][rank <= 20]

sala <- sf::st_as_sf(ala[species%in%d_20$species][,.(x,y,species)], 
                     coords=c("x","y")) 
st_crs(sala) <- st_crs(4326)
st_coordinates(sala)

ala[species==d_20$species[1]][,.(x,y)] %>% as.matrix()
s1 <- raster::rasterize(ala[species==d_20$species[1]][,.(x,y)] %>% as.matrix(),
                        grid,
                        fun=function(x,...)length(x))
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
