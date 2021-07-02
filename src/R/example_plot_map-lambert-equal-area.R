library(tidyverse)
library(sf)
library(rnaturalearth)
library(stars)

nc <- ne_coastline() %>% sf::st_as_sf()
nc <- ne_countries() %>% sf::st_as_sf() %>% select(admin)
st_crs(nc)
st_crs("EPSG:9834")
st_crs(4326)

st_transform(nc, crs="EPSG:9834")

o <- st_transform(nc,
  crs=st_crs("+proj=laea 
             +lat_0=-30 
             +lon_0=140 
             +x_0=432100 
             +y_0=321000 
             +ellps=GRS80 
             +towgs84=0,0,0,0,0,0,0 
             +units=m 
             +no_defs"))




s <- expand_grid(
  x=seq(140,155,1),
  y=seq(-40,-15,1)) %>% 
  mutate(z=rnorm(416))
ss <- st_as_stars(s, dim=c('x','y'))
st_crs(ss) <- st_crs(4326)
sss <- st_transform_proj(ss,
                  crs=st_crs("+proj=laea +lat_0=-30 +lon_0=140 +x_0=432100 +y_0=321000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
plot(sss)


ggplot()+
  geom_sf(data=o)+
  geom_stars(data=sss, aes(fill=z), 
             alpha=0.25)+
  coord_sf(xlim=c(500000,1750000), 
           ylim=c(-1100000,2300000))+
  theme_linedraw()

ggplot()+
  geom_sf(data=o)+
  geom_stars(data=sss, aes(fill=z), 
             alpha=0.25)+
  coord_sf(xlim=c(500000,1750000), 
           ylim=c(-600000,800000))+
  theme_linedraw()



junk <- st_transform_proj(ss,
                         crs=st_crs("+proj=laea +lat_0=-30 +lon_0=140 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
plot(junk)
