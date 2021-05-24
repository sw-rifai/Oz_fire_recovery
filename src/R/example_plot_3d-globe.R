# library(cubs)
library(Ryacas)

Ylm <- function(l,m,phi,theta){
  
  ex <- Ryacas::yac_expr(glue::glue("(-1)^{m}*(1-x^2)^({m}/2)*D(x,{m})OrthoP({l}, x)"))
  
  sqrt((2*l+1)/(4*pi) * factorial(l-m) / factorial(l+m))*
    exp(1i*abs(m)*phi) * eval(ex, list(x = cos(theta)))
  
}


library(threed)
library(rgl) 
library(ggplot2)

N1 <- 90
N2 <- 45
g <- expand.grid(phi=seq(-pi,pi,length=N1),theta=seq(0,pi,length=N2))
g$Y <- Re(Ylm(5,2,phi = g$phi, theta = g$theta))

library(raster)
r <- raster(xmn=-180, ymn=-90, xmx=180, ymx=90, ncols=N1,nrow=N2)
r <- setValues(r, g$Y)

# plot(r)

mesh0 <- quadmesh::quadmesh(r)

## globe versions of longlat
verts <- t(quadmesh::llh2xyz(cbind(mesh0$vb[1, ], mesh0$vb[2,], 0)))
mesh <- mesh3d(vertices = verts, quads = mesh0$ib, meshColor = "faces")
# plot3d(mesh)
# mesh$material$color <- hcl.colors(raster::ncell(r))
# plot3d(mesh)

camera_to_world <- threed::look_at_matrix(eye = c(3, 4, 5)*100, at = c(0, 0, 0))

obj <- mesh %>%
  transform_by(invert_matrix(camera_to_world)) %>%
  orthographic_projection()

fort <- fortify.mesh3d(obj)
fort$value <- rep(g$Y,each=4)

str(fort)
nrow(g)

ggplot(fort) + 
  geom_polygon(aes(x = x, y = y, group = zorder, fill =value), colour = NA, size = 0.2) +
  theme_grey() +
  scale_fill_distiller(palette = 'BrBG') +
  scale_colour_distiller(palette = 'BrBG') +
  theme(
    legend.position = 'none',
    axis.text       = element_blank()
  ) +
  coord_equal() 
