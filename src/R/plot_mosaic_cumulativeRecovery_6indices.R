library(magick)
library(tidyverse)

i_p <- magick::image_read("figures/figure_tile_mean-frac-precip-12mo-anom.png")
i_ppet <- magick::image_read("figures/figure_tile_mean-frac-ppet-12mo-anom.png")
i_vpd <- magick::image_read("figures/figure_tile_mean-frac-vpd15-12mo-anom.png")

i_ndvi <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-ndvi.png")
i_kn <- magick::image_read("figures/figure_cumulativeRecovered_byYear_kndvi-TTR-Def5-ndvi.png")
i_nbr <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-NBR.png")
i_cci <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-cci.png")
i_lsta <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-delta_t.png")
i_tree_cover <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-tree_cover.png")
i_lai <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-lai.png")

p_top <- magick::image_append(c(i_p, i_ppet,i_vpd),stack=F)
p_mid <- magick::image_append(c(i_ndvi, i_lai, i_lsta), stack=F)
p_bottom <- magick::image_append(c(i_kn, i_nbr, i_tree_cover), stack=F)

p_out <- magick::image_append(c(p_top,p_mid,p_bottom), stack=T)
p_out
magick::image_write(p_out, path="figures/figure_mosaic_3-clim_6-vegIndices_cumulativeRecovByYear.png")
