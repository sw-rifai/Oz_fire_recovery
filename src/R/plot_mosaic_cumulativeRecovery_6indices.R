library(magick)
library(tidyverse)

i_ndvi <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-ndvi.png")
i_kn <- magick::image_read("figures/figure_cumulativeRecovered_byYear_kndvi-TTR-Def5-ndvi.png")
i_nbr <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-NBR.png")
i_cci <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-cci.png")
i_lsta <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-delta_t.png")
i_tree_cover <- magick::image_read("figures/figure_cumulativeRecovered_byYear_TTR-Def5-tree_cover.png")


p_left <- magick::image_append(c(i_ndvi, i_kn, i_nbr), stack=TRUE)
p_right <- magick::image_append(c(i_cci, i_lsta, i_tree_cover), stack=TRUE)

p_out <- magick::image_append(c(p_left,p_right), stack=F)
p_out
magick::image_write(p_out, path="figures/figure_mosaic_NDVI-kNDVI-NBR-CCI-deltaT-treeCover_cumulativeRecovByYear.png")
