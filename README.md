**code dump/working analysis on recovery from disturbance**

to do :

\(1\) rm prior model fits

\(2\) update raster extracts with stars::st_extract

\(3\) add 2021 fire data?

\(4\) document h2o mod fitting process

\(5\) re-run log fits w/broader param constraints?

Notes: approximate processing order:

\(1\) extract\_\*.R \# needs 128gb ram

\(2\) extract_linear_time_to_recover_SE_coast.R \# \<30 min

\(3\) fit_logistic \*\*\* .R \# needs a few hours

\(4\) separate_pixel_groups.R

(?) attach_climate_to_vi \# merge into later scripts?

### Processing order for predicting recovery time

\(1\) fit_ttr5-lai-ocat_rf.R \# fit predictive model

\(2\) extract_era5-climate-2020.R \# calibrate ERA5-Land to AWAP

\(3\) prep_postBS-data-for-preds.R \# messy merge of data(s)

Notes on script naming convention: 'extract\_*' denotes extracting data from rasters, and writing output in tabular parquet files. Memory heavy. 'functions\_*' helper functions too complex to define in processing scripts. Needs refactor... 'prep\_*' typically uses already extracted data in parquets, and applies pre-processing to prepare the data for model fitting. Memory heavy. 'fit\_*' fits models to preprocessed data and possibly a figure. Should refactor to be consistent w/naming convention. 'plot\_\*' produces one or multiple figures from (mostly) preprocessed data

Notes on defining "time to recover":  
vi: [ndvi, kndvi, lai, deltaT, cci]

Def 1: vi_anom_3mo \>= 0

Def 2: vi_anom \> 0 during JJA

Def 3: vi_anom_12mo \> mandvi\*0.9

Def 4: vi_anom_12mo \>= mandvi

Def 5: vi_anom_12mo \>= mandvi-0.25\*mandvi_sd
