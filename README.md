## **Tentative repository for "*Burn severity and post-fire weather are key to predicting time-to-recover from Australian forest fires*"** 

## Notes:

Repo needs cleaning of deprecated scripts, and remaining scripts are in need of organization.

Parquet files are compressed using the *snappy* codec, which may cause problems for those who have installed *arrow* without *snappy* support.

## Approximate processing order:

\(1\) extract\_\*.R \# needs 128gb ram

\(2\) extract_linear_time_to_recover_SE_coast.R \# \<30 min

\(3\) fit_logistic \*\*\* .R \# needs a few hours

\(4\) fit_ttr5-lai_h2o-xgb.R \# produces TTR prediction for BS fires

## Processing order for predicting recovery time

\(1\) fit_ttr5-lai-ocat_rf.R \# fit predictive model

\(2\) extract_era5-climate-2020.R \# calibrate ERA5-Land to AWAP

\(3\) prep_postBS-data-for-preds.R \# messy merge of data(s)

Notes on script naming convention: 'extract\_*' denotes extracting data from rasters, and writing output in tabular parquet files. Memory heavy. 'functions\_*' helper functions too complex to define in processing scripts. Needs refactor... 'prep\_*' typically uses already extracted data in parquets, and applies pre-processing to prepare the data for model fitting. Memory heavy. 'fit\_*' fits models to preprocessed data and possibly a figure. Should refactor to be consistent w/naming convention. 'plot\_\*' produces one or multiple figures from (mostly) preprocessed data

Notes on defining "time to recover" for vegetation indices (vi):  
vi: [ndvi, kndvi, lai, deltaT, cci, treecover]

vi_anom_12mo \>= vi\_{mean annual} - 0.25\*vi\_{sd}
