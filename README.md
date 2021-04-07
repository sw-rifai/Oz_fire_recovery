**working analysis on recovery from disturbance**

Main Qs:

1.  How does historical 'greenness' affect time to recover after fire?

2.  Can a historical positive greenness anomaly increase probability of fire, or the size of ndvi reduction if there is fire? \~2012 fires post 2010/11 La Niña

3.  Can we predict the time to recover from the Black Summer fires with only one year of recovery data? nonlinear stat deep dive...

4.  If possible w/historical Landsat record. Has increased CO2 accelerated the time to recover from fire?



Notes: approximate processing order:

(1) extract_*.R \# needs 128gb ram

(2) extract\_linear\_time\_to\_recover\_SE\_coast.R

(3) fit\_weibull\_ttr.R \# needs a few hours

(4) separate\_pixel\_groups.R

(?) attach\_climate\_to\_vi \# merge into later scripts?


Notes on script naming convention:
'extract_*' denotes extracting data from rasters, and writing output in tabular parquet files. Memory heavy. 
'functions_*' helper functions too complex to define in processing scripts. Needs refactor...
'prep_*' typically uses already extracted data in parquets, and applies pre-processing to prepare the data for model fitting. Memory heavy.
'fit_*' fits models to preprocessed data and possibly a figure. Should refactor to be
consistent w/naming convention.
'plot_*' produces one or multiple figures from preprocessed data  


Notes on defining "time to recover": 
* Def 1: ndvi_anom_3mo >= 0
* Def 2: ndvi_anom > 0 during JJA
* Def 3: ndvi_anom_12mo > mandvi*0.9
* Def 4: ndvi_anom_12mo >= mandvi
