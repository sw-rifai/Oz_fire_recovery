**Testing ideas on recovery from disturbance**

Main Qs:

1.  How does historical 'greenness' affect time to recover after fire?

2.  Can a historical positive greenness anomaly increase probability of fire, or the size of ndvi reduction if there is fire? \~2012 fires post 2010/11 La Niña

3.  Can we predict the time to recover from the Black Summer fires with only one year of recovery data? nonlinear stat deep dive...

4.  If possible w/historical Landsat record. Has increased CO2 accelerated the time to recover from fire?



Notes: approximate processing order:

(1) preprocess\_data.R \# needs 128gb ram

(2) extract\_linear\_time\_to\_recover\_SE\_coast.R

(3) fit\_weibull\_ttr.R \# needs a few hours

(4) separate\_pixel\_groups.R

(?) attach\_climate\_to\_vi \# merge into later scripts?
