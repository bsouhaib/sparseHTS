#!/bin/bash


Rscript newresults.R large  all 100
Rscript newresults.R small  all 100
Rscript newresults.R wikipedia all 50
Rscript newresults.R road_traffic all 50

exit

Rscript newresults.R large  arima 100
Rscript newresults.R large  ets 100

Rscript newresults.R small  arima 100
Rscript newresults.R small  ets   100

Rscript newresults.R wikipedia arima 50
Rscript newresults.R wikipedia  ets 40

Rscript newresults.R road_traffic arima 50
Rscript newresults.R road_traffic  ets 50