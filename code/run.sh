#!/bin/bash

rscript="new_main.R"

experiment="large"
#lambdasel="min"
lambdasel="1se"
#njobs=5 #allidjobs=$(seq 1 $njobs )
allidjobs=(200)
nbsimul=50

fmethod_agg="AR1"
fmethod_bot="ARIMA"

for idjob in ${allidjobs[@]}
do

flag="$experiment-$idjob-$fmethod_agg-$fmethod_bot-$lambdasel"
    
  Rscript --vanilla $rscript $experiment $idjob $nbsimul $fmethod_agg $fmethod_bot $lambdasel  > "/home/rstudio/sparseHTS/rout/results-$flag.Rout" 2> "/home/rstudio/sparseHTS/rout/errorFile-$flag.Rout" &

done

