#!/bin/bash

rscript="new_main.R"

experiment="large"
#lambdasel="min"
lambdasel="1se"
#njobs=5 
#allidjobs=$(seq 200 210)
allidjobs=$(seq 2000 2060)
#allidjobs=(200)
nbsimul=10

#fmethods_agg=("ETS" "ARIMA" "ETS")
#fmethods_bot=("ARIMA" "ETS" "ETS")

fmethods_agg=("ARIMA")
fmethods_bot=("ARIMA")

for idjob in ${allidjobs[@]}
do
  for i in "${!fmethods_agg[@]}" 
  do
    fmethod_agg="${fmethods_agg[i]}"
    fmethod_bot="${fmethods_bot[i]}"
    flag="$experiment-$idjob-$fmethod_agg-$fmethod_bot-$lambdasel"
    Rscript --vanilla $rscript $experiment $idjob $nbsimul $fmethod_agg $fmethod_bot $lambdasel  > "/home/rstudio/sparseHTS/rout/results-$flag.Rout" 2> "/home/rstudio/sparseHTS/rout/errorFile-$flag.Rout" &
  done
done

