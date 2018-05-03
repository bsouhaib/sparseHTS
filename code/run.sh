#!/bin/bash

rscript="new_main.R"

experiments=("small-unbiased" "small-biased" "large-unbiased" "large-biased")
#experiments=("small-unbiased" "small-biased")

#lambdasel="min"
lambdasel="1se"

#allidjobs=$(seq 200 210)
#allidjobs=$(seq 2000 2005)
#allidjobs=(200)
allidjobs=$(seq 1 1)

nbsimul=100

#fmethods_agg=("ETS" "ARIMA" "ETS")
#fmethods_bot=("ARIMA" "ETS" "ETS")

fmethods_agg=("ARIMA")
fmethods_bot=("ARIMA")

for experiment in ${experiments[@]}
do
  for idjob in ${allidjobs[@]}
  do
    for i in "${!fmethods_agg[@]}" 
    do
      fmethod_agg="${fmethods_agg[i]}"
      fmethod_bot="${fmethods_bot[i]}"
      flag="$experiment-$idjob-$fmethod_agg-$fmethod_bot-$lambdasel"
      #Rscript --vanilla $rscript $experiment $idjob $nbsimul $fmethod_agg $fmethod_bot $lambdasel  > "/home/rstudio/sparseHTS/rout/results-$flag.Rout" 2> "/home/rstudio/sparseHTS/rout/errorFile-$flag.Rout" &
      Rscript --vanilla $rscript $experiment $idjob $nbsimul $fmethod_agg $fmethod_bot $lambdasel  > "/Users/souhaibt/Dropbox/Code/hts_sparse/rout/results-$flag.Rout" 2> "/Users/souhaibt/Dropbox/Code/hts_sparse/rout/errorFile-$flag.Rout" &
    done
  done
done

