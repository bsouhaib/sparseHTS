#!/bin/bash

rscript="main.R"

#experiments=("small-unbiased" "small-biased" "large-unbiased" "large-biased")
#experiments=("large-unbiased")
experiments=("large-biased")

#lambdasel="min"
lambdasel="1se"

allidjobs=$(seq 1 1)

#allidsimul=$(seq 1 20)
allidsimul=$(seq 21 40)

for experiment in ${experiments[@]}
do
  for idjob in ${allidjobs[@]}
  do
    for idsimul in ${allidsimul[@]}
    do
      flag="$experiment-$idjob-$idsimul-$lambdasel"
      Rscript --vanilla $rscript $experiment $idjob $idsimul $lambdasel  > "/home/rstudio/sparseHTS/rout/results-$flag.Rout" 2> "/home/rstudio/sparseHTS/rout/errorFile-$flag.Rout" &
    done
  done
done

