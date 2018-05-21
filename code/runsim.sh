#!/bin/bash

rscript="main.R"

#experiments=("small-unbiased" "small-biased" "large-unbiased" "large-biased")

#experiments=("small-unbiased" "small-biased")
experiments=("large-unbiased" "large-biased")


#lambdasel="min"
lambdasel="1se"

idjob=1

allidsimul=$(seq 1 100)

#allidsimul=$(seq 21 40)


for idsimul in ${allidsimul[@]}
do
    for experiment in ${experiments[@]}
    do
      flag="$experiment-$idjob-$idsimul-$lambdasel"
      echo "$flag"
      Rscript --vanilla $rscript $experiment $idjob $idsimul $lambdasel  > "/home/rstudio/sparseHTS/rout/results-$flag.Rout" 2> "/home/rstudio/sparseHTS/rout/errorFile-$flag.Rout" &
    done
    if [ `expr $idsimul % 5` -eq 0 ]
    then
      sleep 300
    fi
    

done

