#!/bin/bash

rscript="newexperiments.R"

exp="large"
addbias=FALSE
algobf="ets"
 

allidsimul=$(seq 1 100)

for idsimul in ${allidsimul[@]}
do
    flag="$exp-$idsimul-$algobf"
    echo "$flag"
    
    experiment="$exp-$idsimul"
    Rscript --vanilla $rscript $experiment $addbias $algobf  > "/home/rstudio/sparseHTS/rout/$flag.Rout" 2> "/home/rstudio/sparseHTS/rout/$flag.err" &
    if [ `expr $idsimul % 5` -eq 0 ]
    then
      sleep 160 #240 #10 
    fi
done

