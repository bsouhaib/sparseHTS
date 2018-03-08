#!/bin/bash

rscript="main.R"
njobs=5

towardspbu=TRUE
#towardspbu=FALSE

lambdasel="min"
#lambdasel="1se"

nbsimul=10


allidjobs=$(seq 1 $njobs )

for idjob in ${allidjobs[@]}
do

flag="$idjob-$towardspbu-$lambdasel"
    
  Rscript --vanilla $rscript $idjob $lambdasel $towardspbu $nbsimul > "/home/rstudio/sparseHTS/rout/results-$flag.Rout" 2> "/home/rstudio/sparseHTS/rout/errorFile-$flag.Rout" &

done

