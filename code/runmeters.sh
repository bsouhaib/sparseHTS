#!/bin/bash

rscript="meters_main.R"

hierarchy_names=("UKF23" "UKF16" "UKF21" "UKF13" "UKF15")
for hierarchy_name in ${hierarchy_names[@]}
do
  Rscript --vanilla $rscript  $hierarchy_name > "/home/rstudio/sparseHTS/rout/meters-$hierarchy_name.Rout" 2> "/home/rstudio/sparseHTS/rout/meters-$hierarchy_name.err" &
done

