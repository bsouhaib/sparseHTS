#!/bin/bash

#Rscript newexperiments.R tourism-Bus &
#Rscript newexperiments.R tourism-Vis &
#Rscript newexperiments.R tourism-Hol &
#Rscript newexperiments.R tourism-Oth &
#sleep 10000

#Rscript newexperiments.R road_traffic1 &
#Rscript newexperiments.R road_traffic2 &
#Rscript newexperiments.R road_traffic3 &

#experiment=road_traffic
experiment=wikipedia
#experiment=elec
allidsimul=$(seq 1 5 100)
algobf=ets

#for idsimul in ${allidsimul[@]}
#do
#  Rscript newexperiments.R "$experiment$(($idsimul))" 
#  sleep 10
#done

for idsimul in ${allidsimul[@]}
do
  Rscript newexperiments.R "$experiment-$(($idsimul))" $algobf  &
  Rscript newexperiments.R "$experiment-$(($idsimul+1))" $algobf &
  Rscript newexperiments.R "$experiment-$(($idsimul+2))" $algobf &
  Rscript newexperiments.R "$experiment-$(($idsimul+3))" $algobf &
  Rscript newexperiments.R "$experiment-$(($idsimul+4))" $algobf 
done
