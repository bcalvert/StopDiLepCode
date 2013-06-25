#! /bin/bash

for mdir in `/bin/ls /data/users/tkolberg/new_top/*.root | grep run | sed 's/.root//g' `
  do  
  echo $mdir
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir
done