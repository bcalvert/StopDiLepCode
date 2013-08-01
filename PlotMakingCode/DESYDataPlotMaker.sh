#! /bin/bash
doLocal=$1
for mdir in `/bin/ls /data/users/tkolberg/new_top/*.root | grep run | sed 's/.root//g' `
  do  
  echo $mdir
  if [ $doLocal -gt 0 ]
      then
      ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir      
  else 
      ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir
  fi
done