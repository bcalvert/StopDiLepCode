#! /bin/bash
for mdir in `/bin/ls /data/users/tkolberg/new_top/*.root | grep -v run | sed 's/.root//g' `
  do  
  echo $mdir
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doPURW doHackPURW
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doBookSyst
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doPURW doHackPURW doBookSyst
done