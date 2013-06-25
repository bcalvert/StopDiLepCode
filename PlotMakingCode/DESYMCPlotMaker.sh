#! /bin/bash
for mdir in `/bin/ls /data/users/tkolberg/new_top/*.root | grep -v run | sed 's/.root//g' `
  do  
  echo $mdir
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 doPURW doHackPURW
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 doBookSyst
  ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 doPURW doHackPURW doBookSyst
done