#! /bin/bash
doLocal=$1
for mdir in `/bin/ls /data/users/tkolberg/new_top/*.root | grep -v run | sed 's/.root//g'`
  do  
  echo $mdir
        if [ $doLocal -gt 0 ]
            then
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doPURW doHackPURW
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doBookSyst
            ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doPURW doHackPURW doBookSyst
        else
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doParallel
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doPURW doHackPURW
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doBookSyst
            ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doPURW doHackPURW
        fi
done