#! /bin/bash
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep -v DoubleMu | grep -v DoubleEl | grep -v MuEG | grep SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
    ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW 
    ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doBookSyst
    ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW doBookSyst
done