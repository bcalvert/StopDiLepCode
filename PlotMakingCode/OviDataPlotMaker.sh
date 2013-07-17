#! /bin/bash
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep DoubleMu | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
done
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep DoubleEl | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
done
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep MuEG | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
done