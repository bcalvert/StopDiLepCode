#! /bin/bash
numBPs=$1 #number of break points -- i.e. split the file into $numBPs + 1 different files
numFiles=$[numBPs+1]
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep DoubleMu | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    for ((startPt=1;startPt<=$numFiles;startPt++));
    do
        echo $startPt
        ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doParallel $numBPs $startPt
    done
done
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep DoubleEl | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    for ((startPt=1;startPt<=$numFiles;startPt++));
    do
        echo $startPt
        ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doParallel $numBPs $startPt
    done
done
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep MuEG | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    for ((startPt=1;startPt<=$numFiles;startPt++));
    do
        echo $startPt
        ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doParallel $numBPs $startPt
    done
done