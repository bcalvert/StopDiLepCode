#! /bin/bash
doLocal=$1
grabNTupleVers=$2
baseDir=/data/users/bcalvert/new_oviedo_ntuples_5_14_13/
if [ $grabNTupleVers -gt 0 ]
    then
    baseDir=$(cat SkimmingMacro/outputSavePath.txt)
fi
for mdir in `/bin/ls ${baseDir}*.root | grep -v DoubleMu | grep -v DoubleEl | grep -v MuEG | grep SkimOutput | sed 's/.root//g'`
  do
    echo $mdir    
        if [ $doLocal -gt 0 ]
            then
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doBookSyst
            ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW doBookSyst
        else
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doParallel
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doBookSyst
            ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW
        fi
done