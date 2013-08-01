#! /bin/bash
doLocal=$1
grabNTupleVers=$2
baseDir=/data/users/bcalvert/new_oviedo_ntuples_5_14_13/
if [ $grabNTupleVers -gt 0 ]
    then
    baseDir=$(cat SkimmingMacro/outputSavePath.txt)
fi
for mdir in `/bin/ls ${baseDir}*.root | grep DoubleMu | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    if [ $doLocal -gt 0 ]
	then
	./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
    else
	./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
    fi
done
for mdir in `/bin/ls ${baseDIr}*.root | grep DoubleEl | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    if [ $doLocal -gt 0 ]
	then
	./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
    else
	./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
    fi
done
for mdir in `/bin/ls ${baseDIr}*.root | grep MuEG | grep Oviedo_SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    if [ $doLocal -gt 0 ]
	then
	./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
    else
	./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir
    fi
done