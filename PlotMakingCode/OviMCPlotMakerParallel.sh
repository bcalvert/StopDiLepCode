#! /bin/bash
numBPs=$1 #number of break points -- i.e. split the file into $numBPs + 1 different files
numFiles=$[numBPs+1]
doLocal=$2
grabNTupleVers=$3
baseDir=/data/users/bcalvert/new_oviedo_ntuples_5_14_13/
if [ $grabNTupleVers -gt 0 ]
    then
    baseDir=$(cat SkimmingMacro/outputSavePath.txt)
fi
for ((startPt=1;startPt<=$numFiles;startPt++));
  do
  echo $startPt
  for mdir in `/bin/ls ${baseDir}*.root | grep -v DoubleMu | grep -v DoubleEl | grep -v MuEG | grep SkimOutput | grep -v FineBin | grep -v T2tt | grep -v T2bw | sed 's/.root//g'`
    do
    echo $mdir
    
    if [ $doLocal -gt 0 ]
	then
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doParallel $numBPs $startPt
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW  doParallel $numBPs $startPt
#           ./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doBookSyst doParallel $numBPs $startPt
	./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW doBookSyst doParallel $numBPs $startPt
    else
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doParallel $numBPs $startPt
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW  doParallel $numBPs $startPt
#           ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doBookSyst doParallel $numBPs $startPt
	./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW doBookSyst doParallel $numBPs $startPt
    fi
  done
done