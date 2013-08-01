#! /bin/bash
numBPs=$1 #number of break points -- i.e. split the file into $numBPs + 1 different files
doLocal=$2
numFiles=$[numBPs+1]
for ((startPt=1;startPt<=$numFiles;startPt++));
  do
  echo $startPt  
  for mdir in `/bin/ls /data/users/tkolberg/new_top/*.root | grep run | sed 's/.root//g'`
    do
    echo $mdir    
    if [ $doLocal -gt 0 ]
	then
	./runCondorLocal.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doParallel $numBPs $startPt
    else
	./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doParallel $numBPs $startPt
    fi
  done
done