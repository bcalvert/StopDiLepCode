#! /bin/bash
numBPs=$1 #number of break points -- i.e. split the file into $numBPs + 1 different files
numFiles=$[numBPs+1]
for mdir in `/bin/ls /data/users/tkolberg/new_top/*.root | grep run | sed 's/.root//g'`
  do
    echo $mdir
    for ((startPt=1;startPt<=$numFiles;startPt++));
    do
        echo $startPt
        ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 1 gOutDir doParallel $numBPs $startPt
    done
done