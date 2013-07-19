#! /bin/bash
numBPs=$1 #number of break points -- i.e. split the file into $numBPs + 1 different files
numFiles=$[numBPs+1]
StopMassArray=($(cut -d\| -f2 StopChi0MassRunVals.txt | cut -d ' ' -f2)) #GeV
Chi0MassArray=($(cut -d\| -f3 StopChi0MassRunVals.txt | cut -d ' ' -f2)) #GeV
numStopMassPts=${#StopMassArray[*]}
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/*.root | grep FineBin | grep SkimOutput | sed 's/.root//g'`
  do
    echo $mdir
    for ((startPt=1;startPt<=$numFiles;startPt++));
    do
        echo $startPt
	for (( c=0; c<$numStopMassPts; c++ ))
	do
            StopMass=${StopMassArray[$c]}
            Chi0Mass=${Chi0MassArray[$c]}
            echo RunStopMass $StopMass
            echo RunChi0Mass $Chi0Mass
            ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir isSig $StopMass $Chi0Mass doParallel $numBPs $startPt
            ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW isSig $StopMass $Chi0Mass doParallel $numBPs $startPt
        ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doBookSyst isSig $StopMass $Chi0Mass doParallel $numBPs $startPt
        ./runCondor.sh ./NewOviStopPlotFillerRunOnSkim_wSyst -i $mdir -w 0 gOutDir doPURW doHackPURW doBookSyst isSig $StopMass $Chi0Mass doParallel $numBPs $startPt	
	done
    done
done