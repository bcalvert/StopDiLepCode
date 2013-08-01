#! /bin/bash
../makeRoot.sh plotHaddersStop.C
whichNTuple=$1
maxTTBarGen=8
if [ $whichNTuples -gt 0 ]
    then
    maxTTBarGen=2
fi

for ((whichTTBar=0;whichTTBar<=$maxTTBarGen;whichTTBar++));
do
    for PURW in 0 1 
    do
	for doSyst in 0 1
	do
	    ../PlotMakingCode/runCondor.sh ./StopPlotHadder.sh $whichNTuple $whichTTBar $PURW $doSyst 0
#	    ./StopPlotHadder.sh $whichNTuple $whichTTBar $PURW $doSyst 0
	done
    done
done