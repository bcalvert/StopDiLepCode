#! /bin/bash
../makeRoot.sh plotHaddersStop.C
whichNTuple=$1
for ((whichTTBar=0;whichTTBar<=8;whichTTBar++));
do
    for PURW in 0 1 
    do
	for doSyst in 0 1
	do
	    ../PlotMakingCode/runCondor.sh ./StopPlotHadder.sh $whichNTuple $whichTTBar $PURW $doSyst 0
	done
    done
done