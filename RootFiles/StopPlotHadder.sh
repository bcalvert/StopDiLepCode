#! /bin/bash
whichNTuple=$1
whichTTBar=$2
dPURW=$3
dSyst=$4
bVerb=$5

strPURW=""
strSyst=""
strVerb=""
if [ $dPURW -gt 0 ]
then
    strPURW="doPURW"
    echo $strPURW
fi
if [ $dSyst -lt 1 ]
then
    strSyst="noSyst"
    echo $strSyst
fi
if [ $bVerb -gt 0 ]
then
    strVerb="beVerbose"
    echo $strVerb
fi
./plotHaddersStop wNTuple $whichNTuple wTTbarSyst $whichTTBar $strPURW $strSyst $strVerb 

