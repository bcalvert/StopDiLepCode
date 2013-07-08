#! /bin/bash
whichNTuple=$1
dPURW=$2
dSyst=$3

strNTup=""
strPURW=""
strSyst=""
strVerb=""

if [ $whichNTuple -gt 0 ]
then
    strNTup="_DESY"
    echo $strNTup
else
    strNTup="_Ovi"
    echo $strNTup
fi
if [ $dPURW -gt 0 ]
then
    strPURW="_PURW"
    echo $strPURW
fi
if [ $dSyst -gt 0 ]
then
    strSyst="_wSyst"
    echo $strSyst
fi
for dataset in `/bin/ls *_DESY_PURW_wSyst_Parallel_1.root | sed "s/_DESY_PURW_wSyst_Parallel_1.root//g"`
#for dataset in `/bin/ls ttbarsignalplustau_scaleup_DESY_wSyst_Parallel_1.root | sed "s/_DESY_wSyst_Parallel_1.root//g"`
  do
  echo $dataset
#  for ((numPar=1;numPar<=$numParFiles;numPar++));
#    do
#    echo $numPar
#  ../PlotMakingCode/
#  ../PlotMakingCode/runCondor.sh StopPlotParallelHadder.sh $dataset $strNTup $strPURW $strSyst
  runCondor.sh ./StopPlotParallelHadder.sh $dataset $strNTup $strPURW $strSyst
#  done
#  for 
#hadd -f $}dataset}_DESY
#ttbarsignalplustau_powheg_DESY_PURW_wSyst_Output.root
done
