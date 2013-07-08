#! /bin/bash
dataset=$1
strNTup=$2
strPURW=$3
strSyst=$4
echo $dataset
echo $strNTup
echo $strPURW
echo $strSyst
hadd -f ${dataset}${strNTup}${strPURW}${strSyst}_Output.root ${dataset}${strNTup}${strPURW}${strSyst}_Parallel_*.root

