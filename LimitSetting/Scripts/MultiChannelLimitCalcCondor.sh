#! /bin/bash
doT2TT=$1
doSyst=$2
for chan in 0 1 2 3
  do
  runCondor.sh ./MultiChannelLimitCalc.sh $doT2TT $chan $doSyst
done