#! /bin/bash
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/Tree_*root | grep -v Output | grep -v T2tt | grep -v T2bw | grep -v FineBin | sed 's/.root//g'`
  do  
  echo $mdir
  ../runCondor.sh ./NewOviStopPlotSkimmer -i $mdir gOutDir
done