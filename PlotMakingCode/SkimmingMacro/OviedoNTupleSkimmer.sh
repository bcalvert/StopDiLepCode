#! /bin/bash
for mdir in `/bin/ls /data/users/bcalvert/new_oviedo_ntuples_5_14_13/Tree_*root | grep -v Output | sed 's/.root//g'`
  do  
  echo $mdir
  ../runCondor.sh ./NewOviStopPlotSkimmer -i $mdir
done