#! /bin/bash
for outFile in `/bin/ls -a ._NewOviStopPlotSkimmer_-i_hadoop_store_user_bcalvert_OviedoStopSignal_*_gOutDir_isSig_pEvNum_*.sh.stdout`
  do
  name=$(echo $outFile | sed "s/\(_gOutDir_isSig_pEvNum_[0-9]*.sh.stdout\)//" | sed "s/._NewOviStopPlotSkimmer_-i_hadoop_store_user_bcalvert_OviedoStopSignal_//")
  echo $name
  rm -rf notSortOut${name}.txt
  rm -rf sortOut${name}.txt
  cat .*${name}_*gOutDir_isSig_pEvNum*stdout| grep "in format EventNum:genStopMass:genChi0" | sed "s/in format EventNum:genStopMass:genChi0Mass:genCharginoMass //g" >> notSortOut${name}.txt
  cat notSortOut${name}.txt | sort -s -n -k 1,1 >> sortOut${name}.txt
done