#! /bin/bash
cd ~/CMSSW525HiggsLimit/
cmsenv
cd -
MassArray=($(cut -d\| -f2 StopCrossSection.txt | cut -d ' ' -f2)) #GeV
Method=Asymptotic
for i in ${MassArray[@]}

#for i in 1005
  do 
  echo mass $i
#  combine -S 0 -M HybridNew --frequentist ~/StopSearch/DataCards/StopLimitCardSimple_${i}_GeV.txt
  combine -n SimpleNoSyst -m $i -S 0 -M $Method --run expected ~/StopSearch/DataCards/StopLimitCardSimple_${i}_GeV.txt
  mv  higgsCombineSimpleNoSyst.$Method.mH${i}.root StopLimitsSimpleNoSyst.$Method.mStop_${i}_GeV.root
done
rm -rf ../RootLimitFiles
mkdir -p ../RootLimitFiles/RootLimitIndMassPoints
mv StopLimitsSimple*.root ../RootLimitFiles/RootLimitIndMassPoints/

cd ../RootLimitFiles/
hadd -f StopLimitCompositeSimple${Method}NoSyst.root /RootLimitIndMassPoints/*.root 