#! /bin/bash
rm -rf TrueMassPointsFineBin.txt
ls -lrSh EventNumbersMassPoints/ | grep "[0-9]K" | sed "s/-rw-r--r--.*T/T/" >> TrueMassPointsFineBin.txt
rm -rf SelEffFiles/
mkdir SelEffFiles/
rm -rf Tree_FineBin_EventNumber_Mass*root
rm -rf Tree_FineBin_EventNumber_Mass*Numbers
../../makeRoot.sh SelectionEfficiencyCalc.C
for massPointID in `cat TrueMassPointsFineBin.txt | sed "s/.txt//" | sed "s/Tree_FineBin_EventNumber_MassPt//"`
  do
  echo $massPointID
  massPoint=$(cat EventNumbersMassPoints/Tree_FineBin_MassPt${massPointID}.txt)
  StopMass=$(echo ${massPoint} | cut -d\| -f2 | cut -d ' ' -f2)
  Chi0Mass=$(echo ${massPoint} | cut -d\| -f3 | cut -d ' ' -f2)
  CharginoMass=$(echo ${massPoint} | cut -d\| -f4 | cut -d ' ' -f2)
  echo StopMass $StopMass
  echo Chi0Mass $Chi0Mass
  echo CharginoMass $CharginoMass
  cp EventNumbersMassPoints/Tree_FineBin_MassPt${massPointID}.txt SelEffFiles/Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.txt
  ./SelectionEfficiencyCalc -i EventNumbersMassPoints/Tree_FineBin_EventNumber_MassPt${massPointID} saveRoot >> Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers
  mv Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers SelEffFiles
  mv Tree_FineBin_EventNumber_MassPt${massPointID}.root SelEffFiles
done

