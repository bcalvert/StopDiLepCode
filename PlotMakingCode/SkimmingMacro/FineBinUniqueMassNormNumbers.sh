#! /bin/bash
rm -rf ../FineBinNormFiles
mkdir ../FineBinNormFiles
for massPointID in `cat TrueMassPointsFineBin.txt | sed "s/.txt//" | sed "s/Tree_FineBin_EventNumber_MassPt//"`
  do
  massPoint=$(cat EventNumbersMassPoints/Tree_FineBin_MassPt${massPointID}.txt)
  StopMass=$(echo ${massPoint} | cut -d\| -f2 | cut -d ' ' -f2)
  Chi0Mass=$(echo ${massPoint} | cut -d\| -f3 | cut -d ' ' -f2)
  CharginoMass=$(echo ${massPoint} | cut -d\| -f4 | cut -d ' ' -f2)
  numSurviveSkim=$(cat EventNumbersMassPoints/Tree_FineBin_EventNumber_MassPt${massPointID}.txt | wc -l)
  skimEff1=$(cat SelEffFiles/Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers | grep 1/Mean | sed "s/\(.*\)\([0-9]*\.[0-9]*\) pm \([0-9]*\.[0-9]*\)/0\2/")
  skimEff2=$(cat SelEffFiles/Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers | grep exp | sed "s/\(.*\)\([0-9]*\.[0-9]*\) pm \([0-9]*\.[0-9]*\)/0\2/")
  skimEffErr1=$(cat SelEffFiles/Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers | grep 1/Mean | sed "s/\(.*\)\([0-9]*\.[0-9]*\) pm \([0-9]*\.[0-9]*\)/\3/")
  skimEffErr2=$(cat SelEffFiles/Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers | grep exp | sed "s/\(.*\)\([0-9]*\.[0-9]*\) pm \([0-9]*\.[0-9]*\)/\3/")
  skimEffGeomAve=$(echo "scale=6; sqrt(${skimEff1} * ${skimEff2})" | bc)
  skimEffGeomAveErr=$(echo "scale=6; (1/sqrt(2)) * sqrt((${skimEff1}/${skimEff2}) * (${skimEffErr2} * ${skimEffErr2}) + (${skimEff2}/${skimEff1}) * (${skimEffErr1} * ${skimEffErr1}))" | bc)
  origNum=$(echo "scale=6; ${numSurviveSkim} * (1/${skimEffGeomAve})" | bc)
  origNumErr=$(echo "scale=6; ${numSurviveSkim} * (${skimEffGeomAveErr}/(${skimEffGeomAve} * ${skimEffGeomAve}))" | bc)
#  echo massPointID $massPointID
#  echo StopMass $StopMass
#  echo Chi0Mass $Chi0Mass
#  echo CharginoMass $CharginoMass
#  echo numSurviveSkim ${numSurviveSkim}
#  echo skimEff1 ${skimEff1}
#  echo skimEffErr1 ${skimEffErr1}
#  echo skimEff2 ${skimEff2}
#  echo skimEffErr2 ${skimEffErr2}
  echo skimEffGeomAveEstimate ${skimEffGeomAve} >> ../FineBinNormFiles/NormNumbers_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.txt
  echo skimEffGeomAveEstimateErr ${skimEffGeomAveErr} >> ../FineBinNormFiles/NormNumbers_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.txt
  echo origNumEstimate ${origNum} >> ../FineBinNormFiles/NormNumbers_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.txt
  echo origNumEstimateErr ${origNumErr} >> ../FineBinNormFiles/NormNumbers_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.txt

#  cp EventNumbersMassPoints/Tree_FineBin_MassPt${massPointID}.txt SelEffFiles/Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.txt
#  ./SelectionEfficiencyCalc -i EventNumbersMassPoints/Tree_FineBin_EventNumber_MassPt${massPointID} saveRoot >> Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers
#  mv Tree_FineBin_MassPt_StopMass${StopMass}_Chi0Mass${Chi0Mass}_CharginoMass${CharginoMass}.Numbers SelEffFiles
#  mv Tree_FineBin_EventNumber_MassPt${massPointID}.root SelEffFiles
done

