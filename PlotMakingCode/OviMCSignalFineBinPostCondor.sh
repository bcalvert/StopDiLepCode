#! /bin/bash
baseDir=$(cat SkimmingMacro/outputSavePath.txt)
for massPointFile in `/bin/ls SkimmingMacro/SelEffFiles/Tree_FineBin_MassPt_StopMass*_Chi0Mass*_CharginoMass*.txt | sed "s/.txt//"`
  do
  echo ${massPointFile}
  massPoint=$(cat ${massPointFile}.txt )
  StopMass=$(echo ${massPoint} | cut -d\| -f2 | cut -d ' ' -f2)
  Chi0Mass=$(echo ${massPoint} | cut -d\| -f3 | cut -d ' ' -f2)
  CharginoMass=$(echo ${massPoint} | cut -d\| -f4 | cut -d ' ' -f2)
  hadd -f ../RootFiles/Signal/FineBin_Signal_SignalStop${StopMass}_Chi0_${Chi0Mass}_Oviedo_PURW_wSyst_Output.root ${baseDir}/Signal/FineBin/Tree_FineBin_*_Oviedo_SkimOutput_Oviedo_PURW_wSyst_SignalStop${StopMass}_Chi0${Chi0Mass}_Output.root
done

