#! /bin/bash
rm -r EventNumbersMassPoints
mkdir EventNumbersMassPoints
for uniqMassFile in `/bin/ls *_Comp_UniqMassPts* | sed "s/_Comp_UniqMassPts.txt//" | grep FineBin`
  do
  echo ${uniqMassFile}
  numMassPtCount=0
  for massPoint in `cat ${uniqMassFile}_Comp_UniqMassPts.txt`
    do
    numMassPtCount=$[ numMassPtCount + 1 ]
    StopMass=$(echo ${massPoint} | cut -d\: -f1 | cut -d ' ' -f2)
    Chi0Mass=$(echo ${massPoint} | cut -d\: -f2 | cut -d ' ' -f2)
    CharginoMass=$(echo ${massPoint} | cut -d\: -f3 | cut -d ' ' -f2)
    echo \| ${StopMass} GeV \| ${Chi0Mass} GeV \| ${CharginoMass} GeV \| >> EventNumbersMassPoints/${uniqMassFile}_MassPt${numMassPtCount}.txt
#    echo $massPoint
#    echo $numMassPtCount
#    echo $numMassPtCount
#    echo
#    echo "numberofevents at mass point in sortOut${uniqMassFile}_CompSort.txt"
    cat sortOut${uniqMassFile}_CompSort.txt | grep $massPoint | cut -d\: -f1 >> EventNumbersMassPoints/${uniqMassFile}_EventNumber_MassPt${numMassPtCount}.txt
  done
#  cat sortOut${uniqMassFile}.txt
#  cat sortOut${sortFile}.txt | cut -d\: -f2- | sort -u >> ${sortFile}_UniqMassPts.txt
done
