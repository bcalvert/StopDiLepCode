#! /bin/bash
cd ~/CMSSW525HiggsLimit/
cmsenv
cd -
doT2TT=$1
whichChan=0
doSyst=0
systString=NoSyst
if [ $2 -ge 1 ]; then
    whichChan=$2
fi
if [ $3 -ge 1 ]; then
    doSyst=1
    systString=WSyst
fi
case $whichChan in
  0)
     chanString=Multi
     ;;
  1)
     chanString=_ee
     ;;
  2)
     chanString=_emu
     ;;
  3)
     chanString=_mumu
     ;;
esac
echo $chanString
if [ $doT2TT -ge 1 ]; then
    TS='T2TT'
else
    TS=''
fi
MassArray=($(cut -d\| -f2 StopCrossSection.txt | cut -d ' ' -f2)) #GeV
Method=Asymptotic
rm -rf ../${chanString}${TS}RootLimitFiles/
mkdir ../${chanString}${TS}RootLimitFiles/
for i in ${MassArray[@]}
#for i in 1005 #uncomment this to test a single mass point (comment "for i in ...." line above if you want to do this)
  do 
  echo mass $i
#  combine -S $doSyst -M HybridNew --frequentist ~/StopSearch/DataCards/StopLimitCardSimple_${i}_GeV.txt
  combine -n ${chanString}${TS}${systString} -m $i -S $doSyst -M $Method --run expected ~/StopSearch/${chanString}${TS}DataCards/StopLimitCard${chanString}${TS}_${i}_GeV.txt
  mv  higgsCombine${chanString}${TS}${systString}.$Method.mH${i}.root ../${chanString}${TS}RootLimitFiles/StopLimits${chanString}${TS}${systString}.$Method.mStop_${i}_GeV.root
done
hadd -f ../RootLimitFiles/StopLimitComposite${chanString}${TS}${Method}${systString}.root ../${chanString}${TS}RootLimitFiles/StopLimits${chanString}~*.root