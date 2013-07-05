#! /bin/bash
#Input Parameters
Lumi=19300
LumiFracErr=.01
SigEff=0.5
KinEff=0.02
BkgNumLumi=5300
BkgdNum=100
BkgScaled=$(echo "scale=10; $Lumi * $BkgdNum / $BkgNumLumi" | bc)
#Input Parameters


MassArray=($(cut -d\| -f2 StopCrossSection.txt | cut -d ' ' -f2)) #GeV
XSecArray=($(cut -d\| -f3 StopCrossSection.txt | cut -d ' ' -f2)) #pb
XSecUncertArray=($(cut -d\| -f3 StopCrossSection.txt | cut -d ' ' -f4 | sed 's/%//'))
b=${#MassArray[*]}
#b=1
rm -rf ../DataCards
mkdir ../DataCards
sed "s/BackNum/$BkgScaled/" StopLimitCardTemplateSimple.txt > StopLimitCardSimple.txt
sed "s/.LumiUncert/$LumiFracErr/" StopLimitCardSimple.txt > tmp.txt; mv tmp.txt StopLimitCardSimple.txt
cp StopLimitCardSimple.txt ../DataCards
cd ../DataCards
for (( c=0; c<$b; c++ ))
do
    Mass=${MassArray[$c]}
    XSecSciNot=${XSecArray[$c]}
    XSecUncertP=${XSecUncertArray[$c]}
    XSecUncert=$(echo "scale=6; $XSecUncertP/100.0" | bc)
    XSec=$(echo $XSecSciNot |sed 's/e/\*10\^/')

    ExpEvents=$(echo "scale=12; $XSec*$Lumi*$SigEff*$KinEff" | bc)
#    echo $Mass
 #   echo $XSec
  #  echo $XSecUncert
#    echo $ExpEvents
    sed "s/StopNum/$ExpEvents/" StopLimitCardSimple.txt > StopLimitCardSimple_${Mass}_GeV.txt
    sed "s/STOPMASS/$Mass/" StopLimitCardSimple_${Mass}_GeV.txt > tmp.txt; mv tmp.txt  StopLimitCardSimple_${Mass}_GeV.txt
    sed "s/.XSecUncert/$XSecUncert/" StopLimitCardSimple_${Mass}_GeV.txt > tmp.txt; mv tmp.txt  StopLimitCardSimple_${Mass}_GeV.txt
done
cd ..
#scp -r ../DataCards/ bcalvert@hepcms.umd.edu:~/StopSearch
