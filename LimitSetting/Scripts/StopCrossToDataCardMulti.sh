#! /bin/bash
#Input Parameters
##Run script by calling ./StopCrossToDataCardMulti.sh doT2TT where doT2TT is 1 for yes (Stop->LSP + Top->b + W->LNu) and 0 is for no (Stop->b + chargino-> L + Snu->LSP + Nu | Nu + Slep -> L + LSP)
doT2TT=$1
WENuBR=.1075
WMuNuBR=.1057
if [ $doT2TT -ge 1 ]; then
    TS='T2TT'
    eeBR=$(echo "scale=8; $WENuBR * $WENuBR" | bc)
    emuBR=$(echo "scale=8; 2 * $WENuBR * $WMuNuBR" | bc) 
    mumuBR=$(echo "scale=8; $WMuNuBR * $WMuNuBR" | bc)
else
    TS=''
    eeBR=0.25 ##Democratic Combinatorics
    emuBR=0.5 ##Democratic Combinatorics
    mumuBR=0.25 ##Democratic Combinatorics
fi

Lumi=19300
BkgNumLumi=5300
LumiFracErr=.01
eeKinEff=0.02
eeSigEff=0.3
eeBkgdNum=33
eeBkgScaled=$(echo "scale=10; $Lumi * $eeBkgdNum / $BkgNumLumi" | bc)

emuKinEff=0.02
emuSigEff=0.75
emuBkgdNum=33
emuBkgScaled=$(echo "scale=10; $Lumi * $emuBkgdNum / $BkgNumLumi" | bc)

mumuKinEff=0.02
mumuSigEff=0.3
mumuBkgdNum=33
mumuBkgScaled=$(echo "scale=10; $Lumi * $mumuBkgdNum / $BkgNumLumi" | bc)
#Input Parameters


MassArray=($(cut -d\| -f2 StopCrossSection.txt | cut -d ' ' -f2)) #GeV
XSecArray=($(cut -d\| -f3 StopCrossSection.txt | cut -d ' ' -f2)) #pb
XSecUncertArray=($(cut -d\| -f3 StopCrossSection.txt | cut -d ' ' -f4 | sed 's/%//'))
b=${#MassArray[*]}
#b=1
rm -rf ../Multi${TS}DataCards
rm -rf ../ee${TS}DataCards
rm -rf ../emu${TS}DataCards
rm -rf ../mumu${TS}DataCards
mkdir ../Multi${TS}DataCards
mkdir ../ee${TS}DataCards
mkdir ../emu${TS}DataCards
mkdir ../mumu${TS}DataCards
sed "s/eeBackNum/$eeBkgScaled/" StopLimitCardTemplate_ee.txt > StopLimitCard${TS}_ee.txt
sed "s/.LumiUncert/$LumiFracErr/g" StopLimitCard${TS}_ee.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_ee.txt
sed "s/emuBackNum/$emuBkgScaled/" StopLimitCardTemplate_emu.txt > StopLimitCard${TS}_emu.txt
sed "s/.LumiUncert/$LumiFracErr/g" StopLimitCard${TS}_emu.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_emu.txt
sed "s/mumuBackNum/$mumuBkgScaled/" StopLimitCardTemplate_mumu.txt > StopLimitCard${TS}_mumu.txt
sed "s/.LumiUncert/$LumiFracErr/g" StopLimitCard${TS}_mumu.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_mumu.txt

sed "s/eeBackNum/$eeBkgScaled/" StopLimitMultiChannelCardTemplate.txt > StopLimitCardMulti${TS}.txt
sed "s/emuBackNum/$emuBkgScaled/" StopLimitCardMulti${TS}.txt > tmp.txt; mv tmp.txt StopLimitCardMulti${TS}.txt
sed "s/mumuBackNum/$mumuBkgScaled/" StopLimitCardMulti${TS}.txt > tmp.txt; mv tmp.txt StopLimitCardMulti${TS}.txt
sed "s/.LumiUncert/$LumiFracErr/g" StopLimitCardMulti${TS}.txt > tmp.txt; mv tmp.txt StopLimitCardMulti${TS}.txt

cp StopLimitCard${TS}_ee.txt ../ee${TS}DataCards
cp StopLimitCard${TS}_emu.txt ../emu${TS}DataCards
cp StopLimitCard${TS}_mumu.txt ../mumu${TS}DataCards
cp StopLimitCardMulti${TS}.txt ../Multi${TS}DataCards
for (( c=0; c<$b; c++ ))
do
    Mass=${MassArray[$c]}
    XSecSciNot=${XSecArray[$c]}
    XSecUncertP=${XSecUncertArray[$c]}
    XSecUncert=$(echo "scale=6; $XSecUncertP/100.0" | bc)
    XSec=$(echo $XSecSciNot |sed 's/e/\*10\^/')

    eeExpEvents=$(echo "scale=12; $XSec*$Lumi*$eeSigEff*$eeKinEff * $eeBR" | bc)
    emuExpEvents=$(echo "scale=12; $XSec*$Lumi*$emuSigEff*$emuKinEff * $emuBR" | bc)
    mumuExpEvents=$(echo "scale=12; $XSec*$Lumi*$mumuSigEff*$mumuKinEff * $mumuBR" | bc)

#    echo $Mass
 #   echo $XSec
  #  echo $XSecUncert
#    echo $ExpEvents
    cd ../Multi${TS}DataCards
    sed "s/eeStopNum/$eeExpEvents/" StopLimitCardMulti${TS}.txt > StopLimitCardMulti${TS}_${Mass}_GeV.txt
    sed "s/emuStopNum/$emuExpEvents/" StopLimitCardMulti${TS}_${Mass}_GeV.txt > tmp.txt; mv tmp.txt  StopLimitCardMulti${TS}_${Mass}_GeV.txt
    sed "s/mumuStopNum/$mumuExpEvents/" StopLimitCardMulti${TS}_${Mass}_GeV.txt > tmp.txt; mv tmp.txt  StopLimitCardMulti${TS}_${Mass}_GeV.txt
    sed "s/STOPMASS/$Mass/" StopLimitCardMulti${TS}_${Mass}_GeV.txt > tmp.txt; mv tmp.txt  StopLimitCardMulti${TS}_${Mass}_GeV.txt
    sed "s/.XSecUncert/$XSecUncert/g" StopLimitCardMulti${TS}_${Mass}_GeV.txt > tmp.txt; mv tmp.txt  StopLimitCardMulti${TS}_${Mass}_GeV.txt
    cd ../ee${TS}DataCards
        sed "s/eeStopNum/$eeExpEvents/" StopLimitCard${TS}_ee.txt > StopLimitCard${TS}_ee_${Mass}_GeV.txt
	sed "s/STOPMASS/$Mass/" StopLimitCard${TS}_ee_${Mass}_GeV.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_ee_${Mass}_GeV.txt
	sed "s/.XSecUncert/$XSecUncert/g" StopLimitCard${TS}_ee_${Mass}_GeV.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_ee_${Mass}_GeV.txt
    cd ../emu${TS}DataCards
    sed "s/emuStopNum/$emuExpEvents/" StopLimitCard${TS}_emu.txt > StopLimitCard${TS}_emu_${Mass}_GeV.txt
        sed "s/STOPMASS/$Mass/" StopLimitCard${TS}_emu_${Mass}_GeV.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_emu_${Mass}_GeV.txt
        sed "s/.XSecUncert/$XSecUncert/g" StopLimitCard${TS}_emu_${Mass}_GeV.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_emu_${Mass}_GeV.txt
    cd ../mumu${TS}DataCards
    sed "s/mumuStopNum/$mumuExpEvents/" StopLimitCard${TS}_mumu.txt > StopLimitCard${TS}_mumu_${Mass}_GeV.txt
        sed "s/STOPMASS/$Mass/" StopLimitCard${TS}_mumu_${Mass}_GeV.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_mumu_${Mass}_GeV.txt
        sed "s/.XSecUncert/$XSecUncert/g" StopLimitCard${TS}_mumu_${Mass}_GeV.txt > tmp.txt; mv tmp.txt StopLimitCard${TS}_mumu_${Mass}_GeV.txt
done

#scp -r Multi${TS}DataCards/ bcalvert@hepcms.umd.edu:~/StopSearch
#scp -r ee${TS}DataCards/ bcalvert@hepcms.umd.edu:~/StopSearch
#scp -r emu${TS}DataCards/ bcalvert@hepcms.umd.edu:~/StopSearch
#scp -r mumu${TS}DataCards/ bcalvert@hepcms.umd.edu:~/StopSearch
