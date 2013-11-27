#! /bin/bash
StopMass=$1
Chi0Mass=$2
CharginoMassFrac=$3
MT2llCut=$4
MT2lblbCut=$5
doT2tt=$6
SignalString=""
if [ $doT2tt -gt 0 ]
then
    SignalString="T2tt T2tt_FineBin"
    rm -rf ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt
    rm -rf ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}.txt
else
    SignalString="T2bw T2bw_FineBin"
    rm -rf ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt
    rm -rf ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt
fi
echo "StopMass ${StopMass} Chi0Mass ${Chi0Mass}" >> ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}.txt
for FOMType in 0 1
do
#    echo $StopMass
#    echo $Chi0Mass
#    echo $MT2llCut
#    echo $MT2lblbCut
    if [ $doT2tt -gt 0 ]
    then
	./MT2CutYieldCalc -b -q wChan 4 doPURW wTTbarGen 2 wNTuple 0 doReReco useDDEst doFOM ${FOMType} doOneDeeFOM doSignal T2tt T2tt_FineBin ${StopMass} ${Chi0Mass} ${CharginoMassFrac} noPlots pFOMI MT2AxisCuts ${MT2llCut} ${MT2lblbCut} JsSm 0 >> ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt
	cat ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt | grep -A3 "Printing FOM Info" >> ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}.txt
	rm -rf ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt
#    ./MT2CutYieldCalc -b -q wChan 4 doPURW wTTbarGen 2 wNTuple 0 doReReco useDDEst doFOM ${FOMType} doOneDeeFOM doSignal T2tt T2tt_FineBin ${StopMass} ${Chi0Mass} ${CharginoMassFrac} noPlots pFOMI MT2AxisCuts ${MT2llCut} ${MT2lblbCut} >> ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt
#    cat ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt | grep -A3 "Printing FOM Info" >> ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}.txt
#    rm -rf ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt
	./MT2CutYieldCalc -b -q wChan 4 doPURW wTTbarGen 2 wNTuple 0 doReReco useDDEst doFOM ${FOMType} doTwoDeeFOM doSignal T2tt T2tt_FineBin ${StopMass} ${Chi0Mass} ${CharginoMassFrac} noPlots pFOMI MT2AxisCuts ${MT2llCut} ${MT2lblbCut} JsSm 0 >> ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt
	cat ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt | grep -A4 "Printing FOM Info" >> ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}.txt
	rm -rf ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt
    else
	./MT2CutYieldCalc -b -q wChan 4 doPURW wTTbarGen 2 wNTuple 0 doReReco useDDEst doFOM ${FOMType} doOneDeeFOM doSignal T2tt T2tt_FineBin ${StopMass} ${Chi0Mass} ${CharginoMassFrac} noPlots pFOMI MT2AxisCuts ${MT2llCut} ${MT2lblbCut} JsSm 0 >> ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt
	cat ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}.txt | grep -A3 "Printing FOM Info" >> ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt
	rm -rf ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt


	./MT2CutYieldCalc -b -q wChan 4 doPURW wTTbarGen 2 wNTuple 0 doReReco useDDEst doFOM ${FOMType} doTwoDeeFOM doSignal T2tt T2tt_FineBin ${StopMass} ${Chi0Mass} ${CharginoMassFrac} noPlots pFOMI MT2AxisCuts ${MT2llCut} ${MT2lblbCut} JsSm 0 >> ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt
        cat ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt | grep -A4 "Printing FOM Info" >> ~/Cull_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt
        rm -rf ~/Temp_Stop_${StopMass}_LSP_${Chi0Mass}_Chargino_${CharginoMassFrac}.txt
    fi
done