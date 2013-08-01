#! /bin/bash
whichNTuple=$1
strNTup="Ovi"
if [ $whichNTuple -gt 0 ]
then
    strPURW="DESY"

    hadd -f DataEE_${strPURW}Haddplots.root ee_run2012*MT2Leq80_${strPURW}_Output.root
    hadd -f DataEMu_${strPURW}Haddplots.root emu_run2012*MT2Leq80_${strPURW}_Output.root
    hadd -f DataMuMu_${strPURW}Haddplots.root mumu_run2012*MT2Leq80_${strPURW}_Output.root
    for a in ${strPURW}_PURW ${strPURW}
      do
      for b in _wSystHadd Hadd
	do
        hadd -f QCD_${a}${b}plots.root QCDMu_${a}${b}plots.root QCDEM_${a}${b}plots.root QCDBCEM_${a}${b}plots.root
      done
    done
else
    strPURW="Oviedo"
    hadd -f DataEE_${strPURW}Haddplots.root DoubleEl*_${strPURW}_Output.root
    hadd -f DataMuMu_${strPURW}Haddplots.root DoubleMuMu*_${strPURW}_Output.root
    hadd -f DataEMu_${strPURW}Haddplots.root MuEG*_${strPURW}_Output.root
fi