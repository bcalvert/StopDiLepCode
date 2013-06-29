#! /bin/bash

for whichTTBar in 0 1 2 3 4 5 6 7 8
#for whichTTBar in 1 2 3 4 5 6 7 8
do
  runCondor.sh /data/users/bcalvert/new_top/DESYPlotHadder.sh ${whichTTBar}
done

#hadd -f DataEE_DESYHaddplots.root ee_run2012*MT2Leq80_DESY_Output.root
#hadd -f DataEMu_DESYHaddplots.root emu_run2012*MT2Leq80_DESY_Output.root
#hadd -f DataMuMu_DESYHaddplots.root mumu_run2012*MT2Leq80_DESY_Output.root
#hadd -f QCD_DESY_PURWHaddplots.root QCDMu_DESY_PURWHaddplots.root QCDEM_DESY_PURWHaddplots.root QCDBCEM_DESY_PURWHaddplots.root 
#hadd -f QCD_DESYHaddplots.root QCDMu_DESYHaddplots.root QCDEM_DESYHaddplots.root QCDBCEM_DESYHaddplots.root 
#hadd -f QCD_DESY_PURW_wSystHaddplots.root QCDMu_DESY_PURW_wSystHaddplots.root QCDEM_DESY_PURW_wSystHaddplots.root QCDBCEM_DESY_PURW_wSystHaddplots.root 
#hadd -f QCD_DESY_wSystHaddplots.root QCDMu_DESY_wSystHaddplots.root QCDEM_DESY_wSystHaddplots.root QCDBCEM_DESY_wSystHaddplots.root 



#hadd -f TTBar
#  LocalWords:  hadd
