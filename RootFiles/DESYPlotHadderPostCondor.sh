#! /bin/bash
hadd -f DataEE_DESYHaddplots.root ee_run2012*MT2Leq80_DESY_Output.root
hadd -f DataEMu_DESYHaddplots.root emu_run2012*MT2Leq80_DESY_Output.root
hadd -f DataMuMu_DESYHaddplots.root mumu_run2012*MT2Leq80_DESY_Output.root
hadd -f QCD_DESY_PURWHaddplots.root QCDMu_DESY_PURWHaddplots.root QCDEM_DESY_PURWHaddplots.root QCDBCEM_DESY_PURWHaddplots.root 
hadd -f QCD_DESYHaddplots.root QCDMu_DESYHaddplots.root QCDEM_DESYHaddplots.root QCDBCEM_DESYHaddplots.root 
hadd -f QCD_DESY_PURW_wSystHaddplots.root QCDMu_DESY_PURW_wSystHaddplots.root QCDEM_DESY_PURW_wSystHaddplots.root QCDBCEM_DESY_PURW_wSystHaddplots.root 
hadd -f QCD_DESY_wSystHaddplots.root QCDMu_DESY_wSystHaddplots.root QCDEM_DESY_wSystHaddplots.root QCDBCEM_DESY_wSystHaddplots.root 