If you have successfully run the plot creation in parallel, you need to hadd the files together into a form that then gets hadded into "dataset" files and scaled appropriately

For now, I recommend just using the condor queue to do this.
In order to do this, you should run the following
./StopPlotParallelHadderCondor.sh whichNTuple doPURW doSyst

where whichNTuple is 1 or 0 for DESY or Oviedo
doPURW is 1 or 0 for the PURW version or not
and doSyst is 1 or 0 for the systematic booked version or not
