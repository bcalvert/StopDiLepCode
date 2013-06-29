Ok, so you've gotten files produced, and now you need to hadd and scale those bad boys!
The three key scripts to do so are located in the same directory as this readme (RootFiles in the git dir)
StopPlotHadder.sh
StopPlotHadderCondor.sh
StopPlotHadderPostCondor.sh

The way you utilize these three is, assuming that you've saved the outputs of NewOviStopPlotFillerRunOnSkim_wSyst.C to RootFiles,
1) 
./StopPlotHadderCondor.sh whichNTuple
where whichNTuple is 0 for Oviedo nTuples, or 1 for DESY nTuples
so, for example,
./StopPlotHadderCondor.sh 1
to hadd together the relevant output files for DESY
The way the StopPlotHadderCondor.sh script works is that it first compiles the hadding code (in case you haven't already). It then reads which nTuple to run over from the command line, and then loops over the various input TTBar parameters (essentially which systematic or generator to run over -- this is a bit of a hard code and it's only relevant for the DESY nTuples. I'll be updating as I commision the Oviedo nTuples)

It then also loops over whether or not to do PURW and whether or not to do systematics (so intrinsically, it assumes that you've generated the associated output files for each)

for each of these loop cases, it calls submits condor jobs (using runCondor.sh) of StopPlotHadder.sh for with the various input options
2)
How does StopPlotHadder.sh work? It takes as input 5 command line arguments
which nTuple to run on (0 for Oviedo, 1 for DESY)
which TTBar systematic to run on (documented in the StopPlotHaddInfo.h header file in HeaderFiles)
whether or not to run on pileup reweighted versions of files
whether or not to run on files with systematics booked
whether or not to add in additional verbosity
It then takes these input arguments and calls the compiled version of plotHaddersStop.C (so compile it if you're going to use StopPlotHadder.sh by itself!!!!)
3) 
Once all the condor jobs have finished, you then need to call
./StopPlotHadderPostCondor.sh whichNTuple
where again, whichNTuple is 0 or 1 depending upon which type of nTuple you want to run on

This hadds together the data files and the QCD files into a form that will get called by the plot showing script

#######NOTE#########
6/29/13
For now, do NOT turn off any of the files to hadd in plotHaddersStop.C (i.e. doTTBar or what-not) as I haven't yet made it robust to that, as in right now the way it's set-up it will break unless it runs on all of them...fixing this asap
