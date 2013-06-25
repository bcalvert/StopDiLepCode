ROOTLIBS=`root-config --ldflags --libs` -Wl,--no-as-needed -lGenVector
COMPFLAGS=`root-config --cflags` -Wall -Wextra -pedantic -Wno-vla -O2 -g -std=c++0x -DSTANDALONE_COMPILE
#COMPFLAGS=`root-config --cflags` -Wall -Wextra -pedantic -O2 -g -DSTANDALONE_COMPILE

all: NewOviStopPlotFillerRunOnSkim_wSyst

#StopPlotFillerRunOnSkim_wSyst:  StopDict_ver2.o StopPlotFillerRunOnSkim_wSyst.o
NewOviStopPlotFillerRunOnSkim_wSyst:  StopDict.o NewOviStopPlotFillerRunOnSkim_wSyst.o
	g++ -o $@ $^ $(ROOTLIBS) mt2bisect.cc

#StopDict_ver2.o: StopDict_ver2.cpp
StopDict.o: StopDict.cpp
	g++ -c -o $@ $(COMPFLAGS) $^

NewOviStopPlotFillerRunOnSkim_wSyst.o: NewOviStopPlotFillerRunOnSkim_wSyst.C
	g++ -c -o $@ $(COMPFLAGS) $^

NewOviStopPlotFillerRunOnSkim_wSyst_wSyst:  StopDict.o NewOviStopPlotFillerRunOnSkim_wSyst_wSyst.o
	g++ -o $@ $^ $(ROOTLIBS) mt2bisect.cc

NewOviStopPlotFillerRunOnSkim_wSyst_wSyst.o: NewOviStopPlotFillerRunOnSkim_wSyst_wSyst.C
	g++ -c -o $@ $(COMPFLAGS) $^