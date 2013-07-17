ROOTLIBS=`root-config --ldflags --libs` -Wl,--no-as-needed -lGenVector
COMPFLAGS=`root-config --cflags` -Wall -Wextra -pedantic -Wno-vla -O2 -g -std=c++0x -DSTANDALONE_COMPILE

all: NewOviStopPlotSkimmer

#StopPlotFiller:  StopDict_ver2.o StopPlotSkimmer.o
NewOviStopPlotSkimmer:  StopDict.o NewOviStopPlotSkimmer.o
	g++ -o $@ $^ $(ROOTLIBS)

#StopDict_ver2.o: StopDict_ver2.cpp
StopDict.o: StopDict.cpp
	g++ -c -o $@ $(COMPFLAGS) $^

NewOviStopPlotSkimmer.o: NewOviStopPlotSkimmer.C
	g++ -c -o $@ $(COMPFLAGS) $^