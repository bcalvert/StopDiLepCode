#include <iostream>
#include <vector>
#include "../HeaderFiles/StopFunctionDefinitions_v2.h"
#include "TROOT.h"
#include "TRint.h"

using namespace std;
int main( int argc, char* argv[]) {
    gROOT->ProcessLine("#include <vector>");
    TRint theApp("App", &argc, argv);
    Bool_t retVal = kTRUE;
    vector<SampleT> * subSampVec    = SubSampVec();
    cout << "subsamp size " << subSampVec->size() << endl;
    for (unsigned int i = 0; i < subSampVec->size(); ++i) {
        cout << subSampVec->at(i).histNameSuffix << " is channel " << i << ": " << subSampVec->at(i).descString << endl;
    }
    theApp.Run(retVal);
}
