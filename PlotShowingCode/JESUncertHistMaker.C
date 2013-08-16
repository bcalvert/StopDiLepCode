#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "TH2F.h"
#include "TFile.h"

//using namespace std;
int main( int argc, char* argv[]) {
    gROOT->ProcessLine("#include <vector>");
    using namespace std;
    TString inputFileName;
    TString outputFileName;
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"-i", 5) == 0) {
            inputFileName = TString(argv[k+1]);
        }
        else if (strncmp (argv[k],"outName", 7) == 0) {
            outputFileName = TString(argv[k+1]);
        }
    }
    inputFileName += ".txt";
    outputFileName += ".root";
    vector<vector<float> > numbers;
    double binsJetPt[45] = {9.0, 11.0, 13.5, 16.5, 19.5, 22.5, 26.0, 30.0, 34.5, 40.0, 46.0, 52.5, 60.0, 69.0, 79.0, 90.5, 105.5, 123.5, 143.0, 163.5, 185.0, 208.0, 232.5, 258.5, 286.0, 331.0, 396.0, 468.5, 549.5, 639.0, 738.0, 847.5, 968.5, 1102.0, 1249.5, 1412.0, 1590.5, 1787.0, 1945.0, 2119.0, 2369.0, 2643.5, 2945.0, 3276.5, 4000};
    double binsJetEta[41] = {-5.4, -5.0, -4.4, -4.0, -3.5, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.4};
    float etaBins = 40;
    TFile * outputFile = new TFile(outputFileName, "RECREATE");
    TH2F * JESErrUp = new TH2F("h_JESSystUncertUp", "+1 #sigma Uncertainty on Data/MC (%); Jet p_{T} [GeV]; Jet #eta", 44, binsJetPt, etaBins, binsJetEta); JESErrUp->Sumw2();
    TH2F * JESErrDown = new TH2F("h_JESSystUncertDown", "-1 #sigma Uncertainty on Data/MC (%); Jet p_{T} [GeV]; Jet #eta", 44, binsJetPt, etaBins, binsJetEta); JESErrDown->Sumw2();

    const int numChars = 2000;;
    char Buffer[numChars];
//    ifstream * inputFileStream;
    //    inputFileStream = new ifstream(inputFileName);
    ifstream inputFileStream(inputFileName);
    string lineData;

    getline(inputFileStream, lineData);
    while (getline(inputFileStream, lineData)) {
//        inputFileStream->getline(Buffer, numChars);
        double d;
        std::vector<float> row;
        stringstream lineStream(lineData);
        while (lineStream >> d) {
            cout << "d " << d << endl;
            row.push_back(d);
        }
        numbers.push_back(row);
        //        vector<float> line(istream_iterator<float>(Buffer), 
//        istream_iterator<float>());
        //        numbers.push_back(line);
    }   
    cout << "numbers size " << numbers.size() << endl;
    cout << "numbers first entry " << numbers.at(0).size() << endl;
    float pT, upErr, downErr;
    float etaLowVal, etaHighVal;
    for (unsigned int IEta = 0; IEta < numbers.size(); ++IEta) {
//    for (unsigned int IEta = 0; IEta < 1; ++IEta) {
        vector<float> EtaNumbers = numbers.at(IEta);
        etaLowVal = EtaNumbers.at(0);
        etaHighVal = EtaNumbers.at(1);
        cout << "etaLowVal " << etaLowVal << endl;
        for (unsigned iEtaNumbers = 0; iEtaNumbers < EtaNumbers.size(); iEtaNumbers += 3) {

            if (iEtaNumbers < 3) continue;
            pT = EtaNumbers.at(iEtaNumbers);
            upErr = EtaNumbers.at(iEtaNumbers + 1);
            downErr = EtaNumbers.at(iEtaNumbers + 2);
            cout << "pT " << pT << endl;
            cout << "upErr " << upErr << endl;
            cout << "downErr " << downErr << endl;
            JESErrUp->Fill(pT + .01, etaLowVal + .01, upErr);
            JESErrDown->Fill(pT + .01, etaLowVal + .01, downErr);
        }
    }
    TRint theApp("App", &argc, argv);
    Bool_t retVal = kTRUE;
    outputFile->Write();
    outputFile->Close();
    
    theApp.Run(retVal);
    //    theApp.Terminate(0);
                              // [1-9][0-9]*\.[0-9]*
                            
                               
                               
}
