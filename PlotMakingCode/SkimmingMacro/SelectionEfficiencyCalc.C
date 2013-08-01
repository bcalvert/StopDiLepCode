#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
//#include "TNTuple.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
//#include "TRint.h"
#include <iostream>
#include "TTree.h"
#include "TString.h"
#include <cmath>
#include <sstream>
int main( int argc, const char* argv[] ) {
    //void SelectionEfficiencyCalc(char * inputName, int nBins = 61) {
    TString inputName;
    int nBins = 60;
    bool showPlot = 0;
    bool saveRoot = 0;
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"-i", 2) == 0) {
            inputName = TString(argv[k+1]);
        }
        else if (strncmp (argv[k],"NBins", 5) == 0) {
            nBins = strtol(argv[k+1], NULL, 10 );
        }
        else if (strncmp (argv[k],"showPlot", 8) == 0) {
            showPlot = 1;
            cout << "want to show plot! " << endl;
        }
        else if (strncmp (argv[k],"saveRoot", 8) == 0) {
            saveRoot = 1;
        }
    }
    /*
    if (showPlot) {
        TRint theApp("App", &argc, argv);
        Bool_t retVal = kTRUE;
    }
     */
    //  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
    //  dir.ReplaceAll("basic.C","");
    //  dir.ReplaceAll("/./","/");
    ifstream in;
    TString openName = "./";
    TString outNameBase; outNameBase += inputName;
    TString outName;
    outNameBase.Replace(0, 23, "");
    cout << inputName << endl;
    openName += inputName;
    openName += ".txt";
    cout << "openName " << openName << endl;
    //   in.open((char *) openName);
    in.open(openName);
    outName + "SelEffFiles/";
    outName += outNameBase;
    TString outCanvSaveName = outName;
    outName += ".root";
    cout << "saving to " << outName << endl;
    Float_t eventNum;
    Float_t eventNumNMinOne;
    Int_t nlines = 0;
    TFile * outputRoot;
    if (saveRoot) {
        outputRoot = new TFile(outName,"RECREATE");
    }
    //   int nBins = 61;
    //   float xMin = -0.5;
    float xMin = 0.5;
    float xMax = nBins + 0.5;
    TH1F  * h_NNPlusOneNumDiff = new TH1F("h_NNPlusOneNumDiff","diff in eventNum", nBins,xMin,xMax);
    //   TNtuple * ntuple = new TNtuple("ntuple","data from ascii file","eventNum");
    
    while (1) {
        in >> eventNum;
        //     cout << "eventNum " << eventNum << endl;
        if (!in.good()) break;
        //     if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
        if (nlines > 0) {
            h_NNPlusOneNumDiff->Fill(eventNum - eventNumNMinOne);
        }
        
        //     ntuple->Fill(eventNum);
        nlines++;
        eventNumNMinOne = eventNum;
    }
    printf(" found %d points\n",nlines);
    
    in.close();
    if (saveRoot) {
        outputRoot->Write();
    }
    TFile * inFile = new TFile(outName);
//    TH1F * histRunDiff = (TH1F *) inFile->Get("h_NNPlusOneNumDiff");
    TH1F * histRunDiff = (TH1F *) h_NNPlusOneNumDiff->Clone();
    TF1 * expFit = new TF1("expFit", "expo", 0, nBins - 1);
//    Option_t * fitOption = (showPlot) ? "I" : "IN";
    TCanvas * c1;
    if (showPlot) {
        c1 = new TCanvas("c1", "c1", 600, 600);
        histRunDiff->SetMarkerColor(kBlue);
        histRunDiff->Draw("e1");
        c1->SetLogy();
        histRunDiff->Fit(expFit, "I");
        c1->SaveAs(outCanvSaveName + TString(".pdf"));
    }
    else {
        histRunDiff->Fit(expFit, "IN");   
    }
    float probEst1 = 1/histRunDiff->GetMean();
    float probErrEst1 = histRunDiff->GetMeanError() * (probEst1 * probEst1);
    float probEst2 = 1.0 - TMath::Exp(expFit->GetParameter("Slope"));
    float probErrEst2 = TMath::Exp(expFit->GetParameter("Slope")) * expFit->GetParError(1);
    cout << "First way to estimate Efficiency: 1/Mean -> " << probEst1 << " pm " << probErrEst1 << endl;
    cout << "Second way to estimate Efficiency: 1 - exp(lambda) where lambda is slope -> " << probEst2 << " pm " << probErrEst2 << endl;
  /*
    if (showPlot) {
        theApp.Run(retVal);
        //    theApp.Terminate(0);
    }
    */
}
