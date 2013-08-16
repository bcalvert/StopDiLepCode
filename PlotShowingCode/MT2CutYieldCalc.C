#include <iostream>
#include "../HeaderFiles/BasicFunctions.h"
//#include "../HeaderFiles/HistogramStyleFunctions.h"
//#include "../HeaderFiles/HistogramDataMCPlottingVer2.h"
#include "../HeaderFiles/HistogramSystematics2.h"
#include "../HeaderFiles/StopPlotSetup.h"
#include "../HeaderFiles/StopFunctionDefinitions_v2.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TLegend.h"
//#include "TH3F.h"
//#include "TH3D.h"
//#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRint.h"
//#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <cmath>

//using namespace std;
int main( int argc, char* argv[]) {
    gROOT->ProcessLine("#include <vector>");
    using namespace std;
    int whichChan     = 3;          // which "FullCut" channel to look at -- 0 is MuMu, 1 is EE, 2 is EMu, 3 is all together
    int whichNTuple   = 1;          //as with the plot making code, leave as 1 for now -- 0 is Oviedo, 1 is DESY    
    int whichTTbarGen = 0;          // 0 is Madgraph, 1 is MC@NLO, 2 is Powheg
    bool doExcSamps   = 0;          // For grabbing exclusive (DY + N Jets, TTBar Decay modes) or inclusive samples (As of 8/5/13, only applies to Oviedo)
    bool doNonSig     = 0;          // For whether or not to grab SM background and data
    bool doSignal     = 0;          // For whether or not to grab a signal point
    TString typeSMS   = "";         // Which type of SMS to grab -- either "T2tt" or "T2bw" (as of 8/5/13) only has T2tt FineBin
    TString prefixT2tt   = "";      // prefix for which kind of T2tt to grab
    vector<int> * vecStopMassGrab = new vector<int>;       // vector to hold the list of Stop masses to brab
    vector<int> * vecChi0MassGrab = new vector<int>;       // vector to hold the list of Chi0 masses to brab
    vector<int> * vecCharginoMassGrab = new vector<int>;   // vector to hold the list of Chargino masses to brab
    bool doPURW       = 0;          // grab the nVtx reweighted MC files
    bool doSyst       = 1;          // look at systematics plots
    bool addThings    = 1;          // Add together similar kinds of events (for aesthetic reasons) like VV backgrounds -- (6/25/13) don't turn off for now haven't validated code fully when not adding
    bool useDDEstimate = 0;         // whether or not to use data-driven estimates for appropriate backgrounds -- as of right now it is just the TTBar norm to MT2ll < 80 GeV (7/22/13)
    bool allMT2llSystematic = 0;    // Whether or not to use the MT2ll systematic smearing for all MC or just TTBar
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"wChan", 5) == 0) {
            whichChan = strtol(argv[k+1], NULL, 10 );
        }
        else if (strncmp (argv[k],"wNTuple", 7) == 0) {
            whichNTuple = strtol(argv[k+1], NULL, 10 );
        }
        else if (strncmp (argv[k],"wTTbarGen", 9) == 0) {
            whichTTbarGen = strtol(argv[k+1], NULL, 10 );
        }
        else if (strncmp (argv[k],"doExcSamps", 10) == 0) {
            doExcSamps = 1;
        }
        else if (strncmp (argv[k],"doNonSig", 8) == 0) {
            doNonSig = 1;
        }    
        else if (strncmp (argv[k],"doSignal", 8) == 0) {
            doSignal = 1;
            cout << "NB: Format for post command line arguments is TypeSMS, prefix for the file to grab, StopMassToGrab, Chi0MassToGrab, CharginoMassToGrab " << endl;
            typeSMS = TString(argv[k+1]);
            prefixT2tt = TString(argv[k+2]);
            vecStopMassGrab->push_back(strtol(argv[k+3], NULL, 10 ));
            vecChi0MassGrab->push_back(strtol(argv[k+4], NULL, 10 ));
            vecCharginoMassGrab->push_back(strtol(argv[k+5], NULL, 10 ));
        }      
        else if (strncmp (argv[k],"doPURW", 6) == 0) {
            doPURW = 1;
        }
        else if (strncmp (argv[k],"noSyst", 6) == 0) {
            doSyst = 0;
        }
        else if (strncmp (argv[k],"noAdd", 5) == 0) {
            addThings = 0;
        }
        else if (strncmp (argv[k],"useDDEst", 8) == 0) {
            useDDEstimate = 1;
        }
        else if (strncmp (argv[k],"allMT2ll", 8) == 0) {
            allMT2llSystematic = 1;
        }
    }
    TRint theApp("App", &argc, argv);
    Bool_t retVal = kTRUE;
    
    if (!doExcSamps && whichNTuple == 1 ) {
        doExcSamps = true;
        cout << "setting do Exclusive Samples to true because running on DESY" << endl;
    }
    //Set up the file input
    //    vector<TFile *> * inFiles = new vector<TFile*>;
    vector<TString> * fileInNames = StopFileNames(whichNTuple, whichTTbarGen, doExcSamps);
    vector<TFile *> * inputFiles  = StopFiles(whichNTuple, fileInNames, doExcSamps, whichTTbarGen, doPURW, doSyst);
    
    vector<TString> * sampleAddNames = new vector<TString>;
    vector<int> * sampleStartPositions = new vector<int>;
    sampleStartPositionsNames(whichNTuple, sampleAddNames, sampleStartPositions, doExcSamps);
    
    vector<TFile *> * inputFilesSignal;
    vector<TString> * mcLegendsSignal;
    vector<Color_t> * mcColorsSignal;
    vector<Style_t> * mcStylesSignal;
    if (doSignal) {
        inputFilesSignal              = StopSignalFiles(whichNTuple, typeSMS, prefixT2tt, vecStopMassGrab, vecChi0MassGrab, vecCharginoMassGrab, doPURW, doSyst);
        mcLegendsSignal               = MCSignalLegends(typeSMS, vecStopMassGrab, vecChi0MassGrab, vecCharginoMassGrab);
        mcColorsSignal                = MCSignalColors(vecStopMassGrab->size());
        mcStylesSignal                = MCSignalStyles(vecStopMassGrab->size());
    }
    vector<SampleT> * subSampVec    = SubSampVec();
    cout << "subsamp size " << subSampVec->size() << endl;
    vector<SystT> * systVec         = SystVec();
    
    //some relevant things for saving names        
    TString TTBarGenName[3] = {"_madgraph", "_mcatnlo", "_powheg"};
    TString nameNTuple = (whichNTuple == 0) ? "_Ovi" : "_DESY";
    bool doSymErr = 0;
    TString stringSymmErr = (doSymErr) ? "_wSymmErr" : "";
    TString stringDDEstimate = (useDDEstimate) ? "_wDDEst" : "";
    TString stringExcSamp = (doExcSamps && whichNTuple == 0) ? "_ExcSamps" : "";

    TString canvSuffixSaveName = TTBarGenName[whichTTbarGen];
    canvSuffixSaveName += nameNTuple;
    canvSuffixSaveName += stringSymmErr;
    canvSuffixSaveName += stringDDEstimate;
    canvSuffixSaveName += stringExcSamp;
    
    ///Systematics stuff////
    vector<TH1F *> * dataHist1DVec;
    TH1F * h_DataComp, * h_MCComp;
    
    vector<TH1 *> * dataHistTH1Vec;
    vector<TH1 *> * mcHistTH1Vec;

    vector<TH1 *> * mcHistSystTH1Vec;
    
    vector<TH1F *> * mcIndHist1DCentValVec; //will contain central value histos for individual MC samples
    vector<TH1F *> * mcCompHist1DCentValVec; //will contain the added histos for general categories
    vector<TH1F *> * mcCompHist1DSystVec; //will contain the added histos across all MC for given systematic
    

    vector<float> * nVtxBackScaleVec = ScaleBackVecCalc(inputFiles);
    
    vector<TGraphAsymmErrors *> * errCompSpecSource; // will contain just the systematic error for each specific source -- note, will contain errSystQuadSum as a final dude
    vector<TGraphAsymmErrors *> * errCompSpecSource_pStat; // for each specific source, will contain Stat + respective systematic error -- note, will contain errSystQuadSum_pStat as a final dude
    vector<TString> * systCanvNameVec;
    vector<TGraphAsymmErrors *> * fracRatioSystVec;
    
    
    TH1 * currSignalCentValTH1Hist;
    vector<TH1 *> * vecCurrSignalSystTH1Hists;
    
    TH1F * currSignal1DCentValHist;
    vector<TH1F *> * vecCurrStop1DSystHists;
    vector<TGraphAsymmErrors *> * vecCurrSigCompSpecSource, * vecCurrSigCompSpecSource_pStat;
    
    ///Data-Driven systematics stuff
    float TTBarFullCutSFOvi[3] = {0.738946, 0.776429, 0.832783};
    float TTBarFullCutSFDESY[3] = {0.84936, 0.843207, 0.949907};
    float TTBarSF = (whichNTuple == 0) ? TTBarFullCutSFOvi[whichTTbarGen] : TTBarFullCutSFDESY[whichTTbarGen];
    
    
    //plotting stuff
    float indLumiDESY[4] = {892, 4404, 7032, 7274};
    float indLumiOvi[4] = {892, 4404, 7032, 7274};
    float intLumi = 19602.901;       
    if (whichNTuple == 0) {
        intLumi = indLumiOvi[0] + indLumiOvi[1] + indLumiOvi[2] + indLumiOvi[3];
    }
    else {
        intLumi = indLumiDESY[0] + indLumiDESY[1] + indLumiDESY[2] + indLumiDESY[3];
    }    
    vector<float> * signalSkimScaleVec;
    if (doSignal) signalSkimScaleVec = SignalSkimEfficiencyCalc(typeSMS, prefixT2tt, vecStopMassGrab, vecChi0MassGrab, vecCharginoMassGrab, intLumi);
    
    TString plotGrabBaseName = "h_PassMT2llCut";
    TString plotGrabName;
    TString plotSystGrabName;
    TString plotVarName, subSampName;;
    TString dataLegendComp = "Data";
    float   fracRatioYAxisRange = 0.21;
    
    /*
    const int numMT2llCuts = 5;
    int MT2llCuts[numMT2llCuts] = {80, 90, 100, 110, 120};     
    */
    const int numMT2llCuts = 1;
    int MT2llCuts[numMT2llCuts] = {80};

    int grabChan[4] = {7, 24, 39, 54};

    for (int iMT2Cut = 0; iMT2Cut < numMT2llCuts; ++iMT2Cut) {
        dataHist1DVec = new vector<TH1F *>;
        mcIndHist1DCentValVec = new vector<TH1F *>;
        mcCompHist1DCentValVec = new vector<TH1F *>;
        mcCompHist1DSystVec = new vector<TH1F *>;
        fracRatioSystVec = new vector<TGraphAsymmErrors *>;
        errCompSpecSource = new vector<TGraphAsymmErrors *>;
        errCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
        systCanvNameVec = new vector<TString>;
        dataHistTH1Vec = new vector<TH1 *>;
        mcHistTH1Vec = new vector<TH1 *>;
        mcHistSystTH1Vec = new vector<TH1 *>;
        
        subSampName = "";
        subSampName += subSampVec->at(grabChan[whichChan]).histNameSuffix;
        plotGrabName = plotGrabBaseName;
        plotGrabName += MT2llCuts[iMT2Cut];
        plotSystGrabName = plotGrabName;                
        plotGrabName += subSampVec->at(grabChan[whichChan]).histNameSuffix;
        
        if (doNonSig) {
            HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, plotGrabName, subSampName, useDDEstimate, TTBarSF);
            HistogramAdderData(dataHistTH1Vec, h_DataComp, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");     
            
            cout << "h_DataComp Name " << h_DataComp->GetName() << endl;
            cout << "h_DataComp Bin Content 1 " << h_DataComp->GetBinContent(1) << endl;
            cout << "h_DataComp Bin Error 1 " << h_DataComp->GetBinError(1) << endl;
            cout << "h_DataComp Bin Content 2 " << h_DataComp->GetBinContent(2) << endl;
            cout << "h_DataComp Bin Error 2 " << h_DataComp->GetBinError(2) << endl;
            
            HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, plotGrabName, subSampName, useDDEstimate, TTBarSF);
            HistogramAdderMC(mcHistTH1Vec, mcCompHist1DCentValVec, sampleStartPositions, sampleAddNames, h_MCComp, whichNTuple, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
            
            cout << "h_MCComp Name " << h_MCComp->GetName() << endl;
            cout << "h_MCComp Bin Content 1 " << h_MCComp->GetBinContent(1) << endl;
            cout << "h_MCComp Bin Error 1 " << h_MCComp->GetBinError(1) << endl;
            cout << "h_MCComp Bin Content 2 " << h_MCComp->GetBinContent(2) << endl;
            cout << "h_MCComp Bin Error 2 " << h_MCComp->GetBinError(2) << endl;
            
            if (doSyst) {
                HistogramVecGrabberSystGrab(inputFiles, mcHistTH1Vec, mcHistSystTH1Vec, nVtxBackScaleVec, plotSystGrabName, subSampName, systVec, useDDEstimate, TTBarSF, allMT2llSystematic, whichNTuple);            
                HistogramProjectorSyst(mcHistSystTH1Vec, mcCompHist1DSystVec, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                SystGraphMakers(h_MCComp, mcCompHist1DSystVec, errCompSpecSource, errCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, kGray + 1, plotGrabName, true, fracRatioYAxisRange, doSymErr, false);                
                for (unsigned int iSyst = 0; iSyst < systCanvNameVec->size(); ++iSyst) {
                    cout << "For syst " << systCanvNameVec->at(iSyst) << endl;
                    cout << "ErrGraph Up Err at point 1 " << errCompSpecSource->at(iSyst)->GetErrorYhigh(1) << endl;
                    cout << "ErrGraph Down Err at point 1 " << errCompSpecSource->at(iSyst)->GetErrorYlow(1) << endl;
                    cout << "ErrGraph Up Err at point 2 " << errCompSpecSource->at(iSyst)->GetErrorYhigh(2) << endl;
                    cout << "ErrGraph Down Err at point 2 " << errCompSpecSource->at(iSyst)->GetErrorYlow(2) << endl;
                }
            }
        }
        if (doSignal) {
            vecCurrSignalSystTH1Hists = new vector<TH1 *>;
            vecCurrStop1DSystHists = new vector<TH1F *>;
            vecCurrSigCompSpecSource = new vector<TGraphAsymmErrors *>;
            vecCurrSigCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
            
            
            HistogramVecGrabber_Signal(inputFilesSignal, signalSkimScaleVec, 0, currSignalCentValTH1Hist, plotSystGrabName, vecCurrSignalSystTH1Hists, systVec, subSampName, false, false, doSyst, allMT2llSystematic);
            HistogramAdderSignal(currSignalCentValTH1Hist, currSignal1DCentValHist, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
            cout << "currSignal1DCentValHist Name " << currSignal1DCentValHist->GetName() << endl;
            cout << "currSignal1DCentValHist Bin Content 1 " << currSignal1DCentValHist->GetBinContent(1) << endl;
            cout << "currSignal1DCentValHist Bin Error 1 " << currSignal1DCentValHist->GetBinError(1) << endl;
            cout << "currSignal1DCentValHist Bin Content 2 " << currSignal1DCentValHist->GetBinContent(2) << endl;
            cout << "currSignal1DCentValHist Bin Error 2 " << currSignal1DCentValHist->GetBinError(2) << endl;
            if (doSyst) {
                
                HistogramProjectorSyst(vecCurrSignalSystTH1Hists, vecCurrStop1DSystHists, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                SystGraphMakers(currSignal1DCentValHist, vecCurrStop1DSystHists, vecCurrSigCompSpecSource, vecCurrSigCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, mcColorsSignal->at(0), plotGrabName, true, fracRatioYAxisRange, doSymErr, true);
                for (unsigned int iSyst = 0; iSyst < systCanvNameVec->size(); ++iSyst) {
                    cout << "For syst " << systCanvNameVec->at(iSyst) << endl;
                    cout << "Signal ErrGraph Up Err at point 1 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYhigh(1) << endl;
                    cout << "Signal ErrGraph Down Err at point 1 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYlow(1) << endl;
                    cout << "Signal ErrGraph Up Err at point 2 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYhigh(2) << endl;
                    cout << "Signal ErrGraph Down Err at point 2 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYlow(2) << endl;
                }
            }            
        }        
    }
    theApp.Run(retVal);
    //    theApp.Terminate(0);
}
