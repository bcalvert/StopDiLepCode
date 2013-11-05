#include <iostream>
#include "../HeaderFiles/BasicFunctions.h"
//#include "../HeaderFiles/HistogramStyleFunctions.h"
//#include "../HeaderFiles/HistogramDataMCPlottingVer2.h"
#include "../HeaderFiles/HistogramSystematics2.h"
#include "../HeaderFiles/StopStructDefinitions.h"
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
    int whichChan     = 0;          //which "channel" to run on, this doesn't just mean ee, emu, or mumu, but can mean events in the full cut sequence, events in the ZMass window but otherwise full cut sequence, etc.
    int whichNTuple   = 1;          //as with the plot making code, leave as 1 for now -- 0 is Oviedo, 1 is DESY    
    int whichTTbarGen = 0;          // 0 is Madgraph, 1 is MC@NLO, 2 is Powheg
    bool doExcSamps   = 0;          // For grabbing exclusive (DY + N Jets, TTBar Decay modes) or inclusive samples (As of 8/5/13, only applies to Oviedo)
    bool doSignal     = 0;          // For whether or not to grab a signal point
    bool multSigPts   = 0;          // For whether or not to show multiple different signal points (As of 8/5/13 this functionality hasn't been added in yet)
    TString typeSMS   = "";         // Which type of SMS to grab -- either "T2tt" or "T2bw" (as of 8/5/13) only has T2tt FineBin
    TString prefixT2tt   = "";      // prefix for which kind of T2tt to grab
    vector<int> * vecStopMassGrab = new vector<int>;       // vector to hold the list of Stop masses to brab
    vector<int> * vecChi0MassGrab = new vector<int>;       // vector to hold the list of Chi0 masses to brab
    vector<int> * vecCharginoMassGrab = new vector<int>;   // vector to hold the list of Chargino masses to brab
    bool doPURW       = 0;          // grab the nVtx reweighted MC files
    bool doSyst       = 1;          // look at systematics plots -- (6/25/13) don't turn off for now haven't validated code fully when not making systematicslumi
    bool addThings    = 1;          // Add together similar kinds of events (for aesthetic reasons) like VV backgrounds -- (6/25/13) don't turn off for now haven't validated code fully when not adding
    bool calcTTBarNorm = 0;         // calculate TTBar normalization by utilizing integral to data - (other backgrounds) for MT2ll < 80 in the "Full Cut region"
    bool useDDEstimate = 0;         // whether or not to use data-driven estimates for appropriate backgrounds -- as of right now it is just the TTBar norm to MT2ll < 80 GeV (7/22/13)
    bool allMT2llSystematic = 0;    // Whether or not to use the MT2ll systematic smearing for all MC or just TTBar
    bool doOneDee     = 1;          // Plots 1D histograms, as in the 1D HistTs defined in StopFunctionDefinitions_v2.h
    bool doTwoDee     = 0;          // Plots 2D histograms, as in the 2D HistTs defined in StopFunctionDefinitions_v2.h -- note, this involves making projections and what-not so there (as of 7/10/13) isn't actual plotting of 2D histograms
    bool doThreeDee   = 0;          // Plots 3D histograms, as in the 3D HistTs defined in StopFunctionDefinitions_v2.h -- note, this involves making projections and what-not so there (as of 7/10/13) isn't actual plotting of 3D histograms
    bool doIsoPlots   = 0;          // Plots the electron isolation variables to check what is the deal with them for Oviedo vs. DESY
    bool reportChi2   = 0;          // nominally, reports data/MC chi2 value in the control region
    bool doCDF        = 0;          // nominally, creates Cumulative Distribution Functions for MC and Data -- (6/25/13) haven't tested yet
    bool doStats      = 0;          // show statistics for data and mc spectra on plots -- e.g. mean and RMS -- (6/25/13) haven't tested yet to see how positioning looks
    bool doAbsRatio   = 0;          // Affects what gets shown in fractional ratio plots -- 0: (MC - Data) / Data, 1: Data/MC
    bool doOviDESYFracRatio = 0;    // If on, shows frac ratio plots for both Oviedo and DESY on the given plot (although the top plot will still be whatever you chose)
    bool doSmartFracRatio   = 0;    // If on, makes it so that the ratio plot centers itself to include the ratio plot (i.e. matters when the ratio plot is going out of the nominal axis
    bool saveDotCFile = 0;          // If on, saves a .C version of every plot that is made
    
    bool singSampCompare  = 0;          // For doing comparisons between some set of two single samples (for example, given SM bkgds); can be just a single sample to see what it looks like
    bool showLegend     = 1;
    bool diffFilesSingSampCompare = 0;  // For grabbing same sample across different files
    TString firstSampSetupFile = "";     // file name for the first file containing the information for what sample/histogram to look at
    TString secondSampSetupFile = "";     // file name for the second file containing the information for what sample/histogram to look at
    bool multHistsSingSampCompare = 1;  // if there will be multiple histograms to look at 
    bool multiChannelAdd = 0;
    TString multiChannelSetupFile = "";
    int  versNumber     = 1;
    bool doReReco       = 0;
    bool SmearedPlots   = 0;
    int  JetsSmeared    = 0;
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
        else if (strncmp (argv[k],"calcTTBarNorm", 13) == 0) {
            calcTTBarNorm = 1;
        }
        else if (strncmp (argv[k],"useDDEst", 8) == 0) {
            useDDEstimate = 1;
        }
        else if (strncmp (argv[k],"allMT2ll", 8) == 0) {
            allMT2llSystematic = 1;
        }
        else if (strncmp (argv[k],"noOneDee", 8) == 0) {
            doOneDee = 0;
        }
        else if (strncmp (argv[k],"doTwoDee", 6) == 0) {
            doTwoDee = 1;
        }
        else if (strncmp (argv[k],"doThreeDee", 8) == 0) {
            doThreeDee = 1;
        }
        else if (strncmp (argv[k],"doIsoPlots", 10) == 0) {
            doIsoPlots = 1;
        }
        else if (strncmp (argv[k],"noLeg", 5) == 0) {
            showLegend = 0;
        }
        else if (strncmp (argv[k],"Chi2", 4) == 0) {
            reportChi2 = 1;
        }
        else if (strncmp (argv[k],"doAbsRatio", 10) == 0) {
            doAbsRatio = 1;
        }
        else if (strncmp (argv[k],"doODFR", 6) == 0) {
            doOviDESYFracRatio = 1;
        }
        else if (strncmp (argv[k],"doSFR", 5) == 0) {
            doSmartFracRatio = 1;
        }
        else if (strncmp (argv[k],"doStats", 7) == 0) {
            doStats = 1;
        }
        else if (strncmp (argv[k],"JsSm", 4) == 0) {
            SmearedPlots = 1;
            JetsSmeared = strtol(argv[k+1], NULL, 10);
        }
        else if (strncmp (argv[k],"sDCF", 4) == 0) {
            saveDotCFile = 1;
        }
        else if (strncmp (argv[k],"SSComp", 6) == 0) {
            singSampCompare = 1;
            firstSampSetupFile = TString(argv[k+1]);
            firstSampSetupFile += ".txt";
        }
        else if (strncmp (argv[k],"diffFileSSComp", 14) == 0) {
            diffFilesSingSampCompare = 1;
            secondSampSetupFile = TString(argv[k+1]);
            secondSampSetupFile += ".txt";
        }
        else if (strncmp (argv[k],"multiChan", 9) == 0) {
            multiChannelAdd = 1;
            multiChannelSetupFile = TString(argv[k+1]);
            multiChannelSetupFile += ".txt";
        }
        else if (strncmp (argv[k],"versNum", 7) == 0) {
            versNumber = strtol(argv[k+1], NULL, 10 );
        }       
        else if (strncmp (argv[k],"doReReco", 8) == 0) {
            doReReco = 1;
        }
    }
    TRint theApp("App", &argc, argv);
    Bool_t retVal = kTRUE;
    //void StopPlotsSpectraShowNewDESY_wSyst(int StopMass = 200, int ChiZeroMass = 50, int whichChan = 0, int whichNTuple = 0, bool doPURW = 0, int whichTTbarGen = 0, bool addThings = 1, bool reportChi2 = 0, bool doCDF = 0) {
    float m = 1.3;
    int n_ = 1;
    Double_t e_ = 15*m; //15
    Double_t k_ = 7*m; //15
    Double_t t_ = 30*m; //30
    Double_t b_ = 50*m; //80
    Double_t g_ = 87*m; //87
    Double_t d_ = 30*m; //30
    Double_t h_ = 350*m; //400
    Double_t w_ = 350*m; //350
    Double_t hl_ = 70*m;
    
    //    Double_t ww_ = g_ + w_ + e_ ;
    Double_t W_ = g_ + n_*w_ + 2*(n_-1)*e_ + d_;
    Double_t H_ = t_ + h_ + b_ + hl_ +2*k_ ;
    Int_t wtopx = 100;
    Int_t wtopy = 100;
    
    // Style things
    gStyle->SetErrorX(0.5); // set X uncertainty to 1/2 of a bin
    gStyle->SetTitleFont(42);
    
    
    // Fixing things based on which nTuplevariables
    
    if (!doExcSamps && whichNTuple == 1 ) {
        doExcSamps = true;
        cout << "setting do Exclusive Samples to true because running on DESY" << endl;
    }
    //Set up the file input
    //    vector<TFile *> * inFiles = new vector<TFile*>;
    vector<TString> * fileInNames = StopFileNames(whichNTuple, whichTTbarGen, doExcSamps, doReReco);
    vector<TFile *> * inputFiles  = StopFiles(whichNTuple, fileInNames, doExcSamps, whichTTbarGen, doPURW, doSyst, versNumber);
    vector<TString> * mcLegends   = MCLegends(whichNTuple, addThings);
    vector<Color_t> * mcColors    = MCColors(whichNTuple, addThings);
    
    vector<TString> * fileInNames_Other;
    vector<TFile *> * inputFiles_Other;
    if (doOviDESYFracRatio) {
        fileInNames_Other = StopFileNames((whichNTuple + 1) % 2, whichTTbarGen, doExcSamps, doReReco);
        inputFiles_Other  = StopFiles((whichNTuple + 1) % 2, fileInNames_Other, doExcSamps, whichTTbarGen, doPURW, doSyst, versNumber);
    }
    
    vector<TString> * sampleAddNames = new vector<TString>;
    vector<int> * sampleStartPositions = new vector<int>;
    sampleStartPositionsNames(whichNTuple, whichTTbarGen, sampleAddNames, sampleStartPositions, doExcSamps);    
    vector<int> * sortIndexVec = sortOrder(whichNTuple);
    
    vector<TString> * sampleAddNames_Other;
    vector<int> * sampleStartPositions_Other;
    if (doOviDESYFracRatio) {
        sampleAddNames_Other = new vector<TString>;
        sampleStartPositions_Other = new vector<int>;
        sampleStartPositionsNames((whichNTuple + 1) % 2, whichTTbarGen, sampleAddNames_Other, sampleStartPositions_Other, doExcSamps);    
    }
    
    vector<TFile *> * inputFilesSignal;
    vector<TString> * mcLegendsSignal;
    vector<Color_t> * mcColorsSignal;
    vector<Style_t> * mcStylesSignal;
    if (doSignal) {
        inputFilesSignal              = StopSignalFiles(whichNTuple, typeSMS, prefixT2tt, vecStopMassGrab, vecChi0MassGrab, vecCharginoMassGrab, doPURW, doSyst, doReReco);
        mcLegendsSignal               = MCSignalLegends(typeSMS, vecStopMassGrab, vecChi0MassGrab, vecCharginoMassGrab);
        mcColorsSignal                = MCSignalColors(vecStopMassGrab->size());
        mcStylesSignal                = MCSignalStyles(vecStopMassGrab->size());
    }
    //    unsigned int numFiles = fileInNames->size();
    
    TString fracRatioNumerName = (doAbsRatio) ? "Data" : "(MC-Data)";
    TString fracRatioDenomName = (doAbsRatio) ? "MC" : "Data"; 
    float   fracRatioYAxisRange = (doAbsRatio) ? 0.21 : 0.26;
    
//    vector<HistogramT> * histVec_1D = OneDeeHistTVec();
    vector<HistogramT> * histVec_1D = OneDeeHistTVecStatusReport();
    vector<HistogramT> * histVec_2D = TwoDeeHistTVec();
    vector<HistogramT> * histVec_3D = ThreeDeeHistTVec(); 
    vector<HistogramT> * histVec_1DSmear = OneDeeSmearHistTVec(JetsSmeared);
    vector<HistogramT> * histVec1DtoUse = SmearedPlots ? histVec_1DSmear : histVec_1D;
    vector<SampleT> * subSampVec    = SubSampVec();
    cout << "subsamp size " << subSampVec->size() << endl;
    vector<SystT> * systVec         = SystVec(SmearedPlots);
    vector<TString> * vecIsoPlotNames = IsoPlotNames(whichNTuple);
    
    //some relevant things for saving names        
    TString TTBarGenName[3] = {"_madgraph", "_mcatnlo", "_powheg"};
    TString nameNTuple = (whichNTuple == 0) ? "_Ovi" : "_DESY";
    bool doSymErr = 0;
    TString stringSymmErr = (doSymErr) ? "_wSymmErr" : "";
    TString stringDDEstimate = (useDDEstimate) ? "_wDDEst" : "";
    TString stringExcSamp = (doExcSamps && whichNTuple == 0) ? "_ExcSamps" : "";
    TString stringSignal = "";
    if (doSignal) {
        stringSignal = "_wSignal";
        if (multSigPts) {
            stringSignal += "MultSigPts";
        }
        else {
            stringSignal += "_Stop_";
            stringSignal += vecStopMassGrab->at(0);
            stringSignal += "_Chi0_";
            stringSignal += vecChi0MassGrab->at(0);
            if (typeSMS.Contains("T2bw")) {
                stringSignal += "_Chargino_";
                stringSignal += vecCharginoMassGrab->at(0);
            }
        }
    }   
    TString canvSuffixSaveName = TTBarGenName[whichTTbarGen];
    canvSuffixSaveName += nameNTuple;
    canvSuffixSaveName += stringSymmErr;
    canvSuffixSaveName += stringDDEstimate;
    canvSuffixSaveName += stringExcSamp;
    canvSuffixSaveName += stringSignal;
    if (doReReco) canvSuffixSaveName += "_ReReco";
    if (doAbsRatio) canvSuffixSaveName += "_AbsRatio";
    if (!showLegend) canvSuffixSaveName += "_noLegend";
    if (SmearedPlots) canvSuffixSaveName += "_SmearedMET";
    if (versNumber == 2) canvSuffixSaveName += "_vers2";
    
    ///Systematics stuff////
    vector<TH1F *> * dataHist1DVec;
    vector<TH2F *> * dataHist2DVec;
    vector<TH3F *> * dataHist3DVec;
    TH1 * h_MCCompTH1;
    TH1F * h_DataComp, * h_CDF_Data;
    TH1F * h_MCComp, * h_ErrComp;
    TH1F * h_DataComp_Other, * h_MCComp_Other;
    TH1F * h_FracratioComp, * h_FracratioComp_Other;
    
    TH1F * h_DataComp2D;
    TH1F * h_MCComp2D;
    
    TH1F * h_DataComp3D;
    TH1F * h_MCComp3D;
    
    
    vector<TH1 *> * dataHistTH1Vec, * dataHistTH1Vec_Other;
    vector<TH1 *> * mcHistTH1Vec, * mcHistTH1Vec_Other;
    
    vector<TH1 *> * mcHistSystTH1Vec;
    
    vector<TH1F *> * mcIndHist1DCentValVec; //will contain central value histos for individual MC samples
    vector<TH1F *> * mcCompHist1DCentValVec; //will contain the added histos for general categories
    vector<TH1F *> * mcCompHist1DCentValVec_Other; //will contain the added histos for general categories
    vector<TH1F *> * mcCompHist1DCentValVecSorted; //will contain the added histos for general categories
    vector<TH1F *> * mcCompHist1DSystVec; //will contain the added histos across all MC for given systematic
    
    vector<TH2F *> * mcIndHist2DCentValVec;
    vector<TH2F *> * mcCompHist2DCentValVec;
    vector<TH2F *> * mcCompHist2DSystVec;
    
    vector<TH3F *> * mcIndHist3DCentValVec;
    vector<TH3F *> * mcCompHist3DCentValVec;
    vector<TH3F *> * mcCompHist3DSystVec;    
    
    vector<TGraphAsymmErrors *> * fracRatioSystVec;
    vector<TH1F *> * mcCDFSystVec;
    vector<float> * nVtxBackScaleVec = ScaleBackVecCalc(inputFiles);
    vector<float> * nVtxBackScaleVec_Other;
    if (doOviDESYFracRatio) {
        nVtxBackScaleVec_Other = ScaleBackVecCalc(inputFiles_Other);
    }
    
    TGraphAsymmErrors * errCompStatCentVal, * errSystQuadSum, * errSystQuadSum_pStat;
    TGraphAsymmErrors * errLepEnSc, * errLepEffSF, * errMT2ll, * errGenTopRW;
    TGraphAsymmErrors * errLepEnSc_pStat, * errLepEffSF_pStat, * errMT2ll_pStat, * errGenTopRW_pStat;
    TString stringLepEffSF = "LepEffSF";
    TString stringLepEnSc = "LepES";
    TString stringMT2ll = "MT2ll";  
    TString stringGenTopRW = "genTopRW";
    TString stringStopXSecUncert = "genStopXSec";
    TGraphAsymmErrors * currFracRatioGraph;
    vector<TGraphAsymmErrors *> * errCompSpecSource; // will contain just the systematic error for each specific source -- note, will contain errSystQuadSum as a final dude
    vector<TGraphAsymmErrors *> * errCompSpecSource_pStat; // for each specific source, will contain Stat + respective systematic error -- note, will contain errSystQuadSum_pStat as a final dude
    vector<TString> * systCanvNameVec;
    
    
    TH1 * currSignalCentValTH1Hist;
    vector<TH1 *> * vecCurrSignalSystTH1Hists;
    
    TH1F * currSignal1DCentValHist;
    TH2F * currSignal2DCentValHist;
    TH3F * currSignal3DCentValHist;
    vector<TH1F *> * vecCurrStop1DSystHists;
    vector<TH2F *> * vecCurrStop2DSystHists;
    vector<TH3F *> * vecCurrStop3DSystHists;
    TGraphAsymmErrors * currSignalErrStatCentVal, * currSignalErrSystQuadSum, * currSignalErrSystQuadSum_pStat;
    TGraphAsymmErrors * currSignalErrSystLepEnSc, * currSignalErrSystLepEffSF, * currSignalErrSystMT2ll, * currSignalErrSystStopXSecUncert;
    TGraphAsymmErrors * currSignalErrSystLepEnSc_pStat, * currSignalErrSystLepEffSF_pStat, * currSignalErrSystMT2ll_pStat, * currSignalErrSystStopXSecUncert_pStat;
    vector<TGraphAsymmErrors *> * vecCurrSigCompSpecSource, * vecCurrSigCompSpecSource_pStat;
    
    vector<TH1F *> * vecStop1DCentValHists;
    vector<TH2F *> * vecStop2DCentValHists;
    vector<TH3F *> * vecStop3DCentValHists;
    vector<vector<TH1F *> *> * vecStop1DSystHists;
    vector<vector<TH2F *> *> * vecStop2DSystHists;
    vector<vector<TH3F *> *> * vecStop3DSystHists;
    vector<TGraphAsymmErrors *> * vecErrSigStatCentVal, * vecErrSigSystQuadSum, * vecErrSigSystQuadSum_pStat;
    vector<TGraphAsymmErrors *> * vecErrSigLepEnSc, * vecErrSigLepEffSF, * errSigMT2ll, * vecErrSigStopXSecUncert;
    vector<TGraphAsymmErrors *> * vecErrSigLepEnSc_pStat, * vecErrSigLepEffSF_pStat, * vecErrSigMT2ll_pStat, * vecErrSigStopXSecUncert_pStat;
    vector<vector<TGraphAsymmErrors *> *> * vecErrSigCompVecSpecSource, * vecErrSigCompVecSpecSource_pStat;
    
    
    ////Single sample plot looks
    char Buffer[500];
    bool doAllPlots = false;
    bool scaleTosame = true;
    TString s2sName = (scaleTosame) ? "ON!" : "OFF!";
    cout << "Scaling compared things to same integral is turned " << s2sName << endl;
    float histOneIntegral, histTwoIntegral;
    TString addPath = "../RootFiles/";
    TString singPlot2MakeFirstFile = "";
    TString singPlot2MakeSecondFile = "";
    TString firstSampInFileName = "";
    TString firstSampInFileBaseName = "";
    TFile * firstSampInFile;
    ifstream * firstSampFileStream;
    TString secondSampInFileName = "";
    TString secondSampInFileBaseName = "";
    TFile * secondSampInFile;
    ifstream * secondSampFileStream;
    
    ifstream * multiChannelStream;
    TString multiChannelCompOutputName;
    vector<int> * multiChannelIDs = new vector<int>;
    
    TH1F * histOne1D, * histTwo1D;
    TH2F * histOne2D, * histTwo2D;
    TH3F * histOne3D, * histTwo3D;
    vector<TH1F *> * histOne1DSystVec, * histTwo1DSystVec;
    vector<TH2F *> * histOne2DSystVec, * histTwo2DSystVec;
    vector<TH3F *> * histOne3DSystVec, * histTwo3DSystVec;
    
    vector<TH1F *> * histOneCDF1DSystVec, * histTwoCDF1DSystVec;
    vector<TH2F *> * histOneCDF2DSystVec, * histTwoCDF2DSystVec;
    vector<TH3F *> * histOneCDF3DSystVec, * histTwoCDF3DSystVec;
    
    Color_t histColors[2] = {kBlue, kRed};
    /*
     //Note these will get used in the comparison but don't mean anything otherwise
     fracRatioSystVec = new vector<TGraphAsymmErrors *>;
     errCompSpecSource = new vector<TGraphAsymmErrors *>;
     errCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
     systCanvNameVec = new vector<TString>;
     */
    TString histOneplot, histOneplot_preRW, histTwoplot, histTwoplot_preRW;
    TString addPlotVarName;
    
    if (singSampCompare) {
        firstSampFileStream = new ifstream(firstSampSetupFile);
        if (!(firstSampFileStream->eof())) {
            firstSampFileStream->getline(Buffer, 500);
            firstSampInFileBaseName += TString(string(Buffer));
            firstSampInFileName += TString(string(Buffer));
            firstSampInFileName += (doPURW) ? "_PURW": "";
            firstSampInFileName += (doSyst) ? "_wSyst": "";
            firstSampInFileName += "Haddplots.root";
            firstSampInFile = new TFile(addPath + firstSampInFileName);
            if (!(firstSampFileStream->eof())) {
                firstSampFileStream->getline(Buffer, 500);
                singPlot2MakeFirstFile += TString(string(Buffer));
                cout << "singPlot2MakeFirstFile from second line " << singPlot2MakeFirstFile << endl;
                if (!(firstSampFileStream->eof()) && !diffFilesSingSampCompare) {
                    firstSampFileStream->getline(Buffer, 500);
                    singPlot2MakeSecondFile += TString(string(Buffer));
                    cout << "singPlot2MakeSecondFile from third line " << singPlot2MakeSecondFile << endl;
                }
                if (singPlot2MakeFirstFile == "") doAllPlots = true;
            }
            else {
                doAllPlots = true;
            }
        }
    }
    if (diffFilesSingSampCompare) {
        secondSampFileStream = new ifstream(secondSampSetupFile);
        if (!(secondSampFileStream->eof())) {
            secondSampFileStream->getline(Buffer, 500);
            secondSampInFileBaseName += TString(string(Buffer));
            secondSampInFileName += TString(string(Buffer));
            secondSampInFileName += (doPURW) ? "_PURW": "";
            secondSampInFileName += (doSyst) ? "_wSyst": "";
            secondSampInFileName += "Haddplots.root";
            secondSampInFile = new TFile(addPath + secondSampInFileName);
            if (!(secondSampFileStream->eof())) {
                secondSampFileStream->getline(Buffer, 500);
                singPlot2MakeSecondFile += TString(string(Buffer));
            }
        }
    }
    if (singPlot2MakeFirstFile != singPlot2MakeSecondFile) {
        cout << "###############" << endl;
        cout << "###############" << endl;
        cout << "plots that are being compared are not the same -- do you want this????" << endl;
        cout << "###############" << endl;
        cout << "###############" << endl;
    }
    TString fracRatioNumerSingSampName = (doAbsRatio) ? "1st" : "(2nd - 1st)";
    TString fracRatioDenomSingSampName = (doAbsRatio) ? "2nd" : "1st";
    TString canvNameAdd = "";
    canvNameAdd += firstSampInFileBaseName;
    if (diffFilesSingSampCompare) {
        canvNameAdd += "_vs_";
        canvNameAdd += secondSampInFileBaseName;
    }
    else {
        //nothing in here for now
    }
    
    if (multiChannelAdd) {
        cout << " testing " << endl;
        multiChannelStream = new ifstream(multiChannelSetupFile);
        if (!(multiChannelStream->eof())) {
            multiChannelStream->getline(Buffer, 500);
            cout << " first line " << TString(string(Buffer)) << endl;
            multiChannelCompOutputName = TString(string(Buffer));
            while (!(multiChannelStream->eof())) {
                multiChannelStream->getline(Buffer, 500);
                cout << " next line " << TString(string(Buffer)) << endl;
                if (TString(string(Buffer)).EqualTo("")) {
                    continue;
                }
                else {
                    multiChannelIDs->push_back(ChannelID(subSampVec, TString(Buffer)));
                }
            }
        }
        else {
            cout << "Input MultiChannelAdd file is EMPTY!!! " << endl;
        }
    }
    ///Data-Driven systematics stuff
    //Oviedo EMu: .912598, EE: 0.658692, MuMu: .884796, FullCut: .863585
    float TTBarFullCutSFOviReReco[3] = {1.01276, 0.933961, 1.04784};
    float TTBarFullCutSFOvi[3] = {0.946965, 0.873284, 0.969927};
    float TTBarFullCutSFDESY[3] = {0.84936, 0.843207, 0.949907};
    if (!doExcSamps) {
        TTBarFullCutSFOvi[0] = 0.891653;
        TTBarFullCutSFOviReReco[0] = 0.953711;
    }
    float TTBarSF, TTBarSF_Other;
    if (whichNTuple == 0) {
        if (doReReco) {
            TTBarSF = TTBarFullCutSFOviReReco[whichTTbarGen];
        }
        else {
            TTBarSF = TTBarFullCutSFOvi[whichTTbarGen];            
        }
        TTBarSF_Other = TTBarFullCutSFDESY[whichTTbarGen];
    }
    else {
        TTBarSF = TTBarFullCutSFDESY[whichTTbarGen];
        TTBarSF_Other = TTBarFullCutSFOvi[whichTTbarGen];
    }
    if (useDDEstimate) {
        cout << "TTBarSF used " << TTBarSF << endl;
    }
    
    //plotting stuff        
    //    double KFactor = 1;
    //    bool KFactorScale =      1;
    //    bool DataIntegralScale = 0;
    bool doSystCurrPlot;
    float intLumi;
    float indLumiDESY[4] = {892, 4404, 7032, 7274};
    float indLumiOvi[4] = {892, 4404, 7032, 7274};
    //  intLumi = 19602.901;
    float indLumiOviReReco[4] = {876, 4404, 7016, 7360};  
    float nominalLumi = indLumiOviReReco[0] + indLumiOviReReco[1] + indLumiOviReReco[2] + indLumiOviReReco[3];
    if (whichNTuple == 0) {
        if (!doReReco) {
            intLumi = indLumiOvi[0] + indLumiOvi[1] + indLumiOvi[2] + indLumiOvi[3];
        }
        else {
            intLumi = indLumiOviReReco[0] + indLumiOviReReco[1] + indLumiOviReReco[2] + indLumiOviReReco[3];
        }
    }
    else {
        intLumi = indLumiDESY[0] + indLumiDESY[1] + indLumiDESY[2] + indLumiDESY[3];
    }
    float scaleLumi = intLumi / nominalLumi;
    cout << "scaleLumi " << scaleLumi << endl;
    vector<float> * signalSkimScaleVec;
    if (doSignal) signalSkimScaleVec = SignalSkimEfficiencyCalc(typeSMS, prefixT2tt, vecStopMassGrab, vecChi0MassGrab, vecCharginoMassGrab, intLumi);
    
    TString mcplot, mcplot_preRW, dataplot;
    TString plotVarName, subSampName;;
    TString dataLegendComp = "Data";
    THStack * mcStack;
    TString mcStackName;
    unsigned int num1DPlots = histVec1DtoUse->size();
    unsigned int num2DPlots = histVec_2D->size(); 
    unsigned int num3DPlots = histVec_3D->size();
    unsigned int numIsoPlots = vecIsoPlotNames->size();
    TString saveNameAddition("../Plots/");
    
    bool doOverflow[num1DPlots];
    bool doUnderflow[num1DPlots];
    int  RBNX[num1DPlots];
    int  RBNY[num1DPlots];
    int  RBNZ[num1DPlots];
    float XaxisLegendPos[num1DPlots];
    float YaxisLegendStart[num1DPlots];
    TLegend * leg;
    for (unsigned int i = 0; i < num1DPlots; ++i) {
        doOverflow[i] = false;
        doUnderflow[i] = false;
        RBNX[i] = 1;
        RBNY[i] = 1;
        RBNZ[i] = 1;
        XaxisLegendPos[i] = 0.7;
        YaxisLegendStart[i] = .83;
    }
    float YAxisLB = 0.2;
    float YAxisUBBase = 2E5;
    float YAxisUB, histMax;
    bool  logYPad1;
    //    TCanvas * c_Var[num1DPlots];
    TCanvas * c_Var;
    
    TString canvName, systCanvName;
    logYPad1 = true;
    if (doOneDee) {
        if (singSampCompare) {
            histOne1D = new TH1F();
            histTwo1D = new TH1F();
            histOne1DSystVec = new vector<TH1F *>;
            histTwo1DSystVec = new vector<TH1F *>;
            histOneCDF1DSystVec = new vector<TH1F *>;
            histTwoCDF1DSystVec = new vector<TH1F *>;
            fracRatioSystVec = new vector<TGraphAsymmErrors *>;
            errCompSpecSource = new vector<TGraphAsymmErrors *>;
            errCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
            systCanvNameVec = new vector<TString>;
            if (!multiChannelAdd) {
                subSampName = "";
                subSampName += subSampVec->at(whichChan).histNameSuffix;                
            }
            if (doAllPlots) {
                for (unsigned int k = 0; k < num1DPlots; ++k) {
                    plotVarName = "";
                    plotVarName += histVec1DtoUse->at(k).name;
                    plotVarName.Replace(0, 2, "");                                
                    histOneplot = histVec1DtoUse->at(k).name;
                    if (!multiChannelAdd) histOneplot += subSampName;                
                    if (multHistsSingSampCompare) {
                        histTwoplot = histVec1DtoUse->at(k).name;
                        if (!multiChannelAdd) histTwoplot += subSampName;
                    }
                    else {
                        histTwoplot = "";
                    }
                    canvName = "c_";
                    canvName += plotVarName;
                    if (multiChannelAdd) {
                        canvName += multiChannelCompOutputName;
                    }
                    else {
                        canvName += subSampName;   
                    }
                    canvName += canvNameAdd;
                    c_Var = new TCanvas(canvName, canvName, wtopx, wtopy, W_, H_);
                    doSystCurrPlot = (doSyst && histVec1DtoUse->at(k).doXSyst);
                    SingSampCompHistogramGrabber(firstSampInFile, histOne1D, histOne1DSystVec, histOneplot, secondSampInFile, histTwo1D, histTwo1DSystVec, histTwoplot, subSampName, RBNX[k], RBNY[k], RBNZ[k], doOverflow[k], doUnderflow[k], doSystCurrPlot);
                    histMax = (multHistsSingSampCompare) ? TMath::Max(histOne1D->GetMaximum(), histTwo1D->GetMaximum()) : histOne1D->GetMaximum();
                    YAxisUB = (YAxisUBBase > 1) ? histMax * 2.5E2 : histMax * 2.5;
                    HistMainAttSet(histOne1D, kWhite, 0, histColors[0], 2, histColors[0], 20, 0.9);
                    histOneIntegral = histOne1D->Integral();
                    if (multHistsSingSampCompare) {
                        HistMainAttSet(histTwo1D, histColors[1], 1001, histColors[1], 2, kWhite, 0, 0);
                        histTwoIntegral = histTwo1D->Integral();
                        histTwo1D->Scale(histOneIntegral/histTwoIntegral);
                        h_FracratioComp = FracRatioHist(histOne1D, histTwo1D, fracRatioNumerSingSampName, fracRatioDenomSingSampName, doAbsRatio, histOne1D->GetName() + TString("ratioComp"), fracRatioYAxisRange);
                    }
                    else {
                        h_FracratioComp = NULL;   
                    }                    
                    SpectrumDrawSingSampCompare(c_Var, histOne1D, firstSampInFileBaseName, histTwo1D, secondSampInFileBaseName, h_FracratioComp, XaxisLegendPos[k] - 0.2, YaxisLegendStart[k], YAxisLB, YAxisUB, logYPad1, doStats, "", intLumi);
                }
            }
            else {
                
                /// not commissioned yet! (7/23/13)
                doSystCurrPlot = true;
                
                histOneplot = singPlot2MakeFirstFile;
                if (!multiChannelAdd) histOneplot += subSampName;
                
                plotVarName = "";
                plotVarName += singPlot2MakeFirstFile;
                plotVarName.Replace(0, 2, "");
                if (multHistsSingSampCompare) {
                    addPlotVarName = "";
                    addPlotVarName += singPlot2MakeSecondFile;
                    addPlotVarName.Replace(0, 2, "");
                    if (addPlotVarName != plotVarName) {
                        plotVarName += "_vs_";
                        plotVarName += addPlotVarName;
                    }
                    histTwoplot = singPlot2MakeSecondFile;
                    if (!multiChannelAdd) histTwoplot += subSampName;
                }
                else {
                    histTwoplot = "";
                }
                canvName = "c_";
                canvName += plotVarName;
                if (multiChannelAdd) {
                    canvName += multiChannelCompOutputName;
                }
                else {
                    canvName += subSampName;   
                }
                c_Var = new TCanvas(canvName, canvName, wtopx, wtopy, W_, H_); 
                SingSampCompHistogramGrabber(firstSampInFile, histOne1D, histOne1DSystVec, histOneplot, secondSampInFile, histTwo1D, histTwo1DSystVec, histTwoplot, subSampName, 1, 1, 1, false, false, doSystCurrPlot);
                //                SpectrumDrawSingSampCompare(c_Var, histOne1D, firstSampInFileBaseName, histTwo1D, secondSampInFileBaseName, h_FracratioComp, XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, logYPad1, doStats, "", intLumi);
            }
            //multHistsSingSampCompare
        }
        else {
            for (unsigned int k = 0; k < num1DPlots; ++k) {   
                dataHist1DVec = new vector<TH1F *>;
                mcIndHist1DCentValVec = new vector<TH1F *>;
                mcCompHist1DCentValVec = new vector<TH1F *>;
                mcCompHist1DSystVec = new vector<TH1F *>;
                fracRatioSystVec = new vector<TGraphAsymmErrors *>;
                mcCDFSystVec = new vector<TH1F *>;
                errCompSpecSource = new vector<TGraphAsymmErrors *>;
                errCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
                systCanvNameVec = new vector<TString>;
                dataHistTH1Vec = new vector<TH1 *>;
                mcHistTH1Vec = new vector<TH1 *>;
                mcHistSystTH1Vec = new vector<TH1 *>;
                
                if (doOviDESYFracRatio) {
                    dataHistTH1Vec_Other = new vector<TH1 *>;
                    mcHistTH1Vec_Other = new vector<TH1 *>;
                    mcCompHist1DCentValVec_Other = new vector<TH1F *>;
                }
                
                plotVarName = "";
                plotVarName += histVec1DtoUse->at(k).name;
                plotVarName.Replace(0, 2, "");
                if (multiChannelAdd) {
                    cout << "multiChannelCompOutputName " << multiChannelCompOutputName << endl;
                }
                else {
                    subSampName = "";
                    subSampName += subSampVec->at(whichChan).histNameSuffix;
                    
                }
                cout << "Doing " << plotVarName << endl;
                dataplot = histVec1DtoUse->at(k).name;
                mcplot = histVec1DtoUse->at(k).name;                    
                if (SmearedPlots) {
                    dataplot.Replace(2,5,"");
                }
                canvName = "c_";
                canvName += plotVarName;
                mcStackName = "mcStack_";
                mcStackName += plotVarName;
                if (multiChannelAdd) {
                    canvName += multiChannelCompOutputName;
                    mcStackName += multiChannelCompOutputName;
                }
                else {
                    dataplot += subSampName;
                    canvName += subSampName;
                    mcStackName += subSampName;                    	
                }
                cout << "dataplot " << dataplot << endl;
                canvName += canvSuffixSaveName;
                
                c_Var = new TCanvas(canvName, canvName, wtopx, wtopy, W_, H_);
                mcStack = new THStack(mcStackName, "");
                doSystCurrPlot = (doSyst && histVec1DtoUse->at(k).doXSyst);
                if (multiChannelAdd) {
                    HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, dataplot, subSampVec, multiChannelIDs, useDDEstimate, TTBarSF, scaleLumi); 
                    cout << " test 22" << endl;
                    HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, mcplot, subSampVec, multiChannelIDs, useDDEstimate, TTBarSF, scaleLumi);
                    cout << " test 24" << endl;
                }
                else {
                    cout << "testing 1b " << endl;
                    HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, dataplot, subSampName, useDDEstimate, TTBarSF, scaleLumi);
                    cout << "testing 1b " << endl;
                    HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, mcplot, subSampName, useDDEstimate, TTBarSF, scaleLumi);
                }
                HistogramAdderData(dataHistTH1Vec, h_DataComp, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
//                cout << " grabbed mc central val " << endl;
                HistogramAdderMC(mcHistTH1Vec, mcCompHist1DCentValVec, sampleStartPositions, sampleAddNames, h_MCComp, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
                                cout << "testing 1d" << endl;
                mcCompHist1DCentValVecSorted = sortedVector(mcCompHist1DCentValVec, sortIndexVec);
                                cout << "testing 1e" << endl;
//                cout << " added mc central val " << endl;
                h_FracratioComp = FracRatioHist(h_DataComp, h_MCComp, fracRatioNumerName, fracRatioDenomName, doAbsRatio, "ratioComp", fracRatioYAxisRange);
                
                if (doOviDESYFracRatio) {
                    if (multiChannelAdd) {
                        HistogramVecGrabberCentValGrab(inputFiles_Other, true, dataHistTH1Vec_Other, nVtxBackScaleVec_Other, dataplot, subSampVec, multiChannelIDs, useDDEstimate, TTBarSF_Other, scaleLumi);
                        cout << " test 33" << endl;
                        HistogramVecGrabberCentValGrab(inputFiles_Other, false, mcHistTH1Vec_Other, nVtxBackScaleVec_Other, mcplot, subSampVec, multiChannelIDs, useDDEstimate, TTBarSF_Other, scaleLumi);
                        cout << " test 34" << endl;
                    }
                    else {
                        HistogramVecGrabberCentValGrab(inputFiles_Other, true, dataHistTH1Vec_Other, nVtxBackScaleVec_Other, dataplot, subSampName, useDDEstimate, TTBarSF_Other, scaleLumi);
                        HistogramVecGrabberCentValGrab(inputFiles_Other, false, mcHistTH1Vec_Other, nVtxBackScaleVec_Other, mcplot, subSampName, useDDEstimate, TTBarSF_Other, scaleLumi);
                    }
                    HistogramAdderData(dataHistTH1Vec_Other, h_DataComp_Other, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
                    HistogramAdderMC(mcHistTH1Vec_Other, mcCompHist1DCentValVec_Other, sampleStartPositions_Other, sampleAddNames_Other, h_MCComp_Other, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
                    h_FracratioComp_Other = FracRatioHist(h_DataComp_Other, h_MCComp_Other, fracRatioNumerName, fracRatioDenomName, doAbsRatio, "ratioComp_Other", fracRatioYAxisRange);
                }
                if (doSystCurrPlot) {
                    if (multiChannelAdd) {
                        HistogramVecGrabberSystGrab(inputFiles, mcHistTH1Vec, mcHistSystTH1Vec, nVtxBackScaleVec, mcplot, subSampVec, multiChannelIDs, systVec, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);            
                    }
                    else {
                        HistogramVecGrabberSystGrab(inputFiles, mcHistTH1Vec, mcHistSystTH1Vec, nVtxBackScaleVec, mcplot, subSampName, systVec, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);               
                    }
                    HistogramProjectorSyst(mcHistSystTH1Vec, mcCompHist1DSystVec, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
                }
                if (doSignal) {
                    vecStop1DCentValHists = new vector<TH1F *>;
                    vecErrSigCompVecSpecSource = new vector<vector<TGraphAsymmErrors *> *>;
                    vecErrSigCompVecSpecSource_pStat = new vector<vector<TGraphAsymmErrors *> *>;
                    for (unsigned int iSigPoints = 0; iSigPoints < vecStopMassGrab->size(); ++iSigPoints) {
                        vecCurrSignalSystTH1Hists = new vector<TH1 *>;
                        vecCurrStop1DSystHists = new vector<TH1F *>;
                        vecCurrSigCompSpecSource = new vector<TGraphAsymmErrors *>;
                        vecCurrSigCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
                        HistogramVecGrabber_Signal(inputFilesSignal, signalSkimScaleVec, iSigPoints, currSignalCentValTH1Hist, mcplot, vecCurrSignalSystTH1Hists, systVec, subSampName, doOverflow[k], doUnderflow[k], doSystCurrPlot, allMT2llSystematic);
                        HistogramAdderSignal(currSignalCentValTH1Hist, currSignal1DCentValHist, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
                        HistMainAttSet(currSignal1DCentValHist, kWhite, 0, mcColorsSignal->at(iSigPoints), 2, kWhite, 0, 0, mcStylesSignal->at(iSigPoints));
                        vecStop1DCentValHists->push_back(currSignal1DCentValHist);                    
                        if (doSystCurrPlot) {
                            HistogramProjectorSyst(vecCurrSignalSystTH1Hists, vecCurrStop1DSystHists, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
                            /*
                            cout << "currSignal1DCentValHist " << currSignal1DCentValHist->GetName() << endl;
                            cout << "vecCurrStop1DSystHists size " << vecCurrStop1DSystHists->size() << endl;
                            */
                            SystGraphMakers(currSignal1DCentValHist, vecCurrStop1DSystHists, vecCurrSigCompSpecSource, vecCurrSigCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, mcColorsSignal->at(iSigPoints), plotVarName, doAbsRatio, fracRatioYAxisRange, doSymErr, true, SmearedPlots);
                            vecErrSigCompVecSpecSource->push_back(vecCurrSigCompSpecSource);
                            vecErrSigCompVecSpecSource_pStat->push_back(vecCurrSigCompSpecSource_pStat);
                        }
                        
                    }            
                }
                h_ErrComp = (TH1F *) h_MCComp->Clone();           
                HistMainAttSet(h_DataComp, kWhite, 0, kBlack, 2, kBlack, 20, 0.9);
                HistMainAttSet(h_MCComp, kWhite, 0, kWhite, 2, kWhite, 0, 0);
                HistMainAttSet(h_ErrComp, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0);
                histMax = TMath::Max(h_DataComp->GetMaximum(), h_MCComp->GetMaximum());
                YAxisUB = (YAxisUBBase > 1) ? histMax * 2.5E2 : histMax * 2.5;
                for (unsigned int j = 0; j < mcCompHist1DCentValVecSorted->size(); ++j) {
                    HistMainAttSet(mcCompHist1DCentValVecSorted->at(j), mcColors->at(j), 1001, mcColors->at(j), 2, kWhite, 0, 0);
                    mcStack->Add(mcCompHist1DCentValVecSorted->at(j));
                    cout << "integral for mcCompCentVal " << mcCompHist1DCentValVecSorted->at(j)->GetName() << " is " << mcCompHist1DCentValVecSorted->at(j)->Integral() << endl;
                }   
                cout << "integral for DataComp " << h_DataComp->GetName() << " is " << h_DataComp->Integral() << endl;            
                cout << "integral for MCComp " << h_MCComp->GetName() << " is " << h_MCComp->Integral() << endl;
//                cout << "mcCompHist1DSystVec->size() " << mcCompHist1DSystVec->size() << endl;
                for (unsigned int kMC = 0; kMC < mcCompHist1DSystVec->size(); ++kMC) {
                    cout << "integral for syst ";
                    cout << mcCompHist1DSystVec->at(kMC)->GetName() << " is " << mcCompHist1DSystVec->at(kMC)->Integral() << endl;
                }
                if (mcplot.Contains("h_MT2llControl") && calcTTBarNorm) {
                    TTBarSF = DataDrivenTTBarScaleFactor(h_DataComp, h_MCComp, mcCompHist1DCentValVecSorted); 
                    cout << "For TTBarGen " << whichTTbarGen << " and nTuple " << whichNTuple << ", calculated TTBar SF is " << TTBarSF << endl;
                }
                if (plotVarName == "MT2ll" ) {
                    cout << "DataIntegral for MT2ll > " << h_DataComp->GetBinLowEdge(21) << " is " << h_DataComp->Integral(21, h_DataComp->GetNbinsX()+1) << endl;
                    cout << "MCIntegral for MT2ll > " << h_MCComp->GetBinLowEdge(21) << " is " << h_MCComp->Integral(21, h_DataComp->GetNbinsX()+1) << endl;
                }
                SpectrumDraw(c_Var, h_DataComp, dataLegendComp, h_MCComp, h_FracratioComp, h_ErrComp, mcStack, XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, logYPad1, mcLegends, mcCompHist1DCentValVecSorted, doStats, "", intLumi, leg, doSmartFracRatio, doOviDESYFracRatio, h_FracratioComp_Other, showLegend);
                if (doSignal) {
                    SpectrumDraw_AddSignal(c_Var, vecStop1DCentValHists, mcLegendsSignal, leg, showLegend);
                }
                c_Var->SaveAs(saveNameAddition + canvName + TString(".pdf"));
                if (saveDotCFile) c_Var->SaveAs(saveNameAddition + canvName + TString(".C"));
                
                if (doSystCurrPlot) {
                    
                    SystGraphMakers(h_MCComp, mcCompHist1DSystVec, errCompSpecSource, errCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, kGray + 1, plotVarName, doAbsRatio, fracRatioYAxisRange, doSymErr, false, SmearedPlots);
                    for (unsigned int iSyst = 0; iSyst < systCanvNameVec->size(); ++iSyst) {
                        systCanvName = canvName + systCanvNameVec->at(iSyst);
                        c_Var = new TCanvas(systCanvName, systCanvName, wtopx, wtopy, W_, H_);
                        SpectrumDrawSyst(c_Var, h_DataComp, dataLegendComp, h_MCComp, mcStack, errCompSpecSource_pStat->at(iSyst), errCompSpecSource->at(iSyst), h_FracratioComp, fracRatioSystVec->at(iSyst), XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, logYPad1, mcLegends, mcCompHist1DCentValVecSorted, doStats, "", intLumi, leg, doSmartFracRatio, doOviDESYFracRatio, h_FracratioComp_Other, showLegend);
                        
                        if (doSignal) {
                            SpectrumDrawSyst_AddSignal(c_Var, vecStop1DCentValHists, vecErrSigCompVecSpecSource_pStat, iSyst, mcLegendsSignal, leg, showLegend);
                        }
                        c_Var->SaveAs(saveNameAddition + systCanvName + TString(".pdf"));
                        if (saveDotCFile) c_Var->SaveAs(saveNameAddition + systCanvName + TString(".C"));
                    }
                }
                
                //        cout << "Data Integral " << h_DataComp->Integral() << endl;
                //        HistMainAttSet(h_StopComp, mcColors[5], 1001, mcColors[5], 2, kWhite, 0, 0); 
                /*
                 if (whichNTuple == 0) {
                 DataMCSpectrumDrawSplit(c_Var[k], h_DataComp, mcStack, h_MCComp, h_ErrComp, h_FracratioComp, XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, dataLegendComp, mcLegendsOvi, mcPlotHistsOvi, endJ, doAbsRatio, doStats, "", intLumi[whichNTuple]);
                 }
                 else if (whichNTuple == 1) {
                 DataMCSpectrumDrawSplit(c_Var[k], h_DataComp, mcStack, h_MCComp, h_ErrComp, h_FracratioComp, XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, dataLegendComp, mcLegendsDESY, mcPlotHistsDESY, endJ, doAbsRatio, doStats, "", intLumi[whichNTuple]);
                 }
                 c_Var[k]->SaveAs(canvName + TString(".pdf"));
                 if (plotVar[k] == "MT2llControl" || plotVar[k] == "MT2ll") {
                 doCDF = 1;
                 }
                 else {
                 doCDF = 0;
                 }
                 if (plotVar[k] == "MT2llControl") {
                 reportChi2 = 1;
                 }
                 else {
                 reportChi2 = 0;
                 }
                 if (doCDF) {
                 c_Var[k] = new TCanvas(canvName + TString("_CDF_Data"), canvName + TString("_CDF_Data"), wtopx, wtopy, W_, H_);
                 CDFRight(h_DataComp, h_CDF_Data);
                 h_CDF_Data->Draw("hist");
                 c_Var[k]->SaveAs(canvName + TString("_CDF_Data") + TString(".pdf"));
                 c_Var[k] = new TCanvas(canvName + TString("_CDF_Sim"), canvName + TString("_CDF_Sim"), wtopx, wtopy, W_, H_);
                 CDFRight(h_MCComp, h_CDF_MC);
                 h_CDF_MC->Draw("hist");
                 c_Var[k]->SaveAs(canvName + TString("_CDF_Sim") + TString(".pdf"));
                 }
                 if (reportChi2) {
                 int N = h_DataComp->GetNbinsX();
                 Double_t res[20];
                 cout << "Chi2 WW for histogram " << h_DataComp->GetName() << " is " << h_DataComp->Chi2Test(h_MCComp, "WW P", res) << endl;
                 cout << "Chi2 UW for histogram " << h_DataComp->GetName() << " is " << h_DataComp->Chi2Test(h_MCComp, "UW P", res) << endl;
                 }
                 */
                delete dataHistTH1Vec;
                delete mcHistTH1Vec;
                delete dataHist1DVec;
                delete mcIndHist1DCentValVec;
                delete mcCompHist1DCentValVecSorted;
                delete mcCompHist1DCentValVec;
                delete mcCompHist1DSystVec;
                delete fracRatioSystVec;
                delete mcCDFSystVec;
                delete errCompSpecSource;
                delete errCompSpecSource_pStat;
                delete systCanvNameVec;
                if (doOviDESYFracRatio) {
                    delete dataHistTH1Vec_Other;
                    delete mcHistTH1Vec_Other;
                }                                                
            }
        }
    }
    
    TString cutString, caseCanvName, caseStackName;
    int nCases;
    
    if (doTwoDee) {    
        nCases = 3;
        for (unsigned int k2D = 0; k2D < num2DPlots; ++k2D) {
            dataHist2DVec = new vector<TH2F *>;
            mcIndHist2DCentValVec = new vector<TH2F *>;
            mcCompHist2DCentValVec = new vector<TH2F *>;
            mcCompHist2DSystVec = new vector<TH2F *>;
            
            fracRatioSystVec = new vector<TGraphAsymmErrors *>;
            mcCDFSystVec = new vector<TH1F *>;
            errCompSpecSource = new vector<TGraphAsymmErrors *>;
            errCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
            systCanvNameVec = new vector<TString>;
            
            plotVarName = "";
            plotVarName += histVec_2D->at(k2D).name;
            plotVarName.Replace(0, 2, "");
            subSampName = "";
            subSampName += subSampVec->at(whichChan).histNameSuffix;
            cout << "subSampName " << subSampName << endl;
            cout << "Doing " << plotVarName << endl;
            dataplot = histVec_2D->at(k2D).name;
            dataplot += subSampName;
            cout << "dataplot " << dataplot << endl;
            mcplot = histVec_2D->at(k2D).name;
            canvName = "c_";
            canvName += plotVarName;
            canvName += subSampName;
            canvName += canvSuffixSaveName;
            mcStackName = "mcStack_";
            mcStackName += plotVarName;
            mcStackName += subSampName;
            
            //            doSystCurrPlot = (doSyst && (histVec_2D->at(k2D).doXSyst || histVec_2D->at(k2D).doYSyst));
            doSystCurrPlot = false;
            HistogramVecGrabber(inputFiles, dataHist2DVec, mcIndHist2DCentValVec, mcCompHist2DSystVec, nVtxBackScaleVec, systVec, dataplot, mcplot, subSampName, RBNX[0], RBNY[0], RBNZ[0], doOverflow[0], doUnderflow[0], doSystCurrPlot, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);
            for (int iCase = 0; iCase < nCases; ++iCase) {
                dataHist1DVec = new vector<TH1F *>;
                mcIndHist1DCentValVec = new vector<TH1F *>;
                mcCompHist1DCentValVec = new vector<TH1F *>;                        
                mcCompHist1DSystVec = new vector<TH1F *>;
                
                caseCanvName = canvName;
                caseCanvName += "_case";
                caseCanvName += iCase;   
                caseStackName = mcStackName;
                caseStackName += "_case";
                caseStackName += iCase;
                cutString = HistProjection1D(dataHist2DVec, dataHist1DVec, dataplot, iCase);
                cutString = HistProjection1D(mcIndHist2DCentValVec, mcIndHist1DCentValVec, mcplot, iCase);
                //                cutString = HistProjection1D(mcCompHist2DCentValVec, mcCompHist1DCentValVec, mcplot, iCase);
                cutString = HistProjection1D(mcCompHist2DSystVec, mcCompHist1DSystVec, mcplot, iCase);
                //                cout << "cutString: " << cutString << endl;
                if (cutString.Contains("Nothin")) continue;      
                mcStack = new THStack(caseStackName, "");
                c_Var = new TCanvas(caseCanvName, caseCanvName, wtopx, wtopy, W_, H_);          
                HistogramAdderSyst(dataHist1DVec, mcIndHist1DCentValVec, mcCompHist1DCentValVec, h_DataComp, h_MCComp, h_FracratioComp, whichNTuple, doAbsRatio, fracRatioYAxisRange);
                h_ErrComp = (TH1F *) h_MCComp->Clone();
                errCompStatCentVal = clonePoints(h_ErrComp);
                HistMainAttSet(h_DataComp, kWhite, 0, kBlack, 2, kBlack, 20, 0.9);
                //        HistMainAttSet(h_MCComp, kWhite, 0, kRed, 2, kWhite, 0, 0);
                HistMainAttSet(h_MCComp, kWhite, 0, kWhite, 2, kWhite, 0, 0);
                HistMainAttSet(h_ErrComp, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0);
                histMax = TMath::Max(h_DataComp->GetMaximum(), h_MCComp->GetMaximum());
                YAxisUB = (YAxisUBBase > 1) ? histMax * 2.5E2 : histMax * 2.5;
                
                for (unsigned int j = 0; j < mcCompHist1DCentValVec->size(); ++j) {
                    HistMainAttSet(mcCompHist1DCentValVec->at(j), mcColors->at(j), 1001, mcColors->at(j), 2, kWhite, 0, 0);
                    mcStack->Add(mcCompHist1DCentValVec->at(j));
                    cout << "integral for mcCompCentVal " << mcCompHist1DCentValVec->at(j)->GetName() << " is " << mcCompHist1DCentValVec->at(j)->Integral() << endl;
                }   
                cout << "integral for MCComp " << h_MCComp->GetName() << " is " << h_MCComp->Integral() << endl;
                for (unsigned int kMC = 0; kMC < mcCompHist1DSystVec->size(); ++kMC) {
                    cout << "integral for syst " << mcCompHist1DSystVec->at(kMC)->GetName() << " is " << mcCompHist1DSystVec->at(kMC)->Integral() << endl;
                }
                SpectrumDraw(c_Var, h_DataComp, dataLegendComp, h_MCComp, h_FracratioComp, h_ErrComp, mcStack, XaxisLegendPos[0], YaxisLegendStart[0], YAxisLB, YAxisUB, logYPad1, mcLegends, mcCompHist1DCentValVec, doStats, cutString, intLumi, leg, doSmartFracRatio, doOviDESYFracRatio, h_FracratioComp_Other, showLegend);
                c_Var->SaveAs(saveNameAddition + caseCanvName + TString(".pdf"));
                if (saveDotCFile) c_Var->SaveAs(saveNameAddition + caseCanvName + TString(".C"));   
            }
            
            delete dataHist1DVec;
            delete dataHist2DVec;
            
            delete mcIndHist1DCentValVec;
            delete mcIndHist2DCentValVec;
            
            delete mcCompHist1DCentValVec;
            delete mcCompHist2DCentValVec;
            
            delete mcCompHist1DSystVec;
            delete mcCompHist2DSystVec;
            
            delete fracRatioSystVec;
            delete mcCDFSystVec;
            delete errCompSpecSource;
            delete errCompSpecSource_pStat;
            delete systCanvNameVec;
        }
    }
    if (doThreeDee) {
        nCases = 9;
        for (unsigned int k3D = 0; k3D < num3DPlots; ++k3D) {
            dataHist3DVec = new vector<TH3F *>;
            mcIndHist3DCentValVec = new vector<TH3F *>;
            mcCompHist3DCentValVec = new vector<TH3F *>;
            mcCompHist3DSystVec = new vector<TH3F *>;
            
            fracRatioSystVec = new vector<TGraphAsymmErrors *>;
            mcCDFSystVec = new vector<TH1F *>;
            errCompSpecSource = new vector<TGraphAsymmErrors *>;
            errCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
            systCanvNameVec = new vector<TString>;
            
            plotVarName = "";
            plotVarName += histVec_3D->at(k3D).name;
            plotVarName.Replace(0, 2, "");
            subSampName = "";
            subSampName += subSampVec->at(whichChan).histNameSuffix;
            cout << "subSampName " << subSampName << endl;
            cout << "Doing " << plotVarName << endl;
            dataplot = histVec_3D->at(k3D).name;
            dataplot += subSampName;
            cout << "dataplot " << dataplot << endl;
            mcplot = histVec_3D->at(k3D).name;
            canvName = "c_";
            canvName += plotVarName;
            canvName += subSampName;
            canvName += canvSuffixSaveName;
            mcStackName = "mcStack_";
            mcStackName += plotVarName;
            mcStackName += subSampName;
            
            //            doSystCurrPlot = (doSyst && (histVec_3D->at(k3D).doXSyst || histVec_3D->at(k3D).doYSyst || histVec_3D->at(k3D).doZSyst));
            doSystCurrPlot = false;
            HistogramVecGrabber(inputFiles, dataHist3DVec, mcIndHist3DCentValVec, mcCompHist3DSystVec, nVtxBackScaleVec, systVec, dataplot, mcplot, subSampName, RBNX[0], RBNY[0], RBNZ[0], doOverflow[0], doUnderflow[0], doSystCurrPlot, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);
            for (int iCase = 0; iCase < nCases; ++iCase) {
                dataHist1DVec = new vector<TH1F *>;
                mcIndHist1DCentValVec = new vector<TH1F *>;
                mcCompHist1DCentValVec = new vector<TH1F *>;                        
                mcCompHist1DSystVec = new vector<TH1F *>;
                
                dataHist2DVec = new vector<TH2F *>;
                mcIndHist2DCentValVec = new vector<TH2F *>;
                mcCompHist2DCentValVec = new vector<TH2F *>;
                mcCompHist2DSystVec = new vector<TH2F *>;
                
                caseCanvName = canvName;
                caseCanvName += "_case";
                caseCanvName += iCase;   
                caseStackName = mcStackName;
                caseStackName += "_case";
                caseStackName += iCase;
                cutString = HistProjection1D(dataHist3DVec, dataHist1DVec, dataplot, iCase);
                cutString = HistProjection1D(mcIndHist3DCentValVec, mcIndHist1DCentValVec, mcplot, iCase);
                //                cutString = HistProjection1D(mcCompHist3DCentValVec, mcCompHist1DCentValVec, mcplot, iCase);
                cutString = HistProjection1D(mcCompHist3DSystVec, mcCompHist1DSystVec, mcplot, iCase);
                //                cout << "cutString: " << cutString << endl;
                if (cutString.Contains("Nothin")) continue;      
                mcStack = new THStack(caseStackName, "");
                c_Var = new TCanvas(caseCanvName, caseCanvName, wtopx, wtopy, W_, H_);          
                HistogramAdderSyst(dataHist1DVec, mcIndHist1DCentValVec, mcCompHist1DCentValVec, h_DataComp, h_MCComp, h_FracratioComp, whichNTuple, doAbsRatio, fracRatioYAxisRange);
                h_ErrComp = (TH1F *) h_MCComp->Clone();
                errCompStatCentVal = clonePoints(h_ErrComp);
                HistMainAttSet(h_DataComp, kWhite, 0, kBlack, 2, kBlack, 20, 0.9);
                HistMainAttSet(h_MCComp, kWhite, 0, kWhite, 2, kWhite, 0, 0);
                HistMainAttSet(h_ErrComp, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0);
                histMax = TMath::Max(h_DataComp->GetMaximum(), h_MCComp->GetMaximum());
                YAxisUB = (YAxisUBBase > 1) ? histMax * 2.5E2 : histMax * 2.5;
                
                for (unsigned int j = 0; j < mcCompHist1DCentValVec->size(); ++j) {
                    HistMainAttSet(mcCompHist1DCentValVec->at(j), mcColors->at(j), 1001, mcColors->at(j), 2, kWhite, 0, 0);
                    mcStack->Add(mcCompHist1DCentValVec->at(j));
                    cout << "integral for mcCompCentVal " << mcCompHist1DCentValVec->at(j)->GetName() << " is " << mcCompHist1DCentValVec->at(j)->Integral() << endl;
                }   
                cout << "integral for MCComp " << h_MCComp->GetName() << " is " << h_MCComp->Integral() << endl;
                for (unsigned int kMC = 0; kMC < mcCompHist1DSystVec->size(); ++kMC) {
                    cout << "integral for syst " << mcCompHist1DSystVec->at(kMC)->GetName() << " is " << mcCompHist1DSystVec->at(kMC)->Integral() << endl;
                }
                SpectrumDraw(c_Var, h_DataComp, dataLegendComp, h_MCComp, h_FracratioComp, h_ErrComp, mcStack, XaxisLegendPos[0], YaxisLegendStart[0], YAxisLB, YAxisUB, logYPad1, mcLegends, mcCompHist1DCentValVec, doStats, cutString, intLumi, leg, doSmartFracRatio, doOviDESYFracRatio, h_FracratioComp_Other, showLegend);
                c_Var->SaveAs(saveNameAddition + caseCanvName + TString(".pdf"));
                if (saveDotCFile) c_Var->SaveAs(saveNameAddition + caseCanvName + TString(".C"));
            }
            
            delete dataHist1DVec;
            delete dataHist2DVec;
            delete dataHist3DVec;
            
            delete mcIndHist1DCentValVec;
            delete mcIndHist2DCentValVec;
            delete mcIndHist3DCentValVec;
            
            delete mcCompHist1DCentValVec;
            delete mcCompHist2DCentValVec;
            delete mcCompHist3DCentValVec;
            
            delete mcCompHist1DSystVec;
            delete mcCompHist2DSystVec;
            delete mcCompHist3DSystVec;
            
            delete fracRatioSystVec;
            delete mcCDFSystVec;
            delete errCompSpecSource;
            delete errCompSpecSource_pStat;
            delete systCanvNameVec;
        }   
    }
    
    if (doIsoPlots) { 
        for (unsigned int k = 0; k < numIsoPlots; ++k) {
            dataHist1DVec = new vector<TH1F *>;
            mcIndHist1DCentValVec = new vector<TH1F *>;
            mcCompHist1DCentValVec = new vector<TH1F *>;
            dataHistTH1Vec = new vector<TH1 *>;
            mcHistTH1Vec = new vector<TH1 *>;
            plotVarName = "";
            plotVarName += vecIsoPlotNames->at(k);
            plotVarName.Replace(0, 2, "");
            subSampName = "";
            dataplot = vecIsoPlotNames->at(k);
            mcplot = vecIsoPlotNames->at(k);
//            cout << "test Iso 1 " << dataplot <<  endl;
            canvName = "c_";
            canvName += plotVarName;
            canvName += canvSuffixSaveName;
            mcStackName = "mcStack_";
            mcStackName += plotVarName;
            c_Var = new TCanvas(canvName, canvName, wtopx, wtopy, W_, H_);
            mcStack = new THStack(mcStackName, "");
//            cout << "test Iso 2 " << inputFiles->size() << endl;
            HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, dataplot, subSampName, useDDEstimate, TTBarSF, scaleLumi);
//            cout << "test Iso 3 " << dataHistTH1Vec->size() << endl;
            HistogramAdderData(dataHistTH1Vec, h_DataComp, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
//            cout << "test Iso 4" << endl;
            HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, mcplot, subSampName, useDDEstimate, TTBarSF, scaleLumi);
//            cout << "test Iso 5 " << mcHistTH1Vec->size() << endl;
            HistogramAdderMC(mcHistTH1Vec, mcCompHist1DCentValVec, sampleStartPositions, sampleAddNames, h_MCComp, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
//            cout << "test Iso 6" << endl;
            mcCompHist1DCentValVecSorted = sortedVector(mcCompHist1DCentValVec, sortIndexVec);
//            cout << "test Iso 7" << endl;
            h_FracratioComp = FracRatioHist(h_DataComp, h_MCComp, fracRatioNumerName, fracRatioDenomName, doAbsRatio, "ratioComp", fracRatioYAxisRange);
            if (doSignal) {
                vecStop1DCentValHists = new vector<TH1F *>;
                vecErrSigCompVecSpecSource = new vector<vector<TGraphAsymmErrors *> *>;
                vecErrSigCompVecSpecSource_pStat = new vector<vector<TGraphAsymmErrors *> *>;
                for (unsigned int iSigPoints = 0; iSigPoints < vecStopMassGrab->size(); ++iSigPoints) {
                    vecCurrSignalSystTH1Hists = new vector<TH1 *>;
                    vecCurrStop1DSystHists = new vector<TH1F *>;
                    vecCurrSigCompSpecSource = new vector<TGraphAsymmErrors *>;
                    vecCurrSigCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
                    HistogramVecGrabber_Signal(inputFilesSignal, signalSkimScaleVec, iSigPoints, currSignalCentValTH1Hist, mcplot, vecCurrSignalSystTH1Hists, systVec, subSampName, doOverflow[k], doUnderflow[k], doSystCurrPlot, allMT2llSystematic);
                    HistogramAdderSignal(currSignalCentValTH1Hist, currSignal1DCentValHist, RBNX[k], 1, 1, -1, -1, -1, -1, "", doOverflow[k], doUnderflow[k], "");
                    HistMainAttSet(currSignal1DCentValHist, kWhite, 0, mcColorsSignal->at(iSigPoints), 2, kWhite, 0, 0, mcStylesSignal->at(iSigPoints));
                    vecStop1DCentValHists->push_back(currSignal1DCentValHist);
                }            
            }
            h_ErrComp = (TH1F *) h_MCComp->Clone();           
            HistMainAttSet(h_DataComp, kWhite, 0, kBlack, 2, kBlack, 20, 0.9);
            HistMainAttSet(h_MCComp, kWhite, 0, kWhite, 2, kWhite, 0, 0);
            HistMainAttSet(h_ErrComp, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0);
            histMax = TMath::Max(h_DataComp->GetMaximum(), h_MCComp->GetMaximum());
            YAxisUB = (YAxisUBBase > 1) ? histMax * 2.5E2 : histMax * 2.5;
            for (unsigned int j = 0; j < mcCompHist1DCentValVecSorted->size(); ++j) {
                HistMainAttSet(mcCompHist1DCentValVecSorted->at(j), mcColors->at(j), 1001, mcColors->at(j), 2, kWhite, 0, 0);
                mcStack->Add(mcCompHist1DCentValVecSorted->at(j));
            }
            SpectrumDraw(c_Var, h_DataComp, dataLegendComp, h_MCComp, h_FracratioComp, h_ErrComp, mcStack, XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, logYPad1, mcLegends, mcCompHist1DCentValVecSorted, doStats, "", intLumi, leg, doSmartFracRatio, doOviDESYFracRatio, h_FracratioComp_Other, showLegend);
            if (doSignal) {
                SpectrumDraw_AddSignal(c_Var, vecStop1DCentValHists, mcLegendsSignal, leg, showLegend);
            }
            c_Var->SaveAs(saveNameAddition + canvName + TString(".pdf"));
            if (saveDotCFile) c_Var->SaveAs(saveNameAddition + canvName + TString(".C"));
            delete dataHist1DVec;
            delete mcIndHist1DCentValVec;
            delete mcCompHist1DCentValVecSorted;
            delete mcCompHist1DCentValVec;
        }
    }    
    theApp.Run(retVal);
    //    theApp.Terminate(0);
}
