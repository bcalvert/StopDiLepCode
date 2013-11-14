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
    int whichChan     = 3;          // which "FullCut" channel to look at -- 0 is MuMu, 1 is EE, 2 is EMu, 3 is all together
    int whichNTuple   = 1;          //as with the plot making code, leave as 1 for now -- 0 is Oviedo, 1 is DESY    
    int whichTTbarGen = 0;          // 0 is Madgraph, 1 is MC@NLO, 2 is Powheg
    bool doExcSamps   = 0;          // For grabbing exclusive (DY + N Jets, TTBar Decay modes) or inclusive samples (As of 8/5/13, only applies to Oviedo)
    bool calcTTBarNorm = 0;         // calculate TTBar normalization by utilizing integral to data - (other backgrounds) for MT2ll < 80 in the "Full Cut region"
    bool doNonSig     = 0;          // For whether or not to grab SM background and data
    bool doSignal     = 0;          // For whether or not to grab a signal point
    bool doFOM        = 0;
    bool doOneDeeFOM     = 0;
    bool doTwoDeeFOM     = 0;
    bool doYieldV1    = 1;
    bool doYieldV2    = 0;
    bool doReReco     = 0; 
    bool SmearedPlots   = 0;
    int  JetsSmeared    = 0;   
    TString typeSMS   = "";         // Which type of SMS to grab -- either "T2tt" or "T2bw" (as of 8/5/13) only has T2tt FineBin
    TString prefixT2tt   = "";      // prefix for which kind of T2tt to grab
    vector<int> * vecStopMassGrab = new vector<int>;       // vector to hold the list of Stop masses to brab
    vector<int> * vecChi0MassGrab = new vector<int>;       // vector to hold the list of Chi0 masses to brab
    vector<int> * vecCharginoMassGrab = new vector<int>;   // vector to hold the list of Chargino masses to brab
    bool doPURW       = 0;          // grab the nVtx reweighted MC files
    bool doSyst       = 1;          // look at systematics plots
    bool addThings    = 1;          // Add together similar kinds of events (for aesthetic reasons) like VV backgrounds -- (6/25/13) don't turn off for now haven't validated code fully when not adding
    bool doSymErr     = 0;
    bool saveDotCFile = 0;          // If on, saves a .C version of every plot that is made
    bool useDDEstimate = 0;         // whether or not to use data-driven estimates for appropriate backgrounds -- as of right now it is just the TTBar norm to MT2ll < 80 GeV (7/22/13)
    bool allMT2llSystematic = 0;    // Whether or not to use the MT2ll systematic smearing for all MC or just TTBar
    int  versNumber     = 1;
    bool makePlots      = 1;
    bool dispIndSource = 0;
    bool printFOMInfo  = 0;
    int  typeFOM       = 0;
    TString FOMString = "";
    int  MT2llaxisCut = 80;
    int  MT2lbaxisCut = 170;
    bool UseUnblindedData = 0;
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
        else if (strncmp (argv[k],"DIS", 3) == 0) {
            dispIndSource = 1;
        }
        else if (strncmp (argv[k],"ReleaseTheKraken", 16) == 0) {
            UseUnblindedData = 1;
            cout << "RELEASING THE KRAKEN!!! " << endl;
            cout << "http://www.youtube.com/watch?v=gb2zIR2rvRQ " << endl;
        }
        else if (strncmp (argv[k],"pFOMI", 5) == 0) {
            printFOMInfo = 1;
        }
        else if (strncmp (argv[k],"MT2AxisCuts", 11) == 0) {
            MT2llaxisCut = strtol(argv[k+1], NULL, 10 );
            MT2lbaxisCut = strtol(argv[k+2], NULL, 10 );
        }
        else if (strncmp (argv[k],"noPlots", 7) == 0) {
            makePlots = 0;
        }
        else if (strncmp (argv[k],"versNum", 7) == 0) {
            versNumber = strtol(argv[k+1], NULL, 10 );
        }       
        else if (strncmp (argv[k],"doReReco", 8) == 0) {
            doReReco = 1;
        }
        else if (strncmp (argv[k],"doNonSig", 8) == 0) {
            doNonSig = 1;
        }        
        else if (strncmp (argv[k],"doFOM", 5) == 0) {
            doFOM = 1;
            typeFOM = strtol(argv[k+1], NULL, 10 );
            FOMString = typeFOM == 0 ? "S/sqrt(S+B)" : "S/sqrt(B)";
        }      
        else if (strncmp (argv[k],"doOneDeeFOM", 11) == 0) {
            doOneDeeFOM = 1;
        }       
        else if (strncmp (argv[k],"doTwoDeeFOM", 11) == 0) {
            doTwoDeeFOM = 1;
        }
        else if (strncmp (argv[k],"doYieldV2", 8) == 0) {
            doYieldV1 = 0;
            doYieldV2 = 1;
        }
        else if (strncmp (argv[k],"sDCF", 4) == 0) {
            saveDotCFile = 1;
        }
        else if (strncmp (argv[k],"JsSm", 4) == 0) {
            SmearedPlots = 1;
            JetsSmeared = strtol(argv[k+1], NULL, 10);
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
    }
    TRint theApp("App", &argc, argv);
    Bool_t retVal = kTRUE;
    
    if (!doExcSamps && whichNTuple == 1 ) {
        doExcSamps = true;
        cout << "setting do Exclusive Samples to true because running on DESY" << endl;
    }
    //Set up the file input
    //    vector<TFile *> * inFiles = new vector<TFile*>;
    vector<TString> * fileInNames = StopFileNames(whichNTuple, whichTTbarGen, doExcSamps, doReReco);
    vector<TFile *> * inputFiles  = StopFiles(whichNTuple, fileInNames, doExcSamps, whichTTbarGen, doPURW, doSyst, versNumber, UseUnblindedData);
    
    vector<TString> * sampleAddNames = new vector<TString>;
    vector<int> * sampleStartPositions = new vector<int>;
    sampleStartPositionsNames(whichNTuple, whichTTbarGen, sampleAddNames, sampleStartPositions, doExcSamps);
    
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
    vector<SampleT> * subSampVec    = SubSampVec();
    cout << "subsamp size " << subSampVec->size() << endl;
    vector<SystT> * systVec         = SystVec(SmearedPlots);
    
    //some relevant things for saving names        
    TString TTBarGenName[3] = {"_madgraph", "_mcatnlo", "_powheg"};
    TString nameNTuple = (whichNTuple == 0) ? "_Ovi" : "_DESY";
    TString stringDDEstimate = (useDDEstimate) ? "_wDDEst" : "";
    TString stringExcSamp = (doExcSamps && whichNTuple == 0) ? "_ExcSamps" : "";
    TString stringSignal = "";
    if (doSignal) {
        stringSignal = "_wSignal";
        /*
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
        */        
        stringSignal += "_Stop_";
        stringSignal += vecStopMassGrab->at(0);
        stringSignal += "_Chi0_";
        stringSignal += vecChi0MassGrab->at(0);
        if (typeSMS.Contains("T2bw")) {
            stringSignal += "_Chargino_";
            stringSignal += vecCharginoMassGrab->at(0);
        }
    }   
    TString canvSuffixSaveName = TTBarGenName[whichTTbarGen];
    canvSuffixSaveName += nameNTuple;
    canvSuffixSaveName += stringDDEstimate;
    canvSuffixSaveName += stringExcSamp;
    canvSuffixSaveName += stringSignal;
    if (doReReco) canvSuffixSaveName += "_ReReco";
    if (versNumber == 2) canvSuffixSaveName += "_vers2";
    if (SmearedPlots) canvSuffixSaveName += "_SmearedMET";
    TString saveNameAddition("../Plots/");
    
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
//    float TTBarFullCutSFOviReReco[3] = {1.01276, 0.933961, 1.03732};
 float TTBarFullCutSFOviReReco[3] = {1.01276, 0.933961, 1.00069};
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
    }
    else {
        TTBarSF = TTBarFullCutSFDESY[whichTTbarGen];
    }
    if (useDDEstimate) {
        cout << "TTBarSF used " << TTBarSF << endl;
    }
    
    
    //plotting stuff
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
    vector<float> * signalSkimScaleVec;
    if (doSignal) signalSkimScaleVec = SignalSkimEfficiencyCalc(typeSMS, prefixT2tt, vecStopMassGrab, vecChi0MassGrab, vecCharginoMassGrab, intLumi);
    
    TString plotGrabBaseName = SmearedPlots ? "h_SmearPassMT2llCut" : "h_PassMT2llCut";
    TString plotMCGrabName, plotDataGrabName;
    TString plotSystGrabName;
    TString plotVarName, subSampName;;
    TString dataLegendComp = "Data";
    float   fracRatioYAxisRange = 0.21;
    
    
    const int numMT2llCuts = 5;
    int MT2llCuts[numMT2llCuts] = {80, 90, 100, 110, 120};     
    
    /*
    const int numMT2llCuts = 3;
    int MT2llCuts[numMT2llCuts] = {100, 110, 120};
     */
    if (calcTTBarNorm) {
        MT2llCuts[0] = 80;
    }
    
    int grabChan[5] = {7, 24, 38, 39, 53};
    int grabChanV2[24] = {0, 17, 34, 1, 18, 35, 2, 19, 36, 3, 20, 37, 4, 21, 38, 5, 22, 39, 6, 23, 40, 7, 24, 41};
    if (doYieldV1) {
        for (int iMT2Cut = 0; iMT2Cut < numMT2llCuts; ++iMT2Cut) {
            if (calcTTBarNorm && MT2llCuts[iMT2Cut] != 80) continue;
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
            plotMCGrabName = plotGrabBaseName;
            plotMCGrabName += MT2llCuts[iMT2Cut];
            plotSystGrabName = plotMCGrabName;
            plotDataGrabName = plotMCGrabName;
            plotDataGrabName += subSampVec->at(grabChan[whichChan]).histNameSuffix;
            if (SmearedPlots) {
                plotDataGrabName.Replace(2,5,"");
            }            
            if (doNonSig) {
                HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, plotDataGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
                HistogramAdderData(dataHistTH1Vec, h_DataComp, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");     
                

                
                HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, plotMCGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
                HistogramAdderMC(mcHistTH1Vec, mcCompHist1DCentValVec, sampleStartPositions, sampleAddNames, h_MCComp, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                
                if (calcTTBarNorm) {
                    TTBarSF = DataDrivenTTBarScaleFactor(h_DataComp, h_MCComp, mcCompHist1DCentValVec, 1); 
                    cout << "For TTBarGen " << whichTTbarGen << " and nTuple " << whichNTuple << ", calculated TTBar SF is " << TTBarSF << endl;

                }                
                if (doSyst) {
                    HistogramVecGrabberSystGrab(inputFiles, mcHistTH1Vec, mcHistSystTH1Vec, nVtxBackScaleVec, plotSystGrabName, subSampName, systVec, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);            
                    HistogramProjectorSyst(mcHistSystTH1Vec, mcCompHist1DSystVec, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                    SystGraphMakers(h_MCComp, mcCompHist1DSystVec, errCompSpecSource, errCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, kGray + 1, plotMCGrabName, true, fracRatioYAxisRange, doSymErr, false, SmearedPlots);
                }
                cout << "Cutting on MT2ll equals: " << MT2llCuts[iMT2Cut] << endl;
                cout << "h_DataComp Name " << h_DataComp->GetName() << endl;
                cout << "h_DataComp Bin Content 1 " << h_DataComp->GetBinContent(1) << endl;
                cout << "h_DataComp Bin Error 1 " << h_DataComp->GetBinError(1) << endl;
                cout << "h_DataComp Bin Content 2 " << h_DataComp->GetBinContent(2) << endl;
                cout << "h_DataComp Bin Error 2 " << h_DataComp->GetBinError(2) << endl;
                cout << "h_MCComp Name " << h_MCComp->GetName() << endl;
                cout << "h_MCComp Bin Content 1 " << h_MCComp->GetBinContent(1) << endl;
                cout << "h_MCComp Bin Error 1 " << h_MCComp->GetBinError(1) << endl;
                cout << "h_MCComp Bin Content 2 " << h_MCComp->GetBinContent(2) << endl;
                cout << "h_MCComp Bin Error 2 " << h_MCComp->GetBinError(2) << endl;
                if (dispIndSource) {
                    cout << "Printing Individual Sources yield information: " << endl;
                    for (unsigned int iSource = 0; iSource < mcCompHist1DCentValVec->size(); ++iSource) {
                        cout << "Name of Individual Source: " << mcCompHist1DCentValVec->at(iSource)->GetName() << endl;
                        cout << "Individual Source Bin Content 1 " << mcCompHist1DCentValVec->at(iSource)->GetBinContent(1) << endl;
                        cout << "Individual Source Bin Error 1 " << mcCompHist1DCentValVec->at(iSource)->GetBinError(1) << endl;
                        cout << "Individual Source Bin Content 2 " << mcCompHist1DCentValVec->at(iSource)->GetBinContent(2) << endl;
                        cout << "Individual Source Bin Error 2 " << mcCompHist1DCentValVec->at(iSource)->GetBinError(2) << endl;
                    }
                    if (doSyst) {
                        int grabIndex = -1;                                                
                        vector<float> UpErrBin2, DownErrBin2;
                        
                        for (unsigned int iFile = 0; iFile < inputFiles->size(); ++iFile) {
                            TString fileName = inputFiles->at(iFile)->GetName();
                            cout << "fileName " << fileName << endl;
                            if (fileName.Contains("Data")) continue;                            
                            grabIndex++;
                            
                            vector<TH1F *> * ind1DSystVec = new vector<TH1F *>;
                            vector<TH1 *> * IndSystHistVec = IndividualHistSysts(inputFiles->at(iFile), mcHistTH1Vec->at(grabIndex), plotSystGrabName, subSampName,systVec, whichNTuple);
                            vector<TGraphAsymmErrors *> * errCompSpecSourceSpecFile = new vector<TGraphAsymmErrors *>;
                            vector<TGraphAsymmErrors *> * errCompSpecSourceSpecFile_pStat = new vector<TGraphAsymmErrors *>;
                            HistogramProjectorSyst(IndSystHistVec, ind1DSystVec, 1, 1, 1, -1, -1, -1, -1, "", false, false, fileName);
                            SystGraphMakersIndivSamp((TH1F*) mcHistTH1Vec->at(grabIndex), ind1DSystVec, errCompSpecSourceSpecFile, errCompSpecSourceSpecFile_pStat, doSymErr, SmearedPlots);
                            for (unsigned int iSyst = 0; iSyst < systCanvNameVec->size(); ++iSyst) {
                                cout << "For file: " << fileName << " and syst " << systCanvNameVec->at(iSyst) << endl;
                                cout << "ErrGraph Up Err at point 1 " << errCompSpecSourceSpecFile->at(iSyst)->GetErrorYhigh(1) << endl;
                                cout << "ErrGraph Down Err at point 1 " << errCompSpecSourceSpecFile->at(iSyst)->GetErrorYlow(1) << endl;
                                cout << "ErrGraph Up Err at point 2 " << errCompSpecSourceSpecFile->at(iSyst)->GetErrorYhigh(2) << endl;
                                cout << "ErrGraph Down Err at point 2 " << errCompSpecSourceSpecFile->at(iSyst)->GetErrorYlow(2) << endl;
                                cout << endl;
                                if (systCanvNameVec->at(iSyst).Contains("FullSyst")) {
                                    UpErrBin2.push_back(errCompSpecSourceSpecFile->at(iSyst)->GetErrorYhigh(2));
                                    DownErrBin2.push_back(errCompSpecSourceSpecFile->at(iSyst)->GetErrorYlow(2));
                                }
                            }
                            delete ind1DSystVec;
                            delete IndSystHistVec;
                            delete errCompSpecSourceSpecFile;
                            delete errCompSpecSourceSpecFile_pStat;
                        }
                        PrintSystInfo(&UpErrBin2, &DownErrBin2, mcCompHist1DCentValVec, sampleStartPositions);
                    }
                }
                if (doSyst) {
                    cout << "Printing Comprehensive Systematics Info: " << endl;
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
                if (doSyst) {                    
                    HistogramProjectorSyst(vecCurrSignalSystTH1Hists, vecCurrStop1DSystHists, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                    SystGraphMakers(currSignal1DCentValHist, vecCurrStop1DSystHists, vecCurrSigCompSpecSource, vecCurrSigCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, mcColorsSignal->at(0), plotMCGrabName, true, fracRatioYAxisRange, doSymErr, (doSignal && doNonSig), SmearedPlots);
                }
                cout << "Cutting on MT2ll equals: " << MT2llCuts[iMT2Cut] << endl;
                cout << "currSignal1DCentValHist Name " << currSignal1DCentValHist->GetName() << endl;
                cout << "currSignal1DCentValHist Bin Content 1 " << currSignal1DCentValHist->GetBinContent(1) << endl;
                cout << "currSignal1DCentValHist Bin Error 1 " << currSignal1DCentValHist->GetBinError(1) << endl;
                cout << "currSignal1DCentValHist Bin Content 2 " << currSignal1DCentValHist->GetBinContent(2) << endl;
                cout << "currSignal1DCentValHist Bin Error 2 " << currSignal1DCentValHist->GetBinError(2) << endl;
                if (doSyst) {
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
    }
    
    if (doYieldV2) {
        for (int iChan = 0; iChan < 24; ++iChan) {
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
            subSampName += subSampVec->at(grabChanV2[iChan]).histNameSuffix;
            cout << "grabChanV2[iChan] " << grabChanV2[iChan]  << endl;
            cout << "subSampName " << subSampName << endl;
            plotMCGrabName = plotGrabBaseName;
            plotMCGrabName += 80;
            plotSystGrabName = plotMCGrabName;
            plotDataGrabName = plotMCGrabName;
            plotDataGrabName += subSampVec->at(grabChanV2[iChan]).histNameSuffix;
            
            if (doNonSig) {
                HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, plotDataGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
                HistogramAdderData(dataHistTH1Vec, h_DataComp, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");     
                                
                HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, plotMCGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
                HistogramAdderMC(mcHistTH1Vec, mcCompHist1DCentValVec, sampleStartPositions, sampleAddNames, h_MCComp, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                
                cout << "h_DataComp Name " << h_DataComp->GetName() << endl;
                
                cout << "h_DataComp Yield MT2llCut " << h_DataComp->GetBinContent(1) << endl;
                cout << "h_DataComp Error MT2llCut " << h_DataComp->GetBinError(1) << endl;
                cout << "h_DataComp Yield Total " << h_DataComp->GetBinContent(1) + h_DataComp->GetBinContent(2) << endl;
                cout << "h_DataComp Error Total " << TMath::Sqrt(h_DataComp->GetBinError(2) * h_DataComp->GetBinError(2) + h_DataComp->GetBinError(1) * h_DataComp->GetBinError(1)) << endl;
                
                cout << "h_MCComp Name " << h_MCComp->GetName() << endl;
                cout << "h_MCComp Yield MT2llCut " << h_MCComp->GetBinContent(1) << endl;
                cout << "h_MCComp Error MT2llCut " << h_MCComp->GetBinError(1) << endl;
                cout << "h_MCComp Yield Total " << h_MCComp->GetBinContent(1) + h_MCComp->GetBinContent(2) << endl;
                cout << "h_MCComp Error Total " << TMath::Sqrt(h_MCComp->GetBinError(2) * h_MCComp->GetBinError(2) + h_MCComp->GetBinError(1) * h_MCComp->GetBinError(1)) << endl;                
                for (unsigned int iIndMC = 0; iIndMC < mcCompHist1DCentValVec->size(); ++iIndMC) {
                    cout << "IndSample name " << mcCompHist1DCentValVec->at(iIndMC)->GetName() << endl;
                    cout << "IndSample Yield MT2llCut " << mcCompHist1DCentValVec->at(iIndMC)->GetBinContent(1) << endl;
                    cout << "IndSample Error MT2llCut " << mcCompHist1DCentValVec->at(iIndMC)->GetBinError(1) << endl;
                    cout << "IndSample Yield Total " << mcCompHist1DCentValVec->at(iIndMC)->GetBinContent(1) + mcCompHist1DCentValVec->at(iIndMC)->GetBinContent(2) << endl;
                    cout << "IndSample Error Total " << TMath::Sqrt(mcCompHist1DCentValVec->at(iIndMC)->GetBinError(2) * mcCompHist1DCentValVec->at(iIndMC)->GetBinError(2) + mcCompHist1DCentValVec->at(iIndMC)->GetBinError(1) * mcCompHist1DCentValVec->at(iIndMC)->GetBinError(1)) << endl; 
                }
                
                if (doSyst) {
                    HistogramVecGrabberSystGrab(inputFiles, mcHistTH1Vec, mcHistSystTH1Vec, nVtxBackScaleVec, plotSystGrabName, subSampName, systVec, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);            
                    HistogramProjectorSyst(mcHistSystTH1Vec, mcCompHist1DSystVec, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                    SystGraphMakers(h_MCComp, mcCompHist1DSystVec, errCompSpecSource, errCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, kGray + 1, plotMCGrabName, true, fracRatioYAxisRange, doSymErr, false, SmearedPlots);                
                    for (unsigned int iSyst = 0; iSyst < systCanvNameVec->size(); ++iSyst) {
                        cout << "For hist " << h_MCComp->GetName() << " and syst " << systCanvNameVec->at(iSyst) << endl;
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
                cout << "Signal Name " << currSignal1DCentValHist->GetName() << endl;
                cout << "Signal Yield MT2llCut " << currSignal1DCentValHist->GetBinContent(1) << endl;
                cout << "Signal Error MT2llCut " << currSignal1DCentValHist->GetBinError(1) << endl;
                cout << "Signal Yield Total " << currSignal1DCentValHist->GetBinContent(1) + currSignal1DCentValHist->GetBinContent(2) << endl;
                cout << "Signal Error Total " << TMath::Sqrt(currSignal1DCentValHist->GetBinError(2) * currSignal1DCentValHist->GetBinError(2) + currSignal1DCentValHist->GetBinError(1) * currSignal1DCentValHist->GetBinError(1)) << endl;
                if (doSyst) {
                    
                    HistogramProjectorSyst(vecCurrSignalSystTH1Hists, vecCurrStop1DSystHists, 1, 1, 1, -1, -1, -1, -1, "", false, false, "");
                    SystGraphMakers(currSignal1DCentValHist, vecCurrStop1DSystHists, vecCurrSigCompSpecSource, vecCurrSigCompSpecSource_pStat, fracRatioSystVec, systCanvNameVec, mcColorsSignal->at(0), plotMCGrabName, true, fracRatioYAxisRange, doSymErr, (doSignal && doNonSig), SmearedPlots);
                    for (unsigned int iSyst = 0; iSyst < systCanvNameVec->size(); ++iSyst) {
                        cout << "For hist " << currSignal1DCentValHist->GetName() << " and syst " << systCanvNameVec->at(iSyst) << endl;
                        cout << "Signal ErrGraph Up Err at point 1 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYhigh(1) << endl;
                        cout << "Signal ErrGraph Down Err at point 1 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYlow(1) << endl;
                        cout << "Signal ErrGraph Up Err at point 2 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYhigh(2) << endl;
                        cout << "Signal ErrGraph Down Err at point 2 " << vecCurrSigCompSpecSource->at(iSyst)->GetErrorYlow(2) << endl;
                    }
                }            
            }   
        }
    }                    
    subSampName = subSampVec->at(grabChan[whichChan]).histNameSuffix;
    TCanvas * c_Var;

    
    TString canvNameBase, canvName;
    const Int_t numCanvs = 14;
    TString canvNameAdd[numCanvs] = {"_CentVal", "_LepESShiftUp", "_LepESShiftDown", "_JetESShiftUp", "_JetESShiftDown", "_UncESShiftUp", "_UncESShiftDown", "_LepEffSFShiftUp", "_LepEffSFShiftDown", "_genTopRW", "_genStopXSecShiftUp", "_genStopXSecShiftDown", "_BTagSFShiftUp", "_BTagSFShiftDown"};    
    vector<TH2F *> vecMCComp(numCanvs);
    vector<TH2F *> vecSig(numCanvs);
    vector<TH2F *> vecFOMHists(numCanvs);
    
    vector<TH1F *> vecMCComp1D(numCanvs);
    vector<TH1F *> vecSig1D(numCanvs);
    vector<TH1F *> vecFOMHists1D(numCanvs);
    TString SigStringTitle;
    if (doFOM) {
        SigStringTitle = mcLegendsSignal->at(0);
        if (doOneDeeFOM) {
            plotGrabBaseName = SmearedPlots ? "h_SmearMT2ll" : "h_MT2ll";
            plotMCGrabName = plotGrabBaseName;
            plotSystGrabName = plotMCGrabName;
            plotDataGrabName = plotMCGrabName;
            plotDataGrabName += subSampName;
            if (SmearedPlots) {
                plotDataGrabName.Replace(2,5,"");
            }
            canvNameBase = "c_FOM_OneDee";
            mcHistTH1Vec = new vector<TH1 *>;
            mcHistSystTH1Vec = new vector<TH1 *>;
            vecCurrSignalSystTH1Hists = new vector<TH1 *>;        
            HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, plotMCGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
            vecMCComp1D[0] = HistogramAdderMCOneDee(mcHistTH1Vec);
            HistogramVecGrabber_Signal(inputFilesSignal, signalSkimScaleVec, 0, currSignalCentValTH1Hist, plotSystGrabName, vecCurrSignalSystTH1Hists, systVec, subSampName, false, false, doSyst, allMT2llSystematic);
            vecSig1D[0] = (TH1F*) currSignalCentValTH1Hist;
            //        cout << "vecMCComp[0] ";
            //        cout << vecMCComp[0]->GetName() << " integral " << vecMCComp[0]->Integral() << endl;
            //        cout << "vecSig[0] ";
            //        cout << vecSig[0]->GetName() << " integral " << vecSig[0]->Integral() << endl;
            vecFOMHists1D[0] = FOMHist(vecMCComp1D[0], vecSig1D[0], SigStringTitle, typeFOM, 0);
            if (makePlots) {
                canvName = canvNameBase;
                canvName += canvNameAdd[0];
                canvName += canvSuffixSaveName;
                c_Var = new TCanvas(canvName, canvName, 0, 0, 700, 700);
                //            vecFOMHists1D[0]->Draw("e1");
                vecFOMHists1D[0]->Draw("hist");
                c_Var->SaveAs(saveNameAddition + canvName + TString(".pdf"));
                if (saveDotCFile) c_Var->SaveAs(saveNameAddition + canvName + TString(".C"));
            }
            if (printFOMInfo) {
                
                int XBinMaxVal, YBinMaxVal, ZBinMaxVal;
                TAxis * xAxis = vecFOMHists1D[0]->GetXaxis();
                int MinXBin = xAxis->FindBin(MT2llaxisCut);
                int MaxXBin = xAxis->GetNbins();
                xAxis->SetRange(MinXBin, MaxXBin);
                vecFOMHists1D[0]->GetMaximumBin(XBinMaxVal, YBinMaxVal, ZBinMaxVal);
                cout << "Printing FOM Info " << endl;
                if (SmearedPlots) {
                    cout << "In the range of SMEARED MT2ll = " << xAxis->GetBinLowEdge(MinXBin) << " to " << xAxis->GetBinUpEdge(MaxXBin) << endl;
                    cout << "Max F.O.M., using F.O.M. of " << FOMString << " is " << vecFOMHists1D[0]->GetBinContent(XBinMaxVal, YBinMaxVal) << " and is found cutting at MT2ll = " << xAxis->GetBinLowEdge(XBinMaxVal) << endl;
                }
                else {
                    cout << "In the range of MT2ll = " << xAxis->GetBinLowEdge(MinXBin) << " to " << xAxis->GetBinUpEdge(MaxXBin) << endl;
                    cout << "Max F.O.M., using F.O.M. of " << FOMString << " is " << vecFOMHists1D[0]->GetBinContent(XBinMaxVal, YBinMaxVal) << " and is found cutting at MT2ll = " << xAxis->GetBinLowEdge(XBinMaxVal) << endl;
                }
                /*
                int MinXBin = vecFOMHists1D[0]->FindBin(MT2llaxisCut);
                int MaxXBin = vecFOMHists1D[0]->GetNbinsX();
                TAxis * xAxis = vecFOMHists1D[0]->GetXaxis();
                xAxis->SetRange(MinXBin, MaxXBin);
                cout << "In the range of MT2ll = " << xAxis->GetBinLowEdge(MinXBin) << " to  " << xAxis->GetBinUpEdge(MaxXBin) << endl;
                cout << "Max F.O.M., using F.O.M. of " << FOMString << " is " << vecFOMHists1D[0]->GetMaximum() << " and is found cutting at MT2ll = " << xAxis->GetBinLowEdge(vecFOMHists1D[0]->GetMaximumBin()) << endl;
                */
            }
            else {
                if (doSyst) {
                    HistogramVecGrabberSystGrab(inputFiles, mcHistTH1Vec, mcHistSystTH1Vec, nVtxBackScaleVec, plotSystGrabName, subSampName, systVec, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);
                    for (int iSyst = 0; iSyst < numCanvs - 1; ++iSyst) {
                        if (canvNameAdd[iSyst + 1].Contains("genStop")) {
                            vecMCComp1D[iSyst + 1] = HistogramAdderMCOneDee(mcHistTH1Vec);
                        }
                        else {
                            vecMCComp1D[iSyst + 1] = SystHistFinderOneDee(mcHistSystTH1Vec, canvNameAdd[iSyst + 1]);
                        }
                        cout << "vecMCComp[iSyst + 1] " << vecMCComp1D[iSyst + 1]->GetName() << " Integral " << vecMCComp1D[iSyst + 1]->Integral() << endl;
                        cout << "vecMCComp[iSyst + 1] " << vecMCComp1D[iSyst + 1]->GetName() << " Mean " << vecMCComp1D[iSyst + 1]->GetMean(1) << endl;
                        vecSig1D[iSyst + 1] = SystHistFinderOneDee(vecCurrSignalSystTH1Hists, canvNameAdd[iSyst + 1]);
                        cout << "vecSig[iSyst + 1] " << vecSig1D[iSyst + 1]->GetName() << " Integral " << vecSig1D[iSyst + 1]->Integral() << endl;
                        cout << "vecSig[iSyst + 1] " << vecSig1D[iSyst + 1]->GetName() << " Mean " << vecSig1D[iSyst + 1]->GetMean(1) << endl;
                        vecFOMHists1D[iSyst + 1] = FOMHist(vecMCComp1D[iSyst + 1], vecSig1D[iSyst + 1], SigStringTitle, typeFOM, iSyst + 1);
                        if (makePlots) {
                            canvName = canvNameBase;
                            canvName += canvNameAdd[iSyst + 1];
                            canvName += canvSuffixSaveName;
                            c_Var = new TCanvas(canvName, canvName, 0, 0, 700, 700);
                            //                    vecFOMHists1D[iSyst + 1]->Draw("e1"); 
                            vecFOMHists1D[iSyst + 1]->Draw("hist");
                            c_Var->SaveAs(saveNameAddition + canvName + TString(".pdf"));
                            if (saveDotCFile) c_Var->SaveAs(saveNameAddition + canvName + TString(".C"));
                        }
                    }
                } 
            }
        }
        if (doTwoDeeFOM) {
            plotGrabBaseName = "h_MT2ll_vs_MT2lb";   
            canvNameBase = "c_FOM_TwoDee";
            plotMCGrabName = plotGrabBaseName;
            plotSystGrabName = plotMCGrabName;
            plotDataGrabName = plotMCGrabName;
            plotDataGrabName += subSampName;
            mcHistTH1Vec = new vector<TH1 *>;
            mcHistSystTH1Vec = new vector<TH1 *>;
            vecCurrSignalSystTH1Hists = new vector<TH1 *>;        
            HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, plotMCGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
            vecMCComp[0] = HistogramAdderMCTwoDee(mcHistTH1Vec);
            HistogramVecGrabber_Signal(inputFilesSignal, signalSkimScaleVec, 0, currSignalCentValTH1Hist, plotSystGrabName, vecCurrSignalSystTH1Hists, systVec, subSampName, false, false, doSyst, allMT2llSystematic);
            vecSig[0] = (TH2F*) currSignalCentValTH1Hist;
            //        cout << "vecMCComp[0] ";
            //        cout << vecMCComp[0]->GetName() << " integral " << vecMCComp[0]->Integral() << endl;
            //        cout << "vecSig[0] ";
            //        cout << vecSig[0]->GetName() << " integral " << vecSig[0]->Integral() << endl;
            vecFOMHists[0] = FOMHist(vecMCComp[0], vecSig[0], SigStringTitle, typeFOM, 0);
            if (makePlots) {
                canvName = canvNameBase;
                canvName += canvNameAdd[0];
                canvName += canvSuffixSaveName;
                c_Var = new TCanvas(canvName, canvName, 0, 0, 700, 700);
                vecFOMHists[0]->Draw("colz");
                c_Var->SaveAs(saveNameAddition + canvName + TString(".pdf"));
                if (saveDotCFile) c_Var->SaveAs(saveNameAddition + canvName + TString(".C"));
            }
            if (printFOMInfo) {
                int XBinMaxVal, YBinMaxVal, ZBinMaxVal;
//                cout << " test " << endl;
                TAxis * xAxis = vecFOMHists[0]->GetXaxis();
//                cout << " test2 " << xAxis->GetName() << endl;
                int MinXBin = xAxis->FindBin(MT2llaxisCut);
                int MaxXBin = xAxis->GetNbins();
                TAxis * yAxis = vecFOMHists[0]->GetYaxis();
                int MinYBin = yAxis->FindBin(MT2lbaxisCut);
                int MaxYBin = yAxis->GetNbins();
                xAxis->SetRange(MinXBin, MaxXBin);
                yAxis->SetRange(MinYBin, MaxYBin);
                vecFOMHists[0]->GetMaximumBin(XBinMaxVal, YBinMaxVal, ZBinMaxVal);
                cout << "Printing FOM Info " << endl;
                if (SmearedPlots) {
                    cout << "Currently Smearing doesn't work right with MT2lblb " << endl;
                }
                else {
                    cout << "In the range of MT2ll = " << xAxis->GetBinLowEdge(MinXBin) << " to " << xAxis->GetBinUpEdge(MaxXBin) << endl;
                    cout << "In the range of MT2lb = " << yAxis->GetBinLowEdge(MinYBin) << " to " << yAxis->GetBinUpEdge(MaxYBin) << endl;
                }
                cout << "Max F.O.M., using F.O.M. of " << FOMString << " is " << vecFOMHists[0]->GetBinContent(XBinMaxVal, YBinMaxVal) << " and is found cutting at MT2ll = " << xAxis->GetBinLowEdge(XBinMaxVal) << " and MT2lblb = " << yAxis->GetBinLowEdge(YBinMaxVal) << endl;
            }
            else {
                if (doSyst) {
                    HistogramVecGrabberSystGrab(inputFiles, mcHistTH1Vec, mcHistSystTH1Vec, nVtxBackScaleVec, plotSystGrabName, subSampName, systVec, useDDEstimate, TTBarSF, scaleLumi, allMT2llSystematic, whichNTuple);
                    for (int iSyst = 0; iSyst < numCanvs -1; ++iSyst) {
                        if (canvNameAdd[iSyst + 1].Contains("genStop")) {
                            vecMCComp[iSyst + 1] = HistogramAdderMCTwoDee(mcHistTH1Vec);
                        }
                        else {
                            vecMCComp[iSyst + 1] = SystHistFinderTwoDee(mcHistSystTH1Vec, canvNameAdd[iSyst + 1]);
                        }
                        cout << "vecMCComp[iSyst + 1] " << vecMCComp[iSyst + 1]->GetName() << " Integral " << vecMCComp[iSyst + 1]->Integral() << endl;
                        cout << "vecMCComp[iSyst + 1] " << vecMCComp[iSyst + 1]->GetName() << " Mean " << vecMCComp[iSyst + 1]->GetMean(1) << endl;
                        vecSig[iSyst + 1] = SystHistFinderTwoDee(vecCurrSignalSystTH1Hists, canvNameAdd[iSyst + 1]);
                        cout << "vecSig[iSyst + 1] " << vecSig[iSyst + 1]->GetName() << " Integral " << vecSig[iSyst + 1]->Integral() << endl;
                        cout << "vecSig[iSyst + 1] " << vecSig[iSyst + 1]->GetName() << " Mean " << vecSig[iSyst + 1]->GetMean(1) << endl;
                        vecFOMHists[iSyst + 1] = FOMHist(vecMCComp[iSyst + 1], vecSig[iSyst + 1], SigStringTitle, typeFOM, iSyst + 1);
                        if (makePlots) {
                            canvName = canvNameBase;
                            canvName += canvNameAdd[iSyst + 1];
                            canvName += canvSuffixSaveName;
                            c_Var = new TCanvas(canvName, canvName, 0, 0, 700, 700);
                            vecFOMHists[iSyst + 1]->Draw("colz");
                            c_Var->SaveAs(saveNameAddition + canvName + TString(".pdf"));
                            if (saveDotCFile) c_Var->SaveAs(saveNameAddition + canvName + TString(".C"));
                        }
                    }
                }
            }
        }
    }
    theApp.Run(retVal);
    //    theApp.Terminate(0);
}
