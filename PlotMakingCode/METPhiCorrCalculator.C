//#include "/Users/BrianCalvert/Desktop/LocalRootRunArea/HistogramProjectionFunctions.h"
//#include "/Users/BrianCalvert/Desktop/LocalRootRunArea/HistogramResolutionFunctions.h"
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
#include "TFitResult.h"
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
    vector<int> * vecStopMassGrab = new vector<int>;       // vector to hold the list of Stop masses to brab
    vector<int> * vecChi0MassGrab = new vector<int>;       // vector to hold the list of Chi0 masses to brab
    vector<int> * vecCharginoMassGrab = new vector<int>;   // vector to hold the list of Chargino masses to brab
    bool doPURW       = 0;          // grab the nVtx reweighted MC files
    bool doSyst       = 1;          // look at systematics plots
    bool addThings    = 1;          // Add together similar kinds of events (for aesthetic reasons) like VV backgrounds -- (6/25/13) don't turn off for now haven't validated code fully when not adding
    bool doReReco     = 0;
    bool doSymErr     = 0;
    bool saveDotCFile = 0;          // If on, saves a .C version of every plot that is made
    bool useDDEstimate = 0;         // whether or not to use data-driven estimates for appropriate backgrounds -- as of right now it is just the TTBar norm to MT2ll < 80 GeV (7/22/13)
    bool allMT2llSystematic = 0;    // Whether or not to use the MT2ll systematic smearing for all MC or just TTBar
    int  versNumber     = 1;
    bool doVectorIndDists = 0;
    bool makeRootCopy = 0;
    float c_0d, c_0mc, c_1d, c_1mc;
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
        else if (strncmp (argv[k],"versNum", 7) == 0) {
            versNumber = strtol(argv[k+1], NULL, 10 );
        }        
        else if (strncmp (argv[k],"doReReco", 8) == 0) {
            doReReco = 1;
        }
        else if (strncmp (argv[k],"doNonSig", 8) == 0) {
            doNonSig = 1;
        }
        else if (strncmp (argv[k],"sDCF", 4) == 0) {
            saveDotCFile = 1;
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
    vector<TFile *> * inputFiles  = StopFiles(whichNTuple, fileInNames, doExcSamps, whichTTbarGen, doPURW, doSyst, versNumber);
    
    vector<TString> * sampleAddNames = new vector<TString>;
    vector<int> * sampleStartPositions = new vector<int>;
    sampleStartPositionsNames(whichNTuple, whichTTbarGen, sampleAddNames, sampleStartPositions, doExcSamps);
    
    vector<TFile *> * inputFilesSignal;
    vector<TString> * mcLegendsSignal;
    vector<Color_t> * mcColorsSignal;
    vector<Style_t> * mcStylesSignal;
    vector<SampleT> * subSampVec    = SubSampVec();
    cout << "subsamp size " << subSampVec->size() << endl;
    vector<SystT> * systVec         = SystVec();
    
    //some relevant things for saving names        
    TString TTBarGenName[3] = {"_madgraph", "_mcatnlo", "_powheg"};
    TString nameNTuple = (whichNTuple == 0) ? "_Ovi" : "_DESY";
    TString stringDDEstimate = (useDDEstimate) ? "_wDDEst" : "";
    TString stringExcSamp = (doExcSamps && whichNTuple == 0) ? "_ExcSamps" : "";
    TString stringSignal = "";
    TString canvSuffixSaveName = TTBarGenName[whichTTbarGen];
    canvSuffixSaveName += nameNTuple;
    canvSuffixSaveName += stringDDEstimate;
    canvSuffixSaveName += stringExcSamp;
    canvSuffixSaveName += stringSignal;
    if (doReReco) canvSuffixSaveName += "_ReReco";
    if (versNumber == 2) canvSuffixSaveName += "_vers2"; 
    cout << "canvSave name " << canvSuffixSaveName << endl;
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
    
    
    TLegend * leg;
    
    TLatex * tl = new TLatex();
    tl->SetTextAlign(12);
    tl->SetNDC();
    tl->SetTextSize(0.03);
    char buf[99];
    
    TAxis * YAxis;
    TAxis * XAxis;
    Option_t * FitOption = "IS";    
    float FitLowRange = 2;
    float FitHighRange = 25;
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
    float TTBarFullCutSFOviReReco[3] = {1.01276, 0.933961, 1.03732};
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
        cout << "TTBarSF used " << TTBarSF;
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
    TString plotGrabBaseName = "h_PassMT2llCut";
    TString plotGrabName;
    TString plotSystGrabName;
    TString plotVarName, subSampName;;
    TString dataLegendComp = "Data";
    float   fracRatioYAxisRange = 0.21;
    
    int grabChan[7] = {0, 17, 34, 1, 18, 35, 55};
    
    TString canvNameBase = "h_METX_vs_nVtx_noPhiCorr", canvName;
    TString canvNameAdd[7] = {"_mumu", "_ee", "_emu", "_mumuZVeto", "_eeZVeto", "_emuZVeto", "_inclusive"};
    plotGrabBaseName = "h_METX_vs_nVtx_noPhiCorr";
    TCanvas * c_Var;
    TH2F * DataComp, * MCComp;
    TH1F * DataMean, * MCMean;
    TF1 * f1_resp_Lin = new TF1("fb_Lin","[0] + x*[1]",FitLowRange, FitHighRange);
    TFitResultPtr fitres_Lin;
    vector<TH1F *> * DataVec, * MCVec;
    for (int iChan = 0; iChan < 7; ++iChan) {
        subSampName = subSampVec->at(grabChan[iChan]).histNameSuffix;        
        plotGrabName = plotGrabBaseName;
        plotSystGrabName = plotGrabName;                
        plotGrabName += subSampName; 
        dataHistTH1Vec = new vector<TH1 *>;
        mcHistTH1Vec = new vector<TH1 *>;
        mcHistSystTH1Vec = new vector<TH1 *>;
        DataVec = new vector<TH1F *>;
        MCVec = new vector<TH1F *>;
        HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, plotGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);        
        HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, plotGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
        DataComp = HistogramAdderDataTwoDee(dataHistTH1Vec);
        MCComp = HistogramAdderDataTwoDee(mcHistTH1Vec);
        VectorDistMakerMean(DataMean, DataVec, DataComp, plotGrabName+TString("_mean"), plotGrabName + TString("_Vec"));
        VectorDistMakerMean(MCMean, DataVec, MCComp, plotGrabName+TString("_MC_mean"), plotGrabName + TString("_MC_Vec"));
        HistMainAttSet(DataMean, 0, 0, kBlue, DataMean->GetLineWidth(), kBlue, 20, DataMean->GetMarkerSize());
        HistMainAttSet(MCMean, 0, 0, kRed, MCMean->GetLineWidth(), kRed, 24, MCMean->GetMarkerSize());
        canvName = canvNameBase;
        canvName += canvNameAdd[iChan];
        canvName += canvSuffixSaveName;        
        c_Var = new TCanvas(canvName, canvName, 0, 0, 700, 700);
        DataMean->Draw();
        MCMean->Draw("same");    
        XAxis = DataMean->GetXaxis();
        YAxis = DataMean->GetYaxis();
        HistAxisAttSet(YAxis, "<E_{x}^{miss}>  [GeV/c]", YAxis->GetTitleSize(), YAxis->GetTitleOffset(), YAxis->GetLabelSize(), YAxis->GetLabelOffset(), -30.0, 30.0);
        leg= new TLegend(0.2,0.65,0.60,0.85);
        f1_resp_Lin->SetLineColor(kBlue);
        f1_resp_Lin->SetLineWidth(2);
        f1_resp_Lin->SetLineStyle(1);
        fitres_Lin = DataMean->Fit(f1_resp_Lin,FitOption,"axis same",FitLowRange, FitHighRange);
        c_0d = fitres_Lin->Parameter(0);
        c_1d = fitres_Lin->Parameter(1);
        sprintf(buf,"data A+B+C+D c_{0} = %0.2f #pm %0.2f",fitres_Lin->Parameter(0),fitres_Lin->ParError(0));
        leg->AddEntry(DataMean,buf,"p");
        sprintf(buf,"data A+B+C+D c_{1} = %0.2f #pm %0.2f",fitres_Lin->Parameter(1),fitres_Lin->ParError(1));
        leg->AddEntry(DataMean,buf,"p");
        f1_resp_Lin->SetLineStyle(2);
        f1_resp_Lin->SetLineColor(kRed);
        fitres_Lin = MCMean->Fit(f1_resp_Lin,FitOption,"axis same",FitLowRange, FitHighRange);
        c_0mc = fitres_Lin->Parameter(0);
        c_1mc = fitres_Lin->Parameter(1);
        sprintf(buf,"sim A+B+C+D c_{0} = %0.2f #pm %0.2f",fitres_Lin->Parameter(0),fitres_Lin->ParError(0));
        leg->AddEntry(MCMean,buf,"p");
        sprintf(buf,"sim A+B+C+D c_{1} = %0.2f #pm %0.2f",fitres_Lin->Parameter(1),fitres_Lin->ParError(1));
        leg->AddEntry(MCMean,buf,"p");
        leg->Draw("same");
        c_Var->SaveAs(canvName + TString(".pdf"));
        if (makeRootCopy) {
            c_Var->SaveAs(canvName + TString(".root"));
        }
        cout << "{c_0 X data, c_1 X data, c_0 X MC, c_1 X MC} is: {" << c_0d << ", " << c_1d << ", " << c_0mc << ", " << c_1mc << "}" << endl; 
    }
    
    canvNameBase = "h_METY_vs_nVtx_noPhiCorr";
    plotGrabBaseName = "h_METY_vs_nVtx_noPhiCorr";
    
    for (int iChan = 0; iChan < 7; ++iChan) {
        subSampName = subSampVec->at(grabChan[iChan]).histNameSuffix;        
        plotGrabName = plotGrabBaseName;
        plotSystGrabName = plotGrabName;                
        plotGrabName += subSampName; 
        dataHistTH1Vec = new vector<TH1 *>;
        mcHistTH1Vec = new vector<TH1 *>;
        mcHistSystTH1Vec = new vector<TH1 *>;
        DataVec = new vector<TH1F *>;
        MCVec = new vector<TH1F *>;
        HistogramVecGrabberCentValGrab(inputFiles, true, dataHistTH1Vec, nVtxBackScaleVec, plotGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);        
        HistogramVecGrabberCentValGrab(inputFiles, false, mcHistTH1Vec, nVtxBackScaleVec, plotGrabName, subSampName, useDDEstimate, TTBarSF, scaleLumi);
        DataComp = HistogramAdderDataTwoDee(dataHistTH1Vec);
        MCComp = HistogramAdderDataTwoDee(mcHistTH1Vec);
        VectorDistMakerMean(DataMean, DataVec, DataComp, plotGrabName+TString("_mean"), plotGrabName + TString("_Vec"));
        VectorDistMakerMean(MCMean, DataVec, MCComp, plotGrabName+TString("_MC_mean"), plotGrabName + TString("_MC_Vec"));
        HistMainAttSet(DataMean, 0, 0, kBlue, DataMean->GetLineWidth(), kBlue, 20, DataMean->GetMarkerSize());
        HistMainAttSet(MCMean, 0, 0, kRed, MCMean->GetLineWidth(), kRed, 24, MCMean->GetMarkerSize());
        canvName = canvNameBase;
        canvName += canvNameAdd[iChan];
        canvName += canvSuffixSaveName;        
        c_Var = new TCanvas(canvName, canvName, 0, 0, 700, 700);
        DataMean->Draw();
        MCMean->Draw("same");    
        XAxis = DataMean->GetXaxis();
        YAxis = DataMean->GetYaxis();
        HistAxisAttSet(YAxis, "<E_{y}^{miss}>  [GeV/c]", YAxis->GetTitleSize(), YAxis->GetTitleOffset(), YAxis->GetLabelSize(), YAxis->GetLabelOffset(), -30.0, 30.0);
        leg= new TLegend(0.2,0.65,0.60,0.85);
        f1_resp_Lin->SetLineColor(kBlue);
        f1_resp_Lin->SetLineWidth(2);
        f1_resp_Lin->SetLineStyle(1);
        fitres_Lin = DataMean->Fit(f1_resp_Lin,FitOption,"axis same",FitLowRange, FitHighRange);
        c_0d = fitres_Lin->Parameter(0);
        c_1d = fitres_Lin->Parameter(1);
        sprintf(buf,"data A+B+C+D c_{0} = %0.2f #pm %0.2f",fitres_Lin->Parameter(0),fitres_Lin->ParError(0));
        leg->AddEntry(DataMean,buf,"p");
        sprintf(buf,"data A+B+C+D c_{1} = %0.2f #pm %0.2f",fitres_Lin->Parameter(1),fitres_Lin->ParError(1));
        leg->AddEntry(DataMean,buf,"p");
        f1_resp_Lin->SetLineStyle(2);
        f1_resp_Lin->SetLineColor(kRed);
        fitres_Lin = MCMean->Fit(f1_resp_Lin,FitOption,"axis same",FitLowRange, FitHighRange);
        c_0mc = fitres_Lin->Parameter(0);
        c_1mc = fitres_Lin->Parameter(1);
        sprintf(buf,"sim A+B+C+D c_{0} = %0.2f #pm %0.2f",fitres_Lin->Parameter(0),fitres_Lin->ParError(0));
        leg->AddEntry(MCMean,buf,"p");
        sprintf(buf,"sim A+B+C+D c_{1} = %0.2f #pm %0.2f",fitres_Lin->Parameter(1),fitres_Lin->ParError(1));
        leg->AddEntry(MCMean,buf,"p");
        leg->Draw("same");
        c_Var->SaveAs(canvName + TString(".pdf"));
        if (makeRootCopy) {
            c_Var->SaveAs(canvName + TString(".root"));
        }
        cout << "{c_0 Y data, c_1 Y data, c_0 Y MC, c_1 Y MC} is: {" << c_0d << ", " << c_1d << ", " << c_0mc << ", " << c_1mc << "}" << endl; 
    }
    
    
    
    
    
    theApp.Run(retVal);
    //    theApp.Terminate(0);
}
