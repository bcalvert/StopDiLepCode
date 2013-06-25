#include <iostream>
#include "../HeaderFiles/BasicFunctions.h"
#include "../HeaderFiles/HistogramStyleFunctions.h"
#include "../HeaderFiles/HistogramDataMCPlottingVer2.h"
#include "../HeaderFiles/HistogramSystematics2.h"
#include "../HeaderFiles/StopPlotSetup.h"
#include "../HeaderFiles/StopFunctionDefinitions_v2.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
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
#include <sstream>

using namespace std;
void relErrMaker( TH1F * &inputErrHist) {
    int NBins = inputErrHist->GetNbinsX();
    float binRelErr;
    for (int i = 1; i < NBins+1; ++i) {
        binRelErr = inputErrHist->GetBinError(i)/inputErrHist->GetBinContent(i);    
        inputErrHist->SetBinError(i, binRelErr);
    }
}
void relErrSet(TH1F * inputRelErrHist, TH1F * inputErrHist) {
    int NBins = inputErrHist->GetNbinsX();
    for (int i = 1; i<NBins+1; ++i) {
        inputErrHist->SetBinError(i, inputErrHist->GetBinContent(i)*inputRelErrHist->GetBinError(i));
    }
}
int main( int argc, char* argv[]) {
    using namespace std;
    int whichChan     = 0;
    int whichNTuple   = 1;
    int whichTTbarGen = 0;
    bool doPURW       = 0;
    bool doSyst       = 1;
    bool addThings    = 1;
    bool reportChi2   = 0;
    bool doCDF        = 0;  
    bool doStats      = 0;
    bool doAbsRatio   = 0;
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
        else if (strncmp (argv[k],"doPURW", 6) == 0) {
            doPURW = 1;
        }
        else if (strncmp (argv[k],"noSyst", 6) == 0) {
            doSyst = 0;
        }
        else if (strncmp (argv[k],"noAdd", 5) == 0) {
            addThings = 0;
        }
        else if (strncmp (argv[k],"Chi2", 4) == 0) {
            reportChi2 = 1;
        }
        else if (strncmp (argv[k],"doAbsRatio", 10) == 0) {
            doAbsRatio = 1;
        }
        else if (strncmp (argv[k],"doStats", 7) == 0) {
            doStats = 1;
        }
    }
    gROOT->ProcessLine("#include <vector>");
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
    
    
    //Set up the file input
    //    vector<TFile *> * inFiles = new vector<TFile*>;
    TString fileInNameBase;
    TString fileNameSuffix;
    TString specNTupString;
    
    if (whichNTuple == 0) {
        fileInNameBase = "TreeAnalysisTop_5311pb-1_";
        specNTupString = "_Oviedo";
        fileNameSuffix =  "_Output.root";
    }
    else {
        fileInNameBase = "";
        specNTupString = "_DESY";
        fileNameSuffix =  "Haddplots.root";
    }

    vector<TString> * fileInNames = StopFileNames(whichNTuple);
    vector<TFile *> * inputFiles  = StopFiles(whichNTuple, fileInNames, whichTTbarGen, doPURW, doSyst);
    vector<TString> * mcLegends   = MCLegends(whichNTuple, addThings);
    vector<Color_t> * mcColors    = MCColors(whichNTuple, addThings);
//    unsigned int numFiles = fileInNames->size();
    
    vector<HistogramT> * histVec_1D = OneDeeHistTVec();
    vector<SampleT> * subSampVec    = SubSampVec();
    cout << "subsamp size " << subSampVec->size() << endl;
    vector<SystT> * systVec         = SystVec();

    //some relevant things for saving names
    TString TTBarGenNameDESY[3] = {"_madgraph", "_mcatnlo", "_powheg"};
    bool doSymErr = 0;
    TString AsymmErrString = "AsymmErr.pdf";
    if (doSymErr) {
        AsymmErrString = ".pdf";
    }
    
    ///Systematics stuff////
    vector<TH1F *> * dataHistVec;
    //    TH1F * h_DataEE, * h_DataMuMu, * h_DataEMu;
    TH1F * h_DataComp, * h_CDF_Data;
    TH1F * h_MCComp, * h_ErrComp;
    TH1F * h_FracratioComp;

    vector<TH1F *> * mcIndHistCentValVec; //will contain central value histos for individual MC samples
    vector<TH1F *> * mcCompHistCentValVec; //will contain the added histos for general categories
    vector<TH1F *> * mcCompHistSystVec; //will contain the added histos across all MC for given systematic
    vector<TGraphAsymmErrors *> * fracRatioSystVec;
    vector<TH1F *> * mcCDFSystVec;
    vector<float> * nVtxBackScaleVec = ScaleBackVecCalc(inputFiles);

    TGraphAsymmErrors * errCompCentVal, * errQuadSum, * errQuadSum_pStat;
    TGraphAsymmErrors * errLepEnSc, * errLepEffSF, * errMT21l;
    TGraphAsymmErrors * errLepEnSc_pStat, * errLepEffSF_pStat, * errMT21l_pStat;
    TString stringLepEffSF = "LepEffSF";
    TString stringLepEnSc = "LepES";
    TString stringMT2ll = "MT2ll";  
    TGraphAsymmErrors * currFracRatioGraph;
    vector<TGraphAsymmErrors *> * errCompSpecSource; // will contain just the systematic error for each specific source -- note, will contain errQuadSum as a final dude
    vector<TGraphAsymmErrors *> * errCompSpecSource_pStat; // for each specific source, will contain Stat + respective systematic error -- note, will contain errQuadSum_pStat as a final dude
    vector<TString> * systCanvNameVec;

    //plotting stuff
    
    
//    double KFactor = 1;
//    bool KFactorScale =      1;
//    bool DataIntegralScale = 0;
    bool doSystCurrPlot;
    float intLumi = 19300;
    TString mcplot, mcplot_preRW, dataplot;
    TString plotVarName, subSampName;;
    TString dataLegendComp = "Data";
    THStack * mcStack;
    TString mcStackName;
    float   fracRatioYAxisRange = (doAbsRatio) ? 0.6 : 1.5;
    unsigned int numPlots = histVec_1D->size(); 
    TString saveNameAddition("../Plots/");

    bool doOverflow[numPlots];
    bool doUnderflow[numPlots];
    int  RBNX[numPlots];
    int  RBNY[numPlots];
    int  RBNZ[numPlots];
    float XaxisLegendPos[numPlots];
    float YaxisLegendStart[numPlots];
    for (int i = 0; i < numPlots; ++i) {
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
//    TCanvas * c_Var[numPlots];
    TCanvas * c_Var;
    
    TString canvName, systCanvName;
    logYPad1 = true;
    for (int k = 0; k < numPlots; ++k) {
//    for (int k = 18; k < 19; ++k) {       
        dataHistVec = new vector<TH1F *>;
        mcIndHistCentValVec = new vector<TH1F *>;
        mcCompHistCentValVec = new vector<TH1F *>;
        mcCompHistSystVec = new vector<TH1F *>;
        fracRatioSystVec = new vector<TGraphAsymmErrors *>;
        mcCDFSystVec = new vector<TH1F *>;
        errCompSpecSource = new vector<TGraphAsymmErrors *>;
        errCompSpecSource_pStat = new vector<TGraphAsymmErrors *>;
        systCanvNameVec = new vector<TString>;
        plotVarName = "";
        plotVarName += histVec_1D->at(k).name;
        plotVarName.Replace(0, 2, "");
        subSampName = "";
        subSampName += subSampVec->at(whichChan).histNameSuffix;
        cout << "subSampName " << subSampName << endl;
        cout << "Doing " << plotVarName << endl;
        dataplot = histVec_1D->at(k).name;
        dataplot += subSampName;
        cout << "dataplot " << dataplot << endl;
        mcplot = histVec_1D->at(k).name;
        canvName = "c_";
        canvName += plotVarName;
        canvName += subSampName;
        cout << "test " << endl;
        if (whichNTuple == 1) canvName += TTBarGenNameDESY[whichTTbarGen];
        mcStackName = "mcStack_";
        mcStackName += plotVarName;
        mcStackName += subSampName;
//        c_Var[k] = new TCanvas(canvName, canvName, wtopx, wtopy, W_, H_);
        c_Var = new TCanvas(canvName, canvName, wtopx, wtopy, W_, H_);
        mcStack = new THStack(mcStackName, "");
        cout << "test 2" << endl;
        doSystCurrPlot = (doSyst && histVec_1D->at(k).doXSyst)  ? true : false;
        HistogramVecGrabber(inputFiles, dataHistVec, mcIndHistCentValVec, mcCompHistSystVec, nVtxBackScaleVec, systVec, dataplot, mcplot, subSampName, RBNX[k], RBNY[k], RBNZ[k], doOverflow[k], doUnderflow[k], doSystCurrPlot);
        cout << "test 2a" << endl;
        HistogramAdderSyst(dataHistVec, mcIndHistCentValVec, mcCompHistCentValVec, h_DataComp, h_MCComp, h_FracratioComp, doAbsRatio, fracRatioYAxisRange);
        cout << "test 3" << endl;
//        cout << "QCD " << h_QCDComp->Integral() << endl;
        h_ErrComp = (TH1F *) h_MCComp->Clone();
        errCompCentVal = clonePoints(h_ErrComp);
        cout << "test 3a" << endl;
        // calculate systematics error tgraphs              
        HistMainAttSet(h_DataComp, kWhite, 0, kBlack, 2, kBlack, 20, 0.9);
//        HistMainAttSet(h_MCComp, kWhite, 0, kRed, 2, kWhite, 0, 0);
        HistMainAttSet(h_MCComp, kWhite, 0, kWhite, 2, kWhite, 0, 0);
        HistMainAttSet(h_ErrComp, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0);
        histMax = TMath::Max(h_DataComp->GetMaximum(), h_MCComp->GetMaximum());
        YAxisUB = (YAxisUBBase > 1) ? histMax * 2.5E2 : histMax * 2.5;
        /*
        if (addThings) {
            TH1F * mcPlotHistsOvi[6] = {h_WJComp, h_VVComp, h_ZDYComp, h_SingTopComp, h_TTbarComp};
        }
        else {
            TH1F * mcPlotHistsDESY[9] = {h_WJComp, OneDeeHistsDESY[5], OneDeeHistsDESY[6], OneDeeHistsDESY[7], h_ZDYComp, h_SingTopComp, OneDeeHistsDESY[3], OneDeeHistsDESY[4], h_QCDComp};
        }
        */
        for (unsigned int j = 0; j < mcCompHistCentValVec->size(); ++j) {
            HistMainAttSet(mcCompHistCentValVec->at(j), mcColors->at(j), 1001, mcColors->at(j), 2, kWhite, 0, 0);
            mcStack->Add(mcCompHistCentValVec->at(j));
            cout << "integral for mcCompCentVal " << mcCompHistCentValVec->at(j)->GetName() << " is " << mcCompHistCentValVec->at(j)->Integral() << endl;
        }   
        cout << "integral for MCComp " << h_MCComp->GetName() << " is " << h_MCComp->Integral() << endl;
        for (unsigned int k = 0; k < mcCompHistSystVec->size(); ++k) {
            cout << "integral for syst " << mcCompHistSystVec->at(k)->GetName() << " is " << mcCompHistSystVec->at(k)->Integral() << endl;
        }
//        cout << "Data Integral " << h_DataComp->Integral() << endl;
//        cout << "MC Integral " << h_MCComp->Integral() << endl;
        if (plotVarName == "MT2ll" ) {
            cout << "DataIntegral for MT2ll > " << h_DataComp->GetBinLowEdge(21) << " is " << h_DataComp->Integral(21, h_DataComp->GetNbinsX()+1) << endl;
            cout << "MCIntegral for MT2ll > " << h_MCComp->GetBinLowEdge(21) << " is " << h_MCComp->Integral(21, h_DataComp->GetNbinsX()+1) << endl;
        }
//        SpectrumDraw(c_Var[k], h_DataComp, dataLegendComp, h_MCComp, h_FracratioComp, h_ErrComp, mcStack, XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, mcLegends, mcCompHistCentValVec, doStats, "", intLumi);
        SpectrumDraw(c_Var, h_DataComp, dataLegendComp, h_MCComp, h_FracratioComp, h_ErrComp, mcStack, XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, logYPad1, mcLegends, mcCompHistCentValVec, doStats, "", intLumi);
//        c_Var[k]->SaveAs(canvName + TString(".pdf"));
        c_Var->SaveAs(saveNameAddition + canvName + TString(".pdf"));
        
        if (doSystCurrPlot) {
            errLepEffSF = GraphSystErrorSet_SingleSource(h_MCComp, mcCompHistSystVec, stringLepEffSF + TString("Shift"), doSymErr, 0);
            errLepEffSF_pStat = GraphSystErrorSumErrors(errCompCentVal, errLepEffSF, h_MCComp);
            GraphMainAttSet(errLepEffSF, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
            GraphMainAttSet(errLepEffSF_pStat, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
            errCompSpecSource->push_back(errLepEffSF);
            errCompSpecSource_pStat->push_back(errLepEffSF_pStat);
            systCanvNameVec->push_back(stringLepEffSF); 
            currFracRatioGraph = fracGraph(h_MCComp, errLepEffSF, doAbsRatio, fracRatioYAxisRange);
            fracRatioSystVec->push_back(currFracRatioGraph);
            cout << "test 3b" << endl;
            errLepEnSc = GraphSystErrorSet_SingleSource(h_MCComp, mcCompHistSystVec, stringLepEnSc + TString("Shift"), doSymErr, 0);
            errLepEnSc_pStat = GraphSystErrorSumErrors(errCompCentVal, errLepEnSc, h_MCComp);
            GraphMainAttSet(errLepEnSc, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
            GraphMainAttSet(errLepEnSc_pStat, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
            errCompSpecSource->push_back(errLepEnSc);
            errCompSpecSource_pStat->push_back(errLepEnSc_pStat);
            systCanvNameVec->push_back(stringLepEnSc);
            currFracRatioGraph = fracGraph(h_MCComp, errLepEnSc, doAbsRatio, fracRatioYAxisRange);
            fracRatioSystVec->push_back(currFracRatioGraph);
            cout << "test 3c" << endl;
            if (plotVarName.Contains("MT2ll")) {
                errMT21l = GraphSystErrorSet_SingleSource(h_MCComp, mcCompHistSystVec, stringMT2ll + TString("Shift"), doSymErr, 1);
                errMT21l_pStat = GraphSystErrorSumErrors(errCompCentVal, errMT21l, h_MCComp);
                GraphMainAttSet(errMT21l, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
                GraphMainAttSet(errMT21l_pStat, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
                errCompSpecSource->push_back(errMT21l);
                errCompSpecSource_pStat->push_back(errMT21l_pStat);
                systCanvNameVec->push_back(stringMT2ll);
                currFracRatioGraph = fracGraph(h_MCComp, errMT21l, doAbsRatio, fracRatioYAxisRange);
                fracRatioSystVec->push_back(currFracRatioGraph);
            }        
            else {
                errMT21l = NULL;
            }
            cout << "test 3d" << endl;
            errQuadSum = GraphSystErrorSumErrors(errCompCentVal, errCompSpecSource, false, h_MCComp);
            errQuadSum_pStat = GraphSystErrorSumErrors(errCompCentVal, errCompSpecSource, true, h_MCComp);
            GraphMainAttSet(errQuadSum, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
            GraphMainAttSet(errQuadSum_pStat, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
            errCompSpecSource->push_back(errQuadSum);
            errCompSpecSource_pStat->push_back(errQuadSum_pStat);
            systCanvNameVec->push_back("FullSyst");
            cout << "test 3e" << endl;
            currFracRatioGraph = fracGraph(h_MCComp, errQuadSum, doAbsRatio, fracRatioYAxisRange);
            fracRatioSystVec->push_back(currFracRatioGraph);
            for (unsigned int iSyst = 0; iSyst < systCanvNameVec->size(); ++iSyst) {
                systCanvName = canvName + systCanvNameVec->at(iSyst);
                c_Var = new TCanvas(systCanvName, systCanvName, wtopx, wtopy, W_, H_);
                SpectrumDrawSyst(c_Var, h_DataComp, dataLegendComp, h_MCComp, mcStack, errCompSpecSource_pStat->at(iSyst), errCompSpecSource->at(iSyst), h_FracratioComp, fracRatioSystVec->at(iSyst), XaxisLegendPos[k], YaxisLegendStart[k], YAxisLB, YAxisUB, logYPad1, mcLegends, mcCompHistCentValVec, doStats, "", intLumi);
                c_Var->SaveAs(saveNameAddition + systCanvName + TString(".pdf"));
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
        delete dataHistVec;
        delete mcIndHistCentValVec;
        delete mcCompHistCentValVec;
        delete mcCompHistSystVec;
        delete fracRatioSystVec;
        delete mcCDFSystVec;
        delete errCompSpecSource;
        delete errCompSpecSource_pStat;
        delete systCanvNameVec;
    }
    theApp.Run(retVal);
    //    theApp.Terminate(0);
}
