#include "TROOT.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TString.h"
#include "TProfile.h"
#include <cmath>
#include <sstream>

#include "../HeaderFiles/GraphStyleFunctions.h"

/////Essentially Just LimitPlotGraphMaker_TGraphAsymmErr but changed to include multiple channels (actually implementing them that is) 1/23/13

using namespace std;

void LimitPlotGraphMakerMultiChannel(int whichStat = 0, int StartIndex = 0, bool doT2TT = false, int whichChan = 0, bool doSyst = false) {
    const int NumMassPoints = 381;
    Float_t ObsLim[NumMassPoints], MedExpLim[NumMassPoints], OneSigUpExpLim[NumMassPoints], OneSigDownExpLim[NumMassPoints], TwoSigUpExpLim[NumMassPoints], TwoSigDownExpLim[NumMassPoints];
    
    Float_t StopMass[NumMassPoints];
    Float_t qE;
    Double_t limit;
    int MassIndex;
    double mStop;
    int e;
    TString StatMethods[2] = {"Asymptotic", "HybridNew"};
    /*
     TGraph * ObsLimGraph, * MedExpLimGraph, * OneSigUpExpLimGraph, * OneSigDownExpLimGraph, * TwoSigUpExpLimGraph, * TwoSigDownExpLimGraph;
     TGraph * ObsLim_eeChanGraph, * MedExpLim_eeChanGraph, * OneSigUpExpLim_eeChanGraph, * OneSigDownExpLim_eeChanGraph, * TwoSigUpExpLim_eeChanGraph, * TwoSigDownExpLim_eeChanGraph;
     TGraph * ObsLim_emuChanGraph, * MedExpLim_emuChanGraph, * OneSigUpExpLim_emuChanGraph, * OneSigDownExpLim_emuChanGraph, * TwoSigUpExpLim_emuChanGraph, * TwoSigDownExpLim_emuChanGraph;
     TGraph * ObsLim_mumuChanGraph, * MedExpLim_mumuChanGraph, * OneSigUpExpLim_mumuChanGraph, * OneSigDownExpLim_mumuChanGraph, * TwoSigUpExpLim_mumuChanGraph, * TwoSigDownExpLim_mumuChanGraph;
     */
    TGraph * ObsLimGraph, * MedExpLimGraph, * OneSigUpExpLimGraph, * OneSigDownExpLimGraph, * TwoSigUpExpLimGraph, * TwoSigDownExpLimGraph;
    
    //Trying Asymmetric Errors    
    TGraphAsymmErrors * MedPlusOneSig, * MedPlusTwoSig;
    Double_t MedPlusOneSigErrUp[NumMassPoints],  MedPlusOneSigErrDown[NumMassPoints],  MedPlusTwoSigErrUp[NumMassPoints],  MedPlusTwoSigErrDown[NumMassPoints];
    
    //Trying Asymmetric Errors
    TString T2TTName = "";
    TString doSystName = "";
    TString FileInName = "StopLimitComposite";
    TString SpecChannelName;
    switch (whichChan) {
        case 0:
            SpecChannelName = "Multi"; 
            break;
        case 1:
            SpecChannelName = "_ee"; 
            break;
        case 2:
            SpecChannelName = "_emu"; 
            break;
        case 3:
            SpecChannelName = "_mumu"; 
            break;
        default:
            break;
    }
    FileInName += SpecChannelName;
    if (doT2TT) {
        T2TTName = "T2TT";   
    }
    if (doSyst) {
      doSystName = "WSyst";
    }
    else {
      doSystName = "NoSyst";
    }
    FileInName += T2TTName;
    FileInName += "Asymptotic";
    FileInName += doSystName;
    FileInName += ".root";
    TFile fIn(FileInName);    
    cout << "File In Name " << FileInName << endl;
    TTree * fileTree = (TTree *) fIn.Get("limit");
    fileTree->SetBranchAddress("mh",&mStop);    
    fileTree->SetBranchAddress("limit",&limit);
    fileTree->SetBranchAddress("quantileExpected",&qE);
    
    Float_t LimSigThresholds[6] = {-1, .025, .16, .5, .84, .975};
    for (int i = 0; true; ++i) {
        e = fileTree->GetEntry(i);
        if (e < 1) {
            cout << "i " << i << endl;
            break;   
        }
        
        MassIndex = mStop/5 - 20;
        StopMass[MassIndex] = mStop;
        if (abs(qE - LimSigThresholds[0]) < 1E-5) {
            //            cout << "case -1" << endl;
            ObsLim[MassIndex] = limit; 
        }
        else if (abs(qE - LimSigThresholds[1]) < 1E-5) {
            //            cout << "case .025" << endl;
            TwoSigDownExpLim[MassIndex] = limit;
        }
        else if (abs(qE - LimSigThresholds[2]) < 1E-5) {
            //            cout << "case 0.16" << endl;
            OneSigDownExpLim[MassIndex] = limit;
        }
        else if (abs(qE - LimSigThresholds[3]) < 1E-5) {
            //            cout << "case 0.5" << endl;
            MedExpLim[MassIndex] = limit;
            if (abs(limit - 1) < 5E-2) cout << "limit for mass = " << mStop << " is " << limit <<  endl;
        }
        else if (abs(qE - LimSigThresholds[4]) < 1E-5) {
            //            cout << "case 0.84" << endl;
            OneSigUpExpLim[MassIndex] = limit;
        }
        else if (abs(qE - LimSigThresholds[5]) < 1E-5) {
            //            cout << "case 0.975" << endl;
            TwoSigUpExpLim[MassIndex] = limit;
        } 
    }
    cout << "MedLim 0 " << MedExpLim[0] << endl;
    cout << "MedLim 0 " << MedExpLim[1] << endl;
    cout << "MedLim 200 " << MedExpLim[200] << endl;
    cout << "MedLim 201 " << MedExpLim[201] << endl;
    
    cout << "OneSigUpLim 0 " << OneSigUpExpLim[0] << endl;
    cout << "OneSigUpLim 0 " << OneSigUpExpLim[1] << endl;
    cout << "OneSigUpLim 200 " << OneSigUpExpLim[200] << endl;
    cout << "OneSigUpLim 201 " << OneSigUpExpLim[201] << endl;
    
    cout << "TwoSigUpLim 0 " << TwoSigUpExpLim[0] << endl;
    cout << "TwoSigUpLim 0 " << TwoSigUpExpLim[1] << endl;
    cout << "TwoSigUpLim 200 " << TwoSigUpExpLim[200] << endl;
    cout << "TwoSigUpLim 201 " << TwoSigUpExpLim[201] << endl;
    //Arrays are filled so now make graphs
    TwoSigDownExpLimGraph = new TGraph(NumMassPoints);
    OneSigDownExpLimGraph = new TGraph(NumMassPoints);
    MedExpLimGraph = new TGraph(NumMassPoints);
    OneSigUpExpLimGraph = new TGraph(NumMassPoints);
    TwoSigUpExpLimGraph = new TGraph(NumMassPoints);
    MedPlusOneSig = new TGraphAsymmErrors(NumMassPoints);
    MedPlusTwoSig = new TGraphAsymmErrors(NumMassPoints);
    for (int i = StartIndex; i < NumMassPoints; ++i) {
        MedPlusOneSigErrUp[i] = abs(OneSigUpExpLim[i] - MedExpLim[i]);
        MedPlusOneSigErrDown[i] = abs(OneSigDownExpLim[i] - MedExpLim[i]);
        MedPlusTwoSigErrUp[i] = abs(TwoSigUpExpLim[i] - MedExpLim[i]);
        MedPlusTwoSigErrDown[i] = abs(TwoSigDownExpLim[i] - MedExpLim[i]);
        TwoSigUpExpLimGraph->SetPoint(i, StopMass[i], TwoSigUpExpLim[i]);
        OneSigUpExpLimGraph->SetPoint(i, StopMass[i], OneSigUpExpLim[i]);
        MedExpLimGraph->SetPoint(i, StopMass[i], MedExpLim[i]);
        MedPlusOneSig->SetPoint(i, StopMass[i], MedExpLim[i]);
        MedPlusOneSig->SetPointError(i, 0, 0, MedPlusOneSigErrDown[i], MedPlusOneSigErrUp[i]);
        MedPlusTwoSig->SetPoint(i, StopMass[i], MedExpLim[i]);
        MedPlusTwoSig->SetPointError(i, 0, 0, MedPlusTwoSigErrDown[i], MedPlusTwoSigErrUp[i]);
        OneSigDownExpLimGraph->SetPoint(i, StopMass[i], OneSigDownExpLim[i]);
        TwoSigDownExpLimGraph->SetPoint(i, StopMass[i], TwoSigDownExpLim[i]);
    }
    float GraphLineWidth = 4;
    Style_t FillStyleBase = 1001;
    Color_t LineColors[5] = {kYellow, kGreen, kBlack, kGreen, kYellow};
    Style_t LineStyles[5] = {1, 1, 1, 1, 1};
    Width_t LineWidths[5] = {1, 1, 2, 1, 1};
    Color_t FillColors[5] = {kYellow, kGreen, kGreen, kYellow, kWhite};
    Style_t FillStyles[5] = {FillStyleBase, FillStyleBase, FillStyleBase, FillStyleBase, FillStyleBase};
    Color_t MarkerColors[5] = {kBlue, kRed, kBlack, kRed, kBlue};
    Style_t MarkerStyles[5] = {20, 20, 20, 20, 20};
    Size_t MarkerSizes[5] = {.25, .25, .25, .25, .25};
    GraphMainAttSet(TwoSigUpExpLimGraph, FillColors[0], FillStyles[0], LineColors[0], LineStyles[0], LineWidths[0], MarkerColors[0], MarkerStyles[0], MarkerSizes[0]);
    GraphMainAttSet(OneSigUpExpLimGraph, FillColors[1], FillStyles[1], LineColors[1], LineStyles[1], LineWidths[1], MarkerColors[1], MarkerStyles[1], MarkerSizes[1]);
    GraphMainAttSet(MedExpLimGraph, FillColors[2], FillStyles[2], LineColors[2], LineStyles[2], LineWidths[2], MarkerColors[2], MarkerStyles[2], MarkerSizes[2]);
    GraphMainAttSet(OneSigDownExpLimGraph, FillColors[3], FillStyles[3], LineColors[3], LineStyles[3], LineWidths[3], MarkerColors[3], MarkerStyles[3], MarkerSizes[3]);
    GraphMainAttSet(TwoSigDownExpLimGraph, FillColors[4], FillStyles[4], LineColors[4], LineStyles[4], LineWidths[4], MarkerColors[4], MarkerStyles[4], MarkerSizes[4]);
    
    GraphMainAttSet(MedPlusOneSig, kGreen, FillStyleBase, kBlack, 1, 1, kBlack, 20, 0);
    GraphMainAttSet(MedPlusTwoSig, kYellow, FillStyleBase, kBlack, 1, 1, kBlack, 20, 0);
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(MedPlusTwoSig);
    mg->Add(MedPlusOneSig);
    //Draw the limit vs. Mass
    TString CanvasTitle = TString("LimitVStopMass_");
    CanvasTitle += SpecChannelName;
    CanvasTitle += T2TTName;
    TString SaveTitle = CanvasTitle;
    SaveTitle += ".pdf";
    CanvasTitle += StatMethods[whichStat];
    TCanvas * LimitVMass = new TCanvas(CanvasTitle, CanvasTitle, 0, 0, 700, 700);
    LimitVMass->SetGrid();
    //    LimitVMass->DrawFrame(5 * StartIndex + 100, 0.00001, 2000, 10000);
    //    TPad * Pad1 = LimitVMass->cd(1);
    LimitVMass->SetLogy();
    //    TwoSigUpExpLimGraph->SetFillColor(kYellow);
    mg->Draw("a3");
    //    TwoSigUpExpLimGraph->Draw("P");
    //    OneSigUpExpLimGraph->Draw("P");
    MedExpLimGraph->Draw("L");
    //    OneSigDownExpLimGraph->Draw("P");
    //    TwoSigDownExpLimGraph->Draw("P");
    //    TwoSigUpExpLimGraph->SetTitle();
    //    TAxis * XAxis = TwoSigUpExpLimGraph->GetXaxis();
    //    TAxis * YAxis = TwoSigUpExpLimGraph->GetYaxis();
    mg->SetTitle();
    TAxis * XAxis = mg->GetXaxis();
    TAxis * YAxis = mg->GetYaxis();    
    float TextXPos = 0.4;
    float TextYPos = 0.2;
    XAxis->SetTitle("Mass [GeV]");
    XAxis->SetRangeUser(100, 2000);
    YAxis->SetTitle("Limit on Signal Strength r");
    
    
    TLegend * leg = new TLegend(TextXPos,TextYPos,TextXPos+0.3,TextYPos + .15);
    leg->AddEntry(MedExpLimGraph, "Median expected limit", "l");
    leg->AddEntry(OneSigUpExpLimGraph, "#pm 1 #sigma expected limit", "f");
    leg->AddEntry(TwoSigUpExpLimGraph, "#pm 2 #sigma expected limit", "f");
    leg->Draw("same");
    LimitVMass->SaveAs(SaveTitle);
    //    MedExpLimGraph->Draw("CP");
    
    
    //    OneSigUpExpLimGraph->Draw("ACP same");
    //    OneSigDownExpLimGraph->Draw("ACP same");
}
