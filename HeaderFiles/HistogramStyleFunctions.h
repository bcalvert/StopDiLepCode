/*
#ifndef Users_BrianCalvert_Desktop_LocalRootRunArea_HistogramStyleFunctions_h
#define Users_BrianCalvert_Desktop_LocalRootRunArea_HistogramStyleFunctions_h
#endif
*/
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TProfile.h"
using namespace std;
////////// Histogram Attribute Setting Functions
void GraphMainAttSet(TGraph * inputGraph, Color_t FillColor = kBlack, Style_t FillStyle = 1001, Color_t LineColor = kBlack, Width_t LineWidth = 2, Color_t MarkerColor = kBlack, Style_t MarkerStyle = 20, Size_t MarkerSize = 1) {
    inputGraph->SetFillColor(FillColor);
    inputGraph->SetFillStyle(FillStyle);
    inputGraph->SetLineColor(LineColor);
    inputGraph->SetLineWidth(LineWidth);
    inputGraph->SetMarkerColor(MarkerColor);
    inputGraph->SetMarkerStyle(MarkerStyle);
    inputGraph->SetMarkerSize(MarkerSize);
}
void HistToGraphCopyAttributes(TH1 * inputHist, TGraph * outputGraph) {
    Color_t FillColor = (Color_t) inputHist->GetFillColor();
    Style_t FillStyle = (Style_t) inputHist->GetFillStyle();
    Color_t LineColor = (Color_t) inputHist->GetLineColor();
    Width_t LineWidth = (Width_t) inputHist->GetLineWidth();
    Color_t MarkerColor = (Color_t) inputHist->GetMarkerColor();
    Style_t MarkerStyle = (Style_t) inputHist->GetMarkerStyle();
    Size_t MarkerSize = (Size_t) inputHist->GetMarkerSize();
    outputGraph->SetFillColor(FillColor);
    outputGraph->SetFillStyle(FillStyle);
    outputGraph->SetLineColor(LineColor);
    outputGraph->SetLineWidth(LineWidth);
    outputGraph->SetMarkerColor(MarkerColor);
    outputGraph->SetMarkerStyle(MarkerStyle);
    outputGraph->SetMarkerSize(MarkerSize);
    outputGraph->SetLineStyle((Style_t) inputHist->GetLineStyle());
}
void HistMainAttSet(TH1 * inputHist, Color_t FillColor = kBlack, Style_t FillStyle = 1001, Color_t LineColor = kBlack, Width_t LineWidth = 2, Color_t MarkerColor = kBlack, Style_t MarkerStyle = 20, Size_t MarkerSize = 1, Style_t LineStyle = 1) {
    inputHist->SetFillColor(FillColor);
    inputHist->SetFillStyle(FillStyle);
    inputHist->SetLineColor(LineColor);
    inputHist->SetLineWidth(LineWidth);
    inputHist->SetMarkerColor(MarkerColor);
    inputHist->SetMarkerStyle(MarkerStyle);
    inputHist->SetMarkerSize(MarkerSize);
    inputHist->SetLineStyle(LineStyle);
}
void HistMainAttSet(TProfile * inputHist, Color_t FillColor = kBlack, Style_t FillStyle = 1001, Color_t LineColor = kBlack, Width_t LineWidth = 2, Color_t MarkerColor = kBlack, Style_t MarkerStyle = 20, Size_t MarkerSize = 1) {
    inputHist->SetFillColor(FillColor);
    inputHist->SetFillStyle(FillStyle);
    inputHist->SetLineColor(LineColor);
    inputHist->SetLineWidth(LineWidth);
    inputHist->SetMarkerColor(MarkerColor);
    inputHist->SetMarkerStyle(MarkerStyle);
    inputHist->SetMarkerSize(MarkerSize);
}
void HistAxisAttSet(TAxis * inputAxis, TString AxisTitle = "Axis", Size_t AxisTitleSize = 0, double AxisTitleOffset = 0., Size_t AxisLabelSize = 0, double AxisLabelOffset = 0., double AxisLB = 0.0, double AxisUB = 0.0) {
    if (AxisTitle != "Axis") inputAxis->SetTitle(AxisTitle);
    inputAxis->SetTitleSize(AxisTitleSize);
    inputAxis->SetTitleOffset(AxisTitleOffset);
    inputAxis->SetLabelSize(AxisLabelSize);
    inputAxis->SetLabelOffset(AxisLabelOffset);
    if (AxisLB != 0.0 && AxisUB != 0.0) inputAxis->SetRangeUser(AxisLB, AxisUB);
}
TH1D * HistZeroOutRange(TH1D * inputHist, int XLB, int XUB, TString name) {
    if (name == inputHist->GetName()) name += "_Clone";
    TH1D * OutHist = (TH1D*) inputHist->Clone(name);
    int NBinsX = inputHist->GetNbinsX();
    std::cout << "NBinsX " << NBinsX << endl;
    std::cout << "XLB " << XLB << endl;
    std::cout << "XUB " << XUB << endl;
    for (int i = 0; i < XLB; ++i) {
        //        std::cout << "i " << i << endl;
        OutHist->SetBinContent(i, 0);
        OutHist->SetBinError(i, 0);
    }
    for (int i = XUB; i < NBinsX; ++i) {
        //        std::cout << "i " << i << endl;
        OutHist->SetBinContent(i, 0);
        OutHist->SetBinError(i, 0);
    }
    return OutHist;
}
