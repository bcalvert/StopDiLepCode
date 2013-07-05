/*
#ifndef Users_BrianCalvert_Desktop_LocalRootRunArea_HistogramStyleFunctions_h
#define Users_BrianCalvert_Desktop_LocalRootRunArea_HistogramStyleFunctions_h
#endif
*/
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TProfile.h"
using namespace std;
////////// Histogram Attribute Setting Functions
void GraphMainAttSet(TGraph * inputGraph, Color_t FillColor = kBlack, Style_t FillStyle = 1001, Color_t LineColor = kBlack, Style_t LineStyle = 1,Width_t LineWidth = 2, Color_t MarkerColor = kBlack, Style_t MarkerStyle = 20, Size_t MarkerSize = 1) {
    inputGraph->SetFillColor(FillColor);
    inputGraph->SetFillStyle(FillStyle);
    inputGraph->SetLineColor(LineColor);
    inputGraph->SetLineStyle(LineStyle);
    inputGraph->SetLineWidth(LineWidth);
    inputGraph->SetMarkerColor(MarkerColor);
    inputGraph->SetMarkerStyle(MarkerStyle);
    inputGraph->SetMarkerSize(MarkerSize);
}

void GraphMainAttSet(TGraphAsymmErrors * inputGraph, Color_t FillColor = kBlack, Style_t FillStyle = 1001, Color_t LineColor = kBlack, Style_t LineStyle = 1,Width_t LineWidth = 2, Color_t MarkerColor = kBlack, Style_t MarkerStyle = 20, Size_t MarkerSize = 1) {
    inputGraph->SetFillColor(FillColor);
    inputGraph->SetFillStyle(FillStyle);
    inputGraph->SetLineColor(LineColor);
    inputGraph->SetLineStyle(LineStyle);
    inputGraph->SetLineWidth(LineWidth);
    inputGraph->SetMarkerColor(MarkerColor);
    inputGraph->SetMarkerStyle(MarkerStyle);
    inputGraph->SetMarkerSize(MarkerSize);
}

void GraphAxisAttSet(TAxis * inputAxis, TString AxisTitle = "Axis", Size_t AxisTitleSize = 0, double AxisTitleOffset = 0., Size_t AxisLabelSize = 0, double AxisLabelOffset = 0., double AxisLB = 0.0, double AxisUB = 0.0) {
    if (AxisTitle != "Axis") inputAxis->SetTitle(AxisTitle);
    inputAxis->SetTitleSize(AxisTitleSize);
    inputAxis->SetTitleOffset(AxisTitleOffset);
    inputAxis->SetLabelSize(AxisLabelSize);
    inputAxis->SetLabelOffset(AxisLabelOffset);
    if (AxisLB != 0.0 && AxisUB != 0.0) inputAxis->SetRangeUser(AxisLB, AxisUB);
}

