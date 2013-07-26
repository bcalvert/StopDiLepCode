#include <cmath>
#include "TCanvas.h"
#include "TH1F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "/Users/BrianCalvert/Desktop/LocalRootRunArea/HistogramStyleFunctions.h"
float nDigits(float number, int digits) {
  return round(number * std::pow(10.,digits)) / std::pow(10.,digits);
}

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
TH1F * FracRatioHist(TH1F * hist1, TH1F * hist2, TString numString, TString denomString, bool doAbsRatio, TString fracratioBaseName, float fracRatioYAxisRange) {
    TString fracratioName = hist1->GetName();
    fracratioName += fracratioBaseName;
    TString fracRatioTitle = numString;
    fracRatioTitle += "/";
    fracRatioTitle += denomString;
    TH1F * fracRatioHist;
    if (doAbsRatio) {
        fracRatioHist = (TH1F *) hist1->Clone(fracratioName);
        fracRatioHist->Divide(fracRatioHist, hist2, 1, 1, "");
        HistAxisAttSet(fracRatioHist->GetYaxis(), fracRatioTitle, .15, .54, .14, .011, 1.0 - fracRatioYAxisRange, 1.0 + fracRatioYAxisRange); 
        
    }
    else {
        fracRatioHist = (TH1F *) hist2->Clone(fracratioName);
        fracRatioHist->Add(hist1, -1);
        fracRatioHist->Divide(fracRatioHist, hist1, 1, 1, "");        
        HistAxisAttSet(fracRatioHist->GetYaxis(), fracRatioTitle, .15, .54, .14, .011, -1.0 * fracRatioYAxisRange, 1.0 * fracRatioYAxisRange);        
    }
    return fracRatioHist;
}
void FixPad(TPad * &inputPad, int whichPad, TCanvas * &inputCanvas) {
    //    cout << "inputPad " << inputPad << endl;
    float m = 1.3;
    float n_ = 1;
    //    cout << "m " << m << endl;
    Double_t e_ = 15*m; //15
    Double_t k_ = 7*m; //15
    Double_t t_ = 30*m; //30
    Double_t b_ = 50*m; //80
    Double_t g_ = 87*m; //87
    Double_t d_ = 30*m; //30
    Double_t h_ = 350*m; //400
    Double_t w_ = 350*m; //350
    Double_t hl_ = 70*m;
    
    Double_t ww_ = g_ + w_ + e_ ;
    Double_t W_ = g_ + n_*w_ + 2*(n_-1)*e_ + d_;
    Double_t H_ = t_ + h_ + b_ + hl_ +2*k_ ;
    
    Double_t xlow_= 0;
    Double_t ylow_=0.;
    Double_t xup_=0;
    Double_t yup_=1.;
    inputCanvas->SetWindowSize(W_, H_);
    //    cout << "whichPad == " << whichPad << endl;
    if (whichPad == 1) {
        TString padname_("padhigh_");
        padname_ +=inputCanvas->GetTitle();
        xup_ = xlow_ + ww_ /W_;
        yup_ = 1.; 
        ylow_ = (hl_ + b_ + k_ ) /H_;
        inputPad->SetPad(padname_, padname_, xlow_, ylow_, xup_, yup_, kWhite, 0, 0);
        xlow_ += (w_ + 2*e_)/W_;
        inputPad->SetLeftMargin(  g_/ww_ );
        inputPad->SetRightMargin( e_/ww_ );
        inputPad->SetTopMargin(  t_/H_ );
        inputPad->SetBottomMargin( k_/H_ );
        inputPad->SetFillColor(0);
        inputPad->SetTickx(1);
        inputPad->SetTicky(1);
        inputPad->SetFrameFillStyle(0);
        inputPad->SetFrameLineWidth(2);
        inputPad->SetFrameBorderMode(0);
    }
    else {
        //        cout << "trying some ish" << endl;
        TString padname_("padlow_");
        padname_ +=inputCanvas->GetTitle();
        xup_ = xlow_ + ww_ /W_;
        yup_ = (hl_ + b_ + k_ ) /H_; 
        ylow_ = 0;
        inputPad->SetPad(padname_, padname_, xlow_, ylow_, xup_, yup_, kWhite, 0, 0);
        xlow_ += (w_ + 2*e_)/W_;
        inputPad->SetLeftMargin(  g_/ww_ );
        inputPad->SetRightMargin( e_/ww_ );
        inputPad->SetTopMargin( k_/H_ );
        inputPad->SetBottomMargin( b_ /(hl_ + b_ + k_ ) );
        inputPad->SetFillColor(0);
        inputPad->SetTickx(1);
        inputPad->SetTicky(1);
        inputPad->SetFrameFillStyle(0);
        inputPad->SetFrameLineWidth(2);
        inputPad->SetFrameBorderMode(0);
    }
}
void FixPadSingle(TPad * &inputPad, TCanvas * &inputCanvas) {
    //    cout << "inputPad " << inputPad << endl;
    Double_t padLeftMargin = 0.17;
    Double_t padRightMargin = 0.03;
    Double_t padTopMargin = 0.03;
    Double_t padBottomMargin = 0.11;
    Double_t canvWidth = 800.;
    Double_t canvHeight = 800.;
    
    inputCanvas->SetWindowSize(canvWidth, canvHeight);    
    TString padname_ = "Pad_";
    padname_ += inputCanvas->GetTitle();
    Double_t xLow = 0.;
    Double_t xUp = 1.;
    Double_t yLow = 0.;
    Double_t yUp = 1.;
    Color_t borderColor = kWhite;
    Short_t borderSize = 0;
    Short_t borderMode = 0;
    
    inputPad->SetPad(padname_, padname_, xLow, yLow, xUp, yUp, borderColor, borderSize, borderMode);
    inputPad->SetLeftMargin(padLeftMargin);
    inputPad->SetRightMargin(padRightMargin);
    inputPad->SetBottomMargin(padBottomMargin);
    inputPad->SetTopMargin(padTopMargin);
    inputPad->SetFillColor(0);
    inputPad->SetTickx(1);
    inputPad->SetTicky(1);
    inputPad->SetFrameFillStyle(0);
    inputPad->SetFrameLineWidth(2);
    inputPad->SetFrameBorderMode(0);
}