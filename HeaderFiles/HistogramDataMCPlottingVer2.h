//#include "/Users/BrianCalvert/Desktop/LocalRootRunArea/HistogramStyleFunctions.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
using namespace std;
/*
////Matthieu's Pad Thing////////////////
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
*/
/*
 vector<vector<TPad*> > PreparePadsWithRatio(int npad, TCanvas * InputCanvas) {
 
 float m=1.3;
 
 size_t n_ = npad;
 vector<TPad*> padhigh_;
 vector<TPad*> padlow_;
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
 
 InputCanvas->SetWindowSize(W_, H_);
 //c2=new TCanvas("c2","Test",300,300,W_,H_);
 //  c2=new TCanvas("c2","Test",300,300,1.2*W_,1.2*H_);
 Double_t xlow_= 0;
 Double_t ylow_=0.;
 Double_t xup_=0;
 Double_t yup_=1.;
 
 
 for(int unsigned i=0;i<n_;i++) {
 
 TString padname_("padhigh_");
 padname_ +=InputCanvas->GetTitle();
 padname_+=i;
 xup_ = xlow_ + ww_ /W_;
 yup_ = 1.; 
 ylow_ = (hl_ + b_ + k_ ) /H_;
 
 
 TPad* padtmp = new TPad( padname_, padname_, 
 xlow_, ylow_, xup_, yup_,
 kWhite,0,0);
 xlow_ += (w_ + 2*e_)/W_;
 padtmp->SetLeftMargin(  g_/ww_ );
 padtmp->SetRightMargin( e_/ww_ );
 padtmp->SetTopMargin(  t_/H_ );
 padtmp->SetBottomMargin( k_/H_ );
 padtmp->SetFillColor(0);
 padtmp->SetTickx(1);
 padtmp->SetTicky(1);
 padtmp->SetFrameFillStyle(0);
 padtmp->SetFrameLineWidth(2);
 padtmp->SetFrameBorderMode(0);
 
 padhigh_.push_back(padtmp);
 
 }
 
 xlow_= 0; //reset
 c2->cd(); //reset
 for(int unsigned i=0;i<n_;i++) {
 
 TString padname_("padlow_");
 padname +=InputCanvas->GetTitle();
 padname_+=i;
 xup_ = xlow_ + ww_ /W_;
 yup_ = (hl_ + b_ + k_ ) /H_; 
 ylow_ = 0;
 
 
 TPad* padtmp  = new TPad( padname_, padname_, 
 xlow_, ylow_, xup_, yup_,
 kWhite,0,0);
 xlow_ += (w_ + 2*e_)/W_;
 padtmp->SetLeftMargin(  g_/ww_ );
 padtmp->SetRightMargin( e_/ww_ );
 padtmp->SetTopMargin( k_/H_ );
 padtmp->SetBottomMargin( b_ /(hl_ + b_ + k_ ) );
 padtmp->SetFillColor(0);
 padtmp->SetTickx(1);
 padtmp->SetTicky(1);
 padtmp->SetFrameFillStyle(0);
 padtmp->SetFrameLineWidth(2);
 padtmp->SetFrameBorderMode(0);
 
 padlow_.push_back(padtmp);
 } 
 
 
 vector<vector<TPad*> > pad_;
 pad_.push_back(padhigh_);
 pad_.push_back(padlow_);
 
 return pad_;
 
 }
 */
////Matthieu's Pad Thing//////////
TH1 * HistGrabber(TString name, TFile ** inputTFileArray, int ArrayLength, TString newName) {
    TH1 * OutHist;
    cout << "name: " << name << endl;
    if (ArrayLength != 0) {
        OutHist = (TH1*) inputTFileArray[0]->Get(name)->Clone(newName);   
    }
    else {
        cout << "bad mojo!" << endl;
    }
    for (int i = 1; i < ArrayLength; ++i) {     
        OutHist->Add((TH1*) inputTFileArray[i]->Get(name));        
    }
    return OutHist;
}


void HistogramGrabber(vector<TFile *> * Files,  vector<TH1 *> * inputHistVec, TString dataplot, TString mcplot, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow) {
    int FileSize = Files->size();
    int HistSize = inputHistVec->size();
    int NBinsX, NBinsY, NBinsZ;
    if (FileSize != HistSize) {
        cout << "Ruh Roh, Raggy, FileVector Size and HistVectorSize are not the same!" << endl;
        return;
    }
    for (int k = 0; k < FileSize; ++k) {
        if (k < 3) {
            inputHistVec->at(k) = (TH1*) Files->at(k)->Get(dataplot);
        }
        else {
            inputHistVec->at(k) = (TH1*) Files->at(k)->Get(mcplot);            
        }
        NBinsX = inputHistVec->at(k)->GetNbinsX();
        NBinsY = inputHistVec->at(k)->GetNbinsY();
        NBinsZ = inputHistVec->at(k)->GetNbinsZ();
        cout << "name " << inputHistVec->at(k)->GetName() << endl;
        inputHistVec->at(k)->Rebin(RBNX);
        if (NBinsY > 1) {
            //at least a 2D histogram
            if (NBinsZ > 1) {
                //a 3D histogram
            }
            else {
                //just a 2D histogram
            }
        }
        else {
            //just a 1D histogram
            
        }
        /*
        if (doOverflow) {
            inputHistVec->SetBinContent(
        }
        if (inputHistVec->at(k)->GetNbinsY() > 1) {
            if (doOverflow) {
                inputHistVec->at(k)
            }
        }
         */
    }
}
void HistogramAdder(vector<TH1 *> * inputHistVec, TH1 * &DataComp, TH1 * &MCComp, TH1 * &TTbarSigComp, TH1 * &TTbarBkgComp, TH1 * &VVComp, TH1 * &SingTopComp, TH1 * &WJComp, TH1 * &ZDYComp, TH1 * &QCDComp, TH1 * &StopComp) {
//    cout << " test " << endl;
    DataComp = (TH1 *) inputHistVec->at(0)->Clone(inputHistVec->at(0)->GetName() + TString("_DataComp"));
//    cout << " test 2" << endl;
    MCComp = (TH1 *) inputHistVec->at(3)->Clone(inputHistVec->at(3)->GetName() + TString("_MCComp"));
    TTbarSigComp = (TH1 *) inputHistVec->at(3)->Clone(inputHistVec->at(3)->GetName() + TString("_TTbarSig"));
    TTbarBkgComp = (TH1 *) inputHistVec->at(4)->Clone(inputHistVec->at(4)->GetName() + TString("_TTbarBkg"));
//    cout << " test 3" << endl;
    VVComp = (TH1 *) inputHistVec->at(5)->Clone(inputHistVec->at(5)->GetName() + TString("_VV"));
    SingTopComp = (TH1 *) inputHistVec->at(8)->Clone(inputHistVec->at(8)->GetName() + TString("_SingTop"));
    WJComp = (TH1 *) inputHistVec->at(9)->Clone(inputHistVec->at(10)->GetName() + TString("_WJ"));
//    cout << " test 4" << endl;
    ZDYComp = (TH1 *) inputHistVec->at(10)->Clone(inputHistVec->at(11)->GetName() + TString("_ZDY"));
    QCDComp = (TH1 *) inputHistVec->at(11)->Clone(inputHistVec->at(12)->GetName() + TString("_QCD"));
    StopComp = (TH1 *) inputHistVec->at(12)->Clone(inputHistVec->at(13)->GetName() + TString("_Stop"));
    VVComp->Add(inputHistVec->at(6));
    VVComp->Add(inputHistVec->at(7));
    for (int j = 4; j < 12; ++j) {
        MCComp->Add(inputHistVec->at(j));
    }
}
void HistogramScaler(TH1 ** inputHists, int numInputs, int startPoint, float * weights) {
    for (int i = startPoint; i < numInputs; ++i) {
//        cout << "name " << inputHists[i]->GetName() << endl;
//        cout << "weight: " << weights[i] << endl;
        inputHists[i]->Scale(weights[i]);
    }
}
void HistogramGrabber(TFile ** Files,  TH1 ** inputHists, int numInputs, TString dataplot, TString mcplot, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow) {
    int NBinsX, NBinsY, NBinsZ;
    for (int k = 0; k < numInputs; ++k) {
//        cout << "Files[k] Name: " << Files[k]->GetName() << endl;
        
//        cout << "k " << k << endl;
        if (k < 3) {
            inputHists[k] = (TH1*) Files[k]->Get(dataplot);
        }
        else {
            inputHists[k] = (TH1*) Files[k]->Get(mcplot);
        }
        NBinsX = inputHists[k]->GetNbinsX();
        NBinsY = inputHists[k]->GetNbinsY();
        NBinsZ = inputHists[k]->GetNbinsZ();
//        cout << "file name " << Files[k]->GetName() << endl;
//        cout << "hist name " << inputHists[k]->GetName() << endl;
//        cout << "hist int " << inputHists[k]->Integral() << endl;
        inputHists[k]->Rebin(RBNX);
        if (NBinsY > 1) {
            //at least a 2D histogram
            if (NBinsZ > 1) {
                //a 3D histogram
            }
            else {
                //just a 2D histogram
            }
        }
        else {
            //just a 1D histogram
            
        }
//        cout << "wtf " << endl;
        /*
         if (doOverflow) {
         inputHistVec->SetBinContent(
         }
         if (inputHistVec->at(k)->GetNbinsY() > 1) {
         if (doOverflow) {
         inputHistVec->at(k)
         }
         }
         */
    }
//    cout << "end of grabbing " << endl;
}
float ScaleBackCalc(TH1 * nVtxOrigHist, TH1 * nVtxNewHist) {
    int NBinsX = nVtxOrigHist->GetNbinsX();
    return (float) nVtxOrigHist->Integral(1, NBinsX + 1)/nVtxNewHist->Integral(1, NBinsX + 1);
}
void CDF(TH1 * inputHist, TH1 * &CDFHist) {
    int NBinsX = inputHist->GetNbinsX();
    float TotInt = inputHist->Integral(0, NBinsX+1);
//    cout << "TotInt " << TotInt << endl;
    float currInt;
    TString CDFName = inputHist->GetName();
    CDFName += "_CDF";
    CDFHist = (TH1*) inputHist->Clone(CDFName);
    for (int i = 0; i < NBinsX + 2; ++i) {
        currInt = inputHist->Integral(0, i);
//        cout << "currInt " << currInt << endl;
        CDFHist->SetBinContent(i, currInt/TotInt);
        CDFHist->SetBinError(i, 0); //Currently not setting an error on this badboy        
    }
}
void CDFRight(TH1 * inputHist, TH1 * &CDFHist) {
    int NBinsX = inputHist->GetNbinsX();
    float TotInt = inputHist->Integral(0, NBinsX+1);
//    cout << "TotInt " << TotInt << endl;
    float currInt;
    TString CDFName = inputHist->GetName();
    CDFName += "_CDF";
    CDFHist = (TH1*) inputHist->Clone(CDFName);
    for (int i = 0; i < NBinsX + 2; ++i) {
        currInt = inputHist->Integral(i, NBinsX+1);
//        cout << "currInt " << currInt << endl;
        CDFHist->SetBinContent(i, currInt/TotInt);
        CDFHist->SetBinError(i, 0); //Currently not setting an error on this badboy        
    }
}
void HistogramAdder(TH1 ** inputHists, TH1 * &DataComp, TH1 * &MCComp, TH1 * &TTbarComp, TH1 * &VVComp, TH1 * &SingTopComp, TH1 * &WJComp, TH1 * &ZDYComp, TH1 * &QCDComp, TH1 * &StopComp, int whichNTuple) {
    int endJ = 11;
//    cout << "test A" << endl;
    DataComp = (TH1 *) inputHists[0]->Clone(inputHists[0]->GetName() + TString("_DataComp"));
//    cout << "test B" << endl;
    MCComp = (TH1 *) inputHists[3]->Clone(inputHists[3]->GetName() + TString("_MCComp"));
//    cout << "test C" << endl;
    TTbarComp = (TH1 *) inputHists[3]->Clone(inputHists[3]->GetName() + TString("_TTbar"));
//    cout << "test D" << endl;
    VVComp = (TH1 *) inputHists[5]->Clone(inputHists[5]->GetName() + TString("_VV"));
///    cout << "test F" << endl;
    SingTopComp = (TH1 *) inputHists[8]->Clone(inputHists[8]->GetName() + TString("_SingTop"));
//    cout << "test G" << endl;
    WJComp = (TH1 *) inputHists[9]->Clone(inputHists[9]->GetName() + TString("_WJ"));
//    cout << "test H" << endl;
    ZDYComp = (TH1 *) inputHists[10]->Clone(inputHists[10]->GetName() + TString("_ZDY"));
//    cout << "test I" << endl;
    if (whichNTuple == 0) {
        StopComp = (TH1 *) inputHists[endJ]->Clone(inputHists[endJ]->GetName() + TString("_Stop"));        
//        cout << "test J" << endl;
    }
    else {
        QCDComp = (TH1 *) inputHists[11]->Clone(inputHists[11]->GetName() + TString("_QCD"));
        endJ = 12;
    }


    DataComp->Add(inputHists[1]);
    DataComp->Add(inputHists[2]);
    TTbarComp->Add(inputHists[4]);
    ///////DONT FORGET TO TAKE THIS OUT!!!
//    inputHists[6]->Scale(22.44/33.21);
//    inputHists[7]->Scale(9.03/17.654);
    VVComp->Add(inputHists[6]);
    VVComp->Add(inputHists[7]);
    for (int j = 4; j < endJ; ++j) {
        MCComp->Add(inputHists[j]);
    }
}
void HistogramAdder_v2(TH1 ** inputHists, TH1 * &DataComp, TH1 * &MCComp, TH1 * &TTbarComp, TH1 * &VVComp, TH1 * &SingTopComp, TH1 * &WJComp, TH1 * &ZDYComp, TH1 * &VGComp, TH1 * &HiggsComp, TH1 * &StopComp) {
    int endJ = 12;
    //    cout << "test A" << endl;
    float scaleFactor = (8.30986e+06)/(1.23816e+06);
    //    float scaleFactor = 1;
    DataComp = (TH1 *) inputHists[0]->Clone(inputHists[0]->GetName() + TString("_DataComp"));
    DataComp->Scale(scaleFactor);
    //    cout << "test B" << endl;
    MCComp = (TH1 *) inputHists[3]->Clone(inputHists[3]->GetName() + TString("_MCComp"));
    //    cout << "test C" << endl;
    TTbarComp = (TH1 *) inputHists[3]->Clone(inputHists[3]->GetName() + TString("_TTbarSig"));
    //    cout << "test D" << endl;
    VVComp = (TH1 *) inputHists[4]->Clone(inputHists[4]->GetName() + TString("_VV"));
    ///    cout << "test F" << endl;
    SingTopComp = (TH1 *) inputHists[7]->Clone(inputHists[7]->GetName() + TString("_SingTop"));
//    SingTopComp = (TH1 *) inputHists[11]->Clone(inputHists[11]->GetName() + TString("_SingTop"));
    //    cout << "test G" << endl;
    WJComp = (TH1 *) inputHists[8]->Clone(inputHists[8]->GetName() + TString("_WJ"));
    //    cout << "test H" << endl;
    ZDYComp = (TH1 *) inputHists[9]->Clone(inputHists[9]->GetName() + TString("_ZDY"));
    VGComp = (TH1 *) inputHists[10]->Clone(inputHists[10]->GetName() + TString("_VG"));
    HiggsComp = (TH1 *) inputHists[11]->Clone(inputHists[11]->GetName() + TString("_Higgs"));
    /*
     for (int i = 0; i < 12; ++i) {
         cout << "nbins for input = " << i << " is " << inputHists[i]->GetNbinsX() << endl;
         cout << "integral for input = " << i << " is " << inputHists[i]->Integral() << endl;
     }
     */
    //    cout << "test I" << endl;
    /*
     if (whichNTuple == 0) {
     StopComp = (TH1 *) inputHists[endJ]->Clone(inputHists[endJ]->GetName() + TString("_Stop"));        
     //        cout << "test J" << endl;
     }
    else {
        QCDComp = (TH1 *) inputHists[11]->Clone(inputHists[11]->GetName() + TString("_QCD"));
        endJ = 12;
    }
    */
    /*
     for (int i = 0; i < endJ; ++i) {
     cout << "integral of i = " << i << " " << inputHists[i]->Integral() << endl;
     }
    */
    DataComp->Add(inputHists[1]);
    DataComp->Add(inputHists[2]);
    VVComp->Add(inputHists[5]);
    VVComp->Add(inputHists[6]);
    for (int j = 4; j < endJ; ++j) {
        /*
        if (j == 7) {
            MCComp->Add(inputHists[11]);
        }
        else {
            MCComp->Add(inputHists[j]);
        }
        */
        MCComp->Add(inputHists[j]);
    }
}
void DataMCSpectrumDrawSplit(TCanvas * InputCanvas, TH1F * DataHist, THStack * MCStack, TH1F * MCHist, TH1F * ErrHist, TH1F * FracRatioHist, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB,  TString dataLegendComp, TString * mcLegends, TH1F ** indMCHists, int numIndMC, bool AbsRatio, bool doMean, TString cutUsed, float inputLumi) {
    TLatex * tl = new TLatex();
    TLegend * leg;
    tl->SetTextAlign(12);
    tl->SetNDC();
    tl->SetTextSize(0.03);
    char buf[99];        
    InputCanvas->Divide(1, 2);
    TPad * Pad1 = (TPad *) InputCanvas->cd(1);
    FixPad(Pad1, 1, InputCanvas);
    Pad1->SetLogy();
    ErrHist->Draw("e2");
    MCStack->Draw("hist same");
    ErrHist->Draw("e2 same");
    DataHist->Draw("e1 same");
    DataHist->Draw("axis same");
    int NBins = MCHist->GetNbinsX();
    TAxis * XAxis = ErrHist->GetXaxis();
    float XBinUB = XAxis->GetXmax();
    float XBinLB = XAxis->GetXmin();
    float BinWidthGeV = (XBinUB - XBinLB)/NBins;
    BinWidthGeV *= 1000;
    BinWidthGeV += 0.5;
    BinWidthGeV = floor(BinWidthGeV);
    BinWidthGeV /= 1000.0;
    TString BinWidthString;
    BinWidthString += BinWidthGeV;
//    cout << "BinWidthString " << BinWidthString << endl;
    TAxis * YAxis = ErrHist->GetYaxis();
    TString YAxisTitle = YAxis->GetTitle();    
    int StrPos = YAxisTitle.Index("NUM", 3, TString::kExact);
//    cout << "StrPos " << StrPos << endl;
    if (StrPos != -1) YAxisTitle.Replace(StrPos, 3, BinWidthString); 
//    TString YAxisTitle = "Entries / ";
//    TString YAxisTitle = "Entries / ";
    TString XAxisTitle = XAxis->GetTitle();
//    YAxisTitle += BinWidthGeV;
//    YAxisTitle += " GeV";
    float HistMax = TMath::Max(DataHist->GetMaximum(), MCHist->GetMaximum());
    if (YAxisUB > 1) { YAxisUB = HistMax * 2.5E2; }
    else { YAxisUB = HistMax * 2.5; Pad1->SetLogy(0);}
    HistAxisAttSet(YAxis, YAxisTitle, .06, 1.5, .05, .007, YAxisLB, YAxisUB);
    HistAxisAttSet(XAxis, XAxisTitle, YAxis->GetTitleSize(), 999, YAxis->GetLabelSize(), 999, 0.0, 0.0);
    YAxis->SetNdivisions(5,5,0);
    //    cout << "Y Axis Title " << YAxis->GetTitle() << endl;
    TString DataName = TString(DataHist->GetName());
    float MeanXstart = TextXPos;
    float MeanYstart = TextYStartPos;
    float legXstart, legYstart;
    legXstart = TextXPos;
    legYstart = TextYStartPos;
    if (TextXPos < 0.5) {
        MeanXstart += 0.01;
        MeanYstart -= 0.02;
    }
    if (doMean) {
        legYstart = MeanYstart - 0.25;
        tl->SetTextColor(kBlack);
        sprintf(buf,"Mean = %0.2f #pm %0.2f",DataHist->GetMean(), DataHist->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",DataHist->GetRMS(), DataHist->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.05,buf);
        tl->SetTextColor(kRed);
        sprintf(buf,"Mean = %0.2f #pm %0.2f",MCHist->GetMean(),MCHist->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.1,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",MCHist->GetRMS(),MCHist->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.15,buf);    
    }
//    leg = new TLegend(legXstart,legYstart - 0.05 * (numIndMC + 1),legXstart + 0.05 * (numIndMC + 1),legYstart);
    leg = new TLegend(legXstart,legYstart - 0.05 * (numIndMC + 1),legXstart + 0.40,legYstart);
    leg->AddEntry(DataHist,TString(dataLegendComp),"pl");
    for (int j = 0; j < numIndMC; ++j) {
//        cout << "j " << j << endl;
//        cout << "indMCHist integral " << indMCHists[j]->Integral() << endl;
//        cout << "mcLegends[j] " << mcLegends[j] << endl;
        leg->AddEntry(indMCHists[j], mcLegends[j], "f");
    }
    leg->AddEntry(ErrHist,"Sim unc","f");
    leg->Draw("same");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(31);
    /*
     latex.DrawLatex(0.90,0.96,Form(" %.1f fb^{-1} at #sqrt{s} = 8 TeV",intLumi/1000.));
     latex.SetTextAlign(11); // align left
     latex.DrawLatex(0.14,0.96,"CMS preliminary 2012");
     */
    latex.DrawLatex(0.886,0.862,Form(" %.1f fb^{-1} at #sqrt{s} = 8 TeV",inputLumi/1000.));
    latex.SetTextAlign(11); // align left
    latex.DrawLatex(0.6,0.96,"CMS preliminary 2012");
    TLatex photLatex;
    photLatex.SetNDC();
    photLatex.SetTextSize(0.04);
    photLatex.SetTextAlign(11);
    photLatex.DrawLatex(0.20,0.96,cutUsed);
    Pad1->Update();
    Pad1->Modified();
    //    TPad * Pad2 = new TPad();
    TPad * Pad2 = (TPad *) InputCanvas->cd(2);
    FixPad(Pad2, 2, InputCanvas);
    //    Pad2->SetPad(.005, .005, .995, .2495);
    Pad2->SetGridy(1);
    YAxis = FracRatioHist->GetYaxis();
    XAxis = FracRatioHist->GetXaxis();
    //    cout << XAxis->GetTitle() << endl;
    if (!AbsRatio) {
        HistAxisAttSet(YAxis, TString("(MC-data)/data"), .15, .54, .14, .011, -1.5, 1.5);
    }
    else {
        HistAxisAttSet(YAxis, TString("data/MC"), .15, .54, .14, .011, 0.4, 1.6); 
    }
    HistAxisAttSet(XAxis, XAxisTitle, .17, 1.03, .12,.07, 0.0, 0.0);
    YAxis->SetNdivisions(3,5,0);
    XAxis->SetNdivisions(6,5,0);
    FracRatioHist->SetLineColor(kBlack);
    FracRatioHist->Draw("e1");
    //    ratioCompfitres = fracratioComp->Fit(fit_ratioComp, "S", "axis same", 0., 70.);
}


void OviDESYMCSpectrumDrawSplit(TCanvas * InputCanvas, TH1F * DESYHist, TH1F * OviHist, TH1F * FracRatioHist, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB,  TString DESYLegendComp, TString OviLegendComp, bool AbsRatio, bool doMean, TString cutUsed, float inputLumi) {
    TLatex * tl = new TLatex();
    TLegend * leg;
    tl->SetTextAlign(12);
    tl->SetNDC();
    tl->SetTextSize(0.03);
    char buf[99];        
    InputCanvas->Divide(1, 2);
    TPad * Pad1 = (TPad *) InputCanvas->cd(1);
    FixPad(Pad1, 1, InputCanvas);
    Pad1->SetLogy();
    DESYHist->Draw("e1");
    OviHist->Draw("e1 same");
    int NBins = DESYHist->GetNbinsX();
    TAxis * XAxis = DESYHist->GetXaxis();
    TAxis * YAxis = DESYHist->GetYaxis();
    float XBinUB = XAxis->GetXmax();
    float XBinLB = XAxis->GetXmin();
    float BinWidthGeV = (XBinUB - XBinLB)/NBins;
    BinWidthGeV *= 100;
    BinWidthGeV += 0.5;
    BinWidthGeV = floor(BinWidthGeV);
    BinWidthGeV /= 100;
    TString BinWidthString;
    BinWidthString += BinWidthGeV;
    TString YAxisTitle = YAxis->GetTitle();    
    int StrPos = YAxisTitle.Index("NUM", 3, TString::kExact);
    if (StrPos != -1) YAxisTitle.Replace(StrPos, 3, BinWidthString); 
//    TString YAxisTitle = "Entries / ";
    TString XAxisTitle = XAxis->GetTitle();
//    YAxisTitle += BinWidthGeV;
//    YAxisTitle += " GeV";
    float HistMax = TMath::Max(DESYHist->GetMaximum(), OviHist->GetMaximum());
    if (YAxisUB > 1) { YAxisUB = HistMax * 2.5E2; }
    else { YAxisUB = HistMax * 2.5; Pad1->SetLogy(0);}
    HistAxisAttSet(YAxis, YAxisTitle, .06, 1.5, .05, .007, YAxisLB, YAxisUB);
    HistAxisAttSet(XAxis, XAxisTitle, YAxis->GetTitleSize(), 999, YAxis->GetLabelSize(), 999, 0.0, 0.0);
    YAxis->SetNdivisions(5,5,0);
    //    cout << "Y Axis Title " << YAxis->GetTitle() << endl;
    TString DESYName = TString(DESYHist->GetName());
    float MeanXstart = TextXPos;
    float MeanYstart = TextYStartPos;
    float legXstart, legYstart;
    legXstart = TextXPos;
    legYstart = TextYStartPos;
    if (TextXPos < 0.5) {
        MeanXstart += 0.01;
        MeanYstart -= 0.02;
    }
    if (doMean) {
        legYstart = MeanYstart - 0.25;
        tl->SetTextColor(DESYHist->GetMarkerColor());
        sprintf(buf,"Mean = %0.2f #pm %0.2f",DESYHist->GetMean(), DESYHist->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",DESYHist->GetRMS(), DESYHist->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.05,buf);
        tl->SetTextColor(OviHist->GetMarkerColor());
        sprintf(buf,"Mean = %0.2f #pm %0.2f",OviHist->GetMean(),OviHist->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.1,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",OviHist->GetRMS(),OviHist->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.15,buf);    
    }
    leg = new TLegend(legXstart,legYstart-0.25,legXstart+0.25,legYstart);
    leg->AddEntry(DESYHist,TString(DESYLegendComp),"pl");
    leg->AddEntry(OviHist,TString(OviLegendComp),"pl");
    leg->Draw("same");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.886,0.862,Form("MC comparison: %.1f fb^{-1} at #sqrt{s} = 8 TeV",inputLumi/1000.));
    latex.SetTextAlign(11); // align left
    latex.DrawLatex(0.6,0.96,"CMS preliminary 2012");
    TLatex photLatex;
    photLatex.SetNDC();
    photLatex.SetTextSize(0.04);
    photLatex.SetTextAlign(11);
    photLatex.DrawLatex(0.20,0.96,cutUsed);
    Pad1->Update();
    Pad1->Modified();
    TPad * Pad2 = (TPad *) InputCanvas->cd(2);
    FixPad(Pad2, 2, InputCanvas);
    Pad2->SetGridy(1);
    YAxis = FracRatioHist->GetYaxis();
    XAxis = FracRatioHist->GetXaxis();
    if (!AbsRatio) {
        HistAxisAttSet(YAxis, TString("(O-D)/D"), .15, .54, .14, .011, -1.5, 1.5);
    }
    else {
        HistAxisAttSet(YAxis, TString("DESY/Ovi"), .15, .54, .14, .011, 0.4, 1.6); 
    }
    HistAxisAttSet(XAxis, XAxisTitle, .17, 1.03, .12,.07, 0.0, 0.0);
    YAxis->SetNdivisions(3,5,0);
    XAxis->SetNdivisions(6,5,0);
    FracRatioHist->SetLineColor(kBlack);
    FracRatioHist->Draw("e1");
}
