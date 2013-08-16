#include <vector>
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
using namespace std;

void pointSystErr(Double_t centVal, Double_t shiftUp, Double_t shiftDown, Double_t &errUp, Double_t &errDown, bool doSymError, int justOneSideErr) {
    Double_t diffShiftUp, diffShiftDown;
    Double_t errSum;
    int whichErrCase = 0; 
    // +(-) 3 for if shift up is greater (smaller) than cent Val; likewise for shift down but with +(-) 1 instead;
    // +4 (-4): both more (less) than cent val
    // +2 (-2): upErr more, downErr less (upErr less, downErr more)
    whichErrCase += (shiftUp > centVal) ? 3 : 0;
    whichErrCase += (shiftUp < centVal) ? -3 : 0;
    whichErrCase += (shiftDown > centVal) ? 1 : 0;
    whichErrCase += (shiftDown < centVal) ? -1 : 0;
    diffShiftUp = centVal - shiftUp;
    diffShiftDown = centVal - shiftDown;
    cout << "justOneSideErr " << justOneSideErr << endl;
    cout << "whichErrCase " << whichErrCase << endl;
    switch (justOneSideErr) {
        case 0:
            errSum = TMath::Sqrt(0.5*(diffShiftDown*diffShiftDown + diffShiftUp*diffShiftUp));
            break;
        case 1:
            errSum  = abs(diffShiftUp);
            break;              
        case -1:
            errSum  = abs(diffShiftDown);
            break; 
        default:
            break;
    }
    if (doSymError) {
        errUp = errSum;
        errDown = errSum;
    }
    else {
        switch (whichErrCase) {
            case 0:
                errUp = 0;
                errDown = 0;
                break;
            case 4:
                errUp = errSum;
                errDown = 0;
                break;
            case -4:
                errUp = 0;
                errDown = errSum;
                break;
            case 2:
                errUp = abs(diffShiftUp);
                errDown = abs(diffShiftDown);
                if (justOneSideErr == 1) errDown = errUp;
                if (justOneSideErr == -1) errUp = errDown;
                break;
            case -2:
                errUp = abs(diffShiftDown);
                errDown = abs(diffShiftUp);
                if (justOneSideErr == 1) errDown = errUp;
                if (justOneSideErr == -1) errUp = errDown;
                break;
                // special bit just for MT2ll
            case 3:
                errUp = abs(diffShiftUp);
                errDown = errUp;
                break;
            case -3:
                errUp = abs(diffShiftUp);
                errDown = errUp;
                break;
                // special bit just for MT2ll
            default:
                break;
        }
    }
    //complete hack to handle MT2ll
    if (justOneSideErr > 0) {
        errUp = abs(diffShiftUp);
        errDown = abs(diffShiftUp);
    }
}
TGraphAsymmErrors * GraphSystErrorSet_SingleSource(TH1 * HistCentralVal, TH1 * HistShiftUp, TH1 * HistShiftDown, bool doSymError, int justOneSideErr) {
    int nBinsX = HistCentralVal->GetNbinsX();
    Double_t errUp, errDown;
    Double_t centVal, shiftUp, shiftDown;
    Double_t binWidth;
    Double_t x;
    TString name = HistCentralVal->GetName();
    TGraphAsymmErrors * outGraph = new TGraphAsymmErrors(nBinsX+2); 
    outGraph->SetName(name + "_ErrSingleSource");
    for (int i = 0; i < nBinsX+2; ++i) {

        //coordinate prep
        x = (Double_t) HistCentralVal->GetBinCenter(i);
        if (i == 1) x = (Double_t) HistCentralVal->GetBinLowEdge(i);
        if (i == nBinsX) x = (Double_t) HistCentralVal->GetBinLowEdge(i+1);
        binWidth = (Double_t) HistCentralVal->GetXaxis()->GetBinWidth(i);

        //initial point sets -- y-axis info will contain relevant uncertainties
        outGraph->SetPointEXlow(i, binWidth/2.);
        outGraph->SetPointEXhigh(i, binWidth/2.);
        outGraph->SetPoint(i, x, HistCentralVal->GetBinContent(i));
        
        //get info for setting y-axis
        centVal = (Double_t) HistCentralVal->GetBinContent(i);
        shiftUp = (Double_t) HistShiftUp->GetBinContent(i);
        shiftDown = (Double_t) HistShiftDown->GetBinContent(i);
        pointSystErr(centVal, shiftUp, shiftDown, errUp, errDown, doSymError, justOneSideErr);
        cout << " for i = " << i << " in syst " << HistShiftDown->GetName() << " errors " << endl;
        cout << "upErr: " << errUp << endl;
        cout << "downErr: " << errDown << endl;
        cout << endl;
        if (doSymError) {
            outGraph->SetPointEYlow(i, errDown);
            outGraph->SetPointEYhigh(i, errUp);
        }
        else {
            outGraph->SetPointEYhigh(i, errUp);
            outGraph->SetPointEYlow(i, errDown);
        }
    }
    return outGraph;
}

TGraphAsymmErrors * GraphSystErrorSet_SingleSource(TH1 * HistCentralVal, vector<TH1F *> * systHistVec, TString systString, bool doSymError, int justOneSideErr) {
    // overloaded version of the graph error set used for grabbing specific systematic
    int nBinsX = HistCentralVal->GetNbinsX();
    Double_t errUp, errDown;
    Double_t centVal, shiftUp, shiftDown;
    Double_t binWidth;
    Double_t x;
    bool     filledSystUpHist, filledSystDownHist;
    TString name = HistCentralVal->GetName();
    TGraphAsymmErrors * outGraph = new TGraphAsymmErrors(nBinsX+2); 
    outGraph->SetName(name + systString);// "_ErrSingleSource");
    TH1F * HistShiftUp, * HistShiftDown;
    TString systHistName;
    //picking out correct histogram
    cout << "systVec size " << systHistVec->size() << endl;
    bool foundUp = 0;
    bool foundDown = 0;
    int singleShiftIndex = -99;
    for (unsigned int k = 0; k < systHistVec->size(); ++k) {
        systHistName = systHistVec->at(k)->GetName();
        cout << "histName " << systHistName << endl;
        cout << "systHistName.Contains(systString) for systString " << systString << " is " <<systHistName.Contains(systString) << endl;
        if (!systHistName.Contains(systString)) {
            continue;   
        }
        else {
            singleShiftIndex = k;
        }
        cout << "systHistName.Contains(up) " << systHistName.Contains("Up") << endl;
        cout << "systHistName.Contains(down) " << systHistName.Contains("Down") << endl;
        if (systHistName.Contains("Up")) {
            HistShiftUp = systHistVec->at(k);   
            foundUp = 1;
        }
        else if (systHistName.Contains("Down")) {
            HistShiftDown = systHistVec->at(k);
            foundDown = 1;
        }
    }
    cout << "found up " << foundUp << endl;
    cout << "found down " << foundDown << endl;
//    cout << "singleshift index " << singleShiftIndex << endl;
    if (!foundUp && !foundDown) {
        if (singleShiftIndex >= 0) {
        HistShiftUp = systHistVec->at(singleShiftIndex);
        }
        else {
            HistShiftUp = (TH1F *) HistCentralVal;
        }
        HistShiftDown = (TH1F *) HistCentralVal;//->Clone(HistCentralVal->GetName() + systString + TString("Down"));
    }
    else if (!foundUp) {
        HistShiftUp = (TH1F *) HistCentralVal;//->Clone(HistCentralVal->GetName() + systString + TString("Up"));
    }
    else if (!foundDown) {
        HistShiftDown = (TH1F *) HistCentralVal;//->Clone(HistCentralVal->GetName() + systString + TString("Down"));
    }
    cout << "HACK JOB " << endl;
//    cout << "HistShiftUp " << HistShiftUp->GetName() << " has " << HistShiftUp->GetNbinsX() << endl;
//    cout << "HistShiftDown " << HistShiftDown->GetName() << " has " << HistShiftDown->GetNbinsX() << endl;
    filledSystUpHist = (HistShiftUp->Integral() > 0) ? true : false;
    filledSystDownHist = (HistShiftDown->Integral() > 0) ? true : false;
    for (int i = 0; i < nBinsX+2; ++i) { 

        //coordinate prep
        x = (Double_t) HistCentralVal->GetBinCenter(i);
//        if (i == 1) x = (Double_t) HistCentralVal->GetBinLowEdge(i);
//        if (i == nBinsX) x = (Double_t) HistCentralVal->GetBinLowEdge(i+1);
        binWidth = (Double_t) HistCentralVal->GetXaxis()->GetBinWidth(i);
        //initial point sets -- y-axis info will contain relevant uncertainties
        if (i > 0 && i < nBinsX + 1) {
            cout << "HistCentralVal->GetBinContent(i) for i: " << i << " is " << HistCentralVal->GetBinContent(i) << endl;
            cout << "HistShiftUp->GetBinContent(i) for i: " << i << " is " << HistShiftUp->GetBinContent(i) << endl;
            cout << "HistShiftDown->GetBinContent(i) for i: " << i << " is " << HistShiftDown->GetBinContent(i) << endl;
            cout << "diffShift Up " << HistShiftUp->GetBinContent(i) - HistCentralVal->GetBinContent(i) << endl;
            cout << "diffShift Down " << HistShiftDown->GetBinContent(i) - HistCentralVal->GetBinContent(i) << endl;
        }
        outGraph->SetPointEXlow(i, binWidth/2.);
        outGraph->SetPointEXhigh(i, binWidth/2.);
        outGraph->SetPoint(i, x, HistCentralVal->GetBinContent(i));
        outGraph->SetPointEXlow(i, binWidth/2.);

        //get info for setting y-axis
        centVal = (Double_t) HistCentralVal->GetBinContent(i);
        shiftUp = (Double_t) HistShiftUp->GetBinContent(i);
        shiftDown = (Double_t) HistShiftDown->GetBinContent(i);
//        cout << "within GraphSyst line 150" << endl;
        pointSystErr(centVal, shiftUp, shiftDown, errUp, errDown, doSymError, justOneSideErr);
        if (i > 0 && i < nBinsX + 1) {
            cout << " for i = " << i << " in syst " << systString << " errors " << endl;
            cout << "upErr: " << errUp << endl;
            cout << "downErr: " << errDown << endl;            
            cout << endl;
        }
        if (filledSystUpHist) {
            outGraph->SetPointEYhigh(i, errUp);
            outGraph->SetPointEYlow(i, errDown);
        }
        else {
            outGraph->SetPointEYhigh(i, 0);
            outGraph->SetPointEYlow(i, 0);
        }
        cout << "chose to fill1? " << outGraph->GetErrorYhigh(i) << endl;        
        cout << "chose to fill2? " << outGraph->GetErrorYlow(i) << endl;
    }
    return outGraph;
}

TGraphAsymmErrors * GraphSystErrorSumErrors(TGraphAsymmErrors * centValGraph, TString addName, vector<TGraphAsymmErrors *> * systGraphVec, bool addStatErr, TH1 * inputHist) {
    int nPoints = centValGraph->GetN();
//    cout << "nPoints " << nPoints << endl;
    Double_t upErrSum, downErrSum;
    Double_t upErr, downErr;
    Double_t x, y;
    TString name = centValGraph->GetName();
    TGraphAsymmErrors * outGraph = (TGraphAsymmErrors *) centValGraph->Clone(name + addName + TString("_ErrSummed"));
    for (int i = 0; i < nPoints; ++i) {
        upErrSum = 0;
        downErrSum = 0;
        x = (Double_t) inputHist->GetBinCenter(i);
        y = (Double_t) inputHist->GetBinContent(i);
        if (addStatErr) {            
            upErr = (Double_t) centValGraph->GetErrorYhigh(i);
            downErr = (Double_t) centValGraph->GetErrorYlow(i);
            upErrSum += upErr*upErr;
            downErrSum += downErr*downErr;
        }
        for (unsigned int k = 0; k < systGraphVec->size(); ++k) {
            upErr = (Double_t) systGraphVec->at(k)->GetErrorYhigh(i);
            downErr = (Double_t) systGraphVec->at(k)->GetErrorYlow(i);
            /*
             cout << "k " << i << endl;
             cout << "Graph" << GraphArray[k]->GetName() << endl;
             cout << "upErr " << upErr << endl;
             cout << "downErr " << downErr << endl;
             */
            if (i == 2) {
                cout << "systGraph Name " << systGraphVec->at(k)->GetName() << endl;
                cout << "upErr " << upErr << endl;
                cout << "downErr " << downErr << endl;
            }
            upErrSum += upErr*upErr;
            downErrSum += downErr*downErr;
        }
        outGraph->SetPoint(i, x, y);
        outGraph->SetPointEYlow(i, sqrt(downErrSum));
        outGraph->SetPointEYhigh(i, sqrt(upErrSum));
        outGraph->SetPointEXlow(i, systGraphVec->at(0)->GetErrorXlow(i));
        outGraph->SetPointEXhigh(i, systGraphVec->at(0)->GetErrorXhigh(i));
    }
    return outGraph;
}

TGraphAsymmErrors * GraphSystErrorSumErrors(TGraphAsymmErrors * centValGraph, TString addName, TGraphAsymmErrors * systGraph, TH1 * inputHist) {
    int nPoints = centValGraph->GetN();
    Double_t upErrSum, downErrSum;
    Double_t upErr, downErr;
    Double_t x, y;
    TString name = centValGraph->GetName();
    TGraphAsymmErrors * outGraph = (TGraphAsymmErrors *) centValGraph->Clone(name + addName + TString("_ErrSummed"));
    for (int i = 0; i < nPoints; ++i) {
        upErrSum = 0;
        downErrSum = 0;
        x = (Double_t) inputHist->GetBinCenter(i);
        y = (Double_t) inputHist->GetBinContent(i);  
        upErr = (Double_t) centValGraph->GetErrorYhigh(i);
        downErr = (Double_t) centValGraph->GetErrorYlow(i);
        upErrSum += upErr*upErr;
        downErrSum += downErr*downErr;
        upErr = (Double_t) systGraph->GetErrorYhigh(i);
        downErr = (Double_t) systGraph->GetErrorYlow(i);
        upErrSum += upErr*upErr;
        downErrSum += downErr*downErr;
        outGraph->SetPoint(i, x, y);
        outGraph->SetPointEYlow(i, sqrt(downErrSum));
        outGraph->SetPointEYhigh(i, sqrt(upErrSum));
        outGraph->SetPointEXlow(i, systGraph->GetErrorXlow(i));
        outGraph->SetPointEXhigh(i, systGraph->GetErrorXhigh(i));
    }
    return outGraph;
}
TGraphAsymmErrors * clonePoints(TH1 * inputHist, bool cloneZero = true) {
    int NBinsX = inputHist->GetNbinsX();
    double x, y;
    Double_t binWidth;
    TGraphAsymmErrors * outGraph = new TGraphAsymmErrors(NBinsX + 2);
    outGraph->SetName(inputHist->GetName() + TString("_Graph"));
    for (int i = 0; i < NBinsX + 2; ++i) {
        x = (double) inputHist->GetBinCenter(i);
        y = (double) inputHist->GetBinContent(i);
        if (!cloneZero && y == 0) continue;
        binWidth = (Double_t) inputHist->GetXaxis()->GetBinWidth(i);
        //        cout << "i, x, y: " << i << ", " << x << ", " << y << endl;
        outGraph->SetPoint(i, x, y);
        outGraph->SetPointEYlow(i, inputHist->GetBinError(i));
        outGraph->SetPointEYhigh(i, inputHist->GetBinError(i));
        outGraph->SetPointEXlow(i, binWidth/2.);
        outGraph->SetPointEXhigh(i, binWidth/2.);
    }
    return outGraph;
}

TGraphAsymmErrors * fracGraph(TH1 * centValMCHist, TGraphAsymmErrors * errGraphSystPlusStat, bool doAbsRatio, float yAxisRange) {
    double x, y;
    double rel, relUpErr, relDownErr;
    Double_t binWidth;
    int nBins = centValMCHist->GetNbinsX();
    TGraphAsymmErrors * ratioGraph = new TGraphAsymmErrors(nBins + 2);    
    GraphMainAttSet(ratioGraph, kGray+1, 3001, kGray+1, 2, kWhite, 0, 0); 
    if (doAbsRatio) {
        rel = 1;
        HistAxisAttSet(ratioGraph->GetYaxis(), TString("Data/MC"), .15, .54, .14, .011, 1.0 - yAxisRange, 1.0 + yAxisRange); 
    }
    else {
        rel = 0;
        HistAxisAttSet(ratioGraph->GetYaxis(), TString("(MC-Data)/Data"), .15, .54, .14, .011, -1.0 * yAxisRange, 1.0 * yAxisRange);
    }
    for (int i = 0; i < nBins + 2; ++i) {
        x = centValMCHist->GetBinCenter(i);
        y = centValMCHist->GetBinContent(i);
        ratioGraph->SetPoint(i, x, rel);        
        binWidth = (Double_t) centValMCHist->GetXaxis()->GetBinWidth(i);
        
        //initial point sets -- y-axis info will contain relevant uncertainties
//        cout << " in frac graph " << endl;
//        cout << "x " << x << endl;
//        cout << "y " << y << endl;

        if (!(y > 0.)) {
            relUpErr = 0;
            relDownErr = 0;
        }
        else {
            relUpErr = (doAbsRatio) ? errGraphSystPlusStat->GetErrorYlow(i) :  errGraphSystPlusStat->GetErrorYhigh(i);
            relUpErr /= centValMCHist->GetBinContent(i);
            relDownErr = (doAbsRatio) ? errGraphSystPlusStat->GetErrorYhigh(i) :  errGraphSystPlusStat->GetErrorYlow(i);
            relDownErr /= centValMCHist->GetBinContent(i);
        }
        ratioGraph->SetPointEYhigh(i, relUpErr);
        ratioGraph->SetPointEYlow(i, relDownErr);
//        cout << "error in point i for i = " << i << " is " << ratioGraph->GetErrorYlow(i) << endl;
//        cout << "error in point i for i = " << i << " is " << ratioGraph->GetErrorYhigh(i) << endl;
//        cout << "binWidth / 2 " << binWidth/2 << endl;
        ratioGraph->SetPointEXlow(i, binWidth/2.);
        ratioGraph->SetPointEXhigh(i, binWidth/2.);
        
    }
    return ratioGraph;
}
