#include <vector>
#include "./StopFunctionDefinitions_v2.h"
//#include <boost/format.hpp>
#include <sstream>
using namespace std;


vector<TString> * StopFileNames(int whichNTuple) {
    vector<TString> * outVec = new vector<TString>;
    const int numOviTypes = 12;
    TString fileInNameSpecOviedo[numOviTypes] = {"DataMu12MT2Leq80", "DataEG12MT2Leq80", "DataMuEG12MT2Leq80", "TTbar_MassiveBinDECAY", "TTbar_MassiveBinDECAYBkg", "WWTo2L2Nu_Madgraph", "WZ", "ZZ", "SingTop", "WJets_Madgraph", "ZDY", "Stop"};
        const int numDESYTypes = 12;
    TString fileInNameSpecDESY[numDESYTypes] = {"DataMuMu", "DataEMu", "DataEE", "TTBarSig", "TTBarBkg", "WW", "WZ", "ZZ", "SingleTop", "WToLNu", "ZDY", "QCD"};
    switch (whichNTuple) {
        case 0:
            for (int i = 0; i < numOviTypes; ++i) {
                outVec->push_back(fileInNameSpecOviedo[i]);
            }
            break;
        case 1:
            for (int i = 0; i < numDESYTypes; ++i) {
                outVec->push_back(fileInNameSpecDESY[i]);
            }
            break;
        default:
            break;
    }
    return outVec;
}
vector<TFile *> * StopFiles(int whichNTuple, vector<TString> * fileNames, int whichTTBarGen, bool doPURW, bool doSyst) {
    
    vector<TFile *> * outVec = new vector<TFile *>;
    TFile * outTFile;
    TString addPath, fileInNameBase, specNTupString, fileNameSuffix, fileName, TTBarGenString, PURWString, SystString;
    if (whichNTuple == 0) {
        addPath        = "";
        fileInNameBase = "TreeAnalysisTop_5311pb-1_";
        specNTupString = "_Oviedo";
        fileNameSuffix = "_Output.root";
        TTBarGenString = "";
    }
    else {
        addPath        = "../RootFiles/";
        fileInNameBase = "";
        specNTupString = "_DESY";
        fileNameSuffix =  "Haddplots.root";
        TTBarGenString = "";
        switch (whichTTBarGen) {
            case 1:
                TTBarGenString = "_mcatnlo";
                break;
            case 2:
                TTBarGenString = "_powheg";
                break;                
            default:
                break;
        }
    }
    for (unsigned int i = 0; i < fileNames->size(); ++i) {
        fileName = fileNames->at(i);
        if (fileName.Contains("TTBar")) fileName += TTBarGenString;
        if (!fileName.Contains("Data")) {
            PURWString = doPURW ? "_PURW" : "";
            SystString = doSyst ? "_wSyst" : "";
    
        }
//        outTFile = new TFile(addPath + fileInNameBase + fileName + specNTupString + PURWString + SystString + fileNameSuffix);
        outTFile = TFile::Open(addPath + fileInNameBase + fileName + specNTupString + PURWString + SystString + fileNameSuffix);
        outVec->push_back(outTFile);
    }
    return outVec;
}

vector<TString> * MCLegends(int whichNTuple, bool addThings) {
    vector<TString> * outVec = new vector<TString>;
    TString mcLegendsOviAdd[6] = {"W + Jets", "VV", "Z/#gamma* #rightarrow l^{+}l^{-}", "Single Top", "t#bar{t}", "Stop Mass: "};
    TString mcLegendsDESYAdd[6] = {"W + Jets", "VV", "Z/#gamma* #rightarrow l^{+}l^{-}", "Single Top", "t#bar{t}", "Multijets"};
    TString mcLegendsDESY[9] = {"W + Jets", "WW", "WZ", "ZZ", "Z/#gamma* #rightarrow l^{+}l^{-}", "Single Top", "t#bar{t} sig", "t#bar{t} bkg", "Multijets"};
    switch (whichNTuple) {
        case 0:
            if (addThings) {
                for (int i = 0; i < 6; ++i) {
                    outVec->push_back(mcLegendsOviAdd[i]);
                }
            }
            else {
                /*
                for (int i = 0; i < 6; ++i) {
                    outVec->push_back(mcLegendsOvi[i]);
                }
                */
            }
            break;
        case 1:
            if (addThings) {
                for (int i = 0; i < 6; ++i) {
                    outVec->push_back(mcLegendsDESYAdd[i]);
                }
            }
            else {
                 for (int i = 0; i < 6; ++i) {
                 outVec->push_back(mcLegendsDESY[i]);
                 }
            }
            break;            
        default:
            break;
    }
    return outVec;
}

vector<Color_t> * MCColors(int whichNTuple, bool addThings) {
    vector<Color_t> * outVec = new vector<Color_t>;
    Color_t mcColorsOviAdd[6] = {kGreen + 2, kOrange + 2, kBlue, kRed - 10, kRed, kGreen - 2};
    Color_t mcColorsDESYAdd[6] = {kGreen + 2, kOrange + 2, kBlue, kRed - 10, kRed, kCyan -4};
    Color_t mcColorsDESY[9] = {kGreen + 2, kOrange + 2, kPink+9, kPink - 8, kBlue, kRed - 10, kRed-5, kRed, kCyan -4};
    switch (whichNTuple) {
        case 0:
            if (addThings) {
                for (int i = 0; i < 6; ++i) {
                    outVec->push_back(mcColorsOviAdd[i]);
                }
            }
            else {
                /*
                 for (int i = 0; i < 6; ++i) {
                 outVec->push_back(mcColorsOvi[i]);
                 }
                 */
            }
            break;
        case 1:
            if (addThings) {
                for (int i = 0; i < 6; ++i) {
                    outVec->push_back(mcColorsDESYAdd[i]);
                }
            }
            else {
                for (int i = 0; i < 6; ++i) {
                    outVec->push_back(mcColorsDESY[i]);
                }
            }
            break;            
        default:
            break;
    }
    return outVec;
}
void HistogramVecGrabber(vector<TFile *> * inputFiles, vector<TH1F *> * dataHistVec, vector<TH1F *> * mcIndHistCentValVec, vector<TH1F *> * mcCompHistSystVec, vector<float> * nVtxBackScaleVec, vector<SystT> * systVec, TString dataPlotName, TString mcPlotName, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst) {
    cout << "inside Histogram Vec Grabber line: 144" << endl;
    int NBinsX, NBinsY, NBinsZ;
    TString fileName, systName, mcGrabName;
    TH1F * currHist;
    TH1F * systCompHist = NULL;
    cout << "inside Histogram Vec Grabber line: 149" << endl;
    cout << "inputFiles size " << inputFiles->size() << endl;
    cout << "nVtxBackScale size " << nVtxBackScaleVec->size() << endl;
    TString mcSystPlotName, systCompHistName;
    for (unsigned int k = 0; k < inputFiles->size(); ++k) {
        mcGrabName = "";
        cout << "k " << k << endl;
        fileName = inputFiles->at(k)->GetName();
        cout << "fileName" << " is " << fileName << endl;
        if (fileName.Contains("Data")) {
            cout << "contains data " << endl;
            cout << "trying to grab " << dataPlotName << endl;
            currHist = (TH1F *) inputFiles->at(k)->Get(dataPlotName);
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) {
                if (NBinsZ > 1) {
                }
                else {
                }
            }
            else {
            }
            dataHistVec->push_back(currHist);
        }
        else {
            mcGrabName += mcPlotName;
            mcGrabName += subSampName;
            currHist = (TH1F *) inputFiles->at(k)->Get(mcGrabName);
            currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) {
                if (NBinsZ > 1) {
                }
                else {
                }
            }
            else {
            }
            mcIndHistCentValVec->push_back(currHist);
        }
    }
    if (doSyst) {
        cout << "inside Histogram Vec Grabber line: 187" << endl;
        for (unsigned int j = 0; j < systVec->size(); ++j) {
            systCompHist = NULL;
            mcGrabName = "";
            mcGrabName += mcPlotName;
            cout << "j " << j << endl;
            //        if (!systVec->at(j).doXSyst) continue;
            systName = systVec->at(j).name;
            cout << "systName " << systName << endl;
            cout << "mcGrabName " << mcGrabName << endl;
            if (!mcGrabName.Contains("MT2ll") && systName.Contains("MT2ll")) continue;
            mcSystPlotName = mcGrabName;
            mcSystPlotName += systName;
            mcGrabName += subSampName;
            mcSystPlotName += subSampName;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                if (systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH1F *) inputFiles->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH1F *) inputFiles->at(k)->Get(mcSystPlotName);   
                }
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
                NBinsX = currHist->GetNbinsX();
                NBinsY = currHist->GetNbinsY();
                NBinsZ = currHist->GetNbinsZ();
                currHist->RebinX(RBNX);
                if (NBinsY > 1) {
                    if (NBinsZ > 1) {
                    }
                    else {
                    }
                }
                else {
                }
                cout << "rebinned hist " << endl;
                if (systCompHist == NULL) {
                    cout << "cloning syst " << endl;
                    systCompHistName = currHist->GetName();
                    systCompHistName += "_Comp";
                    systCompHist = (TH1F *) currHist->Clone(systCompHistName);
                }
                else {
                    cout << "adding syst " << endl;
                    cout << "systCompHist nBins " << systCompHist->GetNbinsX() << endl;
                    cout << "currHist nBins " << currHist->GetNbinsX() << endl;
                    systCompHist->Add(currHist);
                }
            }
            mcCompHistSystVec->push_back(systCompHist);
        }
    }
}
void HistogramVecGrabber(vector<TFile *> * inputFiles, vector<TH2F *> * dataHistVec, vector<TH2F *> * mcIndHistCentValVec, vector<TH2F *> * mcCompHistSystVec, vector<float> * nVtxBackScaleVec, vector<SystT> * systVec, TString dataPlotName, TString mcPlotName, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst) {
    cout << "inside Histogram Vec Grabber line: 144" << endl;
    int NBinsX, NBinsY, NBinsZ;
    TString fileName, systName, mcGrabName;
    TH2F * currHist;
    TH2F * systCompHist = NULL;
    cout << "inside Histogram Vec Grabber line: 149" << endl;
    cout << "inputFiles size " << inputFiles->size() << endl;
    cout << "nVtxBackScale size " << nVtxBackScaleVec->size() << endl;
    TString mcSystPlotName, systCompHistName;
    for (unsigned int k = 0; k < inputFiles->size(); ++k) {
        mcGrabName = "";
        cout << "k " << k << endl;
        fileName = inputFiles->at(k)->GetName();
        cout << "fileName" << " is " << fileName << endl;
        if (fileName.Contains("Data")) {
            cout << "contains data " << endl;
            cout << "trying to grab " << dataPlotName << endl;
            currHist = (TH2F *) inputFiles->at(k)->Get(dataPlotName);
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) {
                if (NBinsZ > 1) {
                }
                else {
                }
            }
            else {
            }
            dataHistVec->push_back(currHist);
        }
        else {
            mcGrabName += mcPlotName;
            mcGrabName += subSampName;
            currHist = (TH2F *) inputFiles->at(k)->Get(mcGrabName);
            currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) {
                if (NBinsZ > 1) {
                }
                else {
                }
            }
            else {
            }
            mcIndHistCentValVec->push_back(currHist);
        }
    }
    if (doSyst) {
        cout << "inside Histogram Vec Grabber line: 187" << endl;
        for (unsigned int j = 0; j < systVec->size(); ++j) {
            systCompHist = NULL;
            mcGrabName = "";
            mcGrabName += mcPlotName;
            cout << "j " << j << endl;
            //        if (!systVec->at(j).doXSyst) continue;
            systName = systVec->at(j).name;
            cout << "systName " << systName << endl;
            cout << "mcGrabName " << mcGrabName << endl;
            if (!mcGrabName.Contains("MT2ll") && systName.Contains("MT2ll")) continue;
            mcSystPlotName = mcGrabName;
            mcSystPlotName += systName;
            mcGrabName += subSampName;
            mcSystPlotName += subSampName;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                if (systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH2F *) inputFiles->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH2F *) inputFiles->at(k)->Get(mcSystPlotName);   
                }
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
                NBinsX = currHist->GetNbinsX();
                NBinsY = currHist->GetNbinsY();
                NBinsZ = currHist->GetNbinsZ();
                currHist->RebinX(RBNX);
                if (NBinsY > 1) {
                    if (NBinsZ > 1) {
                    }
                    else {
                    }
                }
                else {
                }
                cout << "rebinned hist " << endl;
                if (systCompHist == NULL) {
                    cout << "cloning syst " << endl;
                    systCompHistName = currHist->GetName();
                    systCompHistName += "_Comp";
                    systCompHist = (TH2F *) currHist->Clone(systCompHistName);
                }
                else {
                    cout << "adding syst " << endl;
                    cout << "systCompHist nBins " << systCompHist->GetNbinsX() << endl;
                    cout << "currHist nBins " << currHist->GetNbinsX() << endl;
                    systCompHist->Add(currHist);
                }
            }
            mcCompHistSystVec->push_back(systCompHist);
        }
    }
}
void HistogramVecGrabber(vector<TFile *> * inputFiles, vector<TH3F *> * dataHistVec, vector<TH3F *> * mcIndHistCentValVec, vector<TH3F *> * mcCompHistSystVec, vector<float> * nVtxBackScaleVec, vector<SystT> * systVec, TString dataPlotName, TString mcPlotName, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst) {
    cout << "inside Histogram Vec Grabber line: 144" << endl;
    int NBinsX, NBinsY, NBinsZ;
    TString fileName, systName, mcGrabName;
    TH3F * currHist;
    TH3F * systCompHist = NULL;
    cout << "inside Histogram Vec Grabber line: 149" << endl;
    cout << "inputFiles size " << inputFiles->size() << endl;
    cout << "nVtxBackScale size " << nVtxBackScaleVec->size() << endl;
    TString mcSystPlotName, systCompHistName;
    for (unsigned int k = 0; k < inputFiles->size(); ++k) {
        mcGrabName = "";
        cout << "k " << k << endl;
        fileName = inputFiles->at(k)->GetName();
        cout << "fileName" << " is " << fileName << endl;
        if (fileName.Contains("Data")) {
            cout << "contains data " << endl;
            cout << "trying to grab " << dataPlotName << endl;
            currHist = (TH3F *) inputFiles->at(k)->Get(dataPlotName);
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) {
                if (NBinsZ > 1) {
                }
                else {
                }
            }
            else {
            }
            dataHistVec->push_back(currHist);
        }
        else {
            mcGrabName += mcPlotName;
            mcGrabName += subSampName;
            currHist = (TH3F *) inputFiles->at(k)->Get(mcGrabName);
            currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) {
                if (NBinsZ > 1) {
                }
                else {
                }
            }
            else {
            }
            mcIndHistCentValVec->push_back(currHist);
        }
    }
    if (doSyst) {
        cout << "inside Histogram Vec Grabber line: 187" << endl;
        for (unsigned int j = 0; j < systVec->size(); ++j) {
            systCompHist = NULL;
            mcGrabName = "";
            mcGrabName += mcPlotName;
            cout << "j " << j << endl;
            //        if (!systVec->at(j).doXSyst) continue;
            systName = systVec->at(j).name;
            cout << "systName " << systName << endl;
            cout << "mcGrabName " << mcGrabName << endl;
            if (!mcGrabName.Contains("MT2ll") && systName.Contains("MT2ll")) continue;
            mcSystPlotName = mcGrabName;
            mcSystPlotName += systName;
            mcGrabName += subSampName;
            mcSystPlotName += subSampName;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                if (systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH3F *) inputFiles->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH3F *) inputFiles->at(k)->Get(mcSystPlotName);   
                }
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
                NBinsX = currHist->GetNbinsX();
                NBinsY = currHist->GetNbinsY();
                NBinsZ = currHist->GetNbinsZ();
                currHist->RebinX(RBNX);
                if (NBinsY > 1) {
                    if (NBinsZ > 1) {
                    }
                    else {
                    }
                }
                else {
                }
                cout << "rebinned hist " << endl;
                if (systCompHist == NULL) {
                    cout << "cloning syst " << endl;
                    systCompHistName = currHist->GetName();
                    systCompHistName += "_Comp";
                    systCompHist = (TH3F *) currHist->Clone(systCompHistName);
                }
                else {
                    cout << "adding syst " << endl;
                    cout << "systCompHist nBins " << systCompHist->GetNbinsX() << endl;
                    cout << "currHist nBins " << currHist->GetNbinsX() << endl;
                    systCompHist->Add(currHist);
                }
            }
            mcCompHistSystVec->push_back(systCompHist);
        }
    }
}
/*
void HistogramVecGrabber(vector<TFile *> * inputFiles, vector<TH1 *> * dataHistVec, vector<TH1 *> * mcIndHistCentValVec, vector<TH1 *> * mcCompHistSystVec, int nDims, vector<float> * nVtxBackScaleVec, vector<SystT> * systVec, TString dataPlotName, TString mcPlotName, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst) {
    cout << "inside Histogram Vec Grabber line: 144" << endl;
    int NBinsX, NBinsY, NBinsZ;
    TString fileName, systName, mcGrabName;
    TH1F * currHist;
    TH1F * systCompHist = NULL;
    cout << "inside Histogram Vec Grabber line: 149" << endl;
    cout << "inputFiles size " << inputFiles->size() << endl;
    cout << "nVtxBackScale size " << nVtxBackScaleVec->size() << endl;
    TString mcSystPlotName, systCompHistName;
    for (unsigned int k = 0; k < inputFiles->size(); ++k) {
        mcGrabName = "";
        cout << "k " << k << endl;
        fileName = inputFiles->at(k)->GetName();
        cout << "fileName" << " is " << fileName << endl;
        if (fileName.Contains("Data")) {
            cout << "contains data " << endl;
            cout << "trying to grab " << dataPlotName << endl;
            switch (nDims) {
                case 1:
                    currHist = (TH1F *) inputFiles->at(k)->Get(dataPlotName);
                    break;
                case 2:
                    currHist = (TH2F *) inputFiles->at(k)->Get(dataPlotName);
                    break;
                case 3:
                    currHist = (TH3F *) inputFiles->at(k)->Get(dataPlotName);
                    break;
                default:
                    break;
            }
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) currHist->RebinY(RBNY);
            if (NBinsZ > 1) currHist->RebinZ(RBNZ);
            dataHistVec->push_back(currHist);
        }
        else {
            mcGrabName += mcPlotName;
            mcGrabName += subSampName;
            switch (nDims) {
                case 1:
                    currHist = (TH1F *) inputFiles->at(k)->Get(mcGrabName);
                    break;
                case 2:
                    currHist = (TH2F *) inputFiles->at(k)->Get(mcGrabName);
                    break;
                case 3:
                    currHist = (TH3F *) inputFiles->at(k)->Get(mcGrabName);
                    break;
                default:
                    break;
            }
            currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
            NBinsX = currHist->GetNbinsX();
            NBinsY = currHist->GetNbinsY();
            NBinsZ = currHist->GetNbinsZ();
            currHist->RebinX(RBNX);
            if (NBinsY > 1) currHist->RebinY(RBNY);
            if (NBinsZ > 1) currHist->RebinZ(RBNZ);
            mcIndHistCentValVec->push_back(currHist);
        }
    }
    if (doSyst) {
        cout << "inside Histogram Vec Grabber line: 187" << endl;
        for (unsigned int j = 0; j < systVec->size(); ++j) {
            systCompHist = NULL;
            mcGrabName = "";
            mcGrabName += mcPlotName;
            cout << "j " << j << endl;
            //        if (!systVec->at(j).doXSyst) continue;
            systName = systVec->at(j).name;
            cout << "systName " << systName << endl;
            cout << "mcGrabName " << mcGrabName << endl;
            if (!mcGrabName.Contains("MT2ll") && systName.Contains("MT2ll")) continue;
            mcSystPlotName = mcGrabName;
            mcSystPlotName += systName;
            mcGrabName += subSampName;
            mcSystPlotName += subSampName;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                switch (nDims) {
                    case 1:                        
                        if (systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                            currHist = (TH1F *) inputFiles->at(k)->Get(mcGrabName);   
                        }
                        else {
                            currHist = (TH1F *) inputFiles->at(k)->Get(mcSystPlotName);   
                        }
                        break;
                    case 2:                        
                        if (systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                            currHist = (TH2F *) inputFiles->at(k)->Get(mcGrabName);   
                        }
                        else {
                            currHist = (TH2F *) inputFiles->at(k)->Get(mcSystPlotName);   
                        }                        break;
                    case 3:                        
                        if (systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                            currHist = (TH3F *) inputFiles->at(k)->Get(mcGrabName);   
                        }
                        else {
                            currHist = (TH3F *) inputFiles->at(k)->Get(mcSystPlotName);   
                        } 
                        break;
                    default:
                        break;
                }
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
                NBinsX = currHist->GetNbinsX();
                NBinsY = currHist->GetNbinsY();
                NBinsZ = currHist->GetNbinsZ();                
                currHist->RebinX(RBNX);
                if (NBinsY > 1) currHist->RebinY(RBNY);
                if (NBinsZ > 1) currHist->RebinZ(RBNZ);
                cout << "rebinned hist " << endl;
                if (systCompHist == NULL) {
                    cout << "cloning syst " << endl;
                    systCompHistName = currHist->GetName();
                    systCompHistName += "_Comp";
                    switch (nDims) {
                        case 1:
                            systCompHist = (TH1F *) currHist->Clone(systCompHistName);
                            break;
                        case 2:
                            systCompHist = (TH2F *) currHist->Clone(systCompHistName);
                            break;
                        case 3:
                            systCompHist = (TH3F *) currHist->Clone(systCompHistName);
                            break;
                        default:
                            break;
                    }
                }
                else {
                    cout << "adding syst " << endl;
                    cout << "systCompHist nBins " << systCompHist->GetNbinsX() << endl;
                    cout << "currHist nBins " << currHist->GetNbinsX() << endl;
                    systCompHist->Add(currHist);
                }
            }
            mcCompHistSystVec->push_back(systCompHist);
        }
    }
}
 */
void HistogramAdderSyst(vector<TH1F *> * dataHistVec, vector<TH1F *> * mcIndHistCentValVec, vector<TH1F *> * mcCompHistCentValVec, TH1F * &DataComp, TH1F * &MCComp, TH1F * &FracComp, bool doAbsRatio, float yAxisRange) {
    cout << "inside Histogram Adder Syst line: 232" << endl;
    /*
    bool cloneTTBar = false;
    bool cloneVV = false;
    bool cloneSingTop = false;
    bool cloneWJ = false;
    bool cloneZDY = false;
    bool cloneStop = false;
    bool cloneQCD = false;
    */
    TString MCName;
    TH1F * TTbarComp, * VVComp, * SingTopComp, * WJComp, * ZDYComp, * QCDComp;
    TH1F * StopComp;
    cout << "inside Histogram Adder Syst line: 242" << endl;
    DataComp = (TH1F *) dataHistVec->at(0)->Clone(dataHistVec->at(0)->GetName() + TString("_DataComp"));
    MCComp = (TH1F *) mcIndHistCentValVec->at(0)->Clone(mcIndHistCentValVec->at(0)->GetName() + TString("_MCComp"));
    cout << "inside Histogram Adder Syst line: 245" << endl;
    cout << "data vec size " <<  dataHistVec->size() << endl;
    for (unsigned int i = 1; i < dataHistVec->size(); ++i) {
        cout << "DataComp nBins " << DataComp->GetNbinsX() << endl;
        cout << "dataHistVec->at(i) " << dataHistVec->at(i)->GetNbinsX() << endl;
        DataComp->Add(dataHistVec->at(i));
        cout << "i " << i << endl;
    }
    cout << "inside Histogram Adder Syst line: 249" << endl;
    cout << "MC vec size " <<  mcIndHistCentValVec->size() << endl;
    TTbarComp = (TH1F *) mcIndHistCentValVec->at(0)->Clone(mcIndHistCentValVec->at(0)->GetName() + TString("_TTbar"));
    TTbarComp->Add(mcIndHistCentValVec->at(1));
    VVComp = (TH1F *) mcIndHistCentValVec->at(2)->Clone(mcIndHistCentValVec->at(2)->GetName() + TString("_VV"));
    VVComp->Add(mcIndHistCentValVec->at(3));
    VVComp->Add(mcIndHistCentValVec->at(4));
    SingTopComp = (TH1F *) mcIndHistCentValVec->at(5)->Clone(mcIndHistCentValVec->at(5)->GetName() + TString("_SingTop"));
    WJComp = (TH1F *) mcIndHistCentValVec->at(6)->Clone(mcIndHistCentValVec->at(6)->GetName() + TString("_WJ"));
    ZDYComp = (TH1F *) mcIndHistCentValVec->at(7)->Clone(mcIndHistCentValVec->at(7)->GetName() + TString("_ZDY"));
    QCDComp = (TH1F *) mcIndHistCentValVec->at(8)->Clone(mcIndHistCentValVec->at(8)->GetName() + TString("_QCD"));
    for (unsigned int j = 0; j < mcIndHistCentValVec->size(); ++j) {
//        MCName = mcIndHistCentValVec->at(j)->GetName();
        if (j != 0) MCComp->Add(mcIndHistCentValVec->at(j));
        /*
        if (MCName.Contains("TTBar")) { 
            if (!cloneTTBar) {
                cloneTTBar = true;
                TTbarComp = (TH1F *) mcIndHistCentValVec->at(j)->Clone(mcIndHistCentValVec->at(j)->GetName() + TString("_TTbar"));
            }
            else {
                TTbarComp->Add(mcIndHistCentValVec->at(j));
            }
        }
        else if (MCName.Contains("SingTop")) { 
            if (!cloneSingTop) {
                cloneSingTop = true;
                SingTopComp = (TH1F *) mcIndHistCentValVec->at(j)->Clone(mcIndHistCentValVec->at(j)->GetName() + TString("_SingTop"));
            }
            else {
                SingTopComp->Add(mcIndHistCentValVec->at(j));
            }
        }
        else if (MCName.Contains("WJ")) { 
            if (!cloneWJ) {
                cloneWJ = true;
                WJComp = (TH1F *) mcIndHistCentValVec->at(j)->Clone(mcIndHistCentValVec->at(j)->GetName() + TString("_WJ"));
            }
            else {
                WJComp->Add(mcIndHistCentValVec->at(j));
            }
        }
        else if (MCName.Contains("ZDY")) {
            if (!cloneZDY) {
                cloneZDY = true;
                ZDYComp = (TH1F *) mcIndHistCentValVec->at(j)->Clone(mcIndHistCentValVec->at(j)->GetName() + TString("_ZDY"));
            }
            else {
                ZDYComp->Add(mcIndHistCentValVec->at(j));
            }
        }
        else if (MCName.Contains("QCD")) {            
            if (!cloneQCD) {
                cloneQCD = true;
                QCDComp = (TH1F *) mcIndHistCentValVec->at(j)->Clone(mcIndHistCentValVec->at(j)->GetName() + TString("_QCD"));
            }
            else {
                QCDComp->Add(mcIndHistCentValVec->at(j));
            }
        }
        else if ((MCName.Contains("WW") || MCName.Contains("WZ") || MCName.Contains("ZZ"))) {
            if (!cloneVV) {
                cloneVV = true;
                VVComp = (TH1F *) mcIndHistCentValVec->at(j)->Clone(mcIndHistCentValVec->at(j)->GetName() + TString("_VV"));
            }
            else {
                VVComp->Add(mcIndHistCentValVec->at(j));                
            }
        }
        else if (MCName.Contains("Stop")) {            
            if (!cloneStop) {
                cloneStop = true;
                StopComp = (TH1F *) mcIndHistCentValVec->at(j)->Clone(mcIndHistCentValVec->at(j)->GetName() + TString("_Stop"));
            }
            else {
                StopComp->Add(mcIndHistCentValVec->at(j));
            }
        }
         */
    }
    cout << "inside Histogram Adder Syst line: 317" << endl;
    if (doAbsRatio) {
        FracComp = (TH1F*) DataComp->Clone("ratioComp");
        FracComp->Divide(FracComp, MCComp, 1, 1, "");
        HistAxisAttSet(FracComp->GetYaxis(), TString("Data/MC"), .15, .54, .14, .011, 1.0 - yAxisRange, 1.0 + yAxisRange); 
    }
    else {
        FracComp = (TH1F*) MCComp->Clone("ratioComp");
        FracComp->Add(DataComp, -1);
        FracComp->Divide(FracComp, DataComp, 1, 1, "");
        HistAxisAttSet(FracComp->GetYaxis(), TString("(MC-Data)/Data"), .15, .54, .14, .011, -1.0 * yAxisRange, 1.0 * yAxisRange);

    }  
    mcCompHistCentValVec->push_back(WJComp);
    mcCompHistCentValVec->push_back(VVComp);
    mcCompHistCentValVec->push_back(ZDYComp);
    mcCompHistCentValVec->push_back(SingTopComp);
    mcCompHistCentValVec->push_back(TTbarComp);
    cout << "inside Histogram Adder Syst line: 340 " << endl;
    cout << "QCDComp->GetName()  " << QCDComp->GetName()  << endl;
    cout << "QCDComp->Integral()  " << QCDComp->Integral()  << endl;
//    if (QCDComp->Integral() > 0) mcCompHistCentValVec->push_back(QCDComp); put this in later
    mcCompHistCentValVec->push_back(QCDComp);
    cout << "inside Histogram Adder Syst line: 342" << endl;
//    if (StopComp->Integral() > 0) mcCompHistCentValVec->push_back(StopComp);
    cout << "inside Histogram Adder Syst line: 344" << endl;
}

vector<float> * ScaleBackVecCalc(vector<TFile *> * inputFiles) {
    vector<float> * outVec = new vector<float>;
    TString mcplot = "h_nVtx_inclusive";
    TString mcplot_preRW = "h_nVtx_preRW_inclusive";
    TString fileName;
    TH1F * nVtxOrigHist;
    TH1F * nVtxNewHist;
    int NBinsX;
    float scaleBack;
    for (unsigned int i = 0; i < inputFiles->size(); ++i) {
        fileName = inputFiles->at(i)->GetName();
//        if (fileName.Contains("Data")) continue;
        nVtxOrigHist = (TH1F*) inputFiles->at(i)->Get(mcplot_preRW);
        nVtxNewHist = (TH1F*) inputFiles->at(i)->Get(mcplot);
        NBinsX = nVtxOrigHist->GetNbinsX();
        scaleBack = (float) nVtxOrigHist->Integral(1, NBinsX + 1) / nVtxNewHist->Integral(1, NBinsX + 1);
        std::cout << "scaleBack " << scaleBack << std::endl;
        outVec->push_back(scaleBack);
    }
    return outVec;
}


void SpectrumDraw(TCanvas * InputCanvas, TH1F * Hist1, TString legHist1, TH1F * Hist2, TH1F * fracRatioHist, TH1F * errHist, THStack * MCStack, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB, bool logYPad1, vector<TString> * mcLegends, vector<TH1F *> * indMCHists, bool doMean, TString cutUsed, float inputLumi) {
    TLatex * tl = new TLatex();
    TLegend * leg;
    tl->SetTextAlign(12);
    tl->SetNDC();
    tl->SetTextSize(0.03);
    char buf[99];        
    InputCanvas->Divide(1, 2);
    TPad * Pad1 = (TPad *) InputCanvas->cd(1);
    FixPad(Pad1, 1, InputCanvas);    
    Pad1->SetLogy(logYPad1);
    errHist->Draw("e2");
    MCStack->Draw("hist same");
    errHist->Draw("e2 same");
    Hist1->Draw("e1 same");
    Hist1->Draw("axis same");
    int NBins = Hist1->GetNbinsX();
    TAxis * XAxis = errHist->GetXaxis();
    float XBinUB = XAxis->GetXmax();
    float XBinLB = XAxis->GetXmin();
    float BinWidthGeVInit = (XBinUB - XBinLB)/NBins;
    float BinWidthGeV = nDigits(BinWidthGeVInit, 3);
    stringstream ss (stringstream::in | stringstream::out);
    ss << BinWidthGeV;   
    string attempt = ss.str();
    TString BinWidthString;
    BinWidthString += attempt;
    TAxis * YAxis = errHist->GetYaxis();
    TString YAxisTitle = YAxis->GetTitle();    
    int StrPos = YAxisTitle.Index("NUM", 3, TString::kExact);
    //    cout << "StrPos " << StrPos << endl;
    if (StrPos != -1) YAxisTitle.Replace(StrPos, 3, BinWidthString); 
    //    TString YAxisTitle = "Entries / ";
    //    TString YAxisTitle = "Entries / ";
    TString XAxisTitle = XAxis->GetTitle();
    //    YAxisTitle += BinWidthGeV;
    //    YAxisTitle += " GeV";
    HistAxisAttSet(YAxis, YAxisTitle, .06, 1.5, .05, .007, YAxisLB, YAxisUB);
    HistAxisAttSet(XAxis, XAxisTitle, YAxis->GetTitleSize(), 999, YAxis->GetLabelSize(), 999, 0.0, 0.0);
    YAxis->SetNdivisions(5,5,0);
    //    cout << "Y Axis Title " << YAxis->GetTitle() << endl;
    TString DataName = TString(Hist1->GetName());
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
        sprintf(buf,"Mean = %0.2f #pm %0.2f",Hist1->GetMean(), Hist1->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",Hist1->GetRMS(), Hist1->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.05,buf);
        tl->SetTextColor(kRed);
        sprintf(buf,"Mean = %0.2f #pm %0.2f",Hist2->GetMean(),Hist2->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.1,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",Hist2->GetRMS(),Hist2->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.15,buf);    
    }
    //    leg = new TLegend(legXstart,legYstart - 0.05 * (numIndMC + 1),legXstart + 0.05 * (numIndMC + 1),legYstart);
    leg = new TLegend(legXstart,legYstart - 0.05 * (mcLegends->size() + 1),legXstart + 0.40,legYstart);
    leg->AddEntry(Hist1,TString(legHist1),"pl");
    for (unsigned int j = 0; j < mcLegends->size(); ++j) {
        leg->AddEntry(indMCHists->at(j), mcLegends->at(j), "f");
    }
    leg->AddEntry(errHist,"Sim unc","f");
    leg->Draw("same");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(31);
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
    YAxis = fracRatioHist->GetYaxis();
    XAxis = fracRatioHist->GetXaxis();
    //    cout << XAxis->GetTitle() << endl;
//    HistAxisAttSet(XAxis, XAxisTitle, .17, 1.03, .12,.07, 0.0, 0.0);
//    HistAxisAttSet(YAxis, YAxis->GetTitle(),);
//    HistAxisAttSet(XAxis, XAxisTitle, .18, 1.03, .18,.07, 0.0, 0.0);
    HistAxisAttSet(XAxis, XAxisTitle, .17, 1.03, .12,.07, 0.0, 0.0);
    YAxis->SetNdivisions(3,5,0);
    XAxis->SetNdivisions(6,5,0);
    fracRatioHist->SetLineColor(kBlack);
    fracRatioHist->Draw("e1");
    Pad2->Update();
    Pad2->Modified();
}


void SpectrumDrawSyst(TCanvas * InputCanvas, TH1F * Hist1, TString legHist1, TH1F * Hist2, THStack * MCStack, TGraphAsymmErrors * errGraph, TGraphAsymmErrors * errGraphjustSyst, TH1F * fracRatioHist, TGraphAsymmErrors * fracRatioGraph, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB, bool logYPad1, vector<TString> * mcLegends, vector<TH1F *> * indMCHists, bool doMean, TString cutUsed, float inputLumi) {
    TLatex * tl = new TLatex();
    TLegend * leg;
    tl->SetTextAlign(12);
    tl->SetNDC();
    tl->SetTextSize(0.03);
    char buf[99]; 
    InputCanvas->Divide(1, 2);
    TPad * Pad1 = (TPad *) InputCanvas->cd(1);
    FixPad(Pad1, 1, InputCanvas);
    Pad1->SetLogy(logYPad1);
//    Hist2->Draw("same");    
    /*
    cout << "errGraph name " << errGraph->GetName() << endl;    
    for (int i = 0; i < errGraph->GetN(); ++i) {
        cout << "i " << i << endl;
        cout << "upErr " << errGraph->GetErrorXhigh(i) << endl;
        cout << "downErr " << errGraph->GetErrorXlow(i) << endl;
        cout << "upErr " << errGraph->GetErrorYhigh(i) << endl;
        cout << "downErr " << errGraph->GetErrorYlow(i) << endl;
        cout << "upErrJustSyst " << errGraphjustSyst->GetErrorYhigh(i) << endl;
        cout << "downErrJustSyst " << errGraphjustSyst->GetErrorYlow(i) << endl;
        cout << "rawStatErr " << Hist2->GetBinError(i) << endl;
        cout << "quad add up " << TMath::Sqrt(Hist2->GetBinError(i) * Hist2->GetBinError(i) + errGraphjustSyst->GetErrorYhigh(i) * errGraphjustSyst->GetErrorYhigh(i)) << endl;
        cout << "quad add down " << TMath::Sqrt(Hist2->GetBinError(i) * Hist2->GetBinError(i) + errGraphjustSyst->GetErrorYlow(i) * errGraphjustSyst->GetErrorYlow(i)) << endl;
    }
    */
    Hist1->Draw("e1");
    MCStack->Draw("hist same");
    errGraph->Draw("2 same");
    Hist1->Draw("e1 same");
    Hist1->Draw("axis same");
    int NBins = Hist1->GetNbinsX();
    TAxis * XAxis = Hist1->GetXaxis();
    float XBinUB = XAxis->GetXmax();
    float XBinLB = XAxis->GetXmin();
    float BinWidthGeVInit = (XBinUB - XBinLB)/NBins;
    float BinWidthGeV = nDigits(BinWidthGeVInit, 3);
    stringstream ss (stringstream::in | stringstream::out);
    ss << BinWidthGeV;
    string attempt = ss.str();
    TString BinWidthString;
    BinWidthString += attempt;
    TAxis * YAxis = Hist1->GetYaxis();
    TString YAxisTitle = YAxis->GetTitle();    
    int StrPos = YAxisTitle.Index("NUM", 3, TString::kExact);
    if (StrPos != -1) YAxisTitle.Replace(StrPos, 3, BinWidthString); 
    TString XAxisTitle = XAxis->GetTitle();
    cout << "YAxis name " << YAxis->GetName() << endl;
    HistAxisAttSet(YAxis, YAxisTitle, .06, 1.5, .05, .007, YAxisLB, YAxisUB);
    HistAxisAttSet(XAxis, XAxisTitle, YAxis->GetTitleSize(), 999, YAxis->GetLabelSize(), 999, 0.0, 0.0);
    YAxis->SetNdivisions(5,5,0);
    TString DataName = TString(Hist1->GetName());
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
        sprintf(buf,"Mean = %0.2f #pm %0.2f",Hist1->GetMean(), Hist1->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",Hist1->GetRMS(), Hist1->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.05,buf);
        tl->SetTextColor(kRed);
        sprintf(buf,"Mean = %0.2f #pm %0.2f",Hist2->GetMean(),Hist2->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.1,buf);
        sprintf(buf,"RMS = %0.2f #pm %0.2f",Hist2->GetRMS(),Hist2->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.15,buf);    
    }
    leg = new TLegend(legXstart,legYstart - 0.05 * (mcLegends->size() + 1),legXstart + 0.40,legYstart);
    leg->AddEntry(Hist1,TString(legHist1),"pl");
    for (unsigned int j = 0; j < mcLegends->size(); ++j) {
        leg->AddEntry(indMCHists->at(j), mcLegends->at(j), "f");
    }
    leg->AddEntry(errGraph,"Sim unc","f");
    leg->Draw("same");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(31);
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
    YAxis = fracRatioHist->GetYaxis();
    XAxis = fracRatioHist->GetXaxis();
    HistAxisAttSet(XAxis, XAxisTitle, .17, 1.03, .12,.07, 0.0, 0.0);
    YAxis->SetNdivisions(3,5,0);
    XAxis->SetNdivisions(6,5,0);
    fracRatioHist->SetLineColor(kBlack);
    fracRatioHist->Draw("e1");
    fracRatioGraph->Draw("3 same");
    fracRatioHist->Draw("e1 same");
    Pad2->Update();
    Pad2->Modified();
}
