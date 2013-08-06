#include <vector>
#include "TLegend.h"
#include "./StopFunctionDefinitions_v2.h"
//#include <boost/format.hpp>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include "TRegexp.h"
using namespace std;

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
float ScaleBackCalc(TFile * inputFile) {
    TString mcplot = "h_nVtx_inclusive";
    TString mcplot_preRW = "h_nVtx_preRW_inclusive";
    TString fileName;
    TH1F * nVtxOrigHist;
    TH1F * nVtxNewHist;
    int NBinsX;
    float scaleBack;
    fileName = inputFile->GetName();
    nVtxOrigHist = (TH1F*) inputFile->Get(mcplot_preRW);
    nVtxNewHist = (TH1F*) inputFile->Get(mcplot);
    NBinsX = nVtxOrigHist->GetNbinsX();
    scaleBack = (float) nVtxOrigHist->Integral(1, NBinsX + 1) / nVtxNewHist->Integral(1, NBinsX + 1);
    return scaleBack;
}
vector<float> * SignalSkimEfficiencyCalc(TString typeSMS, vector<int> * vecStopMassGrab, vector<int> * vecChi0MassGrab, vector<int> * vecCharginoMassGrab, float intLumi) {
    vector<float> * scaleBackVec = new vector<float>;
    char Buffer[500];
    TString basePath = "../PlotMakingCode/FineBinNormFiles/";
    TString nameNormFileSignalBase = "NormNumbers_StopMass";
    TString nameNormFileSignal;
    ifstream * streamNormFileSignal;
    TRegexp fExtractNumber(" [0-9]*\.[0-9]*");
    TString currNumString;
    TString currString;
    float SkimEff, SkimEffErr, estOrigNum, estErrOrigNum;
    for (unsigned int iSkimEff = 0; iSkimEff < vecStopMassGrab->size(); ++iSkimEff) {
        nameNormFileSignal = basePath; nameNormFileSignal += nameNormFileSignalBase;
        nameNormFileSignal += vecStopMassGrab->at(iSkimEff);
        nameNormFileSignal += TString("_Chi0Mass");
        nameNormFileSignal += vecChi0MassGrab->at(iSkimEff);
        nameNormFileSignal += TString("_CharginoMass");
        nameNormFileSignal += vecCharginoMassGrab->at(iSkimEff);
        nameNormFileSignal += TString(".txt");
        streamNormFileSignal = new ifstream(nameNormFileSignal);
        for (int iLine = 0; iLine < 4; ++iLine) {
            streamNormFileSignal->getline(Buffer, 500);
            cout << "line " << Buffer << endl;
            currString = TString(Buffer);
            currNumString = currString(fExtractNumber);
            switch (iLine) {
                case 0:
                    SkimEff = atof(currNumString);
                    cout << "SkimEff " << SkimEff << endl;
                    break;
                case 1:
                    SkimEffErr = atof(currNumString);
                    cout << "SkimEffErr " << SkimEffErr << endl;
                    break;                
                case 2:
                    estOrigNum = atof(currNumString);
                    cout << "estOrigNum " << estOrigNum << endl;
                    break;
                case 3:
                    estErrOrigNum = atof(currNumString);
                    cout << "estErrOrigNum " << estOrigNum << endl;
                    break;
                    
                default:
                    break;
            }
        }
        cout << "intLumi " << intLumi << endl;
        cout << "estOrigNum " << estOrigNum << endl;
        cout << "intLumi / estOrigNum " << (intLumi / estOrigNum) << endl;
        scaleBackVec->push_back(intLumi / estOrigNum);
        streamNormFileSignal->close();
    }
    return scaleBackVec;
}
float DataDrivenTTBarScaleFactor(TH1F * dataHist, TH1F * mcHist, vector<TH1F *> * mcCompHist1DCentValVec) {
    cout << "test of dataHist name " << dataHist->GetName() << endl;
    float dataIntegral = dataHist->Integral();
    float dataIntegralMinNonTTBar = dataIntegral;
    float mcIntegral_v1 = mcHist->Integral();
    float mcIntegral_v2 = 0;
    float integralTTBar;
    float currMCIntegral;
    TString mcName;
    for (unsigned int iMC = 0; iMC < mcCompHist1DCentValVec->size(); ++iMC) {
        mcName = mcCompHist1DCentValVec->at(iMC)->GetName();
        currMCIntegral = mcCompHist1DCentValVec->at(iMC)->Integral();
        cout << "mcName " << mcName << endl;
        if (mcName.Contains("TTbar") || mcName.Contains("TTBar")) {
            integralTTBar = currMCIntegral;
        }
        else {
            dataIntegralMinNonTTBar -= currMCIntegral;
        }
        mcIntegral_v2 += currMCIntegral;
    }
    if (abs(mcIntegral_v2 - mcIntegral_v1) > 10) {
        cout << "something funky...mc comp integral calculated two different but ostensibly equivalent ways isn't the same..." << endl;
        cout << "mcIntegral_v1 " << mcIntegral_v1 << endl;
        cout << "mcIntegral_v2 " << mcIntegral_v2 << endl;
    }
    cout << "test1 " <<  integralTTBar << endl;
    cout << "test2 " << dataIntegralMinNonTTBar << endl;
    return dataIntegralMinNonTTBar/integralTTBar;
}
vector<TString> * StopFileNames(int whichNTuple, int whichTTBarGen, bool doExcSamp) {
    vector<TString> * outVec = new vector<TString>;
    const int numOviTypes = 18;
    TString fileInNameSpecOviedo[numOviTypes] = {"DataMuMu", "DataEMu", "DataEE", "TTBarComp", "TTBarBkg", "WW", "WZ", "ZZ", "SingleTop", "WToLNu", "ZDY", "WG", "ZG", "HiggsWW", "HiggsVBF", "HiggsZZ4L", "TripVecBoson", "TTBarVecBoson"};
    const int numDESYTypes = 12;
    TString fileInNameSpecDESY[numDESYTypes] = {"DataMuMu", "DataEMu", "DataEE", "TTBarSig", "TTBarBkg", "WW", "WZ", "ZZ", "SingleTop", "WToLNu", "ZDY", "QCD"};
    if (doExcSamp && whichTTBarGen == 0) {
        fileInNameSpecOviedo[3] = TString("TTBarSig");
    }
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
vector<TFile *> * StopFiles(int whichNTuple, vector<TString> * fileNames, bool doExcSamp, int whichTTBarGen, bool doPURW, bool doSyst) {
    vector<TFile *> * outVec = new vector<TFile *>;
    TFile * outTFile;
    TString addPath, fileInNameBase, specNTupString, fileNameSuffix, fileName, TTBarGenString, PURWString, SystString;
    TString doExcSampString = (whichNTuple == 0) ? "_Exclusive" : "";
    if (whichNTuple == 0) {
        addPath        = "../RootFiles/";
        fileInNameBase = "";
        specNTupString = "_Oviedo";
        fileNameSuffix = "Haddplots.root";
        TTBarGenString = "";
        switch (whichTTBarGen) {
            case 0:
                TTBarGenString = "_Madgraph";
                break;
            case 1:
                TTBarGenString = "_MCatNLO";
                break;
            case 2:
                TTBarGenString = "_Powheg";
                break;                
            default:
                break;
        }
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
        if (fileName.Contains("TTBar") && !fileName.Contains("VecBoson")) {
            fileName += TTBarGenString;
            fileName += specNTupString;
            if (whichTTBarGen == 0 && doExcSamp) fileName += doExcSampString;
        }
        else {
            fileName += specNTupString;
        }
        if (fileName.Contains("ZDY") && doExcSamp) fileName += doExcSampString;
        if (!fileName.Contains("Data")) {
            PURWString = doPURW ? "_PURW" : "";
            SystString = doSyst ? "_wSyst" : "";
        }
        //        outTFile = new TFile(addPath + fileInNameBase + fileName + PURWString + SystString + fileNameSuffix);
        outTFile = TFile::Open(addPath + fileInNameBase + fileName + PURWString + SystString + fileNameSuffix);
        outVec->push_back(outTFile);
    }
    return outVec;
}
vector<TFile *> * StopSignalFiles(int whichNTuple, TString typeSMS, vector<int> * vecStopMassGrab, vector<int> * vecChi0MassGrab, vector<int> * vecCharginoMassGrab, bool doPURW, bool doSyst) {
    vector<TFile *> * outVec = new vector<TFile *>;
    bool isT2tt = typeSMS.Contains("T2tt");
    TFile * outTFile;
    TString addPath = "../RootFiles/Signal/";
    TString baseString = isT2tt ? "FineBin" : "T2bw";
    TString specNTupString = (whichNTuple == 1) ? "_DESY" : "_Oviedo";
    TString PURWString = doPURW ? "_PURW" : "";
    TString SystString = doSyst ? "_wSyst" : "";
    TString fileNameSuffix = "_Output.root";
    TString fileName;
    for (unsigned int i = 0; i < vecStopMassGrab->size(); ++i) {
        fileName = addPath + baseString;
        fileName += "_Signal_SignalStop";
        fileName += vecStopMassGrab->at(i);
        fileName += "_Chi0_";
        fileName += vecChi0MassGrab->at(i);
        if (!isT2tt) {
            fileName += "_Chargino_";
            fileName += vecCharginoMassGrab->at(i);
        }
        fileName += specNTupString;
        fileName += PURWString;
        fileName += SystString;
        fileName += fileNameSuffix;
        cout << "signal File name to grab " << fileName << endl;
        outTFile = TFile::Open(fileName);
        outVec->push_back(outTFile);
    }
    return outVec;
}


vector<TString> * MCLegends(int whichNTuple, bool addThings) {
    vector<TString> * outVec = new vector<TString>;
    const int numOviTypesAdd = 8;
    //    const int numOviTypes = 7;
    const int numDESYTypesAdd = 6;
    const int numDESYTypes = 9;
    TString mcLegendsOviAdd[numOviTypesAdd] = {"Rare", "Higgs", "W + Jets", "VG", "VV", "Z/#gamma* #rightarrow l^{+}l^{-}", "Single Top", "t#bar{t}"};
    TString mcLegendsDESYAdd[numDESYTypesAdd] = {"W + Jets", "VV", "Z/#gamma* #rightarrow l^{+}l^{-}", "Single Top", "t#bar{t}", "Multijets"};
    TString mcLegendsDESY[numDESYTypes] = {"W + Jets", "WW", "WZ", "ZZ", "Z/#gamma* #rightarrow l^{+}l^{-}", "Single Top", "t#bar{t} sig", "t#bar{t} bkg", "Multijets"};
    switch (whichNTuple) {
        case 0:
            if (addThings) {
                for (int i = 0; i < numOviTypesAdd; ++i) {
                    outVec->push_back(mcLegendsOviAdd[i]);
                }
            }
            else {
                /*
                 for (int i = 0; i < numOviTypes; ++i) {
                 outVec->push_back(mcLegendsOvi[i]);
                 }
                 */
            }
            break;
        case 1:
            if (addThings) {
                for (int i = 0; i < numDESYTypesAdd; ++i) {
                    outVec->push_back(mcLegendsDESYAdd[i]);
                }
            }
            else {
                for (int i = 0; i < numDESYTypes; ++i) {
                    outVec->push_back(mcLegendsDESY[i]);
                }
            }
            break;            
        default:
            break;
    }
    return outVec;
}

vector<TString> * MCSignalLegends(TString typeSMS, vector<int> * vecStopMassGrab, vector<int> * vecChi0MassGrab, vector<int> * vecCharginoMassGrab) {
    vector<TString> * outVec = new vector<TString>;
//    TString stringModelBase = "Signal, ";
    TString stringModelBase = "";
    TString stringStopMass = "M_{#tilde{t}} = ";
    TString stringChi0Mass = "M_{#chi_{0}} = ";
    TString stringCharginoMass = "M_{#chi^{#pm}} = ";
    TString currLegEntry;
    for (unsigned int i = 0; i < vecStopMassGrab->size(); ++i) {
        currLegEntry = stringModelBase + stringStopMass;
        currLegEntry += vecStopMassGrab->at(i);
        currLegEntry += TString(", ");
        currLegEntry += stringChi0Mass;
        currLegEntry += vecChi0MassGrab->at(i);
        if (typeSMS.Contains("T2bw")) {
            currLegEntry += TString(", ");
            currLegEntry += stringCharginoMass;
            currLegEntry += vecCharginoMassGrab->at(i);
        }
        outVec->push_back(currLegEntry);
    }
    return outVec;
}

vector<Color_t> * MCColors(int whichNTuple, bool addThings) {
    vector<Color_t> * outVec = new vector<Color_t>;
    const int numColorsOviAdd = 8;
    //    const int numColorsOvi = 10;    
    const int numColorsDESYAdd = 6;
    const int numColorsDESY= 9;
    Color_t mcColorsOviAdd[numColorsOviAdd] = {kGreen + 3, kCyan - 2, kGreen + 2, kOrange - 5, kOrange + 2, kBlue, kRed - 10, kRed};
    Color_t mcColorsDESYAdd[numColorsDESYAdd] = {kGreen + 2, kOrange + 2, kBlue, kRed - 10, kRed, kCyan -4};
    Color_t mcColorsDESY[numColorsDESY] = {kGreen + 2, kOrange + 2, kPink + 9, kPink - 8, kBlue, kRed - 10, kRed-5, kRed, kCyan - 4};
    switch (whichNTuple) {
        case 0:
            if (addThings) {
                for (int i = 0; i < numColorsOviAdd; ++i) {
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
                for (int i = 0; i <  numColorsDESYAdd; ++i) {
                    outVec->push_back(mcColorsDESYAdd[i]);
                }
            }
            /*
             else {
             for (int i = 0; i < numColorsDESY; ++i) {
             outVec->push_back(mcColorsDESY[i]);
             }
             }
             */
            break;            
        default:
            break;
    }
    return outVec;
}
vector<Color_t> * MCSignalColors(unsigned int numSignalPoints) {
    vector<Color_t> * outVec = new vector<Color_t>;
    Color_t signalColors[6] = {kMagenta, kMagenta + 1, kMagenta + 2, kViolet + 1, kViolet + 5, kViolet + 10};
    if (numSignalPoints > 6) {
        cout << "more signal points than colors for!!" << endl;
    }
    else {
        for (unsigned int iSigColor = 0; iSigColor < numSignalPoints; ++iSigColor) {
            outVec->push_back(signalColors[iSigColor]);
        }
    }
    return outVec;
}
vector<Style_t> * MCSignalStyles(unsigned int numSignalPoints) {
    vector<Style_t> * outVec = new vector<Color_t>;
    Style_t signalStyles[6] = {1, 2, 4, 8, 6, 10};
    if (numSignalPoints > 6) {
        cout << "more signal points than line styles for!!" << endl;
    }
    else {
        for (unsigned int iSigStyle = 0; iSigStyle < numSignalPoints; ++iSigStyle) {
            outVec->push_back(signalStyles[iSigStyle]);
        }
    }
    return outVec;
}
void SingSampCompHistogramGrabber(TFile * inputFile1, TH1F * &inputHist1, vector<TH1F *> * inputSystVec1, TString hist1Name, TFile * inputFile2, TH1F * &inputHist2, vector<TH1F *> * inputSystVec2,  TString hist2Name, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst) {
    int NBinsX, NBinsY, NBinsZ;
    TString fileName1, fileName2, systName, mcGrabName;
    TH1F * currHist;
    TH1F * systCompHist = NULL;
    mcGrabName = "";
    fileName1 = inputFile1->GetName();
    inputHist1 = (TH1F *) inputFile1->Get(hist1Name);
    inputHist1->RebinX(RBNX);
    float nVtxBackScale1 = ScaleBackCalc(inputFile1);
    inputHist1->Scale(nVtxBackScale1);
    if (hist2Name != "") {
        fileName2 = inputFile2->GetName();
        inputHist2 = (TH1F *) inputFile2->Get(hist2Name);
        inputHist2->RebinX(RBNX);
        float nVtxBackScale2 = ScaleBackCalc(inputFile2);
        inputHist2->Scale(nVtxBackScale2);
    }
}
void HistogramVecGrabber(vector<TFile *> * inputFiles, vector<TH1F *> * dataHistVec, vector<TH1F *> * mcIndHistCentValVec, vector<TH1F *> * mcCompHistSystVec, vector<float> * nVtxBackScaleVec, vector<SystT> * systVec, TString dataPlotName, TString mcPlotName, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst, bool useDDEstimate, float scaleFacTTBar, bool allMT2llSystematic, int whichNTuple) {
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
            //            cout << "test " << endl;
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
            fileName = inputFiles->at(k)->GetName();
            if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);
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
            if (systName.Contains("genStop")) continue;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                if (systName.Contains("genTop")) {
                    if (!fileName.Contains("TTBar")) continue;
                    if (whichNTuple == 1 && !fileName.Contains("Sig")) continue;
                }
                if (systName.Contains("genStop") && !fileName.Contains("Stop")) continue;
                
                if (!allMT2llSystematic && systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH1F *) inputFiles->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH1F *) inputFiles->at(k)->Get(mcSystPlotName);   
                }
                
                //currHist = (TH1F *) inputFiles->at(k)->Get(mcSystPlotName);   
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
                if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);
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


void HistogramVecGrabber_Signal(vector<TFile *> * inputFilesSignal, vector<float> * sigSkimEffScaleCalcVec, vector<TH1F *> * vecStop1DCentValHists, TString mcPlotName, vector<vector<TH1F *> *> * vecStop1DSystHists, vector<SystT> * systVec, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst, bool allMT2llSystematic) {
    int NBinsX, NBinsY, NBinsZ;
    TString fileName, systName, mcGrabName;
    TH1F * currHist;
    TH1F * systCompHist = NULL;
    cout << " trying to grab signal scale " << endl;
    vector <float> * nVtxSignalBackScaleVec = ScaleBackVecCalc(inputFilesSignal);
    cout << " trying to grab signal " << endl;
    TString mcSystPlotName, systCompHistName;
    for (unsigned int k = 0; k < inputFilesSignal->size(); ++k) {
        mcGrabName = "";
        fileName = inputFilesSignal->at(k)->GetName();
        cout << "fileName" << " is " << fileName << endl;
        mcGrabName += mcPlotName;
        mcGrabName += subSampName;
        currHist = (TH1F *) inputFilesSignal->at(k)->Get(mcGrabName);
        currHist->Scale(nVtxSignalBackScaleVec->at(k)); // correct for PURW changes to integral            
        currHist->Scale(sigSkimEffScaleCalcVec->at(k)); // correct for PURW changes to integral            
        NBinsX = currHist->GetNbinsX();
        NBinsY = currHist->GetNbinsY();
        NBinsZ = currHist->GetNbinsZ();
        currHist->RebinX(RBNX);
        vecStop1DCentValHists->push_back(currHist);
    }
    if (doSyst) {
        for (unsigned int k = 0; k < inputFilesSignal->size(); ++k) {
            fileName = inputFilesSignal->at(k)->GetName();
            vector<TH1F *> * currStopFileSystVec = new vector<TH1F *>;
            systCompHist = NULL;
            cout << "k " << k << endl;            
            for (unsigned int j = 0; j < systVec->size(); ++j) {
                cout << "j " << j << endl;
                systName = systVec->at(j).name;
                mcGrabName = "";
                mcGrabName += mcPlotName;
                cout << "systName " << systName << endl;
                cout << "mcGrabName " << mcGrabName << endl;
                if (systName.Contains("genTop")) continue;
                if (!mcGrabName.Contains("MT2ll") && systName.Contains("MT2ll")) continue;
                mcSystPlotName = mcGrabName;
                mcSystPlotName += systName;
                mcGrabName += subSampName;
                mcSystPlotName += subSampName;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (!allMT2llSystematic && systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH1F *) inputFilesSignal->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH1F *) inputFilesSignal->at(k)->Get(mcSystPlotName);   
                }   
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxSignalBackScaleVec->at(k)); // correct for PURW changes to integral
                currHist->Scale(sigSkimEffScaleCalcVec->at(k)); // correct for PURW changes to integral
                NBinsX = currHist->GetNbinsX();
                NBinsY = currHist->GetNbinsY();
                NBinsZ = currHist->GetNbinsZ();
                currHist->RebinX(RBNX);
                currStopFileSystVec->push_back(currHist);
                cout << "rebinned hist " << endl;
            }
            vecStop1DSystHists->push_back(currStopFileSystVec);
        }
    }
}




void HistogramVecGrabber(vector<TFile *> * inputFiles, vector<TH2F *> * dataHistVec, vector<TH2F *> * mcIndHistCentValVec, vector<TH2F *> * mcCompHistSystVec, vector<float> * nVtxBackScaleVec, vector<SystT> * systVec, TString dataPlotName, TString mcPlotName, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst, bool useDDEstimate, float scaleFacTTBar, bool allMT2llSystematic, int whichNTuple) {
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
            fileName = inputFiles->at(k)->GetName();
            currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
            if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);
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
            if (systName.Contains("genStop")) continue;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                if (systName.Contains("genTop")) {
                    if (!fileName.Contains("TTBar")) continue;
                    if (whichNTuple == 1 && !fileName.Contains("Sig")) continue;
                }
                if (systName.Contains("genStop") && !fileName.Contains("Stop")) continue;                
                if (!allMT2llSystematic && systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH2F *) inputFiles->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH2F *) inputFiles->at(k)->Get(mcSystPlotName);   
                }                
                //currHist = (TH2F *) inputFiles->at(k)->Get(mcSystPlotName);                   
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
                if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);
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
void HistogramVecGrabber(vector<TFile *> * inputFiles, vector<TH3F *> * dataHistVec, vector<TH3F *> * mcIndHistCentValVec, vector<TH3F *> * mcCompHistSystVec, vector<float> * nVtxBackScaleVec, vector<SystT> * systVec, TString dataPlotName, TString mcPlotName, TString subSampName, int RBNX, int RBNY, int RBNZ, bool doOverflow, bool doUnderflow, bool doSyst, bool useDDEstimate, float scaleFacTTBar, bool allMT2llSystematic, int whichNTuple) {
    int NBinsX, NBinsY, NBinsZ;
    TString fileName, systName, mcGrabName;
    TH3F * currHist;
    TH3F * systCompHist = NULL;
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
            fileName = inputFiles->at(k)->GetName();
            currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
            if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);
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
            if (systName.Contains("genStop")) continue;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                if (systName.Contains("genTop")) {
                    if (!fileName.Contains("TTBar")) continue;
                    if (whichNTuple == 1 && !fileName.Contains("Sig")) continue;
                }                
                if (systName.Contains("genStop") && !fileName.Contains("Stop")) continue;            
                if (systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH3F *) inputFiles->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH3F *) inputFiles->at(k)->Get(mcSystPlotName);   
                }        
                //                currHist = (TH3F *) inputFiles->at(k)->Get(mcSystPlotName);   
                cout << "grabbed hist " << endl;
                currHist->Scale(nVtxBackScaleVec->at(k)); // correct for PURW changes to integral
                if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);
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
void HistogramAdderSyst(vector<TH1F *> * dataHistVec, vector<TH1F *> * mcIndHistCentValVec, vector<TH1F *> * mcCompHistCentValVec, TH1F * &DataComp, TH1F * &MCComp, TH1F * &FracComp, int whichNTuple, bool doAbsRatio, float yAxisRange, bool doExcSamp) {
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
    TH1F * VGComp, * HiggsComp, * RareComp;
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
    if (doExcSamp) {
        TTbarComp->Add(mcIndHistCentValVec->at(1));
    }
    VVComp = (TH1F *) mcIndHistCentValVec->at(2)->Clone(mcIndHistCentValVec->at(2)->GetName() + TString("_VV"));
    VVComp->Add(mcIndHistCentValVec->at(3));
    VVComp->Add(mcIndHistCentValVec->at(4));
    SingTopComp = (TH1F *) mcIndHistCentValVec->at(5)->Clone(mcIndHistCentValVec->at(5)->GetName() + TString("_SingTop"));
    WJComp = (TH1F *) mcIndHistCentValVec->at(6)->Clone(mcIndHistCentValVec->at(6)->GetName() + TString("_WJ"));
    ZDYComp = (TH1F *) mcIndHistCentValVec->at(7)->Clone(mcIndHistCentValVec->at(7)->GetName() + TString("_ZDY"));
    if (whichNTuple == 1) {
        QCDComp = (TH1F *) mcIndHistCentValVec->at(8)->Clone(mcIndHistCentValVec->at(8)->GetName() + TString("_QCD"));
    }
    else {
        VGComp = (TH1F *) mcIndHistCentValVec->at(8)->Clone(mcIndHistCentValVec->at(8)->GetName() + TString("_VG"));
        VGComp->Add(mcIndHistCentValVec->at(9));
        HiggsComp = (TH1F *) mcIndHistCentValVec->at(10)->Clone(mcIndHistCentValVec->at(10)->GetName() + TString("_Higgs"));
        HiggsComp->Add(mcIndHistCentValVec->at(11));
        HiggsComp->Add(mcIndHistCentValVec->at(12));
        RareComp = (TH1F *) mcIndHistCentValVec->at(13)->Clone(mcIndHistCentValVec->at(10)->GetName() + TString("_Rare"));
        RareComp->Add(mcIndHistCentValVec->at(14));
    }
    for (unsigned int j = 0; j < mcIndHistCentValVec->size(); ++j) {
        //        MCName = mcIndHistCentValVec->at(j)->GetName();
        if (j != 0) {
            if (j == 1 && !doExcSamp) continue;
            MCComp->Add(mcIndHistCentValVec->at(j));
        }
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
    
    cout << "inside Histogram Adder Syst line: 322" << endl;
    if (whichNTuple == 0) {
        mcCompHistCentValVec->push_back(RareComp);
        mcCompHistCentValVec->push_back(HiggsComp);
        mcCompHistCentValVec->push_back(WJComp);
        mcCompHistCentValVec->push_back(VGComp);
        mcCompHistCentValVec->push_back(VVComp);
        mcCompHistCentValVec->push_back(ZDYComp);
        mcCompHistCentValVec->push_back(SingTopComp);
        mcCompHistCentValVec->push_back(TTbarComp);
    }
    else {
        mcCompHistCentValVec->push_back(WJComp);
        mcCompHistCentValVec->push_back(VVComp);
        mcCompHistCentValVec->push_back(ZDYComp);
        mcCompHistCentValVec->push_back(SingTopComp);
        mcCompHistCentValVec->push_back(TTbarComp);  
        mcCompHistCentValVec->push_back(QCDComp);
        
        cout << "QCDComp->GetName()  " << QCDComp->GetName()  << endl;
        cout << "QCDComp->Integral()  " << QCDComp->Integral()  << endl;
    }
    //    if (QCDComp->Integral() > 0) mcCompHistCentValVec->push_back(QCDComp); put this in later
    cout << "inside Histogram Adder Syst line: 342" << endl;
    //    if (StopComp->Integral() > 0) mcCompHistCentValVec->push_back(StopComp);
    cout << "inside Histogram Adder Syst line: 344" << endl;
}
void SpectrumDraw(TCanvas * InputCanvas, TH1F * Hist1, TString legHist1, TH1F * Hist2, TH1F * fracRatioHist, TH1F * errHist, THStack * MCStack, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB, bool logYPad1, vector<TString> * mcLegends, vector<TH1F *> * indMCHists, bool doMean, TString cutUsed, float inputLumi, TLegend * &leg) {
    TLatex * tl = new TLatex();
    //    TLegend * leg;
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
    leg->AddEntry(errHist,"Stat. Unc.","f");
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
    TGraphAsymmErrors * fracRatioDrawGraph = clonePoints(fracRatioHist, false);
    HistToGraphCopyAttributes(fracRatioHist, fracRatioDrawGraph);
    TH1F * patsy = (TH1F*) fracRatioHist->Clone("frac patsy");
    patsy->SetLineColor(kWhite);
    patsy->SetMarkerColor(kWhite);
    fracRatioHist->SetLineColor(kBlack);
    patsy->Draw("e1");    
    fracRatioDrawGraph->Draw("p0 same");
    //    fracRatioHist->Draw("e1");
    Pad2->Update();
    Pad2->Modified();
}



void SpectrumDraw_AddSignal(TCanvas * InputCanvas, vector<TH1F *> * vecSignalHists, vector<TString> * vecSignalLegends, TLegend * leg) {
    TPad * Pad1 = (TPad *) InputCanvas->cd(1);
    for (unsigned int iSigPoints = 0; iSigPoints < vecSignalHists->size(); ++iSigPoints) {
        vecSignalHists->at(iSigPoints)->Draw("hist same");
        leg->AddEntry(vecSignalHists->at(iSigPoints), vecSignalLegends->at(iSigPoints), "pl");
        leg->Draw("same");
    }
    Pad1->Update();
    Pad1->Modified();
    /*
     TPad * Pad2 = (TPad *) InputCanvas->cd(2);    
     TGraphAsymmErrors * fracRatioDrawGraph = clonePoints(fracRatioHist, false);
     HistToGraphCopyAttributes(fracRatioHist, fracRatioDrawGraph);
     TH1F * patsy = (TH1F*) fracRatioHist->Clone("frac patsy");
     patsy->SetLineColor(kWhite);
     patsy->SetMarkerColor(kWhite);
     fracRatioHist->SetLineColor(kBlack);
     patsy->Draw("e1");    
     fracRatioDrawGraph->Draw("p0 same");
     //    fracRatioHist->Draw("e1");
     Pad2->Update();
     Pad2->Modified();
     */
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
    leg = new TLegend(legXstart - 0.05,legYstart - 0.05 * (mcLegends->size() + 1),legXstart + 0.40,legYstart);
    leg->AddEntry(Hist1,TString(legHist1),"pl");
    for (unsigned int j = 0; j < mcLegends->size(); ++j) {
        leg->AddEntry(indMCHists->at(j), mcLegends->at(j), "f");
    }
    leg->AddEntry(errGraph,"Stat #oplus syst","f");
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
    TGraphAsymmErrors * fracRatioDrawGraph = clonePoints(fracRatioHist, false);
    HistToGraphCopyAttributes(fracRatioHist, fracRatioDrawGraph);
    TH1F * patsy = (TH1F*) fracRatioHist->Clone("frac patsy");
    patsy->SetLineColor(kWhite);
    patsy->SetMarkerColor(kWhite);
    fracRatioHist->SetLineColor(kBlack);
    patsy->Draw("e1");    
    fracRatioGraph->Draw("2 same");
    //    fracRatioHist->Draw("e1 same");
    fracRatioDrawGraph->Draw("p0 same");
    Pad2->Update();
    Pad2->Modified();
}


void SpectrumDrawSingSampCompare(TCanvas * InputCanvas, TH1F * Hist1, TString legHist1, TH1F * Hist2, TString legHist2, TH1F * fracRatioHist, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB, bool logYPad1, bool doMean, TString cutUsed, float inputLumi) {
    TLatex * tl = new TLatex();
    TLegend * leg;
    tl->SetTextAlign(12);
    tl->SetNDC();
    tl->SetTextSize(0.03);
    char buf[99];        
    TPad * Pad1, * Pad2;
    Color_t hist1Color, hist2Color;
    bool doTwoHists = !(fracRatioHist == NULL || Hist2 == NULL);
    if (!doTwoHists) {
        Pad1 = (TPad *) InputCanvas->cd(1);
        FixPadSingle(Pad1, InputCanvas);
    }
    else {
        InputCanvas->Divide(1, 2);
        Pad1 = (TPad *) InputCanvas->cd(1);
        FixPad(Pad1, 1, InputCanvas);    
    }
    Pad1->SetLogy(logYPad1);
    Hist1->Draw("e1");
    hist1Color = Hist1->GetMarkerColor();
    if (doTwoHists) {
        Hist2->Draw("hist same");
        Hist1->Draw("e1 same");
        Hist1->Draw("axis same");
        hist2Color= Hist2->GetLineColor();
    }
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
    HistAxisAttSet(YAxis, YAxisTitle, .06, 1.5, .05, .007, YAxisLB, YAxisUB);
    HistAxisAttSet(XAxis, XAxisTitle, YAxis->GetTitleSize(), 999, YAxis->GetLabelSize(), 999, 0.0, 0.0);
    YAxis->SetNdivisions(5,5,0);
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
        tl->SetTextColor(hist1Color);
        sprintf(buf, legHist1);
        sprintf(buf + strlen(buf), " Mean = %0.2f #pm %0.2f", Hist1->GetMean(), Hist1->GetMeanError());
        tl->DrawLatex(MeanXstart,MeanYstart,buf);        
        sprintf(buf, legHist1);
        sprintf(buf + strlen(buf)," RMS = %0.2f #pm %0.2f", Hist1->GetRMS(), Hist1->GetRMSError());
        tl->DrawLatex(MeanXstart,MeanYstart-0.05,buf);
        if (doTwoHists) {
            tl->SetTextColor(hist2Color);
            sprintf(buf, legHist2);
            sprintf(buf + strlen(buf)," Mean = %0.2f #pm %0.2f", Hist2->GetMean(),Hist2->GetMeanError());
            tl->DrawLatex(MeanXstart,MeanYstart-0.1,buf);
            sprintf(buf, legHist2);
            sprintf(buf + strlen(buf)," RMS = %0.2f #pm %0.2f", Hist2->GetRMS(),Hist2->GetRMSError());
            tl->DrawLatex(MeanXstart,MeanYstart-0.15,buf);    
        }
    }
    //    leg = new TLegend(legXstart,legYstart - 0.05 * (numIndMC + 1),legXstart + 0.05 * (numIndMC + 1),legYstart);
    leg = new TLegend(legXstart,legYstart - 0.15,legXstart + 0.40,legYstart);
    leg->AddEntry(Hist1,TString(legHist1),"pl");
    if (doTwoHists) {
        leg->AddEntry(Hist2,TString(legHist2),"f");
    }
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
    if (doTwoHists) {
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
        TGraphAsymmErrors * fracRatioDrawGraph = clonePoints(fracRatioHist, false);
        HistToGraphCopyAttributes(fracRatioHist, fracRatioDrawGraph);
        TH1F * patsy = (TH1F*) fracRatioHist->Clone("frac patsy");
        patsy->SetLineColor(kWhite);
        patsy->SetMarkerColor(kWhite);
        fracRatioHist->SetLineColor(kBlack);
        //        fracRatioHist->Draw("e1");
        patsy->Draw("e1");
        fracRatioDrawGraph->Draw("p0 same");
        Pad2->Update();
        Pad2->Modified();
    }
}


