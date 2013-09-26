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

void sampleStartPositionsNames(int whichNTuple, int whichTTBarGen, vector<TString> * sampleAddNames, vector<int> * sampleStartPositions, bool doExcSamp) {
    TString TTbarName       = "_TTBar";
    TString VVName          = "_VV";
    TString SingTopName     = "_SingTop";
    TString WJName          = "_WJ";
    TString ZDYName         = "_ZDY";
    TString QCDName         = "_QCD";
    TString VGName          = "_VG";
    TString HiggsName       = "_Higgs";
    TString RareName        = "_Rare";

    int numTTBar = (doExcSamp && whichTTBarGen == 0) ? 2 : 1;
    if (whichNTuple == 1) numTTBar = 2;
    int numVV               = 3;
    int numSingTop          = 1;
    int numWJ               = 1;
    int numZDY              = 1;
    int numQCD              = 1;
    int numVG               = 2;
    int numHiggs            = 3;
    int numRare             = 2;
    
    int TTbarStart          = 0;
    int VVStart             = TTbarStart + numTTBar;
    int SingTopStart        = VVStart + numVV;
    int WJStart             = SingTopStart + numSingTop;
    int ZDYStart            = WJStart + numWJ;
    int QCDStart            = ZDYStart + numZDY;
    int VGStart             = ZDYStart + numZDY;
    int HiggsStart          = VGStart + numVG;
    int RareStart           = HiggsStart + numHiggs;
    sampleAddNames->push_back(TTbarName);    
    sampleAddNames->push_back(VVName);
    sampleAddNames->push_back(SingTopName);
    sampleAddNames->push_back(WJName);
    sampleAddNames->push_back(ZDYName);
    sampleStartPositions->push_back(TTbarStart);    
    sampleStartPositions->push_back(VVStart);
    sampleStartPositions->push_back(SingTopStart);
    sampleStartPositions->push_back(WJStart);
    sampleStartPositions->push_back(ZDYStart);
    if (whichNTuple == 0) {
        sampleAddNames->push_back(VGName);
        sampleAddNames->push_back(HiggsName);
        sampleAddNames->push_back(RareName);
        sampleStartPositions->push_back(VGStart);
        sampleStartPositions->push_back(HiggsStart);
        sampleStartPositions->push_back(RareStart);
        sampleStartPositions->push_back(RareStart + numRare);
    }
    else {
        sampleAddNames->push_back(QCDName);
        sampleStartPositions->push_back(QCDStart);
        sampleStartPositions->push_back(QCDStart + numQCD);
    }    

    sampleAddNames->push_back("patsy");
}
vector<TH1F *> * sortedVector(vector<TH1F *> * inputVector, vector<int> * whichOrder) {
    int currCounter = 0;
    vector<TH1F *> * outVec = new vector<TH1F *>;
    for (unsigned int iBaseCount = 0; iBaseCount < inputVector->size(); ++iBaseCount) {
        currCounter = iBaseCount;
        for (unsigned int iCount = 0; iCount < inputVector->size(); ++iCount) {
            if (whichOrder->at(iCount) == currCounter) outVec->push_back(inputVector->at(iCount));
        }
    }
    return outVec;
}
vector<int> * sortOrder(int whichNTuple) {
    vector<int> * sortVec = new vector<int>;
    const int numOviTypes = 15;
    TString fileInNameSpecOviedo[numOviTypes] = {"TTBarComp", "TTBarBkg", "WW", "WZ", "ZZ", "SingleTop", "WToLNu", "ZDY", "WG", "ZG", "HiggsWW", "HiggsVBF", "HiggsZZ4L", "TripVecBoson", "TTBarVecBoson"};
    const int numDESYTypes = 9;
    TString fileInNameSpecDESY[numDESYTypes] = {"TTBarSig", "TTBarBkg", "WW", "WZ", "ZZ", "SingleTop", "WToLNu", "ZDY", "QCD"};
    if (whichNTuple == 0) {
        sortVec->push_back(7);
        sortVec->push_back(4);
        sortVec->push_back(6);
        sortVec->push_back(2);
        sortVec->push_back(5);
        sortVec->push_back(3);
        sortVec->push_back(1);
        sortVec->push_back(0);        
        /*
        mcCompHistCentValVec->push_back(RareComp);
        mcCompHistCentValVec->push_back(HiggsComp);
        mcCompHistCentValVec->push_back(WJComp);
        mcCompHistCentValVec->push_back(VGComp);
        mcCompHistCentValVec->push_back(VVComp);
        mcCompHistCentValVec->push_back(ZDYComp);
        mcCompHistCentValVec->push_back(SingTopComp);
        mcCompHistCentValVec->push_back(TTbarComp);
         */
    }
    else {
        sortVec->push_back(4);
        sortVec->push_back(1);
        sortVec->push_back(3);
        sortVec->push_back(0);
        sortVec->push_back(2);
        sortVec->push_back(5);
        /*
        mcCompHistCentValVec->push_back(WJComp);
        mcCompHistCentValVec->push_back(VVComp);
        mcCompHistCentValVec->push_back(ZDYComp);
        mcCompHistCentValVec->push_back(SingTopComp);
        mcCompHistCentValVec->push_back(TTbarComp);  
        mcCompHistCentValVec->push_back(QCDComp);
        */
    }
    return sortVec;
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
vector<float> * SignalSkimEfficiencyCalc(TString typeSMS, TString prefixT2tt, vector<int> * vecStopMassGrab, vector<int> * vecChi0MassGrab, vector<int> * vecCharginoMassGrab, float intLumi) {
    vector<float> * scaleBackVec = new vector<float>;
    char Buffer[500];
    
    TString basePath, nameNormFileSignalBase;
    int specT2ttIndex = -1;
    int whichT2tt = 0;
    int stopMassRef[5] = {250, 400, 500, 600, 750};
    int numEvents[5] = {127789, 126051, 122413, 120570, 118333};
    if (prefixT2tt.Contains("FineBin")) {
        whichT2tt = 0;
        basePath = "../PlotMakingCode/T2ttFineBinNormFiles/";
        nameNormFileSignalBase = "NormNumbers_StopMass";
    }
    else if (prefixT2tt.Contains("225to350")) {
        whichT2tt = 2;
        basePath = "../PlotMakingCode/NonFineBinNonSpecificNormFiles/";
        nameNormFileSignalBase = "NormNumbers_T2tt_225to350LSP25to250_LeptonFilter_StopMass";
    }
    else {
        whichT2tt = 1;
    }

    TString nameNormFileSignal;
    ifstream * streamNormFileSignal;
    TRegexp fExtractNumber(" [0-9]*\.[0-9]*");
    TString currNumString;
    TString currString;
    float SkimEff, SkimEffErr, estOrigNum, estErrOrigNum;
    for (unsigned int iSkimEff = 0; iSkimEff < vecStopMassGrab->size(); ++iSkimEff) {
        if (whichT2tt != 1) {
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
            streamNormFileSignal->close();
        }
        else {
            for (int i = 0; i < 5; ++i) {
                if (abs(vecStopMassGrab->at(iSkimEff) - stopMassRef[i]) < 1) specT2ttIndex = i;                
            }
            if (specT2ttIndex < 0) cout << "something bad with specT2ttIndex" << endl;
            estOrigNum = numEvents[specT2ttIndex];
        }
        cout << "intLumi " << intLumi << endl;
        cout << "estOrigNum " << estOrigNum << endl;
        cout << "intLumi / estOrigNum " << (intLumi / estOrigNum) << endl;
        scaleBackVec->push_back(intLumi / estOrigNum);
    }
    return scaleBackVec;
}
float DataDrivenTTBarScaleFactor(TH1F * dataHist, TH1F * mcHist, vector<TH1F *> * mcCompHist1DCentValVec) {
    cout << "test of dataHist name " << dataHist->GetName() << endl;
    float dataIntegral = dataHist->Integral();
    float dataIntegralMinNonTTBar = dataIntegral;
    float mcIntegral_v1 = mcHist->Integral();
    float mcIntegral_v2 = 0;
    float integralTTBar = 0;
    float currMCIntegral;
    TString mcName;
    for (unsigned int iMC = 0; iMC < mcCompHist1DCentValVec->size(); ++iMC) {
        mcName = mcCompHist1DCentValVec->at(iMC)->GetName();
        currMCIntegral = mcCompHist1DCentValVec->at(iMC)->Integral();
        cout << "mcName " << mcName << endl;
        if (mcName.Contains("TTbar") || mcName.Contains("TTBar")) {
            integralTTBar += currMCIntegral;
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

void HistogramUnderflowOverflow(TH1F * inputHist, bool doUnderflow, bool doOverflow) {
    float newHistErr;
    int NBins = inputHist->GetNbinsX();
    if (doUnderflow) {
        inputHist->SetBinContent(1, inputHist->GetBinContent(1) + inputHist->GetBinContent(0));           
        newHistErr = TMath::Sqrt(inputHist->GetBinError(1)*inputHist->GetBinError(1) + inputHist->GetBinError(0)*inputHist->GetBinError(0));
        inputHist->SetBinError(1, newHistErr);
    }
    if (doOverflow) {
        inputHist->SetBinContent(NBins, inputHist->GetBinContent(NBins) + inputHist->GetBinContent(NBins+1)); 
        newHistErr = TMath::Sqrt(inputHist->GetBinError(NBins)*inputHist->GetBinError(NBins) + inputHist->GetBinError(NBins+1)*inputHist->GetBinError(NBins+1));
        inputHist->SetBinError(NBins, newHistErr);
    }
}


vector<TString> * StopFileNames(int whichNTuple, int whichTTBarGen, bool doExcSamp) {
    vector<TString> * outVec = new vector<TString>;
    const int numOviTypes = 18;
    cout << "test " << endl;
    TString fileInNameSpecOviedo[numOviTypes] = {"DataMuMu", "DataEMu", "DataEE", "TTBarComp", "TTBarBkg", "WW", "WZ", "ZZ", "SingleTop", "WToLNu", "ZDY", "WG", "ZG", "HiggsWW", "HiggsVBF", "HiggsZZ4L", "TripVecBoson", "TTBarVecBoson"};
    const int numDESYTypes = 12;
    TString fileInNameSpecDESY[numDESYTypes] = {"DataMuMu", "DataEMu", "DataEE", "TTBarSig", "TTBarBkg", "WW", "WZ", "ZZ", "SingleTop", "WToLNu", "ZDY", "QCD"};
    if (doExcSamp && whichTTBarGen == 0) {
        fileInNameSpecOviedo[3] = TString("TTBarSig");
    }
    switch (whichNTuple) {
        case 0:
            for (int i = 0; i < numOviTypes; ++i) {
                cout << "file name to nominally push back " << fileInNameSpecOviedo[i] << endl;
                if ((!doExcSamp || whichTTBarGen != 0) && fileInNameSpecOviedo[i].Contains("TTBarBkg")) continue;
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
vector<TFile *> * StopFiles(int whichNTuple, vector<TString> * fileNames, bool doExcSamp, int whichTTBarGen, bool doPURW, bool doSyst, int versNumber) {
    vector<TFile *> * outVec = new vector<TFile *>;
    TFile * outTFile;
    TString addPath, fileInNameBase, specNTupString, fileNameSuffix, fileName, TTBarGenString, PURWString, SystString;
    TString doExcSampString = (whichNTuple == 0) ? "_Exclusive" : "";
    if (whichNTuple == 0) {
        addPath        = "../RootFiles/";
        fileInNameBase = "";
        specNTupString = "_Oviedo";
        if (versNumber == 2) {
            fileNameSuffix = "_subLepPtCut20Haddplots.root";
        }
        else {
            fileNameSuffix = "Haddplots.root";   
        }
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
        specNTupString = versNumber == 2 ? "_DESY_SkimOutput_DESY" : "_DESY";
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
vector<TFile *> * StopSignalFiles(int whichNTuple, TString typeSMS, TString prefixT2tt, vector<int> * vecStopMassGrab, vector<int> * vecChi0MassGrab, vector<int> * vecCharginoMassGrab, bool doPURW, bool doSyst) {
    vector<TFile *> * outVec = new vector<TFile *>;
    bool isT2tt = typeSMS.Contains("T2tt");
    TFile * outTFile;
    TString addPath = "../RootFiles/Signal/";
    TString baseString = isT2tt ? prefixT2tt : "T2bw";
    TString specNTupString = (whichNTuple == 1) ? "_DESY" : "_Oviedo";
    TString PURWString = doPURW ? "_PURW" : "";
    TString SystString = doSyst ? "_wSyst" : "";
    TString fileNameSuffix = "_Output.root";
    TString fileName;
    for (unsigned int i = 0; i < vecStopMassGrab->size(); ++i) {
        fileName = addPath + baseString;
        fileName += "_SignalStop_";
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
vector<TString> * IsoPlotNames(int whichNTuple) {
    vector<TString> * outVec = new vector<TString>;
    outVec->push_back("h_CutFlow");
    outVec->push_back("h_ElecRelIso");
    if (whichNTuple == 0) {
        outVec->push_back("h_ElecCharIso");
        outVec->push_back("h_ElecNeutIso");
        outVec->push_back("h_ElecPhotIso");
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
    TString stringStopMass = "M_{#tilde{t}}: ";
    TString stringChi0Mass = "M_{#chi_{0}}: ";
    TString stringCharginoMass = "M_{#chi^{#pm}}: ";
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
    TString fileName1, fileName2, systName, mcGrabName;
//    TH1F * currHist;
//    TH1F * systCompHist = NULL;
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
            if (!mcGrabName.Contains("MT2ll") && (systName.Contains("MT2ll") || systName.Contains("MT2UncES"))) continue;
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
/*

void HistogramVecGrabber_Signal(vector<TFile *> * inputFilesSignal, vector<float> * sigSkimEffScaleCalcVec, vector<TH1 *> * vecStop1DCentValHists, TString mcPlotName, vector<vector<TH1 *> *> * vecStop1DSystHists, vector<SystT> * systVec, TString subSampName, bool doOverflow, bool doUnderflow, bool doSyst, bool allMT2llSystematic) {
    TString fileName, systName, mcGrabName;
    TH1 * currHist;
    TH1 * systCompHist = NULL;
    vector <float> * nVtxSignalBackScaleVec = ScaleBackVecCalc(inputFilesSignal);
    TString mcSystPlotName, systCompHistName;
    bool grabSyst;
    for (unsigned int k = 0; k < inputFilesSignal->size(); ++k) {
        mcGrabName = "";
        fileName = inputFilesSignal->at(k)->GetName();
        cout << "fileName" << " is " << fileName << endl;
        mcGrabName += mcPlotName;
        mcGrabName += subSampName;
        cout << "nVtx scale factor " << nVtxSignalBackScaleVec->at(k) << endl;
        cout << "SkimEfficiency scale factor " << sigSkimEffScaleCalcVec->at(k) << endl;
        currHist = (TH1 *) inputFilesSignal->at(k)->Get(mcGrabName);
        currHist->Scale(nVtxSignalBackScaleVec->at(k)); // correct for PURW changes to integral            
        currHist->Scale(sigSkimEffScaleCalcVec->at(k)); // correct for PURW changes to integral            
        vecStop1DCentValHists->push_back(currHist);
    }
    if (doSyst) {
        for (unsigned int k = 0; k < inputFilesSignal->size(); ++k) {

            fileName = inputFilesSignal->at(k)->GetName();
            vector<TH1 *> * currStopFileSystVec = new vector<TH1 *>;
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
                if (!mcGrabName.Contains("MT2ll") && (systName.Contains("MT2ll") || systName.Contains("MT2UncES"))) continue;
                mcSystPlotName = mcGrabName;
                mcSystPlotName += systName;
                mcGrabName += subSampName;
                mcSystPlotName += subSampName;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (!allMT2llSystematic && systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    currHist = (TH1 *) inputFilesSignal->at(k)->Get(mcGrabName);   
                }
                else {
                    currHist = (TH1 *) inputFilesSignal->at(k)->Get(mcSystPlotName);   
                }   
                currHist->Scale(nVtxSignalBackScaleVec->at(k)); // correct for PURW changes to integral
                currHist->Scale(sigSkimEffScaleCalcVec->at(k)); // correct for skim efficiencies to integral
                currStopFileSystVec->push_back(currHist);
            }
            vecStop1DSystHists->push_back(currStopFileSystVec);
        }
    }
}
*/


void HistogramVecGrabber_Signal(vector<TFile *> * inputFilesSignal, vector<float> * sigSkimEffScaleCalcVec, unsigned int whichIndex, TH1 * &histStopCentValTH1, TString mcPlotName, vector<TH1 *> * &vecStopTH1SystHists, vector<SystT> * systVec, TString subSampName, bool doOverflow, bool doUnderflow, bool doSyst, bool allMT2llSystematic) {
    TString fileName, systName, mcGrabName;
    TH1 * currHist;
    vector <float> * nVtxSignalBackScaleVec = ScaleBackVecCalc(inputFilesSignal);
    TString mcSystPlotName, systCompHistName;
    mcGrabName = "";
    fileName = inputFilesSignal->at(whichIndex)->GetName();
    cout << "fileName" << " is " << fileName << endl;
    mcGrabName += mcPlotName;
    mcGrabName += subSampName;
    cout << "mcGrabName " << mcGrabName << endl;
    cout << "nVtx scale factor " << nVtxSignalBackScaleVec->at(whichIndex) << endl;
    cout << "SkimEfficiency scale factor " << sigSkimEffScaleCalcVec->at(whichIndex) << endl;
    histStopCentValTH1 = (TH1 *) inputFilesSignal->at(whichIndex)->Get(mcGrabName);
    histStopCentValTH1->Scale(nVtxSignalBackScaleVec->at(whichIndex)); // correct for PURW changes to integral            
    histStopCentValTH1->Scale(sigSkimEffScaleCalcVec->at(whichIndex)); // correct for PURW changes to integral 
    cout << "histStopCentValTH1 bin 1 " << histStopCentValTH1->GetName() << histStopCentValTH1->GetBinContent(1) << endl;
    bool grabSyst;
    if (doSyst) {
        for (unsigned int iSyst = 0; iSyst < systVec->size(); ++iSyst) {
            grabSyst = true;
            cout << "iSyst " << iSyst << endl;
            systName = systVec->at(iSyst).name;
            mcGrabName = "";
            mcGrabName += mcPlotName;
            cout << "systName " << systName << endl;
            cout << "mcGrabName " << mcGrabName << endl;
            mcSystPlotName = mcGrabName;
            mcSystPlotName += systName;
            mcGrabName += subSampName;
            mcSystPlotName += subSampName;
            cout << "mcSystPlotName " << mcSystPlotName << endl;            
            if (systName.Contains("genTop")) {
                grabSyst = false;
            }
            if (!mcGrabName.Contains("MT2ll")) {
                if (systName.Contains("MT2ll")) {
                    if (!allMT2llSystematic && !fileName.Contains("TTBar")) {
                        grabSyst = false;
                    }
                    else {
                        grabSyst = true;
                    }        
                }
                if (systName.Contains("MT2UncES")) {
                    grabSyst = false;
                }
            }
            if (!grabSyst) {
                currHist = (TH1 *) histStopCentValTH1->Clone(histStopCentValTH1->GetName() + systName);
            }
            else {
                currHist = (TH1 *) inputFilesSignal->at(whichIndex)->Get(mcSystPlotName); 
                currHist->Scale(nVtxSignalBackScaleVec->at(whichIndex)); // correct for PURW changes to integral
                currHist->Scale(sigSkimEffScaleCalcVec->at(whichIndex)); // correct for skim efficiencies to integral  
            }
            cout << "curr Syst " << currHist->GetName() << " is bin 1 " << currHist->GetBinContent(1) << endl;
            vecStopTH1SystHists->push_back(currHist);
        }
    }       
}


void HistogramVecGrabberCentValGrab(vector<TFile *> * inputFiles, bool grabData, vector<TH1 *> * vecHist, vector<float> * nVtxBackScaleVec, TString histPlotName, TString subSampName, bool useDDEstimate, float scaleFacTTBar) {
    TString fileName;
    TH1 * currHist;
    for (unsigned int iFile = 0; iFile < inputFiles->size(); ++iFile) {
        fileName = inputFiles->at(iFile)->GetName();
        cout << "fileName" << " is " << fileName << endl;
        cout << "hist plot name " << histPlotName << endl;
        currHist = (TH1 *) inputFiles->at(iFile)->Get(histPlotName);
        cout << "currHist " << currHist->Integral() << endl;
        if (grabData) {
            if (!fileName.Contains("Data")) continue;
        }
        else {
            if (fileName.Contains("Data")) continue;
            currHist->Scale(nVtxBackScaleVec->at(iFile)); // correct for PURW changes to integral
            if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);  
            cout << "currHist bin 1 " << currHist->GetBinContent(1) << endl;
        }
        vecHist->push_back(currHist);
        cout << "pushed hist into vector for iFile " << iFile << endl;
    }
}
void HistogramVecGrabberSystGrab(vector<TFile *> * inputFiles, vector<TH1 *> * vecCentValTH1Hist, vector<TH1 *> * mcCompHistSystVec, vector<float> * nVtxBackScaleVec, TString histPlotName, TString subSampName, vector<SystT> * systVec, bool useDDEstimate, float scaleFacTTBar, bool allMT2llSystematic, int whichNTuple) {
    TString fileName, systName, mcGrabName;
    TH1 * currHist;
    TH1 * systCompHist;
    TString mcSystPlotName, systCompHistName;
    bool grabSyst;
    int grabIndex;
    TH1 * centValTH1Hist;
    for (unsigned int iSyst = 0; iSyst < systVec->size(); ++iSyst) {
        grabSyst = true;
        systCompHist = NULL;
        mcGrabName = "";
        mcGrabName += histPlotName;
        cout << "iSyst " << iSyst << endl;
        systName = systVec->at(iSyst).name;
        cout << "systName " << systName << endl;
        cout << "mcGrabName " << mcGrabName << endl;
        if (!mcGrabName.Contains("MT2ll") && (systName.Contains("MT2ll") || systName.Contains("MT2UncES"))) continue;
        if (systName.Contains("genStop")) grabSyst = false;
        mcSystPlotName = mcGrabName;
        mcSystPlotName += systName;
        mcGrabName += subSampName;
        mcSystPlotName += subSampName;
        grabIndex = -1;
        for (unsigned int iFile = 0; iFile < inputFiles->size(); ++iFile) {
            fileName = inputFiles->at(iFile)->GetName();
            cout << "fileName " << fileName << endl;
            cout << "mcSystPlotName " << mcSystPlotName << endl;
            if (fileName.Contains("Data")) continue;
            ++grabIndex;
            if (systName.Contains("genTop")) {
                if (!fileName.Contains("TTBar")) grabSyst = false;
                if (whichNTuple == 1 && !fileName.Contains("Sig")) grabSyst = false;
            }
            if (systName.Contains("MT2ll")) {
                if (!allMT2llSystematic && !fileName.Contains("TTBar")) {
                    grabSyst = false;
                }
                else {
                    grabSyst = true;
                }        
            }
            if (!grabSyst) {
                centValTH1Hist = vecCentValTH1Hist->at(grabIndex);
                currHist = (TH1 *) centValTH1Hist->Clone(centValTH1Hist->GetName() + systName);
                cout << " currHist name " << currHist->GetName() << endl;
//                patsyHist = (TH1 *) inputFiles->at(iFile)->Get(mcGrabName);
//                currHist = (TH1 *) patsyHist->Clone(patsyHist->GetName() + systName);
            }
            else {
                currHist = (TH1 *) inputFiles->at(iFile)->Get(mcSystPlotName);   
                currHist->Scale(nVtxBackScaleVec->at(iFile)); // correct for PURW changes to integral
                if (useDDEstimate && fileName.Contains("TTBar")) currHist->Scale(scaleFacTTBar);
            }  
            if (systCompHist == NULL) {
                systCompHistName = currHist->GetName();
                systCompHistName += "_Comp";
                systCompHist = (TH1 *) currHist->Clone(systCompHistName);
            }
            else {
                systCompHist->Add(currHist);
            }
            cout << "systCompHist bin 1 " << systCompHist->GetBinContent(1) << endl;
        }
        mcCompHistSystVec->push_back(systCompHist);
    }
}


void HistogramAdderData(vector<TH1 *> * dataHistVec, TH1F * &DataComp, int RBNX, int numDims, int whichAxis, int axis1LB, int axis1UB, int axis2LB, int axis2UB, TString nameProject, bool doOverflow, bool doUnderflow, TString addName) {
    
    TString patsyNamebase = "patsy";
    patsyNamebase += addName;
    TString patsyName;
    
    TH2F * data2DPatsy;
    TH3F * data3DPatsy;
    
    TString dataName = dataHistVec->at(0)->GetName();
    if (numDims > 1) {
        dataName += nameProject;
    }
    dataName += "_DataComp";
    dataName += addName;
    
    switch (numDims) {
        case 1:
            DataComp = (TH1F *) dataHistVec->at(0)->Clone(dataName);
            for (unsigned int i = 1; i < dataHistVec->size(); ++i) {
                DataComp->Add((TH1F*) dataHistVec->at(i));
            }     
            break;
        case 2:
            data2DPatsy = (TH2F *) dataHistVec->at(0);                    
            switch (whichAxis) {
                case 1:
                    DataComp = (TH1F *) data2DPatsy->ProjectionX(dataName, axis1LB, axis1UB, "e");                    
                    for (unsigned int iData = 1; iData < dataHistVec->size(); ++iData) {
                        patsyName = patsyNamebase;
                        patsyName += iData;
                        data2DPatsy = (TH2F *) dataHistVec->at(iData);
                        DataComp->Add((TH1F *) data2DPatsy->ProjectionX(dataName + patsyName, axis1LB, axis1UB, "e"));
                    }            
                    break;
                case 2:
                    DataComp = (TH1F *) data2DPatsy->ProjectionY(dataName, axis1LB, axis1UB, "e");                    
                    for (unsigned int iData = 1; iData < dataHistVec->size(); ++iData) {
                        data2DPatsy = (TH2F *) dataHistVec->at(iData);
                        patsyName = patsyNamebase;
                        patsyName += iData;
                        DataComp->Add((TH1F *) data2DPatsy->ProjectionY(dataName + patsyName, axis1LB, axis1UB, "e"));
                    }           
                    break;                   
                default:
                    break;
            }
            break;
        case 3:
            data3DPatsy = (TH3F *) dataHistVec->at(0);                  
            switch (whichAxis) {
                case 1:
                    DataComp = (TH1F *) data3DPatsy->ProjectionX(dataName, axis1LB, axis1UB, axis2LB, axis2UB, "e");                    
                    for (unsigned int iData = 1; iData < dataHistVec->size(); ++iData) {
                        data3DPatsy = (TH3F *) dataHistVec->at(iData);
                        patsyName = patsyNamebase;
                        patsyName += iData;
                        DataComp->Add((TH1F *) data3DPatsy->ProjectionX(dataName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                    }
                    break;
                case 2:
                    DataComp = (TH1F *) data3DPatsy->ProjectionY(dataName, axis1LB, axis1UB, axis2LB, axis2UB, "e");                    
                    for (unsigned int iData = 1; iData < dataHistVec->size(); ++iData) {
                        data3DPatsy = (TH3F *) dataHistVec->at(iData);
                        patsyName = patsyNamebase;
                        patsyName += iData;
                        DataComp->Add((TH1F *) data3DPatsy->ProjectionY(dataName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                    }
                    break;
                case 3:
                    DataComp = (TH1F *) data3DPatsy->ProjectionZ(dataName, axis1LB, axis1UB, axis2LB, axis2UB, "e");                    
                    for (unsigned int iData = 1; iData < dataHistVec->size(); ++iData) {
                        data3DPatsy = (TH3F *) dataHistVec->at(iData);
                        patsyName = patsyNamebase;
                        patsyName += iData;
                        DataComp->Add((TH1F *) data3DPatsy->ProjectionZ(dataName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                    }           
                    break;
            }
            break;
        default:
            break;
    }
    //Rebin stuff as needed
    DataComp->RebinX(RBNX);
    cout << "DataComp error bin 25 " << DataComp->GetBinError(25) << endl;
    cout << "DataComp Content bin 25 " << DataComp->GetBinContent(25) << endl;
    HistogramUnderflowOverflow(DataComp, doUnderflow, doOverflow);
}
void HistogramAdderMC(vector<TH1 *> * mcIndHistCentValVec, vector<TH1F *> * mcCompHistCentValVec, vector<int> * sampleStartPositions, vector<TString> * sampleAddNames, TH1F * &MCComp, int RBNX, int numDims, int whichAxis, int axis1LB, int axis1UB, int axis2LB, int axis2UB, TString nameProject, bool doOverflow, bool doUnderflow, TString addName) {
    TH1F * currSampHist1D;
    TH2F * currSampHist2DPatsy, * mcComp2DPatsy;
    TH3F * currSampHist3DPatsy, * mcComp3DPatsy;
    TString nameMCComp = mcIndHistCentValVec->at(0)->GetName() + TString("_MCComp");
    TString currSampName;
    
    TString patsyNamebase = "patsy";
    patsyNamebase += addName;
    TString patsyName;
    switch (numDims) {
        case 1:
            MCComp = (TH1F *) mcIndHistCentValVec->at(0)->Clone(nameMCComp);
            for (unsigned int iMC = 1; iMC < mcIndHistCentValVec->size(); ++iMC) {
                MCComp->Add((TH1F*) mcIndHistCentValVec->at(iMC));
            }  
            break;
        case 2:
            mcComp2DPatsy = (TH2F *) mcIndHistCentValVec->at(0);
            switch (whichAxis) {
                case 1:
                    MCComp = (TH1F *) mcComp2DPatsy->ProjectionX(nameMCComp, axis1LB, axis1UB, "e");
                    break;
                case 2:
                    MCComp = (TH1F *) mcComp2DPatsy->ProjectionY(nameMCComp, axis1LB, axis1UB, "e");
                    break;
                default:
                    break;
            }
            for (unsigned int iMC = 1; iMC < mcIndHistCentValVec->size(); ++iMC) {
                patsyName = patsyNamebase;
                patsyName += iMC;
                mcComp2DPatsy = (TH2F *) mcIndHistCentValVec->at(iMC);
                switch (whichAxis) {
                    case 1:
                        MCComp->Add((TH1F*) mcComp2DPatsy->ProjectionX(nameMCComp + patsyName, axis1LB, axis1UB, "e"));
                        break;
                    case 2:
                        MCComp->Add((TH1F*) mcComp2DPatsy->ProjectionY(nameMCComp + patsyName, axis1LB, axis1UB, "e"));
                        break;
                    default:
                        break;
                }
            }
            break;
        case 3:
            mcComp3DPatsy = (TH3F *) mcIndHistCentValVec->at(0);
            switch (whichAxis) {
                case 1:
                    MCComp = (TH1F *) mcComp3DPatsy->ProjectionX(nameMCComp, axis1LB, axis1UB, axis2LB, axis2UB, "e");
                    break;
                case 2:
                    MCComp = (TH1F *) mcComp3DPatsy->ProjectionY(nameMCComp, axis1LB, axis1UB, axis2LB, axis2UB, "e");
                    break;
                case 3:
                    MCComp = (TH1F *) mcComp3DPatsy->ProjectionZ(nameMCComp, axis1LB, axis1UB, axis2LB, axis2UB, "e");
                    break;
                default:
                    break;
            }
            for (unsigned int iMC = 1; iMC < mcIndHistCentValVec->size(); ++iMC) {
                patsyName = patsyNamebase;
                patsyName += iMC;
                mcComp3DPatsy = (TH3F *) mcIndHistCentValVec->at(iMC);
                switch (whichAxis) {
                    case 1:
                        MCComp->Add((TH1F*) mcComp3DPatsy->ProjectionX(nameMCComp + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                        break;
                    case 2:
                        MCComp->Add((TH1F*) mcComp3DPatsy->ProjectionY(nameMCComp + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                        break;
                    case 3:
                        MCComp->Add((TH1F*) mcComp3DPatsy->ProjectionZ(nameMCComp + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                        break;
                    default:
                        break;
                }
            }
            break;
        default:
            break;
    }
    MCComp->RebinX(RBNX);
    cout << "MCComp error bin 25 " << MCComp->GetBinError(25) << endl;
    cout << "MCComp Content bin 25 " << MCComp->GetBinContent(25) << endl;
    HistogramUnderflowOverflow(MCComp, doUnderflow, doOverflow);
/*
    float newHistErr;
    int NBins = MCComp->GetNbinsX();
    if (doUnderflow) {
        MCComp->SetBinContent(1, MCComp->GetBinContent(1) + MCComp->GetBinContent(0));           
        newHistErr = TMath::Sqrt(MCComp->GetBinError(1)*MCComp->GetBinError(1) + MCComp->GetBinError(0)*MCComp->GetBinError(0));
        MCComp->SetBinError(1, newHistErr);
    }
    if (doOverflow) {
        MCComp->SetBinContent(NBins, MCComp->GetBinContent(NBins) + MCComp->GetBinContent(NBins+1)); 
        newHistErr = TMath::Sqrt(MCComp->GetBinError(NBins)*MCComp->GetBinError(NBins) + MCComp->GetBinError(NBins+1)*MCComp->GetBinError(NBins+1));
        MCComp->SetBinError(NBins, newHistErr);
    }
    */
    for (unsigned int iSamp = 0; iSamp < sampleStartPositions->size() - 1; ++iSamp) {
        currSampHist1D = NULL;
        currSampHist2DPatsy = NULL;
        currSampHist3DPatsy = NULL;
        
        currSampName = mcIndHistCentValVec->at(sampleStartPositions->at(iSamp))->GetName();
        if (numDims > 1) currSampName += nameProject;
        currSampName += sampleAddNames->at(iSamp);
        currSampName += addName;
        cout << "stop point " << sampleStartPositions->at(iSamp + 1) << endl;
        for (int iMCSpec = sampleStartPositions->at(iSamp); iMCSpec < sampleStartPositions->at(iSamp + 1); ++iMCSpec) {
            cout << "curr SampName " << currSampName << endl;
            cout << "iMCSpec " << iMCSpec << endl;
            switch (numDims) {
                case 1:
                    if (currSampHist1D == NULL) {
                        currSampHist1D = (TH1F *) mcIndHistCentValVec->at(iMCSpec)->Clone(currSampName);
                    }
                    else {
                        currSampHist1D->Add((TH1F *) mcIndHistCentValVec->at(iMCSpec));
                    }
                    break;
                case 2:
                    if (currSampHist2DPatsy == NULL) {
                        currSampHist2DPatsy = (TH2F *) mcIndHistCentValVec->at(iMCSpec);
                        switch (whichAxis) {
                            case 1:
                                currSampHist1D = (TH1F *) currSampHist2DPatsy->ProjectionX(currSampName, axis1LB, axis1UB, "e");
                                break;
                            case 2:
                                currSampHist1D = (TH1F *) currSampHist2DPatsy->ProjectionY(currSampName, axis1LB, axis1UB, "e");
                                break;
                        }
                    }
                    else {
                        patsyName = patsyNamebase;
                        patsyName += iMCSpec;
                        currSampHist2DPatsy = (TH2F *) mcIndHistCentValVec->at(iMCSpec);
                        switch (whichAxis) {
                            case 1:
                                currSampHist1D->Add((TH1F*) currSampHist2DPatsy->ProjectionX(currSampName + patsyName, axis1LB, axis1UB, "e"));
                                break;
                            case 2:
                                currSampHist1D->Add((TH1F*) currSampHist2DPatsy->ProjectionY(currSampName + patsyName, axis1LB, axis1UB, "e"));
                                break;
                            default:
                                break;
                        }
                    }
                    break;
                case 3:
                    if (currSampHist3DPatsy == NULL) {
                        currSampHist3DPatsy = (TH3F *) mcIndHistCentValVec->at(iMCSpec);
                        switch (whichAxis) {
                            case 1:
                                currSampHist1D = (TH1F *) currSampHist3DPatsy->ProjectionX(currSampName, axis1LB, axis1UB, axis2LB, axis2UB, "e");
                                break;
                            case 2:
                                currSampHist1D = (TH1F *) currSampHist3DPatsy->ProjectionY(currSampName, axis1LB, axis1UB, axis2LB, axis2UB, "e");
                                break;
                            case 3:
                                currSampHist1D = (TH1F *) currSampHist3DPatsy->ProjectionZ(currSampName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e");
                                break;
                        }
                    }
                    else {
                        patsyName = patsyNamebase;
                        patsyName += iMCSpec;
                        currSampHist3DPatsy = (TH3F *) mcIndHistCentValVec->at(iMCSpec);
                        switch (whichAxis) {
                            case 1:
                                currSampHist1D->Add((TH1F*) currSampHist3DPatsy->ProjectionX(currSampName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                                break;
                            case 2:
                                currSampHist1D->Add((TH1F*) currSampHist3DPatsy->ProjectionY(currSampName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                                break;
                            case 3:
                                currSampHist1D->Add((TH1F*) currSampHist3DPatsy->ProjectionZ(currSampName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e"));
                                break;
                            default:
                                break;
                        }
                    }
                    break;
                default:
                    break;
                    
            }
        }
        currSampHist1D->RebinX(RBNX);
        cout << "currSampHist1D " << currSampHist1D->GetName() << " error bin 25 " << currSampHist1D->GetBinError(25) << endl;
        cout << "currSampHist1D " << currSampHist1D->GetName() << " Content bin 25 " << currSampHist1D->GetBinContent(25) << endl;
        HistogramUnderflowOverflow(currSampHist1D, doUnderflow, doOverflow);
        mcCompHistCentValVec->push_back(currSampHist1D);
    }
}



TH2F * HistogramAdderMCTwoDee(vector<TH1 *> * mcIndHistCentValVec) {
    TString cloneName = mcIndHistCentValVec->at(0)->GetName() + TString("_MCComp");;
    TH2F * outHist = (TH2F *) mcIndHistCentValVec->at(0)->Clone(cloneName);
    for (unsigned int iMC = 1; iMC < mcIndHistCentValVec->size(); ++iMC) {
        outHist->Add((TH2F*) mcIndHistCentValVec->at(iMC));
    }  
    return outHist;
}

TH2F * SystHistFinderTwoDee(vector<TH1 *> * vecSystHists, TString searchString) {
    TString name;
    int whichIndex;
    for (unsigned int iSyst = 0; iSyst < vecSystHists->size(); ++iSyst) {
        name = vecSystHists->at(iSyst)->GetName();
        if (name.Contains(searchString)) whichIndex = iSyst;
    }
    return (TH2F *) vecSystHists->at(whichIndex);
}

void HistogramProjectorSyst(vector<TH1 *> * systHistVec, vector<TH1F *> * systProjHistVec, int RBNX, int numDims, int whichAxis, int axis1LB, int axis1UB, int axis2LB, int axis2UB, TString nameProject, bool doOverflow, bool doUnderflow, TString addName) {
    TString systName;
    TString patsyNamebase = "patsy";
    patsyNamebase += addName;
    TString patsyName;
    
    TH1F * currSystHist;
    TH2F * currSystHist2DPatsy;
    TH3F * currSystHist3DPatsy;
    for (unsigned int iSyst = 0; iSyst < systHistVec->size(); ++iSyst) {
        systName = systHistVec->at(iSyst)->GetName();
        if (numDims > 1) systName += nameProject;
        systName += "_SystComp";
        systName += addName;
        switch (numDims) {
            case 1:
                currSystHist = (TH1F *) systHistVec->at(iSyst)->Clone(systName);
                break;                
            case 2:
                currSystHist2DPatsy = (TH2F *) systHistVec->at(iSyst);
                switch (whichAxis) {
                    case 1:
                        currSystHist = (TH1F *) currSystHist2DPatsy->ProjectionX(systName + patsyName, axis1LB, axis1UB, "e");    
                        break;
                    case 2:
                        currSystHist = (TH1F *) currSystHist2DPatsy->ProjectionY(systName + patsyName, axis1LB, axis1UB, "e");    
                        break;
                }
                break;
            case 3:
                currSystHist3DPatsy = (TH3F *) systHistVec->at(iSyst);
                switch (whichAxis) {
                    case 1:
                        currSystHist = (TH1F *) currSystHist3DPatsy->ProjectionX(systName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e");    
                        break;
                    case 2:
                        currSystHist = (TH1F *) currSystHist3DPatsy->ProjectionY(systName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e");    
                        break;
                    case 3:
                        currSystHist = (TH1F *) currSystHist3DPatsy->ProjectionZ(systName + patsyName, axis1LB, axis1UB, axis2LB, axis2UB, "e");    
                        break;
                }
                break;
        }
        currSystHist->RebinX(RBNX);
        float newHistErr;
        int NBins = currSystHist->GetNbinsX();
        if (doUnderflow) {
            currSystHist->SetBinContent(1, currSystHist->GetBinContent(1) + currSystHist->GetBinContent(0));           
            newHistErr = TMath::Sqrt(currSystHist->GetBinError(1)*currSystHist->GetBinError(1) + currSystHist->GetBinError(0)*currSystHist->GetBinError(0));
            currSystHist->SetBinError(1, newHistErr);
        }
        if (doOverflow) {
            currSystHist->SetBinContent(NBins, currSystHist->GetBinContent(NBins) + currSystHist->GetBinContent(NBins+1)); 
            newHistErr = TMath::Sqrt(currSystHist->GetBinError(NBins)*currSystHist->GetBinError(NBins) + currSystHist->GetBinError(NBins+1)*currSystHist->GetBinError(NBins+1));
            currSystHist->SetBinError(NBins, newHistErr);
        }
        systProjHistVec->push_back(currSystHist);
    }
}

void HistogramAdderSignal(TH1 * signalHist, TH1F * &signalProjectionHist, int RBNX, int numDims, int whichAxis, int axis1LB, int axis1UB, int axis2LB, int axis2UB, TString nameProject, bool doOverflow, bool doUnderflow, TString addName) {
    TH2F * signal2DPatsy;
    TH3F * signal3DPatsy;    
    TString signalName = signalHist->GetName();
    if (numDims > 1) {
        signalName += nameProject;
    }
    signalName += "_signalProjectionHist";
    signalName += addName;
    
    switch (numDims) {
        case 1:
            signalProjectionHist = (TH1F *) signalHist;
            break;
        case 2:
            signal2DPatsy = (TH2F *) signalHist;
            switch (whichAxis) {
                case 1:
                    signalProjectionHist = (TH1F *) signal2DPatsy->ProjectionX(signalName, axis1LB, axis1UB, "e");
                    break;
                case 2:
                    signalProjectionHist = (TH1F *) signal2DPatsy->ProjectionY(signalName, axis1LB, axis1UB, "e");
                    break;                   
                default:
                    break;
            }
            break;
        case 3:
            signal3DPatsy = (TH3F *) signalHist;
            switch (whichAxis) {
                case 1:
                    signalProjectionHist = (TH1F *) signal3DPatsy->ProjectionX(signalName, axis1LB, axis1UB, axis2LB, axis2UB, "e");                    
                    break;
                case 2:
                    signalProjectionHist = (TH1F *) signal3DPatsy->ProjectionY(signalName, axis1LB, axis1UB, axis2LB, axis2UB, "e");                    
                    break;
                case 3:
                    signalProjectionHist = (TH1F *) signal3DPatsy->ProjectionZ(signalName, axis1LB, axis1UB, axis2LB, axis2UB, "e");                    
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    signalProjectionHist->RebinX(RBNX);
    HistogramUnderflowOverflow(signalProjectionHist, doUnderflow, doOverflow);
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
    bool grabSyst;
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
            if (!mcGrabName.Contains("MT2ll") && (systName.Contains("MT2ll") || systName.Contains("MT2UncES"))) continue;
            mcSystPlotName = mcGrabName;
            mcSystPlotName += systName;
            mcGrabName += subSampName;
            mcSystPlotName += subSampName;
            if (systName.Contains("genStop")) continue;
            for (unsigned int k = 0; k < inputFiles->size(); ++k) {
                grabSyst = true;
                cout << "k " << k << endl;            
                fileName = inputFiles->at(k)->GetName();
                cout << "fileName " << fileName << endl;
                cout << "mcSystPlotName " << mcSystPlotName << endl;
                if (fileName.Contains("Data")) continue;
                if (systName.Contains("genTop")) {
                    if (!fileName.Contains("TTBar")) grabSyst = false;
                    if (whichNTuple == 1 && !fileName.Contains("Sig")) grabSyst = false;
                }
                if (!allMT2llSystematic && systName.Contains("MT2ll") && !fileName.Contains("TTBar")) {
                    grabSyst = false;
                }
                else {
                    grabSyst = true;
                }    
                if (grabSyst) {
                    currHist = (TH2F *) inputFiles->at(k)->Get(mcSystPlotName);
                }
                else {
                    currHist = (TH2F *) inputFiles->at(k)->Get(mcGrabName);
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



void SystGraphMakers(TH1F * inputBaseMCHist, vector<TH1F *> * inputBaseMCSystHists, vector<TGraphAsymmErrors *> * errCompSpecSource, vector<TGraphAsymmErrors *> * errCompSpecSource_pStat, vector<TGraphAsymmErrors *> * fracRatioSystVec, vector<TString> * systCanvNameVec, Color_t colorErrGraph, TString plotVarName, bool doAbsRatio, float fracRatioYAxisRange, bool doSymErr, bool doSignal) {
    
    TString stringLepEffSF = "LepEffSF";
    TString stringLepEnSc = "LepES";
    TString stringJetEnSc = "JetES";
    TString stringUncEnSc = "UncES";
    TString stringMT2ll = "MT2ll";  
    TString stringGenTopRW = "genTopRW";
    TString stringStopXSecUncert = "genStopXSec";
    
    TGraphAsymmErrors * errCompStatCentVal, * errSystQuadSum, * errSystQuadSum_pStat;
    TGraphAsymmErrors * errLepEnSc, * errLepEffSF, * errJetEnSc, * errUncEnSc, * errMT2ll, * errGenTopRW, * errStopXSecUncert;
    TGraphAsymmErrors * errLepEnSc_pStat, * errLepEffSF_pStat, * errJetEnSc_pStat, * errUncEnSc_pStat, * errMT2ll_pStat, * errGenTopRW_pStat, * errStopXSecUncert_pStat;
    TGraphAsymmErrors * currFracRatioGraph;

    cout << "plotVarName " << plotVarName << endl;
    
    errCompStatCentVal = clonePoints(inputBaseMCHist);
    errLepEnSc = GraphSystErrorSet_SingleSource(inputBaseMCHist, inputBaseMCSystHists, stringLepEnSc + TString("Shift"), doSymErr, 0);
    errLepEnSc_pStat = GraphSystErrorSumErrors(errCompStatCentVal, stringLepEnSc, errLepEnSc, inputBaseMCHist);
    GraphMainAttSet(errLepEnSc, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    GraphMainAttSet(errLepEnSc_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    errCompSpecSource->push_back(errLepEnSc);
    errCompSpecSource_pStat->push_back(errLepEnSc_pStat);
    if (!doSignal) {
        systCanvNameVec->push_back(stringLepEnSc);
        currFracRatioGraph = fracGraph(inputBaseMCHist, errLepEnSc, doAbsRatio, fracRatioYAxisRange);
        fracRatioSystVec->push_back(currFracRatioGraph);
    }
    
    errJetEnSc = GraphSystErrorSet_SingleSource(inputBaseMCHist, inputBaseMCSystHists, stringJetEnSc + TString("Shift"), doSymErr, 0);
    errJetEnSc_pStat = GraphSystErrorSumErrors(errCompStatCentVal, stringJetEnSc, errJetEnSc, inputBaseMCHist);
    GraphMainAttSet(errJetEnSc, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    GraphMainAttSet(errJetEnSc_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    errCompSpecSource->push_back(errJetEnSc);
    errCompSpecSource_pStat->push_back(errJetEnSc_pStat);
    if (!doSignal) {
        systCanvNameVec->push_back(stringJetEnSc);
        currFracRatioGraph = fracGraph(inputBaseMCHist, errJetEnSc, doAbsRatio, fracRatioYAxisRange);
        fracRatioSystVec->push_back(currFracRatioGraph);
    }
    if (plotVarName.Contains("MT2ll")) {
        errUncEnSc = GraphSystErrorSet_SingleSource(inputBaseMCHist, inputBaseMCSystHists, stringUncEnSc + TString("Shift"), doSymErr, 0);
        errUncEnSc_pStat = GraphSystErrorSumErrors(errCompStatCentVal, stringUncEnSc, errUncEnSc, inputBaseMCHist);
        GraphMainAttSet(errUncEnSc, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
        GraphMainAttSet(errUncEnSc_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
        errCompSpecSource->push_back(errUncEnSc);
        errCompSpecSource_pStat->push_back(errUncEnSc_pStat);
        if (!doSignal) {
            systCanvNameVec->push_back(stringUncEnSc);
            currFracRatioGraph = fracGraph(inputBaseMCHist, errUncEnSc, doAbsRatio, fracRatioYAxisRange);
            fracRatioSystVec->push_back(currFracRatioGraph);
        }
    }
    /*
    if (plotVarName.Contains("MT2ll")) {
        errMT2ll = GraphSystErrorSet_SingleSource(inputBaseMCHist, inputBaseMCSystHists, stringMT2ll + TString("Shift"), doSymErr, 1);
        errMT2ll_pStat = GraphSystErrorSumErrors(errCompStatCentVal, stringMT2ll, errMT2ll, inputBaseMCHist);
        GraphMainAttSet(errMT2ll, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
        GraphMainAttSet(errMT2ll_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
        errCompSpecSource->push_back(errMT2ll);
        errCompSpecSource_pStat->push_back(errMT2ll_pStat);
        if (!doSignal) {
            systCanvNameVec->push_back(stringMT2ll);
            currFracRatioGraph = fracGraph(inputBaseMCHist, errMT2ll, doAbsRatio, fracRatioYAxisRange);
            fracRatioSystVec->push_back(currFracRatioGraph);
        }
    }        
    else {
        errMT2ll = NULL;
    }
    */
    
    errLepEffSF = GraphSystErrorSet_SingleSource(inputBaseMCHist, inputBaseMCSystHists, stringLepEffSF + TString("Shift"), doSymErr, 0);
    errLepEffSF_pStat = GraphSystErrorSumErrors(errCompStatCentVal, stringLepEffSF, errLepEffSF, inputBaseMCHist);
    GraphMainAttSet(errLepEffSF, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    GraphMainAttSet(errLepEffSF_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    errCompSpecSource->push_back(errLepEffSF);
    errCompSpecSource_pStat->push_back(errLepEffSF_pStat);
    if (!doSignal) {
        systCanvNameVec->push_back(stringLepEffSF); 
        currFracRatioGraph = fracGraph(inputBaseMCHist, errLepEffSF, doAbsRatio, fracRatioYAxisRange);
        fracRatioSystVec->push_back(currFracRatioGraph);
    }
    
    errGenTopRW = GraphSystErrorSet_SingleSource(inputBaseMCHist, inputBaseMCSystHists, stringGenTopRW, doSymErr, 0);
    errGenTopRW_pStat = GraphSystErrorSumErrors(errCompStatCentVal, stringGenTopRW, errGenTopRW, inputBaseMCHist);
    GraphMainAttSet(errGenTopRW, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    GraphMainAttSet(errGenTopRW_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    errCompSpecSource->push_back(errGenTopRW);
    errCompSpecSource_pStat->push_back(errGenTopRW_pStat);
    if (!doSignal) {
        systCanvNameVec->push_back(stringGenTopRW); 
        currFracRatioGraph = fracGraph(inputBaseMCHist, errGenTopRW, doAbsRatio, fracRatioYAxisRange);
        fracRatioSystVec->push_back(currFracRatioGraph);   
    }
    
    errStopXSecUncert = GraphSystErrorSet_SingleSource(inputBaseMCHist, inputBaseMCSystHists, stringStopXSecUncert + TString("Shift"), doSymErr, 0);
    errStopXSecUncert_pStat = GraphSystErrorSumErrors(errCompStatCentVal, stringStopXSecUncert, errStopXSecUncert, inputBaseMCHist);
    GraphMainAttSet(errStopXSecUncert, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    GraphMainAttSet(errStopXSecUncert_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    errCompSpecSource->push_back(errStopXSecUncert);
    errCompSpecSource_pStat->push_back(errStopXSecUncert_pStat);
    if (!doSignal) {
        systCanvNameVec->push_back(stringStopXSecUncert);
        currFracRatioGraph = fracGraph(inputBaseMCHist, errStopXSecUncert, doAbsRatio, fracRatioYAxisRange);
        fracRatioSystVec->push_back(currFracRatioGraph);
    }
    
    errSystQuadSum = GraphSystErrorSumErrors(errCompStatCentVal, TString("FullSyst"), errCompSpecSource, false, inputBaseMCHist);
    errSystQuadSum_pStat = GraphSystErrorSumErrors(errCompStatCentVal, TString("FullSyst"), errCompSpecSource, true, inputBaseMCHist);
    GraphMainAttSet(errSystQuadSum, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    GraphMainAttSet(errSystQuadSum_pStat, colorErrGraph, 3001, colorErrGraph, 2, kWhite, 0, 0); 
    errCompSpecSource->push_back(errSystQuadSum);
    errCompSpecSource_pStat->push_back(errSystQuadSum_pStat);
    if (!doSignal) {
        systCanvNameVec->push_back("FullSyst");
        currFracRatioGraph = fracGraph(inputBaseMCHist, errSystQuadSum, doAbsRatio, fracRatioYAxisRange);
        fracRatioSystVec->push_back(currFracRatioGraph);
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
            if (!mcGrabName.Contains("MT2ll") && (systName.Contains("MT2ll") || systName.Contains("MT2UncES"))) continue;
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
void HistogramAdderSyst(vector<TH1F *> * dataHistVec, vector<TH1F *> * mcIndHistCentValVec, vector<TH1F *> * mcCompHistCentValVec, TH1F * &DataComp, TH1F * &MCComp, TH1F * &FracComp, int whichNTuple, bool doAbsRatio, float yAxisRange) {
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
void SpectrumDraw(TCanvas * InputCanvas, TH1F * Hist1, TString legHist1, TH1F * Hist2, TH1F * fracRatioHist, TH1F * errHist, THStack * MCStack, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB, bool logYPad1, vector<TString> * mcLegends, vector<TH1F *> * indMCHists, bool doMean, TString cutUsed, float inputLumi, TLegend * &leg, bool doSFR, bool doODFR, TH1F * fracRatioHist_Other) {
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
    cout << "draw test 3d " << mcLegends->size() << endl;
    cout << "draw test 3da " << indMCHists->size() << endl;
    //    leg = new TLegend(legXstart,legYstart - 0.05 * (numIndMC + 1),legXstart + 0.05 * (numIndMC + 1),legYstart);
    leg = new TLegend(legXstart - 0.05,legYstart - 0.05 * (mcLegends->size() + 1),legXstart + 0.40,legYstart);
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
    HistAxisAttSet(XAxis, XAxisTitle, .17, 1.03, .12,.07, 0.0, 0.0);
    YAxis->SetNdivisions(3,5,0);
    XAxis->SetNdivisions(6,5,0);
    float RatioMax;
    float RatioMin;  
    if (doSFR) {
        RatioMax = fracRatioHist->GetBinContent(fracRatioHist->GetMaximumBin());
        RatioMin = fracRatioHist->GetBinContent(fracRatioHist->GetMinimumBin());
        if (abs(RatioMin - 0) < 1E-3) RatioMin = fracRatioHist->GetMinimum(1E-3);
        YAxis->SetRangeUser(RatioMin - 0.1, RatioMax + 0.1);
    }
    TGraphAsymmErrors * fracRatioDrawGraph = clonePoints(fracRatioHist, false);
    fracRatioHist->SetLineColor(kBlack);
    HistToGraphCopyAttributes(fracRatioHist, fracRatioDrawGraph);
    TH1F * patsy = (TH1F*) fracRatioHist->Clone("frac patsy");
    patsy->SetLineColor(kWhite);
    patsy->SetMarkerColor(kWhite);
    patsy->Draw("e1");    
    fracRatioDrawGraph->Draw("p0 same");
    TLegend * fracRatioLegend = new TLegend(0.8, 0.65, 0.95, 0.95);
    TGraphAsymmErrors * fracRatioDrawGraph_Other;
    if (doODFR) {
        fracRatioDrawGraph_Other = clonePoints(fracRatioHist_Other, false);
        fracRatioHist_Other->SetLineColor(kRed);
        fracRatioHist_Other->SetMarkerColor(kRed);
        HistToGraphCopyAttributes(fracRatioHist_Other, fracRatioDrawGraph_Other);
        fracRatioLegend->AddEntry(fracRatioHist, "Selected NTuple Ratio Plot", "pl");
        fracRatioLegend->AddEntry(fracRatioHist_Other, "Other NTuples Ratio Plot", "pl");
        fracRatioDrawGraph_Other->Draw("p0 same"); 
        fracRatioLegend->Draw("same");
    }
    Pad2->Update();
    Pad2->Modified();
}



void SpectrumDraw_AddSignal(TCanvas * InputCanvas, vector<TH1F *> * vecSignalHists, vector<TString> * vecSignalLegends, TLegend * leg) {
    TPad * Pad1 = (TPad *) InputCanvas->cd(1);
    TH1F * currSigHist;
    TH1F * h_SigErr;
    for (unsigned int iSigPoints = 0; iSigPoints < vecSignalHists->size(); ++iSigPoints) {
        currSigHist = vecSignalHists->at(iSigPoints);
        h_SigErr = (TH1F *) currSigHist->Clone();
        HistMainAttSet(h_SigErr, h_SigErr->GetLineColor(), 3001, h_SigErr->GetLineColor(), 2, kWhite, 0, 0);
        h_SigErr->Draw("e2 same");
        currSigHist->Draw("hist same");
        leg->AddEntry(currSigHist, vecSignalLegends->at(iSigPoints), "l");
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
void SpectrumDrawSyst(TCanvas * InputCanvas, TH1F * Hist1, TString legHist1, TH1F * Hist2, THStack * MCStack, TGraphAsymmErrors * errGraph, TGraphAsymmErrors * errGraphjustSyst, TH1F * fracRatioHist, TGraphAsymmErrors * fracRatioGraph, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB, bool logYPad1, vector<TString> * mcLegends, vector<TH1F *> * indMCHists, bool doMean, TString cutUsed, float inputLumi, TLegend * &leg, bool doSFR, bool doODFR, TH1F * fracRatioHist_Other) {
    TLatex * tl = new TLatex();
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
    float RatioMax;
    float RatioMin;
    if (doSFR) {
        RatioMax = fracRatioHist->GetBinContent(fracRatioHist->GetMaximumBin());
        RatioMin = fracRatioHist->GetBinContent(fracRatioHist->GetMinimumBin());
        if (abs(RatioMin - 0) < 1E-3) RatioMin = fracRatioHist->GetMinimum(1E-3);        
        YAxis->SetRangeUser(RatioMin - 0.1, RatioMax + 0.1);
    }
    TGraphAsymmErrors * fracRatioDrawGraph = clonePoints(fracRatioHist, false);
    HistToGraphCopyAttributes(fracRatioHist, fracRatioDrawGraph);
    fracRatioGraph->SetFillColor(kYellow);
    TH1F * patsy = (TH1F*) fracRatioHist->Clone("frac patsy");
    patsy->SetLineColor(kWhite);
    patsy->SetMarkerColor(kWhite);
    fracRatioHist->SetLineColor(kBlack);
    patsy->Draw("e1");    
    fracRatioGraph->Draw("2 same");
    //    fracRatioHist->Draw("e1 same");
    fracRatioDrawGraph->Draw("p0 same");    
    TLegend * fracRatioLegend = new TLegend(0.8, 0.8, 0.95, 0.95);
    fracRatioLegend->AddEntry(fracRatioGraph, "Syst. Uncert.", "f");
    TGraphAsymmErrors * fracRatioDrawGraph_Other;
    if (doODFR) {
        fracRatioDrawGraph_Other = clonePoints(fracRatioHist_Other, false);
        fracRatioHist_Other->SetLineColor(kRed);
        fracRatioHist_Other->SetMarkerColor(kRed);
        HistToGraphCopyAttributes(fracRatioHist_Other, fracRatioDrawGraph_Other);
        fracRatioLegend->AddEntry(fracRatioHist_Other, "Other NTuples Ratio Plot", "pl");
        fracRatioDrawGraph_Other->Draw("p0 same"); 
    }
    fracRatioLegend->Draw("same");
    Pad2->Update();
    Pad2->Modified();
}


void SpectrumDrawSyst_AddSignal(TCanvas * InputCanvas, vector<TH1F *> * vecSignalHists, vector<vector<TGraphAsymmErrors *> *> * vecSigErrGraph, unsigned int whichSyst, vector<TString> * vecSignalLegends, TLegend * leg) {
    TH1F * currSigHist;
    TGraphAsymmErrors * currSigErrGraph;
    
    TPad * Pad1 = (TPad *) InputCanvas->cd(1);
    for (unsigned int iSigPoints = 0; iSigPoints < vecSignalHists->size(); ++iSigPoints) {
        currSigHist = vecSignalHists->at(iSigPoints);
        currSigErrGraph = vecSigErrGraph->at(iSigPoints)->at(whichSyst);
        cout << "SigErrGraph name " << currSigErrGraph->GetName() << endl;
        cout << "SigErrGraph high error at point 2 " << currSigErrGraph->GetErrorYhigh(2) << endl;
        cout << "SigErrGraph low  error at point 2 " << currSigErrGraph->GetErrorYlow(2) << endl;
        currSigErrGraph->Draw("2 same");
        currSigHist->Draw("hist same");
        leg->AddEntry(currSigHist, vecSignalLegends->at(iSigPoints), "l");
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

void SpectrumDrawSingSampCompare(TCanvas * InputCanvas, TH1F * Hist1, TString legHist1, TH1F * Hist2, TString legHist2, TH1F * fracRatioHist, float TextXPos, float TextYStartPos, float YAxisLB, float YAxisUB, bool logYPad1, bool doMean, TString cutUsed, float inputLumi) {
    TLatex * tl = new TLatex();
    TLegend * leg;
    tl->SetTextAlign(12);
    tl->SetNDC();
    tl->SetTextSize(0.03);
    char buf[99];        
    TPad * Pad1;
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

TH2F * FOMHist(TH2F * inputMCCompHist, TH2F * inputMCSigHist, int whichType) {
    int NBinsX = inputMCCompHist->GetNbinsX();
    int XMax = inputMCCompHist->GetXaxis()->GetBinLowEdge(NBinsX + 1);
    int NBinsY = inputMCCompHist->GetNbinsY();
    int YMax = inputMCCompHist->GetYaxis()->GetBinLowEdge(NBinsY + 1);
    float MCInt, SigInt, FOMCentVal, FOMStatErr;
    TString baseName = "h_FOM_forDiffMT2Cuts";
    switch (whichType) {
        case 0:
            baseName += "_CentVal";
            break;
        case 1:
            baseName += "_LepESShiftUp";
            break;
        case 2:
            baseName += "_LepESShiftDown";
            break;
        case 3:
            baseName += "_JetESShiftUp";
            break;
        case 4:
            baseName += "_JetESShiftDown";
            break;
        case 5:
            baseName += "_UncESShiftUp";
            break;
        case 6:
            baseName += "_UncESShiftDown";
            break;
        case 7:
            baseName += "_LepEffSFShiftUp";
            break;
        case 8:
            baseName += "_LepEffSFShiftDown";
            break;            
        case 9:
            baseName += "_genTopRW";
            break;
        case 10:
            baseName += "_StopXSecShiftUp";
            break;
        case 11:
            baseName += "_StopXSecShiftDown";
            break;
        default:
            break;
    }
    TH2F * outHist = new TH2F(baseName, "F.O.M. = #frac{S}{#sqrt{S + B}}; MT2_{ll} Cut [GeV]; MT2_{(lb)(lb)} Cut [GeV]", NBinsX, 0., XMax, NBinsY, 0., YMax); outHist->Sumw2();
    for (int xInd = 1; xInd < NBinsX + 1; ++xInd) {
        for (int yInd = 1; yInd < NBinsY + 1; ++yInd) {
            MCInt = inputMCCompHist->Integral(xInd, NBinsX + 1, yInd, NBinsY + 1);
            SigInt = inputMCSigHist->Integral(xInd, NBinsX + 1, yInd, NBinsY + 1);
            FOMCentVal = SigInt / TMath::Sqrt(SigInt + MCInt);
            if (xInd == 22 && yInd == 3) {
                cout << "MCInt " << MCInt << endl;
                cout << "SigInt " << SigInt << endl;
                cout << "FOMCentVal " << FOMCentVal << endl;
            }
            FOMStatErr = TMath::Sqrt(SigInt) * (SigInt + 2 * MCInt + SigInt * MCInt) / (2 * TMath::Power(SigInt + MCInt, 1.5));
            outHist->SetBinContent(xInd, yInd, FOMCentVal);
            outHist->SetBinError(xInd, yInd, FOMStatErr);
        }
    }
    return outHist;
}

