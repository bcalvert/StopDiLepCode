#include <vector>
#include "./StopFunctionDefinitions_v2.h"
//#include <boost/format.hpp>
#include <sstream>
using namespace std;

TString WhichTTBarString(int whichTTBarSyst) {    
    TString TTBarSystString;
    switch (whichTTBarSyst) {
        case 0:
            TTBarSystString = "_powheg";
            break;
        case 1:
            TTBarSystString = "_mcatnlo";
            break;
        case 2:
            TTBarSystString = "";
            break;
        case 3:
            TTBarSystString = "_scaleup";
            break;
        case 4:
            TTBarSystString = "_scaledown";
            break;
        case 5:
            TTBarSystString = "_massup";
            break;
        case 6:
            TTBarSystString = "_massdown";
            break;
        case 7:
            TTBarSystString = "_matchingup";
            break;
        case 8:
            TTBarSystString = "_matchingdown";
            break;
        default:
            cout << "Not within range!!" << endl;
            return "";
            break;
    }
    return TTBarSystString;
}
void WeightVecFiller(TList * sourcelist, vector<double> * weightVec, TString histName) {
    TFile * first_source = (TFile*) sourcelist->First();
    cout << "first_source name = " << first_source->GetName() << endl;
    TH1F * h_eventCount = (TH1F*) first_source->Get(histName);
    float nEvents = h_eventCount->Integral();
    cout << "nEvents " << nEvents << endl;
    weightVec->push_back(nEvents);
    TFile *nextsource = (TFile*)sourcelist->After( first_source );
    while (nextsource) {
        cout << "nextsource " << nextsource->GetName() << endl;
        h_eventCount = (TH1F*) nextsource->Get(histName);
        nEvents = h_eventCount->Integral();
        cout << "nEvents " << nEvents << endl;
        weightVec->push_back(nEvents);
        nextsource = (TFile*)sourcelist->After( nextsource );
    }
}
vector<double> * WeightBaseVec(vector<TList *> * fileListVec, unsigned int whichFile, TString nEventHistName) {
    vector<double> * outVec = new vector<double>;
    /*
     if (boolVec->at(whichFile)) {
     WeightVecFiller(fileListVec->at(whichFile), outVec, nEventHistName);
     }
     else {
     cout << "file got turned off!! " << endl;
     }    
     */
    WeightVecFiller(fileListVec->at(whichFile), outVec, nEventHistName);
    return outVec;
}
vector<double> * WeightVec(float L_data, vector<double> * baseWeightVec, unsigned int whichFile) {
    vector<double> * outVec = new vector<double>;    
    //Cross Sections mostly taken from: https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    double xsecTTbar                = 244.849;
    double xsecSingTop              = 11.1;
    double xsecWW                   = 54.838;
    double xsecWZ                   = 33.21;
    double xsecZZ                   = 17.654;
    double xsecZDY10to50            = 860.5;
    //double xsecZDY50toInf         = 3532.8;
    double xsecZDY50toInf           = 3503.71;
    double xsecWLNu                 = 36257.2;
    double xsecQCDMu15              = 3.640E8*3.7E-4;
    double xsecQCDMu20to30          = 2.870E8*6.500E-3;
    double xsecQCDMu30to50          = 6.609E7*12.20E-3;
    double xsecQCDMu50to80          = 8.802E6*21.80E-3;
    double xsecQCDMu80to120         = 1.024E6*39.50E-3;
    double xsecQCDMu120to170        = 1.578E5*47.30E-3;
    double xsecQCDEM20to30          = 2.886E8*10.10E-3;
    double xsecQCDEM30to80          = 7.433E7*62.10E-3;
    double xsecQCDEM80to170         = 1.191E6*153.9E-3;
    double xsecQCDBCEM20to30        = 2.886E8*5.800E-4;
    double xsecQCDBCEM30to80        = 7.424E7*2.250E-3;
    double xsecQCDBCEM80to170       = 1.191E6*10.90E-3;
    switch (whichFile) {
        case 0:
            outVec->push_back(L_data * xsecTTbar / baseWeightVec->at(0));
            break;
        case 1:
            outVec->push_back(L_data * xsecTTbar / baseWeightVec->at(0));
            break;
        case 2:
            outVec->push_back(L_data * xsecSingTop / baseWeightVec->at(0));
            outVec->push_back(L_data * xsecSingTop / baseWeightVec->at(1));
            break;
        case 3:
            outVec->push_back(L_data * xsecZDY10to50 / baseWeightVec->at(0));
            outVec->push_back(L_data * xsecZDY50toInf / baseWeightVec->at(1));
            break;
        case 4:
            outVec->push_back(L_data * xsecWW / baseWeightVec->at(0));
            break;
        case 5:
            outVec->push_back(L_data * xsecWZ / baseWeightVec->at(0));
            break;
        case 6:
            outVec->push_back(L_data * xsecZZ / baseWeightVec->at(0));
            break;
        case 7:
            outVec->push_back(L_data * xsecWLNu / baseWeightVec->at(0));
            break;
        case 8:            
            outVec->push_back(L_data * xsecQCDMu15 / baseWeightVec->at(0));
            break;            
        case 9:            
            outVec->push_back(L_data * xsecQCDEM20to30 / baseWeightVec->at(0));
            outVec->push_back(L_data * xsecQCDEM30to80 / baseWeightVec->at(1));
            outVec->push_back(L_data * xsecQCDEM80to170 / baseWeightVec->at(2));
            break; 
        case 10:            
            outVec->push_back(L_data * xsecQCDBCEM20to30 / baseWeightVec->at(0));
            outVec->push_back(L_data * xsecQCDBCEM30to80 / baseWeightVec->at(1));
            outVec->push_back(L_data * xsecQCDBCEM80to170 / baseWeightVec->at(2));
            break; 
        default:
            break;
    }
    return outVec;
}
vector<TFile *> * OutFileVec(vector<TString> * nameVec, vector<bool> * boolVec) {
    vector<TFile *> * outVec = new vector<TFile *>;
    TFile * TargetTTBarSig, * TargetTTBarBkg, * TargetSingTop, * TargetZDY, * TargetWLNu, * TargetWW, * TargetWZ, * TargetZZ, * TargetQCDMu, * TargetQCDEM, * TargetQCDBCEM;
    TString TTBarSystString, whichNTupleString, PURWString, doSystString, suffixString;
    if (nameVec->size() < 5) {
        cout << "nameVec size less than 5!!!" << endl;
    }
    TTBarSystString = nameVec->at(0);
    whichNTupleString = nameVec->at(1);
    PURWString = nameVec->at(2);
    doSystString = nameVec->at(3);
    suffixString = TString("Haddplots.root");
    if (boolVec->at(0)) {
        TargetTTBarSig = TFile::Open(TString("TTBarSig") + TTBarSystString + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetTTBarSig);
        TargetTTBarBkg = TFile::Open(TString("TTBarBkg") + TTBarSystString + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetTTBarBkg);
    }
    if (boolVec->at(2)) {
        TargetSingTop = TFile::Open(TString("SingleTop") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetSingTop);
    }
    if (boolVec->at(3)) {
        TargetZDY = TFile::Open(TString("ZDY") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetZDY);
    }
    if (boolVec->at(4)) {
        TargetWW = TFile::Open(TString("WW") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
        outVec->push_back(TargetWW);
        TargetWZ = TFile::Open(TString("WZ") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
        outVec->push_back(TargetWZ);
        TargetZZ = TFile::Open(TString("ZZ") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE"); 
        outVec->push_back(TargetZZ);
    }
    if (boolVec->at(7)) {
        TargetWLNu = TFile::Open(TString("WToLNu") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetWLNu);
    }
    if (boolVec->at(8)) {
        TargetQCDMu = TFile::Open(TString("QCDMu") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
        outVec->push_back(TargetQCDMu);
        TargetQCDEM = TFile::Open(TString("QCDEM") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
        outVec->push_back(TargetQCDEM);
        TargetQCDBCEM = TFile::Open(TString("QCDBCEM") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetQCDBCEM);
    }
    return outVec;
}

vector<TList *> * FileListVec(int whichNTuple, vector<TString> * nameVec, vector<bool> * boolVec) {
    vector<TList *> * outVec = new vector<TList *>;
    TList * FileListTTBarSig, * FileListTTBarBkg, * FileListSingTop, * FileListZDY, * FileListWLNu, * FileListWW, * FileListWZ, * FileListZZ, * FileListQCDMu, * FileListQCDEM, * FileListQCDBCEM;    
    TString TTBarSystString, whichNTupleString, PURWString, doSystString, suffixString;
    if (nameVec->size() < 5) {
        cout << "nameVec size less than 5!!!" << endl;
    }
    TTBarSystString = nameVec->at(0);
    whichNTupleString = nameVec->at(1);
    PURWString = nameVec->at(2);
    doSystString = nameVec->at(3);
    suffixString = nameVec->at(4);
    switch (whichNTuple) {
        case 0:
            cout << "not currently supported!!" << endl;
            break;
        case 1:
            if (boolVec->at(0)) {
                FileListTTBarSig = new TList();
                FileListTTBarSig->Add(TFile::Open(TString("ttbarsignalplustau") + TTBarSystString + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListTTBarSig);
                FileListTTBarBkg = new TList();
                FileListTTBarBkg->Add(TFile::Open(TString("ttbarbg") + TTBarSystString + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListTTBarBkg);
            }
            if (boolVec->at(2)) {
                FileListSingTop = new TList();
                FileListSingTop->Add(TFile::Open(TString("singletop_tw") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListSingTop->Add(TFile::Open(TString("singleantitop_tw") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListSingTop);
            }
            if (boolVec->at(3)) {
                FileListZDY = new TList();
                FileListZDY->Add(TFile::Open(TString("dy1050") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListZDY->Add(TFile::Open(TString("dy50inf") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListZDY);
            }
            if (boolVec->at(4)) {
                FileListWW = new TList();
                FileListWW->Add(TFile::Open(TString("wwtoall") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListWW);
                FileListWZ = new TList();
                FileListWZ->Add(TFile::Open(TString("wztoall") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListWZ);
                FileListZZ = new TList();
                FileListZZ->Add(TFile::Open(TString("zztoall") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListZZ);
            }
            if (boolVec->at(7)) {
                FileListWLNu = new TList();
                FileListWLNu->Add(TFile::Open(TString("wtolnu") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListWLNu);
            }
            if (boolVec->at(8)) {
                FileListQCDMu = new TList();
                FileListQCDMu->Add(TFile::Open(TString("qcdmu15") + whichNTupleString + PURWString + doSystString + suffixString));
                //        FileListQCDMu->Add(TFile::Open(TString("qcdmu15") + whichNTupleString + PURWString + doSystString + suffixString));                
                outVec->push_back(FileListQCDMu);
                FileListQCDEM = new TList();
                FileListQCDEM->Add(TFile::Open(TString("qcdem2030") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListQCDEM->Add(TFile::Open(TString("qcdem3080") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListQCDEM->Add(TFile::Open(TString("qcdem80170") + whichNTupleString + PURWString + doSystString + suffixString));                
                outVec->push_back(FileListQCDEM);
                FileListQCDBCEM = new TList();
                FileListQCDBCEM->Add(TFile::Open(TString("qcdbcem2030") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListQCDBCEM->Add(TFile::Open(TString("qcdbcem3080") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListQCDBCEM->Add(TFile::Open(TString("qcdbcem80170") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListQCDBCEM);
            }
            break;    
        default:
            break;
    }
    return outVec;
}