#include <vector>
#include "./StopFunctionDefinitions_v2.h"
//#include <boost/format.hpp>
#include <sstream>
using namespace std;

TString WhichTTBarString(int whichTTBarSyst, int whichNTuple) {    
    TString TTBarSystString;
    switch (whichNTuple) {
        case 0:
            switch (whichTTBarSyst) {
                case 0:
                    TTBarSystString = "_Powheg";
                    break;
                case 1:
                    TTBarSystString = "_MCatNLO";
                    break;
                case 2:
                    TTBarSystString = "_Madgraph";
                    break;
                default:
                    break;
            }
            break;
        case 1:
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
            break;
        default:
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
float OviedoNEvents(TString sourceName) {
    // Numbers taken from https://docs.google.com/spreadsheet/ccc?key=0Aq5OAopf_dtsdDJqREJReGQyY21wbERldVFSZVJHbFE&hl=en_US&pli=1#gid=57
    if (sourceName.Contains("TTbar_Powheg")) return 28150723;
    if (sourceName.Contains("TTbar_Madgraph")) return 6923750;
    if (sourceName.Contains("TTbar_MCatNLO")) return 32852589;
    if (sourceName.Contains("TTJetsFullHadrMG")) return 28123821;
    if (sourceName.Contains("TTJetsSemiLeptMG")) return 24424818;
    if (sourceName.Contains("TTJetsFullLeptMGtauola")) return 11673402;
    if (sourceName.Contains("TTJetsFullLeptMG")) return 12114013;
    if (sourceName.Contains("TWDilep")) return 2976510;
    if (sourceName.Contains("TbarWDilep")) return 2971482;
    if (sourceName.Contains("DYJets_Madgraph")) return 7132223;
    if (sourceName.Contains("DY50toLLMad1jet")) return 24045248;
    if (sourceName.Contains("DY50toLLMad2jet")) return 2352304;
    if (sourceName.Contains("DY50toLLMad3jet")) return 11015445;
    if (sourceName.Contains("DY50toLLMad4jet")) return 6402827;
    if (sourceName.Contains("ZJets_Madgraph")) return 30459503;
    if (sourceName.Contains("ggWWto2L")) return 109987;
    if (sourceName.Contains("WWTo2L2Nu_Madgraph")) return 1933235;
    if (sourceName.Contains("WGstarToElNuMad")) return 314653;
    if (sourceName.Contains("WGstarToMuNuMad")) return 299973;
    if (sourceName.Contains("WGstarToTauNuMad")) return 49997;
    if (sourceName.Contains("WZ")) return 10000283;
    if (sourceName.Contains("ZZ")) return 9799908;
    if (sourceName.Contains("WJets_Madgraph")) return 76102995;
    if (sourceName.Contains("ZgammaToLLG")) return 6588161;
    if (sourceName.Contains("WgammaToLNuG")) return 4802358;
    if (sourceName.Contains("HWW125")) return 299975;
    if (sourceName.Contains("VBF125")) return 299687;
    if (sourceName.Contains("HZZ4L")) return 995117;
    if (sourceName.Contains("WWGJets")) return 215121;
    if (sourceName.Contains("WZZJets")) return 219835;
    if (sourceName.Contains("ZZZJets")) return 224904;
    if (sourceName.Contains("WWZJets")) return 222234;
    if (sourceName.Contains("WWWJets")) return 220549;
    if (sourceName.Contains("TTWJets")) return 196046;
    if (sourceName.Contains("TTZJets")) return 210160;
    if (sourceName.Contains("TTWWJets")) return 197820;
    if (sourceName.Contains("TTGJets")) return 71598;
    cout << "didn't MATCH ANY" << endl;
    return 0;
}

float OviedoNEventsSignal(TString sourceName, int genStopMass, int genChi0Mass, int genCharginoMass) {
    // Numbers taken from a combination of: https://docs.google.com/spreadsheet/ccc?key=0Aq5OAopf_dtsdDJqREJReGQyY21wbERldVFSZVJHbFE&hl=en_US&pli=1#gid=57
    // As well as by hand calculations
    
    ///Signal Dudes....should be handled carefully
//    if (sourceName.Contains("T2tt_FineBin")) return 97050784;
    cout << "didn't MATCH ANY" << endl;
    return 0;
}

vector<double> * WeightBaseVec(int whichNTuple, vector<TList *> * fileListVec, unsigned int whichFile, TString nEventHistName, TString nParFileHistName, vector<int> * numParFilesVec) {
    vector<double> * outVec = new vector<double>;
    TH1F * h_eventCount;
    float nEvents;
    TString sourceName;
    TFile * first_source = (TFile*) fileListVec->at(whichFile)->First();
    sourceName = first_source->GetName();
    cout << "first_source name = " << sourceName << endl;
    if (whichNTuple == 1) {
        h_eventCount = (TH1F*) first_source->Get(nEventHistName);
        nEvents = h_eventCount->Integral();
    }
    else {
        nEvents = OviedoNEvents(sourceName);
    }
    TH1F * h_numParFiles = (TH1F*) first_source->Get(nParFileHistName);
    int numParFiles = (h_numParFiles == NULL) ? 1 : (int) h_numParFiles->GetEntries();
    cout << "nEvents " << nEvents << endl;
    outVec->push_back(nEvents);
    numParFilesVec->push_back(numParFiles);
    TFile * nextsource = (TFile*) fileListVec->at(whichFile)->After(first_source);
    while (nextsource) {
        sourceName = nextsource->GetName();
        cout << "nextsource " << sourceName << endl;
        if (whichNTuple == 1) {
            h_eventCount = (TH1F*) nextsource->Get(nEventHistName);
            nEvents = h_eventCount->Integral();
        }
        else {
            nEvents = OviedoNEvents(sourceName);
        }
        cout << "nEvents " << nEvents << endl;
        h_numParFiles = (TH1F*) nextsource->Get(nParFileHistName);
        numParFiles = (h_numParFiles == NULL) ? 1 : (int) h_numParFiles->GetEntries(); 
        outVec->push_back(nEvents);
        numParFilesVec->push_back(numParFiles);
        nextsource = (TFile*)fileListVec->at(whichFile)->After( nextsource );
    }
    return outVec;
}
vector<double> * WeightVec(int whichNTuple, float L_data, vector<double> * baseWeightVec, unsigned int whichFileList, vector<int> * numParFilesVec, bool doExclusiveSample) {
    vector<double> * outVec = new vector<double>;    
    int numParFiles;
    //Cross Sections mostly taken from: https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //Oviedo Cross Sections taken from: https://docs.google.com/spreadsheet/ccc?key=0Aq5OAopf_dtsdDJqREJReGQyY21wbERldVFSZVJHbFE&hl=en#gid=57
    // ttbar cross section taken from 

    
    ////Top and Single Top cross sections
    // Basic TTBar xsec
    //    double xsecTTbar                = 244.849; //old one used...new one is approx 1.003 times larger
    double xsecTTbar                = 245.8;
    // Specific Decay Mode TTBar; Only Oviedo and only Madgraph
    //branching ratios taken from the Oviedo Spreadsheet
    double branchRatioTTBarHadronic = 106.9/235;
    double branchRatioTTBarFullLeptonic = 103./235.;
    double branchRatioTTBarSemiLeptonic = 0.108 * 0.108 * 9;
    double xsecTTbarHadronic        = xsecTTbar * branchRatioTTBarHadronic;
    double xsecTTbarFullLeptonic    = xsecTTbar * branchRatioTTBarFullLeptonic;
    double xsecTTbarSemiLeptonic   = xsecTTbar * branchRatioTTBarSemiLeptonic;
    // Basic Single Top XSec, note that number is quite different for Oviedo...see below
    double xsecSingTop              = 11.1;
    
    //// Di-Weak Vector Boson cross section
    double xsecWW                   = 54.838;
    double xsecWZ                   = 33.21; // Note Oviedo's different versions for WZ and ZZ
    double xsecZZ                   = 17.654;
    // Oviedo version of DiW production
    double xsecWW2L2Nu = 5.8123;
    double xsecWW4L = 0.182852;
    
    ////Drell Yan
    // "Inclusive DY"
    double xsecZDY10to50            = 860.5;
    double xsecZDY50toInf           = 3503.71;
    // DY NumJet Exclusive samples Note: Oviedo only
    double xsecZDY50toInf_1Jets     = 671.830779661017;
    double xsecZDY50toInf_2Jets     = 216.758237288136;
    double xsecZDY50toInf_3Jets     = 61.1952813559322;
    double xsecZDY50toInf_4Jets     = 27.591766779661;

    //// WJets
    double xsecWLNu                 = 36257.2;
    
    //// QCD
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
    
    //// Weak VBoson + Photon
    double xsecZG                   = 132.6;
    double xsecWG                   = 553.92;
    double xsecWGstarToElNu         = 5.873;
    double xsecWGstarToMuNu         = 1.914;
    double xsecWGstarToTauNu        = 0.336;

        
    /////Rare backgrounds
    //Higgs samples
    double xsecHiggsVBF             = 0.0356;
    double xsecHiggsWW              = 0.444;
    double xsecHiggsZZ4L            = 0.052428931968;
    //Triple Vector Boson ish
    double xsecWWG                  = 0.528;
    double xsecWWW                  = 0.0822;
    double xsecWWZ                  = 0.0633;
    double xsecWZZ                  = 0.0192;
    double xsecZZZ                  = 0.00459;
    //TTBar + Vector Boson(s)
    double xsecTTbarW               = 0.232;
    double xsecTTbarZ               = 0.174;
    double xsecTTbarG               = 1.44;
    double xsecTTbarWW              = 0.00204;
    //
    if (whichNTuple == 0) {
        // for xsecs that are different between the two types of nTuples
        cout << "whichNTuple is 0 " << endl;
        xsecZDY50toInf = 3532.8;
        xsecSingTop = 1.17334824;
        xsecWZ = 22.44;
        xsecZZ = 9.03;
        xsecWLNu = 37509;
    }
    switch (whichNTuple) {
        case 0:
            switch (whichFileList) {
                case 0:
                    numParFiles = numParFilesVec->at(0);
                    if (doExclusiveSample) {
                        outVec->push_back(numParFiles * L_data * xsecTTbarFullLeptonic / baseWeightVec->at(0));
                    }
                    else {
                        outVec->push_back(numParFiles * L_data * xsecTTbar / baseWeightVec->at(0));
                    }
                    break;
                case 1:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecTTbarSemiLeptonic / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecTTbarHadronic / baseWeightVec->at(1));
                    break;
                case 2:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecSingTop / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecSingTop / baseWeightVec->at(1));
                    break;
                case 3:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecZDY10to50 / baseWeightVec->at(0));
                    if (doExclusiveSample) {
                        numParFiles = numParFilesVec->at(1);                        
                        outVec->push_back(numParFiles * L_data * xsecZDY50toInf_1Jets / baseWeightVec->at(1));
                        numParFiles = numParFilesVec->at(2);                        
                        outVec->push_back(numParFiles * L_data * xsecZDY50toInf_2Jets / baseWeightVec->at(2));
                        numParFiles = numParFilesVec->at(3);                        
                        outVec->push_back(numParFiles * L_data * xsecZDY50toInf_3Jets / baseWeightVec->at(3));
                        numParFiles = numParFilesVec->at(4);                        
                        outVec->push_back(numParFiles * L_data * xsecZDY50toInf_4Jets / baseWeightVec->at(4));
                    }
                    else {
                        numParFiles = numParFilesVec->at(1);                        
                        outVec->push_back(numParFiles * L_data * xsecZDY50toInf / baseWeightVec->at(1));
                    }
                    break;
                case 4:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWW2L2Nu / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecWW4L / baseWeightVec->at(1));
                    break;
                case 5:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWZ / baseWeightVec->at(0));
                    break;
                case 6:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecZZ / baseWeightVec->at(0));
                    break;
                case 7:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWLNu / baseWeightVec->at(0));
                    break;
                case 8:            
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWG / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecWGstarToElNu / baseWeightVec->at(1));
                    numParFiles = numParFilesVec->at(2);
                    outVec->push_back(numParFiles * L_data * xsecWGstarToMuNu / baseWeightVec->at(2));
                    numParFiles = numParFilesVec->at(3);
                    outVec->push_back(numParFiles * L_data * xsecWGstarToTauNu / baseWeightVec->at(3));
                    break;            
                case 9:            
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecZG / baseWeightVec->at(0));
                    break; 
                case 10:            
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecHiggsWW / baseWeightVec->at(0));
                    break;
                case 11:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecHiggsVBF / baseWeightVec->at(0));
                    break;
                case 12:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecHiggsZZ4L / baseWeightVec->at(0));
                    break;
                case 13:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWWG / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecWWW / baseWeightVec->at(1));
                    numParFiles = numParFilesVec->at(2);
                    outVec->push_back(numParFiles * L_data * xsecWWZ / baseWeightVec->at(2));
                    numParFiles = numParFilesVec->at(3);
                    outVec->push_back(numParFiles * L_data * xsecWZZ / baseWeightVec->at(3));
                    numParFiles = numParFilesVec->at(4);
                    outVec->push_back(numParFiles * L_data * xsecZZZ / baseWeightVec->at(4));
                    break;
                case 14:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecTTbarG / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecTTbarZ / baseWeightVec->at(1));
                    numParFiles = numParFilesVec->at(2);
                    outVec->push_back(numParFiles * L_data * xsecTTbarW / baseWeightVec->at(2));
                    numParFiles = numParFilesVec->at(3);
                    outVec->push_back(numParFiles * L_data * xsecTTbarWW / baseWeightVec->at(3));
                    break;
                default:
                    break;
            }
            break;
        case 1:
            switch (whichFileList) {
                case 0:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecTTbar / baseWeightVec->at(0));
                    break;
                case 1:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecTTbar / baseWeightVec->at(0));
                    break;
                case 2:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecSingTop / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecSingTop / baseWeightVec->at(1));
                    break;
                case 3:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecZDY10to50 / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecZDY50toInf / baseWeightVec->at(1));
                    break;
                case 4:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWW / baseWeightVec->at(0));
                    break;
                case 5:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWZ / baseWeightVec->at(0));
                    break;
                case 6:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecZZ / baseWeightVec->at(0));
                    break;
                case 7:
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecWLNu / baseWeightVec->at(0));
                    break;
                case 8:            
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecQCDMu15 / baseWeightVec->at(0));
                    break;            
                case 9:            
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecQCDEM20to30 / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecQCDEM30to80 / baseWeightVec->at(1));
                    numParFiles = numParFilesVec->at(2);
                    outVec->push_back(numParFiles * L_data * xsecQCDEM80to170 / baseWeightVec->at(2));
                    break; 
                case 10:            
                    numParFiles = numParFilesVec->at(0);
                    outVec->push_back(numParFiles * L_data * xsecQCDBCEM20to30 / baseWeightVec->at(0));
                    numParFiles = numParFilesVec->at(1);
                    outVec->push_back(numParFiles * L_data * xsecQCDBCEM30to80 / baseWeightVec->at(1));
                    numParFiles = numParFilesVec->at(2);
                    outVec->push_back(numParFiles * L_data * xsecQCDBCEM80to170 / baseWeightVec->at(2));
                    break; 
                default:
                    break;
            }
            break;
        default:
            break;
    }
    
    return outVec;
}
vector<TFile *> * OutFileVec(int whichNTuple, vector<TString> * nameVec, vector<bool> * boolVec, bool doExclusiveSample) {
    vector<TFile *> * outVec = new vector<TFile *>;
    TFile * TargetTTBarSig, * TargetTTBarBkg, * TargetSingTop, * TargetZDY, * TargetWLNu, * TargetWW, * TargetWZ, * TargetZZ;
    TFile * TargetQCDMu, * TargetQCDEM, * TargetQCDBCEM;
    TFile * TargetZG, * TargetWG, * TargetHiggsWW, * TargetHiggsVBF, * TargetHiggsZZ4L;
    TFile * TargetTripleVecBoson, * TargetTTBarVecBoson;
    TString TTBarSystString, whichNTupleString, PURWString, doSystString, suffixString;
    if (nameVec->size() < 5) {
        cout << "nameVec size less than 5!!!" << endl;
    }
    TString exclusiveString = (doExclusiveSample) ? "_Exclusive" : "";
    
    TTBarSystString = nameVec->at(0);
    whichNTupleString = nameVec->at(1);
    PURWString = nameVec->at(2);
    doSystString = nameVec->at(3);
    suffixString = TString("Haddplots.root");
    if (boolVec->at(0)) {
        if (doExclusiveSample) {
            TargetTTBarSig = TFile::Open(TString("TTBarSig") + TTBarSystString + whichNTupleString + exclusiveString + PURWString + doSystString + suffixString, "RECREATE");
        }
        else if (whichNTupleString.Contains("DESY")) {
            TargetTTBarSig = TFile::Open(TString("TTBarSig") + TTBarSystString + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        }
        else {
            TargetTTBarSig = TFile::Open(TString("TTBarComp") + TTBarSystString + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        }
        outVec->push_back(TargetTTBarSig);
        TargetTTBarBkg = TFile::Open(TString("TTBarBkg") + TTBarSystString + whichNTupleString + exclusiveString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetTTBarBkg);
    }
    if (boolVec->at(2)) {
        TargetSingTop = TFile::Open(TString("SingleTop") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
        outVec->push_back(TargetSingTop);
    }
    if (boolVec->at(3)) {
        TargetZDY = TFile::Open(TString("ZDY") + whichNTupleString + exclusiveString + PURWString + doSystString + suffixString, "RECREATE");
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
    switch (whichNTuple) {
        case 0:
            if (boolVec->at(8)) {
                TargetWG = TFile::Open(TString("WG") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
                outVec->push_back(TargetWG);
                TargetZG = TFile::Open(TString("ZG") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
                outVec->push_back(TargetZG);
            }        
            if (boolVec->at(10)) {
                TargetHiggsWW = TFile::Open(TString("HiggsWW") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");                
                outVec->push_back(TargetHiggsWW);
                TargetHiggsVBF = TFile::Open(TString("HiggsVBF") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE"); 
                outVec->push_back(TargetHiggsVBF);
                TargetHiggsZZ4L = TFile::Open(TString("HiggsZZ4L") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE"); 
                outVec->push_back(TargetHiggsZZ4L);
            } 
            if (boolVec->at(13)) {
                TargetTripleVecBoson = TFile::Open(TString("TripVecBoson") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");                
                outVec->push_back(TargetTripleVecBoson);
                TargetTTBarVecBoson = TFile::Open(TString("TTBarVecBoson") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");                
                outVec->push_back(TargetTTBarVecBoson);
            }
            break;
        case 1:
            if (boolVec->at(8)) {
                TargetQCDMu = TFile::Open(TString("QCDMu") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
                outVec->push_back(TargetQCDMu);
                TargetQCDEM = TFile::Open(TString("QCDEM") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");   
                outVec->push_back(TargetQCDEM);
                TargetQCDBCEM = TFile::Open(TString("QCDBCEM") + whichNTupleString + PURWString + doSystString + suffixString, "RECREATE");
                outVec->push_back(TargetQCDBCEM);
            }            
            break;
        default:
            break;
    }
    return outVec;
}

vector<TList *> * FileListVec(int whichNTuple, vector<TString> * nameVec, vector<bool> * boolVec, bool doExclusiveSample) {
    vector<TList *> * outVec = new vector<TList *>;
    TList * FileListTTBarSig, * FileListTTBarBkg, * FileListSingTop, * FileListZDY, * FileListWLNu, * FileListWW, * FileListWZ, * FileListZZ, * FileListQCDMu, * FileListQCDEM, * FileListQCDBCEM, * FileListZG, * FileListWG, * FileListHiggsWW, * FileListHiggsVBF, * FileListHiggsZZ4L; 
    TList * FileListTripVecBoson, * FileListTTBarVecBoson;
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
            if (boolVec->at(0)) {
                FileListTTBarSig = new TList();
                if (doExclusiveSample) {
//                    FileListTTBarSig->Add(TFile::Open(TString("TTJetsFullLeptMG") + whichNTupleString + PURWString + doSystString + suffixString));
                    FileListTTBarSig->Add(TFile::Open(TString("TTJetsFullLeptMGtauola") + whichNTupleString + PURWString + doSystString + suffixString));
                }
                else {
                    FileListTTBarSig->Add(TFile::Open(TString("TTbar") + TTBarSystString + whichNTupleString + PURWString + doSystString + suffixString));
                }
                outVec->push_back(FileListTTBarSig);
                FileListTTBarBkg = new TList();
                FileListTTBarBkg->Add(TFile::Open(TString("TTJetsFullHadrMG") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTTBarBkg->Add(TFile::Open(TString("TTJetsSemiLeptMG") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListTTBarBkg);
            }
            if (boolVec->at(2)) {
                FileListSingTop = new TList();
                FileListSingTop->Add(TFile::Open(TString("TWDilep") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListSingTop->Add(TFile::Open(TString("TbarWDilep") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListSingTop);
            }
            if (boolVec->at(3)) {
                FileListZDY = new TList();
                FileListZDY->Add(TFile::Open(TString("DYJets_Madgraph") + whichNTupleString + PURWString + doSystString + suffixString));
                if (doExclusiveSample) {
                    FileListZDY->Add(TFile::Open(TString("DY50toLLMad1jet") + whichNTupleString + PURWString + doSystString + suffixString));
                    FileListZDY->Add(TFile::Open(TString("DY50toLLMad2jet") + whichNTupleString + PURWString + doSystString + suffixString));
                    FileListZDY->Add(TFile::Open(TString("DY50toLLMad3jet") + whichNTupleString + PURWString + doSystString + suffixString));
                    FileListZDY->Add(TFile::Open(TString("DY50toLLMad4jet") + whichNTupleString + PURWString + doSystString + suffixString));
                }
                else {
                    FileListZDY->Add(TFile::Open(TString("ZJets_Madgraph") + whichNTupleString + PURWString + doSystString + suffixString));
                }
                outVec->push_back(FileListZDY);
            }
            if (boolVec->at(4)) {
                FileListWW = new TList();
                FileListWW->Add(TFile::Open(TString("WWTo2L2Nu_Madgraph") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListWW->Add(TFile::Open(TString("ggWWto2L") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListWW);
                FileListWZ = new TList();
                FileListWZ->Add(TFile::Open(TString("WZ") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListWZ);
                FileListZZ = new TList();
                FileListZZ->Add(TFile::Open(TString("ZZ") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListZZ);
            }
            if (boolVec->at(7)) {
                FileListWLNu = new TList();
                FileListWLNu->Add(TFile::Open(TString("WJets_Madgraph") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListWLNu);
            }
            if (boolVec->at(8)) {
                FileListWG = new TList();
                FileListWG->Add(TFile::Open(TString("WgammaToLNuG") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListWG->Add(TFile::Open(TString("WGstarToElNuMad") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListWG->Add(TFile::Open(TString("WGstarToMuNuMad") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListWG->Add(TFile::Open(TString("WGstarToTauNuMad") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListWG);
                FileListZG = new TList();
                FileListZG->Add(TFile::Open(TString("ZgammaToLLG") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListZG);
            }
            if (boolVec->at(10)) {
                FileListHiggsWW = new TList();
                FileListHiggsWW->Add(TFile::Open(TString("HWW125") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListHiggsWW);
                FileListHiggsVBF = new TList();
                FileListHiggsVBF->Add(TFile::Open(TString("VBF125") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListHiggsVBF);
                FileListHiggsZZ4L = new TList();
                FileListHiggsZZ4L->Add(TFile::Open(TString("HZZ4L") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListHiggsZZ4L);
            }
            if (boolVec->at(13)) {
                FileListTripVecBoson = new TList();
                FileListTripVecBoson->Add(TFile::Open(TString("WWG") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTripVecBoson->Add(TFile::Open(TString("WWW") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTripVecBoson->Add(TFile::Open(TString("WWZ") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTripVecBoson->Add(TFile::Open(TString("WZZ") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTripVecBoson->Add(TFile::Open(TString("ZZZ") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListTripVecBoson); 
                FileListTTBarVecBoson = new TList();
                FileListTripVecBoson->Add(TFile::Open(TString("TTGJets") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTripVecBoson->Add(TFile::Open(TString("TTZJets") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTripVecBoson->Add(TFile::Open(TString("TTWJets") + whichNTupleString + PURWString + doSystString + suffixString));
                FileListTripVecBoson->Add(TFile::Open(TString("TTWWJets") + whichNTupleString + PURWString + doSystString + suffixString));
                outVec->push_back(FileListTTBarVecBoson); 
            }
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