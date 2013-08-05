#ifndef HistSampFunc_h_
#define HistSampFunc_h_

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
//#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TMath.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TF1.h"
// Header file for the classes stored in the TTree if any.                                                                            
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <vector>
#include <vector>
#include <Math/GenVector/LorentzVector.h>

#include <iostream>
//#include <vector>
#include <cmath>
#include <sstream>
#include <map>
using namespace std;

typedef struct {
    TString name;
    float   fVar;
    double  dVar;
    string  systVarKey;
    int     whichSystType;   // 0 = universal systematic, 1 = lepton systematic, 2 = jet systematic, 3 = other systematic
} SystT;

typedef struct {
    TString name;
    TString xLabel;
    int xBinN;
    int RBNX;
    float xMin, xMax;
    string xVarKey;
    bool doXSyst;
    
    TString yLabel;
    int yBinN;
    int RBNY;
    float yMin, yMax; 
    string yVarKey;
    bool doYSyst;
    
    TString zLabel;
    int zBinN;
    float zMin, zMax; 
    int RBNZ;
    string zVarKey;
    bool doZSyst;
    
    bool logY1D;
} HistogramT;

typedef struct {
    //    SampleT() : histNameSuffix(""), whichdiLepType(-1), doZVeto(-1), cutNJets(-1), cutNBJets(-1), cutMET(0) {}
    TString histNameSuffix;
    TString histXaxisSuffix;
    TString histYaxisSuffix;
    TString histZaxisSuffix;
    /// variables to store what kinds of
    int     whichdiLepType; // -1: inclusive, 0: MuMu, 1: EE, 2: EMu
    int     doZVeto;   // -1: inclusive, 0: ZMass window, 1: outside ZMass window;
    int     cutNJets;  // -1: inclusive, # > 0: require NJets >= #
    int     cutNBJets;  // -1: inclusive, # > 0: require NBJets >= #
    float   cutMET; // # > 0: require MET >= #
    bool    blindDataChannel;
    
} SampleT;

typedef struct {
    TString HistTXLabel;
    TString HistTYLabel;
    TString HistTZLabel;
    TString SampTLabel;
    
    int newXBinN;
    float newXMin, newXMax;
    
    int newYBinN;
    float newYMin, newYMax;
    
    int newZBinN;
    float newZMin, newZMax;    
    
} SpecHistBinT;


typedef std::pair<HistogramT, SampleT> histKey;
inline bool operator<(const histKey &a, const histKey &b)
{
    return (a.first.name < b.first.name) || (a.first.name == b.first.name && a.second.histNameSuffix < b.second.histNameSuffix);
}
inline bool operator<(const SampleT &a, const SampleT &b)
{
    return (a.histNameSuffix < b.histNameSuffix) || (a.histXaxisSuffix < b.histXaxisSuffix && a.histNameSuffix == b.histNameSuffix);
}
inline SampleT operator+(const SampleT &a, const SampleT &b)
{
    SampleT outSampleT;
    outSampleT.histNameSuffix = a.histNameSuffix + b.histNameSuffix;
    outSampleT.histXaxisSuffix = a.histXaxisSuffix + b.histXaxisSuffix;
    outSampleT.histYaxisSuffix = a.histYaxisSuffix + b.histYaxisSuffix;
    outSampleT.histZaxisSuffix = a.histZaxisSuffix + b.histZaxisSuffix;
    return outSampleT;
}
inline float dPhi(float phi1, float phi2) {
    float result = phi1-phi2;
    while (result >= TMath::Pi()) result -= 2*TMath::Pi();
    while (result < -1*TMath::Pi()) result += 2*TMath::Pi();
    return fabs(result);
}
inline float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float dphi = dPhi(phi1,phi2);
    float deta = eta1-eta2;
    float result = dphi*dphi+deta*deta;
    result = sqrt(result);
    return result;
}

inline double dPhi(double phi1, double phi2) {
    double result = phi1-phi2;
    while (result >= TMath::Pi()) result -= 2*TMath::Pi();
    while (result < -1*TMath::Pi()) result += 2*TMath::Pi();
    return fabs(result);
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dphi = dPhi(phi1,phi2);
    double deta = eta1-eta2;
    double result = dphi*dphi+deta*deta;
    result = sqrt(result);
    return result;
}
/*
 currently not working with this error:
 error: passing 'const TString' as 'this' argument of 'TString& TString::operator=(const TString&)' discards qualifiers [-fpermissive]
 void operator+=(const SampleT &a, const SampleT &b)
 {
 a.histNameSuffix = a.histNameSuffix + b.histNameSuffix;
 a.histXaxisSuffix += a.histXaxisSuffix + b.histXaxisSuffix;
 a.histYaxisSuffix += a.histYaxisSuffix + b.histYaxisSuffix;
 a.histZaxisSuffix += a.histZaxisSuffix + b.histZaxisSuffix;
 }
 */

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
typedef std::vector<LV> VLV;

inline void TreeBranchSet(TTree * analTree, int whichNTuple = 0) {
    
}
inline float PileupRW(TH1 * nVtxSFHist, int nVtx) {
    return (float) nVtxSFHist->GetBinContent(nVtxSFHist->FindBin(nVtx));
}

inline vector<float> * SystWeights(TH1 * origHist, TH1 * rwHist) {
    // origHist is histogram for initial distribution, //rwHist is histogram for distribution
    // you wish to RW to -- both should be 1D (functionality for now) and have same nBins
    // reweighting is just bin by bin for now -- another consideration is to do some sort of shape fit
    // Another consideration is to return a histogram as that facilitates "finding" the weight to use a bit more
    vector<float> * weightVec = new vector<float>;
    float currWeight;
    float rwHistBC, origHistBC;
    int nBins = origHist->GetNbinsX();
    if (nBins != rwHist->GetNbinsX()) {
        cout << "ruh roh!!! nBins between two hists aren't the same" << endl;
        return weightVec;
    }
    for (int i = 0; i < nBins+2; ++i) {
        rwHistBC = rwHist->GetBinContent(i);
        origHistBC = origHist->GetBinContent(i);
        currWeight = origHistBC > 0 ? rwHistBC/origHistBC : 0;
        weightVec->push_back(currWeight);
    }
    return weightVec;
}
inline TLorentzVector LeptonScaleSystShift(TLorentzVector inputLepVec, int inputLepPDGID, float shiftDirection) {
    float barrelEtaEnd = 1.442; float endcapEtaStart = 1.566; 
    float electronESEB = 0.006; float electronESEE = 0.015;
    float muonES = 0.002;
    float scaleToUse = 0;
    TLorentzVector outShiftVec = inputLepVec;
    if (abs(inputLepPDGID) == 11) {
        scaleToUse = ((inputLepVec.Eta() < barrelEtaEnd) ? electronESEB : electronESEE);
    }
    else if (abs(inputLepPDGID) == 13) {
        scaleToUse = muonES;
    }
    else {
        cout << "wtf?! not dealing with electron or muon" << endl;
    }
    outShiftVec *= (1 + shiftDirection * scaleToUse);
    return outShiftVec;
}
/*
 inline TLorentzVector JetScaleSystShift(TLorentzVector inputJetVec) {
 
 }
 */
inline void METSystShift(vector<TLorentzVector> * inputObjVec, vector<TLorentzVector> * shiftObjVec, float &newMET, float &newMETPhi, float origMET, float origMETPhi) {
    TVector3 vecMET;
    TVector3 inputObjTotThreeVec; inputObjTotThreeVec.SetPtEtaPhi(0., 0., 0.);
    TVector3 shiftObjTotThreeVec; shiftObjTotThreeVec.SetPtEtaPhi(0., 0., 0.);
    TVector3 inputObjCurrThreeVec;
    TVector3 shiftObjCurrThreeVec;
    vecMET.SetPtEtaPhi(origMET, 0., origMETPhi);    
    if (inputObjVec->size() != shiftObjVec->size()) {
        cout << "inputObj size " << inputObjVec->size() << endl;
        cout << "shiftObj size " << shiftObjVec->size() << endl;
        cout << "issue: two input object vectors are not same size -- check it out" << endl;
        return;
    }
    for (unsigned int i = 0; i < inputObjVec->size(); ++i) {
        inputObjCurrThreeVec.SetPtEtaPhi(inputObjVec->at(i).Pt(), inputObjVec->at(i).Eta(), inputObjVec->at(i).Phi());
        shiftObjCurrThreeVec.SetPtEtaPhi(shiftObjVec->at(i).Pt(), shiftObjVec->at(i).Eta(), shiftObjVec->at(i).Phi());
        inputObjTotThreeVec = inputObjTotThreeVec + inputObjCurrThreeVec;
        shiftObjTotThreeVec = shiftObjTotThreeVec + shiftObjCurrThreeVec;
    }
    vecMET = vecMET - inputObjTotThreeVec;
    vecMET = vecMET + shiftObjTotThreeVec;
    newMET = vecMET.Pt();
    newMETPhi = vecMET.Phi();
    return;
}

inline float GenLevelTopPtWeight(float pT_Top, float pT_AntiTop) {
    if (pT_Top < 0 || pT_AntiTop < 0) {
        cout << "bad mojo -- either pT is less than 0 " << endl;
        cout << "pT_Top " << pT_Top << endl;
        cout << "pT_AntiTop " << pT_AntiTop << endl;
        return 0.;
    }
    float pT_ToUse = TMath::Sqrt(pT_Top * pT_AntiTop);
    /*
     // Testing speed of using one of these guys
     TF1 * genWeight = new TF1("genTopWeight", "expo");
     genWeight->SetParameter("Constant", .156)
     genWeight->SetParameter("Slope", -.00137)
     return genWeight->Eval(pT_ToUse);
    */
    float expoConst = 0.156;
    float expoSlope = -0.00137;
/*    return TMath::Exp(expoConst) * TMath::Exp(expoSlope * pT_ToUse);*/
//    cout << "TMath::Exp(expoConst) * TMath::Exp(expoSlope * pT_ToUse) " << (TMath::Exp(expoConst) * TMath::Exp(expoSlope * pT_ToUse)) << endl;
    float weightTop = TMath::Exp(expoConst) * TMath::Exp(expoSlope * pT_Top);
    float weightAntiTop = TMath::Exp(expoConst) * TMath::Exp(expoSlope * pT_AntiTop);
//    cout << " weight Top " << weightTop << endl;
//    cout << " weight antiTop " << weightAntiTop << endl;
//    cout << "TMath::Sqrt(weightAntiTop * weightAntiTop + weightTop * weightTop) " << TMath::Sqrt((weightAntiTop * weightAntiTop) + (weightTop * weightTop)) << endl;
//    return TMath::Sqrt(weightAntiTop * weightAntiTop + weightTop * weightTop);
    return TMath::Sqrt(weightAntiTop * weightTop);
}

inline int EventType(vector<int> *lepPdgId, int lep0Index, int lep1Index) {
    int Lep0Type = lepPdgId->at(lep0Index);
    int Lep1Type = lepPdgId->at(lep1Index);
    int Type;
    if (TMath::Abs(Lep0Type) == 11) {
        if (TMath::Abs(Lep1Type) == 11) {
            Type = 1;
        }
        else if (TMath::Abs(Lep1Type) == 13) {
            Type = 2;
        }
    }
    else if (TMath::Abs(Lep0Type) == 13) {
        if (TMath::Abs(Lep1Type) == 11) {
            Type = 2;
        }
        else if (TMath::Abs(Lep1Type) == 13) {
            Type = 0;
        }
    }
    else {
        cout <<"" << endl;
        cout <<"something is wrong" << endl;
        cout <<"Lep0 ID: " << Lep0Type << endl;
        cout <<"Lep1 ID: " << Lep1Type << endl;
        cout <<"" << endl;
    }
    return Type;
}
inline float ScaleFactorMC(int Type, int Syst) {
    
    // Type: 0 MuMu, 1 EE 2 EMu
    // Syst: 0 CentVal, 1 SystShiftUp, 2 SystShiftDown
    //    float SFTrig[3] = {0.994, 0.955, 0.978};
    //    float SFIDIso[3] = {0.993, 0.979, 0.986};    
    //return SFTrig[Type] * SFIDIso[Type];
    float SFTotal[3] = {0.987, 0.957, 0.935};
    float SFSystUp[3] = {0.011, 0.015, 0.013};
    float SFSystDown[3] = {0.011, 0.015, 0.013};
    switch (Syst) {
        case 0:
            return SFTotal[Type];
            break;
        case 1:
            return SFTotal[Type] + SFSystUp[Type];
            break;      
        case -1:
            return SFTotal[Type] - SFSystDown[Type];
            break;  
        default:
            return 1;
            break;
    }
}
inline bool FilterMET(vector<bool> * regFilterVec, vector<bool> * oppFilterVec) {
    for (unsigned int iReg = 0; iReg < regFilterVec->size(); ++iReg) {
        if (!regFilterVec->at(iReg)) return false;
    }
    for (unsigned int iOpp = 0; iOpp < oppFilterVec->size(); ++iOpp) {
        if (oppFilterVec->at(iOpp)) return false;
    }
    return true;
}

inline vector<int> * ElectronPickOvi(vector<float> * ElecPt, vector<float> * ElecEta, vector<int> * ElecCharge, vector<float> * ElecNeutHadIso, vector<float> * ElecCharHadIso, vector<float> * ElecPhotHadIso, vector<bool> * passConvVeto, vector<bool> * isPFElectron) {
    //    cout << "doing electron info" << endl;
    int leadElecIndex = -1; int subElecIndex = -1; int leadPositIndex = -1; int subPositIndex = -1;
    float elecIsoRatioCut = 0.15; float elecEtaCut = 2.5; float leadElecPtCut = 20; float subElecPtCut = 10;
    float barrelEtaEnd = 1.442; float endcapEtaStart = 1.566;
    //    cout << "i max " << ElecPt->size() << endl;
    for (unsigned int i = 0; i < ElecPt->size(); ++i) {
        /*
         cout << "i " << i << endl;
         cout << "passconv " << passConvVeto->at(i) << endl;
         cout << "isPFElec " << isPFElectron->at(i) << endl;
         cout << "neuthadiso " << ElecNeutHadIso->at(i) << endl;
         cout << "charhadiso " << ElecCharHadIso->at(i) << endl;
         cout << "photiso " << ElecPhotHadIso->at(i) << endl;
         cout << "ElecPt " << ElecPt->at(i) << endl;
         cout << "ElecEta " << ElecEta->at(i) << endl;
         cout << "ElecEta " << abs(ElecEta->at(i)) << endl;
         cout << "ElecCharge " << ElecCharge->at(i) << endl;
         cout << "(ElecNeutHadIso->at(i) + ElecCharHadIso->at(i) + ElecPhotHadIso->at(i))/ElecPt->at(i) " << (ElecNeutHadIso->at(i) + ElecCharHadIso->at(i) + ElecPhotHadIso->at(i))/ElecPt->at(i) << endl;
         */
        if (isPFElectron->at(i) && passConvVeto->at(i) && ((ElecNeutHadIso->at(i) + ElecCharHadIso->at(i) + ElecPhotHadIso->at(i))/ElecPt->at(i)) < elecIsoRatioCut && abs(ElecEta->at(i)) < elecEtaCut && !(abs(ElecEta->at(i)) > barrelEtaEnd && abs(ElecEta->at(i)) < endcapEtaStart)) {
            if (ElecCharge->at(i) > 0) {
                if (leadPositIndex < 0) {
                    if (ElecPt->at(i) > leadElecPtCut) leadPositIndex = i;
                }
                else if (subPositIndex < 0) {
                    if (ElecPt->at(i) > subElecPtCut) subPositIndex = i;
                }
                else {
                    continue;
                }
            }
            else {
                if (leadElecIndex < 0) {
                    if (ElecPt->at(i) > leadElecPtCut) leadElecIndex = i;
                }
                else if (subElecIndex < 0) {
                    if (ElecPt->at(i) > subElecPtCut) subElecIndex = i;
                }
                else {
                    continue;
                }
            }
        }
        else {
            continue;
        }
    }
    vector<int> * indices = new vector<int>;
    indices->push_back(leadElecIndex); 
    //    indices->push_back(subElecIndex); 
    indices->push_back(leadPositIndex); 
    //    indices->push_back(subPositIndex);
    return indices;
}


inline vector<int> * MuonPickOvi(vector<float> * MuonPt, vector<float> * MuonEta, vector<int> * MuonCharge, vector<float> * MuonD0, vector<float> *MuonVertZ, vector<float> * VertexZ, vector<float> * MuonNeutHadIso, vector<float> * MuonCharHadIso, vector<float> * MuonPhotHadIso, vector<float> * MuonSumPUPt, vector<bool> * isGMPT, vector<bool> * isPFMuon) {
    int leadMuonIndex = -1; int subMuonIndex = -1; int leadMuBarIndex = -1; int subMuBarIndex = -1;
    float muonIsoRatioCut = 0.15; float muonEtaCut = 2.4; float leadMuonPtCut = 20; float subMuonPtCut = 10;
    float muonD0Cut = 0.2; float muonDZCut = 0.5;
    
    float muonPFIso;
    //    cout << "i max " << MuonPt->size() << endl;
    for (unsigned int i = 0; i < MuonPt->size(); ++i) {
        /*
         cout << "i " << i << endl;
         cout << "isGMPT " << isGMPT->at(i) << endl;
         cout << "isPFMuon " << isPFMuon->at(i) << endl;
         cout << "neuthadiso " << MuonNeutHadIso->at(i) << endl;
         cout << "charhadiso " << MuonCharHadIso->at(i) << endl;
         cout << "photiso " << MuonPhotHadIso->at(i) << endl;
         cout << "MuonPt " << MuonPt->at(i) << endl;
         cout << "MuonEta " << MuonEta->at(i) << endl;
         cout << "VertexZ->at " << VertexZ->at(0) << endl;
         cout << "muonVertexZ->at " << MuonVertZ->at(i) << endl;
         cout << "MuonCharge->at " << MuonCharge->at(i) << endl;
         cout << "abs(VertexZ->at(0) - MuonVertZ->at(i)) " << abs(VertexZ->at(0) - MuonVertZ->at(i)) << endl;
         */
        //         [sumChargedHadronPt+ max(0.,sumNeutralHadronPt+sumPhotonPt-0.5sumPUPt]/pt
        muonPFIso = MuonCharHadIso->at(i) + TMath::Max((float)0., MuonPhotHadIso->at(i) + MuonNeutHadIso->at(i) - MuonSumPUPt->at(i));
        if (isGMPT->at(i) && isPFMuon->at(i) && muonPFIso/MuonPt->at(i) < muonIsoRatioCut && abs(MuonEta->at(i)) < muonEtaCut && MuonD0->at(0) < muonD0Cut && abs(VertexZ->at(0) - MuonVertZ->at(i)) < muonDZCut) {
            if (MuonCharge->at(i) > 0) {
                if (leadMuBarIndex < 0) {
                    if (MuonPt->at(i) > leadMuonPtCut) leadMuBarIndex = i;
                }
                else if (subMuBarIndex < 0) {
                    if (MuonPt->at(i) > subMuonPtCut) subMuBarIndex = i;
                }
                else {
                    continue;
                }
            }
            else {
                if (leadMuonIndex < 0) {
                    if (MuonPt->at(i) > leadMuonPtCut) leadMuonIndex = i;
                }
                else if (subMuonIndex < 0) {
                    if (MuonPt->at(i) > subMuonPtCut) subMuonIndex = i;
                }
                else {
                    continue;
                }
            }
        }
        else {
            continue;
        }
    }
    vector<int> * indices = new vector<int>;
    indices->push_back(leadMuonIndex); 
    //    indices->push_back(subMuonIndex); 
    indices->push_back(leadMuBarIndex); 
    //indices->push_back(subMuBarIndex);
    return indices;
}
inline void LeptonPairOvi(TLorentzVector &Lep0Vec, TLorentzVector &Lep1Vec, int &eventType, bool &doEvent, vector<TLorentzVector> * Leptons, vector<int> * lepPdgId, int &Lep0PdgId,  int &Lep1PdgId) {
    //    cout << "making lepton pair" << endl;
    int lep0Index = -1; int lep1Index = -1;
    //    int currLeadLepIndex, currSubLepIndex;
    //    int currLeadLepPt, currSubLepPt;
    float currDiLepPt; float leadDiLepPt = 0;
    float massCut = 20;
    int productPdgId;
    //Try some cases
    for (unsigned iLep = 0; iLep < Leptons->size(); ++iLep) {
        for (unsigned iLep2 = iLep + 1; iLep2 < Leptons->size(); ++iLep2) {
            if (lepPdgId->at(iLep) * lepPdgId->at(iLep2) > 0) continue;
            if ((Leptons->at(iLep) + Leptons->at(iLep2)).M() < massCut) continue;
            currDiLepPt = Leptons->at(iLep).Pt() + Leptons->at(iLep2).Pt();
            if (currDiLepPt > leadDiLepPt) {
                leadDiLepPt = currDiLepPt;
                if (Leptons->at(iLep).Pt() > Leptons->at(iLep2).Pt()) {
                    lep0Index = iLep;
                    lep1Index = iLep2;
                }
                else {
                    lep0Index = iLep2;
                    lep1Index = iLep;
                }
            }
        }
    }
    if (lep0Index == -1 || lep1Index == -1) {
        doEvent = false;
        return;
    }
    Lep0Vec.SetPxPyPzE(Leptons->at(lep0Index).Px(), Leptons->at(lep0Index).Py(), Leptons->at(lep0Index).Pz(), Leptons->at(lep0Index).E());
    Lep1Vec.SetPxPyPzE(Leptons->at(lep1Index).Px(), Leptons->at(lep1Index).Py(), Leptons->at(lep1Index).Pz(), Leptons->at(lep1Index).E());
    Lep0PdgId = lepPdgId->at(lep0Index);
    Lep1PdgId = lepPdgId->at(lep1Index);
    productPdgId = Lep0PdgId * Lep1PdgId;
    switch (productPdgId) {
        case -169:
            eventType = 0;
            break;
        case -121:
            eventType = 1;
            break;
        case -143:
            eventType = 2;
            break;
        default:
            eventType = -1;
            cout << "something funky with event Type!!" << endl;
            cout << "Lep0PdgId " << Lep0PdgId << endl;
            cout << "Lep1PdgId " << Lep1PdgId << endl;
            break;
    }    
}

inline vector<TLorentzVector> * JetInfo(vector<TLorentzVector> * Leptons, vector<float> * JetEt, vector<float> * JetEta, vector<float> * JetPx, vector<float> * JetPy, vector<float> * JetPz, vector<float> * JetE, vector<float> *JetNHF, vector<float> * JetNEF, vector<float> * JetCHF, vector<float> * JetCEF, vector<float> * JetBTag, vector<int> * JetNDaug, vector<int> * JetCharMult, int &NJets, int &NBJets, vector<int> * BJetIndices, float &HT) {
    NJets = 0;
    NBJets = 0;
    HT = 0;
    vector<TLorentzVector> * Jets = new vector<TLorentzVector>;
    TLorentzVector currJet, patsyJet;
    float BTagWP = 0.679;  //CSV Middle working point, see (remove underscore in address): h_ttps://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    bool leptonJet;
    float dRcut = 0.5;
    //    cout << "iJet Max " << JetEt->size() << endl;
    for (unsigned int iJet = 0; iJet < JetEt->size(); ++iJet) {
        leptonJet = 0;
        /*
         cout << "iJet " << iJet << endl;        
         cout << "JetEta " << JetEta->at(iJet) << endl;
         cout << "JetEt " << JetEt->at(iJet) << endl;
         cout << "JetNDaug " << JetNDaug->at(iJet) << endl;
         cout << "JetNHF " << JetNHF->at(iJet) << endl;
         cout << "JetNEF " << JetNEF->at(iJet) << endl;
         cout << "JetCHF " << JetCHF->at(iJet) << endl;
         cout << "JetCEF " << JetCEF->at(iJet) << endl;
         cout << "JetCharMult " << JetCharMult->at(iJet) << endl;
         cout << "JetBTag " << JetBTag->at(iJet) << endl;
         */
        //        cout << "JetPt " << TMath::Sqrt(JetPx->at(iJet)*JetPx->at(iJet) + JetPy->at(iJet)*JetPy->at(iJet)) << endl;
        patsyJet.SetPxPyPzE(JetPx->at(iJet), JetPy->at(iJet), JetPz->at(iJet), JetE->at(iJet));
        for (unsigned int iLep = 0; iLep < Leptons->size(); ++iLep) {
            if (deltaR(patsyJet.Eta(), patsyJet.Phi(), Leptons->at(iLep).Eta(), Leptons->at(iLep).Phi()) < dRcut) leptonJet = 1;
        }
        if (leptonJet) continue;
        if (abs(JetEta->at(iJet)) > 2.4) continue;
        if (JetEt->at(iJet) < 30) continue;
        //        if (JetNDaug->at(iJet) > 1 && JetNHF->at(iJet) < 0.99 && JetNEF->at(iJet) < 0.99 && JetCEF->at(iJet) < 0.99 && JetCHF->at(iJet) > 0 && JetCharMult->at(iJet) > 0) {
        if (JetNHF->at(iJet) < 0.99 && JetNEF->at(iJet) < 0.99 && JetCEF->at(iJet) < 0.99 && JetCHF->at(iJet) > 0 && JetCharMult->at(iJet) > 0) {
            currJet.SetPxPyPzE(JetPx->at(iJet), JetPy->at(iJet), JetPz->at(iJet), JetE->at(iJet));
            ++NJets;
            HT += JetEt->at(iJet);
            if (JetBTag->at(iJet) > BTagWP) {
                NBJets += 1;
                //                cout << "BTAGG is now " << NBJets << " ocurring at " << iJet << endl;
                BJetIndices->push_back(NJets - 1);
            }
            //            cout << "Jets size pre " << Jets->size() << endl;
            Jets->push_back(currJet);
            /*
             for (int iLep = 0; iLep < Leptons->size(); ++iLep) {
             cout << "deltaR for iLep " << iLep << " is " << deltaR(currJet.Eta(), currJet.Phi(), Leptons->at(iLep).Eta(), Leptons->at(iLep).Phi()) << endl;                
             }
             */
            //            cout << "Jets size post " << Jets->size() << endl;
        }
        else {
            continue;
        }
    }
    return Jets;
}
inline void LeptonInfo(VLV * Leptons, vector<int> *lepPdgId, vector<double> *lepPFIso, int &lep0Index, int &lep1Index, bool &doEvent, int &eventType, int &Lep0PdgId, int &Lep1PdgId) {
    doEvent = true;
    lep0Index = 0;
    lep1Index = 0;
    int vecSize = Leptons->size();
    //    cout << "vecSize " << vecSize << endl;
    if (vecSize < 2) {
        doEvent = false;
        return;
    }
    float muonIsoRatioCut = 0.15; float muonEtaCut = 2.4; float leadMuonPtCut = 20; float subMuonPtCut = 10;
    float elecIsoRatioCut = 0.15; float elecEtaCut = 2.5; float leadElecPtCut = 20; float subElecPtCut = 10;
    float barrelEtaEnd = 1.442; float endcapEtaStart = 1.566;
    float currLeadLepPt, currSubLepPt;
    float currLeadLepEta, currSubLepEta;
    float currLeadLepCharge, currSubLepCharge;
    int productPdgId;
    float currDiLepPt;
    float leadDiLepPt = 0;
    float massCut = 20;
    
    //    cout << "" << endl;
    //    if (vecSize > 2) cout << "vecSize: " << vecSize << endl;
    for (int iLep = 0; iLep < vecSize; ++iLep) {
        currLeadLepPt = Leptons->at(iLep).Pt();
        currLeadLepEta = Leptons->at(iLep).Eta();
        currLeadLepCharge = TMath::Sign(1, lepPdgId->at(iLep));
        //        cout << "currLeadLepPt for iLep = " << iLep << " is " << currLeadLepPt << endl;
        //        cout << "lead lepPDGId" << lepPdgId->at(iLep) << endl;
        if (abs(lepPdgId->at(iLep)) == 13) {
            //            cout << "IsoCut: " << lepPFIso->at(iLep)/currLeadLepPt << endl;
            if (currLeadLepPt < leadMuonPtCut) continue;
            if (lepPFIso->at(iLep)/currLeadLepPt > muonIsoRatioCut) continue;
            if (abs(currLeadLepEta) > muonEtaCut) continue;
        }
        else if (abs(lepPdgId->at(iLep)) == 11) {
            if (currLeadLepPt < leadElecPtCut) continue;
            if (lepPFIso->at(iLep)/currLeadLepPt > elecIsoRatioCut) continue;
            if (abs(currLeadLepEta) > elecEtaCut) continue;
            if (abs(currLeadLepEta) > barrelEtaEnd && abs(currLeadLepEta) < endcapEtaStart) continue;
        }
        for (int iLep2 = iLep+1; iLep2 < vecSize; ++iLep2) {
            currSubLepPt = Leptons->at(iLep2).Pt();
            currSubLepEta = Leptons->at(iLep2).Eta();
            currSubLepCharge = TMath::Sign(1, lepPdgId->at(iLep2));
            if (currLeadLepCharge == currSubLepCharge) continue;    //// opposite sign charge requirement
            //            cout << "currSubLepPt for iLep = " << iLep2 << " is " << currSubLepPt << endl;
            //            cout << "sub lead lepPDGId" << lepPdgId->at(iLep2) << endl;
            if (abs(lepPdgId->at(iLep2)) == 13) {
                if (currSubLepPt < subMuonPtCut) continue;
                if (lepPFIso->at(iLep2)/currSubLepPt > muonIsoRatioCut) continue;
                if (abs(currSubLepEta) > muonEtaCut) continue;
            }
            else if (abs(lepPdgId->at(iLep2)) == 11) {
                if (currSubLepPt < subElecPtCut) continue;
                if (lepPFIso->at(iLep2)/currSubLepPt > elecIsoRatioCut) continue;
                if (abs(currSubLepEta) > elecEtaCut) continue;
                if (abs(currSubLepEta) > barrelEtaEnd && abs(currSubLepEta) < endcapEtaStart) continue;
            }
            if ((Leptons->at(iLep) + Leptons->at(iLep2)).M() < massCut) continue;
            productPdgId = lepPdgId->at(iLep) * lepPdgId->at(iLep2);
            if (productPdgId > 0) continue; //i.e. same sign event
            currDiLepPt = currLeadLepPt + currSubLepPt;
            //            cout << "currDiLep " << currDiLepPt << endl;
            if (currDiLepPt > leadDiLepPt) {
                leadDiLepPt = currDiLepPt;
                lep0Index = iLep;
                lep1Index = iLep2;
            }
        }
    }
    //    cout << "lep0 Index " << lep0Index << endl;
    //    cout << "lep1 Index " << lep1Index << endl;
    //    cout << "" << endl;
    if (lep0Index == lep1Index) {
        doEvent = false;
        /*
         cout << "vecSize: " << vecSize << endl;
         for (int i = 0; i < vecSize; ++i) {
         cout << "i: " << i << endl;
         cout << "lepPDGId: " << lepPdgId->at(i) << endl;
         cout << "lepPt: " << Leptons->at(i).Pt() << endl;
         cout << "lepEta: " << Leptons->at(i).Eta() << endl;
         if (vecSize > 1) cout << "diLepMass: " << (Leptons->at(0) + Leptons->at(1)).M() << endl;
         }
         cout << "no suitable pair found! " << endl;
         cout << "" << endl;
         */
    }
    Lep0PdgId = lepPdgId->at(lep0Index);
    Lep1PdgId = lepPdgId->at(lep1Index);
    productPdgId = Lep0PdgId * Lep1PdgId;
    if (productPdgId == -169) eventType = 0;
    if (productPdgId == -121) eventType = 1;
    if (productPdgId == -143) eventType = 2;
}
inline vector<int> * NumJets(VLV * Jets, float jetPtCut, VLV * Leptons, int lep0Index, int lep1Index,  vector<double> *jetBTagParam, float &HT, float bTagCut) {
    LV CurrentJet;
    LV Lepton0 = Leptons->at(lep0Index);
    LV Lepton1 = Leptons->at(lep1Index);
    int NJets = 0;
    int NBJets = 0;
    HT = 0;
    int jet0Index = -99;
    int jet1Index = -99;
    float CurrentLeadJetPt = 0;
    float CurrentSubJetPt = 0;
    int BJet0Index = -99;
    int BJet1Index = -99;
    float dRcut = 0.5;
    vector<int> * eventJetParams = new vector<int>;
    for (int iJet = 0; iJet < (int) Jets->size(); ++iJet) {        
        CurrentJet = Jets->at(iJet);
        if (CurrentJet.Pt() > jetPtCut) {
            if (deltaR(CurrentJet.Eta(), CurrentJet.Phi(), Lepton0.Eta(), Lepton0.Phi()) < dRcut) {
                //                cout << "delta R failed for lep 0 " << deltaR(CurrentJet.Eta(), CurrentJet.Phi(), Lepton0.Eta(), Lepton0.Phi()) << endl;
                continue; 
            }
            if (deltaR(CurrentJet.Eta(), CurrentJet.Phi(), Lepton1.Eta(), Lepton1.Phi()) < dRcut) continue;
            HT += CurrentJet.Pt();
            ++NJets;
            if (jet0Index < 0 && CurrentJet.Pt() > CurrentLeadJetPt) {
                CurrentLeadJetPt = CurrentJet.Pt();
                jet0Index = iJet;
            }
            else if (jet0Index >= 0 && jet1Index < 0 && CurrentJet.Pt() > CurrentSubJetPt) {
                CurrentSubJetPt = CurrentJet.Pt();
                jet1Index = iJet;
            }
            if (jetBTagParam->at(iJet) > bTagCut) {
                ++NBJets;
                if (BJet0Index < 0) {
                    BJet0Index = iJet;
                }
                else {
                    if (BJet1Index < 0) BJet1Index = iJet;
                }   
            }   
        }
    }
    eventJetParams->push_back(NJets);
    eventJetParams->push_back(NBJets);
    eventJetParams->push_back(jet0Index);
    eventJetParams->push_back(jet1Index);
    eventJetParams->push_back(BJet0Index);
    eventJetParams->push_back(BJet1Index);
    return eventJetParams;
}

inline vector<HistogramT> * AddSystHists(vector<HistogramT> * inputHistTVec, vector<SystT> * inputSystTVec, TString fileInName, bool isSignal) {
    HistogramT H_Curr;
    HistogramT H_CurrNewSyst;
    SystT      Syst_Curr;
    vector<HistogramT> * histTSystVec = new vector<HistogramT>;
    for (unsigned int i = 0; i < inputHistTVec->size(); ++i) {
        H_Curr = inputHistTVec->at(i);
        if (H_Curr.doXSyst) {
            for (unsigned int j = 0; j < inputSystTVec->size(); ++j) {
                Syst_Curr = inputSystTVec->at(j);
                
                switch (Syst_Curr.whichSystType) {
                        /*
                         case 1:
                         if (!(H_Curr.name.Contains("Lep"))) continue;
                         break;
                         case 2:
                         if (!(H_Curr.name.Contains("Jet") || (H_Curr.name.Contains("HT")))) continue;
                         break;
                         */
                    case 3:
                        // temporary one for MT2ll
                        if (Syst_Curr.name.Contains("_MT2llShift")) {
                            if (!(H_Curr.name.Contains("MT2ll"))) continue;
                        }
                        if (Syst_Curr.name.Contains("genTopRW")) {
                            if (!(fileInName.Contains("TT") || fileInName.Contains("ttbar"))) continue;
                        }
                        if (Syst_Curr.name.Contains("genStopXSec")) {
                            if (!isSignal) continue;
                        }
                        break;
                    default:
                        break;
                }
                
                H_CurrNewSyst = H_Curr;
                H_CurrNewSyst.name += Syst_Curr.name;
                H_CurrNewSyst.xVarKey += Syst_Curr.systVarKey;
                histTSystVec->push_back(H_CurrNewSyst);
            }
        }
    }
    return histTSystVec;
}

inline vector<HistogramT> * OneDeeHistTVec() {
    const double PI = 3.14159265;
    int EnergyPtBinN = 50;
    int MassBinN     = 200;
    int EtaBinN      = 50;
    int PhiBinN      = 50;
    int METBinN      = 40;    
    int METXYBinN    = 50;    
    int NJetsBinN    = 11;
    int nVtxBinN     = 60;
    
    float EnergyPtBinLB = 0;
    float EnergyPtBinUB = 400;
    float EtaBinLB      = -6;
    float EtaBinUB      = 6;
    float METBinLB      = 0;
    float METBinUB      = 400;
    float METXYBinLB    = -200;
    float METXYBinUB    = 200;
    float NJetsBinLB    = -0.5;
    float NJetsBinUB    = 10.5;
    float nVtxBinLB     = 0.5;
    float nVtxBinUB     = 60.5;
    
    float numDivs;
    
    HistogramT H_leadLepPt; H_leadLepPt.name = "h_leadLepPt"; 
    H_leadLepPt.xLabel = "Lead Lepton pT [GeV]"; H_leadLepPt.xBinN = EnergyPtBinN; H_leadLepPt.xMin = EnergyPtBinLB; H_leadLepPt.xMax = EnergyPtBinUB;  
    H_leadLepPt.yLabel = "Number of Events / ";
    numDivs = (H_leadLepPt.xMax - H_leadLepPt.xMin) / (float) H_leadLepPt.xBinN;
    H_leadLepPt.yLabel += "NUM"; H_leadLepPt.yLabel += " GeV";
    H_leadLepPt.xVarKey = "leadLepPt";
    H_leadLepPt.doXSyst = true;
    
    HistogramT H_leadLepEta; H_leadLepEta.name = "h_leadLepEta"; 
    H_leadLepEta.xLabel = "Lead Lepton #eta"; H_leadLepEta.xBinN = EtaBinN; H_leadLepEta.xMin = EtaBinLB; H_leadLepEta.xMax = EtaBinUB; 
    H_leadLepEta.yLabel = "Number of Events / ";
    numDivs = (H_leadLepEta.xMax - H_leadLepEta.xMin) / (float) H_leadLepEta.xBinN;
    H_leadLepEta.yLabel += "NUM"; H_leadLepEta.yLabel += " #eta";
    H_leadLepEta.xVarKey = "leadLepEta";
    H_leadLepEta.doXSyst = false;
    
    HistogramT H_subLepPt; H_subLepPt.name = "h_subLepPt"; 
    H_subLepPt.xLabel = "Sub-lead Lepton pT [GeV]"; H_subLepPt.xBinN = EnergyPtBinN; H_subLepPt.xMin = EnergyPtBinLB; H_subLepPt.xMax = EnergyPtBinUB;  
    H_subLepPt.yLabel = "Number of Events / ";
    numDivs = (H_subLepPt.xMax - H_subLepPt.xMin) / (float) H_subLepPt.xBinN;
    H_subLepPt.yLabel += "NUM"; H_subLepPt.yLabel += " GeV";
    H_subLepPt.xVarKey = "subLepPt";
    H_subLepPt.doXSyst = true;
    
    HistogramT H_subLepEta; H_subLepEta.name = "h_subLepEta"; 
    H_subLepEta.xLabel = "Sub-lead Lepton #eta"; H_subLepEta.xBinN = EtaBinN; H_subLepEta.xMin = EtaBinLB; H_subLepEta.xMax = EtaBinUB; 
    H_subLepEta.yLabel = "Number of Events / ";
    numDivs = (H_subLepEta.xMax - H_subLepEta.xMin) / (float) H_subLepEta.xBinN;
    H_subLepEta.yLabel += "NUM"; H_subLepEta.yLabel += " #eta";
    H_subLepEta.xVarKey = "subLepEta";
    H_subLepEta.doXSyst = false;
    
    HistogramT H_leadJetPt; H_leadJetPt.name = "h_leadJetPt"; 
    H_leadJetPt.xLabel = "Lead Jet pT [GeV]"; H_leadJetPt.xBinN = EnergyPtBinN; H_leadJetPt.xMin = EnergyPtBinLB; H_leadJetPt.xMax = EnergyPtBinUB;  
    H_leadJetPt.yLabel = "Number of Events / ";
    numDivs = (H_leadJetPt.xMax - H_leadJetPt.xMin) / (float) H_leadJetPt.xBinN;
    H_leadJetPt.yLabel += "NUM"; H_leadJetPt.yLabel += " GeV";
    H_leadJetPt.xVarKey = "leadJetPt";
    H_leadJetPt.doXSyst = true;
    
    HistogramT H_leadJetEta; H_leadJetEta.name = "h_leadJetEta"; 
    H_leadJetEta.xLabel = "Lead Jet #eta"; H_leadJetEta.xBinN = EtaBinN; H_leadJetEta.xMin = EtaBinLB; H_leadJetEta.xMax = EtaBinUB; 
    H_leadJetEta.yLabel = "Number of Events / ";
    numDivs = (H_leadJetEta.xMax - H_leadJetEta.xMin) / (float) H_leadJetEta.xBinN;
    H_leadJetEta.yLabel += "NUM"; H_leadJetEta.yLabel += " #eta";
    H_leadJetEta.xVarKey = "leadJetEta";
    H_leadJetEta.doXSyst = false;
    
    HistogramT H_subJetPt; H_subJetPt.name = "h_subJetPt"; 
    H_subJetPt.xLabel = "Sub-lead Jet pT [GeV]"; H_subJetPt.xBinN = EnergyPtBinN; H_subJetPt.xMin = EnergyPtBinLB; H_subJetPt.xMax = EnergyPtBinUB;  
    H_subJetPt.yLabel = "Number of Events / ";
    numDivs = (H_subJetPt.xMax - H_subJetPt.xMin) / (float) H_subJetPt.xBinN;
    H_subJetPt.yLabel += "NUM"; H_subJetPt.yLabel += " GeV";
    H_subJetPt.xVarKey = "subJetPt";
    H_subJetPt.doXSyst = true;
    
    HistogramT H_subJetEta; H_subJetEta.name = "h_subJetEta"; 
    H_subJetEta.xLabel = "Sub-lead Jet #eta"; H_subJetEta.xBinN = EtaBinN; H_subJetEta.xMin = EtaBinLB; H_subJetEta.xMax = EtaBinUB; 
    H_subJetEta.yLabel = "Number of Events / ";
    numDivs = (H_subJetEta.xMax - H_subJetEta.xMin) / (float) H_subJetEta.xBinN;
    H_subJetEta.yLabel += "NUM"; H_subJetEta.yLabel += " #eta";
    H_subJetEta.xVarKey = "subJetEta";
    H_subJetEta.doXSyst = false;
    
    HistogramT H_leadBJetPt; H_leadBJetPt.name = "h_leadBJetPt"; 
    H_leadBJetPt.xLabel = "Lead BJet pT [GeV]"; H_leadBJetPt.xBinN = EnergyPtBinN; H_leadBJetPt.xMin = EnergyPtBinLB; H_leadBJetPt.xMax = EnergyPtBinUB;  
    H_leadBJetPt.yLabel = "Number of Events / ";
    numDivs = (H_leadBJetPt.xMax - H_leadBJetPt.xMin) / (float) H_leadBJetPt.xBinN;
    H_leadBJetPt.yLabel += "NUM"; H_leadBJetPt.yLabel += " GeV";
    H_leadBJetPt.xVarKey = "leadBJetPt";
    H_leadBJetPt.doXSyst = true;
    
    HistogramT H_leadBJetEta; H_leadBJetEta.name = "h_leadBJetEta"; 
    H_leadBJetEta.xLabel = "Lead BJet #eta"; H_leadBJetEta.xBinN = EtaBinN; H_leadBJetEta.xMin = EtaBinLB; H_leadBJetEta.xMax = EtaBinUB; 
    H_leadBJetEta.yLabel = "Number of Events / ";
    numDivs = (H_leadBJetEta.xMax - H_leadBJetEta.xMin) / (float) H_leadBJetEta.xBinN;
    H_leadBJetEta.yLabel += "NUM"; H_leadBJetEta.yLabel += " #eta";
    H_leadBJetEta.xVarKey = "leadBJetEta";
    H_leadBJetEta.doXSyst = false;
    
    HistogramT H_leadBJetEn; H_leadBJetEn.name = "h_leadBJetEn"; 
    H_leadBJetEn.xLabel = "Lead BJet Energy"; H_leadBJetEn.xBinN = 3*EnergyPtBinN; H_leadBJetEn.xMin = EnergyPtBinLB; H_leadBJetEn.xMax = 3*EnergyPtBinUB; 
    H_leadBJetEn.yLabel = "Number of Events / ";
    numDivs = (H_leadBJetEn.xMax - H_leadBJetEn.xMin) / (float) H_leadBJetEn.xBinN;
    H_leadBJetEn.yLabel += "NUM"; H_leadBJetEn.yLabel += " GeV";
    H_leadBJetEn.xVarKey = "leadBJetEn";
    H_leadBJetEn.doXSyst = true;
    
    HistogramT H_subBJetPt; H_subBJetPt.name = "h_subBJetPt"; 
    H_subBJetPt.xLabel = "Sub-lead BJet pT [GeV]"; H_subBJetPt.xBinN = EnergyPtBinN; H_subBJetPt.xMin = EnergyPtBinLB; H_subBJetPt.xMax = EnergyPtBinUB;  
    H_subBJetPt.yLabel = "Number of Events / ";
    numDivs = (H_subBJetPt.xMax - H_subBJetPt.xMin) / (float) H_subBJetPt.xBinN;
    H_subBJetPt.yLabel += "NUM"; H_subBJetPt.yLabel += " GeV";
    H_subBJetPt.xVarKey = "subBJetPt";
    H_subBJetPt.doXSyst = true;
    
    HistogramT H_subBJetEta; H_subBJetEta.name = "h_subBJetEta"; 
    H_subBJetEta.xLabel = "Sub-lead BJet #eta"; H_subBJetEta.xBinN = EtaBinN; H_subBJetEta.xMin = EtaBinLB; H_subBJetEta.xMax = EtaBinUB; 
    H_subBJetEta.yLabel = "Number of Events / ";
    numDivs = (H_subBJetEta.xMax - H_subBJetEta.xMin) / (float) H_subBJetEta.xBinN;
    H_subBJetEta.yLabel += "NUM"; H_subBJetEta.yLabel += " #eta";
    H_subBJetEta.xVarKey = "subBJetEta";
    H_subBJetEta.doXSyst = false;
    
    HistogramT H_subBJetEn; H_subBJetEn.name = "h_subBJetEn"; 
    H_subBJetEn.xLabel = "Sub-lead BJet Energy"; H_subBJetEn.xBinN = EnergyPtBinN; H_subBJetEn.xMin = EnergyPtBinLB; H_subBJetEn.xMax = 3*EnergyPtBinUB; 
    H_subBJetEn.yLabel = "Number of Events / ";
    numDivs = (H_subBJetEn.xMax - H_subBJetEn.xMin) / (float) H_subBJetEn.xBinN;
    H_subBJetEn.yLabel += "NUM"; H_subBJetEn.yLabel += " GeV";
    H_subBJetEn.xVarKey = "subBJetEn";
    H_subBJetEn.doXSyst = true;
    
    HistogramT H_diLepPt; H_diLepPt.name = "h_DiLepPt"; 
    H_diLepPt.xLabel = "diLep pT [GeV]"; H_diLepPt.xBinN = EnergyPtBinN; H_diLepPt.xMin = EnergyPtBinLB; H_diLepPt.xMax = EnergyPtBinUB;  
    H_diLepPt.yLabel = "Number of Events / ";
    numDivs = (H_diLepPt.xMax - H_diLepPt.xMin) / (float) H_diLepPt.xBinN;
    H_diLepPt.yLabel += "NUM"; H_diLepPt.yLabel += " GeV";
    H_diLepPt.xVarKey = "diLepPt";
    H_diLepPt.doXSyst = true;
    
    HistogramT H_diLepInvMass; H_diLepInvMass.name = "h_DiLepInvMass"; 
    H_diLepInvMass.xLabel = "M_{ll} [GeV]"; H_diLepInvMass.xBinN = MassBinN; H_diLepInvMass.xMin = EnergyPtBinLB; H_diLepInvMass.xMax = EnergyPtBinUB; 
    H_diLepInvMass.yLabel = "Number of Events / ";
    numDivs = (H_diLepInvMass.xMax - H_diLepInvMass.xMin) / (float) H_diLepInvMass.xBinN;
    H_diLepInvMass.yLabel += "NUM"; H_diLepInvMass.yLabel += " GeV";
    H_diLepInvMass.xVarKey = "diLepInvMass";
    H_diLepInvMass.doXSyst = true;
    
    HistogramT H_diLepEta; H_diLepEta.name = "h_DiLepEta"; 
    H_diLepEta.xLabel = "#eta"; H_diLepEta.xBinN = EtaBinN; H_diLepEta.xMin = EtaBinLB; H_diLepEta.xMax = EtaBinUB; 
    H_diLepEta.yLabel = "Number of Events / ";
    numDivs = (H_diLepEta.xMax - H_diLepEta.xMin) / (float) H_diLepEta.xBinN;
    H_diLepEta.yLabel += "NUM"; H_diLepEta.yLabel += " #eta";
    H_diLepEta.xVarKey = "diLepEta";
    H_diLepEta.doXSyst = true;
    
    HistogramT H_diLepPhi; H_diLepPhi.name = "h_DiLepPhi"; 
    H_diLepPhi.xLabel = "#phi_{ll}"; H_diLepPhi.xBinN = PhiBinN; H_diLepPhi.xMin = -PI; H_diLepPhi.xMax = PI; 
    H_diLepPhi.yLabel = "Number of Events / ";
    numDivs = (H_diLepPhi.xMax - H_diLepPhi.xMin) / (float) H_diLepPhi.xBinN;
    H_diLepPhi.yLabel += "NUM"; H_diLepPhi.yLabel += " radians";
    H_diLepPhi.xVarKey = "diLepPhi";
    H_diLepPhi.doXSyst = true;
    
    ///Met or MET related variable plots
    
    HistogramT H_MT2ll; H_MT2ll.name = "h_MT2ll"; 
    H_MT2ll.xLabel = "MT2_{ll} [GeV]"; H_MT2ll.xBinN = METBinN; H_MT2ll.xMin = METBinLB; H_MT2ll.xMax = METBinUB; 
    H_MT2ll.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll.xMax - H_MT2ll.xMin) / (float) H_MT2ll.xBinN;
    H_MT2ll.yLabel += "NUM"; H_MT2ll.yLabel += " GeV";
    H_MT2ll.xVarKey = "MT2ll";
    H_MT2ll.doXSyst = true;
    
    HistogramT H_MT2llCont; H_MT2llCont.name = "h_MT2llControl"; 
    H_MT2llCont.xLabel = "MT2_{ll} [GeV]"; H_MT2llCont.xBinN = 20; H_MT2llCont.xMin = 0; H_MT2llCont.xMax = 80; 
    H_MT2llCont.yLabel = "Number of Events / ";
    numDivs = (H_MT2llCont.xMax - H_MT2llCont.xMin) / (float) H_MT2llCont.xBinN;
    H_MT2llCont.yLabel += "NUM"; H_MT2llCont.yLabel += " GeV";
    H_MT2llCont.xVarKey = "MT2ll";
    H_MT2llCont.doXSyst = true;
    
    HistogramT H_PassMT2llCut80; H_PassMT2llCut80.name = "h_PassMT2llCut80"; 
    H_PassMT2llCut80.xLabel = "Event MT2_{ll} > 80 [GeV]"; H_PassMT2llCut80.xBinN = 2; H_PassMT2llCut80.xMin = -0.5; H_PassMT2llCut80.xMax = 1.5; 
    H_PassMT2llCut80.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut80.xVarKey = "PassMT2llCut80";
    H_PassMT2llCut80.doXSyst = true;
    
    HistogramT H_PassMT2llCut90; H_PassMT2llCut90.name = "h_PassMT2llCut90"; 
    H_PassMT2llCut90.xLabel = "Event MT2_{ll} > 90 [GeV]"; H_PassMT2llCut90.xBinN = 2; H_PassMT2llCut90.xMin = -0.5; H_PassMT2llCut90.xMax = 1.5; 
    H_PassMT2llCut90.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut90.xVarKey = "PassMT2llCut90";
    H_PassMT2llCut90.doXSyst = true;
    
    HistogramT H_PassMT2llCut100; H_PassMT2llCut100.name = "h_PassMT2llCut100"; 
    H_PassMT2llCut100.xLabel = "Event MT2_{ll} > 100 [GeV]"; H_PassMT2llCut100.xBinN = 2; H_PassMT2llCut100.xMin = -0.5; H_PassMT2llCut100.xMax = 1.5; 
    H_PassMT2llCut100.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut100.xVarKey = "PassMT2llCut100";
    H_PassMT2llCut100.doXSyst = true;
    
    HistogramT H_PassMT2llCut110; H_PassMT2llCut110.name = "h_PassMT2llCut110"; 
    H_PassMT2llCut110.xLabel = "Event MT2_{ll} > 110 [GeV]"; H_PassMT2llCut110.xBinN = 2; H_PassMT2llCut110.xMin = -0.5; H_PassMT2llCut110.xMax = 1.5; 
    H_PassMT2llCut110.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut110.xVarKey = "PassMT2llCut110";
    H_PassMT2llCut110.doXSyst = true;
    
    HistogramT H_PassMT2llCut120; H_PassMT2llCut120.name = "h_PassMT2llCut120"; 
    H_PassMT2llCut120.xLabel = "Event MT2_{ll} > 120 [GeV]"; H_PassMT2llCut120.xBinN = 2; H_PassMT2llCut120.xMin = -0.5; H_PassMT2llCut120.xMax = 1.5; 
    H_PassMT2llCut120.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut120.xVarKey = "PassMT2llCut120";
    H_PassMT2llCut120.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiZMETClose; H_MT2ll_DPhiZMETClose.name = "h_MT2ll_DPhiZMETClose"; 
    H_MT2ll_DPhiZMETClose.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETClose.xBinN = METBinN; H_MT2ll_DPhiZMETClose.xMin = METBinLB; H_MT2ll_DPhiZMETClose.xMax = METBinUB; 
    H_MT2ll_DPhiZMETClose.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETClose.xMax - H_MT2ll_DPhiZMETClose.xMin) / (float) H_MT2ll_DPhiZMETClose.xBinN;
    H_MT2ll_DPhiZMETClose.yLabel += "NUM"; H_MT2ll_DPhiZMETClose.yLabel += " GeV";
    H_MT2ll_DPhiZMETClose.xVarKey = "MT2ll";
    H_MT2ll_DPhiZMETClose.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiZMETMid; H_MT2ll_DPhiZMETMid.name = "h_MT2ll_DPhiZMETMid"; 
    H_MT2ll_DPhiZMETMid.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETMid.xBinN = METBinN; H_MT2ll_DPhiZMETMid.xMin = METBinLB; H_MT2ll_DPhiZMETMid.xMax = METBinUB; 
    H_MT2ll_DPhiZMETMid.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETMid.xMax - H_MT2ll_DPhiZMETMid.xMin) / (float) H_MT2ll_DPhiZMETMid.xBinN;
    H_MT2ll_DPhiZMETMid.yLabel += "NUM"; H_MT2ll_DPhiZMETMid.yLabel += " GeV";
    H_MT2ll_DPhiZMETMid.xVarKey = "MT2ll";
    H_MT2ll_DPhiZMETMid.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiZMETFar; H_MT2ll_DPhiZMETFar.name = "h_MT2ll_DPhiZMETFar"; 
    H_MT2ll_DPhiZMETFar.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETFar.xBinN = METBinN; H_MT2ll_DPhiZMETFar.xMin = METBinLB; H_MT2ll_DPhiZMETFar.xMax = METBinUB; 
    H_MT2ll_DPhiZMETFar.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETFar.xMax - H_MT2ll_DPhiZMETFar.xMin) / (float) H_MT2ll_DPhiZMETFar.xBinN;
    H_MT2ll_DPhiZMETFar.yLabel += "NUM"; H_MT2ll_DPhiZMETFar.yLabel += " GeV";
    H_MT2ll_DPhiZMETFar.xVarKey = "MT2ll";
    H_MT2ll_DPhiZMETFar.doXSyst = true;
    
    HistogramT H_MT2lb; H_MT2lb.name = "h_MT2lb"; 
    H_MT2lb.xLabel = "MT2lb [GeV]"; H_MT2lb.xBinN = METBinN; H_MT2lb.xMin = METBinLB; H_MT2lb.xMax = METBinUB; 
    H_MT2lb.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb.xMax - H_MT2lb.xMin) / (float) H_MT2lb.xBinN;
    H_MT2lb.yLabel += "NUM"; H_MT2lb.yLabel += " GeV";
    H_MT2lb.xVarKey = "MT2lb";
    H_MT2lb.doXSyst = true;
    
    HistogramT H_MT2lbCont; H_MT2lbCont.name = "h_MT2lbControl"; 
    H_MT2lbCont.xLabel = "MT2lb [GeV]"; H_MT2lbCont.xBinN = 43; H_MT2lbCont.xMin = 0; H_MT2lbCont.xMax = 172; 
    H_MT2lbCont.yLabel = "Number of Events / ";
    numDivs = (H_MT2lbCont.xMax - H_MT2lbCont.xMin) / (float) H_MT2lbCont.xBinN;
    H_MT2lbCont.yLabel += "NUM"; H_MT2lbCont.yLabel += " GeV";
    H_MT2lbCont.xVarKey = "MT2lb";
    H_MT2lbCont.doXSyst = true;
    
    HistogramT H_MET; H_MET.name = "h_MET"; 
    H_MET.xLabel = "#slash{E}_{T} [GeV]"; H_MET.xBinN = METBinN; H_MET.xMin = METBinLB; H_MET.xMax = METBinUB; 
    H_MET.yLabel = "Number of Events / ";
    numDivs = (H_MET.xMax - H_MET.xMin) / (float) H_MET.xBinN;
    H_MET.yLabel += "NUM"; H_MET.yLabel += " GeV";
    H_MET.xVarKey = "MET";
    H_MET.doXSyst = true;
    
    HistogramT H_METX; H_METX.name = "h_METX";
    H_METX.xLabel = "#slash{E}_{x} [GeV]"; H_METX.xBinN = METXYBinN; H_METX.xMin = METXYBinLB; H_METX.xMax = METXYBinUB;
    H_METX.yLabel = "Number of Events / ";
    numDivs = (H_METX.xMax - H_METX.xMin) / (float) H_METX.xBinN;
    H_METX.yLabel += "NUM"; H_METX.yLabel += " GeV";
    H_METX.xVarKey = "METX";
    H_METX.doXSyst = true;
    
    HistogramT H_METY; H_METY.name = "h_METY";
    H_METY.xLabel = "#slash{E}_{y} [GeV]"; H_METY.xBinN = METXYBinN; H_METY.xMin = METXYBinLB; H_METY.xMax = METXYBinUB;
    H_METY.yLabel = "Number of Events / ";
    numDivs = (H_METY.xMax - H_METY.xMin) / (float) H_METY.xBinN;
    H_METY.yLabel += "NUM"; H_METY.yLabel += " GeV";
    H_METY.xVarKey = "METY";
    H_METY.doXSyst = true;
    
    HistogramT H_METPhi; H_METPhi.name = "h_METPhi"; 
    H_METPhi.xLabel = "#slash{E}_{T} #phi"; H_METPhi.xBinN = PhiBinN; H_METPhi.xMin = -PI; H_METPhi.xMax = PI; 
    H_METPhi.yLabel = "Number of Events / ";
    numDivs = (H_METPhi.xMax - H_METPhi.xMin) / (float) H_METPhi.xBinN;
    H_METPhi.yLabel += "NUM"; H_METPhi.yLabel += " radians";
    H_METPhi.xVarKey = "METPhi";
    H_METPhi.doXSyst = true;
    
    HistogramT H_METPhi_noCorr; H_METPhi_noCorr.name = "h_METPhi_noPhiCorr"; 
    H_METPhi_noCorr.xLabel = "#slash{E}_{T} #phi"; H_METPhi_noCorr.xBinN = PhiBinN; H_METPhi_noCorr.xMin = -PI; H_METPhi_noCorr.xMax = PI; 
    H_METPhi_noCorr.yLabel = "Number of Events / ";
    numDivs = (H_METPhi_noCorr.xMax - H_METPhi_noCorr.xMin) / (float) H_METPhi_noCorr.xBinN;
    H_METPhi_noCorr.yLabel += "NUM"; H_METPhi_noCorr.yLabel += " radians";
    H_METPhi_noCorr.xVarKey = "METPhi_noPhiCorr";
    H_METPhi_noCorr.doXSyst = false;
    
    HistogramT H_METX_noPhiCorr; H_METX_noPhiCorr.name = "h_METX_noPhiCorr";
    H_METX_noPhiCorr.xLabel = "#slash{E}_{x} [GeV]"; H_METX_noPhiCorr.xBinN = METXYBinN; H_METX_noPhiCorr.xMin = METXYBinLB; H_METX_noPhiCorr.xMax = METXYBinUB;
    H_METX_noPhiCorr.yLabel = "Number of Events / ";
    numDivs = (H_METX_noPhiCorr.xMax - H_METX_noPhiCorr.xMin) / (float) H_METX_noPhiCorr.xBinN;
    H_METX_noPhiCorr.yLabel += "NUM"; H_METX_noPhiCorr.yLabel += " GeV";
    H_METX_noPhiCorr.xVarKey = "METX_noPhiCorr";
    H_METX_noPhiCorr.doXSyst = false;
    
    HistogramT H_METY_noPhiCorr; H_METY_noPhiCorr.name = "h_METY_noPhiCorr";
    H_METY_noPhiCorr.xLabel = "#slash{E}_{y} [GeV]"; H_METY_noPhiCorr.xBinN = METXYBinN; H_METY_noPhiCorr.xMin = METXYBinLB; H_METY_noPhiCorr.xMax = METXYBinUB;
    H_METY_noPhiCorr.yLabel = "Number of Events / ";
    numDivs = (H_METY_noPhiCorr.xMax - H_METY_noPhiCorr.xMin) / (float) H_METY_noPhiCorr.xBinN;
    H_METY_noPhiCorr.yLabel += "NUM"; H_METY_noPhiCorr.yLabel += " GeV";
    H_METY_noPhiCorr.xVarKey = "METY_noPhiCorr";
    H_METY_noPhiCorr.doXSyst = false;
    
    ///Jet plots/////
    
    HistogramT H_NJets; H_NJets.name = "h_NJets"; 
    H_NJets.xLabel = "N_{jets}"; H_NJets.xBinN = NJetsBinN; H_NJets.xMin = NJetsBinLB; H_NJets.xMax = NJetsBinUB; 
    H_NJets.yLabel = "Events / N_{jets}";
    H_NJets.xVarKey = "NJets";
    H_NJets.doXSyst = true;
    
    HistogramT H_NJetswBTag; H_NJetswBTag.name = "h_NJetswBTag"; 
    H_NJetswBTag.xLabel = "N_{b-jets}"; H_NJetswBTag.xBinN = NJetsBinN; H_NJetswBTag.xMin = NJetsBinLB; H_NJetswBTag.xMax = NJetsBinUB; 
    H_NJetswBTag.yLabel = "Events / N_{b-jets}";
    H_NJetswBTag.xVarKey = "NBJets";
    H_NJetswBTag.doXSyst = true;
    
    HistogramT H_diJetPt; H_diJetPt.name = "h_DiJetPt"; 
    H_diJetPt.xLabel = "diJet pT [GeV]"; H_diJetPt.xBinN = EnergyPtBinN; H_diJetPt.xMin = EnergyPtBinLB; H_diJetPt.xMax = EnergyPtBinUB; 
    H_diJetPt.yLabel = "Number of Events / ";
    numDivs = (H_diJetPt.xMax - H_diJetPt.xMin) / (float) H_diJetPt.xBinN;
    H_diJetPt.yLabel += "NUM"; H_diJetPt.yLabel += " GeV";
    H_diJetPt.xVarKey = "diJetPt";
    H_diJetPt.doXSyst = true;
    
    HistogramT H_diJetInvMass; H_diJetInvMass.name = "h_DiJetInvMass"; 
    H_diJetInvMass.xLabel = "M_{jj} [GeV]"; H_diJetInvMass.xBinN = EnergyPtBinN; H_diJetInvMass.xMin = EnergyPtBinLB; H_diJetInvMass.xMax = EnergyPtBinUB; 
    H_diJetInvMass.yLabel = "Number of Events / ";
    numDivs = (H_diJetInvMass.xMax - H_diJetInvMass.xMin) / (float) H_diJetInvMass.xBinN;
    H_diJetInvMass.yLabel += "NUM"; H_diJetInvMass.yLabel += " GeV";
    H_diJetInvMass.xVarKey = "diJetInvMass";
    H_diJetInvMass.doXSyst = true;
    
    HistogramT H_diJetEta; H_diJetEta.name = "h_DiJetEta"; 
    H_diJetEta.xLabel = "#eta"; H_diJetEta.xBinN = EtaBinN; H_diJetEta.xMin = EtaBinLB; H_diJetEta.xMax = EtaBinUB; 
    H_diJetEta.yLabel = "Number of Events / ";
    numDivs = (H_diJetEta.xMax - H_diJetEta.xMin) / (float) H_diJetEta.xBinN;
    H_diJetEta.yLabel += "NUM"; H_diJetEta.yLabel += " #eta";
    H_diJetEta.xVarKey = "diJetEta";
    H_diJetEta.doXSyst = true;
    
    HistogramT H_diJetPhi; H_diJetPhi.name = "h_DiJetPhi"; 
    H_diJetPhi.xLabel = "#phi_{jj}"; H_diJetPhi.xBinN = PhiBinN; H_diJetPhi.xMin = -PI; H_diJetPhi.xMax = PI; 
    H_diJetPhi.yLabel = "Number of Events / ";
    numDivs = (H_diJetPhi.xMax - H_diJetPhi.xMin) / (float) H_diJetPhi.xBinN;
    H_diJetPhi.yLabel += "NUM"; H_diJetPhi.yLabel += " radians";
    H_diJetPhi.xVarKey = "diJetPhi";
    H_diJetPhi.doXSyst = true;
    
    HistogramT H_diBJetPt; H_diBJetPt.name = "h_DiBJetPt"; 
    H_diBJetPt.xLabel = "diBJet pT  [GeV]"; H_diBJetPt.xBinN = EnergyPtBinN; H_diBJetPt.xMin = EnergyPtBinLB; H_diBJetPt.xMax = EnergyPtBinUB; 
    H_diBJetPt.yLabel = "Number of Events / ";
    numDivs = (H_diBJetPt.xMax - H_diBJetPt.xMin) / (float) H_diBJetPt.xBinN;
    H_diBJetPt.yLabel += "NUM"; H_diBJetPt.yLabel += " GeV";
    H_diBJetPt.xVarKey = "diBJetPt";
    H_diBJetPt.doXSyst = true;
    
    HistogramT H_diBJetInvMass; H_diBJetInvMass.name = "h_DiBJetInvMass"; 
    H_diBJetInvMass.xLabel = "M_{bb} [GeV]"; H_diBJetInvMass.xBinN = EnergyPtBinN; H_diBJetInvMass.xMin = EnergyPtBinLB; H_diBJetInvMass.xMax = EnergyPtBinUB; 
    H_diBJetInvMass.yLabel = "Number of Events / ";
    numDivs = (H_diBJetInvMass.xMax - H_diBJetInvMass.xMin) / (float) H_diBJetInvMass.xBinN;
    H_diBJetInvMass.yLabel += "NUM"; H_diBJetInvMass.yLabel += " GeV";
    H_diBJetInvMass.xVarKey = "diBJetInvMass";
    H_diBJetInvMass.doXSyst = true;
    
    HistogramT H_diBJetEta; H_diBJetEta.name = "h_DiBJetEta"; 
    H_diBJetEta.xLabel = "#eta"; H_diBJetEta.xBinN = EtaBinN; H_diBJetEta.xMin = EtaBinLB; H_diBJetEta.xMax = EtaBinUB; 
    H_diBJetEta.yLabel = "Number of Events / ";
    numDivs = (H_diBJetEta.xMax - H_diBJetEta.xMin) / (float) H_diBJetEta.xBinN;
    H_diBJetEta.yLabel += "NUM"; H_diBJetEta.yLabel += " #eta";
    H_diBJetEta.xVarKey = "diBJetEta";
    H_diBJetEta.doXSyst = true;
    
    HistogramT H_diBJetPhi; H_diBJetPhi.name = "h_DiBJetPhi"; 
    H_diBJetPhi.xLabel = "diBJet #phi_{jj}"; H_diBJetPhi.xBinN = PhiBinN; H_diBJetPhi.xMin = -PI; H_diBJetPhi.xMax = PI; 
    H_diBJetPhi.yLabel = "Number of Events / ";
    numDivs = (H_diBJetPhi.xMax - H_diBJetPhi.xMin) / (float) H_diBJetPhi.xBinN;
    H_diBJetPhi.yLabel += "NUM"; H_diBJetPhi.yLabel += " radians";
    H_diBJetPhi.xVarKey = "diBJetPhi";
    H_diBJetPhi.doXSyst = true;
    
    //Angular correlations
    HistogramT H_DeltaPhiLep0Lep1; H_DeltaPhiLep0Lep1.name = "h_DeltaPhiLep0Lep1"; 
    H_DeltaPhiLep0Lep1.xLabel = "#Delta #phi"; H_DeltaPhiLep0Lep1.xBinN = PhiBinN; H_DeltaPhiLep0Lep1.xMin = 0; H_DeltaPhiLep0Lep1.xMax = PI; 
    H_DeltaPhiLep0Lep1.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep0Lep1.xMax - H_DeltaPhiLep0Lep1.xMin) / (float) H_DeltaPhiLep0Lep1.xBinN;
    H_DeltaPhiLep0Lep1.yLabel += "NUM"; H_DeltaPhiLep0Lep1.yLabel += " radians";
    H_DeltaPhiLep0Lep1.xVarKey = "DPhiLep0Lep1";
    H_DeltaPhiLep0Lep1.doXSyst = false;
    
    HistogramT H_DeltaPhiLep0MET; H_DeltaPhiLep0MET.name = "h_DeltaPhiLep0MET"; 
    H_DeltaPhiLep0MET.xLabel = "#Delta #phi"; H_DeltaPhiLep0MET.xBinN = PhiBinN; H_DeltaPhiLep0MET.xMin = 0; H_DeltaPhiLep0MET.xMax = PI; 
    H_DeltaPhiLep0MET.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep0MET.xMax - H_DeltaPhiLep0MET.xMin) / (float) H_DeltaPhiLep0MET.xBinN;
    H_DeltaPhiLep0MET.yLabel += "NUM"; H_DeltaPhiLep0MET.yLabel += " radians";
    H_DeltaPhiLep0MET.xVarKey = "DPhiLep0MET";
    H_DeltaPhiLep0MET.doXSyst = false;
    
    HistogramT H_DeltaPhiLep1MET; H_DeltaPhiLep1MET.name = "h_DeltaPhiLep1MET"; 
    H_DeltaPhiLep1MET.xLabel = "#Delta #phi"; H_DeltaPhiLep1MET.xBinN = PhiBinN; H_DeltaPhiLep1MET.xMin = 0; H_DeltaPhiLep1MET.xMax = PI; 
    H_DeltaPhiLep1MET.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep1MET.xMax - H_DeltaPhiLep1MET.xMin) / (float) H_DeltaPhiLep1MET.xBinN;
    H_DeltaPhiLep1MET.yLabel += "NUM"; H_DeltaPhiLep1MET.yLabel += " radians";
    H_DeltaPhiLep1MET.xVarKey = "DPhiLep1MET";
    H_DeltaPhiLep1MET.doXSyst = false;
    
    HistogramT H_DeltaPhiLep0MET_PreCorr; H_DeltaPhiLep0MET_PreCorr.name = "h_DeltaPhiLep0MET_PreCorr"; 
    H_DeltaPhiLep0MET_PreCorr.xLabel = "#Delta #phi"; H_DeltaPhiLep0MET_PreCorr.xBinN = PhiBinN; H_DeltaPhiLep0MET_PreCorr.xMin = 0; H_DeltaPhiLep0MET_PreCorr.xMax = PI; 
    H_DeltaPhiLep0MET_PreCorr.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep0MET_PreCorr.xMax - H_DeltaPhiLep0MET_PreCorr.xMin) / (float) H_DeltaPhiLep0MET_PreCorr.xBinN;
    H_DeltaPhiLep0MET_PreCorr.yLabel += "NUM"; H_DeltaPhiLep0MET_PreCorr.yLabel += " radians";
    H_DeltaPhiLep0MET_PreCorr.xVarKey = "DPhiLep0MET_PreCorr";
    H_DeltaPhiLep0MET_PreCorr.doXSyst = false;
    
    HistogramT H_DeltaPhiLep1MET_PreCorr; H_DeltaPhiLep1MET_PreCorr.name = "h_DeltaPhiLep1MET_PreCorr"; 
    H_DeltaPhiLep1MET_PreCorr.xLabel = "#Delta #phi"; H_DeltaPhiLep1MET_PreCorr.xBinN = PhiBinN; H_DeltaPhiLep1MET_PreCorr.xMin = 0; H_DeltaPhiLep1MET_PreCorr.xMax = PI; 
    H_DeltaPhiLep1MET_PreCorr.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep1MET_PreCorr.xMax - H_DeltaPhiLep1MET_PreCorr.xMin) / (float) H_DeltaPhiLep1MET_PreCorr.xBinN;
    H_DeltaPhiLep1MET_PreCorr.yLabel += "NUM"; H_DeltaPhiLep1MET_PreCorr.yLabel += " radians";
    H_DeltaPhiLep1MET_PreCorr.xVarKey = "DPhiLep1MET_PreCorr";
    H_DeltaPhiLep1MET_PreCorr.doXSyst = false;

    HistogramT H_DeltaPhiZMET; H_DeltaPhiZMET.name = "h_DeltaPhiZMET"; 
    H_DeltaPhiZMET.xLabel = "#Delta #phi"; H_DeltaPhiZMET.xBinN = PhiBinN; H_DeltaPhiZMET.xMin = 0; H_DeltaPhiZMET.xMax = PI; 
    H_DeltaPhiZMET.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiZMET.xMax - H_DeltaPhiZMET.xMin) / (float) H_DeltaPhiZMET.xBinN;
    H_DeltaPhiZMET.yLabel += "NUM"; H_DeltaPhiZMET.yLabel += " radians";
    H_DeltaPhiZMET.xVarKey = "DPhiZMET";
    H_DeltaPhiZMET.doXSyst = false;
    
    HistogramT H_DeltaPhiZMET_PreCorr; H_DeltaPhiZMET_PreCorr.name = "h_DeltaPhiZMET_PreCorr"; 
    H_DeltaPhiZMET_PreCorr.xLabel = "#Delta #phi"; H_DeltaPhiZMET_PreCorr.xBinN = PhiBinN; H_DeltaPhiZMET_PreCorr.xMin = 0; H_DeltaPhiZMET_PreCorr.xMax = PI; 
    H_DeltaPhiZMET_PreCorr.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiZMET_PreCorr.xMax - H_DeltaPhiZMET_PreCorr.xMin) / (float) H_DeltaPhiZMET_PreCorr.xBinN;
    H_DeltaPhiZMET_PreCorr.yLabel += "NUM"; H_DeltaPhiZMET_PreCorr.yLabel += " radians";
    H_DeltaPhiZMET_PreCorr.xVarKey = "DPhiZMET_PreCorr";
    H_DeltaPhiZMET_PreCorr.doXSyst = false;
    
    
    
    
    HistogramT H_DeltaPhiLep0Jet0; H_DeltaPhiLep0Jet0.name = "h_DeltaPhiLep0Jet0"; 
    H_DeltaPhiLep0Jet0.xLabel = "#Delta #phi Lead Lep:Lead Jet"; H_DeltaPhiLep0Jet0.xBinN = PhiBinN; H_DeltaPhiLep0Jet0.xMin = 0; H_DeltaPhiLep0Jet0.xMax = PI; 
    H_DeltaPhiLep0Jet0.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep0Jet0.xMax - H_DeltaPhiLep0Jet0.xMin) / (float) H_DeltaPhiLep0Jet0.xBinN;
    H_DeltaPhiLep0Jet0.yLabel += "NUM"; H_DeltaPhiLep0Jet0.yLabel += " radians";
    H_DeltaPhiLep0Jet0.xVarKey = "DPhiLep0Jet0";
    H_DeltaPhiLep0Jet0.doXSyst = false;
    
    HistogramT H_DeltaPhiLep0Jet1; H_DeltaPhiLep0Jet1.name = "h_DeltaPhiLep0Jet1"; 
    H_DeltaPhiLep0Jet1.xLabel = "#Delta #phi Lead Lep:subLead Jet"; H_DeltaPhiLep0Jet1.xBinN = PhiBinN; H_DeltaPhiLep0Jet1.xMin = 0; H_DeltaPhiLep0Jet1.xMax = PI; 
    H_DeltaPhiLep0Jet1.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep0Jet1.xMax - H_DeltaPhiLep0Jet1.xMin) / (float) H_DeltaPhiLep0Jet1.xBinN;
    H_DeltaPhiLep0Jet1.yLabel += "NUM"; H_DeltaPhiLep0Jet1.yLabel += " radians";
    H_DeltaPhiLep0Jet1.xVarKey = "DPhiLep0Jet1";
    H_DeltaPhiLep0Jet1.doXSyst = false;
    
    HistogramT H_DeltaPhiLep0BJet0; H_DeltaPhiLep0BJet0.name = "h_DeltaPhiLep0BJet0"; 
    H_DeltaPhiLep0BJet0.xLabel = "#Delta #phi Lead Lep:Lead BJet"; H_DeltaPhiLep0BJet0.xBinN = PhiBinN; H_DeltaPhiLep0BJet0.xMin = 0; H_DeltaPhiLep0BJet0.xMax = PI; 
    H_DeltaPhiLep0BJet0.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep0BJet0.xMax - H_DeltaPhiLep0BJet0.xMin) / (float) H_DeltaPhiLep0BJet0.xBinN;
    H_DeltaPhiLep0BJet0.yLabel += "NUM"; H_DeltaPhiLep0BJet0.yLabel += " radians";
    H_DeltaPhiLep0BJet0.xVarKey = "DPhiLep0BJet0";
    H_DeltaPhiLep0BJet0.doXSyst = false;
    
    HistogramT H_DeltaPhiLep0BJet1; H_DeltaPhiLep0BJet1.name = "h_DeltaPhiLep0BJet1"; 
    H_DeltaPhiLep0BJet1.xLabel = "#Delta #phi Lead Lep:subLead BJet"; H_DeltaPhiLep0BJet1.xBinN = PhiBinN; H_DeltaPhiLep0BJet1.xMin = 0; H_DeltaPhiLep0BJet1.xMax = PI; 
    H_DeltaPhiLep0BJet1.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiLep0BJet1.xMax - H_DeltaPhiLep0BJet1.xMin) / (float) H_DeltaPhiLep0BJet1.xBinN;
    H_DeltaPhiLep0BJet1.yLabel += "NUM"; H_DeltaPhiLep0BJet1.yLabel += " radians";
    H_DeltaPhiLep0BJet1.xVarKey = "DPhiLep0BJet1";
    H_DeltaPhiLep0BJet1.doXSyst = false;
    
    HistogramT H_DeltaPhiJet0BJet0; H_DeltaPhiJet0BJet0.name = "h_DeltaPhiJet0BJet0"; 
    H_DeltaPhiJet0BJet0.xLabel = "#Delta #phi Lead Jet:Lead BJet"; H_DeltaPhiJet0BJet0.xBinN = PhiBinN; H_DeltaPhiJet0BJet0.xMin = 0; H_DeltaPhiJet0BJet0.xMax = PI; 
    H_DeltaPhiJet0BJet0.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiJet0BJet0.xMax - H_DeltaPhiJet0BJet0.xMin) / (float) H_DeltaPhiJet0BJet0.xBinN;
    H_DeltaPhiJet0BJet0.yLabel += "NUM"; H_DeltaPhiJet0BJet0.yLabel += " radians";
    H_DeltaPhiJet0BJet0.xVarKey = "DPhiJet0BJet0";
    H_DeltaPhiJet0BJet0.doXSyst = false;
    
    HistogramT H_DeltaPhiJet1BJet1; H_DeltaPhiJet1BJet1.name = "h_DeltaPhiJet1BJet1"; 
    H_DeltaPhiJet1BJet1.xLabel = "#Delta #phi subLead Jet:subLead BJet"; H_DeltaPhiJet1BJet1.xBinN = PhiBinN; H_DeltaPhiJet1BJet1.xMin = 0; H_DeltaPhiJet1BJet1.xMax = PI; 
    H_DeltaPhiJet1BJet1.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiJet1BJet1.xMax - H_DeltaPhiJet1BJet1.xMin) / (float) H_DeltaPhiJet1BJet1.xBinN;
    H_DeltaPhiJet1BJet1.yLabel += "NUM"; H_DeltaPhiJet1BJet1.yLabel += " radians";
    H_DeltaPhiJet1BJet1.xVarKey = "DPhiJet1BJet1";
    H_DeltaPhiJet1BJet1.doXSyst = false;
    
    HistogramT H_DeltaPhiJet1BJet0; H_DeltaPhiJet1BJet0.name = "h_DeltaPhiJet1BJet0"; 
    H_DeltaPhiJet1BJet0.xLabel = "#Delta #phi subLead Jet:Lead BJet"; H_DeltaPhiJet1BJet0.xBinN = PhiBinN; H_DeltaPhiJet1BJet0.xMin = 0; H_DeltaPhiJet1BJet0.xMax = PI; 
    H_DeltaPhiJet1BJet0.yLabel = "Number of Events / ";
    numDivs = (H_DeltaPhiJet1BJet0.xMax - H_DeltaPhiJet1BJet0.xMin) / (float) H_DeltaPhiJet1BJet0.xBinN;
    H_DeltaPhiJet1BJet0.yLabel += "NUM"; H_DeltaPhiJet1BJet0.yLabel += " radians";
    H_DeltaPhiJet1BJet0.xVarKey = "DPhiJet1BJet0";
    H_DeltaPhiJet1BJet0.doXSyst = false;
    
    //polarization energy "E-malgamation"
    HistogramT H_diLepEMinJetE; H_diLepEMinJetE.name = "h_diLepEMinJetE"; 
    H_diLepEMinJetE.xLabel = "E(l^{+}) + E(l^{-}) - E(J1) - E(J2) [GeV]"; H_diLepEMinJetE.xBinN = EnergyPtBinN; H_diLepEMinJetE.xMin = -EnergyPtBinUB; H_diLepEMinJetE.xMax = EnergyPtBinUB; 
    H_diLepEMinJetE.yLabel = "Number of Events / ";
    numDivs = (H_diLepEMinJetE.xMax - H_diLepEMinJetE.xMin) / (float) H_diLepEMinJetE.xBinN;
    H_diLepEMinJetE.yLabel += "NUM"; H_diLepEMinJetE.yLabel += " GeV";
    H_diLepEMinJetE.xVarKey = "ELepEJet";
    H_diLepEMinJetE.doXSyst = true;
    
    HistogramT H_HT; H_HT.name = "h_HT"; 
    H_HT.xLabel = "H_{T} [GeV]"; H_HT.xBinN = EnergyPtBinN; H_HT.xMin = 0; H_HT.xMax = 3 * EnergyPtBinUB; 
    H_HT.yLabel = "Number of Events / ";
    numDivs = (H_HT.xMax - H_HT.xMin) / (float) H_HT.xBinN;
    H_HT.yLabel += "NUM"; H_HT.yLabel += " GeV";
    H_HT.xVarKey = "HT";
    H_HT.doXSyst = true;
    
    ///nVtx control plot
    HistogramT H_nVtx; H_nVtx.name = "h_nVtx"; 
    H_nVtx.xLabel = "N_{vtx}^{reco}"; H_nVtx.xBinN = nVtxBinN; H_nVtx.xMin = nVtxBinLB; H_nVtx.xMax = nVtxBinUB; 
    H_nVtx.yLabel = "Events / N_{vtx}^{reco}";
    H_nVtx.xVarKey = "nVtx";
    H_nVtx.doXSyst = false;
    
    HistogramT H_nVtx_preRW; H_nVtx_preRW.name = "h_nVtx_preRW"; 
    H_nVtx_preRW.xLabel = "N_{vtx}^{reco}"; H_nVtx_preRW.xBinN = nVtxBinN; H_nVtx_preRW.xMin = nVtxBinLB; H_nVtx_preRW.xMax = nVtxBinUB; 
    H_nVtx_preRW.yLabel = "Events / N_{vtx}^{reco}";
    H_nVtx_preRW.xVarKey = "nVtx"; 
    H_nVtx_preRW.doXSyst = false;
    
    HistogramT H_nVtxTrue; H_nVtxTrue.name = "h_nVtxTrue"; 
    H_nVtxTrue.xLabel = "N_{vtx}^{reco}"; H_nVtxTrue.xBinN = nVtxBinN; H_nVtxTrue.xMin = nVtxBinLB; H_nVtxTrue.xMax = nVtxBinUB; 
    H_nVtxTrue.yLabel = "Events / N_{vtx}^{reco}";
    H_nVtxTrue.xVarKey = "nVtxTrue";
    H_nVtxTrue.doXSyst = false;
    
    HistogramT H_nVtxTrue_preRW; H_nVtxTrue_preRW.name = "h_nVtxTrue_preRW"; 
    H_nVtxTrue_preRW.xLabel = "N_{vtx}^{reco}"; H_nVtxTrue_preRW.xBinN = nVtxBinN; H_nVtxTrue_preRW.xMin = nVtxBinLB; H_nVtxTrue_preRW.xMax = nVtxBinUB; 
    H_nVtxTrue_preRW.yLabel = "Events / N_{vtx}^{reco}";
    H_nVtxTrue_preRW.xVarKey = "nVtxTrue"; 
    H_nVtxTrue_preRW.doXSyst = false;
    
    //push the 1D histograms structures into a vector for eventual use in booking histograms
    vector<HistogramT> * histVec_1D = new vector<HistogramT>;
    histVec_1D->push_back(H_leadLepPt); histVec_1D->push_back(H_leadLepEta); histVec_1D->push_back(H_subLepPt); histVec_1D->push_back(H_subLepEta);
    histVec_1D->push_back(H_leadJetPt); histVec_1D->push_back(H_leadJetEta); histVec_1D->push_back(H_subJetPt); histVec_1D->push_back(H_subJetEta);
    histVec_1D->push_back(H_leadBJetPt); histVec_1D->push_back(H_leadBJetEta); histVec_1D->push_back(H_leadBJetEn); histVec_1D->push_back(H_subBJetPt); histVec_1D->push_back(H_subBJetEta); histVec_1D->push_back(H_subBJetEn);
    histVec_1D->push_back(H_diLepPt); histVec_1D->push_back(H_diLepInvMass); histVec_1D->push_back(H_diLepEta); histVec_1D->push_back(H_diLepPhi); 
    histVec_1D->push_back(H_MT2ll); histVec_1D->push_back(H_MT2lb);
    histVec_1D->push_back(H_MT2llCont); histVec_1D->push_back(H_MT2lbCont);
    histVec_1D->push_back(H_PassMT2llCut80); histVec_1D->push_back(H_PassMT2llCut90); histVec_1D->push_back(H_PassMT2llCut100); histVec_1D->push_back(H_PassMT2llCut110); histVec_1D->push_back(H_PassMT2llCut120);
    histVec_1D->push_back(H_MT2ll_DPhiZMETClose); histVec_1D->push_back(H_MT2ll_DPhiZMETMid); histVec_1D->push_back(H_MT2ll_DPhiZMETFar);    
    histVec_1D->push_back(H_MET); histVec_1D->push_back(H_METX); histVec_1D->push_back(H_METY); histVec_1D->push_back(H_METPhi); histVec_1D->push_back(H_METPhi_noCorr); histVec_1D->push_back(H_METX_noPhiCorr); histVec_1D->push_back(H_METY_noPhiCorr);
    histVec_1D->push_back(H_NJets); histVec_1D->push_back(H_NJetswBTag);
    histVec_1D->push_back(H_diJetPt); histVec_1D->push_back(H_diJetInvMass); histVec_1D->push_back(H_diJetEta); histVec_1D->push_back(H_diJetPhi); 
    histVec_1D->push_back(H_diBJetPt); histVec_1D->push_back(H_diBJetInvMass); histVec_1D->push_back(H_diBJetEta); histVec_1D->push_back(H_diBJetPhi); 
    histVec_1D->push_back(H_DeltaPhiLep0Lep1); 
    histVec_1D->push_back(H_DeltaPhiLep0MET); histVec_1D->push_back(H_DeltaPhiLep1MET); 
    histVec_1D->push_back(H_DeltaPhiLep0MET_PreCorr); histVec_1D->push_back(H_DeltaPhiLep1MET_PreCorr); 
    histVec_1D->push_back(H_DeltaPhiZMET); histVec_1D->push_back(H_DeltaPhiZMET_PreCorr); 
    histVec_1D->push_back(H_DeltaPhiLep0Jet0); histVec_1D->push_back(H_DeltaPhiLep0Jet1);
    histVec_1D->push_back(H_DeltaPhiLep0BJet0); histVec_1D->push_back(H_DeltaPhiLep0BJet1);
    histVec_1D->push_back(H_DeltaPhiJet0BJet0); histVec_1D->push_back(H_DeltaPhiJet1BJet1); histVec_1D->push_back(H_DeltaPhiJet1BJet0);
    histVec_1D->push_back(H_diLepEMinJetE); histVec_1D->push_back(H_HT); 
    histVec_1D->push_back(H_nVtx); histVec_1D->push_back(H_nVtx_preRW);
    histVec_1D->push_back(H_nVtxTrue); histVec_1D->push_back(H_nVtxTrue_preRW);
    
    return histVec_1D;
}

inline vector<HistogramT> * TwoDeeHistTVec() {
    
    vector<HistogramT> * histVec_2D = new vector<HistogramT>;
    
    const double PI = 3.14159265;
    //    int EnergyPtBinN = 400;
    //    int EtaBinN      = 200;
    int PhiBinN      = 100;
    int METBinN      = 40;    
    int METXYBinN    = 100;    
    //    int NJetsBinN    = 11;
    int nVtxBinN     = 35;
    
    //    float EnergyPtBinLB = 0;
    //    float EnergyPtBinUB = 400;
    //    float EtaBinLB      = -6;
    //    float EtaBinUB      = 6;
    float METBinLB      = 0;
    float METBinUB      = 400;
    float METXYBinLB    = -200;
    float METXYBinUB    = 200;
    //    float NJetsBinLB    = -0.5;
    //    float NJetsBinUB    = 10.5;
    float nVtxBinLB     = 0.5;
    float nVtxBinUB     = 35.5;
    
    ///////analogous case for 2D histograms///////
    HistogramT H_METX_vs_nVtx; H_METX_vs_nVtx.name = "h_METX_vs_nVtx";
    H_METX_vs_nVtx.xLabel = "N_{vtx}^{reco}"; H_METX_vs_nVtx.xBinN = nVtxBinN; H_METX_vs_nVtx.xMin = nVtxBinLB; H_METX_vs_nVtx.xMax = nVtxBinUB;
    H_METX_vs_nVtx.yLabel = "#slash{E}_{x} [GeV]"; H_METX_vs_nVtx.yBinN = METXYBinN; H_METX_vs_nVtx.yMin = METXYBinLB; H_METX_vs_nVtx.yMax = METXYBinUB;
    H_METX_vs_nVtx.xVarKey = "nVtx";
    H_METX_vs_nVtx.yVarKey = "METX";
    H_METX_vs_nVtx.doYSyst = true;
    
    HistogramT H_METY_vs_nVtx; H_METY_vs_nVtx.name = "h_METY_vs_nVtx";
    H_METY_vs_nVtx.xLabel = "N_{vtx}^{reco}"; H_METY_vs_nVtx.xBinN = nVtxBinN; H_METY_vs_nVtx.xMin = nVtxBinLB; H_METY_vs_nVtx.xMax = nVtxBinUB;
    H_METY_vs_nVtx.yLabel = "#slash{E}_{y} [GeV]"; H_METY_vs_nVtx.yBinN = METXYBinN; H_METY_vs_nVtx.yMin = METXYBinLB; H_METY_vs_nVtx.yMax = METXYBinUB;
    H_METY_vs_nVtx.xVarKey = "nVtx";
    H_METY_vs_nVtx.yVarKey = "METY";
    H_METY_vs_nVtx.doYSyst = true;
    
    HistogramT H_METX_vs_nVtx_noPhiCorr; H_METX_vs_nVtx_noPhiCorr.name = "h_METX_vs_nVtx_noPhiCorr";
    H_METX_vs_nVtx_noPhiCorr.xLabel = "N_{vtx}^{reco}"; H_METX_vs_nVtx_noPhiCorr.xBinN = nVtxBinN; H_METX_vs_nVtx_noPhiCorr.xMin = nVtxBinLB; H_METX_vs_nVtx_noPhiCorr.xMax = nVtxBinUB;
    H_METX_vs_nVtx_noPhiCorr.yLabel = "#slash{E}_{x} [GeV]"; H_METX_vs_nVtx_noPhiCorr.yBinN = METXYBinN; H_METX_vs_nVtx_noPhiCorr.yMin = METXYBinLB; H_METX_vs_nVtx_noPhiCorr.yMax = METXYBinUB;
    H_METX_vs_nVtx_noPhiCorr.xVarKey = "nVtx";
    H_METX_vs_nVtx_noPhiCorr.yVarKey = "METX_noPhiCorr";
    H_METX_vs_nVtx_noPhiCorr.doYSyst = true;
    
    HistogramT H_METY_vs_nVtx_noPhiCorr; H_METY_vs_nVtx_noPhiCorr.name = "h_METY_vs_nVtx_noPhiCorr";
    H_METY_vs_nVtx_noPhiCorr.xLabel = "N_{vtx}^{reco}"; H_METY_vs_nVtx_noPhiCorr.xBinN = nVtxBinN; H_METY_vs_nVtx_noPhiCorr.xMin = nVtxBinLB; H_METY_vs_nVtx_noPhiCorr.xMax = nVtxBinUB;
    H_METY_vs_nVtx_noPhiCorr.yLabel = "#slash{E}_{y} [GeV]"; H_METY_vs_nVtx_noPhiCorr.yBinN = METXYBinN; H_METY_vs_nVtx_noPhiCorr.yMin = METXYBinLB; H_METY_vs_nVtx_noPhiCorr.yMax = METXYBinUB;
    H_METY_vs_nVtx_noPhiCorr.xVarKey = "nVtx";
    H_METY_vs_nVtx_noPhiCorr.yVarKey = "METY_noPhiCorr";
    H_METY_vs_nVtx_noPhiCorr.doYSyst = true;
    
    HistogramT H_MT2ll_vs_DeltaPhiLep0Lep1; H_MT2ll_vs_DeltaPhiLep0Lep1.name = "h_MT2ll_vs_DeltaPhiLep0Lep1";
    H_MT2ll_vs_DeltaPhiLep0Lep1.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_DeltaPhiLep0Lep1.xBinN = METBinN; H_MT2ll_vs_DeltaPhiLep0Lep1.xMin = METBinLB; H_MT2ll_vs_DeltaPhiLep0Lep1.xMax = METBinUB;
    H_MT2ll_vs_DeltaPhiLep0Lep1.yLabel = "#Delta #phi_{ll}"; H_MT2ll_vs_DeltaPhiLep0Lep1.yBinN = PhiBinN; H_MT2ll_vs_DeltaPhiLep0Lep1.yMin = 0; H_MT2ll_vs_DeltaPhiLep0Lep1.yMax = PI;
    H_MT2ll_vs_DeltaPhiLep0Lep1.xVarKey = "MT2ll";
    H_MT2ll_vs_DeltaPhiLep0Lep1.yVarKey = "DPhiLep0Lep1";
    H_MT2ll_vs_DeltaPhiLep0Lep1.doXSyst = true;
    H_MT2ll_vs_DeltaPhiLep0Lep1.doYSyst = true;
    
    HistogramT H_MT2llControl_vs_DeltaPhiLep0Lep1; H_MT2llControl_vs_DeltaPhiLep0Lep1.name = "h_MT2llControl_vs_DeltaPhiLep0Lep1";
    H_MT2llControl_vs_DeltaPhiLep0Lep1.xLabel = "MT2_{ll} [GeV]"; H_MT2llControl_vs_DeltaPhiLep0Lep1.xBinN = 20; H_MT2llControl_vs_DeltaPhiLep0Lep1.xMin = 0; H_MT2llControl_vs_DeltaPhiLep0Lep1.xMax = 80;
    H_MT2llControl_vs_DeltaPhiLep0Lep1.yLabel = "#Delta #phi_{ll}"; H_MT2llControl_vs_DeltaPhiLep0Lep1.yBinN = PhiBinN; H_MT2llControl_vs_DeltaPhiLep0Lep1.yMin = 0; H_MT2llControl_vs_DeltaPhiLep0Lep1.yMax = PI;
    H_MT2llControl_vs_DeltaPhiLep0Lep1.xVarKey = "MT2ll";
    H_MT2llControl_vs_DeltaPhiLep0Lep1.yVarKey = "DPhiLep0Lep1";
    H_MT2llControl_vs_DeltaPhiLep0Lep1.doXSyst = true;
    H_MT2llControl_vs_DeltaPhiLep0Lep1.doYSyst = true;
    
    HistogramT H_MT2lb_vs_DeltaPhiLepB0LepB1; H_MT2lb_vs_DeltaPhiLepB0LepB1.name = "h_MT2lb_vs_DeltaPhiLepB0LepB1";
    H_MT2lb_vs_DeltaPhiLepB0LepB1.xLabel = "MT2lb [GeV]"; H_MT2lb_vs_DeltaPhiLepB0LepB1.xBinN = METBinN; H_MT2lb_vs_DeltaPhiLepB0LepB1.xMin = METBinLB; H_MT2lb_vs_DeltaPhiLepB0LepB1.xMax = METBinUB;
    H_MT2lb_vs_DeltaPhiLepB0LepB1.yLabel = "#Delta #phi_{lb lb}"; H_MT2lb_vs_DeltaPhiLepB0LepB1.yBinN = PhiBinN; H_MT2lb_vs_DeltaPhiLepB0LepB1.yMin = 0; H_MT2lb_vs_DeltaPhiLepB0LepB1.yMax = PI;
    H_MT2lb_vs_DeltaPhiLepB0LepB1.xVarKey = "MT2lb";
    H_MT2lb_vs_DeltaPhiLepB0LepB1.yVarKey = "DPhiLepB0LepB1";
    H_MT2lb_vs_DeltaPhiLepB0LepB1.doXSyst = true;
    H_MT2lb_vs_DeltaPhiLepB0LepB1.doYSyst = true;
    
    HistogramT H_MT2lbControl_vs_DeltaPhiLepB0LepB1; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.name = "h_MT2lbControl_vs_DeltaPhiLepB0LepB1";
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xLabel = "MT2lb [GeV]"; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xBinN = 43; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xMin = 0; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xMax = 172;
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yLabel = "#Delta #phi_{lb lb}"; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yBinN = PhiBinN; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yMin = 0; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yMax = PI;
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xVarKey = "MT2lb";
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yVarKey = "DPhiLepB0LepB1";
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.doXSyst = true;
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.doYSyst = true;
    //push the 2D histograms structures into a vector for eventual use in booking histograms
    histVec_2D->push_back(H_METX_vs_nVtx); histVec_2D->push_back(H_METY_vs_nVtx);
    histVec_2D->push_back(H_METX_vs_nVtx_noPhiCorr); histVec_2D->push_back(H_METY_vs_nVtx_noPhiCorr);
    histVec_2D->push_back(H_MT2ll_vs_DeltaPhiLep0Lep1); histVec_2D->push_back(H_MT2llControl_vs_DeltaPhiLep0Lep1);
    histVec_2D->push_back(H_MT2lb_vs_DeltaPhiLepB0LepB1); histVec_2D->push_back(H_MT2lbControl_vs_DeltaPhiLepB0LepB1);
    return histVec_2D;
}

inline vector<HistogramT> * ThreeDeeHistTVec() {
    
    vector<HistogramT> * histVec_3D = new vector<HistogramT>;
    
    const double PI = 3.14159265;
    //    int EnergyPtBinN = 400;
    //    int EtaBinN      = 200;
    int PhiBinN      = 100;
    int METBinN      = 40;    
    int METXYBinN    = 200;    
    int NJetsBinN    = 11;
    int nVtxBinN     = 35;
    
    //    float EnergyPtBinLB = 0;
    //    float EnergyPtBinUB = 400;
    //    float EtaBinLB      = -6;
    //    float EtaBinUB      = 6;
    float METBinLB      = 0;
    float METBinUB      = 400;
    float METXYBinLB    = -200;
    float METXYBinUB    = 200;
    float NJetsBinLB    = -0.5;
    float NJetsBinUB    = 10.5;
    float nVtxBinLB     = 0.5;
    float nVtxBinUB     = 35.5;
    
    ///////analogous case for 2D histograms///////
    
    HistogramT H_MT2ll_vs_DeltaPhiZMET_vs_nVtx; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.name = "h_MT2ll_vs_DeltaPhiZMET_vs_nVtx"; 
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.xBinN = METBinN; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.xMin = METBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.xMax = METBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.yLabel = "#Delta #phi_{Z, #slash{E}_{T}}"; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.yBinN = PhiBinN; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.yMin = 0.; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.yMax = PI;
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.zLabel = "N_{vtx}^{reco}"; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.zBinN = nVtxBinN; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.zMin = nVtxBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.zMax = nVtxBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.xVarKey = "MT2ll";
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.yVarKey = "DPhiZMET";
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.zVarKey = "nVtx";
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.doXSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_nVtx.doYSyst = true;
    
    HistogramT H_MT2ll_vs_DeltaPhiZMET_vs_NJets; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.name = "h_MT2ll_vs_DeltaPhiZMET_vs_NJets"; 
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.xBinN = METBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.xMin = METBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.xMax = METBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.yLabel = "#Delta #phi_{Z, #slash{E}_{T}}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.yBinN = PhiBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.yMin = 0.; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.yMax = PI;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.zLabel = "N_{jets}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.zBinN = NJetsBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.zMin = NJetsBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets.zMax = NJetsBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.xVarKey = "MT2ll";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.yVarKey = "DPhiZMET";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.zVarKey = "NJets";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.doXSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.doYSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets.doZSyst = true;
    
    HistogramT H_MT2ll_vs_DeltaPhiZMET_vs_NBJets; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.name = "h_MT2ll_vs_DeltaPhiZMET_vs_NBJets"; 
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.xBinN = METBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.xMin = METBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.xMax = METBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.yLabel = "#Delta #phi_{Z, #slash{E}_{T}}"; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.yBinN = PhiBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.yMin = 0.; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.yMax = PI;
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.zLabel = "N_{b-jets}"; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.zBinN = NJetsBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.zMin = NJetsBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.zMax = NJetsBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.xVarKey = "MT2ll";
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.yVarKey = "DPhiZMET";
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.zVarKey = "NBJets";    
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.doXSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.doYSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NBJets.doZSyst = true;
    
    
    HistogramT H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.name = "h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10"; 
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.xBinN = METBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.xMin = METBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.xMax = METBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.yLabel = "#Delta #phi_{Z, #slash{E}_{T}}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.yBinN = PhiBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.yMin = 0.; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.yMax = PI;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.zLabel = "N_{jets}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.zBinN = NJetsBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.zMin = NJetsBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.zMax = NJetsBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.xVarKey = "MT2ll";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.yVarKey = "DPhiZMET";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.zVarKey = "NJets";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.doXSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.doYSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10.doZSyst = true;
    
    HistogramT H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.name = "h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20"; 
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.xBinN = METBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.xMin = METBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.xMax = METBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.yLabel = "#Delta #phi_{Z, #slash{E}_{T}}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.yBinN = PhiBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.yMin = 0.; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.yMax = PI;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.zLabel = "N_{jets}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.zBinN = NJetsBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.zMin = NJetsBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.zMax = NJetsBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.xVarKey = "MT2ll";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.yVarKey = "DPhiZMET";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.zVarKey = "NJets";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.doXSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.doYSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20.doZSyst = true;
    
    HistogramT H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.name = "h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30"; 
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.xBinN = METBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.xMin = METBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.xMax = METBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.yLabel = "#Delta #phi_{Z, #slash{E}_{T}}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.yBinN = PhiBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.yMin = 0.; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.yMax = PI;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.zLabel = "N_{jets}"; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.zBinN = NJetsBinN; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.zMin = NJetsBinLB; H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.zMax = NJetsBinUB;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.xVarKey = "MT2ll";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.yVarKey = "DPhiZMET";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.zVarKey = "NJets";
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.doXSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.doYSyst = true;
    H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30.doZSyst = true;
    
    //push the 2D histograms structures into a vector for eventual use in booking histograms
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_nVtx);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NJets);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NBJets);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30);
    return histVec_3D;
}


inline vector<SampleT> * SubSampVec() {
    /*
     int     whichdiLepType; // -1: inclusive, 0: MuMu, 1: EE, 2: EMu
     int     doZVeto;   // -1: inclusive, 0: ZMass window, 1: outside ZMass window;
     int     cutNJets;  // -1: inclusive, # > 0: require NJets >= #
     int     cutNBJets;  // -1: inclusive, # > 0: require NBJets >= #
     float   cutMET; // # > 0: require MET >= #
     */
    
    
    /*
     SampleT events_EE; events_EE.histNameSuffix = "_ee"; events_EE.histXaxisSuffix = "_{ee}"; 
     events_EE.histYaxisSuffix = events_EE.histXaxisSuffix; events_EE.histZaxisSuffix = events_EE.histXaxisSuffix;
     events_EE.doZVeto = 0; events_EE.cutNJets = -1; events_EE.cutNBJets = -1; events_EE.cutMET = 0;
     
     SampleT events_MuMu; events_MuMu.histNameSuffix = "_mumu"; events_MuMu.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu.histYaxisSuffix = events_MuMu.histXaxisSuffix; events_MuMu.histZaxisSuffix = events_MuMu.histXaxisSuffix;
     events_MuMu.doZVeto = 0; events_MuMu.cutNJets = -1; events_MuMu.cutNBJets = -1; events_MuMu.cutMET = 0;
     
     SampleT events_EMu; events_EMu.histNameSuffix = "_emu"; events_EMu.histXaxisSuffix = "_{e#mu}"; 
     events_EMu.histYaxisSuffix = events_EMu.histXaxisSuffix; events_EMu.histZaxisSuffix = events_EMu.histXaxisSuffix;
     events_EMu.doZVeto = 0; events_EMu.cutNJets = -1; events_EMu.cutNBJets = -1; events_EMu.cutMET = 0;
     */
    /*
     //Zmasswindow veto sub channel
     SampleT events_EE_ZVeto; events_EE_ZVeto.histNameSuffix = "_ee_ZVeto"; events_EE_ZVeto.histXaxisSuffix = "_{ee}"; 
     events_EE_ZVeto.histYaxisSuffix = events_EE_ZVeto.histXaxisSuffix; events_EE_ZVeto.histZaxisSuffix = events_EE_ZVeto.histXaxisSuffix;
     events_EE_ZVeto.doZVeto = 1; events_EE_ZVeto.cutNJets = -1; events_EE_ZVeto.cutNBJets = -1; events_EE_ZVeto.cutMET = 0;
     
     
     SampleT events_MuMu_ZVeto; events_MuMu_ZVeto.histNameSuffix = "_mumu_ZVeto"; events_MuMu_ZVeto.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu_ZVeto.histYaxisSuffix = events_MuMu_ZVeto.histXaxisSuffix; events_MuMu_ZVeto.histZaxisSuffix = events_MuMu_ZVeto.histXaxisSuffix;
     events_MuMu_ZVeto.doZVeto = 1; events_MuMu_ZVeto.cutNJets = -1; events_MuMu_ZVeto.cutNBJets = -1; events_MuMU_ZVeto.cutMET = 0;
     
     SampleT events_EMu_ZVeto; events_EMu_ZVeto.histNameSuffix = "_emu_ZVeto"; events_EMu_ZVeto.histXaxisSuffix = "_{e#mu}"; 
     events_EMu_ZVeto.histYaxisSuffix = events_EMu_ZVeto.histXaxisSuffix; events_EMu_ZVeto.histZaxisSuffix = events_EMu_ZVeto.histXaxisSuffix;
     events_EMu_ZVeto.doZVeto = 1; events_EMu_ZVeto.cutNJets = -1; events_EMu_ZVeto.cutNBJets = -1; events_EMu_ZVeto.cutMET = 0;
     
     */
    /*
     SampleT events_EE_ZVeto_NJetsGeq2_NBJetsGeq1; events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_ee_ZVeto_Jet2BJet1"; events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{ee}"; 
     events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.doZVeto = 1; events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_EE_ZVeto_NJetsGeq2_NBJetsGeq1.cutMET = 0;
     
     
     SampleT events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1; events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_mumu_ZVeto_Jet2BJet1"; events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.doZVeto = 1; events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_MuMu_ZVeto_NJetsGeq2_NBJetsGeq1.cutMET = 0;
     
     SampleT events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1; events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_emu_ZVeto_Jet2BJet1"; events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{e#mu}"; 
     events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.doZVeto = 1; events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_EMu_ZVeto_NJetsGeq2_NBJetsGeq1.cutMET = 0;
     
     SampleT events_EE_NJetsGeq2_NBJetsGeq1; events_EE_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_ee_Jet2BJet1"; events_EE_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{ee}"; 
     events_EE_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EE_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EE_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EE_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     events_EE_NJetsGeq2_NBJetsGeq1.doZVeto = 0; events_EE_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_EE_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_EE_NJetsGeq2_NBJetsGeq1.cutMET = 0;
     
     SampleT events_MuMu_NJetsGeq2_NBJetsGeq1; events_MuMu_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_mumu_Jet2BJet1"; events_MuMu_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_MuMu_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_MuMu_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_MuMu_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     events_MuMu_NJetsGeq2_NBJetsGeq1.doZVeto = 0; events_MuMu_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_MuMu_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_MuMu_NJetsGeq2_NBJetsGeq1.cutMET = 0;
     
     SampleT events_EMu_NJetsGeq2_NBJetsGeq1; events_EMu_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_emu_Jet2BJet1"; events_EMu_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{e#mu}"; 
     events_EMu_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EMu_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EMu_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EMu_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     events_EMu_NJetsGeq2_NBJetsGeq1.doZVeto = 0; events_EMu_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_EMu_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_EMu_NJetsGeq2_NBJetsGeq1.cutMET = 0;
     
     
     SampleT events_EE_ZVeto_NJetsGeq2; events_EE_ZVeto_NJetsGeq2.histNameSuffix = "_ee_ZVeto_Jet2"; events_EE_ZVeto_NJetsGeq2.histXaxisSuffix = "_{ee}"; 
     events_EE_ZVeto_NJetsGeq2.histYaxisSuffix = events_EE_ZVeto_NJetsGeq2.histXaxisSuffix; events_EE_ZVeto_NJetsGeq2.histZaxisSuffix = events_EE_ZVeto_NJetsGeq2.histXaxisSuffix;
     events_EE_ZVeto_NJetsGeq2.doZVeto = 1; events_EE_ZVeto_NJetsGeq2.cutNJets = 2; events_EE_ZVeto_NJetsGeq2.cutNBJets = -1; events_EE_ZVeto_NJetsGeq2.cutMET = 0;
     
     SampleT events_MuMu_ZVeto_NJetsGeq2; events_MuMu_ZVeto_NJetsGeq2.histNameSuffix = "_mumu_ZVeto_Jet2"; events_MuMu_ZVeto_NJetsGeq2.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu_ZVeto_NJetsGeq2.histYaxisSuffix = events_MuMu_ZVeto_NJetsGeq2.histXaxisSuffix; events_MuMu_ZVeto_NJetsGeq2.histZaxisSuffix = events_MuMu_ZVeto_NJetsGeq2.histXaxisSuffix;
     events_MuMu_ZVeto_NJetsGeq2.doZVeto = 1; events_MuMu_ZVeto_NJetsGeq2.cutNJets = 2; events_MuMu_ZVeto_NJetsGeq2.cutNBJets = -1; events_MuMu_ZVeto_NJetsGeq2.cutMET = 0;
     
     SampleT events_EMu_ZVeto_NJetsGeq2; events_EMu_ZVeto_NJetsGeq2.histNameSuffix = "_emu_ZVeto_Jet2"; events_EMu_ZVeto_NJetsGeq2.histXaxisSuffix = "_{e#mu}"; 
     events_EMu_ZVeto_NJetsGeq2.histYaxisSuffix = events_EMu_ZVeto_NJetsGeq2.histXaxisSuffix; events_EMu_ZVeto_NJetsGeq2.histZaxisSuffix = events_EMu_ZVeto_NJetsGeq2.histXaxisSuffix;
     events_EMu_ZVeto_NJetsGeq2.doZVeto = 1; events_EMu_ZVeto_NJetsGeq2.cutNJets = 2; events_EMu_ZVeto_NJetsGeq2.cutNBJets = -1; events_EMu_ZVeto_NJetsGeq2.cutMET = 0;
     
     SampleT events_EE_NJetsGeq2; events_EE_NJetsGeq2.histNameSuffix = "_ee_Jet2"; events_EE_NJetsGeq2.histXaxisSuffix = "_{ee}"; 
     events_EE_NJetsGeq2.histYaxisSuffix = events_EE_NJetsGeq2.histXaxisSuffix; events_EE_NJetsGeq2.histZaxisSuffix = events_EE_NJetsGeq2.histXaxisSuffix;
     
     SampleT events_MuMu_NJetsGeq2; events_MuMu_NJetsGeq2.histNameSuffix = "_mumu_Jet2"; events_MuMu_NJetsGeq2.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu_NJetsGeq2.histYaxisSuffix = events_MuMu_NJetsGeq2.histXaxisSuffix; events_MuMu_NJetsGeq2.histZaxisSuffix = events_MuMu_NJetsGeq2.histXaxisSuffix;
     
     SampleT events_EMu_NJetsGeq2; events_EMu_NJetsGeq2.histNameSuffix = "_emu_Jet2"; events_EMu_NJetsGeq2.histXaxisSuffix = "_{e#mu}"; 
     events_EMu_NJetsGeq2.histYaxisSuffix = events_EMu_NJetsGeq2.histXaxisSuffix; events_EMu_NJetsGeq2.histZaxisSuffix = events_EMu_NJetsGeq2.histXaxisSuffix;
     
     
     //Combo channels    
     
     SampleT events_EE_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1; events_EE_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_ee_ZVeto_METGeq40_Jet2BJet1"; events_EE_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{ee}"; 
     events_EE_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EE_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EE_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EE_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     
     SampleT events_MuMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1; events_MuMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_mumu_ZVeto_METGeq40_Jet2BJet1"; events_MuMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_MuMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_MuMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_MuMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     
     SampleT events_EMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1; events_EMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_emu_ZVeto_METGeq40_Jet2BJet1"; events_EMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{e#mu}"; 
     events_EMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EMu_ZVeto_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     
     SampleT events_EE_METGeq40_NJetsGeq2_NBJetsGeq1; events_EE_METGeq40_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_ee_METGeq40_Jet2BJet1"; events_EE_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{ee}"; 
     events_EE_METGeq40_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EE_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EE_METGeq40_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EE_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     
     SampleT events_MuMu_METGeq40_NJetsGeq2_NBJetsGeq1; events_MuMu_METGeq40_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_mumu_METGeq40_Jet2BJet1"; events_MuMu_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{#mu#mu}"; 
     events_MuMu_METGeq40_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_MuMu_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_MuMu_METGeq40_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_MuMu_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     
     SampleT events_EMu_METGeq40_NJetsGeq2_NBJetsGeq1; events_EMu_METGeq40_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_emu_METGeq40_Jet2BJet1"; events_EMu_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = "_{e#mu}"; 
     events_EMu_METGeq40_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = events_EMu_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix; events_EMu_METGeq40_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = events_EMu_METGeq40_NJetsGeq2_NBJetsGeq1.histXaxisSuffix;
     */
    
    /*######Define Subchannels###########*/
    TString lepNameSuffix[3] = {"_mumu", "_ee", "_emu"};
    TString lepHistXAxisSuffix[3] = {"_{#mu#mu}", "_{ee}", "_{e#mu}"};    
    TString ZVetoString = "_ZVeto";
    TString METCutString = "_METGeq40";
    TString JetCutString = "_Jet2";
    TString BJetCutString = "_BJet1";
    
    SampleT events_LepInZMass[3], events_LepOutZMass[3];
    SampleT events_LepInZMassJet2[3], events_LepOutZMassJet2[3];
    SampleT events_LepInZMassJet2BJet1[3], events_LepOutZMassJet2BJet1[3];
    SampleT events_LepInZMassJet2BJet1MET40[3], events_LepOutZMassJet2BJet1MET40[3];
    
    SampleT events_LepInZMassBothinBarrel[3];
    SampleT events_LepInZMassOneinBarrel[3];
    SampleT events_LepInZMassBothinEndcap[3];
    
    SampleT events_LepInZMassor0BJetsJet2MET40[3];
    SampleT events_LepInZMass0BJetsJet2MET40[3];
    SampleT events_LepOutZMass0BJetsJet2MET40[3];
    
    SampleT events_LepInZMassor0BJetsJet2[3];
    SampleT events_LepInZMass0BJetsJet2[3];
    SampleT events_LepOutZMass0BJetsJet2[3];
    for (int i = 0; i < 3; ++i) {
        events_LepInZMass[i].histNameSuffix  = lepNameSuffix[i];
        events_LepInZMass[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass[i].whichdiLepType  = i;
        events_LepInZMass[i].blindDataChannel  = 0;
        
        events_LepOutZMass[i].histNameSuffix  = lepNameSuffix[i] + ZVetoString;
        events_LepOutZMass[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass[i].whichdiLepType  = i;
        events_LepOutZMass[i].blindDataChannel  = 0;
        
        events_LepInZMassJet2[i].histNameSuffix  = lepNameSuffix[i] + JetCutString;
        events_LepInZMassJet2[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2[i].whichdiLepType  = i;
        events_LepInZMassJet2[i].blindDataChannel  = 0;
        
        events_LepOutZMassJet2[i].histNameSuffix  = lepNameSuffix[i] + ZVetoString + JetCutString;
        events_LepOutZMassJet2[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2[i].whichdiLepType  = i;
        events_LepOutZMassJet2[i].blindDataChannel  = 0;
        
        events_LepInZMassJet2BJet1[i].histNameSuffix  = lepNameSuffix[i] + JetCutString + BJetCutString;
        events_LepInZMassJet2BJet1[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2BJet1[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2BJet1[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2BJet1[i].whichdiLepType  = i;
        events_LepInZMassJet2BJet1[i].blindDataChannel  = 0;
        
        events_LepOutZMassJet2BJet1[i].histNameSuffix  = lepNameSuffix[i] + ZVetoString + JetCutString + BJetCutString;
        events_LepOutZMassJet2BJet1[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2BJet1[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2BJet1[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2BJet1[i].whichdiLepType  = i;
        events_LepOutZMassJet2BJet1[i].blindDataChannel  = 0;
        
        events_LepInZMassJet2BJet1MET40[i].histNameSuffix  = lepNameSuffix[i] + JetCutString + BJetCutString + METCutString;
        events_LepInZMassJet2BJet1MET40[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2BJet1MET40[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2BJet1MET40[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassJet2BJet1MET40[i].whichdiLepType  = i;
        events_LepInZMassJet2BJet1MET40[i].blindDataChannel  = 0;
        
        events_LepOutZMassJet2BJet1MET40[i].histNameSuffix  = lepNameSuffix[i] + ZVetoString + JetCutString + BJetCutString + METCutString;
        events_LepOutZMassJet2BJet1MET40[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2BJet1MET40[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2BJet1MET40[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMassJet2BJet1MET40[i].whichdiLepType  = i;
        events_LepOutZMassJet2BJet1MET40[i].blindDataChannel  = 0;
        
        events_LepInZMassBothinBarrel[i].histNameSuffix = lepNameSuffix[i] + TString("_BothinBarrel");
        events_LepInZMassBothinBarrel[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassBothinBarrel[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassBothinBarrel[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassBothinBarrel[i].whichdiLepType  = i;
        events_LepInZMassBothinBarrel[i].blindDataChannel  = 0;
        
        events_LepInZMassOneinBarrel[i].histNameSuffix = lepNameSuffix[i] + TString("_OneinBarrel");
        events_LepInZMassOneinBarrel[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassOneinBarrel[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassOneinBarrel[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassOneinBarrel[i].whichdiLepType  = i;
        events_LepInZMassOneinBarrel[i].blindDataChannel  = 0;
        
        events_LepInZMassBothinEndcap[i].histNameSuffix = lepNameSuffix[i] + TString("_BothinEndcap");
        events_LepInZMassBothinEndcap[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassBothinEndcap[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassBothinEndcap[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassBothinEndcap[i].whichdiLepType  = i;
        events_LepInZMassBothinEndcap[i].blindDataChannel  = 0;
        
        
        events_LepInZMassor0BJetsJet2[i].histNameSuffix = lepNameSuffix[i] + JetCutString + TString("_inZMass_or0BJets");
        events_LepInZMassor0BJetsJet2[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassor0BJetsJet2[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassor0BJetsJet2[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassor0BJetsJet2[i].whichdiLepType  = i;
        events_LepInZMassor0BJetsJet2[i].blindDataChannel  = 0;
        
        events_LepInZMass0BJetsJet2[i].histNameSuffix = lepNameSuffix[i] + JetCutString + TString("_0BJets");
        events_LepInZMass0BJetsJet2[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass0BJetsJet2[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass0BJetsJet2[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass0BJetsJet2[i].whichdiLepType  = i;
        events_LepInZMass0BJetsJet2[i].blindDataChannel  = 0;
        
        events_LepOutZMass0BJetsJet2[i].histNameSuffix = lepNameSuffix[i] + JetCutString + ZVetoString + TString("_0BJets");
        events_LepOutZMass0BJetsJet2[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass0BJetsJet2[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass0BJetsJet2[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass0BJetsJet2[i].whichdiLepType  = i;
        events_LepOutZMass0BJetsJet2[i].blindDataChannel  = 0;        
        
        events_LepInZMassor0BJetsJet2MET40[i].histNameSuffix = lepNameSuffix[i] + JetCutString + METCutString + TString("_inZMass_or0BJets");
        events_LepInZMassor0BJetsJet2MET40[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassor0BJetsJet2MET40[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassor0BJetsJet2MET40[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMassor0BJetsJet2MET40[i].whichdiLepType  = i;
        events_LepInZMassor0BJetsJet2MET40[i].blindDataChannel  = 0;
        
        events_LepInZMass0BJetsJet2MET40[i].histNameSuffix = lepNameSuffix[i] + JetCutString + METCutString + TString("_0BJets");
        events_LepInZMass0BJetsJet2MET40[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass0BJetsJet2MET40[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass0BJetsJet2MET40[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepInZMass0BJetsJet2MET40[i].whichdiLepType  = i;
        events_LepInZMass0BJetsJet2MET40[i].blindDataChannel  = 0;
        
        events_LepOutZMass0BJetsJet2MET40[i].histNameSuffix = lepNameSuffix[i] + JetCutString + METCutString + ZVetoString + TString("_0BJets");
        events_LepOutZMass0BJetsJet2MET40[i].histXaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass0BJetsJet2MET40[i].histYaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass0BJetsJet2MET40[i].histZaxisSuffix = lepHistXAxisSuffix[i];
        events_LepOutZMass0BJetsJet2MET40[i].whichdiLepType  = i;
        events_LepOutZMass0BJetsJet2MET40[i].blindDataChannel  = 0;        
    }
    for (int j = 0; j < 3; ++j) {
        
        events_LepInZMass[j].doZVeto         = 0;
        events_LepInZMass[j].cutNJets        = -1;
        events_LepInZMass[j].cutNBJets       = -1;
        events_LepInZMass[j].cutMET          = 0.;
        
        events_LepOutZMass[j].doZVeto         = 1;
        events_LepOutZMass[j].cutNJets        = -1;
        events_LepOutZMass[j].cutNBJets       = -1;
        events_LepOutZMass[j].cutMET          = 0.;
        
        events_LepInZMassJet2[j].doZVeto         = 0;
        events_LepInZMassJet2[j].cutNJets        = 2;
        events_LepInZMassJet2[j].cutNBJets       = -1;
        events_LepInZMassJet2[j].cutMET          = 0.;
        
        events_LepOutZMassJet2[j].doZVeto         = 1;
        events_LepOutZMassJet2[j].cutNJets        = 2;
        events_LepOutZMassJet2[j].cutNBJets       = -1;
        events_LepOutZMassJet2[j].cutMET          = 0.;
        
        events_LepInZMassJet2BJet1[j].doZVeto         = 0;
        events_LepInZMassJet2BJet1[j].cutNJets        = 2;
        events_LepInZMassJet2BJet1[j].cutNBJets       = 1;
        events_LepInZMassJet2BJet1[j].cutMET          = 0.;
        
        events_LepOutZMassJet2BJet1[j].doZVeto         = 1;
        events_LepOutZMassJet2BJet1[j].cutNJets        = 2;
        events_LepOutZMassJet2BJet1[j].cutNBJets       = 1;
        events_LepOutZMassJet2BJet1[j].cutMET          = 0.;
        
        events_LepInZMassJet2BJet1MET40[j].doZVeto         = 0;
        events_LepInZMassJet2BJet1MET40[j].cutNJets        = 2;
        events_LepInZMassJet2BJet1MET40[j].cutNBJets       = 1;
        events_LepInZMassJet2BJet1MET40[j].cutMET          = 40.;
        
        events_LepOutZMassJet2BJet1MET40[j].doZVeto         = 1;
        events_LepOutZMassJet2BJet1MET40[j].cutNJets        = 2;
        events_LepOutZMassJet2BJet1MET40[j].cutNBJets       = 1;
        events_LepOutZMassJet2BJet1MET40[j].cutMET          = 40.;
        events_LepOutZMassJet2BJet1MET40[j].blindDataChannel = 1; // blind with MET + ZVeto + jet cuts            
        
        events_LepInZMassBothinBarrel[j].doZVeto         = 0;
        events_LepInZMassBothinBarrel[j].cutNJets        = -1;
        events_LepInZMassBothinBarrel[j].cutNBJets       = -1;
        events_LepInZMassBothinBarrel[j].cutMET          = 0.;
        
        events_LepInZMassOneinBarrel[j].doZVeto         = 0;
        events_LepInZMassOneinBarrel[j].cutNJets        = -1;
        events_LepInZMassOneinBarrel[j].cutNBJets       = -1;
        events_LepInZMassOneinBarrel[j].cutMET          = 0.;
        
        events_LepInZMassBothinEndcap[j].doZVeto         = 0;
        events_LepInZMassBothinEndcap[j].cutNJets        = -1;
        events_LepInZMassBothinEndcap[j].cutNBJets       = -1;
        events_LepInZMassBothinEndcap[j].cutMET          = 0.;                        
        
        events_LepInZMassor0BJetsJet2[j].doZVeto   = -1;
        events_LepInZMassor0BJetsJet2[j].cutNJets  = 2;
        events_LepInZMassor0BJetsJet2[j].cutNBJets = -1;
        events_LepInZMassor0BJetsJet2[j].cutMET    = 0.;
        
        events_LepInZMass0BJetsJet2[j].doZVeto   = 0;
        events_LepInZMass0BJetsJet2[j].cutNJets  = 2;
        events_LepInZMass0BJetsJet2[j].cutNBJets = -1;
        events_LepInZMass0BJetsJet2[j].cutMET    = 0.;
        
        events_LepOutZMass0BJetsJet2[j].doZVeto   = 1;
        events_LepOutZMass0BJetsJet2[j].cutNJets  = 2;
        events_LepOutZMass0BJetsJet2[j].cutNBJets = -1;
        events_LepOutZMass0BJetsJet2[j].cutMET    = 0.;
        
        events_LepInZMassor0BJetsJet2MET40[j].doZVeto   = -1;
        events_LepInZMassor0BJetsJet2MET40[j].cutNJets  = 2;
        events_LepInZMassor0BJetsJet2MET40[j].cutNBJets = -1;
        events_LepInZMassor0BJetsJet2MET40[j].cutMET    = 40.;
        
        events_LepInZMass0BJetsJet2MET40[j].doZVeto   = 0;
        events_LepInZMass0BJetsJet2MET40[j].cutNJets  = 2;
        events_LepInZMass0BJetsJet2MET40[j].cutNBJets = -1;
        events_LepInZMass0BJetsJet2MET40[j].cutMET    = 40.;
        
        events_LepOutZMass0BJetsJet2MET40[j].doZVeto   = 1;
        events_LepOutZMass0BJetsJet2MET40[j].cutNJets  = 2;
        events_LepOutZMass0BJetsJet2MET40[j].cutNBJets = -1;
        events_LepOutZMass0BJetsJet2MET40[j].cutMET    = 40.;
        
    }    
    events_LepOutZMassJet2BJet1[2].blindDataChannel = 1; // blind emu with ZVeto + jet cuts
    /*
     int     whichdiLepType; // -1: inclusive, 0: MuMu, 1: EE, 2: EMu
     int     doZVeto;   // -1: inclusive, 0: ZMass window, 1: outside ZMass window;
     int     cutNJets;  // -1: inclusive, # > 0: require NJets >= #
     int     cutNBJets;  // -1: inclusive, # > 0: require NBJets >= #
     float   cutMET; // # > 0: require MET >= #
     */
    
    SampleT events_NJetsGeq2; events_NJetsGeq2.histNameSuffix = "_NJetsGeq2";
    events_NJetsGeq2.histXaxisSuffix = ""; events_NJetsGeq2.histYaxisSuffix = "";  events_NJetsGeq2.histZaxisSuffix = "";
    events_NJetsGeq2.whichdiLepType = -1; events_NJetsGeq2.doZVeto = -1; events_NJetsGeq2.cutNJets = 2; events_NJetsGeq2.cutNBJets = -1; events_NJetsGeq2.cutMET = 0.;
    events_NJetsGeq2.blindDataChannel = 0;
    
    SampleT events_NBJetsGeq1; events_NBJetsGeq1.histNameSuffix = "_NBJetsGeq1";
    events_NBJetsGeq1.histXaxisSuffix = ""; events_NBJetsGeq1.histYaxisSuffix = "";  events_NBJetsGeq1.histZaxisSuffix = "";
    events_NBJetsGeq1.whichdiLepType = -1; events_NBJetsGeq1.doZVeto = -1; events_NBJetsGeq1.cutNJets = -1; events_NBJetsGeq1.cutNBJets = 1; events_NBJetsGeq1.cutMET = 0.;
    events_NBJetsGeq1.blindDataChannel = 0;
    
    SampleT events_NBJetsGeq2; events_NBJetsGeq2.histNameSuffix = "_NBJetsGeq2";
    events_NBJetsGeq2.histXaxisSuffix = ""; events_NBJetsGeq2.histYaxisSuffix = "";  events_NBJetsGeq2.histZaxisSuffix = "";
    events_NBJetsGeq2.whichdiLepType = -1; events_NBJetsGeq2.doZVeto = -1; events_NBJetsGeq2.cutNJets = -1; events_NBJetsGeq2.cutNBJets = 2; events_NBJetsGeq2.cutMET = 0.;
    events_NBJetsGeq2.blindDataChannel = 0;
    
    SampleT events_NJetsGeq2_NBJetsGeq1; events_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_Jet2BJet1"; events_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = ""; events_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = ""; events_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = "";
    events_NJetsGeq2_NBJetsGeq1.whichdiLepType = -1; events_NJetsGeq2_NBJetsGeq1.doZVeto = -1; events_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_NJetsGeq2_NBJetsGeq1.cutMET = 0.;
    events_NJetsGeq2_NBJetsGeq1.blindDataChannel = 0;
    
    SampleT events_FullCut; events_FullCut.histNameSuffix = "_FullCut"; events_FullCut.histXaxisSuffix = ""; 
    events_FullCut.histYaxisSuffix = ""; events_FullCut.histZaxisSuffix = "";
    events_FullCut.whichdiLepType = -1; events_FullCut.doZVeto = 1; events_FullCut.cutNJets = 2; events_FullCut.cutNBJets = 1; events_FullCut.cutMET = 40.;
    events_FullCut.blindDataChannel = 1;    
    //"Inclusive subsample"//
    SampleT allEvents; allEvents.histNameSuffix = "_inclusive"; 
    allEvents.histXaxisSuffix = ""; allEvents.histYaxisSuffix = ""; allEvents.histZaxisSuffix = "";
    allEvents.whichdiLepType = -1; allEvents.doZVeto = -1; allEvents.cutNJets = -1; allEvents.cutNBJets = -1; allEvents.cutMET = 0.; allEvents.blindDataChannel = 0;
    
    //inclusive in ZMass window, w/ and w/out jet cuts    
    SampleT events_inZMass; events_inZMass.histNameSuffix = "_allLep_inZ";
    events_inZMass.histXaxisSuffix = ""; events_inZMass.histYaxisSuffix = ""; events_inZMass.histZaxisSuffix = "";
    events_inZMass.whichdiLepType = -1; events_inZMass.doZVeto = 0; events_inZMass.cutNJets = -1; events_inZMass.cutNBJets = -1; events_inZMass.cutMET = 0.;
    events_inZMass.blindDataChannel = 0;
    
    SampleT events_inZMassJet2BJet1; events_inZMassJet2BJet1.histNameSuffix = "_allLep_inZ_w_JetCuts";
    events_inZMassJet2BJet1.histXaxisSuffix = ""; events_inZMassJet2BJet1.histYaxisSuffix = ""; events_inZMassJet2BJet1.histZaxisSuffix = "";
    events_inZMassJet2BJet1.whichdiLepType = -1; events_inZMassJet2BJet1.doZVeto = 0; events_inZMassJet2BJet1.cutNJets = 2; events_inZMassJet2BJet1.cutNBJets = 1; events_inZMassJet2BJet1.cutMET = 0.;
    events_inZMassJet2BJet1.blindDataChannel = 0;
    
    SampleT events_inZMassMET40; events_inZMassMET40.histNameSuffix = "_allLep_inZ_w_METCut";
    events_inZMassMET40.histXaxisSuffix = ""; events_inZMassMET40.histYaxisSuffix = ""; events_inZMassMET40.histZaxisSuffix = "";
    events_inZMassMET40.whichdiLepType = -1; events_inZMassMET40.doZVeto = 0; events_inZMassMET40.cutNJets = -1; events_inZMassMET40.cutNBJets = -1; events_inZMassMET40.cutMET = 40.;
    events_inZMassMET40.blindDataChannel = 0;
    
    SampleT events_inZMassJet2BJet1MET40; events_inZMassJet2BJet1MET40.histNameSuffix = "_allLep_inZ_w_JetCuts_w_METCut";
    events_inZMassJet2BJet1MET40.histXaxisSuffix = ""; events_inZMassJet2BJet1MET40.histYaxisSuffix = ""; events_inZMassJet2BJet1MET40.histZaxisSuffix = "";
    events_inZMassJet2BJet1MET40.whichdiLepType = -1; events_inZMassJet2BJet1MET40.doZVeto = 0; events_inZMassJet2BJet1MET40.cutNJets = 2; events_inZMassJet2BJet1MET40.cutNBJets = 1; events_inZMassJet2BJet1MET40.cutMET = 40.;
    events_inZMassJet2BJet1MET40.blindDataChannel = 0;
    
    ///push Sample_Ts into a vector    
    vector<SampleT> * subSampVec = new vector<SampleT>;
    for (int k = 0; k < 3; ++k) {
        subSampVec->push_back(events_LepInZMass[k]); subSampVec->push_back(events_LepOutZMass[k]);
        subSampVec->push_back(events_LepInZMassJet2[k]); subSampVec->push_back(events_LepOutZMassJet2[k]);
        subSampVec->push_back(events_LepInZMassJet2BJet1[k]); subSampVec->push_back(events_LepOutZMassJet2BJet1[k]);
        subSampVec->push_back(events_LepInZMassJet2BJet1MET40[k]); subSampVec->push_back(events_LepOutZMassJet2BJet1MET40[k]);
        
        subSampVec->push_back(events_LepInZMassBothinBarrel[k]);
        subSampVec->push_back(events_LepInZMassOneinBarrel[k]);
        subSampVec->push_back(events_LepInZMassBothinEndcap[k]);
        
        subSampVec->push_back(events_LepInZMassor0BJetsJet2[k]);
        subSampVec->push_back(events_LepInZMassor0BJetsJet2MET40[k]);
        subSampVec->push_back(events_LepInZMass0BJetsJet2[k]);
        subSampVec->push_back(events_LepInZMass0BJetsJet2MET40[k]);
        subSampVec->push_back(events_LepOutZMass0BJetsJet2[k]);
        subSampVec->push_back(events_LepOutZMass0BJetsJet2MET40[k]);        
    }
    subSampVec->push_back(events_NJetsGeq2); subSampVec->push_back(events_NBJetsGeq1); subSampVec->push_back(events_NBJetsGeq2); subSampVec->push_back(events_NJetsGeq2_NBJetsGeq1);
    subSampVec->push_back(events_FullCut);
    subSampVec->push_back(allEvents);
    subSampVec->push_back(events_inZMass); subSampVec->push_back(events_inZMassJet2BJet1);
    subSampVec->push_back(events_inZMassMET40); subSampVec->push_back(events_inZMassJet2BJet1MET40);
    return subSampVec;
}

inline vector<SystT> * SystVec() {
    //    int     whichSystType;   // 0 = universal systematic, 1 = lepton systematic, 2 = jet systematic, 3 = other systematic
    /*
     systematics considered: MT2ll, Lepton trig/iso SF, Jet ES, Jet ER, Lepton ES, Lepton ER(??)
     */
    // For MT2ll, really should breakdown into components....for now doing agglomerate
    SystT MT2llShiftUp; MT2llShiftUp.name = "_MT2llShiftUp"; MT2llShiftUp.systVarKey = "_MT2llShiftUp";
    MT2llShiftUp.whichSystType = 3;
    /*
    SystT MT2llShiftDown; MT2llShiftDown.name = "_MT2llShiftDown"; MT2llShiftDown.systVarKey = "_MT2llShiftDown";
    MT2llShiftDown.whichSystType = 3;
     */
    
    SystT genTopReweight; genTopReweight.name = "_genTopRW"; genTopReweight.systVarKey = "";
    genTopReweight.whichSystType = 3;
    
    SystT genStopXSecShiftUp; genStopXSecShiftUp.name = "_genStopXSecShiftUp"; genStopXSecShiftUp.systVarKey = "";
    genStopXSecShiftUp.whichSystType = 3;

    SystT genStopXSecShiftDown; genStopXSecShiftDown.name = "_genStopXSecShiftDown"; genStopXSecShiftDown.systVarKey = "";
    genStopXSecShiftDown.whichSystType = 3;
    
    SystT LepEffSFShiftUp; LepEffSFShiftUp.name = "_LepEffSFShiftUp"; LepEffSFShiftUp.systVarKey = ""; //LepEffSFShiftUp.systVarKey = "_LepEffShiftUp";
    LepEffSFShiftUp.whichSystType = 0;
    SystT LepEffSFShiftDown; LepEffSFShiftDown.name = "_LepEffSFShiftDown"; LepEffSFShiftDown.systVarKey = ""; // LepEffSFShiftDown.systVarKey = "_LepEffShiftDown";
    LepEffSFShiftDown.whichSystType = 0;
    SystT LepESShiftUp; LepESShiftUp.name = "_LepESShiftUp"; LepESShiftUp.systVarKey = "_LepESShiftUp";
    LepESShiftUp.whichSystType = 1;
    SystT LepESShiftDown; LepESShiftDown.name = "_LepESShiftDown"; LepESShiftDown.systVarKey = "_LepESShiftDown";
    LepESShiftDown.whichSystType = 1;
    
    /*
     SystT LepERShiftUp; LepERShiftUp.name = "_LepERShiftUp"; LepERShiftUp.systVarKey = "_LepERShiftUp";
     SystT LepERShiftDown; LepERShiftDown.name = "_LepERShiftDown"; LepERShiftDown.systVarKey = "_LepERShiftDown";
     SystT JetESShiftUp; JetESShiftUp.name = "_JetESShiftUp"; JetESShiftUp.systVarKey = "_JetESShiftUp";
     SystT JetESShiftDown; JetESShiftDown.name = "_JetESShiftDown"; JetESShiftDown.systVarKey = "_JetESShiftUpDown";
     SystT JetERShiftUp; JetERShiftUp.name = "_JetERShiftUp"; JetERShiftUp.systVarKey = "_JetERShiftUp";
     SystT JetERShiftDown; JetERShiftDown.name = "_JetERShiftDown"; JetERShiftDown.systVarKey = "_JetERShiftDown";
     */
    vector<SystT> * systVec = new vector<SystT>;
    systVec->push_back(MT2llShiftUp);// systVec->push_back(MT2llShiftDown);
    systVec->push_back(LepEffSFShiftUp); systVec->push_back(LepEffSFShiftDown);
    systVec->push_back(LepESShiftUp); systVec->push_back(LepESShiftDown);
    systVec->push_back(genTopReweight);
    systVec->push_back(genStopXSecShiftUp); systVec->push_back(genStopXSecShiftDown);
    //    systVec->push_back(LepERShiftUp); systVec->push_back(LepERShiftDown);
    //    systVec->push_back(JetESShiftUp); systVec->push_back(JetESShiftDown);
    //    systVec->push_back(JetERShiftUp); systVec->push_back(JetERShiftDown);
    return systVec;
}
inline TString DescriptorString(SampleT inputSubSamp) {
    //descriptor strings
    TString baseDescString[4] = {"All three flavors of lepton events, inclusive other than being/requiring ",  "All mumu events, inclusive other than being/requiring", "All ee events, inclusive other than being/requiring", "All emu events, inclusive other than being/requiring"};
    TString dsZVeto = ", diLepton system invariant mass outside of Z Mass window (76 GeV:106 GeV)";
    TString dsAntiZVeto = ", diLepton system invariant mass inside of Z Mass window (76 GeV:106 GeV) ";
    TString dsJetCutStringP1 = ", at least ";
    TString dsJetCutStringP2 = " jet(s)";
    TString dsBJetCutStringP2 = " b jet(s)";
    TString dsMETCutStringP1 = ", at least ";
    TString dsMETCutStringP2 = " GeV of MET";
    TString etaCutString[3] = {", both leptons in endcap eta range", ", one of the two leptons in barrel eta range", ", both leptons in barrel eta range"};
    //descriptor strings   
    TString outString;
    outString += baseDescString[inputSubSamp.whichdiLepType + 1];
    switch (inputSubSamp.doZVeto) {
        case 0:
            outString += dsAntiZVeto;
            break;
        case 1:
            outString += dsZVeto;
            break;            
        default:
            break;
    }
    if (inputSubSamp.cutNJets > 0) {
        outString += dsJetCutStringP1;
        outString += inputSubSamp.cutNJets;
        outString += dsJetCutStringP2;
    }
    if (inputSubSamp.cutNBJets > 0) {
        outString += dsJetCutStringP1;
        outString += inputSubSamp.cutNBJets;
        outString += dsBJetCutStringP2;
    }
    if (inputSubSamp.cutMET > 0) {
        outString += dsMETCutStringP1;
        outString += inputSubSamp.cutMET;
        outString += dsMETCutStringP2;
    }
    if (inputSubSamp.histNameSuffix.Contains("BothinEndcap")) outString += etaCutString[0];
    if (inputSubSamp.histNameSuffix.Contains("OneinBarrel")) outString += etaCutString[1];
    if (inputSubSamp.histNameSuffix.Contains("BothinBarrel")) outString += etaCutString[2];
    if (inputSubSamp.histNameSuffix.Contains("0BJets")) outString += TString(", with NO BJets");
    if (inputSubSamp.blindDataChannel) outString += TString(", note blinded channel");
    return outString;
}
inline vector<SpecHistBinT> * SpecHistBinVec() {
    SpecHistBinT ZMass; 
    ZMass.HistTXLabel = "M_{ll}";
    ZMass.SampTLabel = "inZ";
    vector<SpecHistBinT> * specHistBinVec = new vector<SpecHistBinT>;
    specHistBinVec->push_back(ZMass); 
    return specHistBinVec;
}


inline TString HistProjection1D(vector<TH2F *> * inputHistVec, vector<TH1F *> * outputHistVec, TString plotName, int whichCase) {
    const double PI = 3.14159265;
    TString outString = "Nothin";
    TH1F * currOutHist;
    TH2F * currInHist;
    TAxis * XAxis, * YAxis;
    int XAxisLB, XAxisUB, YAxisLB, YAxisUB;
    TString outHistName;
    TString outHistPlusName;
    double YAxisLBFind, YAxisUBFind;
    if (plotName.Contains("h_MT2ll_vs_DeltaPhiLep0Lep1")) {
        switch (whichCase) {
            case 0: // plot for DeltaPhiZMET > 2./3. Pi
                outHistPlusName = "_XAxis_DPhiGt2/3_";
                YAxisLBFind = 2./3. * PI;
                YAxisUBFind = PI;
                outString = "#Delta #phi_{ll} #in {2/3#pi:#pi}";
                break;
            case 1: // plot for 1./3. Pi < DeltaPhiZMET < 2/3 Pi
                outHistPlusName = "_XAxis_DPhiGt1/3Lt2/3_";
                YAxisLBFind = 1./3. * PI;
                YAxisUBFind = 2./3. * PI;
                outString = "#Delta #phi_{ll} #in {1/3#pi:2/3#pi}";
                break;   
            case 2: // plot for DeltaPhiZMET < 1/3 Pi
                outHistPlusName = "_XAxis_DPhiLt1/3_";
                YAxisLBFind = 0.;
                YAxisUBFind = 1./3. * PI;
                outString = "#Delta #phi_{ll} #in {0:1/3#pi}";
                break; 
            default:
                break;
        }
    }
    if (outString.Contains("Nothin")) return outString;
    for (unsigned int iHist = 0; iHist < inputHistVec->size(); ++iHist) {
        currInHist = (TH2F*) inputHistVec->at(iHist);
        outHistName = plotName;
        outHistName += outHistPlusName;
        outHistName += iHist;
        YAxis = currInHist->GetYaxis();
        YAxisLB = YAxis->FindBin(YAxisLBFind);
        YAxisUB = YAxis->FindBin(YAxisUBFind);
        cout << "case " << whichCase << endl;
        cout << "YAxisLBFind " << YAxisLBFind << endl;
        cout << "YAxisUBFind " << YAxisUBFind << endl;
        cout << "YAxisLB " << YAxisLB << endl;
        cout << "YAxisUB " << YAxisUB << endl;
        currOutHist = (TH1F*) currInHist->ProjectionX(outHistName, YAxisLB, YAxisUB, "e");
        outputHistVec->push_back(currOutHist);        
    }
    return outString;
}

inline TString HistProjection1D(vector<TH3F *> * inputHistVec, vector<TH1F *> * outputHistVec, TString plotName, int whichCase) {
    const double PI = 3.14159265;
    TString outString = "Nothin";
    TH1F * currOutHist;
    TH3F * currInHist;
    TString addName = "";
    TAxis * XAxis, * YAxis, * ZAxis;
    int XAxisLB, XAxisUB, YAxisLB, YAxisUB, ZAxisLB, ZAxisUB;
    TString outHistName;
    TString outHistPlusName;
    float YAxisLBFind, YAxisUBFind, ZAxisLBFind, ZAxisUBFind;
    if (plotName.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_nVtx")) {
        switch (whichCase) {
            case 0: // plot for DeltaPhiZMET > 2./3. Pi and nVtx = 1-10
                outHistPlusName = "_XAxis_DPhiZMETGt2/3_nVtx1to10_";
                outHistPlusName += addName; 
                YAxisLBFind = 2./3. * PI;
                YAxisUBFind = PI;
                ZAxisLBFind = 1.;
                ZAxisUBFind = 10.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {2/3#pi:#pi}, N_{vtx}^{reco} #in {1:10}";
                break;
            case 1: // plot for 1./3. Pi < DeltaPhiZMET < 2/3 Pi and nVtx = 1-10
                outHistPlusName = "_XAxis_DPhiZMETGt1/3Lt2/3_nVtx1to10_";
                outHistPlusName += addName; 
                YAxisLBFind = 1./3. * PI;
                YAxisUBFind = 2./3. * PI;
                ZAxisLBFind = 1.;
                ZAxisUBFind = 10.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {1/3#pi:2/3 #pi}, N_{vtx}^{reco} #in {1:10}";
                break;   
            case 2: // plot for DeltaPhiZMET < 1/3 Pi and nVtx = 1-10
                outHistPlusName = "_XAxis_DPhiZMETLt1/3_nVtx1to10_";
                outHistPlusName += addName; 
                YAxisLBFind = 0;
                YAxisUBFind = 1./3. * PI;
                ZAxisLBFind = 1.;
                ZAxisUBFind = 10.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {0:1/3#pi}, N_{vtx}^{reco} #in {1:10}";
                break; 
            case 3: // plot for DeltaPhiZMET > 2./3. Pi and nVtx = 11-20
                outHistPlusName = "_XAxis_DPhiZMETGt2/3_nVtx11to20_";
                outHistPlusName += addName; 
                YAxisLBFind = 2./3. * PI;
                YAxisUBFind = PI;
                ZAxisLBFind = 11.;
                ZAxisUBFind = 20.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {2/3#pi:#pi}, N_{vtx}^{reco} #in {11:20}";
                break;
            case 4: // plot for 1./3. Pi < DeltaPhiZMET < 2/3 Pi and nVtx = 11-20
                outHistPlusName = "_XAxis_DPhiZMETGt1/3Lt2/3_nVtx11to20_";
                outHistPlusName += addName; 
                YAxisLBFind = 1./3. * PI;
                YAxisUBFind = 2./3. * PI;
                ZAxisLBFind = 11.;
                ZAxisUBFind = 20.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {1/3#pi:2/3 #pi}, N_{vtx}^{reco} #in {11:20}";
                break;   
            case 5: // plot for DeltaPhiZMET < 1/3 Pi and nVtx = 11-20
                outHistPlusName = "_XAxis_DPhiZMETLt1/3_nVtx11to20_";
                outHistPlusName += addName; 
                YAxisLBFind = 0;
                YAxisUBFind = 1./3. * PI;
                ZAxisLBFind = 11.;
                ZAxisUBFind = 20.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {0:1/3#pi}, N_{vtx}^{reco} #in {11:20}";
                break; 
            case 6: // plot for DeltaPhiZMET > 2./3. Pi and nVtx = 21-30
                outHistPlusName = "_XAxis_DPhiZMETGt2/3_nVtx21to30_";
                outHistPlusName += addName; 
                YAxisLBFind = 2./3. * PI;
                YAxisUBFind = PI;
                ZAxisLBFind = 21.;
                ZAxisUBFind = 30.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {2/3#pi:#pi}, N_{vtx}^{reco} #in {21:30}";
                break;
            case 7: // plot for 1./3. Pi < DeltaPhiZMET < 2/3 Pi and nVtx = 21-30
                outHistPlusName = "_XAxis_DPhiZMETGt1/3Lt2/3_nVtx21to30_";
                outHistPlusName += addName; 
                YAxisLBFind = 1./3. * PI;
                YAxisUBFind = 2./3. * PI;
                ZAxisLBFind = 21.;
                ZAxisUBFind = 30.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {1/3#pi:2/3 #pi}, N_{vtx}^{reco} #in {21:30}";
                break;   
            case 8: // plot for DeltaPhiZMET < 1/3 Pi and nVtx = 21-30
                outHistPlusName = "_XAxis_DPhiZMETLt1/3_nVtx21to30_";
                outHistPlusName += addName; 
                YAxisLBFind = 0;
                YAxisUBFind = 1./3. * PI;
                ZAxisLBFind = 21.;
                ZAxisUBFind = 30.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {0:1/3#pi}, N_{vtx}^{reco} #in {21:30}";
                break;
            default:
                break;
        }
    }    
    else if (plotName.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets")) {
        if (plotName.Contains("nVtx1to10")) {
            addName = "nVtx1to10";
        }
        else if (plotName.Contains("nVtx11to20")) {
            addName = "nVtx11to20";
        }
        else if (plotName.Contains("nVtx21to30")) {
            addName = "nVtx21to30";
        }
        switch (whichCase) {
            case 0: // plot for DeltaPhiZMET > 2./3. Pi and N_{Jets} #in {0}
                outHistPlusName = "_XAxis_DPhiZMETGt2/3_nJets0_";
                outHistPlusName += addName; 
                YAxisLBFind = 2./3. * PI;
                YAxisUBFind = PI;
                ZAxisLBFind = 0.;
                ZAxisUBFind = 0.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {2/3#pi:#pi}, N_{Jets} #in {0}";
                break;
            case 1: // plot for 1./3. Pi < DeltaPhiZMET < 2/3 Pi and N_{Jets} #in {0}
                outHistPlusName = "_XAxis_DPhiZMETGt1/3Lt2/3_nJets0_";
                outHistPlusName += addName; 
                YAxisLBFind = 1./3. * PI;
                YAxisUBFind = 2./3. * PI;
                ZAxisLBFind = 0.;
                ZAxisUBFind = 0.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {1/3#pi:2/3 #pi}, N_{Jets} #in {0}";
                break;   
            case 2: // plot for DeltaPhiZMET < 1/3 Pi and N_{Jets} #in {0}
                outHistPlusName = "_XAxis_DPhiZMETLt1/3_nJets0_";
                outHistPlusName += addName; 
                YAxisLBFind = 0;
                YAxisUBFind = 1./3. * PI;
                ZAxisLBFind = 0.;
                ZAxisUBFind = 0.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {0:1/3#pi}, N_{Jets} #in {0}";
                break; 
            case 3: // plot for DeltaPhiZMET > 2./3. Pi and N_{Jets} #in {1}
                outHistPlusName = "_XAxis_DPhiZMETGt2/3_nJets1_";
                outHistPlusName += addName; 
                YAxisLBFind = 2./3. * PI;
                YAxisUBFind = PI;
                ZAxisLBFind = 1.;
                ZAxisUBFind = 1.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {2/3#pi:#pi}, N_{Jets} #in {1}";
                break;
            case 4: // plot for 1./3. Pi < DeltaPhiZMET < 2/3 Pi and N_{Jets} #in {1}
                outHistPlusName = "_XAxis_DPhiZMETGt1/3Lt2/3_nJets1_";
                outHistPlusName += addName; 
                YAxisLBFind = 1./3. * PI;
                YAxisUBFind = 2./3. * PI;
                ZAxisLBFind = 1.;
                ZAxisUBFind = 1.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {1/3#pi:2/3 #pi}, N_{Jets} #in {1}";
                break;   
            case 5: // plot for DeltaPhiZMET < 1/3 Pi and N_{Jets} #in {1}
                outHistPlusName = "_XAxis_DPhiZMETLt1/3_nJets1_";
                outHistPlusName += addName; 
                YAxisLBFind = 0;
                YAxisUBFind = 1./3. * PI;
                ZAxisLBFind = 1.;
                ZAxisUBFind = 1.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {0:1/3#pi}, N_{Jets} #in {1}";
                break; 
            case 6: // plot for DeltaPhiZMET > 2./3. Pi and N_{Jets} #in {2:#infty}
                outHistPlusName = "_XAxis_DPhiZMETGt2/3_nJetsGt1_";
                outHistPlusName += addName; 
                YAxisLBFind = 2./3. * PI;
                YAxisUBFind = PI;
                ZAxisLBFind = 2.;
                ZAxisUBFind = 20.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {2/3#pi:#pi}, N_{Jets} #in {2:#infty}";
                break;
            case 7: // plot for 1./3. Pi < DeltaPhiZMET < 2/3 Pi and N_{Jets} #in {2:#infty}
                outHistPlusName = "_XAxis_DPhiZMETGt1/3Lt2/3_nJetsGt1_";
                outHistPlusName += addName; 
                YAxisLBFind = 1./3. * PI;
                YAxisUBFind = 2./3. * PI;
                ZAxisLBFind = 2.;
                ZAxisUBFind = 20.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {1/3#pi:2/3 #pi}, N_{Jets} #in {2:#infty}";
                break;   
            case 8: // plot for DeltaPhiZMET < 1/3 Pi and N_{Jets} #in {2:#infty}
                outHistPlusName = "_XAxis_DPhiZMETLt1/3_nJetsGt1_";
                outHistPlusName += addName; 
                YAxisLBFind = 0;
                YAxisUBFind = 1./3. * PI;
                ZAxisLBFind = 2.;
                ZAxisUBFind = 20.;
                outString = "#Delta #phi_{Z, #slash{E}_{T}} #in {0:1/3#pi}, N_{Jets} #in {2:#infty}";
                break;
        }
    }
    if (outString.Contains("Nothin")) return outString;
    for (unsigned int iHist = 0; iHist < inputHistVec->size(); ++iHist) {
        currInHist = (TH3F*) inputHistVec->at(iHist);
        outHistName = plotName;
        outHistName += outHistPlusName;
        outHistName += iHist;
        YAxis = currInHist->GetYaxis();
        YAxisLB = YAxis->FindBin(YAxisLBFind);
        YAxisUB = YAxis->FindBin(YAxisUBFind);
        ZAxis = currInHist->GetZaxis();
        ZAxisLB = ZAxis->FindBin(ZAxisLBFind);
        ZAxisUB = ZAxis->FindBin(ZAxisUBFind);
        currOutHist = (TH1F*) currInHist->ProjectionX(outHistName, YAxisLB, YAxisUB, ZAxisLB, ZAxisUB, "e");
        outputHistVec->push_back(currOutHist);        
    }
    return outString;
}



/*
 float * RebinArray() {
 float
 }
 */

#endif //HistSampFunc_h_
