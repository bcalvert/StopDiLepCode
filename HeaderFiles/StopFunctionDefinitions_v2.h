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
#include <Math/GenVector/LorentzVector.h>
//BTag shenanigans
#include "./BTagSFUtil.C"


#include <iostream>
#include <fstream>
//#include <vector>
#include <cmath>
#include <sstream>
#include <map>
using namespace std;



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
inline vector<TH1F *> * OneDProjectionReturnVec(TH1 * inputHist, int numDims, int whichAxisToProjTo, int whichAxisForDist, int axisProjLB, int axisProjUB, TString nameBase) {
    TString name;
    vector<TH1F *> * outHistVector = new vector<TH1F *>;
    TH1F * projHist;
    TH2F * projPatsy2DHist;
    TH3F * projPatsy3DHist;
    TString patsyNamebase = "patsy";
    TString patsyName;
    int NBins;
    int whichAxisToMix = 6 - (whichAxisForDist + whichAxisToProjTo);
    int axis1LB, axis1UB, axis2LB, axis2UB;
    switch (whichAxisForDist) {
        case 1:
            NBins = inputHist->GetNbinsX();
            break;
        case 2:
            NBins = inputHist->GetNbinsY();
            break;
        case 3:
            NBins = inputHist->GetNbinsZ();
            break;
        default:
            break;
    }    
    for (int ib = 1; ib <= NBins; ++ib) {
        axis1LB = (whichAxisToMix > whichAxisForDist) ? ib : axisProjLB;
        axis1UB = (whichAxisToMix > whichAxisForDist) ? ib : axisProjUB;
        axis2LB = (whichAxisToMix > whichAxisForDist) ? axisProjLB : ib;
        axis2UB = (whichAxisToMix > whichAxisForDist) ? axisProjUB : ib;
        name = nameBase;
        name += ib;
        patsyName = patsyNamebase;
        patsyName += name;
        switch (whichAxisToProjTo) {
            case 1:
                if (numDims > 2) {
                    projPatsy3DHist = (TH3F *) inputHist->Clone(patsyName);
                    projHist = (TH1F *) projPatsy3DHist->ProjectionX(name, axis1LB, axis1UB, axis2LB, axis2UB);   
                }
                else {
                    projPatsy2DHist = (TH2F *) inputHist->Clone(patsyName);
                    
                    projHist = (TH1F *) projPatsy2DHist->ProjectionX(name, axis1LB, axis1UB);   
                }
                break;
            case 2:
                if (numDims > 2) {
                    projPatsy3DHist = (TH3F *) inputHist->Clone(patsyName);
                    projHist = (TH1F *) projPatsy3DHist->ProjectionY(name, axis1LB, axis1UB, axis2LB, axis2UB);   
                }
                else {
                    projPatsy2DHist = (TH2F *) inputHist->Clone(patsyName);
                    
                    projHist = (TH1F *) projPatsy2DHist->ProjectionY(name, axis1LB, axis1UB);   
                }
                break;
            case 3:
                if (numDims > 2) {
                    projPatsy3DHist = (TH3F *) inputHist->Clone(patsyName);
                    projHist = (TH1F *) projPatsy3DHist->ProjectionZ(name, axis1LB, axis1UB, axis2LB, axis2UB);   
                }
                break;
            default:
                break;
        }
        outHistVector->push_back(projHist);
    }
    return outHistVector;
}


inline TLorentzVector LeptonScaleSystShift(TLorentzVector inputLepVec, int inputLepPDGID, float shiftDirection) {
    float barrelEtaEnd = 1.4442; float endcapEtaStart = 1.566; 
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

/*
 inline TLorentzVector JetScaleSystShift(TLorentzVector inputJetVec) {
 
 }
 */

inline bool FilterMET(vector<bool> * regFilterVec, vector<bool> * oppFilterVec) {
    for (unsigned int iReg = 0; iReg < regFilterVec->size(); ++iReg) {
        if (!regFilterVec->at(iReg)) return false;
    }
    for (unsigned int iOpp = 0; iOpp < oppFilterVec->size(); ++iOpp) {
        if (oppFilterVec->at(iOpp)) return false;
    }
    return true;
}
inline Double_t ElectronEffectiveArea(float elecEta) {
    //should check numbers here:
    //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.3&view=markup
    Double_t Aeff = 0.;    
    if (fabs(elecEta) < 1.0 )         Aeff = 0.13; // +/- 0.001
    else if (fabs(elecEta) < 1.479)    Aeff = 0.14; // +/- 0.002
    else if (fabs(elecEta) < 2.0)      Aeff = 0.07; // +/- 0.001
    else if (fabs(elecEta) < 2.2)      Aeff = 0.09; // +/- 0.001
    else if (fabs(elecEta) < 2.3)      Aeff = 0.11; // +/- 0.002
    else if (fabs(elecEta) < 2.4)      Aeff = 0.11; // +/- 0.003
    else if (fabs(elecEta) > 2.4)      Aeff = 0.14; // +/- 0.004
    
    return Aeff;
}
inline bool ElectronPassCutStage1(bool isEB, bool isEE, float SCEta, float deltaPhiIn, float deltaEtaIn, float sigIetaIeta, float HtoERatio) {
    float EBDeltaPhiInCut = 0.06, EEDeltaPhiInCut = 0.03, deltaPhiInCut;
    float EBDeltaEtaInCut = 0.004, EEDeltaEtaInCut = 0.007, deltaEtaInCut;
    float EBSigIetaIetaCut = 0.01, EESigIetaIetaCut = 0.03, sigIetaIetaCut;
    float EBHtoERatio = 0.12, EEHtoERatio = 0.10, HtoERatioCut;
    float SCEtaMax = 2.5, SCEtaBarrelCut = 1.4442, SCEtaEndcapCut = 1.566;
    bool passCut = false;
    if (fabs(SCEta) >= SCEtaMax || ((fabs(SCEta) >= SCEtaBarrelCut ) && (fabs(SCEta) <= SCEtaEndcapCut))) {
        passCut = false;
        return passCut;
    }
    if (isEB) {
        deltaPhiInCut  = EBDeltaPhiInCut;
        deltaEtaInCut  = EBDeltaEtaInCut;
        sigIetaIetaCut = EBSigIetaIetaCut;
        HtoERatioCut   = EBHtoERatio;
    }
    else if (isEE) {
        deltaPhiInCut  = EEDeltaPhiInCut;
        deltaEtaInCut  = EEDeltaEtaInCut;
        sigIetaIetaCut = EESigIetaIetaCut;
        HtoERatioCut   = EEHtoERatio;
    }
    else {
        passCut = false;
        return passCut;
    }
    passCut = ((fabs(deltaPhiIn) < deltaPhiInCut) && (fabs(deltaEtaIn) < deltaEtaInCut) && (sigIetaIeta < sigIetaIetaCut) && (HtoERatio < HtoERatioCut));
    return passCut;
}
inline bool ElectronPassCutStage2(bool isPFElectron, bool passConvVeto, float relPFElecIso, float elecIP, float elecDZ, float elecE, float elecEoverP, float elecMissHits) {
    bool passCut;
    float elecIsoRatioCut = 0.15;
    float elecEPCut = 0.05, elecDZCut = 0.1, elecIPCut = 0.02;
    int elecNumMissHitsCut = 1;
    if (!isPFElectron || !passConvVeto || relPFElecIso >= elecIsoRatioCut) {
        passCut = false;
        return passCut;
    }
    passCut = ((fabs(elecIP) < elecIPCut) && (fabs(elecDZ) < elecDZCut) && (fabs((1 - elecEoverP)/elecE) < elecEPCut) && (elecMissHits <= elecNumMissHitsCut));
    return passCut;
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification    
}
inline bool ElectronPassCutDESY(float relPFElecIso, float elecEta) {
    bool passCut;        
    float elecIsoRatioCut = 0.15; float elecEtaUB = 2.5;
    float barrelEtaEnd = 1.4442; float endcapEtaStart = 1.566;
    if (relPFElecIso > elecIsoRatioCut || fabs(elecEta) > elecEtaUB || (fabs(elecEta) > barrelEtaEnd && fabs(elecEta) < endcapEtaStart)) {
        passCut = false;
    }
    else {
        passCut = true;
    }
    return passCut;
}
//inline void ElectronPickOvi(vector<float> * ElecPx, vector<float> * ElecPy, vector<float> * ElecPz, vector<float> * ElecEn, vector<int> * ElecCharge, vector<float> * PFElecPt, vector<float> * ElecSCEta, vector<float> * ElecDeltaPhiIn, vector<float> * ElecDeltaEtaIn, vector<float> * ElecSigIetaIeta, vector<float> * ElecHtoERatio, vector<float> * ElecIP, vector<float> * ElecDZ, vector<float> * ElecECalE, vector<float> * ElecSCEOverP, vector<int> * ElecNumMissHits, vector<float> * ElecPFNeutHadIso, vector<float> * ElecPFCharHadIso, vector<float> * ElecPFPhotIso, vector<bool> * passConvVeto, vector<bool> * isPFElectron, vector<bool> * ElecisEB, vector<bool> * ElecisEE, float eventRhoIso, vector<Lepton> * vecIsoLeptons, int levelVerbosity) {
inline void ElectronPickOvi(ElectronEventPointers inEEPs, float eventRhoIso, vector<Lepton> * vecIsoLeptons, int levelVerbosity) {
    bool currElecPassCut;
    float currElecPt;
    float currRelPFElecIso;
    float currElecECalE, currElecSCEOverP, currElecSCEta;
    float currElecIP, currElecDZ;
    float currElecDeltaPhiIn, currElecDeltaEtaIn, currElecSigIetaIeta, currElecHtoERatio;
    int currElecNumMissHits;
    Lepton patsyLep;
    TLorentzVector patsyVec;
//    for (unsigned int i = 0; i < ElecPx->size(); ++i) {
    for (unsigned int i = 0; i < inEEPs.numElectrons; ++i) {
        patsyVec.SetPxPyPzE(inEEPs.ElecPx->at(i), inEEPs.ElecPy->at(i), inEEPs.ElecPz->at(i), inEEPs.ElecEn->at(i));        
        if (levelVerbosity > 0) {
            cout << "electron " << i << " has:" << endl;
            cout << "pT is " << patsyVec.Pt() << endl;
            cout << "Eta is " << patsyVec.Eta() << endl;
            cout << "Phi is " << patsyVec.Phi() << endl;
            cout << "PF pT is " << inEEPs.PFElecPt->at(i) << endl;
        }
        currElecPt = patsyVec.Pt();
        if (fabs(currElecPt - inEEPs.PFElecPt->at(i)) >= 10) continue;
        currElecECalE = inEEPs.ElecECalE->at(i); currElecSCEOverP = inEEPs.ElecSCEOverP->at(i); currElecSCEta = inEEPs.ElecSCEta->at(i);
        currRelPFElecIso = (inEEPs.ElecPFCharHadIso->at(i) + TMath::Max(0., inEEPs.ElecPFNeutHadIso->at(i) +  inEEPs.ElecPFPhotIso->at(i) - (eventRhoIso * ElectronEffectiveArea(currElecSCEta))))/currElecPt;
        currElecDZ = inEEPs.ElecDZ->at(i); currElecIP = inEEPs.ElecIP->at(i);
        currElecDeltaPhiIn = inEEPs.ElecDeltaPhiIn->at(i); currElecDeltaEtaIn = inEEPs.ElecDeltaEtaIn->at(i);
        currElecSigIetaIeta = inEEPs.ElecSigIetaIeta->at(i); currElecHtoERatio = inEEPs.ElecHtoERatio->at(i);
        currElecNumMissHits = inEEPs.ElecNumMissHits->at(i);        
        if (levelVerbosity > 0) {
            cout << "continuing for electron " << i << endl;
            cout << "currElecECalE " << currElecECalE << endl;
            cout << "currElecSCEOverP " << currElecSCEOverP << endl;
            cout << "currElecSCEta " << currElecSCEta << endl;
            cout << "currRelPFElecIso " << currRelPFElecIso << endl;
            cout << "currElecDZ " << currElecDZ << endl;
            cout << "currElecDeltaPhiIn " << currElecDeltaPhiIn << endl;
            cout << "currElecDeltaEtaIn " << currElecDeltaEtaIn << endl;
            cout << "currElecSigIetaIeta " << currElecSigIetaIeta << endl;
            cout << "currElecHtoERatio " << currElecHtoERatio << endl;
            cout << "currElecNumMissHits " << currElecNumMissHits << endl;
            cout << "ElecisEB->at(i) " << inEEPs.ElecisEB->at(i) << endl;
            cout << "ElecisEE->at(i) " << inEEPs.ElecisEE->at(i) << endl;
            cout << "isPFElectron->at(i) " << inEEPs.isPFElectron->at(i) << endl;
            cout << "passConvVeto->at(i) " << inEEPs.passConvVeto->at(i) << endl;
        }
        currElecPassCut = ElectronPassCutStage1(inEEPs.ElecisEB->at(i), inEEPs.ElecisEE->at(i), currElecSCEta, currElecDeltaPhiIn, currElecDeltaEtaIn, currElecSigIetaIeta, currElecHtoERatio);    
        if (levelVerbosity > 0) {
            cout << "continuing for electron " << i << endl;
            cout << "pass stage 1? " << currElecPassCut << endl;
        }
        currElecPassCut &= ElectronPassCutStage2(inEEPs.isPFElectron->at(i), inEEPs.passConvVeto->at(i), currRelPFElecIso, currElecIP, currElecDZ, currElecECalE, currElecSCEOverP, currElecNumMissHits);
        if (levelVerbosity > 0) {
            cout << "continuing for electron " << i << endl;
            cout << "pass stage 1 & 2? " << currElecPassCut << endl;
        }        
        if (currElecPassCut) {
            patsyLep.P4 = patsyVec;
            patsyLep.PDGID = inEEPs.ElecCharge->at(i) > 0 ? -11 : 11;
            patsyLep.relPFLepIso = currRelPFElecIso;
            vecIsoLeptons->push_back(patsyLep);
        }
        else {
            continue;
        }
    }
}
//inline void ElectronPickOviSyst(vector<float> * ElecPx, vector<float> * ElecPy, vector<float> * ElecPz, vector<float> * ElecEn, vector<int> * ElecCharge, vector<float> * PFElecPt, vector<float> * ElecSCEta, vector<float> * ElecDeltaPhiIn, vector<float> * ElecDeltaEtaIn, vector<float> * ElecSigIetaIeta, vector<float> * ElecHtoERatio, vector<float> * ElecIP, vector<float> * ElecDZ, vector<float> * ElecECalE, vector<float> * ElecSCEOverP, vector<int> * ElecNumMissHits, vector<float> * ElecPFNeutHadIso, vector<float> * ElecPFCharHadIso, vector<float> * ElecPFPhotIso, vector<bool> * passConvVeto, vector<bool> * isPFElectron, vector<bool> * ElecisEB, vector<bool> * ElecisEE, float eventRhoIso, float whichSystCase, vector<Lepton> * vecIsoLeptons, int levelVerbosity, vector<Lepton> * vecIsoLeptons_CentVal) {
inline void ElectronPickOviSyst(ElectronEventPointers inEEPs, float eventRhoIso, float whichSystCase, vector<Lepton> * vecIsoLeptons, int levelVerbosity, vector<Lepton> * vecIsoLeptons_CentVal) {
    bool currElecPassCut;
    float currElecPt;
    float currRelPFElecIso;
    float currElecECalE, currElecSCEOverP, currElecSCEta;
    float currElecIP, currElecDZ;
    float currElecDeltaPhiIn, currElecDeltaEtaIn, currElecSigIetaIeta, currElecHtoERatio;
    int currElecNumMissHits;
    Lepton patsyLep; 
    TLorentzVector patsyVec, patsyVec2;
    for (unsigned int i = 0; i < inEEPs.numElectrons; ++i) {
        patsyVec.SetPxPyPzE(inEEPs.ElecPx->at(i), inEEPs.ElecPy->at(i), inEEPs.ElecPz->at(i), inEEPs.ElecEn->at(i));
        currElecPt = patsyVec.Pt();
        if (fabs(currElecPt - inEEPs.PFElecPt->at(i)) >= 10) continue;
        patsyVec2.SetPxPyPzE(inEEPs.ElecPx->at(i), inEEPs.ElecPy->at(i), inEEPs.ElecPz->at(i), inEEPs.ElecEn->at(i));
        if (whichSystCase != 0) {
            patsyVec = LeptonScaleSystShift(patsyVec, 11, whichSystCase);
        }
        currElecECalE = inEEPs.ElecECalE->at(i); currElecSCEOverP = inEEPs.ElecSCEOverP->at(i); currElecSCEta = inEEPs.ElecSCEta->at(i);
        currRelPFElecIso = (inEEPs.ElecPFCharHadIso->at(i) + TMath::Max(0., inEEPs.ElecPFNeutHadIso->at(i) +  inEEPs.ElecPFPhotIso->at(i) - (eventRhoIso * ElectronEffectiveArea(currElecSCEta))))/currElecPt;
        currElecDZ = inEEPs.ElecDZ->at(i); currElecIP = inEEPs.ElecIP->at(i);
        currElecDeltaPhiIn = inEEPs.ElecDeltaPhiIn->at(i); currElecDeltaEtaIn = inEEPs.ElecDeltaEtaIn->at(i); 
        currElecSigIetaIeta = inEEPs.ElecSigIetaIeta->at(i); currElecHtoERatio = inEEPs.ElecHtoERatio->at(i);
        currElecNumMissHits = inEEPs.ElecNumMissHits->at(i);
        currElecPassCut = ElectronPassCutStage1(inEEPs.ElecisEB->at(i), inEEPs.ElecisEE->at(i), currElecSCEta, currElecDeltaPhiIn, currElecDeltaEtaIn, currElecSigIetaIeta, currElecHtoERatio);
        currElecPassCut &= ElectronPassCutStage2(inEEPs.isPFElectron->at(i), inEEPs.passConvVeto->at(i), currRelPFElecIso, currElecIP, currElecDZ, currElecECalE, currElecSCEOverP, currElecNumMissHits);
        
        if (currElecPassCut) {
            patsyLep.P4 = patsyVec;
            patsyLep.PDGID = inEEPs.ElecCharge->at(i) > 0 ? -11 : 11;
            patsyLep.relPFLepIso = currRelPFElecIso;
            vecIsoLeptons->push_back(patsyLep);
            patsyLep.P4 = patsyVec2;
            vecIsoLeptons_CentVal->push_back(patsyLep);
        }
        else {
            continue;
        }
    }
}
inline float goodVertexSelection(vector<float> * VertZ, vector<float> * VertRho, vector<float> * VertNDOF, vector<bool> * VertIsFake, int &nGoodVertices) {
    float firstGoodVertZ = -99;
    nGoodVertices = 0;
    for (unsigned int iVert = 0; iVert < VertZ->size(); ++iVert) {
        /*
        cout << "VertZ->at(iVert) " << VertZ->at(iVert) << endl;
        cout << "fabs(VertZ->at(iVert)) < 24 " << (fabs(VertZ->at(iVert)) < 24) << endl;
        cout << "VertRho   ->at(iVert) " << VertRho   ->at(iVert) << endl;
        cout << "VertRho   ->at(iVert)  <  2 " << (VertRho   ->at(iVert)  <  2) << endl;
        cout << "VertIsFake ->at(iVert) " << VertIsFake ->at(iVert) << endl;
        cout << "!VertIsFake ->at(iVert) " << (!VertIsFake ->at(iVert)) << endl;
        cout << " VertNDOF  ->at(iVert) " << VertNDOF  ->at(iVert) << endl;
        cout << " VertNDOF  ->at(iVert)  >  4 " << (VertNDOF  ->at(iVert)  >  4) << endl;
         */
        if (fabs(VertZ->at(iVert)) < 24 &&
            VertRho   ->at(iVert)  <  2 &&
            VertNDOF  ->at(iVert)  >  4 &&
            !VertIsFake ->at(iVert)) {
            //cout << " passed vertex cut! " << endl;
            ++nGoodVertices;
            if (nGoodVertices == 1) firstGoodVertZ = VertZ->at(iVert);
        }
    }
    //cout << "number of good vertices " << nGoodVertices << endl;
    return firstGoodVertZ;
}
inline bool MuonPassCutStage1(bool isGMPTMuon, bool isPFMuon, bool isGlobMuon, float relPFMuonIso, float muonEta) {
    bool muonPassCut;
    float muonIsoRatioCut = 0.15; float muonEtaCut = 2.4;
//    if (!isGMPTMuon || !isPFMuon || (!isGlobMuon && !isTrackArbitMuon) || !(relPFMuonIso < muonIsoRatioCut) || !(fabs(muonEta) < muonEtaCut)) {
    if (!isGMPTMuon || !isPFMuon || !isGlobMuon || !(relPFMuonIso < muonIsoRatioCut) || !(fabs(muonEta) < muonEtaCut)) {
        muonPassCut = false;
    }
    else {
        muonPassCut = true;
    }
    return muonPassCut;
}
inline bool MuonPassCutStage2(float muonD0, float muonDZ, int muonNumMatchStations, int muonNumValidPixHitsinTrack, int muonNumLayers) {
    float muonD0Cut = 0.2; float muonDZCut = 0.5;
    int muonNumMatchStationsCut = 1, muonNumValidPixHitsinTrackCut = 0, muonNumLayersCut = 5;
    bool muonPassCut;
    if (!(muonD0 < muonD0Cut) || !(muonDZ < muonDZCut) || !(muonNumMatchStations > muonNumMatchStationsCut) || !(muonNumLayers > muonNumLayersCut) || !(muonNumValidPixHitsinTrack > muonNumValidPixHitsinTrackCut)) {
        muonPassCut = false;
    }
    else {
        muonPassCut = true;
    }
    return muonPassCut;
}
inline bool MuonPassCutDESY(float relPFMuonIso, float muonEta) {
    bool muonPassCut;
    float muonIsoRatioCut = 0.15; float muonEtaCut = 2.4;
    if (!(relPFMuonIso < muonIsoRatioCut) || !(fabs(muonEta) < muonEtaCut)) {
        muonPassCut = false;
    }
    else {
        muonPassCut = true;
    }
    return muonPassCut;
}
//inline void MuonPickOvi(vector<float> * MuonPt, vector<float> * MuonPx, vector<float> * MuonPy, vector<float> * MuonPz, vector<float> * MuonEn, vector<int> * MuonCharge, vector<float> * PFMuonPt, vector<float> * MuonD0, vector<float> * MuonVertZ, float firstGoodVertZ, vector<float> * MuonPFNeutHadIso, vector<float> * MuonPFCharHadIso, vector<float> * MuonPFPhotIso, vector<float> * MuonSumPUPt, vector<int> * MuonNumLayers, vector<int> * MuonNumMatchStations, vector<int> * MuonNumValidPixHitsinTrack, vector<bool> * isGMPTMuon, vector<bool> * isPFMuon, vector<bool> * isGlobMuon, vector<Lepton> * vecIsoLeptons, int levelVerbosity) {
inline void MuonPickOvi(MuonEventPointers inMEPs, float firstGoodVertZ, vector<Lepton> * vecIsoLeptons, int levelVerbosity) {

    int currMuonNumMatchStations, currMuonNumValidPixHitsinTrack, currMuonNumLayers;
    float currMuonPt, currMuonEta;
    float currMuonRelPFIso;
    float currMuonDZ, currMuonD0;
    bool  currMuonPassCut;
    Lepton patsyLep;
    
    TLorentzVector patsyVec;
    for (unsigned int i = 0; i < inMEPs.numMuons; ++i) {
        //          https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleLepton2012#Muon_Selection
        //          https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
        patsyVec.SetPxPyPzE(inMEPs.MuonPx->at(i), inMEPs.MuonPy->at(i), inMEPs.MuonPz->at(i), inMEPs.MuonEn->at(i));
        currMuonPt = patsyVec.Pt();
        if (levelVerbosity > 0) {
            cout << "muon " << i << " has:" << endl;
            cout << "pT is " << currMuonPt << endl;
            cout << "Eta is " << patsyVec.Eta() << endl;
            cout << "Phi is " << patsyVec.Phi() << endl;
            cout << "PF pT is " << inMEPs.PFMuonPt->at(i) << endl;
        }
        //        if (fabs(currMuonPt - PFMuonPt->at(i)) >= 5.) continue;
        if (fabs(inMEPs.MuonPt->at(i) - inMEPs.PFMuonPt->at(i)) >= 5.) continue;
        currMuonEta = patsyVec.Eta();
        currMuonRelPFIso = (inMEPs.MuonPFCharHadIso->at(i) + TMath::Max(0., inMEPs.MuonPFPhotIso->at(i) + inMEPs.MuonPFNeutHadIso->at(i) - 0.5 * inMEPs.MuonSumPUPt->at(i)))/currMuonPt;
//        currMuonRelPFIso /= currMuonPt;
        currMuonDZ = fabs(firstGoodVertZ - inMEPs.MuonVertZ->at(i));        
        currMuonD0 = inMEPs.MuonD0->at(i);
        currMuonNumMatchStations = inMEPs.MuonNumMatchStations->at(i);
        currMuonNumLayers = inMEPs.MuonNumLayers->at(i);
        currMuonNumValidPixHitsinTrack = inMEPs.MuonNumValidPixHitsinTrack->at(i);
        
        if (levelVerbosity > 0) {
            cout << "continuing for muon " << i << endl;
            cout << "currMuonRelPFIso " << currMuonRelPFIso << endl;
            cout << "currMuonDZ " << currMuonDZ << endl;
            cout << "^^ coming from firstGoodVertZ = " << firstGoodVertZ << " and MuonVertZ->at(i) = " << inMEPs.MuonVertZ->at(i) << endl;
            cout << "currMuonD0 " << currMuonD0 << endl;
            cout << "currMuonNumMatchStations " << currMuonNumMatchStations << endl;
            cout << "currMuonNumLayers " << currMuonNumLayers << endl;
            cout << "currMuonNumValidPixHitsinTrack " << currMuonNumValidPixHitsinTrack << endl;
            cout << "isGMPTMuon->at(i) " << inMEPs.isGMPTMuon->at(i) << endl;
            cout << "isPFMuon->at(i) " << inMEPs.isPFMuon->at(i) << endl;
            cout << "isGlobMuon->at(i) " << inMEPs.isGlobMuon->at(i) << endl;
        }
        
        currMuonPassCut = MuonPassCutStage1(inMEPs.isGMPTMuon->at(i), inMEPs.isPFMuon->at(i), inMEPs.isGlobMuon->at(i), currMuonRelPFIso, currMuonEta);
        if (levelVerbosity > 0) {
            cout << "continuing for muon " << i << endl;
            cout << "pass stage 1? " << currMuonPassCut << endl;
        }
        currMuonPassCut &= MuonPassCutStage2(currMuonD0, currMuonDZ, currMuonNumMatchStations, currMuonNumValidPixHitsinTrack, currMuonNumLayers);
        if (levelVerbosity > 0) {
            cout << "continuing for muon " << i << endl;
            cout << "pass stage 1 & 2? " << currMuonPassCut << endl;
        }
        if (currMuonPassCut){
            patsyLep.P4 = patsyVec;
            patsyLep.PDGID = inMEPs.MuonCharge->at(i) > 0 ? -13 : 13;
            patsyLep.relPFLepIso = currMuonRelPFIso;
            vecIsoLeptons->push_back(patsyLep);
        }
        else {
            continue;
        }
    }
}

//inline void MuonPickOviSyst(vector<float> * MuonPt, vector<float> * MuonPx, vector<float> * MuonPy, vector<float> * MuonPz, vector<float> * MuonEn, vector<int> * MuonCharge, vector<float> * PFMuonPt, vector<float> * MuonD0, vector<float> * MuonVertZ, float firstGoodVertZ, vector<float> * MuonPFNeutHadIso, vector<float> * MuonPFCharHadIso, vector<float> * MuonPFPhotIso, vector<float> * MuonSumPUPt, vector<int> * MuonNumLayers, vector<int> * MuonNumMatchStations, vector<int> * MuonNumValidPixHitsinTrack, vector<bool> * isGMPTMuon, vector<bool> * isPFMuon, vector<bool> * isGlobMuon, float whichSystCase, vector<Lepton> * vecIsoLeptons, vector<Lepton> * vecIsoLeptons_CentVal) {    
inline void MuonPickOviSyst(MuonEventPointers inMEPs, float firstGoodVertZ, float whichSystCase, vector<Lepton> * vecIsoLeptons, int levelVerbosity, vector<Lepton> * vecIsoLeptons_CentVal) {    
    
    int currMuonNumMatchStations, currMuonNumValidPixHitsinTrack, currMuonNumLayers;
    float currMuonPt, currMuonEta;
    float currMuonRelPFIso;
    float currMuonDZ, currMuonD0;
    bool  currMuonPassCut;  
    Lepton patsyLep;

    TLorentzVector patsyVec, patsyVec2;        
    for (unsigned int i = 0; i < inMEPs.numMuons; ++i) {
        patsyVec.SetPxPyPzE(inMEPs.MuonPx->at(i), inMEPs.MuonPy->at(i), inMEPs.MuonPz->at(i), inMEPs.MuonEn->at(i));
        patsyVec2.SetPxPyPzE(inMEPs.MuonPx->at(i), inMEPs.MuonPy->at(i), inMEPs.MuonPz->at(i), inMEPs.MuonEn->at(i));
        //        if (fabs(patsyVec.Pt() - PFMuonPt->at(i)) >= 5.) continue;
        
        if (levelVerbosity > 0) {
            cout << "muon " << i << " has:" << endl;
            cout << "pT is " << patsyVec.Pt() << endl;
            cout << "Eta is " << patsyVec.Eta() << endl;
            cout << "Phi is " << patsyVec.Phi() << endl;
            cout << "PF pT is " << inMEPs.PFMuonPt->at(i) << endl;
        }
        
        if (fabs(inMEPs.MuonPt->at(i) - inMEPs.PFMuonPt->at(i)) >= 5.) continue; // want to check before shifting energy scale around
        if (whichSystCase != 0) {
            patsyVec = LeptonScaleSystShift(patsyVec, 13, whichSystCase);
        }    
        currMuonPt = patsyVec.Pt();
        currMuonEta = patsyVec.Eta();
        currMuonRelPFIso = (inMEPs.MuonPFCharHadIso->at(i) + TMath::Max(0., inMEPs.MuonPFPhotIso->at(i) + inMEPs.MuonPFNeutHadIso->at(i) - 0.5 * inMEPs.MuonSumPUPt->at(i)))/currMuonPt;
//        currMuonRelPFIso /= currMuonPt;
        currMuonDZ = fabs(firstGoodVertZ - inMEPs.MuonVertZ->at(i));
        currMuonD0 = inMEPs.MuonD0->at(i);
        currMuonNumMatchStations = inMEPs.MuonNumMatchStations->at(i);
        currMuonNumLayers = inMEPs.MuonNumLayers->at(i);
        currMuonNumValidPixHitsinTrack = inMEPs.MuonNumValidPixHitsinTrack->at(i);
        if (levelVerbosity > 0) {
            cout << "continuing for muon " << i << endl;
            cout << "currMuonRelPFIso " << currMuonRelPFIso << endl;
            cout << "currMuonDZ " << currMuonDZ << endl;
            cout << "^^ coming from firstGoodVertZ = " << firstGoodVertZ << " and MuonVertZ->at(i) = " << inMEPs.MuonVertZ->at(i) << endl;
            cout << "currMuonD0 " << currMuonD0 << endl;
            cout << "currMuonNumMatchStations " << currMuonNumMatchStations << endl;
            cout << "currMuonNumLayers " << currMuonNumLayers << endl;
            cout << "currMuonNumValidPixHitsinTrack " << currMuonNumValidPixHitsinTrack << endl;
            cout << "isGMPTMuon->at(i) " << inMEPs.isGMPTMuon->at(i) << endl;
            cout << "isPFMuon->at(i) " << inMEPs.isPFMuon->at(i) << endl;
            cout << "isGlobMuon->at(i) " << inMEPs.isGlobMuon->at(i) << endl;
        }        
        currMuonPassCut = MuonPassCutStage1(inMEPs.isGMPTMuon->at(i), inMEPs.isPFMuon->at(i), inMEPs.isGlobMuon->at(i), currMuonRelPFIso, currMuonEta);
        if (levelVerbosity > 0) {
            cout << "continuing for muon " << i << endl;
            cout << "pass stage 1? " << currMuonPassCut << endl;
        }
        currMuonPassCut &= MuonPassCutStage2(currMuonD0, currMuonDZ, currMuonNumMatchStations, currMuonNumValidPixHitsinTrack, currMuonNumLayers);
        if (levelVerbosity > 0) {
            cout << "continuing for muon " << i << endl;
            cout << "pass stage 1 & 2? " << currMuonPassCut << endl;
        }
        if (currMuonPassCut){
            patsyLep.P4 = patsyVec;
            patsyLep.PDGID = inMEPs.MuonCharge->at(i) > 0 ? -13 : 13;
            patsyLep.relPFLepIso = currMuonRelPFIso;
            vecIsoLeptons->push_back(patsyLep);
            //            vecIsoLeptons->push_back(patsyVec);
            patsyLep.P4 = patsyVec2;
            vecIsoLeptons_CentVal->push_back(patsyLep);
        }
        else {
            continue;
        }
    }
}

inline void IsoLeptonsPickDESY(VLV * Leptons, vector<int> *lepPdgId, vector<double> *lepPFIso, float whichSystCase, vector<Lepton> * vecIsoLeptons, vector<Lepton> * vecIsoLeptons_CentVal) {    
    unsigned int vecSize = Leptons->size();
    if (vecSize < 2) {
        return;
    }
    TLorentzVector patsyVec, patsyVec2;
    float currLepEta, currLepRelPFIso;
    int   currLepPDGID;
    Lepton patsyLep;
    bool  currLepPassCut = false;    
    for (unsigned int iLep = 0; iLep < vecSize; ++iLep) {
        patsyVec.SetPxPyPzE(Leptons->at(iLep).Px(), Leptons->at(iLep).Py(), Leptons->at(iLep).Pz(), Leptons->at(iLep).E());
        patsyVec2.SetPxPyPzE(Leptons->at(iLep).Px(), Leptons->at(iLep).Py(), Leptons->at(iLep).Pz(), Leptons->at(iLep).E());
        currLepEta = patsyVec.Eta();
        currLepRelPFIso = lepPFIso->at(iLep);
        currLepPDGID = lepPdgId->at(iLep);        
        if (whichSystCase != 0) {
            patsyVec = LeptonScaleSystShift(patsyVec, abs(currLepPDGID), whichSystCase);
        }
        if (abs(currLepPDGID) == 13) {
            currLepPassCut = MuonPassCutDESY(currLepRelPFIso, currLepEta);
        }
        else if (abs(currLepPDGID) == 11) {
            currLepPassCut = ElectronPassCutDESY(currLepRelPFIso, currLepEta);
        }
        if (currLepPassCut) {            
            patsyLep.P4 = patsyVec;
            patsyLep.PDGID = currLepPDGID;
            patsyLep.relPFLepIso = currLepRelPFIso;
            vecIsoLeptons->push_back(patsyLep);
            patsyLep.P4 = patsyVec2;
            vecIsoLeptons_CentVal->push_back(patsyLep);
        }
        else {
            continue;
        }
    }
}

inline void IsoLeptonsPickDESY(VLV * Leptons, vector<int> *lepPdgId, vector<double> *lepPFIso, vector<Lepton> * vecIsoLeptons) {    
    unsigned int vecSize = Leptons->size();
    if (vecSize < 2) {
        return;
    }
    TLorentzVector patsyVec;
    float currLepEta, currLepRelPFIso;
    int   currLepPDGID;
    bool  currLepPassCut = false;  
    Lepton patsyLep;
    for (unsigned int iLep = 0; iLep < vecSize; ++iLep) {
        patsyVec.SetPxPyPzE(Leptons->at(iLep).Px(), Leptons->at(iLep).Py(), Leptons->at(iLep).Pz(), Leptons->at(iLep).E());
        currLepEta = patsyVec.Eta();
        currLepRelPFIso = lepPFIso->at(iLep);
        currLepPDGID = lepPdgId->at(iLep);
        if (abs(currLepPDGID) == 13) {
            currLepPassCut = MuonPassCutDESY(currLepRelPFIso, currLepEta);
        }
        else if (abs(currLepPDGID) == 11) {
            currLepPassCut = ElectronPassCutDESY(currLepRelPFIso, currLepEta);
        }
        if (currLepPassCut) {
            patsyLep.P4 = patsyVec;
            patsyLep.PDGID = currLepPDGID;
            patsyLep.relPFLepIso = currLepRelPFIso;
            vecIsoLeptons->push_back(patsyLep);
        }
        else {
            continue;
        }
    }
}
inline EventLepInfo LeptonPair(vector<Lepton> * vecIsoLeptons, int levelVerbosity) {
    EventLepInfo outELI; outELI.ELIDefaultVarVals();
    int lep0Index = -1; int lep1Index = -1;
    int vecSize = vecIsoLeptons->size();
    if (levelVerbosity > 0) {
        cout << endl;
        cout << "starting on trying to pair up leptons " << endl;
        cout << "vecSize " << vecSize << endl;   
    }
    if (vecSize < 2) {
        outELI.EventFails();
        return outELI;
    }
    float leadLepPtCut = 20., subLepPtCut = 10.;
    float currLeadLepPt, currSubLepPt;
    float currLeadLepEta, currSubLepEta;
    float currLeadLepCharge, currSubLepCharge;
    int   currLeadLepPDGID, currSubLepPDGID;
    int   productPdgId;
    float currDiLepPt;
    float leadDiLepPt = 0;
    float diLepmassCut = 20.;
    for (int iLep = 0; iLep < vecSize; ++iLep) {
        currLeadLepPt = vecIsoLeptons->at(iLep).P4.Pt();
        currLeadLepEta = vecIsoLeptons->at(iLep).P4.Eta();
        currLeadLepPDGID = vecIsoLeptons->at(iLep).PDGID;
        currLeadLepCharge = TMath::Sign(1, currLeadLepPDGID);
        if (currLeadLepPDGID == -13) {
            if (currLeadLepPt > leadLepPtCut) {
                outELI.EventNIsoMubars_pT20++;
            }
            else if (currLeadLepPt > subLepPtCut) {
                outELI.EventNIsoMubars_pT10to20++;
            }
        }
        else if (currLeadLepPDGID == 13) {
            if (currLeadLepPt > leadLepPtCut) {
                outELI.EventNIsoMuons_pT20++;
            }
            else if (currLeadLepPt > subLepPtCut) {
                outELI.EventNIsoMuons_pT10to20++;
            }
        }
        else if (currLeadLepPDGID == 11) {
            if (currLeadLepPt > leadLepPtCut) {
                outELI.EventNIsoElecs_pT20++;
            }
            else if (currLeadLepPt > subLepPtCut) {
                outELI.EventNIsoElecs_pT10to20++;
            }
        }
        else if (currLeadLepPDGID == -11) {
            if (currLeadLepPt > leadLepPtCut) {
                outELI.EventNIsoPosits_pT20++;
            }
            else if (currLeadLepPt > subLepPtCut) {
                outELI.EventNIsoPosits_pT10to20++;
            }
        }
        if (levelVerbosity > 0) {
            cout << "For attempted leading lepton, vector position " << iLep << endl;
            cout << "currLeadLepPt = " << currLeadLepPt << endl;
            cout << "currLeadLepEta = " << currLeadLepEta << endl;
            cout << "currLeadLepPhi = " << vecIsoLeptons->at(iLep).P4.Phi() << endl;
            cout << "currLeadLepPDGID = " << currLeadLepPDGID << endl;
            cout << "currLeadLepCharge = " << currLeadLepCharge << endl;
        }
        if (currLeadLepPt < leadLepPtCut) continue;
        for (int iLep2 = iLep+1; iLep2 < vecSize; ++iLep2) {
            currSubLepPt = vecIsoLeptons->at(iLep2).P4.Pt();
            currSubLepEta = vecIsoLeptons->at(iLep2).P4.Eta();
            currSubLepPDGID = vecIsoLeptons->at(iLep2).PDGID;
            currSubLepCharge = TMath::Sign(1, currSubLepPDGID);
            if (levelVerbosity > 0) {
                cout << "For attempted leading lepton, vector position " << iLep2 << endl;
                cout << "currSubLepPt = " << currSubLepPt << endl;
                cout << "currSubLepEta = " << currSubLepEta << endl;
                cout << "currSubLepPhi = " << vecIsoLeptons->at(iLep2).P4.Phi() << endl;
                cout << "currSubLepPDGID = " << currSubLepPDGID << endl;
                cout << "currSubLepCharge = " << currSubLepCharge << endl;
            }
            if (currSubLepPt < subLepPtCut) continue;
            outELI.EventNViableLepPairsPreMassCut++;
            if (levelVerbosity > 0) {
                productPdgId = currLeadLepPDGID * currSubLepPDGID;
                float currPairMass = (vecIsoLeptons->at(iLep).P4 + vecIsoLeptons->at(iLep2).P4).M();
                currDiLepPt = currLeadLepPt + currSubLepPt;
                cout << "for iLep = " << iLep << " and iLep2 = " << iLep2 << endl;
                cout << "productPdgId = " << productPdgId << endl;
                cout << "currPairMass = " << currPairMass << endl;
                cout << "currDiLepPt = " << currDiLepPt << endl;
            }
            if ((vecIsoLeptons->at(iLep).P4 + vecIsoLeptons->at(iLep2).P4).M() < diLepmassCut) continue;
            if (currLeadLepCharge * currSubLepCharge > 0) continue;//i.e. same sign event
            currDiLepPt = currLeadLepPt + currSubLepPt;            
            if (currDiLepPt > leadDiLepPt) {
                leadDiLepPt = currDiLepPt;
                lep0Index = iLep;
                lep1Index = iLep2;
            }
        }
    }
    if (levelVerbosity) {
        cout << "lep0index in isolated leptons vector " << lep0Index << endl;
        cout << "lep1index in isolated leptons vector " << lep1Index << endl;
    }
    if (lep0Index < 0 || lep1Index < 0) {
        outELI.EventFails();
    }
    else {
        outELI.EventPasses(lep0Index, lep1Index, vecIsoLeptons);
    }
    if (levelVerbosity) {
        cout << "Lep0PdgId " << outELI.EventLep0PDGID << endl;
        cout << "Lep1PdgId " << outELI.EventLep1PDGID << endl;
        cout << "Dilepton Event Type " << outELI.EventDiLepType << endl;
    }
    return outELI;
}

inline float JESUncertFactor(TH2F * histJES, float shiftDirection, float JetPt, float JetEta) {
    int NXbins = histJES->GetNbinsX();
    TAxis * XAxis = histJES->GetXaxis();
    TAxis * YAxis = histJES->GetYaxis();
    int binX = XAxis->FindBin(JetPt);
    int binY = YAxis->FindBin(JetEta);
    float uncertFactor;
    float linYPointOne, linYPointTwo, linXPointOne, linXPointTwo, rise, run, slope;
    if (binX >= NXbins) uncertFactor = histJES->GetBinContent(NXbins, binY);
    else if (binX < 1) uncertFactor = histJES->GetBinContent(1, binY);
    else {
        linYPointOne = histJES->GetBinContent(binX, binY);
        linYPointTwo = histJES->GetBinContent(binX + 1, binY);
        linXPointOne = XAxis->GetBinLowEdge(binX);
        linXPointTwo = XAxis->GetBinUpEdge(binX);
        rise = linYPointTwo - linYPointOne;
        run = linXPointTwo - linXPointOne;
        slope = rise / run;
        uncertFactor = linYPointOne + slope * (JetPt - linXPointOne);
    }
    return (1 + shiftDirection * uncertFactor);
}
inline bool isLeptonJet(vector<Lepton> * isoLeptons, TLorentzVector * inputJet) {
    float dRcut = 0.4;
    float lepPtCut;// = 10.0;
    float jetPtCutForCutLep = 20.0;
    for (unsigned int iLep = 0; iLep < isoLeptons->size(); ++iLep) {
        lepPtCut = (inputJet->Pt() > jetPtCutForCutLep) ? 10.0: 0.0; // for jets > 20 GeV, check that isoLepton's pT is > 10 GeV, otherwise don't care about pT
        if (isoLeptons->at(iLep).P4.Pt() > lepPtCut) {
            if (deltaR(inputJet->Eta(), inputJet->Phi(), isoLeptons->at(iLep).P4.Eta(), isoLeptons->at(iLep).P4.Phi()) < dRcut) return true;
        }
    }
    return false;
}
/*
bool jetPassPFID(float jetNHF, float jetNEF, float jetEta, int jetNConst, float jetCHF, float jetCEF, int jetChMult, int whichCut, bool beVerbose) {
    //taken from this twiki https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    // whichCut: 0 = Loose, 1 = Medium, 2 = Tight    
    float etaSwitch = 2.4;
    float cutNHF[3] = {0.99, 0.95, 0.90};
    float cutNEF[3] = {0.99, 0.95, 0.90};
    float cutCHF[3] = {0., 0., 0.};
    float cutCEF[3] = {0.99, 0.99, 0.99};
    int cutNConst[3] = {1, 1, 1};
    int cutChMult[3] = {0, 0, 0};
    
    bool failNHF = (jetNHF >= cutNHF[whichCut]);
    bool failNEF = (jetNEF >= cutNEF[whichCut]);
    bool failNConst = (jetNConst <= cutNConst[whichCut]);
    bool failCHF = false;
    bool failCEF = false;
    bool failChMult = false;
    if (abs(jetEta) < etaSwitch) {
        failCHF = (jetCHF <= cutCHF[whichCut]);
        failCEF = (jetCEF >= cutCEF[whichCut]);
        failChMult = (jetChMult <= cutChMult[whichCut]);
    }
    bool failAny = (failNHF || failNEF || failNConst || failChMult || failCHF || failCEF);
    
    if (beVerbose) {
        cout << "jetNHF " << jetNHF << " is " << failNHF << endl;
        cout << "jetNEF " << jetNEF << " is " << failNEF <<endl;
        cout << "jetNConst " << jetNConst << " is " << failNConst <<endl;
        cout << "jetEta " << jetEta << endl;
        cout << "jetCHF " << jetCHF << " is " << failCHF <<endl;
        cout << "jetCEF " << jetCEF << " is " << failCEF <<endl;
        cout << "jetChMult " << jetChMult << " is " << failChMult <<endl;
        cout << "failAny? " << failAny << endl;
    }
    return !failAny;
}
*/
/*
inline bool passLoosePFJetID(float currJetEta, int currJetIndex, PFJetEventPointers inPFJEPs, int strengthCut) {
    bool passID = (inPFJEPs.JetNHF->at(iJet) < 0.99 && inPFJEPs.JetNEF->at(iJet) < 0.99);
    if (fabs(currJetEta < 2.4)) {
        passID &= (inPFJEPs.JetCEF->at(iJet) < 0.99 && inPFJEPs.JetCHF->at(iJet) > 0 && inPFJEPs.JetCharMult->at(iJet) > 0);
    }
    return passID;
}
*/
inline vector<PFJet> * JetInfoDESY(vector<Lepton> * isoLeptons, VLV * Jets, vector<double> * JetBTag, float whichSystCase, TH2F * shiftHist) {
    vector<PFJet> * outVecJets = new vector<PFJet>;
    TLorentzVector patsyJetLV;
    PFJet currPFJet;
    bool leptonJet;
    float dRcut = 0.4;
    float JetFactor;
    Double_t currJetPt, currJetEta, currJetPhi;
    for (unsigned int iJet = 0; iJet < Jets->size(); ++iJet) {
        leptonJet = 0;        
        patsyJetLV.SetPxPyPzE(Jets->at(iJet).Px(), Jets->at(iJet).Py(), Jets->at(iJet).Pz(), Jets->at(iJet).E());
        currJetPt = patsyJetLV.Pt();
        currJetEta = patsyJetLV.Eta();
        currJetPhi = patsyJetLV.Phi();
        if (isLeptonJet(isoLeptons, &patsyJetLV)) continue;
        if (currJetPt < 9.0) continue;
        if (whichSystCase != 0) {
            JetFactor = JESUncertFactor(shiftHist, whichSystCase, (float) currJetPt, (float) currJetEta);
            patsyJetLV *= JetFactor;
        }
        currPFJet.P4 = patsyJetLV;
        currPFJet.valBTagDisc = (float) JetBTag->at(iJet);
        currPFJet.partonFlavor = -999999;
        outVecJets->push_back(currPFJet);
    }
    return outVecJets;
}
//inline vector<PFJet> * JetInfo(vector<Lepton> * Leptons, vector<float> * JetPx, vector<float> * JetPy, vector<float> * JetPz, vector<float> * JetE, vector<float> * JetNHF, vector<float> * JetNEF, vector<float> * JetCHF, vector<float> * JetCEF, vector<int> * JetNDaug, vector<int> * JetCharMult, vector<float> * JetBTag, vector<int> * JetPartFlav, float whichSystCase, TH2F * shiftHist) {
inline vector<PFJet> * JetInfo(vector<Lepton> * isoLeptons, PFJetEventPointers inPFJEPs, JetCutInfo * inJCI, float whichSystCase, TH2F * shiftHist) {
    vector<PFJet> * Jets = new vector<PFJet>;
    TLorentzVector patsyJetLV;
    PFJet currPFJet;
    bool leptonJet;
    float dRcut = 0.4;
    float JetFactor;
    Double_t currJetPt, currJetEta, currJetPhi;
    for (unsigned int iJet = 0; iJet < inPFJEPs.numPFJets; ++iJet) {
        leptonJet = 0;
        currPFJet.SetJetVars(inPFJEPs, iJet);
        currJetPt = currPFJet.P4.Pt();
        currJetEta = currPFJet.P4.Eta();
        currJetPhi = currPFJet.P4.Phi();
        if (isLeptonJet(isoLeptons, &currPFJet.P4)) continue;
        if (currJetPt < 9.0) continue;
/*
        currPFJet.PrintVals();
        cout <<endl;
        */
        if (currPFJet.PassesID(inJCI)) {
            /*
            cout <<endl;
            cout << "passedID " << endl;
            cout <<endl;
            cout <<endl;
            */
            if (whichSystCase != 0) {
                JetFactor = JESUncertFactor(shiftHist, whichSystCase, (float) currJetPt, (float) currJetEta);
                currPFJet.P4 *= JetFactor;
            }
            Jets->push_back(currPFJet);
        }
        else {
            continue;
        }
    }
    return Jets;
}
inline float JetResolutionHighPt(TLorentzVector * inputJet, vector<TF1> * inputTF1Vec) {
    float currJetAbsEta = fabs(inputJet->Eta());
    float EtaBinEnds[10] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 9.9}; // in reality 9 TF1 vecs
    int   currJetEtaIndex = -1;
    bool  cornerCaseBinEdge = false;
    float outResolution;
    for (unsigned int iEta = 0; iEta < 9; ++iEta) {
        if (currJetAbsEta > EtaBinEnds[iEta] && currJetAbsEta < EtaBinEnds[iEta + 1]) {
            currJetEtaIndex = iEta;
            break;
        }
        else if (currJetAbsEta == EtaBinEnds[iEta] && currJetAbsEta > 0) {
            currJetEtaIndex = iEta;
            cornerCaseBinEdge = true;
            cout << "corner case came up, jet eta is " << inputJet->Eta() << endl;
            // you've somehow managed to fall right on the bin edge -- floating point precision -- so take average of two different bins -- currently rolling with straight up mean
        }
    }
    if (currJetEtaIndex < 0) {          
        cout << "can't do this jet!!...technically jet eta is: " << inputJet->Eta() << endl; 
        currJetEtaIndex = 8;        
    }
    cout << "currJetEtaIndex " << currJetEtaIndex << endl;;
    cout << "currJetPt " << inputJet->Pt() << endl;
    outResolution = inputTF1Vec->at(currJetEtaIndex).Eval(inputJet->Pt());
    outResolution *= inputJet->Pt(); // assuming the outResolution is a fractional resolution
    /*
    if (cornerCaseBinEdge) {
        outResolution += inputTF1Vec->at(currJetEtaIndex - 1).Eval(inputJet->Pt());
        outResolution /= 2;
    }
    */
    return outResolution;
}
inline bool GenJetMatch(TLorentzVector * inputJet, TLorentzVector * genJet, float &dRGenReco) {
    bool isGenJetMatched;
    float dRmaxGenJetMatch = TMath::Min(0.5, 0.1 + 0.3 * TMath::Exp(-0.05 * (genJet->Pt() - 10)));
    dRGenReco = deltaR(inputJet->Eta(), inputJet->Phi(), genJet->Eta(), genJet->Phi());
    isGenJetMatched = (dRGenReco < dRmaxGenJetMatch);
    return isGenJetMatched;
}
inline vector<GenJet> * GenJetsNonZero(GenJetEventPointers inGJEPs) {
    GenJet currGenJet;
    vector<GenJet> * outGenJetVec = new vector<GenJet>;
    for (unsigned int iGenJet = 0; iGenJet < inGJEPs.numGenJets; ++iGenJet) {
        if (inGJEPs.genJetIsGenJet->at(iGenJet)) {
            currGenJet.P4.SetPxPyPzE(inGJEPs.genJetPx->at(iGenJet), inGJEPs.genJetPy->at(iGenJet), inGJEPs.genJetPz->at(iGenJet), inGJEPs.genJetEn->at(iGenJet));
            currGenJet.invisE = inGJEPs.genJetInvisE->at(iGenJet);
            currGenJet.seedPDGID = 0;
            currGenJet.isBJet = false;
            outGenJetVec->push_back(currGenJet);
        }
    }
    return outGenJetVec;
}
inline TLorentzVector SmearJet(TLorentzVector * inputJet, vector<GenJet> * vecGoodGenJets, bool &isGenMatched, float &dEnRecoGen, TH2F * GenJetSmearHist, TH2F * RecoJetLowPtSmearHist, vector<TF1> * inputTF1Vec, float smearBy = 1.0, float shiftSmearSyst = 0.0, int levelVerbosity = 0) {
    /// Function to smear the jets according to C.V.'s prescription. Look at the produce method for the SmearedJet EDProducer shown here:
    /// http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#207
    /// As of 10/10/13 I am missing the Reco vs. Corr Jet Max function, but otherwise things are in order with what he does
    /// This function basically comes about from tracing down associated configurations and EDProducers from running the runMEtUncertainties.py tool defined here:
    /// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatUtils/python/tools/runType1PFMEtUncertainties.py?view=log
    TRandom3 rnd_;
    float x = fabs(inputJet->Eta()), y = inputJet->Pt(); // set up the "x" and "y" to look for in the resolution histograms
    float inEn = inputJet->E();
    TLorentzVector smearJet = *inputJet;    // initializing incase I don't smear jet for some reasons
    if (levelVerbosity > 0) {
        cout << endl;
        cout << "TRYING TO SMEAR JET " << endl;
        cout << "smear jet initial px " << smearJet.Px() << endl;
        cout << "smear jet initial py " << smearJet.Py() << endl;
        cout << "smear jet initial pz " << smearJet.Pz() << endl;
        cout << "smear jet initial E " << smearJet.E() << endl;
        cout << "smear jet initial Pt " << smearJet.Pt() << endl;
    }
    float smearJetEnergy = inEn;            // initializing incase I don't smear jet for some reasons
    
    TLorentzVector currGenJetLV;
    float dRGenReco;
    float inputJetBaseRes, inputJetScaleRes;
    float smearFactor = 1.0, smearFactorErr = 0.0;
    float pTResolThreshold = 10.0, sigmaMaxGenJetMatch = 3.0; 
    float skipRawJetPtThreshold = 10.0, skipCorrJetPtThreshold = 0.01, minSmearJetEn = 0.01;
    int bestGenMatchIndex = -1;
    float bestdRGenReco = 1.0;
    isGenMatched = false;
    
    if (x > GenJetSmearHist->GetXaxis()->GetXmin() && x < GenJetSmearHist->GetXaxis()->GetXmax() && 
        y > GenJetSmearHist->GetYaxis()->GetXmin() && y < GenJetSmearHist->GetYaxis()->GetXmax() ) {
        int binIndexGenJetMatch = GenJetSmearHist->FindBin(x, y);
        smearFactor += smearBy * (GenJetSmearHist->GetBinContent(binIndexGenJetMatch) - 1);
        smearFactorErr = GenJetSmearHist->GetBinError(binIndexGenJetMatch);   
    }
    int binIndexNoGenJetMatch = RecoJetLowPtSmearHist->FindBin(x, y);
    if (levelVerbosity > 0) {
        cout << " smear factor pre shift " << smearFactor << endl;
        cout << " smearFactorErr " << smearFactorErr << endl; 
    }

    if (shiftSmearSyst != 0) smearFactor += shiftSmearSyst * smearFactorErr;
    if (smearFactor < 1.) smearFactor = 1.;
    if (levelVerbosity > 0) {
        cout << " smear factor post shift " << smearFactor << endl;        
        cout << "multiply factor " << TMath::Sqrt(smearFactor * smearFactor - 1.) << endl;
    }
    if (y < pTResolThreshold && y < 20.) {
        inputJetBaseRes = RecoJetLowPtSmearHist->GetBinContent(binIndexNoGenJetMatch);
    }
    else {
        inputJetBaseRes = JetResolutionHighPt(inputJet, inputTF1Vec);
    }
    if (levelVerbosity > 0) {
        cout << "inputJetBaseRes pre multiply " << inputJetBaseRes << endl;   
    }
    inputJetScaleRes = inputJetBaseRes * TMath::Sqrt(smearFactor * smearFactor - 1.);
    if (levelVerbosity > 0) {
        cout << "inputJetBaseRes post multiply (i.e. inputJetScaleRes) " << inputJetScaleRes << endl;   
    }
    for (unsigned int iGenJet = 0; iGenJet < vecGoodGenJets->size(); ++iGenJet) {
//        cout << " test " << endl;        
        currGenJetLV = vecGoodGenJets->at(iGenJet).P4;
//        cout << " test 2 " << currGenJetLV.Pt() <<        endl;
        if (GenJetMatch(inputJet, &currGenJetLV, dRGenReco)) {
            if (bestdRGenReco > dRGenReco) {
                if (!isGenMatched) isGenMatched = true;
                bestGenMatchIndex = iGenJet;
                bestdRGenReco = dRGenReco;
            }
        }
        else {
            continue;
        }
    }
    if (levelVerbosity > 0) {
        cout << "isGenMatched? " << isGenMatched << endl;   
    }
    if (isGenMatched) {
        dEnRecoGen = inEn - vecGoodGenJets->at(bestGenMatchIndex).P4.E();
        if (levelVerbosity > 0) {
            cout << " dEnRecoGen " << dEnRecoGen << endl;
            cout << " sigmaMaxGenJetMatch * inputJetBaseRes " << sigmaMaxGenJetMatch * inputJetBaseRes << endl;
            cout << "(smearFactor - 1.) * dEnRecoGen / inEn " << (smearFactor - 1.) * dEnRecoGen / inEn << endl;
        }
        if (fabs(dEnRecoGen) < (sigmaMaxGenJetMatch * inputJetBaseRes)) {
            smearJetEnergy = inEn * (1. + (smearFactor - 1.) * dEnRecoGen / inEn); 
            //^^ technically this isn't matching the analogous line from the smeared jet producer (line 271) but I don't have access to raw jet info as of right now
        }
        else {
            isGenMatched = false;
        }
    }
    if (!isGenMatched) {
        // gaussian smearing when no gen jet matches
        if (smearFactor > 1. ) {
            smearJetEnergy = inEn * (1. + rnd_.Gaus(0., inputJetScaleRes)/inEn);
            cout << "relative difference: smearJetEnergy - normal E " << (smearJetEnergy - inEn)/inEn << endl;
            //^^ technically this isn't matching the analogous line from the smeared jet producer (line 291) but I don't have access to raw jet info as of right now
        }
    }
    if (smearJetEnergy < minSmearJetEn) smearJetEnergy = minSmearJetEn;
    if (!(y < skipCorrJetPtThreshold )) {
        smearJet *= (smearJetEnergy/inEn);
    }
    if (levelVerbosity > 0) {
        if (isGenMatched ) {
            cout << " jet was Gen Jet matched! " << endl;
        }
        else {
            cout << " jet was NOT Gen Jet matched! " << endl;
        }
        cout << "smearJetEnergy (i.e. post smearing energy) " << smearJetEnergy << endl;
        cout << "smear jet post px " << smearJet.Px() << endl;
        cout << "smear jet post py " << smearJet.Py() << endl;
        cout << "smear jet post pz " << smearJet.Pz() << endl;
        cout << "smear jet post E " << smearJet.E() << endl;
        cout << "smear jet post Pt " << smearJet.Pt() << endl;
        cout << endl;
    }
    return smearJet;
    // 
    // Line 254 of the EDProducer shown above has the following
    // -- double sigmaEn = jetResolutionExtractor_(jet)*TMath::Sqrt(smearFactor*smearFactor - 1.);
    // Looking at line 341 of this file
    // -- SmearedJetProducer_namespace::JetResolutionExtractorT<T> jetResolutionExtractor_;
    // Look here for the definition of what JetResolutionExtractorT is:
    // http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatUtils/plugins/SmearedJetProducer.cc#015
    // You'll see this salient line (line 26):
    // -- metsig::SigInputObj pfJetResolution = jetResolutions_.evalPFJet(&jet);
    // Looking at a combination of 
    // http://cmslxr.fnal.gov/lxr/source/DataFormats/METReco/interface/SigInputObj.h#047
    // and 
    // http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/SignAlgoResolutions.cc#109
    // you can recognize that line 134 of the second of the aforementioned two files
    // -- metsig::SigInputObj obj_jet(inputtype,jpt,jphi,jdeltapt,jdeltapphi);
    // shows that jdeltapt is the sigmaEn used in line 254 of the SmearedJetEDProducer
    // This jdeltapt is different for low and high pT jets
    // low pT jets are based off of the jdpt vectors (numbering for different eta bins)
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoMET/METProducers/python/METSigParams_cfi.py?revision=1.3&view=markup
    // Binning for |eta| is steps of 0.5 from 0.0 to 5.0
    // while binning for pT is steps of 2 GeV from 1 to 21 GeV
    // high pT jets are based off of a double sided crystal ball function, specifically the sigma of the gaussian core
    // as can basically be surmised from a combination of 
    // http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/SignAlgoResolutions.cc#129
    // http://cmslxr.fnal.gov/lxr/source/CondFormats/JetMETObjects/src/JetResolution.cc#082
    // and 
    // http://cmslxr.fnal.gov/lxr/source/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt
    // pT/Eta dependent functions were fit for each of the parameters of the DSCB function, and you can see more discussion for the 
    // "sigma" parameter in the function NoGenJetMatchHighPtResolutionTF1 defined in this header file
}
inline TH2F * ResolutionHistMaker(TString inputFile) {
    // Following Discussion in SmearJet Function this function tries to load in the resolutions defined here:
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoMET/METProducers/python/METSigParams_cfi.py?revision=1.3&view=markup
    // (look for jdpt# (where # runs from 0 to 9)
    // the jdpt# refers to different Eta bins, as can be ascertained from 
    // http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/SignAlgoResolutions.cc#121
    // and 
    // http://cmslxr.fnal.gov/lxr/source/RecoMET/METAlgorithms/src/SignAlgoResolutions.cc#151
    // The function reads in a text file where these parameters have been put in, 
    // with the first two lines contain pT eta information (number of eta/Pt bins and bin ranges)
    
    TH2F * outHist;
    ifstream inputFileStream(inputFile);
    vector<unsigned int> pTParams, EtaParams;
    //    int nPtBins, PtMin, PtMax, nEtaBins, EtaMin, EtaMax;
    string line, field;
    int pTCounter;
    if (!(inputFileStream.eof())) {
        getline(inputFileStream, line);
        stringstream ss( line );
        while (std::getline(ss, field, ',' )) {
            stringstream fs ( field );
            int currField = 0;
            fs >> currField;
            pTParams.push_back(currField);
        }
    }
    else {
        cout << "issue with inputFile: prematurely reached end of file " << endl;
    }
    if (!(inputFileStream.eof())) {
        getline(inputFileStream, line);
        stringstream ss( (line) );
        while (std::getline(ss, field, ',' )) {
            stringstream fs ( field );
            int currField = 0;
            fs >> currField;
            EtaParams.push_back(currField);
        }
    }
    else {
        cout << "issue with inputFile: prematurely reached end of file " << endl;
    }
    outHist = new TH2F("h_MCJetResol", ";Jet #eta;Jet p_{T}", EtaParams[0], EtaParams[1], EtaParams[2], pTParams[0], pTParams[1], pTParams[2]);
    for (unsigned int iEta = 0; iEta < EtaParams[0]; ++iEta) {
        pTCounter = 0;
        if (!(inputFileStream.eof())) {
            getline(inputFileStream, line);
            stringstream ss( (line) );
            while (std::getline(ss, field, ',' )) {
                stringstream fs ( field );
                float currField = 0.0;
                fs >> currField;
                if (pTCounter == 0) outHist->SetBinContent(iEta + 1, pTCounter, currField);
                    outHist->SetBinContent(iEta + 1, pTCounter + 1, currField);
                ++pTCounter;
            }
        }
        else {
            cout << "issue with inputFile: prematurely reached end of file " << endl;
        }
    }
    for (int i = 1; i < outHist->GetNbinsX() + 1; ++i) {
        for (int j = 1; j < outHist->GetNbinsY() + 1; ++j) {
            cout << "bin (" << i << ", " << j << ") is (eta = " << outHist->GetXaxis()->GetBinLowEdge(i) << ", pT = " << outHist->GetYaxis()->GetBinLowEdge(j) << endl;
            cout << "bin content in this bin: " << outHist->GetBinContent(i, j) << endl; 
        }        
    }
    return outHist;
}


inline vector<PFJet> * JetSmear(vector<PFJet> * InputJets, vector<GenJet> * vecGoodGenJets, float smearBy, int whichSystToDo, float whichSystCase, TH2F * shiftHistJetES, TH2F * shiftHistLowPt, TH2F * shiftHistSmearFactor, vector<TF1> * inputTF1Vec, int levelVerbosity, TH1F * DeltaEnergyHist) {
    vector<PFJet> * outSmearJets = new vector<PFJet>;
    TLorentzVector patsyJetLV, patsyJetLVSmear;
    PFJet currPFJet;
    float JetFactor;
    float currJetPt, currJetEta;
    float smearSystShift = (whichSystToDo == 2) ? whichSystCase : 0.0;
    for (unsigned int iJet = 0; iJet < InputJets->size(); ++iJet) {
        currPFJet = InputJets->at(iJet);
        if (levelVerbosity > 0) {
            cout << "currently on Jet " << iJet << endl;
            cout << " which has: " << endl;
            cout << "Pt: " << currPFJet.P4.Pt() << endl;
            cout << "Eta: " << currPFJet.P4.Eta() << endl;
            cout << "Phi: " << currPFJet.P4.Phi() << endl;
            cout << "Energy: " << currPFJet.P4.E() << endl;
        }
        patsyJetLVSmear = SmearJet(&currPFJet.P4, vecGoodGenJets, currPFJet.isGenJetMatched, currPFJet.dEnRecoGen, shiftHistSmearFactor, shiftHistLowPt, inputTF1Vec, smearBy, smearSystShift, levelVerbosity);
        /*
        cout << "patsyJetLVSmear px " << patsyJetLVSmear.Px() << endl;
        cout << "patsyJetLVSmear py " << patsyJetLVSmear.Py() << endl;
        cout << "patsyJetLVSmear pz " << patsyJetLVSmear.Pz() << endl;
        cout << "patsyJetLVSmear E " << patsyJetLVSmear.E() << endl;    
         */
        currJetPt = patsyJetLVSmear.Pt();
        currJetEta = patsyJetLVSmear.Eta();
        if (whichSystToDo == 1 && currJetPt > 9.0) {
            JetFactor = JESUncertFactor(shiftHistJetES, whichSystCase, currJetPt, currJetEta);
            patsyJetLVSmear *= JetFactor;
        }
        DeltaEnergyHist->Fill(patsyJetLVSmear.E() - currPFJet.P4.E());
        currPFJet.P4 = patsyJetLVSmear;
        outSmearJets->push_back(currPFJet);
    }
    return outSmearJets;
}

inline EventJetInfo JetKinematicsCut(vector<PFJet> * vecGoodPFJets) {
    int NJets = 0;
    int NBJets = 0;
    float HT = 0;
    int indexJet0 = -1; int indexJet1 = -1;
    int indexBtagJet0 = -1; int indexBtagJet1 = -1;
    int subIndexBtagJet0 = -1; int subIndexBtagJet1 = -1;
    float currJetPt, currJetEta;
    float BTagWP = 0.679;  //CSV Middle working point, see (remove underscore in address): h_ttps://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    float JetPtCut = 30, JetEtaCut = 2.4;
    unsigned int numGoodJets = vecGoodPFJets->size();
    for (unsigned int iJet = 0; iJet < numGoodJets; ++iJet) {
        currJetPt = vecGoodPFJets->at(iJet).P4.Pt();
        currJetEta = vecGoodPFJets->at(iJet).P4.Eta();
        if (currJetPt < JetPtCut) continue;
        if (fabs(currJetEta) > JetEtaCut) continue;
        ++NJets;
        HT += currJetPt;
        if (indexJet0 < 0) {
            indexJet0 = iJet;
        }
        else if (indexJet1 < 0) {
            indexJet1 = iJet;
        }
        if (vecGoodPFJets->at(iJet).valBTagDisc > BTagWP) {
            vecGoodPFJets->at(iJet).isBJet = true;
            NBJets++;
            if (indexBtagJet0 < 0) {
                indexBtagJet0 = iJet;
                subIndexBtagJet0 = NJets - 1;
            }
            else if (indexBtagJet1 < 0) {
                indexBtagJet1 = iJet;
                subIndexBtagJet1 = NJets - 1;
            }
        }
    }
    EventJetInfo outEJI; outEJI.EJISetValsInput(HT, NJets, indexJet0, indexJet1, NBJets, indexBtagJet0, subIndexBtagJet0, indexBtagJet1, subIndexBtagJet1, vecGoodPFJets);
    return outEJI;
}
inline EventJetInfo JetKinematicsCut(vector<PFJet> * vecGoodPFJets, BTagSFUtil * inputBTagSFUtil, bool isData, bool setGenJetInfo = false) {
    int NJets = 0;
    int NBJets = 0;
    float HT = 0;
    int indexJet0 = -1; int indexJet1 = -1;
    int indexBtagJet0 = -1; int indexBtagJet1 = -1;
    int subIndexBtagJet0 = -1; int subIndexBtagJet1 = -1;
    float currJetPt, currJetEta, currJetValBTagDisc;
    int currJetPartFlav;
    float JetPtCut = 30, JetEtaCut = 2.4;
    unsigned int numGoodJets = vecGoodPFJets->size();
    for (unsigned int iJet = 0; iJet < numGoodJets; ++iJet) {
        currJetPt = vecGoodPFJets->at(iJet).P4.Pt();
        currJetEta = vecGoodPFJets->at(iJet).P4.Eta();
        if (currJetPt < JetPtCut) continue;
        if (fabs(currJetEta) > JetEtaCut) continue;
        ++NJets;
        HT += currJetPt;
        if (indexJet0 < 0) {
            indexJet0 = iJet;
        }
        else if (indexJet1 < 0) {
            indexJet1 = iJet;
        }
        currJetPartFlav = (isData) ? -999999 : vecGoodPFJets->at(iJet).partonFlavor;
        currJetValBTagDisc = vecGoodPFJets->at(iJet).valBTagDisc;
        if (inputBTagSFUtil->IsTagged(currJetValBTagDisc, currJetPartFlav, currJetPt, currJetEta, 0, 0, false)) {
            ++NBJets;
            if (indexBtagJet0 < 0) {
                indexBtagJet0 = iJet;
                subIndexBtagJet0 = NJets - 1;
            }
            else if (indexBtagJet1 < 0) {
                indexBtagJet1 = iJet;
                subIndexBtagJet1 = NJets - 1;
            }
        }       
    }
    EventJetInfo outEJI; outEJI.EJISetValsInput(HT, NJets, indexJet0, indexJet1, NBJets, indexBtagJet0, subIndexBtagJet0, indexBtagJet1, subIndexBtagJet1, vecGoodPFJets);
    if (setGenJetInfo) outEJI.EJISetGenMatchValsInput(NJets, indexJet0, indexJet1, vecGoodPFJets);
    return outEJI;
}
inline EventJetInfo JetKinematicsCutBTagSyst(vector<PFJet> * vecGoodPFJets, BTagSFUtil * inputBTagSFUtil, int SystVar, bool isSignal) {
    int NJets = 0;
    int NBJets = 0;
    float HT = 0;
    int indexJet0 = -1; int indexJet1 = -1;
    int indexBtagJet0 = -1; int indexBtagJet1 = -1;
    int subIndexBtagJet0 = -1; int subIndexBtagJet1 = -1;
    float currJetPt, currJetEta, currJetValBTagDisc;
    int currJetPartFlav;
    float JetPtCut = 30, JetEtaCut = 2.4;
        
    int FastSimSystVar = (isSignal) ? SystVar : 0;        
    unsigned int numGoodJets = vecGoodPFJets->size();
    for (unsigned int iJet = 0; iJet < numGoodJets; ++iJet) {
        currJetPt = vecGoodPFJets->at(iJet).P4.Pt();
        currJetEta = vecGoodPFJets->at(iJet).P4.Eta();
        if (currJetPt < JetPtCut) continue;
        if (fabs(currJetEta) > JetEtaCut) continue;
        ++NJets;
        HT += currJetPt;
        if (indexJet0 < 0) {
            indexJet0 = iJet;
        }
        else if (indexJet1 < 0) {
            indexJet1 = iJet;
        }
        currJetPartFlav = vecGoodPFJets->at(iJet).partonFlavor;
        currJetValBTagDisc = vecGoodPFJets->at(iJet).valBTagDisc;
        if (inputBTagSFUtil->IsTagged(currJetValBTagDisc, currJetPartFlav, currJetPt, currJetEta, SystVar, FastSimSystVar, false)) {
            ++NBJets;
            if (indexBtagJet0 < 0) {
                indexBtagJet0 = iJet;
                subIndexBtagJet0 = NJets - 1;
            }
            else if (indexBtagJet1 < 0) {
                indexBtagJet1 = iJet;
                subIndexBtagJet1 = NJets - 1;
            }
        }       
    }
    EventJetInfo outEJI; outEJI.EJISetValsInput(HT, NJets, indexJet0, indexJet1, NBJets, indexBtagJet0, subIndexBtagJet0, indexBtagJet1, subIndexBtagJet1, vecGoodPFJets);
    return outEJI;
}
/*
inline vector<HistogramT> * AddSystHists(vector<HistogramT> * inputHistTVec, vector<SystT> * inputSystTVec, TString fileInName, bool isSignal) {
    HistogramT H_Curr;
    HistogramT H_CurrNewSyst;
    SystT      Syst_Curr;
    vector<HistogramT> * histTSystVec = new vector<HistogramT>;
    bool isLepXAxis = false, isJetXAxis = false, isOtherXAxis = false;
    bool isLepYAxis = false, isJetYAxis = false, isOtherYAxis = false;
    bool isLepZAxis = false, isJetZAxis = false, isOtherZAxis = false;
    for (unsigned int i = 0; i < inputHistTVec->size(); ++i) {
        H_Curr = inputHistTVec->at(i);
        if (H_Curr.doXSyst || H_Curr.doYSyst || H_Curr.doZSyst) {
            if (H_Curr.doXSyst) {
                isLepXAxis = (TString(H_Curr.xVarKey).Contains("Lep"));
                isJetXAxis = (TString(H_Curr.xVarKey).Contains("Jet"));
                isOtherXAxis = (TString(H_Curr.xVarKey).Contains("MT2") || TString(H_Curr.xVarKey).Contains("MET"));
            }
            if (H_Curr.doYSyst) {
                isLepYAxis = (TString(H_Curr.yVarKey).Contains("Lep"));
                isJetYAxis = (TString(H_Curr.yVarKey).Contains("Jet"));
                isOtherYAxis = (TString(H_Curr.yVarKey).Contains("MT2") || TString(H_Curr.yVarKey).Contains("MET"));            
            }
            if (H_Curr.doZSyst) {
                isLepZAxis = (TString(H_Curr.zVarKey).Contains("Lep"));
                isJetZAxis = (TString(H_Curr.zVarKey).Contains("Jet"));
                isOtherZAxis = (TString(H_Curr.zVarKey).Contains("MT2") || TString(H_Curr.zVarKey).Contains("MET"));            
            }
            for (unsigned int j = 0; j < inputSystTVec->size(); ++j) {
                Syst_Curr = inputSystTVec->at(j);
                switch (Syst_Curr.whichSystType) {
                    case 3:
                        // temporary one for MT2ll
                        if (Syst_Curr.name.Contains("_MT2llShift")) {
                            if (!(H_Curr.name.Contains("MT2ll"))) continue;
                        }
                        if (Syst_Curr.name.Contains("_MT2UncESSpec")) {
                            if (!(H_Curr.name.Contains("MT2"))) continue;
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
                if (Syst_Curr.name.Contains("_BTagSF")) {
                    H_CurrNewSyst.xVarKey += Syst_Curr.systVarKey;
                    H_CurrNewSyst.yVarKey += Syst_Curr.systVarKey;
                    H_CurrNewSyst.zVarKey += Syst_Curr.systVarKey;
                    cout << "H_CurrNewSyst.xVarKey " << H_CurrNewSyst.xVarKey << endl;
                }
                else {
                    if (isOtherXAxis) {
                        H_CurrNewSyst.xVarKey += Syst_Curr.systVarKey;
                    }
                    else if (isLepXAxis && Syst_Curr.whichSystType == 1) {
                        H_CurrNewSyst.xVarKey += Syst_Curr.systVarKey;
                    }                
                    else if (isJetXAxis && Syst_Curr.whichSystType == 2) {
                        H_CurrNewSyst.xVarKey += Syst_Curr.systVarKey;
                    }                    
                    if (isOtherYAxis) {
                        H_CurrNewSyst.yVarKey += Syst_Curr.systVarKey;
                    }
                    else if (isLepYAxis && Syst_Curr.whichSystType == 1) {
                        H_CurrNewSyst.yVarKey += Syst_Curr.systVarKey;
                    }                
                    else if (isJetYAxis && Syst_Curr.whichSystType == 2) {
                        H_CurrNewSyst.yVarKey += Syst_Curr.systVarKey;
                    }
                    if (isOtherZAxis) {
                        H_CurrNewSyst.zVarKey += Syst_Curr.systVarKey;
                    }
                    else if (isLepZAxis && Syst_Curr.whichSystType == 1) {
                        H_CurrNewSyst.zVarKey += Syst_Curr.systVarKey;
                    }                
                    else if (isJetZAxis && Syst_Curr.whichSystType == 2) {
                        H_CurrNewSyst.zVarKey += Syst_Curr.systVarKey;
                    }
                }                
                histTSystVec->push_back(H_CurrNewSyst);
            }
        }
    }
    return histTSystVec;
}
*/
inline vector<HistogramT> * AddSystHists(vector<HistogramT> * inputHistTVec, vector<SystT> * inputSystTVec, TString fileInName, bool isSignal) {
    HistogramT H_Curr;
    HistogramT H_CurrNewSyst;
    SystT      Syst_Curr;
    vector<HistogramT> * histTSystVec = new vector<HistogramT>;
    for (unsigned int i = 0; i < inputHistTVec->size(); ++i) {
        H_Curr = inputHistTVec->at(i);
        if (H_Curr.doXSyst || H_Curr.doYSyst || H_Curr.doZSyst) {
            for (unsigned int j = 0; j < inputSystTVec->size(); ++j) {
                Syst_Curr = inputSystTVec->at(j);
                H_CurrNewSyst = H_Curr;
                if (Syst_Curr.name.Contains("Smear") && !(H_CurrNewSyst.name.Contains("Smear"))) continue;
                if (Syst_Curr.name.Contains("_MT2llShift") && !(H_Curr.name.Contains("MT2ll"))) continue;
                if (Syst_Curr.name.Contains("_MT2UncESSpec") && !(H_Curr.name.Contains("MT2"))) continue;
                if (Syst_Curr.name.Contains("genTopRW") && !(fileInName.Contains("TT") || fileInName.Contains("ttbar"))) continue;
                if (Syst_Curr.name.Contains("genStopXSec") && !isSignal) continue;
                H_CurrNewSyst.name += Syst_Curr.name;
                if (H_Curr.doXSyst) H_CurrNewSyst.xVarKey += Syst_Curr.systVarKey;
                if (H_Curr.doYSyst) H_CurrNewSyst.yVarKey += Syst_Curr.systVarKey;
                if (H_Curr.doZSyst) H_CurrNewSyst.zVarKey += Syst_Curr.systVarKey;            
                histTSystVec->push_back(H_CurrNewSyst);
            }
        }
    }
    return histTSystVec;
}
inline vector<HistogramT> * OneDeeHistTVec() {
    const double PI = 3.14159265;
    int EnergyPtBinN = 40;
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
    
    HistogramT H_RelLeadLepPFIso; H_RelLeadLepPFIso.name = "h_RelLeadLepPFIso";
    H_RelLeadLepPFIso.xLabel = "Lead Lepton Relative PF Isolation"; H_RelLeadLepPFIso.xBinN = 32; H_RelLeadLepPFIso.xMin = 0.; H_RelLeadLepPFIso.xMax = 0.16;
    H_RelLeadLepPFIso.yLabel = "Number of Events / ";
    H_RelLeadLepPFIso.yLabel += "NUM"; H_RelLeadLepPFIso.yLabel += " of Iso.";
    H_RelLeadLepPFIso.xVarKey = "lep0RelPFIso";
    H_RelLeadLepPFIso.doXSyst = false;
    
    HistogramT H_RelSubLepPFIso; H_RelSubLepPFIso.name = "h_RelSubLepPFIso";
    H_RelSubLepPFIso.xLabel = "Sub-Lead Lepton Relative PF Isolation"; H_RelSubLepPFIso.xBinN = 32; H_RelSubLepPFIso.xMin = 0.; H_RelSubLepPFIso.xMax = 0.16;
    H_RelSubLepPFIso.yLabel = "Number of Events / ";
    H_RelSubLepPFIso.yLabel += "NUM"; H_RelSubLepPFIso.yLabel += " of Iso.";
    H_RelSubLepPFIso.xVarKey = "lep1RelPFIso";
    H_RelSubLepPFIso.doXSyst = false;

    HistogramT H_CutFlow; H_CutFlow.name = "h_ChannelCutFlow";
    H_CutFlow.xLabel = "Cut Stage"; H_CutFlow.xBinN = 4; H_CutFlow.xMin = 0.5; H_CutFlow.xMax = 4.5;
    H_CutFlow.yLabel = "Number of Events passing CutStage";
    H_CutFlow.xVarKey = "CutFlowEntry";
    H_CutFlow.doXSyst = false;
    
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
    
    
    HistogramT H_MT2ll_DPhiLep0Lep1Close; H_MT2ll_DPhiLep0Lep1Close.name = "h_MT2ll_DPhiLep0Lep1Close"; 
    H_MT2ll_DPhiLep0Lep1Close.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Close.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Close.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Close.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Close.xMax - H_MT2ll_DPhiLep0Lep1Close.xMin) / (float) H_MT2ll_DPhiLep0Lep1Close.xBinN;
    H_MT2ll_DPhiLep0Lep1Close.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Close.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Close.xVarKey = "MT2ll";
    H_MT2ll_DPhiLep0Lep1Close.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Mid; H_MT2ll_DPhiLep0Lep1Mid.name = "h_MT2ll_DPhiLep0Lep1Mid"; 
    H_MT2ll_DPhiLep0Lep1Mid.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Mid.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Mid.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Mid.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Mid.xMax - H_MT2ll_DPhiLep0Lep1Mid.xMin) / (float) H_MT2ll_DPhiLep0Lep1Mid.xBinN;
    H_MT2ll_DPhiLep0Lep1Mid.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Mid.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Mid.xVarKey = "MT2ll";
    H_MT2ll_DPhiLep0Lep1Mid.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Far; H_MT2ll_DPhiLep0Lep1Far.name = "h_MT2ll_DPhiLep0Lep1Far"; 
    H_MT2ll_DPhiLep0Lep1Far.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Far.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Far.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Far.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Far.xMax - H_MT2ll_DPhiLep0Lep1Far.xMin) / (float) H_MT2ll_DPhiLep0Lep1Far.xBinN;
    H_MT2ll_DPhiLep0Lep1Far.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Far.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Far.xVarKey = "MT2ll";
    H_MT2ll_DPhiLep0Lep1Far.doXSyst = true;
    
    
    HistogramT H_MT2lb; H_MT2lb.name = "h_MT2lb"; 
    H_MT2lb.xLabel = "MT2lb [GeV]"; H_MT2lb.xBinN = METBinN; H_MT2lb.xMin = METBinLB; H_MT2lb.xMax = METBinUB; 
    H_MT2lb.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb.xMax - H_MT2lb.xMin) / (float) H_MT2lb.xBinN;
    H_MT2lb.yLabel += "NUM"; H_MT2lb.yLabel += " GeV";
    H_MT2lb.xVarKey = "MT2lb";
    H_MT2lb.doXSyst = true;
    HistogramT H_MT2lb_DPhiBLep0BLep1Close; H_MT2lb_DPhiBLep0BLep1Close.name = "h_MT2lb_DPhiBLep0BLep1Close"; 
    H_MT2lb_DPhiBLep0BLep1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Close.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Close.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Close.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Close.xMax - H_MT2lb_DPhiBLep0BLep1Close.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Close.xBinN;
    H_MT2lb_DPhiBLep0BLep1Close.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Close.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Close.xVarKey = "MT2lb";
    H_MT2lb_DPhiBLep0BLep1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Mid; H_MT2lb_DPhiBLep0BLep1Mid.name = "h_MT2lb_DPhiBLep0BLep1Mid"; 
    H_MT2lb_DPhiBLep0BLep1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Mid.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Mid.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Mid.xMax - H_MT2lb_DPhiBLep0BLep1Mid.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Mid.xBinN;
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Mid.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Mid.xVarKey = "MT2lb";
    H_MT2lb_DPhiBLep0BLep1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Far; H_MT2lb_DPhiBLep0BLep1Far.name = "h_MT2lb_DPhiBLep0BLep1Far"; 
    H_MT2lb_DPhiBLep0BLep1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Far.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Far.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Far.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Far.xMax - H_MT2lb_DPhiBLep0BLep1Far.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Far.xBinN;
    H_MT2lb_DPhiBLep0BLep1Far.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Far.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Far.xVarKey = "MT2lb";
    H_MT2lb_DPhiBLep0BLep1Far.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Close; H_MT2lb_DPhiJet0Jet1Close.name = "h_MT2lb_DPhiJet0Jet1Close"; 
    H_MT2lb_DPhiJet0Jet1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Close.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Close.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Close.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Close.xMax - H_MT2lb_DPhiJet0Jet1Close.xMin) / (float) H_MT2lb_DPhiJet0Jet1Close.xBinN;
    H_MT2lb_DPhiJet0Jet1Close.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Close.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Close.xVarKey = "MT2lb";
    H_MT2lb_DPhiJet0Jet1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Mid; H_MT2lb_DPhiJet0Jet1Mid.name = "h_MT2lb_DPhiJet0Jet1Mid"; 
    H_MT2lb_DPhiJet0Jet1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Mid.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Mid.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Mid.xMax - H_MT2lb_DPhiJet0Jet1Mid.xMin) / (float) H_MT2lb_DPhiJet0Jet1Mid.xBinN;
    H_MT2lb_DPhiJet0Jet1Mid.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Mid.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Mid.xVarKey = "MT2lb";
    H_MT2lb_DPhiJet0Jet1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Far; H_MT2lb_DPhiJet0Jet1Far.name = "h_MT2lb_DPhiJet0Jet1Far"; 
    H_MT2lb_DPhiJet0Jet1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Far.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Far.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Far.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Far.xMax - H_MT2lb_DPhiJet0Jet1Far.xMin) / (float) H_MT2lb_DPhiJet0Jet1Far.xBinN;
    H_MT2lb_DPhiJet0Jet1Far.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Far.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Far.xVarKey = "MT2lb";
    H_MT2lb_DPhiJet0Jet1Far.doXSyst = true;
    
    
    HistogramT H_MT2lbCont; H_MT2lbCont.name = "h_MT2lbControl"; 
    H_MT2lbCont.xLabel = "MT2lb [GeV]"; H_MT2lbCont.xBinN = 34; H_MT2lbCont.xMin = 0; H_MT2lbCont.xMax = 170; 
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
    
    HistogramT H_MET_div_Meff; H_MET_div_Meff.name = "h_MET_div_Meff";
    H_MET_div_Meff.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff.xBinN = 50; H_MET_div_Meff.xMin = 0.; H_MET_div_Meff.xMax = 1.;
    H_MET_div_Meff.yLabel = "Number of Events / "; H_MET_div_Meff.yLabel += "NUM";
    H_MET_div_Meff.xVarKey = "METdivMeff";
    H_MET_div_Meff.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut80; H_MET_div_Meff_PassMT2llCut80.name = "h_MET_div_Meff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut80.xBinN = 50; H_MET_div_Meff_PassMT2llCut80.xMin = 0.; H_MET_div_Meff_PassMT2llCut80.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut80.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut80.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut80.xVarKey = "METdivMeff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut90; H_MET_div_Meff_PassMT2llCut90.name = "h_MET_div_Meff_PassMT2llCut90";
    H_MET_div_Meff_PassMT2llCut90.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut90.xBinN = 50; H_MET_div_Meff_PassMT2llCut90.xMin = 0.; H_MET_div_Meff_PassMT2llCut90.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut90.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut90.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut90.xVarKey = "METdivMeff_PassMT2llCut90";
    H_MET_div_Meff_PassMT2llCut90.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut100; H_MET_div_Meff_PassMT2llCut100.name = "h_MET_div_Meff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut100.xBinN = 50; H_MET_div_Meff_PassMT2llCut100.xMin = 0.; H_MET_div_Meff_PassMT2llCut100.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut100.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut100.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut100.xVarKey = "METdivMeff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut110; H_MET_div_Meff_PassMT2llCut110.name = "h_MET_div_Meff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut110.xBinN = 50; H_MET_div_Meff_PassMT2llCut110.xMin = 0.; H_MET_div_Meff_PassMT2llCut110.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut110.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut110.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut110.xVarKey = "METdivMeff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut120; H_MET_div_Meff_PassMT2llCut120.name = "h_MET_div_Meff_PassMT2llCut120";
    H_MET_div_Meff_PassMT2llCut120.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut120.xBinN = 50; H_MET_div_Meff_PassMT2llCut120.xMin = 0.; H_MET_div_Meff_PassMT2llCut120.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut120.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut120.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut120.xVarKey = "METdivMeff_PassMT2llCut120";
    H_MET_div_Meff_PassMT2llCut120.doXSyst = true;
    
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
    histVec_1D->push_back(H_RelLeadLepPFIso); histVec_1D->push_back(H_RelSubLepPFIso); histVec_1D->push_back(H_CutFlow);
    histVec_1D->push_back(H_leadLepPt); histVec_1D->push_back(H_leadLepEta); histVec_1D->push_back(H_subLepPt); histVec_1D->push_back(H_subLepEta);
    histVec_1D->push_back(H_leadJetPt); histVec_1D->push_back(H_leadJetEta); histVec_1D->push_back(H_subJetPt); histVec_1D->push_back(H_subJetEta);
    histVec_1D->push_back(H_leadBJetPt); histVec_1D->push_back(H_leadBJetEta); histVec_1D->push_back(H_leadBJetEn); histVec_1D->push_back(H_subBJetPt); histVec_1D->push_back(H_subBJetEta); histVec_1D->push_back(H_subBJetEn);
    histVec_1D->push_back(H_diLepPt); histVec_1D->push_back(H_diLepInvMass); //histVec_1D->push_back(H_diLepEta); histVec_1D->push_back(H_diLepPhi); 
    histVec_1D->push_back(H_MT2ll); histVec_1D->push_back(H_MT2lb);
    histVec_1D->push_back(H_MT2llCont); histVec_1D->push_back(H_MT2lbCont);
    histVec_1D->push_back(H_PassMT2llCut80); histVec_1D->push_back(H_PassMT2llCut90); histVec_1D->push_back(H_PassMT2llCut100); histVec_1D->push_back(H_PassMT2llCut110); histVec_1D->push_back(H_PassMT2llCut120);
    histVec_1D->push_back(H_MT2ll_DPhiZMETClose); histVec_1D->push_back(H_MT2ll_DPhiZMETMid); histVec_1D->push_back(H_MT2ll_DPhiZMETFar);    
    histVec_1D->push_back(H_MT2ll_DPhiLep0Lep1Close); histVec_1D->push_back(H_MT2ll_DPhiLep0Lep1Mid); histVec_1D->push_back(H_MT2ll_DPhiLep0Lep1Far);    
    histVec_1D->push_back(H_MT2lb_DPhiBLep0BLep1Close); histVec_1D->push_back(H_MT2lb_DPhiBLep0BLep1Mid); histVec_1D->push_back(H_MT2lb_DPhiBLep0BLep1Far);
    histVec_1D->push_back(H_MT2lb_DPhiJet0Jet1Close); histVec_1D->push_back(H_MT2lb_DPhiJet0Jet1Mid); histVec_1D->push_back(H_MT2lb_DPhiJet0Jet1Far);
    histVec_1D->push_back(H_MET); histVec_1D->push_back(H_METX); histVec_1D->push_back(H_METY); histVec_1D->push_back(H_METPhi); histVec_1D->push_back(H_METPhi_noCorr); histVec_1D->push_back(H_METX_noPhiCorr); histVec_1D->push_back(H_METY_noPhiCorr);
    histVec_1D->push_back(H_NJets); histVec_1D->push_back(H_NJetswBTag);
    histVec_1D->push_back(H_diJetPt); histVec_1D->push_back(H_diJetInvMass); //histVec_1D->push_back(H_diJetEta); histVec_1D->push_back(H_diJetPhi); 
    histVec_1D->push_back(H_diBJetPt); histVec_1D->push_back(H_diBJetInvMass); //histVec_1D->push_back(H_diBJetEta); histVec_1D->push_back(H_diBJetPhi); 
    histVec_1D->push_back(H_DeltaPhiLep0Lep1); 
    histVec_1D->push_back(H_DeltaPhiLep0MET); histVec_1D->push_back(H_DeltaPhiLep1MET); 
    histVec_1D->push_back(H_DeltaPhiLep0MET_PreCorr); histVec_1D->push_back(H_DeltaPhiLep1MET_PreCorr); 
    histVec_1D->push_back(H_DeltaPhiZMET); histVec_1D->push_back(H_DeltaPhiZMET_PreCorr); 
    histVec_1D->push_back(H_DeltaPhiLep0Jet0); histVec_1D->push_back(H_DeltaPhiLep0Jet1);
    histVec_1D->push_back(H_DeltaPhiLep0BJet0); histVec_1D->push_back(H_DeltaPhiLep0BJet1);
    histVec_1D->push_back(H_DeltaPhiJet0BJet0); histVec_1D->push_back(H_DeltaPhiJet1BJet1); histVec_1D->push_back(H_DeltaPhiJet1BJet0);
    histVec_1D->push_back(H_diLepEMinJetE); histVec_1D->push_back(H_HT); 
    histVec_1D->push_back(H_MET_div_Meff); 
    histVec_1D->push_back(H_MET_div_Meff_PassMT2llCut80); histVec_1D->push_back(H_MET_div_Meff_PassMT2llCut90); histVec_1D->push_back(H_MET_div_Meff_PassMT2llCut100); histVec_1D->push_back(H_MET_div_Meff_PassMT2llCut110); histVec_1D->push_back(H_MET_div_Meff_PassMT2llCut120); 
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
    
    HistogramT H_MT2lb_vs_DeltaPhiLepB0LepB1; H_MT2lb_vs_DeltaPhiLepB0LepB1.name = "h_MT2lb_vs_DeltaPhiLepB0LepB1";
    H_MT2lb_vs_DeltaPhiLepB0LepB1.xLabel = "MT2lb [GeV]"; H_MT2lb_vs_DeltaPhiLepB0LepB1.xBinN = METBinN; H_MT2lb_vs_DeltaPhiLepB0LepB1.xMin = METBinLB; H_MT2lb_vs_DeltaPhiLepB0LepB1.xMax = METBinUB;
    H_MT2lb_vs_DeltaPhiLepB0LepB1.yLabel = "#Delta #phi_{lb lb}"; H_MT2lb_vs_DeltaPhiLepB0LepB1.yBinN = PhiBinN; H_MT2lb_vs_DeltaPhiLepB0LepB1.yMin = 0; H_MT2lb_vs_DeltaPhiLepB0LepB1.yMax = PI;
    H_MT2lb_vs_DeltaPhiLepB0LepB1.xVarKey = "MT2lb";
    H_MT2lb_vs_DeltaPhiLepB0LepB1.yVarKey = "DPhiLepB0LepB1";
    H_MT2lb_vs_DeltaPhiLepB0LepB1.doXSyst = true;
    H_MT2lb_vs_DeltaPhiLepB0LepB1.doYSyst = true;
    
    
    
    HistogramT H_MT2ll_vs_MT2lb; H_MT2ll_vs_MT2lb.name = "h_MT2ll_vs_MT2lb";
    H_MT2ll_vs_MT2lb.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_MT2lb.xBinN = METBinN; H_MT2ll_vs_MT2lb.xMin = METBinLB; H_MT2ll_vs_MT2lb.xMax = METBinUB;
    H_MT2ll_vs_MT2lb.yLabel = "MT2lb [GeV]"; H_MT2ll_vs_MT2lb.yBinN = METBinN; H_MT2ll_vs_MT2lb.yMin = METBinLB; H_MT2ll_vs_MT2lb.yMax = METBinUB;
    H_MT2ll_vs_MT2lb.xVarKey = "MT2ll";
    H_MT2ll_vs_MT2lb.yVarKey = "MT2lb";
    H_MT2ll_vs_MT2lb.doXSyst = true;
    H_MT2ll_vs_MT2lb.doYSyst = true;
    /*
    HistogramT H_MT2llControl_vs_DeltaPhiLep0Lep1; H_MT2llControl_vs_DeltaPhiLep0Lep1.name = "h_MT2llControl_vs_DeltaPhiLep0Lep1";
    H_MT2llControl_vs_DeltaPhiLep0Lep1.xLabel = "MT2_{ll} [GeV]"; H_MT2llControl_vs_DeltaPhiLep0Lep1.xBinN = 20; H_MT2llControl_vs_DeltaPhiLep0Lep1.xMin = 0; H_MT2llControl_vs_DeltaPhiLep0Lep1.xMax = 80;
    H_MT2llControl_vs_DeltaPhiLep0Lep1.yLabel = "#Delta #phi_{ll}"; H_MT2llControl_vs_DeltaPhiLep0Lep1.yBinN = PhiBinN; H_MT2llControl_vs_DeltaPhiLep0Lep1.yMin = 0; H_MT2llControl_vs_DeltaPhiLep0Lep1.yMax = PI;
    H_MT2llControl_vs_DeltaPhiLep0Lep1.xVarKey = "MT2ll";
    H_MT2llControl_vs_DeltaPhiLep0Lep1.yVarKey = "DPhiLep0Lep1";
    H_MT2llControl_vs_DeltaPhiLep0Lep1.doXSyst = true;
    H_MT2llControl_vs_DeltaPhiLep0Lep1.doYSyst = true;
    */
    /*    
    HistogramT H_MT2lbControl_vs_DeltaPhiLepB0LepB1; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.name = "h_MT2lbControl_vs_DeltaPhiLepB0LepB1";
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xLabel = "MT2lb [GeV]"; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xBinN = 43; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xMin = 0; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xMax = 172;
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yLabel = "#Delta #phi_{lb lb}"; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yBinN = PhiBinN; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yMin = 0; H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yMax = PI;
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.xVarKey = "MT2lb";
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.yVarKey = "DPhiLepB0LepB1";
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.doXSyst = true;
    H_MT2lbControl_vs_DeltaPhiLepB0LepB1.doYSyst = true;
     */
    //push the 2D histograms structures into a vector for eventual use in booking histograms
    histVec_2D->push_back(H_METX_vs_nVtx); histVec_2D->push_back(H_METY_vs_nVtx);
    histVec_2D->push_back(H_METX_vs_nVtx_noPhiCorr); histVec_2D->push_back(H_METY_vs_nVtx_noPhiCorr);
    histVec_2D->push_back(H_MT2ll_vs_DeltaPhiLep0Lep1); //histVec_2D->push_back(H_MT2llControl_vs_DeltaPhiLep0Lep1);
    histVec_2D->push_back(H_MT2lb_vs_DeltaPhiLepB0LepB1); // histVec_2D->push_back(H_MT2lbControl_vs_DeltaPhiLepB0LepB1);
    histVec_2D->push_back(H_MT2ll_vs_MT2lb);
    return histVec_2D;
}

inline vector<HistogramT> * ThreeDeeHistTVec() {
    
    vector<HistogramT> * histVec_3D = new vector<HistogramT>;
    
    const double PI = 3.14159265;
    //    int EnergyPtBinN = 400;
    //    int EtaBinN      = 200;
    int PhiBinN      = 100;
    int METBinN      = 40;    
//    int METXYBinN    = 200;    
    int NJetsBinN    = 11;
    int nVtxBinN     = 35;
    
    //    float EnergyPtBinLB = 0;
    //    float EnergyPtBinUB = 400;
    //    float EtaBinLB      = -6;
    //    float EtaBinUB      = 6;
    float METBinLB      = 0;
    float METBinUB      = 400;
//    float METXYBinLB    = -200;
//    float METXYBinUB    = 200;
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
    
    /*
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
    */
    
    
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
//    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NBJets);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx1to10);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx11to20);
    histVec_3D->push_back(H_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx21to30);
    return histVec_3D;
}

inline vector<HistogramT> * OneDeeHistTVecStatusReport() {
    const double PI = 3.14159265;
    int EnergyPtBinN = 40;
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
    
    HistogramT H_RelLeadLepPFIso; H_RelLeadLepPFIso.name = "h_RelLeadLepPFIso";
    H_RelLeadLepPFIso.xLabel = "Lead Lepton Relative PF Isolation"; H_RelLeadLepPFIso.xBinN = 32; H_RelLeadLepPFIso.xMin = 0.; H_RelLeadLepPFIso.xMax = 0.16;
    H_RelLeadLepPFIso.yLabel = "Number of Events / ";
    H_RelLeadLepPFIso.yLabel += "NUM"; H_RelLeadLepPFIso.yLabel += " of Iso.";
    H_RelLeadLepPFIso.xVarKey = "lep0RelPFIso";
    H_RelLeadLepPFIso.doXSyst = false;
    
    HistogramT H_RelSubLepPFIso; H_RelSubLepPFIso.name = "h_RelSubLepPFIso";
    H_RelSubLepPFIso.xLabel = "Sub-Lead Lepton Relative PF Isolation"; H_RelSubLepPFIso.xBinN = 32; H_RelSubLepPFIso.xMin = 0.; H_RelSubLepPFIso.xMax = 0.16;
    H_RelSubLepPFIso.yLabel = "Number of Events / ";
    H_RelSubLepPFIso.yLabel += "NUM"; H_RelSubLepPFIso.yLabel += " of Iso.";
    H_RelSubLepPFIso.xVarKey = "lep1RelPFIso";
    H_RelSubLepPFIso.doXSyst = false;
    
    HistogramT H_CutFlow; H_CutFlow.name = "h_ChannelCutFlow";
    H_CutFlow.xLabel = "Cut Stage"; H_CutFlow.xBinN = 4; H_CutFlow.xMin = 0.5; H_CutFlow.xMax = 4.5;
    H_CutFlow.yLabel = "Number of Events passing CutStage";
    H_CutFlow.xVarKey = "CutFlowEntry";
    H_CutFlow.doXSyst = false;
    
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
    
    
    HistogramT H_MT2ll_DPhiLep0Lep1Close; H_MT2ll_DPhiLep0Lep1Close.name = "h_MT2ll_DPhiLep0Lep1Close"; 
    H_MT2ll_DPhiLep0Lep1Close.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Close.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Close.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Close.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Close.xMax - H_MT2ll_DPhiLep0Lep1Close.xMin) / (float) H_MT2ll_DPhiLep0Lep1Close.xBinN;
    H_MT2ll_DPhiLep0Lep1Close.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Close.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Close.xVarKey = "MT2ll";
    H_MT2ll_DPhiLep0Lep1Close.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Mid; H_MT2ll_DPhiLep0Lep1Mid.name = "h_MT2ll_DPhiLep0Lep1Mid"; 
    H_MT2ll_DPhiLep0Lep1Mid.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Mid.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Mid.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Mid.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Mid.xMax - H_MT2ll_DPhiLep0Lep1Mid.xMin) / (float) H_MT2ll_DPhiLep0Lep1Mid.xBinN;
    H_MT2ll_DPhiLep0Lep1Mid.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Mid.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Mid.xVarKey = "MT2ll";
    H_MT2ll_DPhiLep0Lep1Mid.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Far; H_MT2ll_DPhiLep0Lep1Far.name = "h_MT2ll_DPhiLep0Lep1Far"; 
    H_MT2ll_DPhiLep0Lep1Far.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Far.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Far.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Far.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Far.xMax - H_MT2ll_DPhiLep0Lep1Far.xMin) / (float) H_MT2ll_DPhiLep0Lep1Far.xBinN;
    H_MT2ll_DPhiLep0Lep1Far.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Far.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Far.xVarKey = "MT2ll";
    H_MT2ll_DPhiLep0Lep1Far.doXSyst = true;
    
    
    HistogramT H_MT2lb; H_MT2lb.name = "h_MT2lb"; 
    H_MT2lb.xLabel = "MT2lb [GeV]"; H_MT2lb.xBinN = METBinN; H_MT2lb.xMin = METBinLB; H_MT2lb.xMax = METBinUB; 
    H_MT2lb.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb.xMax - H_MT2lb.xMin) / (float) H_MT2lb.xBinN;
    H_MT2lb.yLabel += "NUM"; H_MT2lb.yLabel += " GeV";
    H_MT2lb.xVarKey = "MT2lb";
    H_MT2lb.doXSyst = true;
    HistogramT H_MT2lb_DPhiBLep0BLep1Close; H_MT2lb_DPhiBLep0BLep1Close.name = "h_MT2lb_DPhiBLep0BLep1Close"; 
    H_MT2lb_DPhiBLep0BLep1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Close.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Close.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Close.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Close.xMax - H_MT2lb_DPhiBLep0BLep1Close.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Close.xBinN;
    H_MT2lb_DPhiBLep0BLep1Close.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Close.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Close.xVarKey = "MT2lb";
    H_MT2lb_DPhiBLep0BLep1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Mid; H_MT2lb_DPhiBLep0BLep1Mid.name = "h_MT2lb_DPhiBLep0BLep1Mid"; 
    H_MT2lb_DPhiBLep0BLep1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Mid.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Mid.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Mid.xMax - H_MT2lb_DPhiBLep0BLep1Mid.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Mid.xBinN;
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Mid.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Mid.xVarKey = "MT2lb";
    H_MT2lb_DPhiBLep0BLep1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Far; H_MT2lb_DPhiBLep0BLep1Far.name = "h_MT2lb_DPhiBLep0BLep1Far"; 
    H_MT2lb_DPhiBLep0BLep1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Far.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Far.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Far.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Far.xMax - H_MT2lb_DPhiBLep0BLep1Far.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Far.xBinN;
    H_MT2lb_DPhiBLep0BLep1Far.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Far.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Far.xVarKey = "MT2lb";
    H_MT2lb_DPhiBLep0BLep1Far.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Close; H_MT2lb_DPhiJet0Jet1Close.name = "h_MT2lb_DPhiJet0Jet1Close"; 
    H_MT2lb_DPhiJet0Jet1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Close.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Close.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Close.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Close.xMax - H_MT2lb_DPhiJet0Jet1Close.xMin) / (float) H_MT2lb_DPhiJet0Jet1Close.xBinN;
    H_MT2lb_DPhiJet0Jet1Close.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Close.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Close.xVarKey = "MT2lb";
    H_MT2lb_DPhiJet0Jet1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Mid; H_MT2lb_DPhiJet0Jet1Mid.name = "h_MT2lb_DPhiJet0Jet1Mid"; 
    H_MT2lb_DPhiJet0Jet1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Mid.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Mid.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Mid.xMax - H_MT2lb_DPhiJet0Jet1Mid.xMin) / (float) H_MT2lb_DPhiJet0Jet1Mid.xBinN;
    H_MT2lb_DPhiJet0Jet1Mid.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Mid.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Mid.xVarKey = "MT2lb";
    H_MT2lb_DPhiJet0Jet1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Far; H_MT2lb_DPhiJet0Jet1Far.name = "h_MT2lb_DPhiJet0Jet1Far"; 
    H_MT2lb_DPhiJet0Jet1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Far.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Far.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Far.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Far.xMax - H_MT2lb_DPhiJet0Jet1Far.xMin) / (float) H_MT2lb_DPhiJet0Jet1Far.xBinN;
    H_MT2lb_DPhiJet0Jet1Far.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Far.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Far.xVarKey = "MT2lb";
    H_MT2lb_DPhiJet0Jet1Far.doXSyst = true;
    
    
    HistogramT H_MT2lbCont; H_MT2lbCont.name = "h_MT2lbControl"; 
    H_MT2lbCont.xLabel = "MT2lb [GeV]"; H_MT2lbCont.xBinN = 34; H_MT2lbCont.xMin = 0; H_MT2lbCont.xMax = 170; 
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
    
    HistogramT H_MET_div_Meff; H_MET_div_Meff.name = "h_MET_div_Meff";
    H_MET_div_Meff.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff.xBinN = 50; H_MET_div_Meff.xMin = 0.; H_MET_div_Meff.xMax = 1.;
    H_MET_div_Meff.yLabel = "Number of Events / "; H_MET_div_Meff.yLabel += "NUM";
    H_MET_div_Meff.xVarKey = "METdivMeff";
    H_MET_div_Meff.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut80; H_MET_div_Meff_PassMT2llCut80.name = "h_MET_div_Meff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut80.xBinN = 50; H_MET_div_Meff_PassMT2llCut80.xMin = 0.; H_MET_div_Meff_PassMT2llCut80.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut80.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut80.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut80.xVarKey = "METdivMeff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut90; H_MET_div_Meff_PassMT2llCut90.name = "h_MET_div_Meff_PassMT2llCut90";
    H_MET_div_Meff_PassMT2llCut90.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut90.xBinN = 50; H_MET_div_Meff_PassMT2llCut90.xMin = 0.; H_MET_div_Meff_PassMT2llCut90.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut90.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut90.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut90.xVarKey = "METdivMeff_PassMT2llCut90";
    H_MET_div_Meff_PassMT2llCut90.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut100; H_MET_div_Meff_PassMT2llCut100.name = "h_MET_div_Meff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut100.xBinN = 50; H_MET_div_Meff_PassMT2llCut100.xMin = 0.; H_MET_div_Meff_PassMT2llCut100.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut100.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut100.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut100.xVarKey = "METdivMeff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut110; H_MET_div_Meff_PassMT2llCut110.name = "h_MET_div_Meff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut110.xBinN = 50; H_MET_div_Meff_PassMT2llCut110.xMin = 0.; H_MET_div_Meff_PassMT2llCut110.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut110.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut110.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut110.xVarKey = "METdivMeff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut120; H_MET_div_Meff_PassMT2llCut120.name = "h_MET_div_Meff_PassMT2llCut120";
    H_MET_div_Meff_PassMT2llCut120.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut120.xBinN = 50; H_MET_div_Meff_PassMT2llCut120.xMin = 0.; H_MET_div_Meff_PassMT2llCut120.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut120.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut120.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut120.xVarKey = "METdivMeff_PassMT2llCut120";
    H_MET_div_Meff_PassMT2llCut120.doXSyst = true;
    
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
    histVec_1D->push_back(H_diLepPt); histVec_1D->push_back(H_diLepInvMass); //histVec_1D->push_back(H_diLepEta); histVec_1D->push_back(H_diLepPhi); 
    histVec_1D->push_back(H_MT2ll); histVec_1D->push_back(H_MT2lb);
    histVec_1D->push_back(H_MT2llCont); histVec_1D->push_back(H_MT2lbCont);
    histVec_1D->push_back(H_PassMT2llCut80); histVec_1D->push_back(H_PassMT2llCut90); histVec_1D->push_back(H_PassMT2llCut100); histVec_1D->push_back(H_PassMT2llCut110); histVec_1D->push_back(H_PassMT2llCut120);
    histVec_1D->push_back(H_MT2ll_DPhiZMETClose); histVec_1D->push_back(H_MT2ll_DPhiZMETMid); histVec_1D->push_back(H_MT2ll_DPhiZMETFar);    
    histVec_1D->push_back(H_MT2lb);
    histVec_1D->push_back(H_MT2lb_DPhiBLep0BLep1Close); histVec_1D->push_back(H_MT2lb_DPhiBLep0BLep1Mid); histVec_1D->push_back(H_MT2lb_DPhiBLep0BLep1Far);
    histVec_1D->push_back(H_MET); histVec_1D->push_back(H_METX); histVec_1D->push_back(H_METY); histVec_1D->push_back(H_METPhi); histVec_1D->push_back(H_METPhi_noCorr); histVec_1D->push_back(H_METX_noPhiCorr); histVec_1D->push_back(H_METY_noPhiCorr);
    histVec_1D->push_back(H_NJets); histVec_1D->push_back(H_NJetswBTag);
    histVec_1D->push_back(H_DeltaPhiLep0Lep1); 
    histVec_1D->push_back(H_DeltaPhiZMET); histVec_1D->push_back(H_DeltaPhiZMET_PreCorr); 
    histVec_1D->push_back(H_nVtx); histVec_1D->push_back(H_nVtx_preRW);
    
    return histVec_1D;
}

inline vector<HistogramT> * OneDeeSmearHistTVec(bool fileContainsTrueSmearedJets) {
    const double PI = 3.14159265;
    int EnergyPtBinN = 40;
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
    
    HistogramT H_leadJetPt; H_leadJetPt.name = "h_SmearleadJetPt"; 
    H_leadJetPt.xLabel = "Lead Jet pT [GeV]"; H_leadJetPt.xBinN = EnergyPtBinN; H_leadJetPt.xMin = EnergyPtBinLB; H_leadJetPt.xMax = EnergyPtBinUB;  
    H_leadJetPt.yLabel = "Number of Events / ";
    numDivs = (H_leadJetPt.xMax - H_leadJetPt.xMin) / (float) H_leadJetPt.xBinN;
    H_leadJetPt.yLabel += "NUM"; H_leadJetPt.yLabel += " GeV";
    H_leadJetPt.xVarKey = "SmearleadJetPt";
    H_leadJetPt.doXSyst = true;
    
    HistogramT H_subJetPt; H_subJetPt.name = "h_SmearsubJetPt"; 
    H_subJetPt.xLabel = "Sub-lead Jet pT [GeV]"; H_subJetPt.xBinN = EnergyPtBinN; H_subJetPt.xMin = EnergyPtBinLB; H_subJetPt.xMax = EnergyPtBinUB;  
    H_subJetPt.yLabel = "Number of Events / ";
    numDivs = (H_subJetPt.xMax - H_subJetPt.xMin) / (float) H_subJetPt.xBinN;
    H_subJetPt.yLabel += "NUM"; H_subJetPt.yLabel += " GeV";
    H_subJetPt.xVarKey = "SmearsubJetPt";
    H_subJetPt.doXSyst = true;
    
    HistogramT H_leadBJetPt; H_leadBJetPt.name = "h_SmearleadBJetPt"; 
    H_leadBJetPt.xLabel = "Lead BJet pT [GeV]"; H_leadBJetPt.xBinN = EnergyPtBinN; H_leadBJetPt.xMin = EnergyPtBinLB; H_leadBJetPt.xMax = EnergyPtBinUB;  
    H_leadBJetPt.yLabel = "Number of Events / ";
    numDivs = (H_leadBJetPt.xMax - H_leadBJetPt.xMin) / (float) H_leadBJetPt.xBinN;
    H_leadBJetPt.yLabel += "NUM"; H_leadBJetPt.yLabel += " GeV";
    H_leadBJetPt.xVarKey = "SmearleadBJetPt";
    H_leadBJetPt.doXSyst = true;
    
    HistogramT H_leadBJetEn; H_leadBJetEn.name = "h_SmearleadBJetEn"; 
    H_leadBJetEn.xLabel = "Lead BJet Energy"; H_leadBJetEn.xBinN = 3*EnergyPtBinN; H_leadBJetEn.xMin = EnergyPtBinLB; H_leadBJetEn.xMax = 3*EnergyPtBinUB; 
    H_leadBJetEn.yLabel = "Number of Events / ";
    numDivs = (H_leadBJetEn.xMax - H_leadBJetEn.xMin) / (float) H_leadBJetEn.xBinN;
    H_leadBJetEn.yLabel += "NUM"; H_leadBJetEn.yLabel += " GeV";
    H_leadBJetEn.xVarKey = "SmearleadBJetEn";
    H_leadBJetEn.doXSyst = true;
    
    HistogramT H_subBJetPt; H_subBJetPt.name = "h_SmearsubBJetPt"; 
    H_subBJetPt.xLabel = "Sub-lead BJet pT [GeV]"; H_subBJetPt.xBinN = EnergyPtBinN; H_subBJetPt.xMin = EnergyPtBinLB; H_subBJetPt.xMax = EnergyPtBinUB;  
    H_subBJetPt.yLabel = "Number of Events / ";
    numDivs = (H_subBJetPt.xMax - H_subBJetPt.xMin) / (float) H_subBJetPt.xBinN;
    H_subBJetPt.yLabel += "NUM"; H_subBJetPt.yLabel += " GeV";
    H_subBJetPt.xVarKey = "SmearsubBJetPt";
    H_subBJetPt.doXSyst = true;
    
    HistogramT H_subBJetEn; H_subBJetEn.name = "h_SmearsubBJetEn"; 
    H_subBJetEn.xLabel = "Sub-lead BJet Energy"; H_subBJetEn.xBinN = EnergyPtBinN; H_subBJetEn.xMin = EnergyPtBinLB; H_subBJetEn.xMax = 3*EnergyPtBinUB; 
    H_subBJetEn.yLabel = "Number of Events / ";
    numDivs = (H_subBJetEn.xMax - H_subBJetEn.xMin) / (float) H_subBJetEn.xBinN;
    H_subBJetEn.yLabel += "NUM"; H_subBJetEn.yLabel += " GeV";
    H_subBJetEn.xVarKey = "SmearsubBJetEn";
    H_subBJetEn.doXSyst = true;
    
    ///Met or MET related variable plots
    
    HistogramT H_MT2ll; H_MT2ll.name = "h_SmearMT2ll"; 
    H_MT2ll.xLabel = "MT2_{ll} [GeV]"; H_MT2ll.xBinN = METBinN; H_MT2ll.xMin = METBinLB; H_MT2ll.xMax = METBinUB; 
    H_MT2ll.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll.xMax - H_MT2ll.xMin) / (float) H_MT2ll.xBinN;
    H_MT2ll.yLabel += "NUM"; H_MT2ll.yLabel += " GeV";
    H_MT2ll.xVarKey = "SmearMT2ll";
    H_MT2ll.doXSyst = true;
    
    HistogramT H_PassMT2llCut80; H_PassMT2llCut80.name = "h_SmearPassMT2llCut80"; 
    H_PassMT2llCut80.xLabel = "Event MT2_{ll} > 80 [GeV]"; H_PassMT2llCut80.xBinN = 2; H_PassMT2llCut80.xMin = -0.5; H_PassMT2llCut80.xMax = 1.5; 
    H_PassMT2llCut80.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut80.xVarKey = "SmearPassMT2llCut80";
    H_PassMT2llCut80.doXSyst = true;
    
    HistogramT H_PassMT2llCut90; H_PassMT2llCut90.name = "h_SmearPassMT2llCut90"; 
    H_PassMT2llCut90.xLabel = "Event MT2_{ll} > 90 [GeV]"; H_PassMT2llCut90.xBinN = 2; H_PassMT2llCut90.xMin = -0.5; H_PassMT2llCut90.xMax = 1.5; 
    H_PassMT2llCut90.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut90.xVarKey = "SmearPassMT2llCut90";
    H_PassMT2llCut90.doXSyst = true;
    
    HistogramT H_PassMT2llCut100; H_PassMT2llCut100.name = "h_SmearPassMT2llCut100"; 
    H_PassMT2llCut100.xLabel = "Event MT2_{ll} > 100 [GeV]"; H_PassMT2llCut100.xBinN = 2; H_PassMT2llCut100.xMin = -0.5; H_PassMT2llCut100.xMax = 1.5; 
    H_PassMT2llCut100.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut100.xVarKey = "SmearPassMT2llCut100";
    H_PassMT2llCut100.doXSyst = true;
    
    HistogramT H_PassMT2llCut110; H_PassMT2llCut110.name = "h_SmearPassMT2llCut110"; 
    H_PassMT2llCut110.xLabel = "Event MT2_{ll} > 110 [GeV]"; H_PassMT2llCut110.xBinN = 2; H_PassMT2llCut110.xMin = -0.5; H_PassMT2llCut110.xMax = 1.5; 
    H_PassMT2llCut110.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut110.xVarKey = "SmearPassMT2llCut110";
    H_PassMT2llCut110.doXSyst = true;
    
    HistogramT H_PassMT2llCut120; H_PassMT2llCut120.name = "h_SmearPassMT2llCut120"; 
    H_PassMT2llCut120.xLabel = "Event MT2_{ll} > 120 [GeV]"; H_PassMT2llCut120.xBinN = 2; H_PassMT2llCut120.xMin = -0.5; H_PassMT2llCut120.xMax = 1.5; 
    H_PassMT2llCut120.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut120.xVarKey = "SmearPassMT2llCut120";
    H_PassMT2llCut120.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiZMETClose; H_MT2ll_DPhiZMETClose.name = "h_SmearMT2ll_DPhiZMETClose"; 
    H_MT2ll_DPhiZMETClose.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETClose.xBinN = METBinN; H_MT2ll_DPhiZMETClose.xMin = METBinLB; H_MT2ll_DPhiZMETClose.xMax = METBinUB; 
    H_MT2ll_DPhiZMETClose.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETClose.xMax - H_MT2ll_DPhiZMETClose.xMin) / (float) H_MT2ll_DPhiZMETClose.xBinN;
    H_MT2ll_DPhiZMETClose.yLabel += "NUM"; H_MT2ll_DPhiZMETClose.yLabel += " GeV";
    H_MT2ll_DPhiZMETClose.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiZMETClose.doXSyst = true;
    
    
    HistogramT H_MT2ll_DPhiZMETMid; H_MT2ll_DPhiZMETMid.name = "h_SmearMT2ll_DPhiZMETMid"; 
    H_MT2ll_DPhiZMETMid.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETMid.xBinN = METBinN; H_MT2ll_DPhiZMETMid.xMin = METBinLB; H_MT2ll_DPhiZMETMid.xMax = METBinUB; 
    H_MT2ll_DPhiZMETMid.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETMid.xMax - H_MT2ll_DPhiZMETMid.xMin) / (float) H_MT2ll_DPhiZMETMid.xBinN;
    H_MT2ll_DPhiZMETMid.yLabel += "NUM"; H_MT2ll_DPhiZMETMid.yLabel += " GeV";
    H_MT2ll_DPhiZMETMid.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiZMETMid.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiZMETFar; H_MT2ll_DPhiZMETFar.name = "h_SmearMT2ll_DPhiZMETFar"; 
    H_MT2ll_DPhiZMETFar.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETFar.xBinN = METBinN; H_MT2ll_DPhiZMETFar.xMin = METBinLB; H_MT2ll_DPhiZMETFar.xMax = METBinUB; 
    H_MT2ll_DPhiZMETFar.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETFar.xMax - H_MT2ll_DPhiZMETFar.xMin) / (float) H_MT2ll_DPhiZMETFar.xBinN;
    H_MT2ll_DPhiZMETFar.yLabel += "NUM"; H_MT2ll_DPhiZMETFar.yLabel += " GeV";
    H_MT2ll_DPhiZMETFar.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiZMETFar.doXSyst = true;
    
    
    HistogramT H_MT2ll_DPhiLep0Lep1Close; H_MT2ll_DPhiLep0Lep1Close.name = "h_SmearMT2ll_DPhiLep0Lep1Close"; 
    H_MT2ll_DPhiLep0Lep1Close.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Close.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Close.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Close.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Close.xMax - H_MT2ll_DPhiLep0Lep1Close.xMin) / (float) H_MT2ll_DPhiLep0Lep1Close.xBinN;
    H_MT2ll_DPhiLep0Lep1Close.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Close.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Close.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiLep0Lep1Close.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Mid; H_MT2ll_DPhiLep0Lep1Mid.name = "h_SmearMT2ll_DPhiLep0Lep1Mid"; 
    H_MT2ll_DPhiLep0Lep1Mid.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Mid.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Mid.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Mid.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Mid.xMax - H_MT2ll_DPhiLep0Lep1Mid.xMin) / (float) H_MT2ll_DPhiLep0Lep1Mid.xBinN;
    H_MT2ll_DPhiLep0Lep1Mid.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Mid.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Mid.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiLep0Lep1Mid.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Far; H_MT2ll_DPhiLep0Lep1Far.name = "h_SmearMT2ll_DPhiLep0Lep1Far"; 
    H_MT2ll_DPhiLep0Lep1Far.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Far.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Far.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Far.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Far.xMax - H_MT2ll_DPhiLep0Lep1Far.xMin) / (float) H_MT2ll_DPhiLep0Lep1Far.xBinN;
    H_MT2ll_DPhiLep0Lep1Far.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Far.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Far.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiLep0Lep1Far.doXSyst = true;
    
    
    HistogramT H_MT2lb; H_MT2lb.name = "h_SmearMT2lb"; 
    H_MT2lb.xLabel = "MT2lb [GeV]"; H_MT2lb.xBinN = METBinN; H_MT2lb.xMin = METBinLB; H_MT2lb.xMax = METBinUB; 
    H_MT2lb.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb.xMax - H_MT2lb.xMin) / (float) H_MT2lb.xBinN;
    H_MT2lb.yLabel += "NUM"; H_MT2lb.yLabel += " GeV";
    H_MT2lb.xVarKey = "SmearMT2lb";
    H_MT2lb.doXSyst = true;
    HistogramT H_MT2lb_DPhiBLep0BLep1Close; H_MT2lb_DPhiBLep0BLep1Close.name = "h_SmearMT2lb_DPhiBLep0BLep1Close"; 
    H_MT2lb_DPhiBLep0BLep1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Close.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Close.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Close.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Close.xMax - H_MT2lb_DPhiBLep0BLep1Close.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Close.xBinN;
    H_MT2lb_DPhiBLep0BLep1Close.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Close.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Close.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiBLep0BLep1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Mid; H_MT2lb_DPhiBLep0BLep1Mid.name = "h_SmearMT2lb_DPhiBLep0BLep1Mid"; 
    H_MT2lb_DPhiBLep0BLep1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Mid.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Mid.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Mid.xMax - H_MT2lb_DPhiBLep0BLep1Mid.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Mid.xBinN;
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Mid.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Mid.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiBLep0BLep1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Far; H_MT2lb_DPhiBLep0BLep1Far.name = "h_SmearMT2lb_DPhiBLep0BLep1Far"; 
    H_MT2lb_DPhiBLep0BLep1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Far.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Far.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Far.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Far.xMax - H_MT2lb_DPhiBLep0BLep1Far.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Far.xBinN;
    H_MT2lb_DPhiBLep0BLep1Far.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Far.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Far.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiBLep0BLep1Far.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Close; H_MT2lb_DPhiJet0Jet1Close.name = "h_SmearMT2lb_DPhiJet0Jet1Close"; 
    H_MT2lb_DPhiJet0Jet1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Close.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Close.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Close.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Close.xMax - H_MT2lb_DPhiJet0Jet1Close.xMin) / (float) H_MT2lb_DPhiJet0Jet1Close.xBinN;
    H_MT2lb_DPhiJet0Jet1Close.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Close.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Close.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiJet0Jet1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Mid; H_MT2lb_DPhiJet0Jet1Mid.name = "h_SmearMT2lb_DPhiJet0Jet1Mid"; 
    H_MT2lb_DPhiJet0Jet1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Mid.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Mid.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Mid.xMax - H_MT2lb_DPhiJet0Jet1Mid.xMin) / (float) H_MT2lb_DPhiJet0Jet1Mid.xBinN;
    H_MT2lb_DPhiJet0Jet1Mid.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Mid.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Mid.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiJet0Jet1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Far; H_MT2lb_DPhiJet0Jet1Far.name = "h_SmearMT2lb_DPhiJet0Jet1Far"; 
    H_MT2lb_DPhiJet0Jet1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Far.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Far.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Far.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Far.xMax - H_MT2lb_DPhiJet0Jet1Far.xMin) / (float) H_MT2lb_DPhiJet0Jet1Far.xBinN;
    H_MT2lb_DPhiJet0Jet1Far.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Far.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Far.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiJet0Jet1Far.doXSyst = true;
    
    HistogramT H_MET; H_MET.name = "h_SmearMET"; 
    H_MET.xLabel = "#slash{E}_{T} [GeV]"; H_MET.xBinN = METBinN; H_MET.xMin = METBinLB; H_MET.xMax = METBinUB; 
    H_MET.yLabel = "Number of Events / ";
    numDivs = (H_MET.xMax - H_MET.xMin) / (float) H_MET.xBinN;
    H_MET.yLabel += "NUM"; H_MET.yLabel += " GeV";
    H_MET.xVarKey = "SmearMET";
    H_MET.doXSyst = true;
        
    HistogramT H_METPhi; H_METPhi.name = "h_SmearMETPhi"; 
    H_METPhi.xLabel = "#slash{E}_{T} #phi"; H_METPhi.xBinN = PhiBinN; H_METPhi.xMin = -PI; H_METPhi.xMax = PI; 
    H_METPhi.yLabel = "Number of Events / ";
    numDivs = (H_METPhi.xMax - H_METPhi.xMin) / (float) H_METPhi.xBinN;
    H_METPhi.yLabel += "NUM"; H_METPhi.yLabel += " radians";
    H_METPhi.xVarKey = "SmearMETPhi";
    H_METPhi.doXSyst = true;   
    
    HistogramT H_METX; H_METX.name = "h_SmearMETX";
    H_METX.xLabel = "#slash{E}_{x} [GeV]"; H_METX.xBinN = METXYBinN; H_METX.xMin = METXYBinLB; H_METX.xMax = METXYBinUB;
    H_METX.yLabel = "Number of Events / ";
    numDivs = (H_METX.xMax - H_METX.xMin) / (float) H_METX.xBinN;
    H_METX.yLabel += "NUM"; H_METX.yLabel += " GeV";
    H_METX.xVarKey = "SmearMETX";
    H_METX.doXSyst = true;
    
    HistogramT H_METY; H_METY.name = "h_SmearMETY";
    H_METY.xLabel = "#slash{E}_{y} [GeV]"; H_METY.xBinN = METXYBinN; H_METY.xMin = METXYBinLB; H_METY.xMax = METXYBinUB;
    H_METY.yLabel = "Number of Events / ";
    numDivs = (H_METY.xMax - H_METY.xMin) / (float) H_METY.xBinN;
    H_METY.yLabel += "NUM"; H_METY.yLabel += " GeV";
    H_METY.xVarKey = "SmearMETY";
    H_METY.doXSyst = true;
    
    ///Jet plots/////
    
    HistogramT H_NJets; H_NJets.name = "h_SmearNJets"; 
    H_NJets.xLabel = "N_{jets}"; H_NJets.xBinN = NJetsBinN; H_NJets.xMin = NJetsBinLB; H_NJets.xMax = NJetsBinUB; 
    H_NJets.yLabel = "Events / N_{jets}";
    H_NJets.xVarKey = "SmearNJets";
    H_NJets.doXSyst = true;
    
    HistogramT H_NJetswBTag; H_NJetswBTag.name = "h_SmearNJetswBTag"; 
    H_NJetswBTag.xLabel = "N_{b-jets}"; H_NJetswBTag.xBinN = NJetsBinN; H_NJetswBTag.xMin = NJetsBinLB; H_NJetswBTag.xMax = NJetsBinUB; 
    H_NJetswBTag.yLabel = "Events / N_{b-jets}";
    H_NJetswBTag.xVarKey = "SmearNBJets";
    H_NJetswBTag.doXSyst = true;
    
    //polarization energy "E-malgamation"
    HistogramT H_diLepEMinJetE; H_diLepEMinJetE.name = "h_SmeardiLepEMinJetE"; 
    H_diLepEMinJetE.xLabel = "E(l^{+}) + E(l^{-}) - E(J1) - E(J2) [GeV]"; H_diLepEMinJetE.xBinN = EnergyPtBinN; H_diLepEMinJetE.xMin = -EnergyPtBinUB; H_diLepEMinJetE.xMax = EnergyPtBinUB; 
    H_diLepEMinJetE.yLabel = "Number of Events / ";
    numDivs = (H_diLepEMinJetE.xMax - H_diLepEMinJetE.xMin) / (float) H_diLepEMinJetE.xBinN;
    H_diLepEMinJetE.yLabel += "NUM"; H_diLepEMinJetE.yLabel += " GeV";
    H_diLepEMinJetE.xVarKey = "SmearELepEJet";
    H_diLepEMinJetE.doXSyst = true;
    
    HistogramT H_MET_div_Meff; H_MET_div_Meff.name = "h_SmearMET_div_Meff";
    H_MET_div_Meff.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff.xBinN = 50; H_MET_div_Meff.xMin = 0.; H_MET_div_Meff.xMax = 1.;
    H_MET_div_Meff.yLabel = "Number of Events / "; H_MET_div_Meff.yLabel += "NUM";
    H_MET_div_Meff.xVarKey = "SmearMETdivMeff";
    H_MET_div_Meff.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut80; H_MET_div_Meff_PassMT2llCut80.name = "h_SmearMET_div_Meff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut80.xBinN = 50; H_MET_div_Meff_PassMT2llCut80.xMin = 0.; H_MET_div_Meff_PassMT2llCut80.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut80.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut80.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut80.xVarKey = "SmearMETdivMeff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut100; H_MET_div_Meff_PassMT2llCut100.name = "h_SmearMET_div_Meff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut100.xBinN = 50; H_MET_div_Meff_PassMT2llCut100.xMin = 0.; H_MET_div_Meff_PassMT2llCut100.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut100.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut100.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut100.xVarKey = "SmearMETdivMeff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut110; H_MET_div_Meff_PassMT2llCut110.name = "h_SmearMET_div_Meff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut110.xBinN = 50; H_MET_div_Meff_PassMT2llCut110.xMin = 0.; H_MET_div_Meff_PassMT2llCut110.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut110.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut110.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut110.xVarKey = "SmearMETdivMeff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.doXSyst = true;
    
    HistogramT H_HT; H_HT.name = "h_SmearHT"; 
    H_HT.xLabel = "H_{T} [GeV]"; H_HT.xBinN = EnergyPtBinN; H_HT.xMin = 0; H_HT.xMax = 3 * EnergyPtBinUB; 
    H_HT.yLabel = "Number of Events / ";
    numDivs = (H_HT.xMax - H_HT.xMin) / (float) H_HT.xBinN;
    H_HT.yLabel += "NUM"; H_HT.yLabel += " GeV";
    H_HT.xVarKey = "SmearHT";
    H_HT.doXSyst = true;
        
    //push the 1D histograms structures into a vector for eventual use in booking histograms
    vector<HistogramT> * histVecSmear_1D = new vector<HistogramT>;
    if (fileContainsTrueSmearedJets) {
        histVecSmear_1D->push_back(H_leadJetPt); histVecSmear_1D->push_back(H_subJetPt);
        histVecSmear_1D->push_back(H_leadBJetPt); histVecSmear_1D->push_back(H_leadBJetEn); histVecSmear_1D->push_back(H_subBJetPt); histVecSmear_1D->push_back(H_subBJetEn);
        histVecSmear_1D->push_back(H_NJets); histVecSmear_1D->push_back(H_NJetswBTag);
        histVecSmear_1D->push_back(H_diLepEMinJetE); histVecSmear_1D->push_back(H_HT);
    }
    histVecSmear_1D->push_back(H_MT2ll); 
    histVecSmear_1D->push_back(H_PassMT2llCut80); histVecSmear_1D->push_back(H_PassMT2llCut90); histVecSmear_1D->push_back(H_PassMT2llCut100); histVecSmear_1D->push_back(H_PassMT2llCut110); histVecSmear_1D->push_back(H_PassMT2llCut120);
    histVecSmear_1D->push_back(H_MT2ll_DPhiZMETClose); histVecSmear_1D->push_back(H_MT2ll_DPhiZMETMid); histVecSmear_1D->push_back(H_MT2ll_DPhiZMETFar);    
    histVecSmear_1D->push_back(H_MT2ll_DPhiLep0Lep1Close); histVecSmear_1D->push_back(H_MT2ll_DPhiLep0Lep1Mid); histVecSmear_1D->push_back(H_MT2ll_DPhiLep0Lep1Far);
    histVecSmear_1D->push_back(H_MET); histVecSmear_1D->push_back(H_METX); histVecSmear_1D->push_back(H_METY); histVecSmear_1D->push_back(H_METPhi);  
    histVecSmear_1D->push_back(H_MT2lb);            
    histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Close); histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Mid); histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Far);
    histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Close); histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Mid); histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Far);
    histVecSmear_1D->push_back(H_MET_div_Meff); 
    histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut80); histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut100); histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut110); 
    return histVecSmear_1D;
}


inline vector<HistogramT> * OneDeeHistTVecSpecial(bool doSmearedPlots, bool fileContainsTrueSmearedJets = 0) {
    const double PI = 3.14159265;
    int EnergyPtBinN = 40;
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
    ///Met or MET related variable plots
    
    HistogramT H_MT2ll; H_MT2ll.name = "h_SmearMT2ll"; 
    H_MT2ll.xLabel = "MT2_{ll} [GeV]"; H_MT2ll.xBinN = METBinN; H_MT2ll.xMin = METBinLB; H_MT2ll.xMax = METBinUB; 
    H_MT2ll.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll.xMax - H_MT2ll.xMin) / (float) H_MT2ll.xBinN;
    H_MT2ll.yLabel += "NUM"; H_MT2ll.yLabel += " GeV";
    H_MT2ll.xVarKey = "SmearMT2ll";
    H_MT2ll.doXSyst = true;
    
    HistogramT H_PassMT2llCut80; H_PassMT2llCut80.name = "h_SmearPassMT2llCut80"; 
    H_PassMT2llCut80.xLabel = "Event MT2_{ll} > 80 [GeV]"; H_PassMT2llCut80.xBinN = 2; H_PassMT2llCut80.xMin = -0.5; H_PassMT2llCut80.xMax = 1.5; 
    H_PassMT2llCut80.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut80.xVarKey = "SmearPassMT2llCut80";
    H_PassMT2llCut80.doXSyst = true;
    
    HistogramT H_PassMT2llCut90; H_PassMT2llCut90.name = "h_SmearPassMT2llCut90"; 
    H_PassMT2llCut90.xLabel = "Event MT2_{ll} > 90 [GeV]"; H_PassMT2llCut90.xBinN = 2; H_PassMT2llCut90.xMin = -0.5; H_PassMT2llCut90.xMax = 1.5; 
    H_PassMT2llCut90.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut90.xVarKey = "SmearPassMT2llCut90";
    H_PassMT2llCut90.doXSyst = true;
    
    HistogramT H_PassMT2llCut100; H_PassMT2llCut100.name = "h_SmearPassMT2llCut100"; 
    H_PassMT2llCut100.xLabel = "Event MT2_{ll} > 100 [GeV]"; H_PassMT2llCut100.xBinN = 2; H_PassMT2llCut100.xMin = -0.5; H_PassMT2llCut100.xMax = 1.5; 
    H_PassMT2llCut100.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut100.xVarKey = "SmearPassMT2llCut100";
    H_PassMT2llCut100.doXSyst = true;
    
    HistogramT H_PassMT2llCut110; H_PassMT2llCut110.name = "h_SmearPassMT2llCut110"; 
    H_PassMT2llCut110.xLabel = "Event MT2_{ll} > 110 [GeV]"; H_PassMT2llCut110.xBinN = 2; H_PassMT2llCut110.xMin = -0.5; H_PassMT2llCut110.xMax = 1.5; 
    H_PassMT2llCut110.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut110.xVarKey = "SmearPassMT2llCut110";
    H_PassMT2llCut110.doXSyst = true;
    
    HistogramT H_PassMT2llCut120; H_PassMT2llCut120.name = "h_SmearPassMT2llCut120"; 
    H_PassMT2llCut120.xLabel = "Event MT2_{ll} > 120 [GeV]"; H_PassMT2llCut120.xBinN = 2; H_PassMT2llCut120.xMin = -0.5; H_PassMT2llCut120.xMax = 1.5; 
    H_PassMT2llCut120.yLabel = "Events Passing/Failing MT2ll Cut";
    H_PassMT2llCut120.xVarKey = "SmearPassMT2llCut120";
    H_PassMT2llCut120.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiZMETClose; H_MT2ll_DPhiZMETClose.name = "h_SmearMT2ll_DPhiZMETClose"; 
    H_MT2ll_DPhiZMETClose.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETClose.xBinN = METBinN; H_MT2ll_DPhiZMETClose.xMin = METBinLB; H_MT2ll_DPhiZMETClose.xMax = METBinUB; 
    H_MT2ll_DPhiZMETClose.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETClose.xMax - H_MT2ll_DPhiZMETClose.xMin) / (float) H_MT2ll_DPhiZMETClose.xBinN;
    H_MT2ll_DPhiZMETClose.yLabel += "NUM"; H_MT2ll_DPhiZMETClose.yLabel += " GeV";
    H_MT2ll_DPhiZMETClose.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiZMETClose.doXSyst = true;
    
    
    HistogramT H_MT2ll_DPhiZMETMid; H_MT2ll_DPhiZMETMid.name = "h_SmearMT2ll_DPhiZMETMid"; 
    H_MT2ll_DPhiZMETMid.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETMid.xBinN = METBinN; H_MT2ll_DPhiZMETMid.xMin = METBinLB; H_MT2ll_DPhiZMETMid.xMax = METBinUB; 
    H_MT2ll_DPhiZMETMid.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETMid.xMax - H_MT2ll_DPhiZMETMid.xMin) / (float) H_MT2ll_DPhiZMETMid.xBinN;
    H_MT2ll_DPhiZMETMid.yLabel += "NUM"; H_MT2ll_DPhiZMETMid.yLabel += " GeV";
    H_MT2ll_DPhiZMETMid.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiZMETMid.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiZMETFar; H_MT2ll_DPhiZMETFar.name = "h_SmearMT2ll_DPhiZMETFar"; 
    H_MT2ll_DPhiZMETFar.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiZMETFar.xBinN = METBinN; H_MT2ll_DPhiZMETFar.xMin = METBinLB; H_MT2ll_DPhiZMETFar.xMax = METBinUB; 
    H_MT2ll_DPhiZMETFar.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiZMETFar.xMax - H_MT2ll_DPhiZMETFar.xMin) / (float) H_MT2ll_DPhiZMETFar.xBinN;
    H_MT2ll_DPhiZMETFar.yLabel += "NUM"; H_MT2ll_DPhiZMETFar.yLabel += " GeV";
    H_MT2ll_DPhiZMETFar.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiZMETFar.doXSyst = true;
    
    
    HistogramT H_MT2ll_DPhiLep0Lep1Close; H_MT2ll_DPhiLep0Lep1Close.name = "h_SmearMT2ll_DPhiLep0Lep1Close"; 
    H_MT2ll_DPhiLep0Lep1Close.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Close.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Close.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Close.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Close.xMax - H_MT2ll_DPhiLep0Lep1Close.xMin) / (float) H_MT2ll_DPhiLep0Lep1Close.xBinN;
    H_MT2ll_DPhiLep0Lep1Close.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Close.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Close.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiLep0Lep1Close.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Mid; H_MT2ll_DPhiLep0Lep1Mid.name = "h_SmearMT2ll_DPhiLep0Lep1Mid"; 
    H_MT2ll_DPhiLep0Lep1Mid.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Mid.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Mid.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Mid.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Mid.xMax - H_MT2ll_DPhiLep0Lep1Mid.xMin) / (float) H_MT2ll_DPhiLep0Lep1Mid.xBinN;
    H_MT2ll_DPhiLep0Lep1Mid.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Mid.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Mid.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiLep0Lep1Mid.doXSyst = true;
    
    HistogramT H_MT2ll_DPhiLep0Lep1Far; H_MT2ll_DPhiLep0Lep1Far.name = "h_SmearMT2ll_DPhiLep0Lep1Far"; 
    H_MT2ll_DPhiLep0Lep1Far.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_DPhiLep0Lep1Far.xBinN = METBinN; H_MT2ll_DPhiLep0Lep1Far.xMin = METBinLB; H_MT2ll_DPhiLep0Lep1Far.xMax = METBinUB; 
    H_MT2ll_DPhiLep0Lep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2ll_DPhiLep0Lep1Far.xMax - H_MT2ll_DPhiLep0Lep1Far.xMin) / (float) H_MT2ll_DPhiLep0Lep1Far.xBinN;
    H_MT2ll_DPhiLep0Lep1Far.yLabel += "NUM"; H_MT2ll_DPhiLep0Lep1Far.yLabel += " GeV";
    H_MT2ll_DPhiLep0Lep1Far.xVarKey = "SmearMT2ll";
    H_MT2ll_DPhiLep0Lep1Far.doXSyst = true;
    
    
    HistogramT H_MT2lb; H_MT2lb.name = "h_SmearMT2lb"; 
    H_MT2lb.xLabel = "MT2lb [GeV]"; H_MT2lb.xBinN = METBinN; H_MT2lb.xMin = METBinLB; H_MT2lb.xMax = METBinUB; 
    H_MT2lb.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb.xMax - H_MT2lb.xMin) / (float) H_MT2lb.xBinN;
    H_MT2lb.yLabel += "NUM"; H_MT2lb.yLabel += " GeV";
    H_MT2lb.xVarKey = "SmearMT2lb";
    H_MT2lb.doXSyst = true;
    HistogramT H_MT2lb_DPhiBLep0BLep1Close; H_MT2lb_DPhiBLep0BLep1Close.name = "h_SmearMT2lb_DPhiBLep0BLep1Close"; 
    H_MT2lb_DPhiBLep0BLep1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Close.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Close.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Close.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Close.xMax - H_MT2lb_DPhiBLep0BLep1Close.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Close.xBinN;
    H_MT2lb_DPhiBLep0BLep1Close.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Close.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Close.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiBLep0BLep1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Mid; H_MT2lb_DPhiBLep0BLep1Mid.name = "h_SmearMT2lb_DPhiBLep0BLep1Mid"; 
    H_MT2lb_DPhiBLep0BLep1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Mid.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Mid.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Mid.xMax - H_MT2lb_DPhiBLep0BLep1Mid.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Mid.xBinN;
    H_MT2lb_DPhiBLep0BLep1Mid.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Mid.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Mid.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiBLep0BLep1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiBLep0BLep1Far; H_MT2lb_DPhiBLep0BLep1Far.name = "h_SmearMT2lb_DPhiBLep0BLep1Far"; 
    H_MT2lb_DPhiBLep0BLep1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiBLep0BLep1Far.xBinN = METBinN; H_MT2lb_DPhiBLep0BLep1Far.xMin = METBinLB; H_MT2lb_DPhiBLep0BLep1Far.xMax = METBinUB; 
    H_MT2lb_DPhiBLep0BLep1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiBLep0BLep1Far.xMax - H_MT2lb_DPhiBLep0BLep1Far.xMin) / (float) H_MT2lb_DPhiBLep0BLep1Far.xBinN;
    H_MT2lb_DPhiBLep0BLep1Far.yLabel += "NUM"; H_MT2lb_DPhiBLep0BLep1Far.yLabel += " GeV";
    H_MT2lb_DPhiBLep0BLep1Far.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiBLep0BLep1Far.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Close; H_MT2lb_DPhiJet0Jet1Close.name = "h_SmearMT2lb_DPhiJet0Jet1Close"; 
    H_MT2lb_DPhiJet0Jet1Close.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Close.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Close.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Close.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Close.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Close.xMax - H_MT2lb_DPhiJet0Jet1Close.xMin) / (float) H_MT2lb_DPhiJet0Jet1Close.xBinN;
    H_MT2lb_DPhiJet0Jet1Close.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Close.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Close.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiJet0Jet1Close.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Mid; H_MT2lb_DPhiJet0Jet1Mid.name = "h_SmearMT2lb_DPhiJet0Jet1Mid"; 
    H_MT2lb_DPhiJet0Jet1Mid.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Mid.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Mid.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Mid.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Mid.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Mid.xMax - H_MT2lb_DPhiJet0Jet1Mid.xMin) / (float) H_MT2lb_DPhiJet0Jet1Mid.xBinN;
    H_MT2lb_DPhiJet0Jet1Mid.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Mid.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Mid.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiJet0Jet1Mid.doXSyst = true;
    
    HistogramT H_MT2lb_DPhiJet0Jet1Far; H_MT2lb_DPhiJet0Jet1Far.name = "h_SmearMT2lb_DPhiJet0Jet1Far"; 
    H_MT2lb_DPhiJet0Jet1Far.xLabel = "MT2lb [GeV] [GeV]"; H_MT2lb_DPhiJet0Jet1Far.xBinN = METBinN; H_MT2lb_DPhiJet0Jet1Far.xMin = METBinLB; H_MT2lb_DPhiJet0Jet1Far.xMax = METBinUB; 
    H_MT2lb_DPhiJet0Jet1Far.yLabel = "Number of Events / ";
    numDivs = (H_MT2lb_DPhiJet0Jet1Far.xMax - H_MT2lb_DPhiJet0Jet1Far.xMin) / (float) H_MT2lb_DPhiJet0Jet1Far.xBinN;
    H_MT2lb_DPhiJet0Jet1Far.yLabel += "NUM"; H_MT2lb_DPhiJet0Jet1Far.yLabel += " GeV";
    H_MT2lb_DPhiJet0Jet1Far.xVarKey = "SmearMT2lb";
    H_MT2lb_DPhiJet0Jet1Far.doXSyst = true;
        
    HistogramT H_MET_div_Meff_PassMT2llCut80; H_MET_div_Meff_PassMT2llCut80.name = "h_SmearMET_div_Meff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut80.xBinN = 50; H_MET_div_Meff_PassMT2llCut80.xMin = 0.; H_MET_div_Meff_PassMT2llCut80.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut80.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut80.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut80.xVarKey = "SmearMETdivMeff_PassMT2llCut80";
    H_MET_div_Meff_PassMT2llCut80.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut90; H_MET_div_Meff_PassMT2llCut90.name = "h_SmearMET_div_Meff_PassMT2llCut90";
    H_MET_div_Meff_PassMT2llCut90.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut90.xBinN = 50; H_MET_div_Meff_PassMT2llCut90.xMin = 0.; H_MET_div_Meff_PassMT2llCut90.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut90.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut90.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut90.xVarKey = "SmearMETdivMeff_PassMT2llCut90";
    H_MET_div_Meff_PassMT2llCut90.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut100; H_MET_div_Meff_PassMT2llCut100.name = "h_SmearMET_div_Meff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut100.xBinN = 50; H_MET_div_Meff_PassMT2llCut100.xMin = 0.; H_MET_div_Meff_PassMT2llCut100.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut100.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut100.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut100.xVarKey = "SmearMETdivMeff_PassMT2llCut100";
    H_MET_div_Meff_PassMT2llCut100.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut110; H_MET_div_Meff_PassMT2llCut110.name = "h_SmearMET_div_Meff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut110.xBinN = 50; H_MET_div_Meff_PassMT2llCut110.xMin = 0.; H_MET_div_Meff_PassMT2llCut110.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut110.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut110.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut110.xVarKey = "SmearMETdivMeff_PassMT2llCut110";
    H_MET_div_Meff_PassMT2llCut110.doXSyst = true;
    
    HistogramT H_MET_div_Meff_PassMT2llCut120; H_MET_div_Meff_PassMT2llCut120.name = "h_SmearMET_div_Meff_PassMT2llCut120";
    H_MET_div_Meff_PassMT2llCut120.xLabel = "#slash{E}_{T} / M_{eff}"; H_MET_div_Meff_PassMT2llCut120.xBinN = 50; H_MET_div_Meff_PassMT2llCut120.xMin = 0.; H_MET_div_Meff_PassMT2llCut120.xMax = 1.;
    H_MET_div_Meff_PassMT2llCut120.yLabel = "Number of Events / "; H_MET_div_Meff_PassMT2llCut120.yLabel += "NUM";
    H_MET_div_Meff_PassMT2llCut120.xVarKey = "SmearMETdivMeff_PassMT2llCut120";
    H_MET_div_Meff_PassMT2llCut120.doXSyst = true;
        
    //push the 1D histograms structures into a vector for eventual use in booking histograms
    vector<HistogramT> * histVecSmear_1D = new vector<HistogramT>;
    if (fileContainsTrueSmearedJets) {
        histVecSmear_1D->push_back(H_MT2lb);            
        histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Close); histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Mid); histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Far);
        histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Close); histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Mid); histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Far);
        histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut80); histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut90); histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut100); histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut110); histVecSmear_1D->push_back(H_MET_div_Meff_PassMT2llCut120);
    }
    else if (!doSmearedPlots) {
        histVecSmear_1D->push_back(H_MT2lb);            
        histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Close); histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Mid); histVecSmear_1D->push_back(H_MT2lb_DPhiBLep0BLep1Far);
        histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Close); histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Mid); histVecSmear_1D->push_back(H_MT2lb_DPhiJet0Jet1Far);
    }       
    histVecSmear_1D->push_back(H_MT2ll); 
    histVecSmear_1D->push_back(H_PassMT2llCut80); histVecSmear_1D->push_back(H_PassMT2llCut90); histVecSmear_1D->push_back(H_PassMT2llCut100); histVecSmear_1D->push_back(H_PassMT2llCut110); histVecSmear_1D->push_back(H_PassMT2llCut120);
    histVecSmear_1D->push_back(H_MT2ll_DPhiZMETClose); histVecSmear_1D->push_back(H_MT2ll_DPhiZMETMid); histVecSmear_1D->push_back(H_MT2ll_DPhiZMETFar);    
    histVecSmear_1D->push_back(H_MT2ll_DPhiLep0Lep1Close); histVecSmear_1D->push_back(H_MT2ll_DPhiLep0Lep1Mid); histVecSmear_1D->push_back(H_MT2ll_DPhiLep0Lep1Far);
    if (!doSmearedPlots) {
        for (unsigned int iHistT = 0; iHistT < histVecSmear_1D->size(); ++iHistT) {
            histVecSmear_1D->at(iHistT).name.Replace(2, 5, "");
        }
    }
    return histVecSmear_1D;
}


inline vector<HistogramT> * TwoDeeSmearHistTVec() {
    
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
    
    HistogramT H_MT2ll_vs_MT2lb; H_MT2ll_vs_MT2lb.name = "h_SmearMT2ll_vs_MT2lb";
    H_MT2ll_vs_MT2lb.xLabel = "MT2_{ll} [GeV]"; H_MT2ll_vs_MT2lb.xBinN = METBinN; H_MT2ll_vs_MT2lb.xMin = METBinLB; H_MT2ll_vs_MT2lb.xMax = METBinUB;
    H_MT2ll_vs_MT2lb.yLabel = "MT2lb [GeV]"; H_MT2ll_vs_MT2lb.yBinN = METBinN; H_MT2ll_vs_MT2lb.yMin = METBinLB; H_MT2ll_vs_MT2lb.yMax = METBinUB;
    H_MT2ll_vs_MT2lb.xVarKey = "SmearMT2ll";
    H_MT2ll_vs_MT2lb.yVarKey = "SmearMT2lb";
    H_MT2ll_vs_MT2lb.doXSyst = true;
    H_MT2ll_vs_MT2lb.doYSyst = true;
    
    histVec_2D->push_back(H_MT2ll_vs_MT2lb);
    return histVec_2D;
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
    events_LepInZMassJet2BJet1[2].blindDataChannel = 1; // blind emu without ZVeto + jet cuts
    events_LepInZMassJet2BJet1MET40[2].blindDataChannel = 1; // blind with MET + ZVeto + jet cuts            

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
    /*
    SampleT events_NBJetsGeq1; events_NBJetsGeq1.histNameSuffix = "_NBJetsGeq1";
    events_NBJetsGeq1.histXaxisSuffix = ""; events_NBJetsGeq1.histYaxisSuffix = "";  events_NBJetsGeq1.histZaxisSuffix = "";
    events_NBJetsGeq1.whichdiLepType = -1; events_NBJetsGeq1.doZVeto = -1; events_NBJetsGeq1.cutNJets = -1; events_NBJetsGeq1.cutNBJets = 1; events_NBJetsGeq1.cutMET = 0.;
    events_NBJetsGeq1.blindDataChannel = 0;
    */
    /*
    SampleT events_NBJetsGeq2; events_NBJetsGeq2.histNameSuffix = "_NBJetsGeq2";
    events_NBJetsGeq2.histXaxisSuffix = ""; events_NBJetsGeq2.histYaxisSuffix = "";  events_NBJetsGeq2.histZaxisSuffix = "";
    events_NBJetsGeq2.whichdiLepType = -1; events_NBJetsGeq2.doZVeto = -1; events_NBJetsGeq2.cutNJets = -1; events_NBJetsGeq2.cutNBJets = 2; events_NBJetsGeq2.cutMET = 0.;
    events_NBJetsGeq2.blindDataChannel = 0;
    */
    
    SampleT events_NJetsGeq2_NBJetsGeq1; events_NJetsGeq2_NBJetsGeq1.histNameSuffix = "_Jet2BJet1"; events_NJetsGeq2_NBJetsGeq1.histXaxisSuffix = ""; events_NJetsGeq2_NBJetsGeq1.histYaxisSuffix = ""; events_NJetsGeq2_NBJetsGeq1.histZaxisSuffix = "";
    events_NJetsGeq2_NBJetsGeq1.whichdiLepType = -1; events_NJetsGeq2_NBJetsGeq1.doZVeto = -1; events_NJetsGeq2_NBJetsGeq1.cutNJets = 2; events_NJetsGeq2_NBJetsGeq1.cutNBJets = 1; events_NJetsGeq2_NBJetsGeq1.cutMET = 0.;
    events_NJetsGeq2_NBJetsGeq1.blindDataChannel = 0;
    
    SampleT events_FullCut; events_FullCut.histNameSuffix = "_FullCut"; events_FullCut.histXaxisSuffix = ""; 
    events_FullCut.histYaxisSuffix = ""; events_FullCut.histZaxisSuffix = "";
    events_FullCut.whichdiLepType = -1; events_FullCut.doZVeto = 1; events_FullCut.cutNJets = 2; events_FullCut.cutNBJets = 1; events_FullCut.cutMET = 40.;
    events_FullCut.blindDataChannel = 1;
    
    SampleT events_FullCutBlind; events_FullCutBlind.histNameSuffix = "_FullCutBlind"; events_FullCutBlind.histXaxisSuffix = ""; 
    events_FullCutBlind.histYaxisSuffix = ""; events_FullCutBlind.histZaxisSuffix = "";
    events_FullCutBlind.whichdiLepType = -1; events_FullCutBlind.doZVeto = 1; events_FullCutBlind.cutNJets = 2; events_FullCutBlind.cutNBJets = 1; events_FullCutBlind.cutMET = 40.;
    events_FullCutBlind.blindDataChannel = 1; 
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
    
    SampleT events_inZMass0Jets; events_inZMass0Jets.histNameSuffix = "_allLep_inZ_0Jets";
    events_inZMass0Jets.histXaxisSuffix = ""; events_inZMass0Jets.histYaxisSuffix = ""; events_inZMass0Jets.histZaxisSuffix = "";
    events_inZMass0Jets.whichdiLepType = -1; events_inZMass0Jets.doZVeto = 0; events_inZMass0Jets.cutNJets = -1; events_inZMass0Jets.cutNBJets = 0; events_inZMass0Jets.cutMET = 0.;
    events_inZMass0Jets.blindDataChannel = 0;
    
    SampleT events_inZMass1Jet; events_inZMass1Jet.histNameSuffix = "_allLep_inZ_1Jet";
    events_inZMass1Jet.histXaxisSuffix = ""; events_inZMass1Jet.histYaxisSuffix = ""; events_inZMass1Jet.histZaxisSuffix = "";
    events_inZMass1Jet.whichdiLepType = -1; events_inZMass1Jet.doZVeto = 0; events_inZMass1Jet.cutNJets = -1; events_inZMass1Jet.cutNBJets = 0; events_inZMass1Jet.cutMET = 0.;
    events_inZMass1Jet.blindDataChannel = 0;
    
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
    subSampVec->push_back(events_NJetsGeq2);// subSampVec->push_back(events_NBJetsGeq1); 
//    subSampVec->push_back(events_NBJetsGeq2); 
    subSampVec->push_back(events_NJetsGeq2_NBJetsGeq1);
    subSampVec->push_back(events_FullCut); subSampVec->push_back(events_FullCutBlind);
    subSampVec->push_back(allEvents);
    subSampVec->push_back(events_inZMass); subSampVec->push_back(events_inZMassJet2BJet1);
    subSampVec->push_back(events_inZMassMET40); subSampVec->push_back(events_inZMassJet2BJet1MET40);
    subSampVec->push_back(events_inZMass0Jets); subSampVec->push_back(events_inZMass1Jet);
    return subSampVec;
}

inline vector<SystT> * SystVec(bool doJetSmearSyst) {
    // int whichSystType is used to demarcate between different systematics, with positive values designating a "shift up" systematic and negative values designating a "shift down" systematic -- note that systematics with no "shift" to them, i.e. GenTopRW which is unidirectional, are marked with a positive shift type by default
    // whichSystType number meanings:
    // 1: Lepton Energy Scale
    // 2: Jet Energy Scale
    // 3: BTag SF Shift
    // 4: Lepton Efficiency SF Shift
    // 5: Unclustered Energy Scale
    // 6: Stop XSec Shift
    // 7: Gen TopPt RW
    // 8: MT2ll Smear Shift
    
    SystT LepESShiftUp; LepESShiftUp.name = "_LepESShiftUp"; LepESShiftUp.systVarKey = "_LepESShiftUp";
    LepESShiftUp.whichSystType = 1;
    SystT LepESShiftDown; LepESShiftDown.name = "_LepESShiftDown"; LepESShiftDown.systVarKey = "_LepESShiftDown";
    LepESShiftDown.whichSystType = -1;
    
    SystT JetESShiftUp; JetESShiftUp.name = "_JetESShiftUp"; JetESShiftUp.systVarKey = "_JetESShiftUp";
    JetESShiftUp.whichSystType = 2;
    SystT JetESShiftDown; JetESShiftDown.name = "_JetESShiftDown"; JetESShiftDown.systVarKey = "_JetESShiftDown";
    JetESShiftDown.whichSystType = -2;
    
    SystT BTagSFShiftUp; BTagSFShiftUp.name = "_BTagSFShiftUp"; BTagSFShiftUp.systVarKey = "_BTagSFShiftUp";
    BTagSFShiftUp.whichSystType = 3;
    SystT BTagSFShiftDown; BTagSFShiftDown.name = "_BTagSFShiftDown"; BTagSFShiftDown.systVarKey = "_BTagSFShiftDown";
    BTagSFShiftDown.whichSystType = -3;
    
    SystT LepEffSFShiftUp; LepEffSFShiftUp.name = "_LepEffSFShiftUp"; LepEffSFShiftUp.systVarKey = "";
    LepEffSFShiftUp.whichSystType = 4;
    SystT LepEffSFShiftDown; LepEffSFShiftDown.name = "_LepEffSFShiftDown"; LepEffSFShiftDown.systVarKey = "";
    LepEffSFShiftDown.whichSystType = -4;
    
    SystT MT2UncEnShiftUp; MT2UncEnShiftUp.name = "_MT2UncESSpecShiftUp"; MT2UncEnShiftUp.systVarKey = "_MT2UncESSpecShiftUp";
    MT2UncEnShiftUp.whichSystType = 5;
    SystT MT2UncEnShiftDown; MT2UncEnShiftDown.name = "_MT2UncESSpecShiftDown"; MT2UncEnShiftDown.systVarKey = "_MT2UncESSpecShiftDown";
    MT2UncEnShiftDown.whichSystType = -5;
    
    SystT genStopXSecShiftUp; genStopXSecShiftUp.name = "_genStopXSecShiftUp"; genStopXSecShiftUp.systVarKey = "";
    genStopXSecShiftUp.whichSystType = 6;
    
    SystT genStopXSecShiftDown; genStopXSecShiftDown.name = "_genStopXSecShiftDown"; genStopXSecShiftDown.systVarKey = "";
    genStopXSecShiftDown.whichSystType =-6;
    
    SystT genTopReweight; genTopReweight.name = "_genTopRW"; genTopReweight.systVarKey = "";
    genTopReweight.whichSystType = 7;
    
    SystT MT2llShiftUp; MT2llShiftUp.name = "_MT2llShiftUp"; MT2llShiftUp.systVarKey = "_MT2llShiftUp";
    MT2llShiftUp.whichSystType = 8;
    
    SystT UncESShiftUp; UncESShiftUp.name = "_UncESShiftUp"; UncESShiftUp.systVarKey = "_UncESShiftUp";
    UncESShiftUp.whichSystType = 9;   
    SystT UncESShiftDown; UncESShiftDown.name = "_UncESShiftDown"; UncESShiftDown.systVarKey = "_UncESShiftDown";
    UncESShiftDown.whichSystType = -9;
    
    SystT JetSmearShiftUp; JetSmearShiftUp.name = "_JetSmearShiftUp"; JetSmearShiftUp.systVarKey = "_JetSmearShiftUp";
    JetSmearShiftUp.whichSystType = 10;   
    SystT JetSmearShiftDown; JetSmearShiftDown.name = "_JetSmearShiftDown"; JetSmearShiftDown.systVarKey = "_JetSmearShiftDown";
    JetSmearShiftDown.whichSystType = -10;
    
    
    vector<SystT> * systVec = new vector<SystT>;
    systVec->push_back(MT2llShiftUp);// systVec->push_back(MT2llShiftDown);
    systVec->push_back(LepEffSFShiftUp); systVec->push_back(LepEffSFShiftDown);
    systVec->push_back(BTagSFShiftUp); systVec->push_back(BTagSFShiftDown);
    systVec->push_back(LepESShiftUp); systVec->push_back(LepESShiftDown);
    systVec->push_back(JetESShiftUp); systVec->push_back(JetESShiftDown);
    systVec->push_back(MT2UncEnShiftUp), systVec->push_back(MT2UncEnShiftDown);
    systVec->push_back(genTopReweight);
    systVec->push_back(genStopXSecShiftUp); systVec->push_back(genStopXSecShiftDown);
    if (doJetSmearSyst) {
        systVec->push_back(JetSmearShiftUp); systVec->push_back(JetSmearShiftDown);
    }
    systVec->push_back(UncESShiftUp); systVec->push_back(UncESShiftDown);
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
//    TAxis * XAxis;
    TAxis * YAxis;
//    int XAxisLB, XAxisUB, 
    int YAxisLB, YAxisUB;
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
//    TAxis * XAxis;
    TAxis * YAxis, * ZAxis;
//    int XAxisLB, XAxisUB;
    int YAxisLB, YAxisUB, ZAxisLB, ZAxisUB;
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



inline vector<TF1> * VecJetHighPtResolutionTF1() {
    // Formula and parameter values for "sigma" taken from http://cmslxr.fnal.gov/lxr/source/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt
    // sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))
    /*
     EtaMin EtaMax  NMinuitParams   FuncMin FuncMax Par0        Par1     Par2   Par3
     0      0.5     6               0.      99999.  -0.349206   0.297831 0      0.471121
     0.5    1       6               0.      99999.  -0.499735   0.336391 0      0.430689
     1      1.5     6               0.      99999.  -0.561649   0.420293 0      0.392398
     1.5    2       6               0.      99999.  -1.12329    0.657891 0      0.139595
     2      2.5     6               0.      99999.  1.04792     0.466763 0      0.193137
     2.5    3       6               0.      99999.  2.56933     0.305802 0      0.398929
     3      3.5     6               0.      99999.  2.81978     0.272373 0      0.579396
     3.5    4       6               0.      99999.  1.65966     0.223683 0      0.60873
     4      9.9     6               0.      99999.  1.41584     0.209477 0      0.588872
     */
    
    TString baseTF1Name = "SigmaResolTF1", TF1Name;
    TString function = "sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))";
    vector<TF1> * outTF1Vec = new vector<TF1>;
    const int numTF1s = 9;
    float Par0[numTF1s] = {-0.349206, -0.499735, -0.561649, -1.12329, 1.04792, 2.56933, 2.81978, 1.65966, 1.41584};
    float Par1[numTF1s] = {0.297831, 0.336391, 0.420293, 0.657891, 0.466763, 0.305802, 0.272373, 0.223683, 0.209477};
    float Par2[numTF1s] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    float Par3[numTF1s] = {0.471121, 0.430689, 0.392398, 0.139595, 0.193137, 0.398929, 0.579396, 0.60873, 0.588872};
    for (int iFunc = 0; iFunc < numTF1s; ++iFunc) {
        TF1Name = baseTF1Name;
        TF1Name += "_EtaBin"; TF1Name += iFunc;
        TF1 currTF1(TF1Name, function, 0., 99999.);
        currTF1.SetParameter(0, Par0[iFunc]);
        currTF1.SetParameter(1, Par1[iFunc]);
        currTF1.SetParameter(2, Par2[iFunc]);
        currTF1.SetParameter(3, Par3[iFunc]);
        outTF1Vec->push_back(currTF1);
    }
    return outTF1Vec;
}
inline TH2F * SmearHist(TString inputFileName, std::map<string, int> * inputHistSmearMap, vector<TH2F *> * inputVecMETSmearHists) {
    TH2F * outHist;
    for(map<string, int>::iterator xIter = inputHistSmearMap->begin(); xIter != inputHistSmearMap->end(); xIter++) {
//        cout << "inputFileName " << inputFileName << endl;
//        cout << "xIter->first " << xIter->first << endl;
        if (inputFileName.Contains(xIter->first)) {
            if ((unsigned int) xIter->second > inputVecMETSmearHists->size() - 1) {
                outHist = new TH2F;
                break;
            }
            cout << inputVecMETSmearHists->at(xIter->second)->GetName() << endl;
            outHist = inputVecMETSmearHists->at(xIter->second);
            break;
        }
        else {
            outHist = new TH2F;
        }
        // iterator->first = key
        // iterator->second = value
        // Repeat if you also want to iterate through the second map.
    }
    return outHist;
}
inline void SetMETSmearHistMap(std::map<string, int> &inputHistSmearMap) {
    int TTBarIndex  = 0;
    int DYIndex     = 1;
    int WWIndex     = 2;
    int WZIndex     = 3;
    int ZZIndex     = 4;
    int WJIndex     = 5;
    int ZGIndex     = 6;
    int WGIndex     = 7;
    TH2F currHist;
    //Signal -- use TTBar for now...
    inputHistSmearMap["T2tt"] = TTBarIndex;
    inputHistSmearMap["T2bw"] = TTBarIndex;
    
    inputHistSmearMap["TTbar_Powheg"]           = TTBarIndex;
    inputHistSmearMap["TTbar_Madgraph"]         = TTBarIndex;
    inputHistSmearMap["TTbar_MCatNLO"]          = TTBarIndex;
    inputHistSmearMap["TTJetsFullHadrMG"]       = TTBarIndex;
    inputHistSmearMap["TTJetsSemiLeptMG"]       = TTBarIndex;
    inputHistSmearMap["TTJetsFullLeptMGtauola"] = TTBarIndex;
    inputHistSmearMap["TTJetsFullLeptMG"]       = TTBarIndex;
    inputHistSmearMap["TWDilep"]                = TTBarIndex;
    inputHistSmearMap["TbarWDilep"]             = TTBarIndex;
    inputHistSmearMap["TTWJets"]                = TTBarIndex;
    inputHistSmearMap["TTZJets"]                = TTBarIndex;
    inputHistSmearMap["TTWWJets"]               = TTBarIndex;
    inputHistSmearMap["TTGJets"]                = TTBarIndex;
    
    inputHistSmearMap["DYJets_Madgraph"]    = DYIndex;
    inputHistSmearMap["DY50toLLMad1jet"]    = DYIndex;
    inputHistSmearMap["DY50toLLMad2jet"]    = DYIndex;
    inputHistSmearMap["DY50toLLMad3jet"]    = DYIndex;
    inputHistSmearMap["DY50toLLMad4jet"]    = DYIndex;
    inputHistSmearMap["ZJets_Madgraph"]     = DYIndex;
        
    inputHistSmearMap["SingleTop_Tchannel"]     = WJIndex;
    inputHistSmearMap["SingleTopBar_Tchannel"]  = WJIndex;
    inputHistSmearMap["SingleTop_Schannel"]     = WJIndex;
    inputHistSmearMap["SingleTopBar_Schannel"]  = WJIndex;
    inputHistSmearMap["WJets_Madgraph"]         = WJIndex;
    
    inputHistSmearMap["WGstarToElNuMad"]    = WGIndex;
    inputHistSmearMap["WGstarToMuNuMad"]    = WGIndex;
    inputHistSmearMap["WGstarToTauNuMad"]   = WGIndex;
    inputHistSmearMap["WgammaToLNuG"]       = WGIndex;        

    inputHistSmearMap["ZgammaToLLG"]        = ZGIndex;
    
    inputHistSmearMap["ZZ_"]                = ZZIndex;    
    inputHistSmearMap["ZZZJets"]            = ZZIndex;
    inputHistSmearMap["HZZ4L"]              = ZZIndex;    
    
    inputHistSmearMap["WZZJets"]            = WZIndex;
    inputHistSmearMap["WZ_"]                = WZIndex;

    inputHistSmearMap["WWincl"]             = WWIndex;
    inputHistSmearMap["ggWWto2L"]           = WWIndex;
    inputHistSmearMap["WWTo2L2Nu_Madgraph"] = WWIndex;
    inputHistSmearMap["HWW125"]             = WWIndex;
    inputHistSmearMap["VBF125"]             = WWIndex;
    inputHistSmearMap["WWGJets"]            = WWIndex;
    inputHistSmearMap["WWZJets"]            = WWIndex;
    inputHistSmearMap["WWWJets"]            = WWIndex;                                                                   
}
inline vector<TH2F *> * SetMETSmearHistVec(int METorMETPhi, vector<TFile *> * inputVecMETSmearFiles, int whichSystCase = 0, int whichSystDirection = 0) {    
    TString METSmearString = "METSmearMinMET_vs_MET";
    TString METPhiSmearString = "METPhiSmearMinMETPhi_vs_MET";
    TString METSmearSystString[4] = {"METSmearJetERUpMinMET_vs_MET", "METSmearJetERDownMinMET_vs_MET", "METSmearUncESUpMinMET_vs_MET", "METSmearUncESDownMinMET_vs_MET"};
    TString METPhiSmearSystString[4] = {"METPhiSmearJetERUpMinMETPhi_vs_MET", "METPhiSmearJetERDownMinMETPhi_vs_MET", "METPhiSmearUncESUpMinMETPhi_vs_MET", "METPhiSmearUncESDownMinMETPhi_vs_MET"};
    TString TTBarHistString = "h_";
    int systGrabIndex = 0;
    if (whichSystCase > 0) {
        if (whichSystDirection == 0) cout << "error with Systs: they're turned on, but no direction has been specified " << whichSystDirection << endl;
        if (whichSystCase == 1) {
            systGrabIndex = whichSystDirection > 0 ? 0 : 1;
        }
        else if (whichSystCase == 2) {
            systGrabIndex = whichSystDirection > 0 ? 2 : 3;
        }
        TTBarHistString += METorMETPhi > 0 ? METSmearSystString[systGrabIndex] : METPhiSmearSystString[systGrabIndex];
    }
    else if (whichSystCase == 0) {
        TTBarHistString += METorMETPhi > 0 ? METSmearString : METPhiSmearString;
    }
    else {
        cout << "bug in whichSystCase! " << whichSystCase << endl;
    }
    TTBarHistString += "_inclusive";
    vector<TH2F *> * outVec = new vector<TH2F *>;
    TString OtherHistString = METorMETPhi > 0 ? "h_SmearMetMinUnsmearMet_vs_UnsmearMet" : "h_SmearMetPhiMinUnsmearMetPhi_vs_UnsmearMet";
    TString grabString;
    TH2F * currHist;
    for (unsigned int iFile = 0; iFile < inputVecMETSmearFiles->size(); ++iFile) {
        if (whichSystCase != 0 && iFile > 4) continue; // no JER or UncEnSmearing for nonTTBar for now
        grabString = iFile < 5 ? TTBarHistString : OtherHistString;
/*
        if (iFile == 1) {
            if (whichSystCase != 0) continue;
            grabString = OtherHistString;
        }
        */
        //        cout << "fileName " << inputVecMETSmearFiles->at(iFile)->GetName();
        currHist = (TH2F *) inputVecMETSmearFiles->at(iFile)->Get(grabString);
        outVec->push_back(currHist);

    }
    return outVec;
}

inline float DeltaMET(std::vector<TH1F *> * vecOneDeeHists, TH2F * TwoDeeHist, float inputMETValue, bool isMETPhi, float SF = 1.0) {
    //returns float which is a random number drawn from the distribution of MET (Phi) Smeared - MET (Phi) Unsmeared
    // Also works for the "Jet Smear systematics version"
    // the "isMETPhi" boolean turns on an additional random number that is +/- 1 (since dPhi for the METPhi p.d.f. is a magnitude thing)
    double sign = 1.0;
    
    if (isMETPhi) {
        TRandom3 rand(0);
        if (rand.Integer(2)) sign = -1;
//        cout << "sign " << sign << endl;
    }
    int whichOneDeeHist = TwoDeeHist->GetXaxis()->FindBin(inputMETValue) - 1;
    unsigned int numOneDeeHists = vecOneDeeHists->size();
    if (whichOneDeeHist < 0) {
        std::cout << inputMETValue << std::endl;
        std::cout << "ERROR with which one dee hist!" << std::endl;
        return 0.;
    }
    else if (whichOneDeeHist >= (int) numOneDeeHists) {
        whichOneDeeHist = (int) numOneDeeHists - 1;
    }
    if (vecOneDeeHists->at(whichOneDeeHist)->Integral() == 0) {
        return 0;
    }
    else {
        return SF * sign * vecOneDeeHists->at(whichOneDeeHist)->GetRandom();   
    }
}





inline vector<TH2F *> * SetMETXYSmearHistVec(int METXorMETY, vector<TFile *> * inputVecMETSmearFiles, int whichSystCase = 0, int whichSystDirection = 0) {
    TString METXSmearString = "METXSmearMinMETX_vs_MET";
    TString METYSmearString = "METYSmearMinMETY_vs_MET";
    TString METXSmearSystString[4] = {"METXSmearJetERUpMinMETX_vs_MET", "METXSmearJetERDownMinMETX_vs_MET", "METXSmearUncESUpMinMETX_vs_MET", "METXSmearUncESDownMinMETX_vs_MET"};
    TString METYSmearSystString[4] = {"METYSmearJetERUpMinMETY_vs_MET", "METYSmearJetERDownMinMETY_vs_MET", "METYSmearUncESUpMinMETY_vs_MET", "METYSmearUncESDownMinMETY_vs_MET"};
    TString TTBarHistString = "h_";
    int systGrabIndex = 0;
    if (whichSystCase > 0) {
        if (whichSystDirection == 0) cout << "error with Systs: they're turned on, but no direction has been specified " << whichSystDirection << endl;
        if (whichSystCase == 1) {
            systGrabIndex = whichSystDirection > 0 ? 0 : 1;
        }
        else if (whichSystCase == 2) {
            systGrabIndex = whichSystDirection > 0 ? 2 : 3;
        }
        TTBarHistString += METXorMETY > 0 ? METXSmearSystString[systGrabIndex] : METYSmearSystString[systGrabIndex];
    }
    else if (whichSystCase == 0) {
        TTBarHistString += METXorMETY > 0 ? METXSmearString : METYSmearString;
    }
    else {
        cout << "bug in whichSystCase! " << whichSystCase << endl;
    }
    TTBarHistString += "_inclusive";
    vector<TH2F *> * outVec = new vector<TH2F *>;
    TString OtherHistString = METXorMETY > 0 ? "h_SmearMetXMinUnsmearMetX_vs_UnsmearMet" : "h_SmearMetYMinUnsmearMetY_vs_UnsmearMet";
    TString grabString;
    TH2F * currHist;
    for (unsigned int iFile = 0; iFile < inputVecMETSmearFiles->size(); ++iFile) {
        if (whichSystCase != 0 && iFile > 4) continue; // no JER or UncEnSmearing for nonTTBar for now
        grabString = iFile < 5 ? TTBarHistString : OtherHistString;
        //        cout << "fileName " << inputVecMETSmearFiles->at(iFile)->GetName();
        currHist = (TH2F *) inputVecMETSmearFiles->at(iFile)->Get(grabString);
        outVec->push_back(currHist);        
    }
    return outVec;
}

inline float DeltaMETXY(std::vector<TH1F *> * vecOneDeeHists, TH2F * TwoDeeHist, float inputMETValue) {
    //returns float which is a random number drawn from the distribution of MET (Phi) Smeared - MET (Phi) Unsmeared
    // Also works for the "Jet Smear systematics version"
    // the "isMETPhi" boolean turns on an additional random number that is +/- 1 (since dPhi for the METPhi p.d.f. is a magnitude thing)
    int whichOneDeeHist = TwoDeeHist->GetXaxis()->FindBin(inputMETValue) - 1;
    unsigned int numOneDeeHists = vecOneDeeHists->size();
    if (whichOneDeeHist < 0) {
        std::cout << inputMETValue << std::endl;
        std::cout << "ERROR with which one dee hist!" << std::endl;
        return 0.;
    }
    else if (whichOneDeeHist >= (int) numOneDeeHists) {
        whichOneDeeHist = (int) numOneDeeHists - 1;
    }
    if (vecOneDeeHists->at(whichOneDeeHist)->Integral() == 0) {
        return 0;
    }
    else {
        return vecOneDeeHists->at(whichOneDeeHist)->GetRandom();   
    }
}




#endif //HistSampFunc_h_
