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
#include "../PlotMakingCode/mt2bisect.h"

#include <iostream>
//#include <vector>
#include <cmath>
#include <sstream>
#include <map>
using namespace std;
inline float getMT2(TLorentzVector lept1, TLorentzVector lept2, float theMET, float theMETphi){
    // Calculate MT2 variable for two leptons and missing energy, assuming zero testmass                                                                       
    double pa[3];
    double pb[3];
    double pmiss[3];
    
    TLorentzVector pmet;
    pmet.SetPtEtaPhiM(theMET, 0., theMETphi, 0.);
    pmiss[0] = 0.; // irrelevant                                                                                                                               
    pmiss[1] = pmet.Px();
    pmiss[2] = pmet.Py();
    
    pa[0] = 0.;
    pa[1] = lept1.Px();
    pa[2] = lept1.Py();
    
    pb[0] = 0.;
    pb[1] = lept2.Px();
    pb[2] = lept2.Py();
    
    mt2bisect* MT2bisect = new mt2bisect();
    MT2bisect->set_verbose(0);
    MT2bisect->set_momenta(pa, pb, pmiss);
    MT2bisect->set_mn(0.); // testmass                                                                                                                         
    double MT2 = MT2bisect->get_mt2();
    delete MT2bisect;
    return MT2;
}
inline void MetPhiCorrect(bool doData, float &MetX, float &MetY, int nVtx, bool doReReco) {
    int c0 = 0;
    int c1 = 1;
    if (!doData) {
        c0 = 2;
        c1 = 3;
    }
//    float CorrType1PFMETX[4] = {-0.165731, 0.393398, 0.0637551, 0.000364516};
//    float CorrType1PFMETY[4] = {-0.00360751, -0.253457, 0.430063, -0.228208};
    float CorrType1PFMETX[4] = {-0.067, 0.387, 0.110, -0.005};
    float CorrType1PFMETY[4] = {-0.336, -0.247, 0.321, -0.214};
    
    float CorrType1PFMETX_ReReco[4] = {-0.254, 0.309, 0.110, -0.005};
    float CorrType1PFMETY_ReReco[4] = {-0.158, -0.160, 0.321, -0.214};
    if (doReReco) {
        MetX = MetX - CorrType1PFMETX_ReReco[c0] - nVtx * CorrType1PFMETX_ReReco[c1];
        MetY = MetY - CorrType1PFMETY_ReReco[c0] - nVtx * CorrType1PFMETY_ReReco[c1];
    }
    else {
        MetX = MetX - CorrType1PFMETX[c0] - nVtx * CorrType1PFMETX[c1];
        MetY = MetY - CorrType1PFMETY[c0] - nVtx * CorrType1PFMETY[c1];
    }
}
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

inline void METSystShift(vector<Lepton> * inputLepVec, vector<Lepton> * shiftLepVec, float &newMET, float &newMETPhi, float origMET, float origMETPhi) {
    TVector3 vecMET;
    TVector3 inputLepTotThreeVec; inputLepTotThreeVec.SetPtEtaPhi(0., 0., 0.);
    TVector3 shiftLepTotThreeVec; shiftLepTotThreeVec.SetPtEtaPhi(0., 0., 0.);
    TVector3 inputLepCurrThreeVec;
    TVector3 shiftLepCurrThreeVec;
    vecMET.SetPtEtaPhi(origMET, 0., origMETPhi);    
    if (inputLepVec->size() != shiftLepVec->size()) {
        cout << "inputLep size " << inputLepVec->size() << endl;
        cout << "shiftLep size " << shiftLepVec->size() << endl;
        cout << "issue: two input object vectors are not same size -- check it out" << endl;
        return;
    }
    for (unsigned int i = 0; i < inputLepVec->size(); ++i) {
        inputLepCurrThreeVec.SetPtEtaPhi(inputLepVec->at(i).P4.Pt(), inputLepVec->at(i).P4.Eta(), inputLepVec->at(i).P4.Phi());
        shiftLepCurrThreeVec.SetPtEtaPhi(shiftLepVec->at(i).P4.Pt(), shiftLepVec->at(i).P4.Eta(), shiftLepVec->at(i).P4.Phi());
        inputLepTotThreeVec = inputLepTotThreeVec + inputLepCurrThreeVec;
        shiftLepTotThreeVec = shiftLepTotThreeVec + shiftLepCurrThreeVec;
    }
    vecMET = vecMET - inputLepTotThreeVec;
    vecMET = vecMET + shiftLepTotThreeVec;
    newMET = vecMET.Pt();
    newMETPhi = vecMET.Phi();
    return;
}

inline void METSystShift(vector<PFJet> * inputJetVec, vector<PFJet> * shiftJetVec, float &newMET, float &newMETPhi, float origMET, float origMETPhi) {
    TVector3 vecMET;
    TVector3 inputJetTotThreeVec; inputJetTotThreeVec.SetPtEtaPhi(0., 0., 0.);
    TVector3 shiftJetTotThreeVec; shiftJetTotThreeVec.SetPtEtaPhi(0., 0., 0.);
    TVector3 inputJetCurrThreeVec;
    TVector3 shiftJetCurrThreeVec;
    vecMET.SetPtEtaPhi(origMET, 0., origMETPhi);    
    if (inputJetVec->size() != shiftJetVec->size()) {
        cout << "inputJet size " << inputJetVec->size() << endl;
        cout << "shiftJet size " << shiftJetVec->size() << endl;
        cout << "issue: two input object vectors are not same size -- check it out" << endl;
        return;
    }
    for (unsigned int i = 0; i < inputJetVec->size(); ++i) {
        inputJetCurrThreeVec.SetPtEtaPhi(inputJetVec->at(i).P4.Pt(), inputJetVec->at(i).P4.Eta(), inputJetVec->at(i).P4.Phi());
        shiftJetCurrThreeVec.SetPtEtaPhi(shiftJetVec->at(i).P4.Pt(), shiftJetVec->at(i).P4.Eta(), shiftJetVec->at(i).P4.Phi());
        inputJetTotThreeVec = inputJetTotThreeVec + inputJetCurrThreeVec;
        shiftJetTotThreeVec = shiftJetTotThreeVec + shiftJetCurrThreeVec;
    }
    vecMET = vecMET - inputJetTotThreeVec;
    vecMET = vecMET + shiftJetTotThreeVec;
    newMET = vecMET.Pt();
    newMETPhi = vecMET.Phi();
    return;
}

inline float MT2lbCalculator(vector<TLorentzVector> &vecLeps, vector<TLorentzVector> &vecJets, float MET, float METPhi, vector<TLorentzVector> &vecBLeps) {
    float MT2lbPair1;
    float MT2lbPair2;
    TLorentzVector vecBLeadLep, vecBSubLep;
    if (vecLeps.size() < 2 || vecJets.size() < 2) {
        cout << "Houston, we've had a problem here: one of the two vectors is less than 2!" << endl;
        cout << "vecLep size " << vecLeps.size() << endl;
        cout << "vecJet size " << vecJets.size() << endl;
    }
    MT2lbPair1 = getMT2(vecLeps[0] + vecJets[0], vecLeps[1] + vecJets[1], MET, METPhi);
    MT2lbPair2 = getMT2(vecLeps[0] + vecJets[1], vecLeps[1] + vecJets[0], MET, METPhi);
    if (MT2lbPair1 > MT2lbPair2) {
        vecBLeadLep = vecLeps[0] + vecJets[1];
        vecBSubLep = vecLeps[1] + vecJets[0];
    }
    else {
        vecBLeadLep = vecLeps[0] + vecJets[0];
        vecBSubLep = vecLeps[1] + vecJets[1];               
    }
    vecBLeps[0] = vecBLeadLep;
    vecBLeps[1] = vecBSubLep;
    return TMath::Min(MT2lbPair1, MT2lbPair2);
}
inline float GenLevelTopPtWeight(float pT_Top, float pT_AntiTop) {
    if (pT_Top < 0 || pT_AntiTop < 0) {
        cout << "bad mojo -- either pT is less than 0 " << endl;
        cout << "pT_Top " << pT_Top << endl;
        cout << "pT_AntiTop " << pT_AntiTop << endl;
        return 0.;
    }
    //    float pT_ToUse = TMath::Sqrt(pT_Top * pT_AntiTop);
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
inline float ScaleFactorMC(int Type, int Syst) {
    
    // Type: 0 MuMu, 1 EE 2 EMu
    // Syst: 0 CentVal, 1 SystShiftUp, 2 SystShiftDown
    //    float SFTrig[3] = {0.994, 0.955, 0.978};
    //    float SFIDIso[3] = {0.993, 0.979, 0.986};    
    //return SFTrig[Type] * SFIDIso[Type];
    // Note! As of 8/10/13, trying an additional SF just for funsies, ok not really just for funsies, basically we find that the DD TTBar normalization is different for the three separate channels, which is bad news bears because it is indicative of different lepton reconstruction efficiency scale factors for data/MC for the different leptons
    /// 
    /*
     float SF_trigger_mumu=0.965;// +/- 0.0102;
     float SF_trigger_ee  =0.962;// +/- 0.0130;
     float SF_trigger_mue =0.943;// +/- 0.0120;
     //02-11-2012
     float SF_IDISO_mumu=0.997;// +/- 0.00085;
     float SF_IDISO_ee  =0.975;// +/- 0.0006;
     float SF_IDISO_mue =0.986;// +/- 0.0007;
     *////
    float SFTotal[3] = {0.962105, 0.9379500, 0.9297980};
    float SFSystUp[3] = {0.0102024, 0.0126881, 0.0185040};
    float SFSystDown[3] = {0.0102024, 0.0126881, 0.0185040};
    /*
     float SFTotal[3] = {0.987, 0.957, 0.935};
     float SFSystUp[3] = {0.011, 0.015, 0.013};
     float SFSystDown[3] = {0.011, 0.015, 0.013};
     float MuMuAdditionalSF = 1.02169370264668;
     float EEAddtionalSF    = 0.977225138882017;
     float EMuAdditionalSF  = 0.999212074818846;
     
     float SFAdditional[3] = {MuMuAdditionalSF, EEAddtionalSF, EMuAdditionalSF};
     switch (Syst) {
     case 0:
     return SFTotal[Type] * SFAdditional[Type];
     break;
     case 1:
     return (SFTotal[Type] * SFAdditional[Type]) + SFSystUp[Type];
     break;      
     case -1:
     return (SFTotal[Type] * SFAdditional[Type]) - SFSystDown[Type];
     break;  
     default:
     return 1;
     break;
     }
     */
    
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
/*
float FastSimScaleFactor(Lepton lep0, Lepton lep1) {
// scale factors taken from slideshttps://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=275431
}
*/