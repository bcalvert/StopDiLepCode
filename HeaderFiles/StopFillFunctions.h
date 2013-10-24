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


typedef std::pair<HistogramT, SampleT> histKey;
typedef std::map<histKey, TH1 *>      HMap_1D;
typedef std::map<histKey, TH2 *>      HMap_2D;
typedef std::map<histKey, TH3 *>      HMap_3D;

float    getMT2(TLorentzVector, TLorentzVector, float, float);

inline float PileupRW(TH1 * nVtxSFHist, int nVtx) {
    int nVtxBin = nVtxSFHist->GetXaxis()->FindBin(nVtx);
    if (nVtxBin > nVtxSFHist->GetNbinsX()) nVtxBin =  nVtxSFHist->GetNbinsX();
    return (float) nVtxSFHist->GetBinContent(nVtxBin);
}
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

inline float MT2lbCalculator(vector<TLorentzVector> * vecLeps, vector<TLorentzVector> * vecJets, float MET, float METPhi, vector<TLorentzVector> &vecBLeps) {
    float MT2lbPair1;
    float MT2lbPair2;
    TLorentzVector vecBLeadLep, vecBSubLep;
    if (vecLeps->size() < 2 || vecJets->size() < 2) {
        cout << "Houston, we've had a problem here: one of the two vectors is less than 2!" << endl;
        cout << "vecLep size " << vecLeps->size() << endl;
        cout << "vecJet size " << vecJets->size() << endl;
    }
    MT2lbPair1 = getMT2(vecLeps->at(0) + vecJets->at(0), vecLeps->at(1) + vecJets->at(1), MET, METPhi);
    MT2lbPair2 = getMT2(vecLeps->at(0) + vecJets->at(1), vecLeps->at(1) + vecJets->at(0), MET, METPhi);
    if (MT2lbPair1 > MT2lbPair2) {
        vecBLeadLep = vecLeps->at(0) + vecJets->at(1);
        vecBSubLep = vecLeps->at(1) + vecJets->at(0);
    }
    else {
        vecBLeadLep = vecLeps->at(0) + vecJets->at(0);
        vecBSubLep = vecLeps->at(1) + vecJets->at(1);               
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
inline void CorrectPlusSetMET(BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo &inputEMI) {
    inputEMI.EventMETX_preCorr = inputEMI.EventMET_preCorr * TMath::Cos(inputEMI.EventMETPhi_preCorr);
    inputEMI.EventMETY_preCorr = inputEMI.EventMET_preCorr * TMath::Sin(inputEMI.EventMETPhi_preCorr);
    inputEMI.EventMETX = inputEMI.EventMETX_preCorr;
    inputEMI.EventMETY = inputEMI.EventMETY_preCorr;
    if (inBEI->doPhiCorr) MetPhiCorrect(inBEI->doData, inputEMI.EventMETX, inputEMI.EventMETY, inBEI->nVtx, inBEI->doReReco);    
    /*
     TRandom3 rand;
     if (doMETSmear) {
     inputEMI.EventMETX *= rand.Gaus(1, METSF * inputEMI.EventMETX);   
     inputEMI.EventMETY *= rand.Gaus(1, METSF * inputEMI.EventMETY);   
     }
     */
    inputEMI.EventMETPhi = TMath::ATan2(inputEMI.EventMETY, inputEMI.EventMETX);
    inputEMI.EventMET = TMath::Sqrt(inputEMI.EventMETX * inputEMI.EventMETX + inputEMI.EventMETY * inputEMI.EventMETY);
    inputEMI.EventMETdivMeff = inputEMI.EventMET / (inputEMI.EventMET + inEJI->EventJetST + inELI->EventLepST);    
}
inline void CalcMT2(EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo &inputEMI) {
    vector<TLorentzVector> vecLepMT2lb(2), vecJetMT2lb(2);
    
    // Set the events M_{T2}(ll)
    inputEMI.EventMT2ll = getMT2(inELI->Lep0.P4, inELI->Lep1.P4, inputEMI.EventMET, inputEMI.EventMETPhi);
    // Set the events M_{T2}(lb)(lb)
    if (inEJI->EventNJets > 1) {
        vecLepMT2lb[0] = inELI->Lep0.P4;
        vecLepMT2lb[1] = inELI->Lep1.P4;
        if (inEJI->EventNBtagJets > 1) {
            vecJetMT2lb[0] = inEJI->BtagJet0.P4;
            vecJetMT2lb[1] = inEJI->BtagJet1.P4;
            inputEMI.caseMT2lb = 0;
        }
        else if (inEJI->EventNBtagJets == 1) {
            vecJetMT2lb[0] = inEJI->BtagJet0.P4;
            if (inEJI->EventBtagJet0Index == 0) {
                vecJetMT2lb[1] = inEJI->Jet1.P4;
                inputEMI.caseMT2lb = 1;
            }
            else {
                vecJetMT2lb[1] = inEJI->Jet0.P4;
                inputEMI.caseMT2lb = 2;
            }            
        }
        else {
            vecJetMT2lb[0] = inEJI->Jet0.P4;
            vecJetMT2lb[1] = inEJI->Jet1.P4;
            inputEMI.caseMT2lb = 3;
        }
        inputEMI.EventMT2lb = MT2lbCalculator(&vecLepMT2lb, &vecJetMT2lb, inputEMI.EventMET, inputEMI.EventMETPhi, inputEMI.EventVecBLepsMT2lb);
        inputEMI.EventDeltaPhiMT2lb_JetsUsed = dPhi(vecJetMT2lb[0].Phi(), vecJetMT2lb[1].Phi());
        inputEMI.EventDeltaPhiMT2lb_BLepsUsed = dPhi(inputEMI.EventVecBLepsMT2lb[0].Phi(), inputEMI.EventVecBLepsMT2lb[1].Phi());
    }
    else {
        inputEMI.EventMT2lb = -99.;
        inputEMI.EventDeltaPhiMT2lb_JetsUsed = -99.;
        inputEMI.EventDeltaPhiMT2lb_BLepsUsed = -99.;
    }
}

inline bool EventPassTrigger(BasicEventInfo * inBEI, EventLepInfo * inELI) {
    bool stillDoEvent;
    if (inELI->doEvent) {
        if (inELI->EventDiLepType == 0) stillDoEvent = inBEI->passTrigDoubleMu;
        else if (inELI->EventDiLepType == 1) stillDoEvent = inBEI->passTrigDoubleEl;
        else if (inELI->EventDiLepType == 2) stillDoEvent = inBEI->passTrigElMu;
        else {
            std::cout << "oddity with EventDiLepType for input ELI" << std::endl;
            stillDoEvent = false;
        }
    }
    return stillDoEvent;
}
inline float nVtxWeight(BasicEventInfo * inBEI, TH2F * nVtxSFHistOviToDESY, TH2F * nVtxSFHist_v2, TH2F * h_S7toS10RWHist, TH2F * nVtxSFHist) {
    float outNVtxWeight = inBEI->weight;
    if (!inBEI->doData) {
        if (inBEI->doPURWOviToDESY) outNVtxWeight *= PileupRW(nVtxSFHistOviToDESY, inBEI->nVtx);
        if (inBEI->doHackPURW) {
            if (inBEI->isSignal) {
                outNVtxWeight = PileupRW(nVtxSFHist_v2, inBEI->nVtx);
                outNVtxWeight *= PileupRW(h_S7toS10RWHist, inBEI->nVtx);
            }
            else {
                outNVtxWeight = PileupRW(nVtxSFHist, inBEI->nVtx);      
            }
        }
    }
    return outNVtxWeight;
}

inline float nVtxWeight(BasicEventInfo * inBEI, TH1F * nVtxSFHistOviToDESY, TH1F * nVtxSFHist_v2, TH1F * h_S7toS10RWHist, TH1F * nVtxSFHist) {
    float outNVtxWeight = inBEI->weight;
    if (!inBEI->doData) {
        if (inBEI->doPURWOviToDESY) outNVtxWeight *= PileupRW(nVtxSFHistOviToDESY, inBEI->nVtx);
        if (inBEI->doHackPURW) {
            if (inBEI->isSignal) {
                outNVtxWeight = PileupRW(nVtxSFHist_v2, inBEI->nVtx);
                outNVtxWeight *= PileupRW(h_S7toS10RWHist, inBEI->nVtx);
            }
            else {
                outNVtxWeight = PileupRW(nVtxSFHist, inBEI->nVtx);      
            }
        }
    }
    return outNVtxWeight;
}

inline void SetEventInformation(BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo &inputEMI, EventDiStructureInfo &inputEDSI) {
    inELI->GrabbedFromTree();
    inEJI->GrabbedFromTree();
    CorrectPlusSetMET(inBEI, inELI, inEJI, inputEMI);
    CalcMT2(inELI, inEJI, inputEMI);
    inputEDSI.SetVars(inELI, inEJI, &inputEMI);
}

inline void SetStringKeyMap(map<string, float> &stringKeyToVar, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI, EventDiStructureInfo * inEDSI, vector<SystT> * systVec, int whichSyst = 0) {
    ///Set up the mapping of string keys to the appropriate event variables
    /*
     list of keys needed can be found via the command:
     cat StopFunctionDefinitions_v2.h | grep VarKey
     _______________________________
     */
    TString SystStringAppend = "";
    if (whichSyst != 0) {
        for (unsigned int iSyst = 0; iSyst < systVec->size(); ++iSyst) {
            if (systVec->at(iSyst).whichSystType == whichSyst) SystStringAppend = TString(systVec->at(iSyst).systVarKey);
        }
    }
    stringKeyToVar[string(TString("leadLepPt") + SystStringAppend)]         = inELI->Lep0.P4.Pt();
    stringKeyToVar[string(TString("leadLepEta") + SystStringAppend)]        = inELI->Lep0.P4.Eta();
    stringKeyToVar[string(TString("lep0RelPFIso") + SystStringAppend)]      = inELI->Lep0.relPFLepIso;
    stringKeyToVar[string(TString("subLepPt") + SystStringAppend)]          = inELI->Lep1.P4.Pt();
    stringKeyToVar[string(TString("subLepEta") + SystStringAppend)]         = inELI->Lep1.P4.Eta();
    stringKeyToVar[string(TString("lep1RelPFIso") + SystStringAppend)]      = inELI->Lep1.relPFLepIso;
    
    stringKeyToVar[string(TString("MT2ll") + SystStringAppend)]                         = inEMI->EventMT2ll;
    stringKeyToVar[string(TString("PassMT2llCut80") + SystStringAppend)]                = (inEMI->EventMT2ll > 80.);
    stringKeyToVar[string(TString("PassMT2llCut90") + SystStringAppend)]                = (inEMI->EventMT2ll > 90.);
    stringKeyToVar[string(TString("PassMT2llCut100") + SystStringAppend)]               = (inEMI->EventMT2ll > 100.);
    stringKeyToVar[string(TString("PassMT2llCut110") + SystStringAppend)]               = (inEMI->EventMT2ll > 110.);
    stringKeyToVar[string(TString("PassMT2llCut120") + SystStringAppend)]               = (inEMI->EventMT2ll > 120.);
    stringKeyToVar[string(TString("MT2lb") + SystStringAppend)]                         = inEMI->EventMT2lb;
    stringKeyToVar[string(TString("MET") + SystStringAppend)]                           = inEMI->EventMET;
    stringKeyToVar[string(TString("METPhi") + SystStringAppend)]                        = inEMI->EventMETPhi;
    stringKeyToVar[string(TString("METPhi_noPhiCorr") + SystStringAppend)]              = inEMI->EventMETPhi_preCorr;    
    stringKeyToVar[string(TString("METX") + SystStringAppend)]                          = inEMI->EventMETX;
    stringKeyToVar[string(TString("METY") + SystStringAppend)]                          = inEMI->EventMETY;
    stringKeyToVar[string(TString("METX_noPhiCorr") + SystStringAppend)]                = inEMI->EventMETX_preCorr;
    stringKeyToVar[string(TString("METY_noPhiCorr") + SystStringAppend)]                = inEMI->EventMETY_preCorr;
    stringKeyToVar[string(TString("METdivMeff") + SystStringAppend)]                    = inEMI->EventMETdivMeff;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut80") + SystStringAppend)]     = inEMI->EventMT2ll > 80 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut90") + SystStringAppend)]     = inEMI->EventMT2ll > 90 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut190") + SystStringAppend)]    = inEMI->EventMT2ll > 100 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut110") + SystStringAppend)]    = inEMI->EventMT2ll > 110 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut120") + SystStringAppend)]    = inEMI->EventMT2ll > 120 ? inEMI->EventMETdivMeff : -99;
    
    stringKeyToVar[string(TString("diLepPt") + SystStringAppend)]               = inEDSI->diLepPt;
    stringKeyToVar[string(TString("diLepInvMass") + SystStringAppend)]          = inEDSI->diLepInvMass;
    stringKeyToVar[string(TString("diLepEta") + SystStringAppend)]              = inEDSI->diLepEta;
    stringKeyToVar[string(TString("diLepPhi") + SystStringAppend)]              = inEDSI->diLepPhi;
    stringKeyToVar[string(TString("DPhiLep0Lep1") + SystStringAppend)]          = inEDSI->DPhiLep0Lep1;
    stringKeyToVar[string(TString("DPhiLep0MET") + SystStringAppend)]           = inEDSI->DPhiLep0MET;
    stringKeyToVar[string(TString("DPhiLep1MET") + SystStringAppend)]           = inEDSI->DPhiLep1MET;
    stringKeyToVar[string(TString("DPhiLep0MET_PreCorr") + SystStringAppend)]   = inEDSI->DPhiLep0MET_PreCorr;
    stringKeyToVar[string(TString("DPhiLep1MET_PreCorr") + SystStringAppend)]   = inEDSI->DPhiLep1MET_PreCorr;
    stringKeyToVar[string(TString("DPhiZMET") + SystStringAppend)]              = inEDSI->DPhiZMET;
    stringKeyToVar[string(TString("DPhiZMET_PreCorr") + SystStringAppend)]      = inEDSI->DPhiZMET_PreCorr;
    
    stringKeyToVar[string(TString("NJets") + SystStringAppend)]     = inEJI->EventNJets;
    stringKeyToVar[string(TString("NBJets") + SystStringAppend)]    = inEJI->EventNBtagJets;
    stringKeyToVar[string(TString("HT") + SystStringAppend)]        = inEJI->EventHT;    
    if (inEJI->EventNJets > 0) {        
        stringKeyToVar[string(TString("leadJetPt") + SystStringAppend)]     = inEJI->Jet0.P4.Pt();
        stringKeyToVar[string(TString("leadJetEta") + SystStringAppend)]    = inEJI->Jet0.P4.Eta();
        stringKeyToVar[string(TString("DPhiLep0Jet0") + SystStringAppend)]  = inEDSI->DPhiLep0Jet0;        
        if (inEJI->EventNBtagJets > 0) {
            stringKeyToVar[string(TString("leadBJetPt") + SystStringAppend)]    = inEJI->BtagJet0.P4.Pt();
            stringKeyToVar[string(TString("leadBJetEta") + SystStringAppend)]   = inEJI->BtagJet0.P4.Eta();
            stringKeyToVar[string(TString("leadBJetEn") + SystStringAppend)]    = inEJI->BtagJet0.P4.E();
            stringKeyToVar[string(TString("DPhiLep0BJet0") + SystStringAppend)] = inEDSI->DPhiLep0BJet0;
            stringKeyToVar[string(TString("DPhiJet0BJet0") + SystStringAppend)] = inEDSI->DPhiJet0BJet0;
            stringKeyToVar[string(TString("DPhiJet1BJet0") + SystStringAppend)] = inEDSI->DPhiJet1BJet0;
        }
        if (inEJI->EventNJets > 1) {
            stringKeyToVar[string(TString("subJetPt") + SystStringAppend)]       = inEJI->Jet1.P4.Pt();
            stringKeyToVar[string(TString("subJetEta") + SystStringAppend)]      = inEJI->Jet1.P4.Eta();                
            stringKeyToVar[string(TString("diJetPt") + SystStringAppend)]        = inEDSI->diJetPt;
            stringKeyToVar[string(TString("diJetInvMass") + SystStringAppend)]   = inEDSI->diJetInvMass;
            stringKeyToVar[string(TString("diJetEta") + SystStringAppend)]       = inEDSI->diJetEta;
            stringKeyToVar[string(TString("diJetPhi") + SystStringAppend)]       = inEDSI->diJetPhi;
            stringKeyToVar[string(TString("DPhiLep0Jet1") + SystStringAppend)]   = inEDSI->DPhiLep0Jet1;
            stringKeyToVar[string(TString("ELepEJet") + SystStringAppend)]       = inEDSI->ELepEJet;
            stringKeyToVar[string(TString("DPhiLepB0LepB1") + SystStringAppend)] = inEDSI->DPhiBLep0BLep1;            
            if (inEJI->EventNBtagJets > 1) {
                stringKeyToVar[string(TString("subBJetPt") + SystStringAppend)]     = inEJI->BtagJet1.P4.Pt();
                stringKeyToVar[string(TString("subBJetEta") + SystStringAppend)]    = inEJI->BtagJet1.P4.Eta(); 
                stringKeyToVar[string(TString("subBJetEn") + SystStringAppend)]     = inEJI->BtagJet1.P4.E();
                stringKeyToVar[string(TString("diBJetPt") + SystStringAppend)]      = inEDSI->diBJetPt;
                stringKeyToVar[string(TString("diBJetInvMass") + SystStringAppend)] = inEDSI->diBJetInvMass;
                stringKeyToVar[string(TString("diBJetEta") + SystStringAppend)]     = inEDSI->diBJetEta;
                stringKeyToVar[string(TString("diBJetPhi") + SystStringAppend)]     = inEDSI->diBJetPhi;                    
                stringKeyToVar[string(TString("DPhiLep0BJet1") + SystStringAppend)] = inEDSI->DPhiLep0BJet1;
                stringKeyToVar[string(TString("DPhiJet1BJet1") + SystStringAppend)] = inEDSI->DPhiJet1BJet1;
            }  
        } 
    }
}


inline void SetStringKeyMapSpecial(map<string, float> &stringKeyToVar, EventMETInfo * inEMI, EventSpecialMT2Info * inESMT2I, vector<SystT> * systVec, int whichSyst = 0) {
    ///Set up the mapping of string keys to the appropriate event variables
    /*
     list of keys needed can be found via the command:
     cat StopFunctionDefinitions_v2.h | grep VarKey
     _______________________________
     */
    float MT2llToUse, MT2lbToUse;
    TString SystStringAppend = "";
    if (whichSyst != 0) {
        for (unsigned int iSyst = 0; iSyst < systVec->size(); ++iSyst) {
            if (systVec->at(iSyst).whichSystType == whichSyst) SystStringAppend = TString(systVec->at(iSyst).systVarKey);
        }
    }
    if (SystStringAppend.Contains("UncESShiftUp")) {
        MT2llToUse = inESMT2I->EventMT2ll_UncESUp;
        MT2lbToUse = inESMT2I->EventMT2lb_UncESUp;
    }
    else if (SystStringAppend.Contains("UncESShiftDown")) {
        MT2llToUse = inESMT2I->EventMT2ll_UncESDown;
        MT2lbToUse = inESMT2I->EventMT2lb_UncESDown;
    }
    else if (SystStringAppend.Contains("MT2llShift")) {
        MT2llToUse = inESMT2I->EventMT2ll_ShiftUp;
        MT2lbToUse = -99.;
    }
    stringKeyToVar[string(TString("MT2ll") + SystStringAppend)]                         = MT2llToUse;
    stringKeyToVar[string(TString("PassMT2llCut80") + SystStringAppend)]                = (MT2llToUse > 80.);
    stringKeyToVar[string(TString("PassMT2llCut90") + SystStringAppend)]                = (MT2llToUse > 90.);
    stringKeyToVar[string(TString("PassMT2llCut100") + SystStringAppend)]               = (MT2llToUse > 100.);
    stringKeyToVar[string(TString("PassMT2llCut110") + SystStringAppend)]               = (MT2llToUse > 110.);
    stringKeyToVar[string(TString("PassMT2llCut120") + SystStringAppend)]               = (MT2llToUse > 120.);
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut80") + SystStringAppend)]     = MT2llToUse > 80 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut90") + SystStringAppend)]     = MT2llToUse > 90 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut190") + SystStringAppend)]    = MT2llToUse > 100 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut110") + SystStringAppend)]    = MT2llToUse > 110 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("METdivMeff_PassMT2llCut120") + SystStringAppend)]    = MT2llToUse > 120 ? inEMI->EventMETdivMeff : -99;
    stringKeyToVar[string(TString("MT2lb") + SystStringAppend)]                         = MT2lbToUse;
}
inline bool LepInBarrel(Lepton * inLep) {
    float barrelEtaEnd = 1.4442;
    return (fabs(inLep->P4.Eta()) < barrelEtaEnd);
}
inline bool LepInEndcap(Lepton * inLep) {
    float endcapEtaStart = 1.566;
    return (fabs(inLep->P4.Eta()) > endcapEtaStart);
}
inline bool DiLeptonGeometry(EventLepInfo * inELI, int whichCase) {
    // whichCase -- which case for individual lepton eta location we care about
    // 0 for both leptons in Barrel
    // 1 for one lepton in Barrel, one lepton not in Barrel (bot not necessarily in the Endcap)
    // 2 for both leptons in Endcap
    bool  Lep0InBarrel = LepInBarrel(&inELI->Lep0);
    bool  Lep0InEndcap = LepInEndcap(&inELI->Lep0);
    bool  Lep1InBarrel = LepInBarrel(&inELI->Lep1);
    bool  Lep1InEndcap = LepInEndcap(&inELI->Lep1);
    bool PassesGeometryCut = true;
    if (whichCase == 0) {
        PassesGeometryCut = Lep0InBarrel && Lep1InBarrel;
    }
    else if (whichCase == 1) {
        PassesGeometryCut = Lep0InBarrel || Lep1InBarrel;
        PassesGeometryCut &= !(Lep0InBarrel && Lep1InBarrel);
    }
    else if (whichCase == 2) {
        PassesGeometryCut = Lep0InEndcap && Lep1InEndcap;
    }
    else {
        cout << "issue with whichCase in DiLeptonGeometry function: " << whichCase << endl;
    }
    return PassesGeometryCut;    
}

inline void SetPassCutMap(std::map<SampleT, bool> &inputCutMap, vector<SampleT> * subSampVec, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI) {
    SampleT S_Current;
    for (unsigned int iSS = 0; iSS < subSampVec->size(); ++iSS) {
        S_Current = subSampVec->at(iSS);
        inputCutMap[S_Current] = false;
        if (!inELI->doEvent) continue;
        if (S_Current.whichdiLepType >= 0 && inELI->EventDiLepType != S_Current.whichdiLepType) continue;
        if (inEJI->EventNJets < S_Current.cutNJets) continue;
        if (inEJI->EventNBtagJets < S_Current.cutNBJets) continue;
        if (!(inELI->EventDiLepType == 2 && S_Current.histNameSuffix.Contains("FullCut"))) {
            if (inEMI->EventMET < S_Current.cutMET) continue;
//            if (S_Current.doZVeto >= 0 && inELI->EventDiLepinZMass == S_Current.doZVeto) continue;
            if (!(S_Current.doZVeto < 0 || inELI->EventDiLepinZMass != S_Current.doZVeto)) continue;
        }
        if (S_Current.histNameSuffix.Contains("BothinBarrel")) {
            if (!DiLeptonGeometry(inELI, 0)) continue;
        }
        if (S_Current.histNameSuffix.Contains("OneinBarrel")) { 
            if (!DiLeptonGeometry(inELI, 1)) continue;
        }   
        if (S_Current.histNameSuffix.Contains("BothinEndcap")) {
            if (!DiLeptonGeometry(inELI, 2)) continue;
        }
        if (S_Current.histNameSuffix.Contains("0BJets")) {
            if (inEJI->EventNBtagJets > 0) {
                if (!(S_Current.histNameSuffix.Contains("inZMass") && inELI->EventDiLepinZMass)) continue;
            }
        }
        if (S_Current.histNameSuffix.Contains("_0Jets") && inEJI->EventNJets != 0) continue;
        if (S_Current.histNameSuffix.Contains("_1Jet") && inEJI->EventNJets != 1) continue;
        if (S_Current.histNameSuffix.Contains("FullCutBlind")) {
            if (inEMI->EventMT2ll > 80) continue;
        }
        inputCutMap[S_Current] = true;
//        cout << "set inputCut Map to be true for S_Current = " << S_Current.histNameSuffix << endl;
    }    
}
inline void HistogramFillOneDee(std::map<SampleT, bool> * inputCutMap, map<string, float> * stringKeyToVar, vector<SampleT> * subSampVec, vector<HistogramT> * HistTVec, HMap_1D * inputHMap1D, BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI, EventDiStructureInfo * inEDSI, bool doVerbosity) {
    SampleT S_Current;
    HistogramT H_Current;
    float fillWeight;
    map<string, float>::iterator xIter;
//    map<string, float>::iterator yIter;
//    map<string, float>::iterator zIter;
    float MT2llCut = 80;
    float MT2lbCut = 170;
    const double PI = 3.14159265;
    for (unsigned int iSS = 0; iSS < subSampVec->size(); ++iSS) {
        S_Current = subSampVec->at(iSS);
//        cout << "S_Current = " << S_Current.histNameSuffix << endl;
        if ((*inputCutMap)[S_Current]) {
//            cout << "passed inputCut Map for S_Current = " << S_Current.histNameSuffix << endl;
            for (unsigned int iHT = 0; iHT < HistTVec->size(); ++iHT) {
                H_Current = HistTVec->at(iHT);
                xIter = stringKeyToVar->find(H_Current.xVarKey);
                /*
                if (doVerbosity) {
                    cout << "" << endl;
                    cout << "ievt " << ievt << endl;
                    cout << "i " << i << endl;
                    cout << "j " << j << endl;
                    cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                    cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                }
                */
                if (xIter != stringKeyToVar->end()) {
                    if (doVerbosity) {
                        cout << "xIter first " << xIter->first << endl;
                        cout << "xIter second " << xIter->second << endl;
                    }
                    ///Some necessary continue checks
                    if (S_Current.blindDataChannel && TString(H_Current.xVarKey).Contains("MT2ll")) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2ll > MT2llCut) continue;
                    }
                    if (S_Current.blindDataChannel && TString(H_Current.xVarKey).Contains("MT2lb")) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2lb > MT2lbCut) continue;
                    }
                    if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                        if (inEDSI->DPhiZMET > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                        if (inEDSI->DPhiZMET < 1./3. * PI || inEDSI->DPhiZMET > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                        if (inEDSI->DPhiZMET < 2./3. * PI) continue;
                    }                        
                    if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                        if (inEDSI->DPhiLep0Lep1 > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                        if (inEDSI->DPhiLep0Lep1 < 1./3. * PI || inEDSI->DPhiLep0Lep1 > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                        if (inEDSI->DPhiLep0Lep1 < 2./3. * PI) continue;
                    }
                    if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                        if (inEMI->EventDeltaPhiMT2lb_JetsUsed > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                        if (inEMI->EventDeltaPhiMT2lb_JetsUsed < 1./3. * PI || inEMI->EventDeltaPhiMT2lb_JetsUsed > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                        if (inEMI->EventDeltaPhiMT2lb_JetsUsed < 2./3. * PI) continue;
                    }
                    if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                        if (inEDSI->DPhiBLep0BLep1 > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                        if (inEDSI->DPhiBLep0BLep1 < 1./3. * PI || inEDSI->DPhiBLep0BLep1 > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                        if (inEDSI->DPhiBLep0BLep1 < 2./3. * PI) continue;
                    }
                    fillWeight = H_Current.name.Contains("preRW") ? inBEI->preNVtxRWWeight : inBEI->weight;
                    if (H_Current.name.Contains("h_ChannelCutFlow")) fillWeight = 1.;
                    (*inputHMap1D)[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                }
                if (doVerbosity) cout << "" << endl;
            }
        }
    }
}

inline void HistogramFillOneDeeSyst(std::map<SampleT, bool> * inputCutMap, map<string, float> * stringKeyToVar, vector<SampleT> * subSampVec, vector<HistogramT> * HistTVec, HMap_1D * inputHMap1D, BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI, EventDiStructureInfo * inEDSI, TString SystCheck, bool doVerbosity) {
    SampleT S_Current;
    HistogramT H_Current;
    float fillWeight;
    map<string, float>::iterator xIter;
    //    map<string, float>::iterator yIter;
    //    map<string, float>::iterator zIter;
    float MT2llCut = 80;
    float MT2lbCut = 170;
    const double PI = 3.14159265;
    for (unsigned int iSS = 0; iSS < subSampVec->size(); ++iSS) {
        S_Current = subSampVec->at(iSS);
        if ((*inputCutMap)[S_Current]) {
            for (unsigned int iHT = 0; iHT < HistTVec->size(); ++iHT) {
                H_Current = HistTVec->at(iHT);
                if (!H_Current.name.Contains(SystCheck)) continue;                
                xIter = stringKeyToVar->find(H_Current.xVarKey);
                /*
                if (doVerbosity) {
                    cout << "" << endl;
                    cout << "ievt " << ievt << endl;
                    cout << "i " << i << endl;
                    cout << "j " << j << endl;
                    cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                    cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                }
                */
                if (xIter != stringKeyToVar->end()) {
                    if (doVerbosity) {
                        cout << "xIter first " << xIter->first << endl;
                        cout << "xIter second " << xIter->second << endl;
                    }
                    ///Some necessary continue checks
                    if (S_Current.blindDataChannel && TString(H_Current.xVarKey).Contains("MT2ll")) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2ll > MT2llCut) continue;
                    }
                    if (S_Current.blindDataChannel && TString(H_Current.xVarKey).Contains("MT2lb")) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2lb > MT2lbCut) continue;
                    }
                    if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                        if (inEDSI->DPhiZMET > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                        if (inEDSI->DPhiZMET < 1./3. * PI || inEDSI->DPhiZMET > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                        if (inEDSI->DPhiZMET < 2./3. * PI) continue;
                    }                        
                    if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                        if (inEDSI->DPhiLep0Lep1 > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                        if (inEDSI->DPhiLep0Lep1 < 1./3. * PI || inEDSI->DPhiLep0Lep1 > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                        if (inEDSI->DPhiLep0Lep1 < 2./3. * PI) continue;
                    }
                    if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                        if (inEMI->EventDeltaPhiMT2lb_JetsUsed > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                        if (inEMI->EventDeltaPhiMT2lb_JetsUsed < 1./3. * PI || inEMI->EventDeltaPhiMT2lb_JetsUsed > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                        if (inEMI->EventDeltaPhiMT2lb_JetsUsed < 2./3. * PI) continue;
                    }
                    if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                        if (inEDSI->DPhiBLep0BLep1 > 1./3. * PI) continue;
                    }
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                        if (inEDSI->DPhiBLep0BLep1 < 1./3. * PI || inEDSI->DPhiBLep0BLep1 > 2./3. * PI) continue;                         
                    }                        
                    else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                        if (inEDSI->DPhiBLep0BLep1 < 2./3. * PI) continue;
                    }
                    fillWeight = H_Current.name.Contains("preRW") ? inBEI->preNVtxRWWeight : inBEI->weight;
                    if (H_Current.name.Contains("h_ChannelCutFlow")) fillWeight = 1.;
                    (*inputHMap1D)[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                }
                if (doVerbosity) cout << "" << endl;
            }
        }
    }
}
inline void HistogramFillTwoDee(std::map<SampleT, bool> * inputCutMap, map<string, float> * stringKeyToVar, vector<SampleT> * subSampVec, vector<HistogramT> * HistTVec, HMap_2D * inputHMap2D, BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI, EventDiStructureInfo * inEDSI, bool doVerbosity) {
    SampleT S_Current;
    HistogramT H_Current;
    float fillWeight;
    map<string, float>::iterator xIter;
    map<string, float>::iterator yIter;
    //    map<string, float>::iterator zIter;
    float MT2llCut = 80;
    float MT2lbCut = 170;
    const double PI = 3.14159265;
    for (unsigned int iSS = 0; iSS < subSampVec->size(); ++iSS) {
        S_Current = subSampVec->at(iSS);
        if ((*inputCutMap)[S_Current]) {
            
            for (unsigned int iHT = 0; iHT < HistTVec->size(); ++iHT) {
                H_Current = HistTVec->at(iHT);
                xIter = stringKeyToVar->find(H_Current.xVarKey);
                yIter = stringKeyToVar->find(H_Current.yVarKey);
                if (xIter != stringKeyToVar->end() && yIter != stringKeyToVar->end()) { 
                    //                        if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2ll") || TString(H_Current.yVarKey).Contains("MT2ll"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2ll > MT2llCut) continue;
                    }
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2lb") || TString(H_Current.yVarKey).Contains("MT2lb"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2lb > MT2lbCut) continue;
                    }
                    fillWeight = H_Current.name.Contains("preRW") ? inBEI->preNVtxRWWeight : inBEI->weight;
                    (*inputHMap2D)[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, fillWeight);
                }
            }
        }
    }
}
inline void HistogramFillTwoDeeSyst(std::map<SampleT, bool> * inputCutMap, map<string, float> * stringKeyToVar, vector<SampleT> * subSampVec, vector<HistogramT> * HistTVec, HMap_2D * inputHMap2D, BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI, EventDiStructureInfo * inEDSI, TString SystCheck, bool doVerbosity) {
    SampleT S_Current;
    HistogramT H_Current;
    float fillWeight;
    map<string, float>::iterator xIter;
    map<string, float>::iterator yIter;
    //    map<string, float>::iterator zIter;
    float MT2llCut = 80;
    float MT2lbCut = 170;
    const double PI = 3.14159265;
    for (unsigned int iSS = 0; iSS < subSampVec->size(); ++iSS) {
        S_Current = subSampVec->at(iSS);
        if ((*inputCutMap)[S_Current]) {            
            for (unsigned int iHT = 0; iHT < HistTVec->size(); ++iHT) {
                H_Current = HistTVec->at(iHT);
                if (!H_Current.name.Contains(SystCheck)) continue;
                xIter = stringKeyToVar->find(H_Current.xVarKey);
                yIter = stringKeyToVar->find(H_Current.yVarKey);
                if (xIter != stringKeyToVar->end() && yIter != stringKeyToVar->end()) { 
                    //                        if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2ll") || TString(H_Current.yVarKey).Contains("MT2ll"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2ll > MT2llCut) continue;
                    }
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2lb") || TString(H_Current.yVarKey).Contains("MT2lb"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2lb > MT2lbCut) continue;
                    }
                    fillWeight = H_Current.name.Contains("preRW") ? inBEI->preNVtxRWWeight : inBEI->weight;
                    (*inputHMap2D)[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, fillWeight);
                }
            }
        }
    }
}

inline void HistogramFillThreeDee(std::map<SampleT, bool> * inputCutMap, map<string, float> * stringKeyToVar, vector<SampleT> * subSampVec, vector<HistogramT> * HistTVec, HMap_3D * inputHMap3D, BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI, EventDiStructureInfo * inEDSI, bool doVerbosity) {
    SampleT S_Current;
    HistogramT H_Current;
    float fillWeight;
    map<string, float>::iterator xIter;
    map<string, float>::iterator yIter;
    map<string, float>::iterator zIter;
    float MT2llCut = 80;
    float MT2lbCut = 170;
    const double PI = 3.14159265;
    for (unsigned int iSS = 0; iSS < subSampVec->size(); ++iSS) {
        S_Current = subSampVec->at(iSS);
        if ((*inputCutMap)[S_Current]) {            
            for (unsigned int iHT = 0; iHT < HistTVec->size(); ++iHT) {
                H_Current = HistTVec->at(iHT);
                xIter = stringKeyToVar->find(H_Current.xVarKey);
                yIter = stringKeyToVar->find(H_Current.yVarKey);
                zIter = stringKeyToVar->find(H_Current.zVarKey);
                if (xIter != stringKeyToVar->end() && yIter != stringKeyToVar->end() && zIter != stringKeyToVar->end()) {
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2ll") || TString(H_Current.yVarKey).Contains("MT2ll") || TString(H_Current.zVarKey).Contains("MT2ll"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2ll > MT2llCut) continue;
                    }
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2lb") || TString(H_Current.yVarKey).Contains("MT2lb") || TString(H_Current.zVarKey).Contains("MT2lb"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2lb > MT2lbCut) continue;
                    }
                    if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                        if (H_Current.name.Contains("21to30")) {
                            if ((inBEI->nVtx > 30 || inBEI->nVtx < 21)) continue;
                        }
                        else if (H_Current.name.Contains("11to20")) {
                            if ((inBEI->nVtx > 20 || inBEI->nVtx < 11)) continue;
                        }
                        else if (H_Current.name.Contains("1to10")) {
                            if ((inBEI->nVtx > 10 || inBEI->nVtx < 1)) continue;
                        }
                    }
                    fillWeight = H_Current.name.Contains("preRW") ? inBEI->preNVtxRWWeight : inBEI->weight;
                    (*inputHMap3D)[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, fillWeight);
                }
            }
        }
    }
}


inline void HistogramFillThreeDeeSyst(std::map<SampleT, bool> * inputCutMap, map<string, float> * stringKeyToVar, vector<SampleT> * subSampVec, vector<HistogramT> * HistTVec, HMap_3D * inputHMap3D, BasicEventInfo * inBEI, EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI, EventDiStructureInfo * inEDSI, TString SystCheck, bool doVerbosity) {
    SampleT S_Current;
    HistogramT H_Current;
    float fillWeight;
    map<string, float>::iterator xIter;
    map<string, float>::iterator yIter;
    map<string, float>::iterator zIter;
    float MT2llCut = 80;
    float MT2lbCut = 170;
    const double PI = 3.14159265;
    for (unsigned int iSS = 0; iSS < subSampVec->size(); ++iSS) {
        S_Current = subSampVec->at(iSS);
        if ((*inputCutMap)[S_Current]) {            
            for (unsigned int iHT = 0; iHT < HistTVec->size(); ++iHT) {
                H_Current = HistTVec->at(iHT);
                if (!H_Current.name.Contains(SystCheck)) continue;
                xIter = stringKeyToVar->find(H_Current.xVarKey);
                yIter = stringKeyToVar->find(H_Current.yVarKey);
                zIter = stringKeyToVar->find(H_Current.zVarKey);
                if (xIter != stringKeyToVar->end() && yIter != stringKeyToVar->end() && zIter != stringKeyToVar->end()) {
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2ll") || TString(H_Current.yVarKey).Contains("MT2ll") || TString(H_Current.zVarKey).Contains("MT2ll"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2ll > MT2llCut) continue;
                    }
                    if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2lb") || TString(H_Current.yVarKey).Contains("MT2lb") || TString(H_Current.zVarKey).Contains("MT2lb"))) {
                        if (inBEI->blindData && inBEI->doData && inEMI->EventMT2lb > MT2lbCut) continue;
                    }
                    if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                        if (H_Current.name.Contains("21to30")) {
                            if ((inBEI->nVtx > 30 || inBEI->nVtx < 21)) continue;
                        }
                        else if (H_Current.name.Contains("11to20")) {
                            if ((inBEI->nVtx > 20 || inBEI->nVtx < 11)) continue;
                        }
                        else if (H_Current.name.Contains("1to10")) {
                            if ((inBEI->nVtx > 10 || inBEI->nVtx < 1)) continue;
                        }
                    }
                    fillWeight = H_Current.name.Contains("preRW") ? inBEI->preNVtxRWWeight : inBEI->weight;
                    (*inputHMap3D)[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, fillWeight);
                }
            }
        }
    }
}