
#include "TChain.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TF1.h"
#include "TTree.h"
#include "TString.h"
#include "TRegexp.h"
#include "TProfile.h"

//#include <exception>                                                                                                      
#include <sys/stat.h>

#include "TCut.h"
//#include "mt2bisect.h"
//#include "StopDict_ver2.h"
#include "../../HeaderFiles/StopFunctionDefinitions_v2.h"

//#include "PileUpMC.h"

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <map>

using namespace std;
void MetPhiCorrect(bool doData, float &MetX, float &MetY, int nVtx) {
    int c0 = 0;
    int c1 = 1;
    if (!doData) {
        c0 = 2;
        c1 = 3;
    }
    float CorrType1PFMETX[4] = {-0.165731, 0.393398, 0.0637551, 0.000364516};
    float CorrType1PFMETY[4] = {-0.00360751, -0.253457, 0.430063, -0.228208};
    MetX = MetX - CorrType1PFMETX[c0] - nVtx * CorrType1PFMETX[c1];
    MetY = MetY - CorrType1PFMETY[c0] - nVtx * CorrType1PFMETY[c1];
}
int main( int argc, const char* argv[] ) {
    /////////Variable initializations/////////////////
    /////Organization Variables//////////////   
    TFile * inputWeightFile = new TFile("StopABCPUWeight_vers2.root");
    TH3F * weights = (TH3F*) inputWeightFile->Get("WHist");
    TH1F * OneDPUDist = (TH1F*) weights->ProjectionX("OneDWeight");
    TFile * inputPUFileMC = new TFile(TString("OneDMCPURW.root"));
    TFile * inputPUFileMCOviToDESY = new TFile(TString("OneDMCPURW_OviToDESY.root"));
    //    TH1F * truePUDistMC = (TH1F*) inputPUFileMC->Get("MCPU");
    //    TFile * inputPUFileData = new TFile("RunABCPUDist.root");
    //    TH1F * truePUDistData = (TH1F*) inputPUFileData->Get("pileup");
    //    float intMCPU = truePUDistMC->
    TH1F * nVtxSFHist = (TH1F*) inputPUFileMC->Get("nVtxSF_preRW");
    TH1F * nVtxSFHistOviToDESY = (TH1F*) inputPUFileMCOviToDESY->Get("normFrac");
    TString fileTreeName;    
    TString fInName;
    TString fOutName;    
    TH1F * h_eventCount = new TH1F("h_eventCount", "; numEvents (# of Entries);", 2, -0.5, 1.5);
    
    /////Event Variables/////////////////////
    
    //////bunch of vectors
    //Muons
    vector<float> * MuonPt, * MuonEta, * MuonPx, * MuonPy, * MuonPz, * MuonE, * MuonPFCharIso, * MuonPFNeutIso, * MuonPFPhotIso, * MuonSumPUPt, * MuonD0, * MuonVertZ; //ambiguity in which Muon Pt...rolling with PFMuonPt
    vector<bool> * isPFMuon, * isGMPTMuons;
    //Electrons
    vector<float> * ElecPt, * ElecEta, * ElecPx, * ElecPy, * ElecPz, * ElecE, * ElecPFCharIso, * ElecPFNeutIso, * ElecPFPhotIso; //ambiguity in which Electron Pt...rolling with PFElecPt
    vector<bool> * isPFElectron, *passConvVeto;
    //Jets
    vector<float> * JetEt, *JetEta, * JetPx, * JetPy, * JetPz, * JetE, * JetNHF, * JetNEF, * JetCHF, * JetCEF, * JetBTag;
    vector<int> * JetNDaug, * JetCharMult;
    
    vector<int> * MuonCharge, * ElecCharge;
    
    bool filterECalDeadCell, filterLogErrorTooManyClusters, filterTrackingFailure, filterHCalLaser, filterECalLaserCorr, filterTooManyStripClus, filterManyStripClus, filterEEBadSC;    
    vector<bool> * regularMETFilterVec = new vector<bool>;
    vector<bool> * oppositeMETFilterVec = new vector<bool>;
    
    //Vertex info (both for nVtx info and also muon dZ)
    vector<float> * VertZ;    
    
    MuonPt = new vector<float>;
    MuonEta = new vector<float>;
    MuonPx = new vector<float>;
    MuonPy = new vector<float>;
    MuonPz = new vector<float>;
    MuonE = new vector<float>;
    MuonCharge = new vector<int>;
    MuonPFCharIso = new vector<float>;
    MuonPFNeutIso = new vector<float>;
    MuonPFPhotIso = new vector<float>;
    MuonSumPUPt = new vector<float>;
    MuonD0 = new vector<float>;
    MuonVertZ = new vector<float>; 
    
    isPFMuon = new vector<bool>;
    isGMPTMuons = new vector<bool>;
    
    ElecPt = new vector<float>;
    ElecEta = new vector<float>;
    ElecPx = new vector<float>;
    ElecPy = new vector<float>;
    ElecPz = new vector<float>;
    ElecE = new vector<float>;
    ElecCharge = new vector<int>;
    ElecPFCharIso = new vector<float>;
    ElecPFNeutIso = new vector<float>;
    ElecPFPhotIso = new vector<float>;
    
    isPFElectron = new vector<bool>;    
    passConvVeto = new vector<bool>;
    
    JetEt = new vector<float>;
    JetEta = new vector<float>;
    JetPx = new vector<float>;
    JetPy = new vector<float>;
    JetPz = new vector<float>;
    JetE = new vector<float>;
    JetNHF = new vector<float>;
    JetNEF = new vector<float>;
    JetCHF = new vector<float>;
    JetCEF = new vector<float>;
    JetBTag = new vector<float>;
    
    JetNDaug = new vector<int>;
    JetCharMult = new vector<int>;
    
    VertZ = new vector<float>;
    
    int   Type, nVtx, nVtxTrue;
    int   RunNum, EventNum, LumiBlock;
    float MET,MET_Phi,METSig;
    float METX, METY;
    float Lep0Px,Lep0Py,Lep0Pz,Lep0E,Lep1Px,Lep1Py,Lep1Pz,Lep1E;
    int   Lep0PDGID, Lep1PDGID;
    int   NJets, NBtagJets;
    float HT, Btag_j0,Btag_j1;
    float Jet0Px,Jet0Py,Jet0Pz,Jet0E,Jet0Et,Jet1Px,Jet1Py,Jet1Pz,Jet1E,Jet1Et;
    float BtagJet0Px,BtagJet0Py,BtagJet0Pz,BtagJet0E,BtagJet0Et,BtagJet1Px,BtagJet1Py,BtagJet1Pz,BtagJet1E,BtagJet1Et;
    int   BtagJet0Index, BtagJet1Index;
    float BDT,BDTDis;
    
    int genStopMass0, genStopMass1, genChi0Mass0, genChi0Mass1, genCharginoMass0, genCharginoMass1;
    float genStopMassCut, genChi0MassCut, genCharginoMassCut;
    genStopMassCut = 0;
    genChi0MassCut = 0;
    genCharginoMassCut = 0;
    // Gen info
    vector<float> * genStopMass, * genChi0Mass, * genCharginoMass;
    vector<double> * genPolWeights;
    genStopMass = new vector<float>;
    genChi0Mass = new vector<float>;
    genCharginoMass = new vector<float>;
    genPolWeights = new vector<double>;  
    
    float StopMass0, StopMass1, Chi0Mass0, Chi0Mass1;
    // Leading Status 3 particle characteristics
    float genTopSt3_0_Energy, genTopSt3_0_Pt, genTopSt3_0_Eta, genTopSt3_0_Phi;
    int   genTopSt3_0_PID, genTopSt3_0_Index, genTopSt3_0_FirstMom;
    float genBSt3_0_Energy, genBSt3_0_Pt, genBSt3_0_Eta, genBSt3_0_Phi;
    int   genBSt3_0_PID, genBSt3_0_Index, genBSt3_0_FirstMom;
    float genMuonSt3_0_Energy, genMuonSt3_0_Pt, genMuonSt3_0_Eta, genMuonSt3_0_Phi;
    int   genMuonSt3_0_PID, genMuonSt3_0_Index, genMuonSt3_0_FirstMom;
    float genElecSt3_0_Energy, genElecSt3_0_Pt, genElecSt3_0_Eta, genElecSt3_0_Phi;
    int   genElecSt3_0_PID, genElecSt3_0_Index, genElecSt3_0_FirstMom;
    
    // Sub-Leading Status 3 particle characteristics
    float genTopSt3_1_Energy, genTopSt3_1_Pt, genTopSt3_1_Eta, genTopSt3_1_Phi;
    int   genTopSt3_1_PID, genTopSt3_1_Index, genTopSt3_1_FirstMom;
    float genBSt3_1_Energy, genBSt3_1_Pt, genBSt3_1_Eta, genBSt3_1_Phi;
    int   genBSt3_1_PID, genBSt3_1_Index, genBSt3_1_FirstMom;
    float genMuonSt3_1_Energy, genMuonSt3_1_Pt, genMuonSt3_1_Eta, genMuonSt3_1_Phi;
    int   genMuonSt3_1_PID, genMuonSt3_1_Index, genMuonSt3_1_FirstMom;
    float genElecSt3_1_Energy, genElecSt3_1_Pt, genElecSt3_1_Eta, genElecSt3_1_Phi;
    int   genElecSt3_1_PID, genElecSt3_1_Index, genElecSt3_1_FirstMom;
    
    // Leading Status 1 particle characteristics
    float genBSt1_0_Energy, genBSt1_0_Px, genBSt1_0_Py, genBSt1_0_Pz;
    float genBSt1_0_MomEnergy, genBSt1_0_MomPx, genBSt1_0_MomPy, genBSt1_0_MomPz;
    int   genBSt1_0_PID, genBSt1_0_MomPID, genBSt1_0_MomStatus;
    float genMuonSt1_0_Energy, genMuonSt1_0_Px, genMuonSt1_0_Py, genMuonSt1_0_Pz;
    float genMuonSt1_0_MomEnergy, genMuonSt1_0_MomPx, genMuonSt1_0_MomPy, genMuonSt1_0_MomPz;
    int   genMuonSt1_0_PID, genMuonSt1_0_MomPID, genMuonSt1_0_MomStatus;
    float genElecSt1_0_Energy, genElecSt1_0_Px, genElecSt1_0_Py, genElecSt1_0_Pz;
    float genElecSt1_0_MomEnergy, genElecSt1_0_MomPx, genElecSt1_0_MomPy, genElecSt1_0_MomPz;
    int   genElecSt1_0_PID, genElecSt1_0_MomPID, genElecSt1_0_MomStatus;
    
    // Sub-Leading Status 1 particle characteristics
    float genBSt1_1_Energy, genBSt1_1_Px, genBSt1_1_Py, genBSt1_1_Pz;
    float genBSt1_1_MomEnergy, genBSt1_1_MomPx, genBSt1_1_MomPy, genBSt1_1_MomPz;
    int   genBSt1_1_PID, genBSt1_1_MomPID, genBSt1_1_MomStatus;
    float genMuonSt1_1_Energy, genMuonSt1_1_Px, genMuonSt1_1_Py, genMuonSt1_1_Pz;
    float genMuonSt1_1_MomEnergy, genMuonSt1_1_MomPx, genMuonSt1_1_MomPy, genMuonSt1_1_MomPz;
    int   genMuonSt1_1_PID, genMuonSt1_1_MomPID, genMuonSt1_1_MomStatus;
    float genElecSt1_1_Energy, genElecSt1_1_Px, genElecSt1_1_Py, genElecSt1_1_Pz;
    float genElecSt1_1_MomEnergy, genElecSt1_1_MomPx, genElecSt1_1_MomPy, genElecSt1_1_MomPz;
    int   genElecSt1_1_PID, genElecSt1_1_MomPID, genElecSt1_1_MomStatus;
    
    
    vector<float> * genTopSt3En, * genTopSt3Pt, * genTopSt3Eta, * genTopSt3Phi;
    vector<int> * genTopSt3_i, * genTopSt3_firstMom, * genTopSt3_pdgId;
    genTopSt3En = new vector<float>;
    genTopSt3Eta = new vector<float>;
    genTopSt3Phi = new vector<float>;
    genTopSt3Pt = new vector<float>;
    genTopSt3_i = new vector<int>;
    genTopSt3_firstMom = new vector<int>;
    genTopSt3_pdgId = new vector<int>;
    
    vector<float> * genBSt3En, * genBSt3Pt, * genBSt3Eta, * genBSt3Phi;
    vector<int> * genBSt3_i, * genBSt3_firstMom, * genBSt3_pdgId;
    genBSt3En = new vector<float>;
    genBSt3Eta = new vector<float>;
    genBSt3Phi = new vector<float>;
    genBSt3Pt = new vector<float>;
    genBSt3_i = new vector<int>;
    genBSt3_firstMom = new vector<int>;
    genBSt3_pdgId = new vector<int>; 
    
    vector<float> * genMuonSt3En, * genMuonSt3Pt, * genMuonSt3Eta, * genMuonSt3Phi;
    vector<int> * genMuonSt3_i, * genMuonSt3_firstMom, * genMuonSt3_pdgId;
    genMuonSt3En = new vector<float>;
    genMuonSt3Eta = new vector<float>;
    genMuonSt3Phi = new vector<float>;
    genMuonSt3Pt = new vector<float>;
    genMuonSt3_i = new vector<int>;
    genMuonSt3_firstMom = new vector<int>;
    genMuonSt3_pdgId = new vector<int>;
    
    vector<float> * genElecSt3En, * genElecSt3Pt, * genElecSt3Eta, * genElecSt3Phi;
    vector<int> * genElecSt3_i, * genElecSt3_firstMom, * genElecSt3_pdgId;
    genElecSt3En = new vector<float>;
    genElecSt3Eta = new vector<float>;
    genElecSt3Phi = new vector<float>;
    genElecSt3Pt = new vector<float>;
    genElecSt3_i = new vector<int>;
    genElecSt3_firstMom = new vector<int>;
    genElecSt3_pdgId = new vector<int>;
    
    vector<float> * genBSt1En, * genBSt1Px, * genBSt1Py, * genBSt1Pz;
    vector<float> * genBSt1MomEn, * genBSt1MomPx, * genBSt1MomPy, * genBSt1MomPz;
    vector<int> * genBSt1MomPID, * genBSt1_pdgId, * genBSt1MomStatus;
    genBSt1En = new vector<float>;
    genBSt1Px = new vector<float>;
    genBSt1Py = new vector<float>;
    genBSt1Pz = new vector<float>;
    genBSt1_pdgId = new vector<int>;
    genBSt1MomEn = new vector<float>;
    genBSt1MomPx = new vector<float>;
    genBSt1MomPy = new vector<float>;
    genBSt1MomPz = new vector<float>;
    genBSt1MomPID = new vector<int>;
    genBSt1MomStatus = new vector<int>;
    
    vector<float> * genMuonSt1En, * genMuonSt1Px, * genMuonSt1Py, * genMuonSt1Pz;
    vector<float> * genMuonSt1MomEn, * genMuonSt1MomPx, * genMuonSt1MomPy, * genMuonSt1MomPz;
    vector<int> * genMuonSt1MomPID, * genMuonSt1_pdgId, * genMuonSt1MomStatus;
    genMuonSt1En = new vector<float>;
    genMuonSt1Px = new vector<float>;
    genMuonSt1Py = new vector<float>;
    genMuonSt1Pz = new vector<float>;
    genMuonSt1_pdgId = new vector<int>;
    genMuonSt1MomEn = new vector<float>;
    genMuonSt1MomPx = new vector<float>;
    genMuonSt1MomPy = new vector<float>;
    genMuonSt1MomPz = new vector<float>;
    genMuonSt1MomPID = new vector<int>;
    genMuonSt1MomStatus = new vector<int>;
    
    vector<float> * genElecSt1En, * genElecSt1Px, * genElecSt1Py, * genElecSt1Pz;
    vector<float> * genElecSt1MomEn, * genElecSt1MomPx, * genElecSt1MomPy, * genElecSt1MomPz;
    vector<int> * genElecSt1MomPID, * genElecSt1_pdgId, * genElecSt1MomStatus;
    genElecSt1En = new vector<float>;
    genElecSt1Px = new vector<float>;
    genElecSt1Py = new vector<float>;
    genElecSt1Pz = new vector<float>;
    genElecSt1_pdgId = new vector<int>;
    genElecSt1MomEn = new vector<float>;
    genElecSt1MomPx = new vector<float>;
    genElecSt1MomPy = new vector<float>;
    genElecSt1MomPz = new vector<float>;
    genElecSt1MomPID = new vector<int>;
    genElecSt1MomStatus = new vector<int>;
    
    /*
     vector<float> * genStopEn, * genStopPt, * genStopEta, * genStopPhi;
     vector<int> * genStop_i, * genStop_firstMom, * genStop_pdgId;
     genStopEn = new vector<float>;
     genStopEta = new vector<float>;
     genStopPhi = new vector<float>;
     genStopPt = new vector<float>; 
     
     vector<float> * genChi0En, * genChi0Pt, * genChi0Eta, * genChi0Phi;
     vector<int> * genChi0_i, * genChi0_firstMom, * genChi0_pdgId;
     genChi0En = new vector<float>;
     genChi0Eta = new vector<float>;
     genChi0Phi = new vector<float>;
     genChi0Pt = new vector<float>; 
     
     vector<float> * genChiPMEn, * genChiPMPt, * genChiPMEta, * genChiPMPhi;
     vector<int> * genChiPM_i, * genChiPM_firstMom, * genChiPM_pdgId;
     genChiPMEn = new vector<float>;
     genChiPMEta = new vector<float>;
     genChiPMPhi = new vector<float>;
     genChiPMPt = new vector<float>; 
     
     vector<float> * genBEn, * genBPt, * genBEta, * genBPhi;
     vector<int> * genB_i, * genB_firstMom, * genB_pdgId;
     genBEn = new vector<float>;
     genBEta = new vector<float>;
     genBPhi = new vector<float>;
     genBPt = new vector<float>;  
     
     vector<float> * genWEn, * genWPt, * genWEta, * genWPhi;
     vector<int> * genW_i, * genW_firstMom, * genW_pdgId;
     genWEn = new vector<float>;
     genWEta = new vector<float>;
     genWPhi = new vector<float>;
     genWPt = new vector<float>;
     
     genW_i = new vector<int>;
     genW_firstMom = new vector<int>;
     genW_pdgId = new vector<int>;
     
     vector<float> * genTauLepPx, * genTauLepPy, * genTauLepPz, * genTauLepEn;
     vector<int> * genTauLepPID;
     vector<bool> * genTauLepisLepDecay;
     genTauLepPx = new vector<float>;
     genTauLepPy = new vector<float>;
     genTauLepPz = new vector<float>;
     genTauLepEn = new vector<float>;
     genTauLepPID = new vector<int>;
     genTauLepisLepDecay = new vector<bool>;
     
     
     vector<float> * genTauEn, * genTauPt, * genTauEta, * genTauPhi;
     vector<int> * genTau_i, * genTau_firstMom, * genTau_pdgId;
     genTauEn = new vector<float>;
     genTauEta = new vector<float>;
     genTauPhi = new vector<float>;
     genTauPt = new vector<float>; 
     
     vector<float> * genElecEn, * genElecPt, * genElecEta, * genElecPhi;
     vector<int> * genElec_i, * genElec_firstMom, * genElec_pdgId;
     genElecEn = new vector<float>;
     genElecEta = new vector<float>;
     genElecPhi = new vector<float>;
     genElecPt = new vector<float>; 
     
     vector<float> * genMuonEn, * genMuonPt, * genMuonEta, * genMuonPhi;
     vector<int> * genMuon_i, * genMuon_firstMom, * genMuon_pdgId;
     genMuonEn = new vector<float>;
     genMuonEta = new vector<float>;
     genMuonPhi = new vector<float>;
     genMuonPt = new vector<float>; 
     
     vector<float> * genJetPx, * genJetPy, * genJetPz, * genJetEt, * genJetEta, * genJetEn, * genJetInvisE, 
     vector<bool> * genJetIsGenJet;
     genJetPx = new vector<float>;
     genJetPy = new vector<float>;
     genJetPz = new vector<float>;
     genJetEn = new vector<float>;
     genJetEt = new vector<float>;
     genJetEta = new vector<float>;
     genJetInvisE = new vector<float>;
     genJetIsGenJet = new vector<bool>;
     */
    
    float genMET, genMETPhi;
    
    float SysVar;
    
    TLorentzVector Lep0Vec, Lep1Vec;
    TLorentzVector Jet0Vec, Jet1Vec, BtagJet0Vec, BtagJet1Vec;
    float weight;
    
    float lumi = 19603.691; //5296.3; // ipb                                                                                     
    const float genStopMassMin = 295, genStopMassMax = 355, genDeltaM_stopChi0_Min = 195, genDeltaM_stopChi0_Max = 255; 
    // add 5 GeV safety margin (deltaM = 10 GeV in the FineBin sample)  
    float Nevt_stop_oneMassPoint = 50000 * ( (genStopMassMax-genStopMassMin)/10. ) * ( (genDeltaM_stopChi0_Max-genDeltaM_stopChi0_Min)/10. );  
    // 50k evts per point x Npoints
    ////input cuts/commands    
    //    const double PI = 3.14159265;
    bool doPFElecMu      = 0;
    bool grabOutDir      = 0;      // whether or not to use the file: "outputSavePath.txt" for where to save output
    bool doData          = 0;
    bool doVerbosity     = 0;
    int  whichNTupleType = 0; //0 IFCA Oviedo; 1 DESY
    bool doPURW          = 0;
    bool doHackPURW      = 0;
    bool doPURWOviToDESY = 0;
    bool doSignal        = 0;
    bool printEventNum   = 0;
    bool doMassCut       = 0;
    /////loop over inputs
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"-i",2) == 0) { 
            fInName = TString(argv[k+1]);
        }
        else if (strncmp (argv[k],"-w",2) == 0) {
            whichNTupleType = strtol(argv[k+1], NULL, 10);   
        }
        else if (strncmp (argv[k],"doVerbosity",11) == 0) {
            doVerbosity = 1;   
        }
        else if (strncmp (argv[k],"doPURW",6) == 0) {
            doPURW = 1;           
        }
        else if (strncmp (argv[k],"doHackPURW",10) == 0) {
            doHackPURW = 1;   
        }
        else if (strncmp (argv[k],"doPURWOviToDESY",15) == 0) {
            doPURWOviToDESY = 1;   
        }
        else if (strncmp (argv[k],"gOutDir", 7) == 0) {
            grabOutDir = 1;
        }
        else if (strncmp (argv[k],"isSig", 5) == 0) {
            doSignal = 1;
        }
        else if (strncmp (argv[k],"pEvNum", 6) == 0) {
            printEventNum = 1;
        }
        else if (strncmp (argv[k],"doMassCut", 9) == 0) {
            doMassCut = 1;
            genStopMassCut = strtol(argv[k+1], NULL, 10);   
            genChi0MassCut = strtol(argv[k+2], NULL, 10);   
            genCharginoMassCut = strtol(argv[k+3], NULL, 10);   
        }
    }
    char Buffer[500];
    char MyRootFile[2000];
    ifstream * outDirFile;
    TRegexp fCutSlash("[^/]+$");
    fOutName = "";
    if (grabOutDir) {
        outDirFile = new ifstream(TString("outputSavePath.txt"));
        if (!(outDirFile->eof())) {
            outDirFile->getline(Buffer,500);
            fOutName += TString(string(Buffer));
            fOutName += "/"; //in case user forgot a slash
        }
    }
    fOutName += fInName(fCutSlash);
    
    if (fInName.Contains("MuEG") || fInName.Contains("DoubleMu") || fInName.Contains("DoubleEl")) {
        cout << "Running on Data" << endl;
        doData = 1;
    }
    if (whichNTupleType == 0) {
        fOutName += "_Oviedo";        
    }
    else {
        fOutName += "_DESY";
    }
    if (doPURW && !doData) fOutName += "_PURW";
    if (doPURWOviToDESY && !doData) fOutName += "OviToDESY";
    fOutName += "_SkimOutput.root";
    cout << "saving to " << fOutName << endl;
    TFile * outputFile;
    outputFile = new TFile(fOutName,"RECREATE");
    TTree * outTree = new TTree("OviSkimTree", "OviSkimTree");
    TChain fileTree("Tree");
    TFile inputFile(fInName + TString(".root"));
    //////////////////////////
    fileTree.Add(fInName + TString(".root"));
    if (doPFElecMu) {
        fileTree.SetBranchAddress("T_Elec_PFElecPt", &ElecPt);
        fileTree.SetBranchAddress("T_Elec_PFElecPx", &ElecPx);
        fileTree.SetBranchAddress("T_Elec_PFElecPy", &ElecPy);
        fileTree.SetBranchAddress("T_Elec_PFElecPz", &ElecPz);
        fileTree.SetBranchAddress("T_Elec_PFElecE",  &ElecE);
        
        fileTree.SetBranchAddress("T_Muon_PFMuonPt", &MuonPt);
        fileTree.SetBranchAddress("T_Muon_PFMuonPx", &MuonPx);
        fileTree.SetBranchAddress("T_Muon_PFMuonPy", &MuonPy);
        fileTree.SetBranchAddress("T_Muon_PFMuonPz", &MuonPz);
        fileTree.SetBranchAddress("T_Muon_PFMuonE", &MuonE);
    }
    else {
        fileTree.SetBranchAddress("T_Elec_Pt", &ElecPt);
        fileTree.SetBranchAddress("T_Elec_Px", &ElecPx);
        fileTree.SetBranchAddress("T_Elec_Py", &ElecPy);
        fileTree.SetBranchAddress("T_Elec_Pz", &ElecPz);
        fileTree.SetBranchAddress("T_Elec_Energy", &ElecE);
        
        fileTree.SetBranchAddress("T_Muon_Pt", &MuonPt);
        fileTree.SetBranchAddress("T_Muon_Px", &MuonPx);
        fileTree.SetBranchAddress("T_Muon_Py", &MuonPy);
        fileTree.SetBranchAddress("T_Muon_Pz", &MuonPz);
        fileTree.SetBranchAddress("T_Muon_Energy", &MuonE);                    
    }    
    fileTree.SetBranchAddress("T_Elec_Eta", &ElecEta);
    fileTree.SetBranchAddress("T_Elec_Charge", &ElecCharge);    
    fileTree.SetBranchAddress("T_Elec_chargedHadronIso", &ElecPFCharIso);
    fileTree.SetBranchAddress("T_Elec_neutralHadronIso", &ElecPFNeutIso);
    fileTree.SetBranchAddress("T_Elec_photonIso", &ElecPFPhotIso);
    fileTree.SetBranchAddress("T_Elec_passConversionVeto", &passConvVeto);
    fileTree.SetBranchAddress("T_Elec_isPF", &isPFElectron);
    
    
    fileTree.SetBranchAddress("T_Muon_Eta", &MuonEta);
    fileTree.SetBranchAddress("T_Muon_Charge", &MuonCharge);    
    fileTree.SetBranchAddress("T_Muon_chargedHadronIsoR04", &MuonPFCharIso);
    fileTree.SetBranchAddress("T_Muon_neutralHadronIsoR04", &MuonPFNeutIso);
    fileTree.SetBranchAddress("T_Muon_photonIsoR04", &MuonPFPhotIso);
    fileTree.SetBranchAddress("T_Muon_sumPUPtR04", &MuonSumPUPt);
    fileTree.SetBranchAddress("T_Muon_IsGMPTMuons", &isGMPTMuons);
    fileTree.SetBranchAddress("T_Muon_isPFMuon", &isPFMuon);
    fileTree.SetBranchAddress("T_Muon_IPwrtAveBSInTrack", &MuonD0);
    fileTree.SetBranchAddress("T_Muon_vz", &MuonVertZ);
    
    
    fileTree.SetBranchAddress("T_JetAKCHS_Et", &JetEt);
    fileTree.SetBranchAddress("T_JetAKCHS_Eta", &JetEta);  
    fileTree.SetBranchAddress("T_JetAKCHS_Px", &JetPx);
    fileTree.SetBranchAddress("T_JetAKCHS_Py", &JetPy);
    fileTree.SetBranchAddress("T_JetAKCHS_Pz", &JetPz);
    fileTree.SetBranchAddress("T_JetAKCHS_Energy", &JetE);
    fileTree.SetBranchAddress("T_JetAKCHS_NeutHadEnergyFrac", &JetNHF);
    fileTree.SetBranchAddress("T_JetAKCHS_NeutEmEnergyFrac", &JetNEF);
    fileTree.SetBranchAddress("T_JetAKCHS_CharHadEnergyFrac", &JetCHF);
    fileTree.SetBranchAddress("T_JetAKCHS_CharEmEnergyFrac", &JetCEF);
    fileTree.SetBranchAddress("T_JetAKCHS_Tag_CombSVtx", &JetBTag);
    fileTree.SetBranchAddress("T_JetAKCHS_nDaughters", &JetNDaug);
    fileTree.SetBranchAddress("T_JetAKCHS_ChargedMultiplicity", &JetCharMult);
    fileTree.SetBranchAddress("T_Vertex_z", &VertZ);
    
    fileTree.SetBranchAddress( "T_METPFTypeI_ET",     &MET );
    fileTree.SetBranchAddress( "T_METPFTypeI_Phi", &MET_Phi );
    fileTree.SetBranchAddress( "T_METPFTypeI_Sig",  &METSig );
    
    fileTree.SetBranchAddress( "T_EventF_EcalDeadCell", &filterECalDeadCell );
    fileTree.SetBranchAddress( "T_EventF_logErrorTooManyClusters", &filterLogErrorTooManyClusters );
    fileTree.SetBranchAddress( "T_EventF_trackingFailure", &filterTrackingFailure );
    fileTree.SetBranchAddress( "T_EventF_hcalLaser", &filterHCalLaser );
    fileTree.SetBranchAddress( "T_EventF_ecalLaserCorr", &filterECalLaserCorr );
    fileTree.SetBranchAddress( "T_EventF_toomanystripclus", &filterTooManyStripClus );
    fileTree.SetBranchAddress( "T_EventF_manystripclus", &filterManyStripClus );
    fileTree.SetBranchAddress( "T_EventF_eeBadSc", &filterEEBadSC );
    
    //// Generator information
    fileTree.SetBranchAddress( "T_METgen_ET", &genMET );
    fileTree.SetBranchAddress( "T_METgen_Phi", &genMETPhi );
    /// Status 3 particles
    fileTree.SetBranchAddress( "T_Gen_tSt3_pdgId", &genTopSt3_pdgId );
    fileTree.SetBranchAddress( "T_Gen_tSt3_firstMother", &genTopSt3_firstMom );
    fileTree.SetBranchAddress( "T_Gen_tSt3_i", &genTopSt3_i  );
    fileTree.SetBranchAddress( "T_Gen_tSt3_energy", &genTopSt3En );
    fileTree.SetBranchAddress( "T_Gen_tSt3_pt", &genTopSt3Pt );
    fileTree.SetBranchAddress( "T_Gen_tSt3_eta", &genTopSt3Eta );
    fileTree.SetBranchAddress( "T_Gen_tSt3_phi", &genTopSt3Phi );
    
    fileTree.SetBranchAddress( "T_Gen_bSt3_pdgId", &genBSt3_pdgId );
    fileTree.SetBranchAddress( "T_Gen_bSt3_firstMother", &genBSt3_firstMom );
    fileTree.SetBranchAddress( "T_Gen_bSt3_i", &genBSt3_i  );
    fileTree.SetBranchAddress( "T_Gen_bSt3_energy", &genBSt3En );
    fileTree.SetBranchAddress( "T_Gen_bSt3_pt", &genBSt3Pt );
    fileTree.SetBranchAddress( "T_Gen_bSt3_eta", &genBSt3Eta );
    fileTree.SetBranchAddress( "T_Gen_bSt3_phi", &genBSt3Phi );
    
    fileTree.SetBranchAddress( "T_Gen_MuonSt3_pdgId", &genMuonSt3_pdgId );
    fileTree.SetBranchAddress( "T_Gen_MuonSt3_firstMother", &genMuonSt3_firstMom );
    fileTree.SetBranchAddress( "T_Gen_MuonSt3_i", &genMuonSt3_i  );
    fileTree.SetBranchAddress( "T_Gen_MuonSt3_energy", &genMuonSt3En );
    fileTree.SetBranchAddress( "T_Gen_MuonSt3_pt", &genMuonSt3Pt );
    fileTree.SetBranchAddress( "T_Gen_MuonSt3_eta", &genMuonSt3Eta );
    fileTree.SetBranchAddress( "T_Gen_MuonSt3_phi", &genMuonSt3Phi );
    
    fileTree.SetBranchAddress( "T_Gen_ElecSt3_pdgId", &genElecSt3_pdgId );
    fileTree.SetBranchAddress( "T_Gen_ElecSt3_firstMother", &genElecSt3_firstMom );
    fileTree.SetBranchAddress( "T_Gen_ElecSt3_i", &genElecSt3_i  );
    fileTree.SetBranchAddress( "T_Gen_ElecSt3_energy", &genElecSt3En );
    fileTree.SetBranchAddress( "T_Gen_ElecSt3_pt", &genElecSt3Pt );
    fileTree.SetBranchAddress( "T_Gen_ElecSt3_eta", &genElecSt3Eta );
    fileTree.SetBranchAddress( "T_Gen_ElecSt3_phi", &genElecSt3Phi );
    
    /// Status 1 particles (technically, b-quark is a status 2 -- intermediate particle)
    fileTree.SetBranchAddress( "T_Gen_b_PID", &genBSt1_pdgId );
    fileTree.SetBranchAddress( "T_Gen_b_Px", &genBSt1Px );
    fileTree.SetBranchAddress( "T_Gen_b_Py", &genBSt1Py );
    fileTree.SetBranchAddress( "T_Gen_b_Pz", &genBSt1Pz );
    fileTree.SetBranchAddress( "T_Gen_b_Energy", &genBSt1En );
    fileTree.SetBranchAddress( "T_Gen_b_MPID", &genBSt1MomPID );
    fileTree.SetBranchAddress( "T_Gen_b_MPx", &genBSt1MomPx );
    fileTree.SetBranchAddress( "T_Gen_b_MPy", &genBSt1MomPy );
    fileTree.SetBranchAddress( "T_Gen_b_MPz", &genBSt1MomPz );
    fileTree.SetBranchAddress( "T_Gen_b_MEnergy", &genBSt1MomEn );
    fileTree.SetBranchAddress( "T_Gen_b_MSt", &genBSt1MomStatus );
    
    fileTree.SetBranchAddress( "T_Gen_Muon_PID", &genMuonSt1_pdgId );
    fileTree.SetBranchAddress( "T_Gen_Muon_Px", &genMuonSt1Px );
    fileTree.SetBranchAddress( "T_Gen_Muon_Py", &genMuonSt1Py );
    fileTree.SetBranchAddress( "T_Gen_Muon_Pz", &genMuonSt1Pz );
    fileTree.SetBranchAddress( "T_Gen_Muon_Energy", &genMuonSt1En );
    fileTree.SetBranchAddress( "T_Gen_Muon_MPID", &genMuonSt1MomPID );
    fileTree.SetBranchAddress( "T_Gen_Muon_MPx", &genMuonSt1MomPx );
    fileTree.SetBranchAddress( "T_Gen_Muon_MPy", &genMuonSt1MomPy );
    fileTree.SetBranchAddress( "T_Gen_Muon_MPz", &genMuonSt1MomPz );
    fileTree.SetBranchAddress( "T_Gen_Muon_MEnergy", &genMuonSt1MomEn );
    fileTree.SetBranchAddress( "T_Gen_Muon_MSt", &genMuonSt1MomStatus );
    
    fileTree.SetBranchAddress( "T_Gen_Elec_PID", &genElecSt1_pdgId );
    fileTree.SetBranchAddress( "T_Gen_Elec_Px", &genElecSt1Px );
    fileTree.SetBranchAddress( "T_Gen_Elec_Py", &genElecSt1Py );
    fileTree.SetBranchAddress( "T_Gen_Elec_Pz", &genElecSt1Pz );
    fileTree.SetBranchAddress( "T_Gen_Elec_Energy", &genElecSt1En );
    fileTree.SetBranchAddress( "T_Gen_Elec_MPID", &genElecSt1MomPID );
    fileTree.SetBranchAddress( "T_Gen_Elec_MPx", &genElecSt1MomPx );
    fileTree.SetBranchAddress( "T_Gen_Elec_MPy", &genElecSt1MomPy );
    fileTree.SetBranchAddress( "T_Gen_Elec_MPz", &genElecSt1MomPz );
    fileTree.SetBranchAddress( "T_Gen_Elec_MEnergy", &genElecSt1MomEn );
    fileTree.SetBranchAddress( "T_Gen_Elec_MSt", &genElecSt1MomStatus );
    
    //Signal stuff
    fileTree.SetBranchAddress( "T_Gen_StopMass", &genStopMass );
    fileTree.SetBranchAddress( "T_Gen_Chi0Mass", &genChi0Mass );
    fileTree.SetBranchAddress( "T_Gen_CharginoMass", &genCharginoMass );
    
    //run information
    fileTree.SetBranchAddress( "T_Event_RunNumber", &RunNum );
    fileTree.SetBranchAddress( "T_Event_EventNumber", &EventNum );
    fileTree.SetBranchAddress( "T_Event_LuminosityBlock", &LumiBlock );
    
    
    ///Out tree information
    
    outTree->Branch("TWeight",   &weight);
    outTree->Branch("TChannel",  &Type);
    outTree->Branch("TNPV",      &nVtx);
    outTree->Branch("TRunNum",   &RunNum);
    outTree->Branch("TEventNum", &EventNum);
    outTree->Branch("TLumiBlock",&LumiBlock);
    outTree->Branch("TMET",      &MET);
    outTree->Branch("TMET_Phi",  &MET_Phi);
    outTree->Branch("TMETSig",   &METSig);
    
    outTree->Branch("TLep0Px",   &Lep0Px); 
    outTree->Branch("TLep0Py",   &Lep0Py); 
    outTree->Branch("TLep0Pz",   &Lep0Pz); 
    outTree->Branch("TLep0E",    &Lep0E);
    outTree->Branch("TLep0PdgId",&Lep0PDGID);
    outTree->Branch("TLep1Px",   &Lep1Px); 
    outTree->Branch("TLep1Py",   &Lep1Py); 
    outTree->Branch("TLep1Pz",   &Lep1Pz); 
    outTree->Branch("TLep1E",    &Lep1E);
    outTree->Branch("TLep1PdgId",&Lep0PDGID);
    
    outTree->Branch("TNJets",    &NJets); 
    //    outTree->Branch("THT",     &Jet0Py); 
    outTree->Branch("TNJetsBtag",&NBtagJets); 
    
    outTree->Branch("HT",        &HT);
    outTree->Branch("TJet0Px",   &Jet0Px); 
    outTree->Branch("TJet0Py",   &Jet0Py); 
    outTree->Branch("TJet0Pz",   &Jet0Pz); 
    outTree->Branch("TJet0E",    &Jet0E);
    outTree->Branch("TJet1Px",   &Jet1Px); 
    outTree->Branch("TJet1Py",   &Jet1Py); 
    outTree->Branch("TJet1Pz",   &Jet1Pz); 
    outTree->Branch("TJet1E",    &Jet1E);
    
    outTree->Branch("TBtagJet0Px",   &BtagJet0Px); 
    outTree->Branch("TBtagJet0Py",   &BtagJet0Py); 
    outTree->Branch("TBtagJet0Pz",   &BtagJet0Pz); 
    outTree->Branch("TBtagJet0E",    &BtagJet0E);
    outTree->Branch("TBtagJet0Index",    &BtagJet0Index);    
    outTree->Branch("TBtagJet1Px",   &BtagJet1Px); 
    outTree->Branch("TBtagJet1Py",   &BtagJet1Py); 
    outTree->Branch("TBtagJet1Pz",   &BtagJet1Pz); 
    outTree->Branch("TBtagJet1E",    &BtagJet1E);
    outTree->Branch("TBtagJet1Index",    &BtagJet1Index);
    
    if (doSignal) {
        outTree->Branch("TGenStopMass0", &genStopMass0);
        outTree->Branch("TGenStopMass1", &genStopMass1);
        outTree->Branch("TGenChi0Mass0", &genChi0Mass0);
        outTree->Branch("TGenChi0Mass1", &genChi0Mass1);
        outTree->Branch("TGenCharginoMass0", &genCharginoMass0);
        outTree->Branch("TGenCharginoMass1", &genCharginoMass1);
    }
    
    if (!doData) {
        outTree->Branch("TGenMET", &genMET);
        outTree->Branch("TGenMETPhi", &genMETPhi);
        
        /// Status 3 lead particle branches
        outTree->Branch("TGenTopSt3_0_PDGID", &genTopSt3_0_PID );
        outTree->Branch("TGenTopSt3_0_FirstMom", &genTopSt3_0_FirstMom );
        outTree->Branch("TGenTopSt3_0_Index", &genTopSt3_0_Index );
        outTree->Branch("TGenTopSt3_0_Energy", &genTopSt3_0_Energy );
        outTree->Branch("TGenTopSt3_0_Pt", &genTopSt3_0_Pt );
        outTree->Branch("TGenTopSt3_0_Eta", &genTopSt3_0_Eta );
        outTree->Branch("TGenTopSt3_0_Phi", &genTopSt3_0_Phi );
        
        outTree->Branch("TGenBSt3_0_PDGID", &genBSt3_0_PID );
        outTree->Branch("TGenBSt3_0_FirstMom", &genBSt3_0_FirstMom );
        outTree->Branch("TGenBSt3_0_Index", &genBSt3_0_Index );
        outTree->Branch("TGenBSt3_0_Energy", &genBSt3_0_Energy );
        outTree->Branch("TGenBSt3_0_Pt", &genBSt3_0_Pt );
        outTree->Branch("TGenBSt3_0_Eta", &genBSt3_0_Eta );
        outTree->Branch("TGenBSt3_0_Phi", &genBSt3_0_Phi );
        
        outTree->Branch("TGenMuonSt3_0_PDGID", &genMuonSt3_0_PID );
        outTree->Branch("TGenMuonSt3_0_FirstMom", &genMuonSt3_0_FirstMom );
        outTree->Branch("TGenMuonSt3_0_Index", &genMuonSt3_0_Index );
        outTree->Branch("TGenMuonSt3_0_Energy", &genMuonSt3_0_Energy );
        outTree->Branch("TGenMuonSt3_0_Pt", &genMuonSt3_0_Pt );
        outTree->Branch("TGenMuonSt3_0_Eta", &genMuonSt3_0_Eta );
        outTree->Branch("TGenMuonSt3_0_Phi", &genMuonSt3_0_Phi );
        
        outTree->Branch("TGenElecSt3_0_PDGID", &genElecSt3_0_PID );
        outTree->Branch("TGenElecSt3_0_FirstMom", &genElecSt3_0_FirstMom );
        outTree->Branch("TGenElecSt3_0_Index", &genElecSt3_0_Index );
        outTree->Branch("TGenElecSt3_0_Energy", &genElecSt3_0_Energy );
        outTree->Branch("TGenElecSt3_0_Pt", &genElecSt3_0_Pt );
        outTree->Branch("TGenElecSt3_0_Eta", &genElecSt3_0_Eta );
        outTree->Branch("TGenElecSt3_0_Phi", &genElecSt3_0_Phi );
        
        /// Status 3 sub-lead particle branches
        outTree->Branch("TGenTopSt3_1_PDGID", &genTopSt3_1_PID );
        outTree->Branch("TGenTopSt3_1_FirstMom", &genTopSt3_1_FirstMom );
        outTree->Branch("TGenTopSt3_1_Index", &genTopSt3_1_Index );
        outTree->Branch("TGenTopSt3_1_Energy", &genTopSt3_1_Energy );
        outTree->Branch("TGenTopSt3_1_Pt", &genTopSt3_1_Pt );
        outTree->Branch("TGenTopSt3_1_Eta", &genTopSt3_1_Eta );
        outTree->Branch("TGenTopSt3_1_Phi", &genTopSt3_1_Phi );
        
        outTree->Branch("TGenBSt3_1_PDGID", &genBSt3_1_PID );
        outTree->Branch("TGenBSt3_1_FirstMom", &genBSt3_1_FirstMom );
        outTree->Branch("TGenBSt3_1_Index", &genBSt3_1_Index );
        outTree->Branch("TGenBSt3_1_Energy", &genBSt3_1_Energy );
        outTree->Branch("TGenBSt3_1_Pt", &genBSt3_1_Pt );
        outTree->Branch("TGenBSt3_1_Eta", &genBSt3_1_Eta );
        outTree->Branch("TGenBSt3_1_Phi", &genBSt3_1_Phi );
        
        outTree->Branch("TGenMuonSt3_1_PDGID", &genMuonSt3_1_PID );
        outTree->Branch("TGenMuonSt3_1_FirstMom", &genMuonSt3_1_FirstMom );
        outTree->Branch("TGenMuonSt3_1_Index", &genMuonSt3_1_Index );
        outTree->Branch("TGenMuonSt3_1_Energy", &genMuonSt3_1_Energy );
        outTree->Branch("TGenMuonSt3_1_Pt", &genMuonSt3_1_Pt );
        outTree->Branch("TGenMuonSt3_1_Eta", &genMuonSt3_1_Eta );
        outTree->Branch("TGenMuonSt3_1_Phi", &genMuonSt3_1_Phi );
        
        outTree->Branch("TGenElecSt3_1_PDGID", &genElecSt3_1_PID );
        outTree->Branch("TGenElecSt3_1_FirstMom", &genElecSt3_1_FirstMom );
        outTree->Branch("TGenElecSt3_1_Index", &genElecSt3_1_Index );
        outTree->Branch("TGenElecSt3_1_Energy", &genElecSt3_1_Energy );
        outTree->Branch("TGenElecSt3_1_Pt", &genElecSt3_1_Pt );
        outTree->Branch("TGenElecSt3_1_Eta", &genElecSt3_1_Eta );
        outTree->Branch("TGenElecSt3_1_Phi", &genElecSt3_1_Phi );
        
        /// Status 1 lead particle branches
        outTree->Branch( "T_Gen_BSt1_0_PID", &genBSt1_0_PID );
        outTree->Branch( "T_Gen_BSt1_0_Px", &genBSt1_0_Px );
        outTree->Branch( "T_Gen_BSt1_0_Py", &genBSt1_0_Py );
        outTree->Branch( "T_Gen_BSt1_0_Pz", &genBSt1_0_Pz );
        outTree->Branch( "T_Gen_BSt1_0_Energy", &genBSt1_0_Energy );
        outTree->Branch( "T_Gen_BSt1_0_MomPID", &genBSt1_0_MomPID );
        outTree->Branch( "T_Gen_BSt1_0_MomPx", &genBSt1_0_MomPx );
        outTree->Branch( "T_Gen_BSt1_0_MomPy", &genBSt1_0_MomPy );
        outTree->Branch( "T_Gen_BSt1_0_MomPz", &genBSt1_0_MomPz );
        outTree->Branch( "T_Gen_BSt1_0_MomEnergy", &genBSt1_0_MomEnergy );
        outTree->Branch( "T_Gen_BSt1_0_MomStat", &genBSt1_0_MomStatus );
        
        outTree->Branch( "T_Gen_MuonSt1_0_PID", &genMuonSt1_0_PID );
        outTree->Branch( "T_Gen_MuonSt1_0_Px", &genMuonSt1_0_Px );
        outTree->Branch( "T_Gen_MuonSt1_0_Py", &genMuonSt1_0_Py );
        outTree->Branch( "T_Gen_MuonSt1_0_Pz", &genMuonSt1_0_Pz );
        outTree->Branch( "T_Gen_MuonSt1_0_Energy", &genMuonSt1_0_Energy );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPID", &genMuonSt1_0_MomPID );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPx", &genMuonSt1_0_MomPx );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPy", &genMuonSt1_0_MomPy );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPz", &genMuonSt1_0_MomPz );
        outTree->Branch( "T_Gen_MuonSt1_0_MomEnergy", &genMuonSt1_0_MomEnergy );
        outTree->Branch( "T_Gen_MuonSt1_0_MomStat", &genMuonSt1_0_MomStatus );
        
        outTree->Branch( "T_Gen_ElecSt1_0_PID", &genElecSt1_0_PID );
        outTree->Branch( "T_Gen_ElecSt1_0_Px", &genElecSt1_0_Px );
        outTree->Branch( "T_Gen_ElecSt1_0_Py", &genElecSt1_0_Py );
        outTree->Branch( "T_Gen_ElecSt1_0_Pz", &genElecSt1_0_Pz );
        outTree->Branch( "T_Gen_ElecSt1_0_Energy", &genElecSt1_0_Energy );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPID", &genElecSt1_0_MomPID );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPx", &genElecSt1_0_MomPx );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPy", &genElecSt1_0_MomPy );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPz", &genElecSt1_0_MomPz );
        outTree->Branch( "T_Gen_ElecSt1_0_MomEnergy", &genElecSt1_0_MomEnergy );
        outTree->Branch( "T_Gen_ElecSt1_0_MomStat", &genElecSt1_0_MomStatus );  
        
        /// Status 1 sub-lead particle branches
        outTree->Branch( "T_Gen_BSt1_1_PID", &genBSt1_1_PID );
        outTree->Branch( "T_Gen_BSt1_1_Px", &genBSt1_1_Px );
        outTree->Branch( "T_Gen_BSt1_1_Py", &genBSt1_1_Py );
        outTree->Branch( "T_Gen_BSt1_1_Pz", &genBSt1_1_Pz );
        outTree->Branch( "T_Gen_BSt1_1_Energy", &genBSt1_1_Energy );
        outTree->Branch( "T_Gen_BSt1_1_MomPID", &genBSt1_1_MomPID );
        outTree->Branch( "T_Gen_BSt1_1_MomPx", &genBSt1_1_MomPx );
        outTree->Branch( "T_Gen_BSt1_1_MomPy", &genBSt1_1_MomPy );
        outTree->Branch( "T_Gen_BSt1_1_MomPz", &genBSt1_1_MomPz );
        outTree->Branch( "T_Gen_BSt1_1_MomEnergy", &genBSt1_1_MomEnergy );
        outTree->Branch( "T_Gen_BSt1_1_MomStat", &genBSt1_1_MomStatus );
        
        outTree->Branch( "T_Gen_MuonSt1_1_PID", &genMuonSt1_1_PID );
        outTree->Branch( "T_Gen_MuonSt1_1_Px", &genMuonSt1_1_Px );
        outTree->Branch( "T_Gen_MuonSt1_1_Py", &genMuonSt1_1_Py );
        outTree->Branch( "T_Gen_MuonSt1_1_Pz", &genMuonSt1_1_Pz );
        outTree->Branch( "T_Gen_MuonSt1_1_Energy", &genMuonSt1_1_Energy );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPID", &genMuonSt1_1_MomPID );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPx", &genMuonSt1_1_MomPx );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPy", &genMuonSt1_1_MomPy );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPz", &genMuonSt1_1_MomPz );
        outTree->Branch( "T_Gen_MuonSt1_1_MomEnergy", &genMuonSt1_1_MomEnergy );
        outTree->Branch( "T_Gen_MuonSt1_1_MomStat", &genMuonSt1_1_MomStatus );
        
        outTree->Branch( "T_Gen_ElecSt1_1_PID", &genElecSt1_1_PID );
        outTree->Branch( "T_Gen_ElecSt1_1_Px", &genElecSt1_1_Px );
        outTree->Branch( "T_Gen_ElecSt1_1_Py", &genElecSt1_1_Py );
        outTree->Branch( "T_Gen_ElecSt1_1_Pz", &genElecSt1_1_Pz );
        outTree->Branch( "T_Gen_ElecSt1_1_Energy", &genElecSt1_1_Energy );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPID", &genElecSt1_1_MomPID );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPx", &genElecSt1_1_MomPx );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPy", &genElecSt1_1_MomPy );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPz", &genElecSt1_1_MomPz );
        outTree->Branch( "T_Gen_ElecSt1_1_MomEnergy", &genElecSt1_1_MomEnergy );
        outTree->Branch( "T_Gen_ElecSt1_1_MomStat", &genElecSt1_1_MomStatus );
    }
    
    cout << "--- Processing: " << fileTree.GetEntries() << " events" << endl;
    h_eventCount->Fill(1);
    h_eventCount->SetEntries(fileTree.GetEntries());
    outputFile->cd();
    //    int lep0Index, lep1Index;
    bool doEvent;
    // cout << "test2 " << endl;
    
    
    vector<TLorentzVector> * Jets = new vector<TLorentzVector>;
    vector<int> * BJetIndices = new vector<int>;
    vector<int> * ElecIndices = new vector<int>;
    vector<int> * MuonIndices = new vector<int>;
    vector<int> * pdgId = new vector<int>;
    vector<TLorentzVector> * Leptons = new vector<TLorentzVector>;
    TLorentzVector patsyVec;
    /////Iterate over events    
    
    
    float roundNum = 1.0;
    int roundMult = 1;
    if (doSignal) {
        roundNum = (fInName.Contains("to") || fInName.Contains("FineBin")) ? 10.0 : 1.0;
        roundMult = (fInName.Contains("to") || fInName.Contains("FineBin")) ? 10 : 1;
    }
    
    for (Long64_t ievt=0; ievt<fileTree.GetEntries();ievt++) {
        //    for (Long64_t ievt=0; ievt<100;ievt++) {
        Leptons = new vector<TLorentzVector>;
        pdgId = new vector<int>;
        Jets = new vector<TLorentzVector>;
        BJetIndices = new vector<int>;
        ElecIndices = new vector<int>;
        MuonIndices = new vector<int>;
        //        lep0Index = 0; lep1Index = 1;
        doEvent = true;
        map<string, float> stringKeyToVar;
        fileTree.GetEntry(ievt);
        genStopMass0 = -1.;
        genStopMass1 = -1.;
        genChi0Mass0 = -1.;
        genChi0Mass1 = -1.;
        genCharginoMass0 = -1.;
        genCharginoMass1 = -1.;
        if (doSignal) {
            if (genStopMass->size() > 1) {
                genStopMass0 = TMath::Nint(genStopMass->at(0)/roundNum) * roundMult;
                genStopMass1 = TMath::Nint(genStopMass->at(1)/roundNum) * roundMult;
            }
            if (genChi0Mass->size() > 1) {
                genChi0Mass0 = TMath::Nint(genChi0Mass->at(0)/roundNum) * roundMult;
                genChi0Mass1 = TMath::Nint(genChi0Mass->at(1)/roundNum) * roundMult;
            }
            if (genCharginoMass->size() > 1) {
                genCharginoMass0 = TMath::Nint(genCharginoMass->at(0)/roundNum) * roundMult;
                genCharginoMass0 = TMath::Nint(genCharginoMass->at(1)/roundNum) * roundMult;
            }
            if (doMassCut) {
                if (genStopMass0 >=0 && abs(genStopMassCut - genStopMass0) < 2.5) continue;
                if (genChi0Mass0 >=0 && abs(genChi0MassCut - genChi0Mass0) < 2.5) continue;
                if (genCharginoMass0 >=0 && abs(genCharginoMassCut - genCharginoMass0) < 2.5) continue;
            }
        }
        if (printEventNum) {
            cout << "in format EventNum:genStopMass:genChi0Mass:genCharginoMass ";
            cout << EventNum << ":" << genStopMass0 << ":" << genChi0Mass0 << ":" <<  genCharginoMass0 << endl;
        }
        /// Status 3 lead particle info
        genTopSt3_0_Energy = ((int) genTopSt3En->size() > 0) ? genTopSt3En->at(0) : -1;
        genTopSt3_0_Pt = ((int) genTopSt3Pt->size() > 0) ? genTopSt3Pt->at(0) : -1;
        genTopSt3_0_Eta = ((int) genTopSt3Eta->size() > 0) ? genTopSt3Eta->at(0) : -1;
        genTopSt3_0_Phi = ((int) genTopSt3Phi->size() > 0) ? genTopSt3Phi->at(0) : -1;
        genTopSt3_0_PID = ((int) genTopSt3_pdgId->size() > 0) ? genTopSt3_pdgId->at(0) : -1;
        genTopSt3_0_Index = ((int) genTopSt3_i->size() > 0) ? genTopSt3_i->at(0) : -1;
        genTopSt3_0_FirstMom = ((int) genTopSt3_firstMom->size() > 0) ? genTopSt3_firstMom->at(0) : -1;
        
        genBSt3_0_Energy = ((int) genBSt3En->size() > 0) ? genBSt3En->at(0) : -1;
        genBSt3_0_Pt = ((int) genBSt3Pt->size() > 0) ? genBSt3Pt->at(0) : -1;
        genBSt3_0_Eta = ((int) genBSt3Eta->size() > 0) ? genBSt3Eta->at(0) : -1;
        genBSt3_0_Phi = ((int) genBSt3Phi->size() > 0) ? genBSt3Phi->at(0) : -1;
        genBSt3_0_PID = ((int) genBSt3_pdgId->size() > 0) ? genBSt3_pdgId->at(0) : -1;
        genBSt3_0_Index = ((int) genBSt3_i->size() > 0) ? genBSt3_i->at(0) : -1;
        genBSt3_0_FirstMom = ((int) genBSt3_firstMom->size() > 0) ? genBSt3_firstMom->at(0) : -1;
        
        genMuonSt3_0_Energy = ((int) genMuonSt3En->size() > 0) ? genMuonSt3En->at(0) : -1;
        genMuonSt3_0_Pt = ((int) genMuonSt3Pt->size() > 0) ? genMuonSt3Pt->at(0) : -1;
        genMuonSt3_0_Eta = ((int) genMuonSt3Eta->size() > 0) ? genMuonSt3Eta->at(0) : -1;
        genMuonSt3_0_Phi = ((int) genMuonSt3Phi->size() > 0) ? genMuonSt3Phi->at(0) : -1;
        genMuonSt3_0_PID = ((int) genMuonSt3_pdgId->size() > 0) ? genMuonSt3_pdgId->at(0) : -1;
        genMuonSt3_0_Index = ((int) genMuonSt3_i->size() > 0) ? genMuonSt3_i->at(0) : -1;
        genMuonSt3_0_FirstMom = ((int) genMuonSt3_firstMom->size() > 0) ? genMuonSt3_firstMom->at(0) : -1;
        
        genElecSt3_0_Energy = ((int) genElecSt3En->size() > 0) ? genElecSt3En->at(0) : -1;
        genElecSt3_0_Pt = ((int) genElecSt3Pt->size() > 0) ? genElecSt3Pt->at(0) : -1;
        genElecSt3_0_Eta = ((int) genElecSt3Eta->size() > 0) ? genElecSt3Eta->at(0) : -1;
        genElecSt3_0_Phi = ((int) genElecSt3Phi->size() > 0) ? genElecSt3Phi->at(0) : -1;
        genElecSt3_0_PID = ((int) genElecSt3_pdgId->size() > 0) ? genElecSt3_pdgId->at(0) : -1;
        genElecSt3_0_Index = ((int) genElecSt3_i->size() > 0) ? genElecSt3_i->at(0) : -1;
        genElecSt3_0_FirstMom = ((int) genElecSt3_firstMom->size() > 0) ? genElecSt3_firstMom->at(0) : -1;
        
        /// Status 3 sub-lead particle info
        genTopSt3_1_Energy = ((int) genTopSt3En->size() > 1) ? genTopSt3En->at(1) : -1;
        genTopSt3_1_Pt = ((int) genTopSt3Pt->size() > 1) ? genTopSt3Pt->at(1) : -1;
        genTopSt3_1_Eta = ((int) genTopSt3Eta->size() > 1) ? genTopSt3Eta->at(1) : -1;
        genTopSt3_1_Phi = ((int) genTopSt3Phi->size() > 1) ? genTopSt3Phi->at(1) : -1;
        genTopSt3_1_PID = ((int) genTopSt3_pdgId->size() > 1) ? genTopSt3_pdgId->at(1) : -1;
        genTopSt3_1_Index = ((int) genTopSt3_i->size() > 1) ? genTopSt3_i->at(1) : -1;
        genTopSt3_1_FirstMom = ((int) genTopSt3_firstMom->size() > 1) ? genTopSt3_firstMom->at(1) : -1;
        
        genBSt3_1_Energy = ((int) genBSt3En->size() > 1) ? genBSt3En->at(1) : -1;
        genBSt3_1_Pt = ((int) genBSt3Pt->size() > 1) ? genBSt3Pt->at(1) : -1;
        genBSt3_1_Eta = ((int) genBSt3Eta->size() > 1) ? genBSt3Eta->at(1) : -1;
        genBSt3_1_Phi = ((int) genBSt3Phi->size() > 1) ? genBSt3Phi->at(1) : -1;
        genBSt3_1_PID = ((int) genBSt3_pdgId->size() > 1) ? genBSt3_pdgId->at(1) : -1;
        genBSt3_1_Index = ((int) genBSt3_i->size() > 1) ? genBSt3_i->at(1) : -1;
        genBSt3_1_FirstMom = ((int) genBSt3_firstMom->size() > 1) ? genBSt3_firstMom->at(1) : -1;
        
        genMuonSt3_1_Energy = ((int) genMuonSt3En->size() > 1) ? genMuonSt3En->at(1) : -1;
        genMuonSt3_1_Pt = ((int) genMuonSt3Pt->size() > 1) ? genMuonSt3Pt->at(1) : -1;
        genMuonSt3_1_Eta = ((int) genMuonSt3Eta->size() > 1) ? genMuonSt3Eta->at(1) : -1;
        genMuonSt3_1_Phi = ((int) genMuonSt3Phi->size() > 1) ? genMuonSt3Phi->at(1) : -1;
        genMuonSt3_1_PID = ((int) genMuonSt3_pdgId->size() > 1) ? genMuonSt3_pdgId->at(1) : -1;
        genMuonSt3_1_Index = ((int) genMuonSt3_i->size() > 1) ? genMuonSt3_i->at(1) : -1;
        genMuonSt3_1_FirstMom = ((int) genMuonSt3_firstMom->size() > 1) ? genMuonSt3_firstMom->at(1) : -1;
        
        genElecSt3_1_Energy = ((int) genElecSt3En->size() > 1) ? genElecSt3En->at(1) : -1;
        genElecSt3_1_Pt = ((int) genElecSt3Pt->size() > 1) ? genElecSt3Pt->at(1) : -1;
        genElecSt3_1_Eta = ((int) genElecSt3Eta->size() > 1) ? genElecSt3Eta->at(1) : -1;
        genElecSt3_1_Phi = ((int) genElecSt3Phi->size() > 1) ? genElecSt3Phi->at(1) : -1;
        genElecSt3_1_PID = ((int) genElecSt3_pdgId->size() > 1) ? genElecSt3_pdgId->at(1) : -1;
        genElecSt3_1_Index = ((int) genElecSt3_i->size() > 1) ? genElecSt3_i->at(1) : -1;
        genElecSt3_1_FirstMom = ((int) genElecSt3_firstMom->size() > 1) ? genElecSt3_firstMom->at(1) : -1;
        
        /// Status 1 lead-particle info
        genBSt1_0_Energy = ((int) genBSt1En->size() > 0) ? genBSt1En->at(0) : -1;
        genBSt1_0_MomEnergy = ((int) genBSt1MomEn->size() > 0) ? genBSt1MomEn->at(0) : -1;
        genBSt1_0_Px = ((int) genBSt1Px->size() > 0) ? genBSt1Px->at(0) : -1;
        genBSt1_0_MomPx = ((int) genBSt1MomPx->size() > 0) ? genBSt1MomPx->at(0) : -1;
        genBSt1_0_Py = ((int) genBSt1Py->size() > 0) ? genBSt1Py->at(0) : -1;
        genBSt1_0_MomPy = ((int) genBSt1MomPy->size() > 0) ? genBSt1MomPy->at(0) : -1;
        genBSt1_0_Pz = ((int) genBSt1Pz->size() > 0) ? genBSt1Pz->at(0) : -1;
        genBSt1_0_MomPz = ((int) genBSt1MomPz->size() > 0) ? genBSt1MomPz->at(0) : -1;
        genBSt1_0_PID = ((int) genBSt1_pdgId->size() > 0) ? genBSt1_pdgId->at(0) : -1;
        genBSt1_0_MomPID = ((int) genBSt1MomPID->size() > 0) ? genBSt1MomPID->at(0) : -1;
        genBSt1_0_MomStatus = ((int) genBSt1MomStatus->size() > 0) ? genBSt1MomStatus->at(0) : -1;
        
        genMuonSt1_0_Energy = ((int) genMuonSt1En->size() > 0) ? genMuonSt1En->at(0) : -1;
        genMuonSt1_0_MomEnergy = ((int) genMuonSt1MomEn->size() > 0) ? genMuonSt1MomEn->at(0) : -1;
        genMuonSt1_0_Px = ((int) genMuonSt1Px->size() > 0) ? genMuonSt1Px->at(0) : -1;
        genMuonSt1_0_MomPx = ((int) genMuonSt1MomPx->size() > 0) ? genMuonSt1MomPx->at(0) : -1;
        genMuonSt1_0_Py = ((int) genMuonSt1Py->size() > 0) ? genMuonSt1Py->at(0) : -1;
        genMuonSt1_0_MomPy = ((int) genMuonSt1MomPy->size() > 0) ? genMuonSt1MomPy->at(0) : -1;
        genMuonSt1_0_Pz = ((int) genMuonSt1Pz->size() > 0) ? genMuonSt1Pz->at(0) : -1;
        genMuonSt1_0_MomPz = ((int) genMuonSt1MomPz->size() > 0) ? genMuonSt1MomPz->at(0) : -1;
        genMuonSt1_0_PID = ((int) genMuonSt1_pdgId->size() > 0) ? genMuonSt1_pdgId->at(0) : -1;
        genMuonSt1_0_MomPID = ((int) genMuonSt1MomPID->size() > 0) ? genMuonSt1MomPID->at(0) : -1;
        genMuonSt1_0_MomStatus = ((int) genMuonSt1MomStatus->size() > 0) ? genMuonSt1MomStatus->at(0) : -1;
        
        genElecSt1_0_Energy = ((int) genElecSt1En->size() > 0) ? genElecSt1En->at(0) : -1;
        genElecSt1_0_MomEnergy = ((int) genElecSt1MomEn->size() > 0) ? genElecSt1MomEn->at(0) : -1;
        genElecSt1_0_Px = ((int) genElecSt1Px->size() > 0) ? genElecSt1Px->at(0) : -1;
        genElecSt1_0_MomPx = ((int) genElecSt1MomPx->size() > 0) ? genElecSt1MomPx->at(0) : -1;
        genElecSt1_0_Py = ((int) genElecSt1Py->size() > 0) ? genElecSt1Py->at(0) : -1;
        genElecSt1_0_MomPy = ((int) genElecSt1MomPy->size() > 0) ? genElecSt1MomPy->at(0) : -1;
        genElecSt1_0_Pz = ((int) genElecSt1Pz->size() > 0) ? genElecSt1Pz->at(0) : -1;
        genElecSt1_0_MomPz = ((int) genElecSt1MomPz->size() > 0) ? genElecSt1MomPz->at(0) : -1;
        genElecSt1_0_PID = ((int) genElecSt1_pdgId->size() > 0) ? genElecSt1_pdgId->at(0) : -1;
        genElecSt1_0_MomPID = ((int) genElecSt1MomPID->size() > 0) ? genElecSt1MomPID->at(0) : -1;
        genElecSt1_0_MomStatus = ((int) genElecSt1MomStatus->size() > 0) ? genElecSt1MomStatus->at(0) : -1; 
        
        /// Status 1 sub-lead-particle info
        genBSt1_1_Energy = ((int) genBSt1En->size() > 1) ? genBSt1En->at(1) : -1;
        genBSt1_1_MomEnergy = ((int) genBSt1MomEn->size() > 1) ? genBSt1MomEn->at(1) : -1;
        genBSt1_1_Px = ((int) genBSt1Px->size() > 1) ? genBSt1Px->at(1) : -1;
        genBSt1_1_MomPx = ((int) genBSt1MomPx->size() > 1) ? genBSt1MomPx->at(1) : -1;
        genBSt1_1_Py = ((int) genBSt1Py->size() > 1) ? genBSt1Py->at(1) : -1;
        genBSt1_1_MomPy = ((int) genBSt1MomPy->size() > 1) ? genBSt1MomPy->at(1) : -1;
        genBSt1_1_Pz = ((int) genBSt1Pz->size() > 1) ? genBSt1Pz->at(1) : -1;
        genBSt1_1_MomPz = ((int) genBSt1MomPz->size() > 1) ? genBSt1MomPz->at(1) : -1;
        genBSt1_1_PID = ((int) genBSt1_pdgId->size() > 1) ? genBSt1_pdgId->at(1) : -1;
        genBSt1_1_MomPID = ((int) genBSt1MomPID->size() > 1) ? genBSt1MomPID->at(1) : -1;
        genBSt1_1_MomStatus = ((int) genBSt1MomStatus->size() > 1) ? genBSt1MomStatus->at(1) : -1;
        
        genMuonSt1_1_Energy = ((int) genMuonSt1En->size() > 1) ? genMuonSt1En->at(1) : -1;
        genMuonSt1_1_MomEnergy = ((int) genMuonSt1MomEn->size() > 1) ? genMuonSt1MomEn->at(1) : -1;
        genMuonSt1_1_Px = ((int) genMuonSt1Px->size() > 1) ? genMuonSt1Px->at(1) : -1;
        genMuonSt1_1_MomPx = ((int) genMuonSt1MomPx->size() > 1) ? genMuonSt1MomPx->at(1) : -1;
        genMuonSt1_1_Py = ((int) genMuonSt1Py->size() > 1) ? genMuonSt1Py->at(1) : -1;
        genMuonSt1_1_MomPy = ((int) genMuonSt1MomPy->size() > 1) ? genMuonSt1MomPy->at(1) : -1;
        genMuonSt1_1_Pz = ((int) genMuonSt1Pz->size() > 1) ? genMuonSt1Pz->at(1) : -1;
        genMuonSt1_1_MomPz = ((int) genMuonSt1MomPz->size() > 1) ? genMuonSt1MomPz->at(1) : -1;
        genMuonSt1_1_PID = ((int) genMuonSt1_pdgId->size() > 1) ? genMuonSt1_pdgId->at(1) : -1;
        genMuonSt1_1_MomPID = ((int) genMuonSt1MomPID->size() > 1) ? genMuonSt1MomPID->at(1) : -1;
        genMuonSt1_1_MomStatus = ((int) genMuonSt1MomStatus->size() > 1) ? genMuonSt1MomStatus->at(1) : -1;
        
        genElecSt1_1_Energy = ((int) genElecSt1En->size() > 1) ? genElecSt1En->at(1) : -1;
        genElecSt1_1_MomEnergy = ((int) genElecSt1MomEn->size() > 1) ? genElecSt1MomEn->at(1) : -1;
        genElecSt1_1_Px = ((int) genElecSt1Px->size() > 1) ? genElecSt1Px->at(1) : -1;
        genElecSt1_1_MomPx = ((int) genElecSt1MomPx->size() > 1) ? genElecSt1MomPx->at(1) : -1;
        genElecSt1_1_Py = ((int) genElecSt1Py->size() > 1) ? genElecSt1Py->at(1) : -1;
        genElecSt1_1_MomPy = ((int) genElecSt1MomPy->size() > 1) ? genElecSt1MomPy->at(1) : -1;
        genElecSt1_1_Pz = ((int) genElecSt1Pz->size() > 1) ? genElecSt1Pz->at(1) : -1;
        genElecSt1_1_MomPz = ((int) genElecSt1MomPz->size() > 1) ? genElecSt1MomPz->at(1) : -1;
        genElecSt1_1_PID = ((int) genElecSt1_pdgId->size() > 1) ? genElecSt1_pdgId->at(1) : -1;
        genElecSt1_1_MomPID = ((int) genElecSt1MomPID->size() > 1) ? genElecSt1MomPID->at(1) : -1;
        genElecSt1_1_MomStatus = ((int) genElecSt1MomStatus->size() > 1) ? genElecSt1MomStatus->at(1) : -1;
        
        weight = 1.;
        float preNVtxRWweight = 1;  
        
        if (doData) {
            regularMETFilterVec = new vector<bool>;
            regularMETFilterVec->push_back(filterECalDeadCell);
            regularMETFilterVec->push_back(filterTrackingFailure);
            regularMETFilterVec->push_back(filterHCalLaser);
            regularMETFilterVec->push_back(filterECalLaserCorr);
            regularMETFilterVec->push_back(filterEEBadSC);
            
            oppositeMETFilterVec = new vector<bool>;
            oppositeMETFilterVec->push_back(filterLogErrorTooManyClusters);
            oppositeMETFilterVec->push_back(filterTooManyStripClus);
            oppositeMETFilterVec->push_back(filterManyStripClus);
            doEvent = FilterMET(regularMETFilterVec, oppositeMETFilterVec);
        } 
        if (!doEvent) continue;
        ElecIndices = ElectronPickOvi(ElecPt, ElecEta, ElecCharge, ElecPFNeutIso, ElecPFCharIso, ElecPFPhotIso, passConvVeto, isPFElectron);
        if (doVerbosity) {
            cout << "iEvent " << ievt << endl;
            cout << "E0 Index " << ElecIndices->at(0) << endl;
            cout << "E1 Index " << ElecIndices->at(1) << endl;
        }
        MuonIndices = MuonPickOvi(MuonPt, MuonEta, MuonCharge, MuonD0, MuonVertZ, VertZ, MuonPFNeutIso, MuonPFCharIso, MuonPFPhotIso, MuonSumPUPt, isGMPTMuons, isPFMuon);
        if (doVerbosity) {
            cout << "Mu0 Index " << MuonIndices->at(0) << endl;
            cout << "Mu1 Index " << MuonIndices->at(1) << endl;
        }
        int E0I = ElecIndices->at(0);
        int E1I = ElecIndices->at(1);
        int M0I = MuonIndices->at(0);
        int M1I = MuonIndices->at(1);
        if (!(M0I < 0)) {
            //            cout << "muon E 0 " << MuonE->at(MuonIndices->at(0)) << endl;
            //            cout << "muon Px " << MuonPx->at(MuonIndices->at(0)) << endl;
            //            cout << "muon Py " << MuonPy->at(MuonIndices->at(0)) << endl;
            //            cout << "muon Pz " << MuonPz->at(MuonIndices->at(0)) << endl;
            patsyVec.SetPxPyPzE(MuonPx->at(MuonIndices->at(0)), MuonPy->at(MuonIndices->at(0)), MuonPz->at(MuonIndices->at(0)), MuonE->at(MuonIndices->at(0)));
            //            cout << " patsyVec Eta " << patsyVec.Eta() << endl;
            //            cout << "MuonEtavector " << MuonEta->at(MuonIndices->at(0)) << endl;
            //            patsyVec.SetPtEtaPhiE(MuonPt->at(MuonIndices->at(0)), MuonEta->at(MuonIndices->at(0)), patsyVec.Phi(), patsyVec.E());
            //            cout << "MuonEtapatsy " << patsyVec.Eta() << endl;
            //            cout << "MuonE " << patsyVec.E() << endl;
            Leptons->push_back(patsyVec);
            pdgId->push_back(-13);
        }
        if (!(M1I < 0)) {            
            patsyVec.SetPxPyPzE(MuonPx->at(MuonIndices->at(1)), MuonPy->at(MuonIndices->at(1)), MuonPz->at(MuonIndices->at(1)), MuonE->at(MuonIndices->at(1)));
            Leptons->push_back(patsyVec);
            pdgId->push_back(13);
        }
        if (!(E0I < 0)) {      
            //            cout << "elec E 0" << ElecE->at(ElecIndices->at(0)) << endl;
            //            cout << "elec Px " << ElecPx->at(ElecIndices->at(0)) << endl;
            //            cout << "elec Py " << ElecPy->at(ElecIndices->at(0)) << endl;
            //            cout << "elec Pz " << ElecPz->at(ElecIndices->at(0)) << endl;
            patsyVec.SetPxPyPzE(ElecPx->at(ElecIndices->at(0)), ElecPy->at(ElecIndices->at(0)), ElecPz->at(ElecIndices->at(0)), ElecE->at(ElecIndices->at(0)));
            //            cout << "elec E 1 " << patsyVec.E() << endl;
            //            cout << "elec E 2 " << ElecE->at(ElecIndices->at(0)) << endl;
            Leptons->push_back(patsyVec);
            pdgId->push_back(-11); 
        }
        if (!(E1I < 0)) {            
            patsyVec.SetPxPyPzE(ElecPx->at(ElecIndices->at(1)), ElecPy->at(ElecIndices->at(1)), ElecPz->at(ElecIndices->at(1)), ElecE->at(ElecIndices->at(1)));
            Leptons->push_back(patsyVec);
            pdgId->push_back(11);
        }
        if (Leptons->size() < 2) continue;
        LeptonPairOvi(Lep0Vec, Lep1Vec, Type, doEvent, Leptons, pdgId, Lep0PDGID, Lep1PDGID);
        if (!doEvent) continue;   
        Lep0Px = Lep0Vec.Px();
        Lep0Py = Lep0Vec.Py();
        Lep0Pz = Lep0Vec.Pz();
        Lep0E  = Lep0Vec.E();
        Lep1Px = Lep1Vec.Px();
        Lep1Py = Lep1Vec.Py();
        Lep1Pz = Lep1Vec.Pz();
        Lep1E  = Lep1Vec.E();
        /*
         cout << "Lep0Px " << Lep0Px << endl;
         cout << "Lep0Py " << Lep0Py << endl;
         cout << "Lep0Pz " << Lep0Pz << endl;
         cout << "Lep1Px " << Lep1Px << endl;
         cout << "Lep1Py " << Lep1Py << endl;
         cout << "Lep1Pz " << Lep1Pz << endl;
         */
        
        Jets = JetInfo(Leptons, JetEt, JetEta, JetPx, JetPy, JetPz, JetE, JetNHF, JetNEF, JetCHF, JetCEF, JetBTag, JetNDaug, JetCharMult, NJets, NBtagJets, BJetIndices, HT);
        BtagJet0Index = (NBtagJets > 0) ? BJetIndices->at(0) : -1;
        BtagJet1Index = (NBtagJets > 1) ? BJetIndices->at(1) : -1;
        if (NJets > 0) {
            Jet0Vec = Jets->at(0);
            Jet0Px = Jet0Vec.Px();
            Jet0Py = Jet0Vec.Py();
            Jet0Pz = Jet0Vec.Pz();
            Jet0E  = Jet0Vec.E();
            if (NJets > 1) {
                Jet1Vec = Jets->at(1); 
                Jet1Px = Jet1Vec.Px();
                Jet1Py = Jet1Vec.Py();
                Jet1Pz = Jet1Vec.Pz();
                Jet1E  = Jet1Vec.E();
            }  
            else {
                Jet1Px = -99999;
                Jet1Py = -99999;
                Jet1Pz = -99999;
                Jet1E  = -99999;
            }
        }
        else {
            Jet0Px = -99999;
            Jet0Py = -99999;
            Jet0Pz = -99999;
            Jet0E  = -99999;
            Jet1Px = -99999;
            Jet1Py = -99999;
            Jet1Pz = -99999;
            Jet1E  = -99999;
            BtagJet0Px = -99999;
            BtagJet0Py = -99999;
            BtagJet0Pz = -99999;
            BtagJet0E  = -99999;
            BtagJet1Px = -99999;
            BtagJet1Py = -99999;
            BtagJet1Pz = -99999;
            BtagJet1E  = -99999;            
        }
        if (NBtagJets > 0) {
            BtagJet0Vec = Jets->at(BJetIndices->at(0));
            BtagJet0Px = BtagJet0Vec.Px();
            BtagJet0Py = BtagJet0Vec.Py();
            BtagJet0Pz = BtagJet0Vec.Pz();
            BtagJet0E  = BtagJet0Vec.E();
            if (NBtagJets > 1) {
                BtagJet1Vec = Jets->at(BJetIndices->at(1));
                BtagJet1Px = BtagJet1Vec.Px();
                BtagJet1Py = BtagJet1Vec.Py();
                BtagJet1Pz = BtagJet1Vec.Pz();
                BtagJet1E  = BtagJet1Vec.E();
            }
            else {
                BtagJet1Px = -99999;
                BtagJet1Py = -99999;
                BtagJet1Pz = -99999;
                BtagJet1E  = -99999; 
            }
        }
        nVtx = VertZ->size();
        if (doPURW && !doData) {
            if (doHackPURW) weight = PileupRW(nVtxSFHist, nVtx);
        }
        
        if(Type==-2) Type=2;
        if (doData) {
            if (Type == 0) {
                if (!(fInName.Contains("DoubleMu"))) continue;
            }
            else if (Type == 1) {
                if (!(fInName.Contains("DoubleEl"))) continue;
            }
            else if (Type == 2) {
                if (!(fInName.Contains("MuEG"))) continue;
            }
        }         
        outTree->Fill();
    }
    cout << "All events done" << endl;
    outputFile->cd();
    cout << "cd-ing to output directory" << endl;
    outputFile->Write();
    h_eventCount->Write();
    cout << "Writing of output file done" << endl;
    outputFile->Close();
    cout << "end of code" << endl;
}
