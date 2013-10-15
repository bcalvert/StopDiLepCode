
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

//#include "../../HeaderFiles/BTagSFUtil.C"
#include "../../HeaderFiles/StopStructDefinitions.h"
#include "../../HeaderFiles/StopFunctionDefinitions_v2.h"
#include "../../HeaderFiles/StopFillFunctions.h"

//#include "PileUpMC.h"

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <map>

using namespace std;
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
    TFile * JetESFile = new TFile("../../PlotShowingCode/Summer13JESUncert.root");
    TH2F * h_JetESUp = (TH2F *) JetESFile->Get("h_JESSystUncertUp");
    TH2F * h_JetESDown = (TH2F *) JetESFile->Get("h_JESSystUncertDown");
    TFile * GenJetSmearFile = new TFile("pfJetResolutionMCtoDataCorrLUT.root");
    TH2F * h_GenJetSmearHist = (TH2F*) GenJetSmearFile->Get("pfJetResolutionMCtoDataCorrLUT");
    TH2F * h_RecoJetLowPtSmearHist = ResolutionHistMaker(TString("JetResolution.txt"));
    vector<TF1> * JetResolutionTF1Vec = VecJetHighPtResolutionTF1();
    
    TString fileInTreeName, fileOutTreeName;
    TString fInName;
    TString fOutName;
    TString outputSavePathString = "outputSavePath";
    //    Double_t relIsoBins[] = {-1.01, 0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16};
    TH1F * h_eventCount = new TH1F("h_eventCount", "; numEvents (# of Entries);", 2, -0.5, 1.5);
    TH1F * h_CutFlow = new TH1F("h_CutFlow", "; Cut Stage; N_{evts} passing Cut Stage", 3, 0.5, 3.5);
    TH1F * h_CutFlow_LepESUp, * h_CutFlow_LepESDown;
    TH1F * h_ElecCharIso = new TH1F("h_ElecCharIso", "; Electron Char. Had. Iso.; N_{evts} / bin", 51, -2, 100);
    TH1F * h_ElecNeutIso = new TH1F("h_ElecNeutIso", "; Electron Neut. Had. Iso.; N_{evts} / bin", 51, -2, 100);
    TH1F * h_ElecPhotIso = new TH1F("h_ElecPhotIso", "; Electron Photon Iso.; N_{evts} / bin", 51, -2, 100);
    
    //    TH1F * h_ElecRelIso = new TH1F("h_ElecRelIso", "; Electron Relative Iso.; N_{evts} / bin", 32, 0., 0.16);
    /////Event Variables/////////////////////
    
    //////bunch of vectors
    //Muons
    MuonEventPointers MEPs;
    /*
    vector<float> * MuonPt, * MuonPx, * MuonPy, * MuonPz, * MuonE, * MuonPFCharHadIso, * MuonPFNeutHadIso, * MuonPFPhotIso, * MuonSumPUPt, * MuonD0, * MuonVertZ;
    vector<float> * PFMuonPt;
    
    vector<bool> * isPFMuon, * isGMPTMuon, * isGlobMuon; //, * isTrackArbitMuon;
    vector<int> * MuonNumMatchStations, * MuonNumLayers, * MuonNumValidPixHitsinTrack;
    */
    //Electrons
    ElectronEventPointers EEPs;
    /*
    vector<float> * ElecPt, * ElecPx, * ElecPy, * ElecPz, * ElecE, * ElecPFCharHadIso, * ElecPFNeutHadIso, * ElecPFPhotIso; 
    vector<float> * PFElecPt;
    vector<float> * ElecSCEta, * ElecDeltaPhiIn, * ElecDeltaEtaIn, * ElecSigIetaIeta, * ElecHtoERatio, * ElecIP, * ElecDZ, * ElecECalE, * ElecSCEOverP;
    vector<int> * ElecNumMissHits;
    vector<bool> * ElecisEB, * ElecisEE;
    //ambiguity in which Muon Pt...rolling with PFMuonPt...never mind ambiguity resolved see https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleLepton2012#Muon_Selection
    //ambiguity in which Electron Pt...rolling with PFElecPt...never mind ambiguity resolved    
    vector<bool> * isPFElectron, * passConvVeto;
    */
    
    
    
    
    //Jets
//    vector<float> * JetPx, * JetPy, * JetPz, * JetE, * JetNHF, * JetNEF, * JetCHF, * JetCEF, * JetBTag;
//    vector<int> * JetNDaug, * JetCharMult, * JetPartFlav;
    PFJetEventPointers PFJEPs;
    
    BTagSFUtil * BTagSFUtilBase, * BTagSFUtilSignal, * BTagSFUtilToUse;
    TString BTagSFAlgorithm = "CSVM", BTagSFDataPeriod = "ABCD", BTagSFSignalString;
    
//    vector<int> * MuonCharge, * ElecCharge;
    
    bool filterECalDeadCell, filterLogErrorTooManyClusters, filterTrackingFailure, filterHCalLaser, filterECalLaserCorr, filterTooManyStripClus, filterManyStripClus, filterEEBadSC;    
    vector<bool> * regularMETFilterVec = new vector<bool>;
    vector<bool> * oppositeMETFilterVec = new vector<bool>;
    bool passTrigDoubleEl, passTrigDoubleMu, passTrigElMu;
    
    //Vertex info (both for nVtx info and also muon dZ)
    vector<float> * VertZ, * VertNDOF, * VertRho; 
    vector<bool> * VertIsFake;
    float firstGoodVertZ;
    
    MEPs.MuonPt = new vector<float>;
    MEPs.MuonPx = new vector<float>; MEPs.MuonPy = new vector<float>; 
    MEPs.MuonPz = new vector<float>; MEPs.MuonEn = new vector<float>;
    MEPs.MuonCharge = new vector<int>;
    MEPs.MuonPFCharHadIso = new vector<float>; MEPs.MuonPFNeutHadIso = new vector<float>; 
    MEPs.MuonPFPhotIso = new vector<float>; MEPs.MuonSumPUPt = new vector<float>;
    MEPs.MuonD0 = new vector<float>; MEPs.MuonVertZ = new vector<float>;
    MEPs.PFMuonPt = new vector<float>;
    
    MEPs.isPFMuon = new vector<bool>; MEPs.isGMPTMuon = new vector<bool>; 
    MEPs.isGlobMuon = new vector<bool>; //isTrackArbitMuon = new vector<bool>;
    
    MEPs.MuonNumMatchStations = new vector<int>; MEPs.MuonNumLayers = new vector<int>;
    MEPs.MuonNumValidPixHitsinTrack = new vector<int>;
    
    EEPs.ElecPt = new vector<float>;
    EEPs.ElecPx = new vector<float>; EEPs.ElecPy = new vector<float>;
    EEPs.ElecPz = new vector<float>; EEPs.ElecEn = new vector<float>;
    EEPs.ElecCharge = new vector<int>;
    EEPs.ElecPFCharHadIso = new vector<float>; EEPs.ElecPFNeutHadIso = new vector<float>;
    EEPs.ElecPFPhotIso = new vector<float>;
    EEPs.PFElecPt = new vector<float>;
    
    EEPs.ElecSCEta = new vector<float>; EEPs.ElecSCEOverP = new vector<float>; EEPs.ElecECalE = new vector<float>;
    EEPs.ElecDeltaPhiIn = new vector<float>; EEPs.ElecDeltaEtaIn = new vector<float>;
    EEPs.ElecSigIetaIeta = new vector<float>; EEPs.ElecHtoERatio = new vector<float>;
    EEPs.ElecIP = new vector<float>; EEPs.ElecDZ = new vector<float>;
    EEPs.ElecNumMissHits = new vector<int>;
    EEPs.isPFElectron = new vector<bool>; EEPs.passConvVeto = new vector<bool>;
    EEPs.ElecisEB = new vector<bool>; EEPs.ElecisEE = new vector<bool>;
    
    PFJEPs.JetPx = new vector<float>; PFJEPs.JetPy = new vector<float>;
    PFJEPs.JetPz = new vector<float>; PFJEPs.JetE = new vector<float>;
    PFJEPs.JetNHF = new vector<float>; PFJEPs.JetNEF = new vector<float>;
    PFJEPs.JetCHF = new vector<float>; PFJEPs.JetCEF = new vector<float>;
    PFJEPs.JetBTag = new vector<float>;
    PFJEPs.JetPartFlav = new vector<int>;        
    PFJEPs.JetNDaug = new vector<int>; PFJEPs.JetCharMult = new vector<int>;
    
    VertZ = new vector<float>; VertNDOF = new vector<float>;
    VertRho = new vector<float>; 
    VertIsFake = new vector<bool>;
    
    int   nVtx, nVtxTrue;
    unsigned int RunNum, EventNum, LumiBlock;
    float MET,MET_Phi,METSig;
    //    float METX, METY;
    float eventRhoIso;
    
    EventLepInfo ELI;
    
    EventJetInfo EJI;
    
    bool  hasTopInfo = 0, hasMETInfo = 0;
    
    //variables for MC systematics
    float MET_LepESUp, MET_Phi_LepESUp;
    float MET_LepESDown, MET_Phi_LepESDown;
    float MET_JetESUp, MET_Phi_JetESUp;
    float MET_JetESDown, MET_Phi_JetESDown;
    //    float METX_LepESUp, METY_LepESUp;
    //    float METX_LepESDown, METY_LepESDown;
    //    float METX_JetESUp, METY_JetESUp;
    //    float METX_JetESDown, METY_JetESDown;
    
    EventLepInfo ELI_LepESUp, ELI_LepESDown;
    
    EventJetInfo EJI_JetESUp, EJI_JetESDown;   
    EventJetInfo EJI_BTagSFUp, EJI_BTagSFDown;
    
    
    //// Smear Jet branches (MET at the end)
    EventJetInfo EJISmear, EJISmear_JetESUp, EJISmear_JetESDown;
    /*
    int NSmearJets, NSmearJets_JetESUp, NSmearJets_JetESDown;
    int NSmearBtagJets, NSmearBtagJets_JetESUp, NSmearBtagJets_JetESDown;
    float SmearJet0Px, SmearJet0Px_JetESUp, SmearJet0Px_JetESDown;
    float SmearJet0Py, SmearJet0Py_JetESUp, SmearJet0Py_JetESDown;
    float SmearJet0Pz, SmearJet0Pz_JetESUp, SmearJet0Pz_JetESDown;
    float SmearJet0E,  SmearJet0E_JetESUp,  SmearJet0E_JetESDown;
    float SmearJet1Px, SmearJet1Px_JetESUp, SmearJet1Px_JetESDown;
    float SmearJet1Py, SmearJet1Py_JetESUp, SmearJet1Py_JetESDown;
    float SmearJet1Pz, SmearJet1Pz_JetESUp, SmearJet1Pz_JetESDown;
    float SmearJet1E,  SmearJet1E_JetESUp,  SmearJet1E_JetESDown;
    
    float SmearBtagJet0Px, SmearBtagJet0Px_JetESUp, SmearBtagJet0Px_JetESDown;
    float SmearBtagJet0Py, SmearBtagJet0Py_JetESUp, SmearBtagJet0Py_JetESDown;
    float SmearBtagJet0Pz, SmearBtagJet0Pz_JetESUp, SmearBtagJet0Pz_JetESDown;
    float SmearBtagJet0E,  SmearBtagJet0E_JetESUp,  SmearBtagJet0E_JetESDown;
    int   SmearBtagJet0Index, SmearBtagJet0Index_JetESUp, SmearBtagJet0Index_JetESDown;
    float SmearBtagJet1Px, SmearBtagJet1Px_JetESUp, SmearBtagJet1Px_JetESDown;
    float SmearBtagJet1Py, SmearBtagJet1Py_JetESUp, SmearBtagJet1Py_JetESDown;
    float SmearBtagJet1Pz, SmearBtagJet1Pz_JetESUp, SmearBtagJet1Pz_JetESDown;
    float SmearBtagJet1E,  SmearBtagJet1E_JetESUp,  SmearBtagJet1E_JetESDown;
    int   SmearBtagJet1Index, SmearBtagJet1Index_JetESUp, SmearBtagJet1Index_JetESDown;
    */
    EventJetInfo EJISmear_BTagSFUp, EJISmear_BTagSFDown;
    /*
    int NSmearBtagJets_BTagSFUp, NSmearBtagJets_BTagSFDown;
    float SmearBtagJet0Px_BTagSFUp, SmearBtagJet0Px_BTagSFDown;
    float SmearBtagJet0Py_BTagSFUp, SmearBtagJet0Py_BTagSFDown;
    float SmearBtagJet0Pz_BTagSFUp, SmearBtagJet0Pz_BTagSFDown;
    float SmearBtagJet0E_BTagSFUp,  SmearBtagJet0E_BTagSFDown;
    int   SmearBtagJet0Index_BTagSFUp, SmearBtagJet0Index_BTagSFDown;
    float SmearBtagJet1Px_BTagSFUp, SmearBtagJet1Px_BTagSFDown;
    float SmearBtagJet1Py_BTagSFUp, SmearBtagJet1Py_BTagSFDown;
    float SmearBtagJet1Pz_BTagSFUp, SmearBtagJet1Pz_BTagSFDown;
    float SmearBtagJet1E_BTagSFUp,  SmearBtagJet1E_BTagSFDown;
    int   SmearBtagJet1Index_BTagSFUp, SmearBtagJet1Index_BTagSFDown;
    */
    EventJetInfo EJISmear_JetSmearUp, EJISmear_JetSmearDown;

    float SmearMET, SmearMET_Phi;
    float SmearMET_LepESUp, SmearMET_LepESDown, SmearMET_Phi_LepESUp, SmearMET_Phi_LepESDown;
    float SmearMET_JetESUp, SmearMET_JetESDown, SmearMET_Phi_JetESUp, SmearMET_Phi_JetESDown;                
    float SmearMET_JetSmearUp, SmearMET_JetSmearDown, SmearMET_Phi_JetSmearUp, SmearMET_Phi_JetSmearDown;

    
    //SUSY particle Gen Mass stuff
    
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
    
    GenJetEventPointers GJEPs;
    GJEPs.genJetPx = new vector<float>;
    GJEPs.genJetPy = new vector<float>;
    GJEPs.genJetPz = new vector<float>;
    GJEPs.genJetEn = new vector<float>;
    GJEPs.genJetInvisE = new vector<float>;
    GJEPs.genJetIsGenJet = new vector<bool>;
    /*
    vector<float> * genJetPx, * genJetPy, * genJetPz, * genJetEn, * genJetInvisE; //* genJetEt, * genJetEta;
    vector<bool> * genJetIsGenJet;
    genJetPx = new vector<float>;
    genJetPy = new vector<float>;
    genJetPz = new vector<float>;
    genJetEn = new vector<float>;
//    genJetEt = new vector<float>;
//    genJetEta = new vector<float>;
    genJetInvisE = new vector<float>;
    genJetIsGenJet = new vector<bool>;
    */
    
    //    float StopMass0, StopMass1, Chi0Mass0, Chi0Mass1;
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
     

     */
    
    float genMET, genMETPhi;
    
    bool doEvent, doEvent_LepESDown, doEvent_LepESUp;    
    
    float weight;
    
    //DESY stuff
    TBranch * b_runNumber, * b_lumiBlock, * b_eventNumber;
    TBranch * b_vertMulti, * b_vertMultiTrue;
    TBranch * b_lepton, * b_lepPdgId, * b_lepPFIso;
    VLV            * leptons = 0;
    vector<int>    * lepPdgId = 0;
    vector<double> * lepPFIso = 0;
    
    TBranch * b_jet, * b_jetBTagCSV;
    VLV            * jets = 0;
    vector<double> * jetBTagCSV = 0;
    
    TBranch * b_GenTop, * b_GenAntiTop;
    LV * genTop = 0, * genAntiTop = 0;
    /// MET branches
    TBranch        * b_met;
    LV             * met = 0;
    
    /// GenMET branch
    TBranch        * b_genMET;
    LV             * genMETVec = 0;
    
    TBranch        * b_GenWeight;
    double         GenWeight;        
    
    ////input cuts/commands    
    //    const double PI = 3.14159265;
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
    bool doSpecRun       = 0;
    bool doSpecRunEvent  = 0;
    int  whichRun        = -1;
    int  whichEvent      = -1;
    int  levelVerbosity  = 0;
    /////loop over inputs
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"-i",2) == 0) { 
            fInName = TString(argv[k+1]);
        }
        else if (strncmp (argv[k],"-w",2) == 0) {
            whichNTupleType = strtol(argv[k+1], NULL, 10);   
        }
        else if (strncmp (argv[k],"doVerbosity", 11) == 0) {
            doVerbosity = 1;
            levelVerbosity = strtol(argv[k+1], NULL, 10);
        }
        else if (strncmp (argv[k],"doPURW", 6) == 0) {
            doPURW = 1;           
        }
        else if (strncmp (argv[k],"doHackPURW", 10) == 0) {
            doHackPURW = 1;   
        }
        else if (strncmp (argv[k],"doPURWOviToDESY", 15) == 0) {
            doPURWOviToDESY = 1;   
        }
        else if (strncmp (argv[k],"gOutDir", 7) == 0) {
            grabOutDir = 1;
        }
        else if (strncmp (argv[k], "OutDirName", 10) == 0) {
            outputSavePathString = TString(argv[k+1]);
        }
        else if (strncmp (argv[k],"isSig", 5) == 0) {
            doSignal = 1;
        }
        else if (strncmp (argv[k],"pEvNum", 6) == 0) {
            printEventNum = 1;
        }    
        else if (strncmp (argv[k],"doSpecRun", 9) == 0) {
            doSpecRun = 1;
            whichRun = strtol(argv[k+1], NULL, 10);
        }        
        else if (strncmp (argv[k],"doSpecEventRun", 14) == 0) {
            doSpecRunEvent = 1;
            doSpecRun = 0;
            whichRun = strtol(argv[k+1], NULL, 10);
            whichEvent = strtol(argv[k+2], NULL, 10);
        }
        else if (strncmp (argv[k],"doMassCut", 9) == 0) {
            doMassCut = 1;
            genStopMassCut = strtol(argv[k+1], NULL, 10);   
            genChi0MassCut = strtol(argv[k+2], NULL, 10);   
            genCharginoMassCut = strtol(argv[k+3], NULL, 10);   
        }
    }
    char Buffer[500];
    ifstream * outDirFile;
    TRegexp fCutSlash("[^/]+$");
    fOutName = "";
    if (grabOutDir) {
        outDirFile = new ifstream(outputSavePathString + TString(".txt"));
        if (!(outDirFile->eof())) {
            outDirFile->getline(Buffer,500);
            fOutName += TString(string(Buffer));
            fOutName += "/"; //in case user forgot a slash
        }
    }
    fOutName += fInName(fCutSlash);
    if (fInName.Contains("MuEG") || fInName.Contains("DoubleMu") || fInName.Contains("DoubleEl") || fInName.Contains("run")) {
        cout << "Running on Data" << endl;
        doData = 1;
    }    
    if (!doData) {
        h_CutFlow_LepESUp = new TH1F("h_CutFlow_LepESUp", "; Cut Stage; N_{evts} passing Cut Stage", 3, 0.5, 3.5);
        h_CutFlow_LepESDown = new TH1F("h_CutFlow_LepESDown", "; Cut Stage; N_{evts} passing Cut Stage", 3, 0.5, 3.5);
    }
    
    if (whichNTupleType == 0) {
        fOutName += "_Oviedo";        
    }
    else {
        fOutName += "_DESY";
    }
    if (doPURW && !doData) fOutName += "_PURW";
    if (doPURWOviToDESY && !doData) fOutName += "OviToDESY";
    if (doSpecRun) {
        fOutName += "_specRun";
    }
    else if (doSpecRunEvent) {
        fOutName += "_specRunEvent";
    }
    fOutName += "_SkimOutput.root";
    cout << "saving to " << fOutName << endl;
    TFile * outputFile;
    outputFile = new TFile(fOutName, "RECREATE");
    fileOutTreeName = (whichNTupleType == 1) ? "DESY" : "Ovi";
    fileOutTreeName += "SkimTree";
    TTree * outTree = new TTree(fileOutTreeName, fileOutTreeName);
    
    fileInTreeName = (whichNTupleType == 1) ? "writeNTuple/NTuple" : "Tree";
    TChain fileTree(fileInTreeName);
    TFile inputFile(fInName + TString(".root"));
    //////////////////////////
    fileTree.Add(fInName + TString(".root"));
    if (whichNTupleType == 0) {
        BTagSFUtilBase = new BTagSFUtil("CSVM");
        
        ////Electron Branches////
        //Electron Kinematic Parameters
        fileTree.SetBranchAddress("T_Elec_PFElecPt",            &EEPs.PFElecPt);
        fileTree.SetBranchAddress("T_Elec_Pt",                  &EEPs.ElecPt);
        fileTree.SetBranchAddress("T_Elec_Px",                  &EEPs.ElecPx);
        fileTree.SetBranchAddress("T_Elec_Py",                  &EEPs.ElecPy);
        fileTree.SetBranchAddress("T_Elec_Pz",                  &EEPs.ElecPz);
        fileTree.SetBranchAddress("T_Elec_Energy",              &EEPs.ElecEn);
        fileTree.SetBranchAddress("T_Elec_Charge",              &EEPs.ElecCharge);        
        
        //Electron Supercluster/Shower parameters
        fileTree.SetBranchAddress("T_Elec_SC_Eta",              &EEPs.ElecSCEta);
        fileTree.SetBranchAddress("T_Elec_deltaPhiIn",          &EEPs.ElecDeltaPhiIn);
        fileTree.SetBranchAddress("T_Elec_deltaEtaIn",          &EEPs.ElecDeltaEtaIn);
        fileTree.SetBranchAddress("T_Elec_HtoE",                &EEPs.ElecHtoERatio);
        fileTree.SetBranchAddress("T_Elec_sigmaIetaIeta",       &EEPs.ElecSigIetaIeta);
        fileTree.SetBranchAddress("T_Elec_eSuperClusterOverP",  &EEPs.ElecSCEOverP);
        fileTree.SetBranchAddress("T_Elec_ecalEnergy",          &EEPs.ElecECalE);        
        
        //Electron Vertex Geometry
        fileTree.SetBranchAddress("T_Elec_IPwrtPV",             &EEPs.ElecIP);
        fileTree.SetBranchAddress("T_Elec_dzwrtPV",             &EEPs.ElecDZ);
        fileTree.SetBranchAddress("T_Elec_nHits",               &EEPs.ElecNumMissHits);
        
        //Electron Isolation parameters
        fileTree.SetBranchAddress("T_Elec_chargedHadronIso",    &EEPs.ElecPFCharHadIso);
        fileTree.SetBranchAddress("T_Elec_neutralHadronIso",    &EEPs.ElecPFNeutHadIso);
        fileTree.SetBranchAddress("T_Elec_photonIso",           &EEPs.ElecPFPhotIso);
        
        //Electron booleans
        fileTree.SetBranchAddress("T_Elec_passConversionVeto",  &EEPs.passConvVeto);
        fileTree.SetBranchAddress("T_Elec_isPF",                &EEPs.isPFElectron);
        fileTree.SetBranchAddress("T_Elec_isEB",                &EEPs.ElecisEB);
        fileTree.SetBranchAddress("T_Elec_isEE",                &EEPs.ElecisEE);
        ////Muon Branches////
        //Muon Kinematics
        fileTree.SetBranchAddress("T_Muon_PFMuonPt",                &MEPs.PFMuonPt);
        fileTree.SetBranchAddress("T_Muon_Pt",                      &MEPs.MuonPt);
        fileTree.SetBranchAddress("T_Muon_Px",                      &MEPs.MuonPx);
        fileTree.SetBranchAddress("T_Muon_Py",                      &MEPs.MuonPy);
        fileTree.SetBranchAddress("T_Muon_Pz",                      &MEPs.MuonPz);
        fileTree.SetBranchAddress("T_Muon_Energy",                  &MEPs.MuonEn);
        fileTree.SetBranchAddress("T_Muon_Charge",                  &MEPs.MuonCharge); 
        
        //Muon Isolation parameters
        fileTree.SetBranchAddress("T_Muon_chargedHadronIsoR03",     &MEPs.MuonPFCharHadIso);
        fileTree.SetBranchAddress("T_Muon_neutralHadronIsoR03",     &MEPs.MuonPFNeutHadIso);
        fileTree.SetBranchAddress("T_Muon_photonIsoR03",            &MEPs.MuonPFPhotIso);
        fileTree.SetBranchAddress("T_Muon_sumPUPtR03",              &MEPs.MuonSumPUPt);
        
        //Muon vertex geometry 
        fileTree.SetBranchAddress("T_Muon_IPwrtAveBSInTrack",       &MEPs.MuonD0);
        fileTree.SetBranchAddress("T_Muon_vz",                      &MEPs.MuonVertZ);
        
        //Muon track/detector parameters
        fileTree.SetBranchAddress("T_Muon_NumOfMatchedStations",    &MEPs.MuonNumMatchStations);
        fileTree.SetBranchAddress("T_Muon_NValidPixelHitsInTrk",    &MEPs.MuonNumValidPixHitsinTrack);
        fileTree.SetBranchAddress("T_Muon_NLayers",                 &MEPs.MuonNumLayers);
        
        //Muon booleans
        fileTree.SetBranchAddress("T_Muon_IsGMPTMuons",             &MEPs.isGMPTMuon);
        fileTree.SetBranchAddress("T_Muon_isPFMuon",                &MEPs.isPFMuon);
        fileTree.SetBranchAddress("T_Muon_IsGlobalMuon",            &MEPs.isGlobMuon);
//        fileTree.SetBranchAddress("T_Muon_IsTrackerMuonArbitrated", &MEPs.isTrackArbitMuon);
        
        ////Jet Branches
        //        fileTree.SetBranchAddress("T_JetAKCHS_Et",        &PFJEPs.JetEt);
        //        fileTree.SetBranchAddress("T_JetAKCHS_Eta",       &PFJEPs.JetEta);  
        fileTree.SetBranchAddress("T_JetAKCHS_Px",                  &PFJEPs.JetPx);
        fileTree.SetBranchAddress("T_JetAKCHS_Py",                  &PFJEPs.JetPy);
        fileTree.SetBranchAddress("T_JetAKCHS_Pz",                  &PFJEPs.JetPz);
        fileTree.SetBranchAddress("T_JetAKCHS_Energy",              &PFJEPs.JetE);
        fileTree.SetBranchAddress("T_JetAKCHS_NeutHadEnergyFrac",   &PFJEPs.JetNHF);
        fileTree.SetBranchAddress("T_JetAKCHS_NeutEmEnergyFrac",    &PFJEPs.JetNEF);
        fileTree.SetBranchAddress("T_JetAKCHS_CharHadEnergyFrac",   &PFJEPs.JetCHF);
        fileTree.SetBranchAddress("T_JetAKCHS_CharEmEnergyFrac",    &PFJEPs.JetCEF);
        fileTree.SetBranchAddress("T_JetAKCHS_Tag_CombSVtx",        &PFJEPs.JetBTag);
        fileTree.SetBranchAddress("T_JetAKCHS_Parton_Flavour",      &PFJEPs.JetPartFlav);
        fileTree.SetBranchAddress("T_JetAKCHS_nDaughters",          &PFJEPs.JetNDaug);
        fileTree.SetBranchAddress("T_JetAKCHS_ChargedMultiplicity", &PFJEPs.JetCharMult);
        
        //Vertex information
        fileTree.SetBranchAddress("T_Vertex_z",      &VertZ);
        fileTree.SetBranchAddress("T_Vertex_ndof",   &VertNDOF);
        fileTree.SetBranchAddress("T_Vertex_rho",    &VertRho);
        fileTree.SetBranchAddress("T_Vertex_isFake", &VertIsFake);
        
        //MET information
        fileTree.SetBranchAddress( "T_METPFTypeI_ET",     &MET );
        fileTree.SetBranchAddress( "T_METPFTypeI_Phi", &MET_Phi );
        fileTree.SetBranchAddress( "T_METPFTypeI_Sig",  &METSig );
        
        //run information
        fileTree.SetBranchAddress( "T_Event_RunNumber", &RunNum );
        fileTree.SetBranchAddress( "T_Event_EventNumber", &EventNum );
        fileTree.SetBranchAddress( "T_Event_LuminosityBlock", &LumiBlock );
        
        //event Rho variable (used for electron Isolation)        
        fileTree.SetBranchAddress( "T_Event_RhoIso", &eventRhoIso );
        
        //Trigger information
        fileTree.SetBranchAddress( "T_passTriggerDoubleMu", &passTrigDoubleMu );
        fileTree.SetBranchAddress( "T_passTriggerDoubleEl", &passTrigDoubleEl );
        fileTree.SetBranchAddress( "T_passTriggerElMu", &passTrigElMu );
        
        //MET Filter information        
        fileTree.SetBranchAddress( "T_EventF_EcalDeadCell", &filterECalDeadCell );
        fileTree.SetBranchAddress( "T_EventF_logErrorTooManyClusters", &filterLogErrorTooManyClusters );
        fileTree.SetBranchAddress( "T_EventF_trackingFailure", &filterTrackingFailure );
        fileTree.SetBranchAddress( "T_EventF_hcalLaser", &filterHCalLaser );
        fileTree.SetBranchAddress( "T_EventF_ecalLaserCorr", &filterECalLaserCorr );
        fileTree.SetBranchAddress( "T_EventF_toomanystripclus", &filterTooManyStripClus );
        fileTree.SetBranchAddress( "T_EventF_manystripclus", &filterManyStripClus );
        fileTree.SetBranchAddress( "T_EventF_eeBadSc", &filterEEBadSC );
        
        //// Generator information
        /// Gen MET info
        fileTree.SetBranchAddress( "T_METgen_ET",                   &genMET );
        fileTree.SetBranchAddress( "T_METgen_Phi",                  &genMETPhi );
        
        /// Gen Jet info
        fileTree.SetBranchAddress("T_JetAKCHS_GenJet_Px",           &GJEPs.genJetPx);
        fileTree.SetBranchAddress("T_JetAKCHS_GenJet_Py",           &GJEPs.genJetPy);
        fileTree.SetBranchAddress("T_JetAKCHS_GenJet_Pz",           &GJEPs.genJetPz);
        fileTree.SetBranchAddress("T_JetAKCHS_GenJet_Energy",       &GJEPs.genJetEn);
        fileTree.SetBranchAddress("T_JetAKCHS_GenJet_InvisibleE",   &GJEPs.genJetInvisE);
        fileTree.SetBranchAddress("T_JetAKCHS_IsGenJet",            &GJEPs.genJetIsGenJet);
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
    }
    else if (whichNTupleType == 1) {
        fileTree.SetBranchAddress("runNumber", &RunNum, &b_runNumber);
        fileTree.SetBranchAddress("lumiBlock", &LumiBlock, &b_lumiBlock);
        fileTree.SetBranchAddress("eventNumber", &EventNum, &b_eventNumber);
        fileTree.SetBranchAddress("vertMulti", &nVtx, &b_vertMulti);
        fileTree.SetBranchAddress("vertMultiTrue", &nVtxTrue, &b_vertMultiTrue);
        fileTree.SetBranchAddress("lepPdgId", &lepPdgId, &b_lepPdgId );	
        fileTree.SetBranchAddress("lepPfIso", &lepPFIso, &b_lepPFIso );	
        fileTree.SetBranchAddress("jets", &jets, &b_jet );
        fileTree.SetBranchAddress("leptons", &leptons, &b_lepton );
        fileTree.SetBranchAddress("jetBTagCSV", &jetBTagCSV, &b_jetBTagCSV );
        fileTree.SetBranchAddress("met", &met, &b_met );
        fileTree.SetBranchAddress("weightGenerator", &GenWeight, &b_GenWeight);
        if (fileTree.GetBranch("GenTop")) {
            hasTopInfo = 1;
            fileTree.SetBranchAddress("GenTop", &genTop, &b_GenTop );   
        }
        if (fileTree.GetBranch("GenAntiTop")) {
            fileTree.SetBranchAddress("GenAntiTop", &genAntiTop, &b_GenAntiTop );   
        }
        if (fileTree.GetBranch("GenMET")) {
            hasMETInfo = 1;
            fileTree.SetBranchAddress("GenMET", &genMETVec, &b_genMET );   
        }
    }
    
    ///Out tree information
    if (whichNTupleType == 0) {
        outTree->Branch("TPassDoubleEl",  &passTrigDoubleEl);
        outTree->Branch("TPassDoubleMu",  &passTrigDoubleMu);
        outTree->Branch("TPassElMu",  &passTrigElMu);
    }
    
    outTree->Branch("TWeight",      &weight);
    outTree->Branch("TNPV",         &nVtx);
    outTree->Branch("TNPV_True",    &nVtxTrue);
    outTree->Branch("TRunNum",      &RunNum);
    outTree->Branch("TEventNum",    &EventNum);
    outTree->Branch("TLumiBlock",   &LumiBlock);
    outTree->Branch("TMET",         &MET);
    outTree->Branch("TMET_Phi",     &MET_Phi);
    outTree->Branch("TMETSig",      &METSig);
    
    outTree->Branch("TDoEvent",                     &ELI.doEvent);
    outTree->Branch("TChannel",                     &ELI.EventDiLepType);
    outTree->Branch("TNIsoElecs_pT20",              &ELI.EventNIsoElecs_pT20);
    outTree->Branch("TNIsoMuons_pT20",              &ELI.EventNIsoMuons_pT20);
    outTree->Branch("TNIsoPosits_pT20",             &ELI.EventNIsoPosits_pT20);
    outTree->Branch("TNIsoMubars_pT20",             &ELI.EventNIsoMubars_pT20);
    outTree->Branch("TNIsoElecs_pT10to20",          &ELI.EventNIsoElecs_pT10to20);
    outTree->Branch("TNIsoMuons_pT10to20",          &ELI.EventNIsoMuons_pT10to20);
    outTree->Branch("TNIsoPosits_pT10to20",         &ELI.EventNIsoPosits_pT10to20);
    outTree->Branch("TNIsoMubars_pT10to20",         &ELI.EventNIsoMubars_pT10to20);
    outTree->Branch("TNViableLepPairsPreMassCut",   &ELI.EventNViableLepPairsPreMassCut);
    
    outTree->Branch("TDiLepMass",       &ELI.EventDiLepMass);
    outTree->Branch("TLep0Px",          &ELI.EventLep0Px); 
    outTree->Branch("TLep0Py",          &ELI.EventLep0Py); 
    outTree->Branch("TLep0Pz",          &ELI.EventLep0Pz); 
    outTree->Branch("TLep0E",           &ELI.EventLep0E);
    outTree->Branch("TLep0PdgId",       &ELI.EventLep0PDGID);
    outTree->Branch("TLep0RelPFIso",    &ELI.EventLep0RelPFIso);
    outTree->Branch("TLep1Px",          &ELI.EventLep1Px); 
    outTree->Branch("TLep1Py",          &ELI.EventLep1Py); 
    outTree->Branch("TLep1Pz",          &ELI.EventLep1Pz); 
    outTree->Branch("TLep1E",           &ELI.EventLep1E);
    outTree->Branch("TLep1PdgId",       &ELI.EventLep1PDGID);
    outTree->Branch("TLep1RelPFIso",    &ELI.EventLep1RelPFIso);    
    
    outTree->Branch("TNJets",           &EJI.EventNJets); 
    outTree->Branch("TNJetsBtag",       &EJI.EventNBtagJets);     
    outTree->Branch("THT",              &EJI.EventHT);
    outTree->Branch("TJet0Px",          &EJI.EventJet0Px); 
    outTree->Branch("TJet0Py",          &EJI.EventJet0Py); 
    outTree->Branch("TJet0Pz",          &EJI.EventJet0Pz); 
    outTree->Branch("TJet0E",           &EJI.EventJet0E);
    outTree->Branch("TJet1Px",          &EJI.EventJet1Px); 
    outTree->Branch("TJet1Py",          &EJI.EventJet1Py); 
    outTree->Branch("TJet1Pz",          &EJI.EventJet1Pz); 
    outTree->Branch("TJet1E",           &EJI.EventJet1E);    
    outTree->Branch("TBtagJet0Px",      &EJI.EventBtagJet0Px); 
    outTree->Branch("TBtagJet0Py",      &EJI.EventBtagJet0Py); 
    outTree->Branch("TBtagJet0Pz",      &EJI.EventBtagJet0Pz); 
    outTree->Branch("TBtagJet0E",       &EJI.EventBtagJet0E);
    outTree->Branch("TBtagJet0Index",   &EJI.EventBtagJet0Index);    
    outTree->Branch("TBtagJet1Px",      &EJI.EventBtagJet1Px); 
    outTree->Branch("TBtagJet1Py",      &EJI.EventBtagJet1Py); 
    outTree->Branch("TBtagJet1Pz",      &EJI.EventBtagJet1Pz); 
    outTree->Branch("TBtagJet1E",       &EJI.EventBtagJet1E);
    outTree->Branch("TBtagJet1Index",   &EJI.EventBtagJet1Index);
    
    if (doSignal) {
        if (fInName.Contains("T2tt")) {
            BTagSFSignalString = "T2tt";
        }
        else if (fInName.Contains("T2bw")) {
            BTagSFSignalString = "T2bw";
        }
        else {
            BTagSFSignalString = "";
        }
        BTagSFUtilSignal = new BTagSFUtil(BTagSFAlgorithm, BTagSFDataPeriod, BTagSFSignalString);
        BTagSFUtilToUse = BTagSFUtilSignal;
        outTree->Branch("TGenStopMass0",     &genStopMass0);
        outTree->Branch("TGenStopMass1",     &genStopMass1);
        outTree->Branch("TGenChi0Mass0",     &genChi0Mass0);
        outTree->Branch("TGenChi0Mass1",     &genChi0Mass1);
        outTree->Branch("TGenCharginoMass0", &genCharginoMass0);
        outTree->Branch("TGenCharginoMass1", &genCharginoMass1);
    }
    else {
        BTagSFUtilToUse = BTagSFUtilBase;
    }
    if (!doData) {
        outTree->Branch("TDoEvent_LepESUp",                     &ELI_LepESUp.doEvent);
        outTree->Branch("TChannel_LepESUp",                     &ELI_LepESUp.EventDiLepType);
        outTree->Branch("TNIsoElecs_pT20_LepESUp",              &ELI_LepESUp.EventNIsoElecs_pT20);
        outTree->Branch("TNIsoMuons_pT20_LepESUp",              &ELI_LepESUp.EventNIsoMuons_pT20);
        outTree->Branch("TNIsoPosits_pT20_LepESUp",             &ELI_LepESUp.EventNIsoPosits_pT20);
        outTree->Branch("TNIsoMubars_pT20_LepESUp",             &ELI_LepESUp.EventNIsoMubars_pT20);
        outTree->Branch("TNIsoElecs_pT10to20_LepESUp",          &ELI_LepESUp.EventNIsoElecs_pT10to20);
        outTree->Branch("TNIsoMuons_pT10to20_LepESUp",          &ELI_LepESUp.EventNIsoMuons_pT10to20);
        outTree->Branch("TNIsoPosits_pT10to20_LepESUp",         &ELI_LepESUp.EventNIsoPosits_pT10to20);
        outTree->Branch("TNIsoMubars_pT10to20_LepESUp",         &ELI_LepESUp.EventNIsoMubars_pT10to20);
        outTree->Branch("TNViableLepPairsPreMassCut_LepESUp",   &ELI_LepESUp.EventNViableLepPairsPreMassCut);
        
        outTree->Branch("TDiLepMass_LepESUp",                   &ELI_LepESUp.EventDiLepMass);
        outTree->Branch("TLep0Px_LepESUp",                      &ELI_LepESUp.EventLep0Px); 
        outTree->Branch("TLep0Py_LepESUp",                      &ELI_LepESUp.EventLep0Py); 
        outTree->Branch("TLep0Pz_LepESUp",                      &ELI_LepESUp.EventLep0Pz); 
        outTree->Branch("TLep0E_LepESUp",                       &ELI_LepESUp.EventLep0E);
        outTree->Branch("TLep0PdgId_LepESUp",                   &ELI_LepESUp.EventLep0PDGID);
        outTree->Branch("TLep0RelPFIso_LepESUp",                &ELI_LepESUp.EventLep0RelPFIso);
        outTree->Branch("TLep1Px_LepESUp",                      &ELI_LepESUp.EventLep1Px); 
        outTree->Branch("TLep1Py_LepESUp",                      &ELI_LepESUp.EventLep1Py); 
        outTree->Branch("TLep1Pz_LepESUp",                      &ELI_LepESUp.EventLep1Pz); 
        outTree->Branch("TLep1E_LepESUp",                       &ELI_LepESUp.EventLep1E);
        outTree->Branch("TLep1PdgId_LepESUp",                   &ELI_LepESUp.EventLep1PDGID);
        outTree->Branch("TLep1RelPFIso_LepESUp",                &ELI_LepESUp.EventLep1RelPFIso);        
        
        outTree->Branch("TDoEvent_LepESDown",                   &ELI_LepESDown.doEvent);
        outTree->Branch("TChannel_LepESDown",                   &ELI_LepESDown.EventDiLepType);
        outTree->Branch("TNIsoElecs_pT20_LepESDown",            &ELI_LepESDown.EventNIsoElecs_pT20);
        outTree->Branch("TNIsoMuons_pT20_LepESDown",            &ELI_LepESDown.EventNIsoMuons_pT20);
        outTree->Branch("TNIsoPosits_pT20_LepESDown",           &ELI_LepESDown.EventNIsoPosits_pT20);
        outTree->Branch("TNIsoMubars_pT20_LepESDown",           &ELI_LepESDown.EventNIsoMubars_pT20);
        outTree->Branch("TNIsoElecs_pT10to20_LepESDown",        &ELI_LepESDown.EventNIsoElecs_pT10to20);
        outTree->Branch("TNIsoMuons_pT10to20_LepESDown",        &ELI_LepESDown.EventNIsoMuons_pT10to20);
        outTree->Branch("TNIsoPosits_pT10to20_LepESDown",       &ELI_LepESDown.EventNIsoPosits_pT10to20);
        outTree->Branch("TNIsoMubars_pT10to20_LepESDown",       &ELI_LepESDown.EventNIsoMubars_pT10to20);
        outTree->Branch("TNViableLepPairsPreMassCut_LepESDown", &ELI_LepESDown.EventNViableLepPairsPreMassCut);
        
        outTree->Branch("TDiLepMass_LepESDown",                 &ELI_LepESDown.EventDiLepMass);
        outTree->Branch("TLep0Px_LepESDown",                    &ELI_LepESDown.EventLep0Px); 
        outTree->Branch("TLep0Py_LepESDown",                    &ELI_LepESDown.EventLep0Py); 
        outTree->Branch("TLep0Pz_LepESDown",                    &ELI_LepESDown.EventLep0Pz); 
        outTree->Branch("TLep0E_LepESDown",                     &ELI_LepESDown.EventLep0E);
        outTree->Branch("TLep0PdgId_LepESDown",                 &ELI_LepESDown.EventLep0PDGID);
        outTree->Branch("TLep0RelPFIso_LepESDown",              &ELI_LepESDown.EventLep0RelPFIso);
        outTree->Branch("TLep1Px_LepESDown",                    &ELI_LepESDown.EventLep1Px); 
        outTree->Branch("TLep1Py_LepESDown",                    &ELI_LepESDown.EventLep1Py); 
        outTree->Branch("TLep1Pz_LepESDown",                    &ELI_LepESDown.EventLep1Pz); 
        outTree->Branch("TLep1E_LepESDown",                     &ELI_LepESDown.EventLep1E);
        outTree->Branch("TLep1PdgId_LepESDown",                 &ELI_LepESDown.EventLep1PDGID);
        outTree->Branch("TLep1RelPFIso_LepESDown",              &ELI_LepESDown.EventLep1RelPFIso);
              
        outTree->Branch("TNJets_JetESUp",           &EJI_JetESUp.EventNJets); 
        outTree->Branch("TNJetsBtag_JetESUp",       &EJI_JetESUp.EventNBtagJets);     
        outTree->Branch("THT_JetESUp",              &EJI_JetESUp.EventHT);
        outTree->Branch("TJet0Px_JetESUp",          &EJI_JetESUp.EventJet0Px); 
        outTree->Branch("TJet0Py_JetESUp",          &EJI_JetESUp.EventJet0Py); 
        outTree->Branch("TJet0Pz_JetESUp",          &EJI_JetESUp.EventJet0Pz); 
        outTree->Branch("TJet0E_JetESUp",           &EJI_JetESUp.EventJet0E);
        outTree->Branch("TJet1Px_JetESUp",          &EJI_JetESUp.EventJet1Px); 
        outTree->Branch("TJet1Py_JetESUp",          &EJI_JetESUp.EventJet1Py); 
        outTree->Branch("TJet1Pz_JetESUp",          &EJI_JetESUp.EventJet1Pz); 
        outTree->Branch("TJet1E_JetESUp",           &EJI_JetESUp.EventJet1E);    
        outTree->Branch("TBtagJet0Px_JetESUp",      &EJI_JetESUp.EventBtagJet0Px); 
        outTree->Branch("TBtagJet0Py_JetESUp",      &EJI_JetESUp.EventBtagJet0Py); 
        outTree->Branch("TBtagJet0Pz_JetESUp",      &EJI_JetESUp.EventBtagJet0Pz); 
        outTree->Branch("TBtagJet0E_JetESUp",       &EJI_JetESUp.EventBtagJet0E);
        outTree->Branch("TBtagJet0Index_JetESUp",   &EJI_JetESUp.EventBtagJet0Index);    
        outTree->Branch("TBtagJet1Px_JetESUp",      &EJI_JetESUp.EventBtagJet1Px); 
        outTree->Branch("TBtagJet1Py_JetESUp",      &EJI_JetESUp.EventBtagJet1Py); 
        outTree->Branch("TBtagJet1Pz_JetESUp",      &EJI_JetESUp.EventBtagJet1Pz); 
        outTree->Branch("TBtagJet1E_JetESUp",       &EJI_JetESUp.EventBtagJet1E);
        outTree->Branch("TBtagJet1Index_JetESUp",   &EJI_JetESUp.EventBtagJet1Index);        
        
        outTree->Branch("TNJets_JetESDown",           &EJI_JetESDown.EventNJets); 
        outTree->Branch("TNJetsBtag_JetESDown",       &EJI_JetESDown.EventNBtagJets);     
        outTree->Branch("THT_JetESDown",               &EJI_JetESDown.EventHT);
        outTree->Branch("TJet0Px_JetESDown",          &EJI_JetESDown.EventJet0Px); 
        outTree->Branch("TJet0Py_JetESDown",          &EJI_JetESDown.EventJet0Py); 
        outTree->Branch("TJet0Pz_JetESDown",          &EJI_JetESDown.EventJet0Pz); 
        outTree->Branch("TJet0E_JetESDown",           &EJI_JetESDown.EventJet0E);
        outTree->Branch("TJet1Px_JetESDown",          &EJI_JetESDown.EventJet1Px); 
        outTree->Branch("TJet1Py_JetESDown",          &EJI_JetESDown.EventJet1Py); 
        outTree->Branch("TJet1Pz_JetESDown",          &EJI_JetESDown.EventJet1Pz); 
        outTree->Branch("TJet1E_JetESDown",           &EJI_JetESDown.EventJet1E);    
        outTree->Branch("TBtagJet0Px_JetESDown",      &EJI_JetESDown.EventBtagJet0Px); 
        outTree->Branch("TBtagJet0Py_JetESDown",      &EJI_JetESDown.EventBtagJet0Py); 
        outTree->Branch("TBtagJet0Pz_JetESDown",      &EJI_JetESDown.EventBtagJet0Pz); 
        outTree->Branch("TBtagJet0E_JetESDown",       &EJI_JetESDown.EventBtagJet0E);
        outTree->Branch("TBtagJet0Index_JetESDown",   &EJI_JetESDown.EventBtagJet0Index);    
        outTree->Branch("TBtagJet1Px_JetESDown",      &EJI_JetESDown.EventBtagJet1Px); 
        outTree->Branch("TBtagJet1Py_JetESDown",      &EJI_JetESDown.EventBtagJet1Py); 
        outTree->Branch("TBtagJet1Pz_JetESDown",      &EJI_JetESDown.EventBtagJet1Pz); 
        outTree->Branch("TBtagJet1E_JetESDown",       &EJI_JetESDown.EventBtagJet1E);
        outTree->Branch("TBtagJet1Index_JetESDown",   &EJI_JetESDown.EventBtagJet1Index);
        
        
                
        outTree->Branch("TNJetsBtag_BTagSFUp",       &EJI_BTagSFUp.EventNBtagJets);   
        outTree->Branch("TBtagJet0Px_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet0Px); 
        outTree->Branch("TBtagJet0Py_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet0Py); 
        outTree->Branch("TBtagJet0Pz_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet0Pz); 
        outTree->Branch("TBtagJet0E_BTagSFUp",       &EJI_BTagSFUp.EventBtagJet0E);
        outTree->Branch("TBtagJet0Index_BTagSFUp",   &EJI_BTagSFUp.EventBtagJet0Index);    
        outTree->Branch("TBtagJet1Px_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet1Px); 
        outTree->Branch("TBtagJet1Py_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet1Py); 
        outTree->Branch("TBtagJet1Pz_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet1Pz); 
        outTree->Branch("TBtagJet1E_BTagSFUp",       &EJI_BTagSFUp.EventBtagJet1E);
        outTree->Branch("TBtagJet1Index_BTagSFUp",   &EJI_BTagSFUp.EventBtagJet1Index);
                
        outTree->Branch("TNJetsBtag_BTagSFDown",       &EJI_BTagSFDown.EventNBtagJets);
        outTree->Branch("TBtagJet0Px_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet0Px); 
        outTree->Branch("TBtagJet0Py_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet0Py); 
        outTree->Branch("TBtagJet0Pz_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet0Pz); 
        outTree->Branch("TBtagJet0E_BTagSFDown",       &EJI_BTagSFDown.EventBtagJet0E);
        outTree->Branch("TBtagJet0Index_BTagSFDown",   &EJI_BTagSFDown.EventBtagJet0Index);    
        outTree->Branch("TBtagJet1Px_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet1Px); 
        outTree->Branch("TBtagJet1Py_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet1Py); 
        outTree->Branch("TBtagJet1Pz_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet1Pz); 
        outTree->Branch("TBtagJet1E_BTagSFDown",       &EJI_BTagSFDown.EventBtagJet1E);
        outTree->Branch("TBtagJet1Index_BTagSFDown",   &EJI_BTagSFDown.EventBtagJet1Index);        
        
        outTree->Branch("TMET_LepESUp",             &MET_LepESUp);
        outTree->Branch("TMET_LepESDown",           &MET_LepESDown);
        outTree->Branch("TMET_JetESUp",             &MET_JetESUp);
        outTree->Branch("TMET_JetESDown",           &MET_JetESDown);
        outTree->Branch("TMET_Phi_LepESUp",         &MET_Phi_LepESUp);
        outTree->Branch("TMET_Phi_LepESDown",       &MET_Phi_LepESDown);
        outTree->Branch("TMET_Phi_JetESUp",         &MET_Phi_JetESUp);
        outTree->Branch("TMET_Phi_JetESDown",       &MET_Phi_JetESDown); 
        
        ////####SMEAR JET BRANCHES#####///////
        outTree->Branch("TSmearMET",                &SmearMET);
        outTree->Branch("TSmearMET_Phi",            &SmearMET_Phi);
        
        outTree->Branch("TNSmearJets",              &EJISmear.EventNJets); 
        outTree->Branch("TNSmearJetsBtag",          &EJISmear.EventNBtagJets);
        outTree->Branch("THTSmear",                 &EJISmear.EventHT);
        
        outTree->Branch("TSmearJet0Px",             &EJISmear.EventJet0Px); 
        outTree->Branch("TSmearJet0Py",             &EJISmear.EventJet0Py); 
        outTree->Branch("TSmearJet0Pz",             &EJISmear.EventJet0Pz); 
        outTree->Branch("TSmearJet0E",              &EJISmear.EventJet0E);
        outTree->Branch("TSmearJet1Px",             &EJISmear.EventJet1Px); 
        outTree->Branch("TSmearJet1Py",             &EJISmear.EventJet1Py); 
        outTree->Branch("TSmearJet1Pz",             &EJISmear.EventJet1Pz); 
        outTree->Branch("TSmearJet1E",              &EJISmear.EventJet1E);
        
        outTree->Branch("TSmearBtagJet0Px",         &EJISmear.EventBtagJet0Px); 
        outTree->Branch("TSmearBtagJet0Py",         &EJISmear.EventBtagJet0Py); 
        outTree->Branch("TSmearBtagJet0Pz",         &EJISmear.EventBtagJet0Pz); 
        outTree->Branch("TSmearBtagJet0E",          &EJISmear.EventBtagJet0E);
        outTree->Branch("TSmearBtagJet0Index",      &EJISmear.EventBtagJet0Index);    
        outTree->Branch("TSmearBtagJet1Px",         &EJISmear.EventBtagJet1Px); 
        outTree->Branch("TSmearBtagJet1Py",         &EJISmear.EventBtagJet1Py); 
        outTree->Branch("TSmearBtagJet1Pz",         &EJISmear.EventBtagJet1Pz); 
        outTree->Branch("TSmearBtagJet1E",          &EJISmear.EventBtagJet1E);
        outTree->Branch("TSmearBtagJet1Index",      &EJISmear.EventBtagJet1Index);
        
        ////####SMEAR JET Systematics Branches#####///////
        outTree->Branch("TSmearMET_LepESUp",        &SmearMET_LepESUp);
        outTree->Branch("TSmearMET_LepESDown",      &SmearMET_LepESDown);
        outTree->Branch("TSmearMET_Phi_LepESUp",    &SmearMET_Phi_LepESUp);
        outTree->Branch("TSmearMET_Phi_LepESDown",  &SmearMET_Phi_LepESDown);
        outTree->Branch("TSmearMET_JetESUp",        &SmearMET_JetESUp);
        outTree->Branch("TSmearMET_JetESDown",      &SmearMET_JetESDown);
        outTree->Branch("TSmearMET_Phi_JetESUp",    &SmearMET_Phi_JetESUp);
        outTree->Branch("TSmearMET_Phi_JetESDown",  &SmearMET_Phi_JetESDown);        
        outTree->Branch("TSmearMET_JetSmearUp",        &SmearMET_JetSmearUp);
        outTree->Branch("TSmearMET_JetSmearDown",      &SmearMET_JetSmearDown);
        outTree->Branch("TSmearMET_Phi_JetSmearUp",    &SmearMET_Phi_JetSmearUp);
        outTree->Branch("TSmearMET_Phi_JetSmearDown",  &SmearMET_Phi_JetSmearDown);
        
        outTree->Branch("TNSmearJets_JetESUp",              &EJISmear_JetESUp.EventNJets); 
        outTree->Branch("TNSmearJetsBtag_JetESUp",          &EJISmear_JetESUp.EventNBtagJets);
        outTree->Branch("THTSmear_JetESUp",                 &EJISmear_JetESUp.EventHT);
        
        outTree->Branch("TNSmearJets_JetESDown",            &EJISmear_JetESDown.EventNJets); 
        outTree->Branch("TNSmearJetsBtag_JetESDown",        &EJISmear_JetESDown.EventNBtagJets);
        outTree->Branch("THTSmear_JetESDown",               &EJISmear_JetESDown.EventHT);
        
        outTree->Branch("TSmearJet0Px_JetESUp",             &EJISmear_JetESUp.EventJet0Px); 
        outTree->Branch("TSmearJet0Py_JetESUp",             &EJISmear_JetESUp.EventJet0Py); 
        outTree->Branch("TSmearJet0Pz_JetESUp",             &EJISmear_JetESUp.EventJet0Pz); 
        outTree->Branch("TSmearJet0E_JetESUp",              &EJISmear_JetESUp.EventJet0E);
        outTree->Branch("TSmearJet1Px_JetESUp",             &EJISmear_JetESUp.EventJet1Px); 
        outTree->Branch("TSmearJet1Py_JetESUp",             &EJISmear_JetESUp.EventJet1Py); 
        outTree->Branch("TSmearJet1Pz_JetESUp",             &EJISmear_JetESUp.EventJet1Pz); 
        outTree->Branch("TSmearJet1E_JetESUp",              &EJISmear_JetESUp.EventJet1E);       
        
        outTree->Branch("TSmearJet0Px_JetESDown",           &EJISmear_JetESDown.EventJet0Px); 
        outTree->Branch("TSmearJet0Py_JetESDown",           &EJISmear_JetESDown.EventJet0Py); 
        outTree->Branch("TSmearJet0Pz_JetESDown",           &EJISmear_JetESDown.EventJet0Pz); 
        outTree->Branch("TSmearJet0E_JetESDown",            &EJISmear_JetESDown.EventJet0E);
        outTree->Branch("TSmearJet1Px_JetESDown",           &EJISmear_JetESDown.EventJet1Px); 
        outTree->Branch("TSmearJet1Py_JetESDown",           &EJISmear_JetESDown.EventJet1Py); 
        outTree->Branch("TSmearJet1Pz_JetESDown",           &EJISmear_JetESDown.EventJet1Pz); 
        outTree->Branch("TSmearJet1E_JetESDown",            &EJISmear_JetESDown.EventJet1E);
        
        outTree->Branch("TSmearBtagJet0Px_JetESUp",         &EJISmear_JetESUp.EventBtagJet0Px); 
        outTree->Branch("TSmearBtagJet0Py_JetESUp",         &EJISmear_JetESUp.EventBtagJet0Py); 
        outTree->Branch("TSmearBtagJet0Pz_JetESUp",         &EJISmear_JetESUp.EventBtagJet0Pz); 
        outTree->Branch("TSmearBtagJet0E_JetESUp",          &EJISmear_JetESUp.EventBtagJet0E);
        outTree->Branch("TSmearBtagJet0Index_JetESUp",      &EJISmear_JetESUp.EventBtagJet0Index);    
        outTree->Branch("TSmearBtagJet1Px_JetESUp",         &EJISmear_JetESUp.EventBtagJet1Px); 
        outTree->Branch("TSmearBtagJet1Py_JetESUp",         &EJISmear_JetESUp.EventBtagJet1Py); 
        outTree->Branch("TSmearBtagJet1Pz_JetESUp",         &EJISmear_JetESUp.EventBtagJet1Pz); 
        outTree->Branch("TSmearBtagJet1E_JetESUp",          &EJISmear_JetESUp.EventBtagJet1E);
        outTree->Branch("TSmearBtagJet1Index_JetESUp",      &EJISmear_JetESUp.EventBtagJet1Index); 
        
        outTree->Branch("TSmearBtagJet0Px_JetESDown",       &EJISmear_JetESDown.EventBtagJet0Px); 
        outTree->Branch("TSmearBtagJet0Py_JetESDown",       &EJISmear_JetESDown.EventBtagJet0Py); 
        outTree->Branch("TSmearBtagJet0Pz_JetESDown",       &EJISmear_JetESDown.EventBtagJet0Pz); 
        outTree->Branch("TSmearBtagJet0E_JetESDown",        &EJISmear_JetESDown.EventBtagJet0E);
        outTree->Branch("TSmearBtagJet0Index_JetESDown",    &EJISmear_JetESDown.EventBtagJet0Index);    
        outTree->Branch("TSmearBtagJet1Px_JetESDown",       &EJISmear_JetESDown.EventBtagJet1Px); 
        outTree->Branch("TSmearBtagJet1Py_JetESDown",       &EJISmear_JetESDown.EventBtagJet1Py); 
        outTree->Branch("TSmearBtagJet1Pz_JetESDown",       &EJISmear_JetESDown.EventBtagJet1Pz); 
        outTree->Branch("TSmearBtagJet1E_JetESDown",        &EJISmear_JetESDown.EventBtagJet1E);
        outTree->Branch("TSmearBtagJet1Index_JetESDown",    &EJISmear_JetESDown.EventBtagJet1Index);
        
        outTree->Branch("TNSmearJetsBtag_BTagSFUp",         &EJISmear_BTagSFUp.EventNBtagJets);        
        outTree->Branch("TSmearBtagJet0Px_BTagSFUp",        &EJISmear_BTagSFUp.EventBtagJet0Px); 
        outTree->Branch("TSmearBtagJet0Py_BTagSFUp",        &EJISmear_BTagSFUp.EventBtagJet0Py); 
        outTree->Branch("TSmearBtagJet0Pz_BTagSFUp",        &EJISmear_BTagSFUp.EventBtagJet0Pz); 
        outTree->Branch("TSmearBtagJet0E_BTagSFUp",         &EJISmear_BTagSFUp.EventBtagJet0E);
        outTree->Branch("TSmearBtagJet0Index_BTagSFUp",     &EJISmear_BTagSFUp.EventBtagJet0Index);    
        outTree->Branch("TSmearBtagJet1Px_BTagSFUp",        &EJISmear_BTagSFUp.EventBtagJet1Px); 
        outTree->Branch("TSmearBtagJet1Py_BTagSFUp",        &EJISmear_BTagSFUp.EventBtagJet1Py); 
        outTree->Branch("TSmearBtagJet1Pz_BTagSFUp",        &EJISmear_BTagSFUp.EventBtagJet1Pz); 
        outTree->Branch("TSmearBtagJet1E_BTagSFUp",         &EJISmear_BTagSFUp.EventBtagJet1E);
        outTree->Branch("TSmearBtagJet1Index_BTagSFUp",     &EJISmear_BTagSFUp.EventBtagJet1Index);
        
        outTree->Branch("TNSmearJetsBtag_BTagSFDown",       &EJISmear_BTagSFDown.EventNBtagJets);        
        outTree->Branch("TSmearBtagJet0Px_BTagSFDown",      &EJISmear_BTagSFDown.EventBtagJet0Px); 
        outTree->Branch("TSmearBtagJet0Py_BTagSFDown",      &EJISmear_BTagSFDown.EventBtagJet0Py); 
        outTree->Branch("TSmearBtagJet0Pz_BTagSFDown",      &EJISmear_BTagSFDown.EventBtagJet0Pz); 
        outTree->Branch("TSmearBtagJet0E_BTagSFDown",       &EJISmear_BTagSFDown.EventBtagJet0E);
        outTree->Branch("TSmearBtagJet0Index_BTagSFDown",   &EJISmear_BTagSFDown.EventBtagJet0Index);    
        outTree->Branch("TSmearBtagJet1Px_BTagSFDown",      &EJISmear_BTagSFDown.EventBtagJet1Px); 
        outTree->Branch("TSmearBtagJet1Py_BTagSFDown",      &EJISmear_BTagSFDown.EventBtagJet1Py); 
        outTree->Branch("TSmearBtagJet1Pz_BTagSFDown",      &EJISmear_BTagSFDown.EventBtagJet1Pz); 
        outTree->Branch("TSmearBtagJet1E_BTagSFDown",       &EJISmear_BTagSFDown.EventBtagJet1E);
        outTree->Branch("TSmearBtagJet1Index_BTagSFDown",   &EJISmear_BTagSFDown.EventBtagJet1Index);
        
        
        
        
        outTree->Branch("TNSmearJets_JetSmearUp",           &EJISmear_JetSmearUp.EventNJets); 
        outTree->Branch("TNSmearJetsBtag_JetSmearUp",       &EJISmear_JetSmearUp.EventNBtagJets);
        outTree->Branch("THTSmear_JetSmearUp",              &EJISmear_JetSmearUp.EventHT);
        
        outTree->Branch("TNSmearJets_JetSmearDown",         &EJISmear_JetSmearDown.EventNJets); 
        outTree->Branch("TNSmearJetsBtag_JetSmearDown",     &EJISmear_JetSmearDown.EventNBtagJets);
        outTree->Branch("THTSmear_JetSmearDown",            &EJISmear_JetSmearDown.EventHT);
        
        outTree->Branch("TSmearJet0Px_JetSmearUp",          &EJISmear_JetSmearUp.EventJet0Px); 
        outTree->Branch("TSmearJet0Py_JetSmearUp",          &EJISmear_JetSmearUp.EventJet0Py); 
        outTree->Branch("TSmearJet0Pz_JetSmearUp",          &EJISmear_JetSmearUp.EventJet0Pz); 
        outTree->Branch("TSmearJet0E_JetSmearUp",           &EJISmear_JetSmearUp.EventJet0E);
        outTree->Branch("TSmearJet1Px_JetSmearUp",          &EJISmear_JetSmearUp.EventJet1Px); 
        outTree->Branch("TSmearJet1Py_JetSmearUp",          &EJISmear_JetSmearUp.EventJet1Py); 
        outTree->Branch("TSmearJet1Pz_JetSmearUp",          &EJISmear_JetSmearUp.EventJet1Pz); 
        outTree->Branch("TSmearJet1E_JetSmearUp",           &EJISmear_JetSmearUp.EventJet1E);       
        
        outTree->Branch("TSmearJet0Px_JetSmearDown",        &EJISmear_JetSmearDown.EventJet0Px); 
        outTree->Branch("TSmearJet0Py_JetSmearDown",        &EJISmear_JetSmearDown.EventJet0Py); 
        outTree->Branch("TSmearJet0Pz_JetSmearDown",        &EJISmear_JetSmearDown.EventJet0Pz); 
        outTree->Branch("TSmearJet0E_JetSmearDown",         &EJISmear_JetSmearDown.EventJet0E);
        outTree->Branch("TSmearJet1Px_JetSmearDown",        &EJISmear_JetSmearDown.EventJet1Px); 
        outTree->Branch("TSmearJet1Py_JetSmearDown",        &EJISmear_JetSmearDown.EventJet1Py); 
        outTree->Branch("TSmearJet1Pz_JetSmearDown",        &EJISmear_JetSmearDown.EventJet1Pz); 
        outTree->Branch("TSmearJet1E_JetSmearDown",         &EJISmear_JetSmearDown.EventJet1E);
        
        outTree->Branch("TSmearBtagJet0Px_JetSmearUp",      &EJISmear_JetSmearUp.EventBtagJet0Px); 
        outTree->Branch("TSmearBtagJet0Py_JetSmearUp",      &EJISmear_JetSmearUp.EventBtagJet0Py); 
        outTree->Branch("TSmearBtagJet0Pz_JetSmearUp",      &EJISmear_JetSmearUp.EventBtagJet0Pz); 
        outTree->Branch("TSmearBtagJet0E_JetSmearUp",       &EJISmear_JetSmearUp.EventBtagJet0E);
        outTree->Branch("TSmearBtagJet0Index_JetSmearUp",   &EJISmear_JetSmearUp.EventBtagJet0Index);    
        outTree->Branch("TSmearBtagJet1Px_JetSmearUp",      &EJISmear_JetSmearUp.EventBtagJet1Px); 
        outTree->Branch("TSmearBtagJet1Py_JetSmearUp",      &EJISmear_JetSmearUp.EventBtagJet1Py); 
        outTree->Branch("TSmearBtagJet1Pz_JetSmearUp",      &EJISmear_JetSmearUp.EventBtagJet1Pz); 
        outTree->Branch("TSmearBtagJet1E_JetSmearUp",       &EJISmear_JetSmearUp.EventBtagJet1E);
        outTree->Branch("TSmearBtagJet1Index_JetSmearUp",   &EJISmear_JetSmearUp.EventBtagJet1Index); 
        
        outTree->Branch("TSmearBtagJet0Px_JetSmearDown",    &EJISmear_JetSmearDown.EventBtagJet0Px); 
        outTree->Branch("TSmearBtagJet0Py_JetSmearDown",    &EJISmear_JetSmearDown.EventBtagJet0Py); 
        outTree->Branch("TSmearBtagJet0Pz_JetSmearDown",    &EJISmear_JetSmearDown.EventBtagJet0Pz); 
        outTree->Branch("TSmearBtagJet0E_JetSmearDown",     &EJISmear_JetSmearDown.EventBtagJet0E);
        outTree->Branch("TSmearBtagJet0Index_JetSmearDown", &EJISmear_JetSmearDown.EventBtagJet0Index);    
        outTree->Branch("TSmearBtagJet1Px_JetSmearDown",    &EJISmear_JetSmearDown.EventBtagJet1Px); 
        outTree->Branch("TSmearBtagJet1Py_JetSmearDown",    &EJISmear_JetSmearDown.EventBtagJet1Py); 
        outTree->Branch("TSmearBtagJet1Pz_JetSmearDown",    &EJISmear_JetSmearDown.EventBtagJet1Pz); 
        outTree->Branch("TSmearBtagJet1E_JetSmearDown",     &EJISmear_JetSmearDown.EventBtagJet1E);
        outTree->Branch("TSmearBtagJet1Index_JetSmearDown", &EJISmear_JetSmearDown.EventBtagJet1Index);
        
        
        ///###GENMET Branches###/////
        outTree->Branch("TGenMET",                      &genMET);
        outTree->Branch("TGenMETPhi",                   &genMETPhi);                
        
        /// Status 3 lead particle branches
        outTree->Branch("TGenTopSt3_0_PDGID",           &genTopSt3_0_PID );
        outTree->Branch("TGenTopSt3_0_FirstMom",        &genTopSt3_0_FirstMom );
        outTree->Branch("TGenTopSt3_0_Index",           &genTopSt3_0_Index );
        outTree->Branch("TGenTopSt3_0_Energy",          &genTopSt3_0_Energy );
        outTree->Branch("TGenTopSt3_0_Pt",              &genTopSt3_0_Pt );
        outTree->Branch("TGenTopSt3_0_Eta",             &genTopSt3_0_Eta );
        outTree->Branch("TGenTopSt3_0_Phi",             &genTopSt3_0_Phi );
        
        outTree->Branch("TGenBSt3_0_PDGID",             &genBSt3_0_PID );
        outTree->Branch("TGenBSt3_0_FirstMom",          &genBSt3_0_FirstMom );
        outTree->Branch("TGenBSt3_0_Index",             &genBSt3_0_Index );
        outTree->Branch("TGenBSt3_0_Energy",            &genBSt3_0_Energy );
        outTree->Branch("TGenBSt3_0_Pt",                &genBSt3_0_Pt );
        outTree->Branch("TGenBSt3_0_Eta",               &genBSt3_0_Eta );
        outTree->Branch("TGenBSt3_0_Phi",               &genBSt3_0_Phi );
        
        outTree->Branch("TGenMuonSt3_0_PDGID",          &genMuonSt3_0_PID );
        outTree->Branch("TGenMuonSt3_0_FirstMom",       &genMuonSt3_0_FirstMom );
        outTree->Branch("TGenMuonSt3_0_Index",          &genMuonSt3_0_Index );
        outTree->Branch("TGenMuonSt3_0_Energy",         &genMuonSt3_0_Energy );
        outTree->Branch("TGenMuonSt3_0_Pt",             &genMuonSt3_0_Pt );
        outTree->Branch("TGenMuonSt3_0_Eta",            &genMuonSt3_0_Eta );
        outTree->Branch("TGenMuonSt3_0_Phi",            &genMuonSt3_0_Phi );
        
        outTree->Branch("TGenElecSt3_0_PDGID",          &genElecSt3_0_PID );
        outTree->Branch("TGenElecSt3_0_FirstMom",       &genElecSt3_0_FirstMom );
        outTree->Branch("TGenElecSt3_0_Index",          &genElecSt3_0_Index );
        outTree->Branch("TGenElecSt3_0_Energy",         &genElecSt3_0_Energy );
        outTree->Branch("TGenElecSt3_0_Pt",             &genElecSt3_0_Pt );
        outTree->Branch("TGenElecSt3_0_Eta",            &genElecSt3_0_Eta );
        outTree->Branch("TGenElecSt3_0_Phi",            &genElecSt3_0_Phi );
        
        /// Status 3 sub-lead particle branches
        outTree->Branch("TGenTopSt3_1_PDGID",           &genTopSt3_1_PID );
        outTree->Branch("TGenTopSt3_1_FirstMom",        &genTopSt3_1_FirstMom );
        outTree->Branch("TGenTopSt3_1_Index",           &genTopSt3_1_Index );
        outTree->Branch("TGenTopSt3_1_Energy",          &genTopSt3_1_Energy );
        outTree->Branch("TGenTopSt3_1_Pt",              &genTopSt3_1_Pt );
        outTree->Branch("TGenTopSt3_1_Eta",             &genTopSt3_1_Eta );
        outTree->Branch("TGenTopSt3_1_Phi",             &genTopSt3_1_Phi );
        
        outTree->Branch("TGenBSt3_1_PDGID",             &genBSt3_1_PID );
        outTree->Branch("TGenBSt3_1_FirstMom",          &genBSt3_1_FirstMom );
        outTree->Branch("TGenBSt3_1_Index",             &genBSt3_1_Index );
        outTree->Branch("TGenBSt3_1_Energy",            &genBSt3_1_Energy );
        outTree->Branch("TGenBSt3_1_Pt",                &genBSt3_1_Pt );
        outTree->Branch("TGenBSt3_1_Eta",               &genBSt3_1_Eta );
        outTree->Branch("TGenBSt3_1_Phi",               &genBSt3_1_Phi );
        
        outTree->Branch("TGenMuonSt3_1_PDGID",          &genMuonSt3_1_PID );
        outTree->Branch("TGenMuonSt3_1_FirstMom",       &genMuonSt3_1_FirstMom );
        outTree->Branch("TGenMuonSt3_1_Index",          &genMuonSt3_1_Index );
        outTree->Branch("TGenMuonSt3_1_Energy",         &genMuonSt3_1_Energy );
        outTree->Branch("TGenMuonSt3_1_Pt",             &genMuonSt3_1_Pt );
        outTree->Branch("TGenMuonSt3_1_Eta",            &genMuonSt3_1_Eta );
        outTree->Branch("TGenMuonSt3_1_Phi",            &genMuonSt3_1_Phi );
        
        outTree->Branch("TGenElecSt3_1_PDGID",          &genElecSt3_1_PID );
        outTree->Branch("TGenElecSt3_1_FirstMom",       &genElecSt3_1_FirstMom );
        outTree->Branch("TGenElecSt3_1_Index",          &genElecSt3_1_Index );
        outTree->Branch("TGenElecSt3_1_Energy",         &genElecSt3_1_Energy );
        outTree->Branch("TGenElecSt3_1_Pt",             &genElecSt3_1_Pt );
        outTree->Branch("TGenElecSt3_1_Eta",            &genElecSt3_1_Eta );
        outTree->Branch("TGenElecSt3_1_Phi",            &genElecSt3_1_Phi );
        
        /// Status 1 lead particle branches
        outTree->Branch( "T_Gen_BSt1_0_PID",            &genBSt1_0_PID );
        outTree->Branch( "T_Gen_BSt1_0_Px",             &genBSt1_0_Px );
        outTree->Branch( "T_Gen_BSt1_0_Py",             &genBSt1_0_Py );
        outTree->Branch( "T_Gen_BSt1_0_Pz",             &genBSt1_0_Pz );
        outTree->Branch( "T_Gen_BSt1_0_Energy",         &genBSt1_0_Energy );
        outTree->Branch( "T_Gen_BSt1_0_MomPID",         &genBSt1_0_MomPID );
        outTree->Branch( "T_Gen_BSt1_0_MomPx",          &genBSt1_0_MomPx );
        outTree->Branch( "T_Gen_BSt1_0_MomPy",          &genBSt1_0_MomPy );
        outTree->Branch( "T_Gen_BSt1_0_MomPz",          &genBSt1_0_MomPz );
        outTree->Branch( "T_Gen_BSt1_0_MomEnergy",      &genBSt1_0_MomEnergy );
        outTree->Branch( "T_Gen_BSt1_0_MomStat",        &genBSt1_0_MomStatus );
        
        outTree->Branch( "T_Gen_MuonSt1_0_PID",         &genMuonSt1_0_PID );
        outTree->Branch( "T_Gen_MuonSt1_0_Px",          &genMuonSt1_0_Px );
        outTree->Branch( "T_Gen_MuonSt1_0_Py",          &genMuonSt1_0_Py );
        outTree->Branch( "T_Gen_MuonSt1_0_Pz",          &genMuonSt1_0_Pz );
        outTree->Branch( "T_Gen_MuonSt1_0_Energy",      &genMuonSt1_0_Energy );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPID",      &genMuonSt1_0_MomPID );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPx",       &genMuonSt1_0_MomPx );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPy",       &genMuonSt1_0_MomPy );
        outTree->Branch( "T_Gen_MuonSt1_0_MomPz",       &genMuonSt1_0_MomPz );
        outTree->Branch( "T_Gen_MuonSt1_0_MomEnergy",   &genMuonSt1_0_MomEnergy );
        outTree->Branch( "T_Gen_MuonSt1_0_MomStat",     &genMuonSt1_0_MomStatus );
        
        outTree->Branch( "T_Gen_ElecSt1_0_PID",         &genElecSt1_0_PID );
        outTree->Branch( "T_Gen_ElecSt1_0_Px",          &genElecSt1_0_Px );
        outTree->Branch( "T_Gen_ElecSt1_0_Py",          &genElecSt1_0_Py );
        outTree->Branch( "T_Gen_ElecSt1_0_Pz",          &genElecSt1_0_Pz );
        outTree->Branch( "T_Gen_ElecSt1_0_Energy",      &genElecSt1_0_Energy );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPID",      &genElecSt1_0_MomPID );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPx",       &genElecSt1_0_MomPx );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPy",       &genElecSt1_0_MomPy );
        outTree->Branch( "T_Gen_ElecSt1_0_MomPz",       &genElecSt1_0_MomPz );
        outTree->Branch( "T_Gen_ElecSt1_0_MomEnergy",   &genElecSt1_0_MomEnergy );
        outTree->Branch( "T_Gen_ElecSt1_0_MomStat",     &genElecSt1_0_MomStatus );  
        
        /// Status 1 sub-lead particle branches
        outTree->Branch( "T_Gen_BSt1_1_PID",            &genBSt1_1_PID );
        outTree->Branch( "T_Gen_BSt1_1_Px",             &genBSt1_1_Px );
        outTree->Branch( "T_Gen_BSt1_1_Py",             &genBSt1_1_Py );
        outTree->Branch( "T_Gen_BSt1_1_Pz",             &genBSt1_1_Pz );
        outTree->Branch( "T_Gen_BSt1_1_Energy",         &genBSt1_1_Energy );
        outTree->Branch( "T_Gen_BSt1_1_MomPID",         &genBSt1_1_MomPID );
        outTree->Branch( "T_Gen_BSt1_1_MomPx",          &genBSt1_1_MomPx );
        outTree->Branch( "T_Gen_BSt1_1_MomPy",          &genBSt1_1_MomPy );
        outTree->Branch( "T_Gen_BSt1_1_MomPz",          &genBSt1_1_MomPz );
        outTree->Branch( "T_Gen_BSt1_1_MomEnergy",      &genBSt1_1_MomEnergy );
        outTree->Branch( "T_Gen_BSt1_1_MomStat",        &genBSt1_1_MomStatus );
        
        outTree->Branch( "T_Gen_MuonSt1_1_PID",         &genMuonSt1_1_PID );
        outTree->Branch( "T_Gen_MuonSt1_1_Px",          &genMuonSt1_1_Px );
        outTree->Branch( "T_Gen_MuonSt1_1_Py",          &genMuonSt1_1_Py );
        outTree->Branch( "T_Gen_MuonSt1_1_Pz",          &genMuonSt1_1_Pz );
        outTree->Branch( "T_Gen_MuonSt1_1_Energy",      &genMuonSt1_1_Energy );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPID",      &genMuonSt1_1_MomPID );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPx",       &genMuonSt1_1_MomPx );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPy",       &genMuonSt1_1_MomPy );
        outTree->Branch( "T_Gen_MuonSt1_1_MomPz",       &genMuonSt1_1_MomPz );
        outTree->Branch( "T_Gen_MuonSt1_1_MomEnergy",   &genMuonSt1_1_MomEnergy );
        outTree->Branch( "T_Gen_MuonSt1_1_MomStat",     &genMuonSt1_1_MomStatus );
        
        outTree->Branch( "T_Gen_ElecSt1_1_PID",         &genElecSt1_1_PID );
        outTree->Branch( "T_Gen_ElecSt1_1_Px",          &genElecSt1_1_Px );
        outTree->Branch( "T_Gen_ElecSt1_1_Py",          &genElecSt1_1_Py );
        outTree->Branch( "T_Gen_ElecSt1_1_Pz",          &genElecSt1_1_Pz );
        outTree->Branch( "T_Gen_ElecSt1_1_Energy",      &genElecSt1_1_Energy );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPID",      &genElecSt1_1_MomPID );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPx",       &genElecSt1_1_MomPx );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPy",       &genElecSt1_1_MomPy );
        outTree->Branch( "T_Gen_ElecSt1_1_MomPz",       &genElecSt1_1_MomPz );
        outTree->Branch( "T_Gen_ElecSt1_1_MomEnergy",   &genElecSt1_1_MomEnergy );
        outTree->Branch( "T_Gen_ElecSt1_1_MomStat",     &genElecSt1_1_MomStatus );
    }
    
    cout << "--- Processing: " << fileTree.GetEntries() << " events" << endl;
    h_eventCount->Fill(1);
    h_eventCount->SetEntries(fileTree.GetEntries());
    outputFile->cd();    
    
    vector<PFJet> * Jets;     // Will contain all jets not overlapping with isolated electrons or muons
    vector<Lepton> * vecIsoLeptons;
    TLorentzVector patsyVec;
    
    vector<PFJet> * Jets_JetESUp, * Jets_JetESDown;
    
    vector<Lepton> * vecIsoLeptons_LepESUp, * vecIsoLeptons_LepESDown, * vecIsoLeptonsCentValMETPatsy_LepESUp, * vecIsoLeptonsCentValMETPatsy_LepESDown;
    
    ///####Smeared Jet Vector####//////
    vector<PFJet> * SmearJets;
    vector<PFJet> * SmearJets_JetESUp, * SmearJets_JetESDown;
    vector<PFJet> * SmearJets_JetSmearUp, * SmearJets_JetSmearDown;
    
    vector<GenJet> * vecGoodGenJets;
    
    float roundNum = 1.0;
    int roundMult = 1;
    if (doSignal) {
        roundNum = (fInName.Contains("to") || fInName.Contains("FineBin")) ? 10.0 : 1.0;
        roundMult = (fInName.Contains("to") || fInName.Contains("FineBin")) ? 10 : 1;
        if (fInName.Contains("T2bw")) {
            roundNum = 25.0;
            roundMult = 25; //Check this for non FineBin T2bw samples            
        }
    }
    /////Iterate over events  
    //    for (Long64_t ievt=0; ievt < 10;ievt++) {
    for (Long64_t ievt=0; ievt<fileTree.GetEntries();ievt++) {
        if (ievt%10000 == 0) cout << ievt << endl;
        //    for (Long64_t ievt=0; ievt<100;ievt++) {
        vecIsoLeptons = new vector<Lepton>;
        Jets = new vector<PFJet>;
        doEvent = true;        
        if (!doData) {    
            Jets_JetESUp = new vector<PFJet>;
            Jets_JetESDown = new vector<PFJet>;
            
            SmearJets = new vector<PFJet>;
            SmearJets_JetESUp = new vector<PFJet>;
            SmearJets_JetESDown = new vector<PFJet>;
            SmearJets_JetSmearUp = new vector<PFJet>;
            SmearJets_JetSmearDown = new vector<PFJet>;
            
            vecIsoLeptonsCentValMETPatsy_LepESUp = new vector<Lepton>;
            vecIsoLeptonsCentValMETPatsy_LepESDown = new vector<Lepton>;
            vecIsoLeptons_LepESUp = new vector<Lepton>;
            vecIsoLeptons_LepESDown = new vector<Lepton>;
            
            doEvent_LepESUp = true;
            doEvent_LepESDown = true;
        }
        
        
        //        map<string, float> stringKeyToVar;
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
                // Check this for T2bw samples
            }
            if (doMassCut) {
                if (genStopMass0 >=0 && abs(genStopMassCut - genStopMass0) < 2.5) continue;
                if (genChi0Mass0 >=0 && abs(genChi0MassCut - genChi0Mass0) < 2.5) continue;
                if (genCharginoMass0 >=0 && abs(genCharginoMassCut - genCharginoMass0) < 2.5) continue;
            }
        }
        if (printEventNum && doSignal) {
            cout << "in format EventNum:genStopMass:genChi0Mass:genCharginoMass ";
            cout << EventNum << ":" << genStopMass0 << ":" << genChi0Mass0 << ":" <<  genCharginoMass0 << endl;
        }
        else if (printEventNum) {
            cout << "Event Prior to Selections: RunNum:EventNum:LumiBlock " << RunNum  << ":" << EventNum << ":" << LumiBlock << endl;
        }
        //        cout << " test " << endl;
        if (whichNTupleType == 0) {
            firstGoodVertZ = goodVertexSelection(VertZ, VertRho, VertNDOF,VertIsFake, nVtx);
            if (nVtx < 1) {
                cout << "failed vertex cut " << endl; 
                continue;
            }
            weight = 1.;                        
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
        }
        else {
            weight = (float) GenWeight;                      
        }
        h_CutFlow->Fill(1.);
        if (!doData) {
            h_CutFlow_LepESUp->Fill(1);
            h_CutFlow_LepESDown->Fill(1);
        }
        if (!doEvent) continue;
        if (doSpecRun) {
            if (RunNum == (unsigned int) whichRun) {
                cout << "Pre LepCut Run:Event:Lumi" << RunNum << ":" << EventNum << ":" << LumiBlock << endl;
            }
            else {
                continue;
            }
        }
        else if (doSpecRunEvent) {
            if (!(RunNum == (unsigned int) whichRun &&EventNum == (unsigned int) whichEvent)) {
                continue;
            }
            else {
                cout << "Pre LepCut Run:Event:Lumi" << RunNum << ":" << EventNum << ":" << LumiBlock << endl;
            }
        }
        h_CutFlow->Fill(2.);
        if (!doData) {
            h_CutFlow_LepESUp->Fill(2.);
            h_CutFlow_LepESDown->Fill(2.);
        }
        if (whichNTupleType == 0) {
            EEPs.numElectrons = EEPs.ElecPx->size();
            MEPs.numMuons = MEPs.MuonPx->size();
            PFJEPs.numPFJets = PFJEPs.JetPx->size();
            //            cout << "size of vecIsoLeptons pre electron pick " << vecIsoLeptons->size() << endl;
//            ElectronPickOvi(ElecPx, ElecPy, ElecPz, ElecE, ElecCharge, PFElecPt, ElecSCEta, ElecDeltaPhiIn, ElecDeltaEtaIn, ElecSigIetaIeta, ElecHtoERatio, ElecIP, ElecDZ, ElecECalE, ElecSCEOverP, ElecNumMissHits, ElecPFNeutHadIso, ElecPFCharHadIso, ElecPFPhotIso, passConvVeto, isPFElectron, ElecisEB, ElecisEE, eventRhoIso, vecIsoLeptons, levelVerbosity);
            ElectronPickOvi(EEPs, eventRhoIso, vecIsoLeptons, levelVerbosity);
            //            cout << "size of vecIsoLeptons pre muon/post elec pick " << vecIsoLeptons->size() << endl;
//            MuonPickOvi(MuonPt, MuonPx, MuonPy, MuonPz, MuonE, MuonCharge, PFMuonPt, MuonD0, MuonVertZ, firstGoodVertZ, MuonPFNeutHadIso, MuonPFCharHadIso, MuonPFPhotIso, MuonSumPUPt, MuonNumLayers, MuonNumMatchStations, MuonNumValidPixHitsinTrack,  isGMPTMuon, isPFMuon, isGlobMuon, vecIsoLeptons, levelVerbosity);
            MuonPickOvi(MEPs, firstGoodVertZ, vecIsoLeptons, levelVerbosity);
            sort(vecIsoLeptons->begin(), vecIsoLeptons->end(), greater<Lepton>());
            ELI = LeptonPair(vecIsoLeptons, levelVerbosity);
            doEvent = ELI.doEvent;
            //            cout << "size of vecIsoLeptons post muon/electron pick " << vecIsoLeptons->size() << endl;
            if (!doData) {
//                ElectronPickOviSyst(ElecPx, ElecPy, ElecPz, ElecE, ElecCharge, PFElecPt, ElecSCEta, ElecDeltaPhiIn, ElecDeltaEtaIn, ElecSigIetaIeta, ElecHtoERatio, ElecIP, ElecDZ, ElecECalE, ElecSCEOverP, ElecNumMissHits, ElecPFNeutHadIso, ElecPFCharHadIso, ElecPFPhotIso, passConvVeto, isPFElectron, ElecisEB, ElecisEE, eventRhoIso, 1., vecIsoLeptons_LepESUp, levelVerbosity, vecIsoLeptonsCentValMETPatsy_LepESUp);                
                ElectronPickOviSyst(EEPs, eventRhoIso, 1., vecIsoLeptons_LepESUp, levelVerbosity, vecIsoLeptonsCentValMETPatsy_LepESUp);                
//                MuonPickOviSyst(MuonPt, MuonPx, MuonPy, MuonPz, MuonE, MuonCharge, PFMuonPt, MuonD0, MuonVertZ, firstGoodVertZ, MuonPFNeutHadIso, MuonPFCharHadIso, MuonPFPhotIso, MuonSumPUPt, MuonNumLayers, MuonNumMatchStations, MuonNumValidPixHitsinTrack, isGMPTMuon, isPFMuon, isGlobMuon, 1., vecIsoLeptons_LepESUp, vecIsoLeptonsCentValMETPatsy_LepESUp);
                MuonPickOviSyst(MEPs, firstGoodVertZ, 1., vecIsoLeptons_LepESUp, vecIsoLeptonsCentValMETPatsy_LepESUp);
                sort(vecIsoLeptons_LepESUp->begin(), vecIsoLeptons_LepESUp->end(), greater<Lepton>());
                ELI_LepESUp = LeptonPair(vecIsoLeptons_LepESUp, levelVerbosity);
                MET_LepESUp = MET; MET_Phi_LepESUp = MET_Phi;
                METSystShift(vecIsoLeptonsCentValMETPatsy_LepESUp, vecIsoLeptons_LepESUp, MET_LepESUp, MET_Phi_LepESUp, MET, MET_Phi);
                doEvent_LepESUp = ELI_LepESUp.doEvent;
                
//                ElectronPickOviSyst(ElecPx, ElecPy, ElecPz, ElecE, ElecCharge, PFElecPt, ElecSCEta, ElecDeltaPhiIn, ElecDeltaEtaIn, ElecSigIetaIeta, ElecHtoERatio, ElecIP, ElecDZ, ElecECalE, ElecSCEOverP, ElecNumMissHits, ElecPFNeutHadIso, ElecPFCharHadIso, ElecPFPhotIso, passConvVeto, isPFElectron, ElecisEB, ElecisEE, eventRhoIso, -1., vecIsoLeptons_LepESDown, levelVerbosity, vecIsoLeptonsCentValMETPatsy_LepESDown);
                ElectronPickOviSyst(EEPs, eventRhoIso, -1., vecIsoLeptons_LepESUp, levelVerbosity, vecIsoLeptonsCentValMETPatsy_LepESUp);                
//                MuonPickOviSyst(MuonPt, MuonPx, MuonPy, MuonPz, MuonE, MuonCharge, PFMuonPt, MuonD0, MuonVertZ, firstGoodVertZ, MuonPFNeutHadIso, MuonPFCharHadIso, MuonPFPhotIso, MuonSumPUPt, MuonNumLayers, MuonNumMatchStations, MuonNumValidPixHitsinTrack, isGMPTMuon, isPFMuon, isGlobMuon, -1., vecIsoLeptons_LepESDown, vecIsoLeptonsCentValMETPatsy_LepESDown);
                MuonPickOviSyst(MEPs, firstGoodVertZ, -1., vecIsoLeptons_LepESUp, vecIsoLeptonsCentValMETPatsy_LepESUp);
                sort(vecIsoLeptons_LepESDown->begin(), vecIsoLeptons_LepESDown->end(), greater<Lepton>());
                ELI_LepESDown = LeptonPair(vecIsoLeptons_LepESDown, levelVerbosity);
                MET_LepESDown = MET; MET_Phi_LepESDown = MET_Phi;
                METSystShift(vecIsoLeptonsCentValMETPatsy_LepESDown, vecIsoLeptons_LepESDown, MET_LepESDown, MET_Phi_LepESDown, MET, MET_Phi); 
                doEvent_LepESDown = ELI_LepESDown.doEvent;
            }
        }
        else {;
            MET = met->Pt();
            MET_Phi = met->Phi();
            IsoLeptonsPickDESY(leptons, lepPdgId, lepPFIso, vecIsoLeptons);
            sort(vecIsoLeptons->begin(), vecIsoLeptons->end(), greater<Lepton>());
            ELI = LeptonPair(vecIsoLeptons, levelVerbosity);
            doEvent = ELI.doEvent;
            if (!doData) {
                //                cout << "here 3.0" << endl;
                IsoLeptonsPickDESY(leptons, lepPdgId, lepPFIso, 1., vecIsoLeptons_LepESUp, vecIsoLeptonsCentValMETPatsy_LepESUp);
                sort(vecIsoLeptons_LepESUp->begin(), vecIsoLeptons_LepESUp->end(), greater<Lepton>());
                ELI_LepESUp = LeptonPair(vecIsoLeptons_LepESUp, levelVerbosity);
                MET_LepESUp = MET; MET_Phi_LepESUp = MET_Phi;
                METSystShift(vecIsoLeptonsCentValMETPatsy_LepESUp, vecIsoLeptons_LepESUp, MET_LepESUp, MET_Phi_LepESUp, MET, MET_Phi);
                doEvent_LepESUp = ELI_LepESDown.doEvent;
                
                IsoLeptonsPickDESY(leptons, lepPdgId, lepPFIso, -1., vecIsoLeptons_LepESDown, vecIsoLeptonsCentValMETPatsy_LepESDown);
                sort(vecIsoLeptons_LepESDown->begin(), vecIsoLeptons_LepESDown->end(), greater<Lepton>());
                ELI_LepESDown = LeptonPair(vecIsoLeptons_LepESDown, levelVerbosity);
                MET_LepESDown = MET; MET_Phi_LepESDown = MET_Phi;
                METSystShift(vecIsoLeptonsCentValMETPatsy_LepESDown, vecIsoLeptons_LepESDown, MET_LepESDown, MET_Phi_LepESDown, MET, MET_Phi);
                doEvent_LepESDown = ELI_LepESDown.doEvent;
            }
        }
        if (!doEvent) {
            if (doData) { 
                continue;
            }
            else if (!doEvent_LepESUp && !doEvent_LepESDown) {
                continue;
            }
        }
        if (whichNTupleType == 0) {
//            Jets = JetInfo(vecIsoLeptons, JetPx, JetPy, JetPz, JetE, JetNHF, JetNEF, JetCHF, JetCEF, JetNDaug, JetCharMult, JetBTag, JetPartFlav,  0., h_JetESUp);            
            Jets = JetInfo(vecIsoLeptons, PFJEPs,  0., h_JetESUp);            
            sort(Jets->begin(), Jets->end(), greater<PFJet>());
            EJI = JetKinematicsCut(Jets, BTagSFUtilToUse, doData); 
        }
        else {
            Jets = JetInfoDESY(vecIsoLeptons, jets, jetBTagCSV, 0., h_JetESUp);
            sort(Jets->begin(), Jets->end(), greater<PFJet>());
            EJI = JetKinematicsCut(Jets);
        }
        if (!doData) {
            if (whichNTupleType == 0) {
                GJEPs.numGenJets = GJEPs.genJetPx->size();
                vecGoodGenJets = GenJetsNonZero(GJEPs);
                SmearJets = JetSmear(Jets, vecGoodGenJets, 0, 0.0, h_JetESUp, h_RecoJetLowPtSmearHist, h_GenJetSmearHist, JetResolutionTF1Vec, levelVerbosity);
                SmearJets_JetESUp = JetSmear(Jets, vecGoodGenJets, 1, 1.0, h_JetESUp, h_RecoJetLowPtSmearHist, h_GenJetSmearHist, JetResolutionTF1Vec, levelVerbosity);
                SmearJets_JetESDown = JetSmear(Jets, vecGoodGenJets, 1, -1.0, h_JetESDown, h_RecoJetLowPtSmearHist, h_GenJetSmearHist, JetResolutionTF1Vec, levelVerbosity);
                SmearJets_JetSmearUp = JetSmear(Jets, vecGoodGenJets, 2, 1.0, h_JetESUp, h_RecoJetLowPtSmearHist, h_GenJetSmearHist, JetResolutionTF1Vec, levelVerbosity);
                SmearJets_JetSmearDown = JetSmear(Jets, vecGoodGenJets, 2, -1.0, h_JetESDown, h_RecoJetLowPtSmearHist, h_GenJetSmearHist, JetResolutionTF1Vec, levelVerbosity);
                
                sort(SmearJets->begin(), SmearJets->end(), greater<PFJet>());
                sort(SmearJets_JetESUp->begin(), SmearJets_JetESUp->end(), greater<PFJet>());
                sort(SmearJets_JetESDown->begin(), SmearJets_JetESDown->end(), greater<PFJet>());
                sort(SmearJets_JetSmearUp->begin(), SmearJets_JetSmearUp->end(), greater<PFJet>());
                sort(SmearJets_JetSmearDown->begin(), SmearJets_JetSmearDown->end(), greater<PFJet>());

                EJISmear = JetKinematicsCut(SmearJets, BTagSFUtilToUse, doData);
                EJISmear_JetESUp = JetKinematicsCut(SmearJets_JetESUp, BTagSFUtilToUse, doData);
                EJISmear_JetESDown = JetKinematicsCut(SmearJets_JetESDown, BTagSFUtilToUse, doData);
                EJISmear_BTagSFUp = JetKinematicsCutBTagSyst(SmearJets, BTagSFUtilToUse, 1, doSignal);
                EJISmear_BTagSFDown = JetKinematicsCutBTagSyst(SmearJets, BTagSFUtilToUse, -1, doSignal);
                
                SmearMET = MET; SmearMET_Phi = MET_Phi;
                SmearMET_JetESUp = MET; SmearMET_Phi_JetESUp = MET_Phi;
                SmearMET_JetESDown = MET; SmearMET_Phi_JetESDown = MET_Phi;
                SmearMET_JetSmearUp = MET; SmearMET_Phi_JetSmearUp = MET_Phi;
                SmearMET_JetSmearDown = MET; SmearMET_Phi_JetSmearDown = MET_Phi;
                
                METSystShift(Jets, SmearJets, SmearMET, SmearMET_Phi, MET, MET_Phi);
                METSystShift(Jets, SmearJets_JetESUp, SmearMET_JetESUp, SmearMET_Phi_JetESUp, MET, MET_Phi);
                METSystShift(Jets, SmearJets_JetESDown, SmearMET_JetESDown, SmearMET_Phi_JetESDown, MET, MET_Phi);
                METSystShift(Jets, SmearJets_JetSmearUp, SmearMET_JetSmearUp, SmearMET_Phi_JetSmearUp, MET, MET_Phi);
                METSystShift(Jets, SmearJets_JetSmearDown, SmearMET_JetSmearDown, SmearMET_Phi_JetSmearDown, MET, MET_Phi);
                
                EJI_BTagSFUp = JetKinematicsCutBTagSyst(Jets, BTagSFUtilToUse, 1, doSignal);
                EJI_BTagSFDown = JetKinematicsCutBTagSyst(Jets, BTagSFUtilToUse, -1, doSignal);
                
//                Jets_JetESUp = JetInfo(vecIsoLeptons, JetPx, JetPy, JetPz, JetE, JetNHF, JetNEF, JetCHF, JetCEF, JetNDaug, JetCharMult, JetBTag, JetPartFlav, 1.0, h_JetESUp);
                Jets_JetESUp = JetInfo(vecIsoLeptons, PFJEPs, 1.0, h_JetESUp);
                sort(Jets_JetESUp->begin(), Jets_JetESUp->end(), greater<PFJet>());
                EJI_JetESUp = JetKinematicsCut(Jets_JetESUp, BTagSFUtilToUse, doData);
                
//                Jets_JetESDown = JetInfo(vecIsoLeptons, JetPx, JetPy, JetPz, JetE, JetNHF, JetNEF, JetCHF, JetCEF, JetNDaug, JetCharMult, JetBTag, JetPartFlav, -1.0, h_JetESDown);
                Jets_JetESUp = JetInfo(vecIsoLeptons, PFJEPs, -1.0, h_JetESUp);
                sort(Jets_JetESDown->begin(), Jets_JetESDown->end(), greater<PFJet>());
                EJI_JetESDown = JetKinematicsCut(Jets_JetESDown, BTagSFUtilToUse, doData);
                
            }
            else {
                Jets_JetESUp = JetInfoDESY(vecIsoLeptons, jets, jetBTagCSV, 1.0, h_JetESUp);
                sort(Jets_JetESUp->begin(), Jets_JetESUp->end(), greater<PFJet>());
                EJI_JetESUp = JetKinematicsCut(Jets_JetESUp);                
                
                Jets_JetESDown = JetInfoDESY(vecIsoLeptons, jets, jetBTagCSV, -1.0, h_JetESUp);
                sort(Jets_JetESDown->begin(), Jets_JetESDown->end(), greater<PFJet>());
                EJI_JetESDown = JetKinematicsCut(Jets_JetESDown);
            }
            MET_JetESUp = MET; MET_Phi_JetESUp = MET_Phi;
            METSystShift(Jets, Jets_JetESUp, MET_JetESUp, MET_Phi_JetESUp, MET, MET_Phi);
            
            MET_JetESDown = MET; MET_Phi_JetESDown = MET_Phi;
            METSystShift(Jets, Jets_JetESDown, MET_JetESDown, MET_Phi_JetESDown, MET, MET_Phi);
        }
        if (whichNTupleType == 0) {
            nVtxTrue = nVtx;
        }
        if (doPURW && !doData) {
            if (doHackPURW) weight = PileupRW(nVtxSFHist, nVtx);
        }
        if (ELI.EventDiLepType == -2) ELI.EventDiLepType = 2;
        if (ELI_LepESUp.EventDiLepType == -2) ELI_LepESUp.EventDiLepType = 2;
        if (ELI_LepESDown.EventDiLepType == -2) ELI_LepESDown.EventDiLepType = 2;        
        
        if (printEventNum && !doSignal) {
            if (ELI.EventDiLepType == 0 && passTrigDoubleMu) {
                cout << "MuMu -- Event Passed: RunNum:EventNum:LumiBlock " << RunNum  << ":" << EventNum << ":" << LumiBlock << endl;                
            }
            else if (ELI.EventDiLepType == 1 && passTrigDoubleEl) {
                cout << "EE -- Event Passed: RunNum:EventNum:LumiBlock " << RunNum  << ":" << EventNum << ":" << LumiBlock << endl;                
            }
            else if (ELI.EventDiLepType == 2 && passTrigElMu) {
                cout << "EMu -- Event Passed: RunNum:EventNum:LumiBlock " << RunNum  << ":" << EventNum << ":" << LumiBlock << endl;                
            }
        }        
        
        if (doData) {
            if (ELI.EventDiLepType == 0) {
                if (!(fInName.Contains("DoubleMu") || fInName.Contains("mumu_run2012"))) continue;
            }
            else if (ELI.EventDiLepType == 1) {
                if (!(fInName.Contains("DoubleEl") || fInName.Contains("ee_run2012"))) continue;
            }
            else if (ELI.EventDiLepType == 2) {
                if (!(fInName.Contains("MuEG") || fInName.Contains("emu_run2012"))) continue;
            }
        }
        if (doEvent) {
            h_CutFlow->Fill(3.);    
        }
        if (!doData) {
            if (doEvent_LepESUp) {
                h_CutFlow_LepESUp->Fill(3.);
            }
            if (doEvent_LepESDown) {
                h_CutFlow_LepESDown->Fill(3.);
            }
        }
        ////Grab Generator Info
        // grab generator MET info if there
        if (whichNTupleType == 1) {            
            if (hasMETInfo) {
                genMET = genMETVec->Pt();
                genMETPhi = genMETVec->Phi();
            }
            //grab Generator Top info if there
            if (fInName.Contains("ttbar") && hasTopInfo) {
                genTopSt3_0_Energy = genTop->E();
                genTopSt3_0_Pt = genTop->Pt();
                genTopSt3_0_Phi = genTop->Phi();
                genTopSt3_0_Eta = genTop->Eta();
                genTopSt3_0_PID = 6;
                
                genTopSt3_1_Energy = genAntiTop->E();
                genTopSt3_1_Pt = genAntiTop->Pt();
                genTopSt3_1_Phi = genAntiTop->Phi();
                genTopSt3_1_Eta = genAntiTop->Eta();
                genTopSt3_1_PID = -6;                
            }
        }        
        if (doSpecRun) {
            if (RunNum == (unsigned int) whichRun) {
                cout << "Post LepCut Run:Event:Lumi" << RunNum << ":" << EventNum << ":" << LumiBlock << endl;
            }
            else {
                continue;
            }
        }
        outTree->Fill();        
        delete vecIsoLeptons;
        delete Jets;
        
    }
    cout << "All events done" << endl;
    outputFile->cd();
    cout << "cd-ing to output directory" << endl;
    outputFile->Write();
    if (whichNTupleType == 1) {
        h_eventCount = (TH1F*) inputFile.Get("EventsBeforeSelection/weightedEvents");
    }
    h_eventCount->Write();
    h_CutFlow->Write();
    if (!doData) {
        h_CutFlow_LepESUp->Write();
        h_CutFlow_LepESDown->Write();
    }
    h_ElecCharIso->Write();
    h_ElecNeutIso->Write();
    h_ElecPhotIso->Write();
    cout << "Writing of output file done" << endl;
    outputFile->Close();
    cout << "end of code" << endl;
}