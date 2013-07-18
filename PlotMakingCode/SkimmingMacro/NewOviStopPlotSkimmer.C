
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
#include "TProfile.h"

//#include <exception>                                                                                                      
#include <sys/stat.h>

#include "TCut.h"
//#include "mt2bisect.h"
//#include "StopDict_ver2.h"
#include "../../HeaderFiles/StopFunctionDefinitions_v2.h"

//#include "PileUpMC.h"
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
    TH1F * eventCount;
    
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
    
    float genStopMass0,genStopMass1,genChi0Mass0,genChi0Mass1;
    
// Gen info
    vector<float> * genStopMass, * genChi0Mass, * genCharginoMass;
    vector<double> * genPolWeights;
    genStopMass = new vector<float>;
    genChi0Mass = new vector<float>;
    genCharginoMass = new vector<float>;
    genPolWeights = new vector<double>;  
    
    float StopMass0, StopMass1, Chi0Mass0, Chi0Mass1;
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
    
    vector<float> * genTopEn, * genTopPt, * genTopEta, * genTopPhi;
    vector<int> * genTop_i, * genTop_firstMom, * genTop_pdgId;
    genTopEn = new vector<float>;
    genTopEta = new vector<float>;
    genTopPhi = new vector<float>;
    genTopPt = new vector<float>;    
    
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
    bool doData          = 0;
    bool doVerbosity     = 0;
    int  whichNTupleType = 0; //0 IFCA Oviedo; 1 DESY
    bool doPURW          = 0;
    bool doHackPURW      = 0;
    bool doPURWOviToDESY = 0;
    /////loop over inputs
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"-i",2) == 0) fInName = TString(argv[k+1]);
        if (strncmp (argv[k],"-w",2) == 0) whichNTupleType = strtol(argv[k+1], NULL, 10);
        if (strncmp (argv[k],"doVerbosity",11) == 0) doVerbosity = 1;
        if (strncmp (argv[k],"doPURW",6) == 0) doPURW = 1;        
        if (strncmp (argv[k],"doHackPURW",10) == 0) doHackPURW = 1;
        if (strncmp (argv[k],"doPURWOviToDESY",15) == 0) doPURWOviToDESY = 1;
    }
    if (fInName.Contains("tkolberg")) {
        fOutName = fInName;
        fOutName.Replace(12, 8, "bcalvert");   
    }
    else {
        fOutName = fInName;
    }
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
    
    //generator information
    fileTree.SetBranchAddress( "T_Gen_StopMass", &genStopMass );
    fileTree.SetBranchAddress( "T_Gen_Chi0Mass", &genChi0Mass );
    fileTree.SetBranchAddress( "T_METgen_ET", &genMET );
    fileTree.SetBranchAddress( "T_METgen_Phi", &genMETPhi );
    
    outTree->Branch("TWeight",   &weight);
    outTree->Branch("TChannel",  &Type);
    outTree->Branch("TNPV",      &nVtx);
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
    
    if (fInName.Contains("FineBin") {
        outTree->Branch("TGenStopMass0", &genStopMass0);
        outTree->Branch("TGenStopMass1", &genStopMass1);
        outTree->Branch("TGenChi0Mass0", &genChi0Mass0);
        outTree->Branch("TGenChi0Mass1", &genChi0Mass1);
    }
    outTree->Branch("TGenMET", &genMET);
    outTree->Branch("TGenMETPhi", &genMETPhi);

    cout << "--- Processing: " << fileTree.GetEntries() << " events" << endl;
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
        if (fInName.Contains("FineBin")) {
            genStopMass0 = genStopMass->at(0);
            genStopMass1 = genStopMass->at(1);
            genChi0Mass0 = genChi0Mass->at(0);
            genChi0Mass1 = genChi0Mass->at(1);
        }
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
    cout << "Writing of output file done" << endl;
    outputFile->Close();
    cout << "end of code" << endl;
}
