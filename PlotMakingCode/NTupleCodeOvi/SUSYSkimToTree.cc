// 
// Package:    SUSYSkimToTree
// Class:      SUSYSkimToTree
// 
/**\class SUSYSkimToTree SUSYSkimToTree.cc TopAnalysis/SUSYSkimToTree/src/SUSYSkimToTree.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  
//         Created:  Fri Jan  9 11:55:23 CET 2009
// $Id$
//
// 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include  "WWAnalysis/AnalysisStep/interface/TriggerBitChecker.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
//#include "TrackingTools/IPTools/interface/ImpactParameterComputer.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "SUSYMCANDXPAG/Stop_reweight/Stop_TopChi0_Reweighting.C"  

#include "TH1D.h"
#include <map>
#include "TFile.h"
#include <math.h> 
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"
#include "TTree.h"


//OOOOOOOOOOOOOOOOJOOOOOOOOOOOOOOOOOOO
//if I want to add an additional jet collection I should modify this variable
#define NumberOfJetCollections 1

//#define NumberOfHLTLabels 10  

using namespace std;
using namespace reco;
using namespace edm;


//
// class declaration
//



class SUSYSkimToTree : public edm::EDAnalyzer {
public:
  explicit SUSYSkimToTree(const edm::ParameterSet&);
  ~SUSYSkimToTree();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  /*  bool passTriggerSingleMu(size_t, bool) const;
      bool passTriggerDoubleMu(size_t, bool) const;
      bool passTriggerDoubleEl(size_t, bool) const;
      bool passTriggerElMu(size_t, bool) const;
  */
  //  edm::OwnVector<reco::RecoCandidate> leps_;
  bool IsRealData;  
  string theHistosFileName;
  TFile* theHistosFile;
  // const edm::ggerNames_TriggerNames &triggerNames_;
  //  const edm::TriggerNames & triggerNames = event.triggerNames(*triggerResults);

 
  TTree *Tree;
  
  void SetJetInfo(int idx, edm::View<pat::Jet> JET, const reco::VertexCollection& vtxs, bool calojet);
  void SetJetBranchAddress(int idx, TString namecol, bool calojet);
  
  
  void LeptonicTauDecay(const reco::Candidate& tau, bool& elecdec, bool& muondec,
			int &pid, float &px, float &py, float &pz, float &energy);
  
  //Input tags
  edm::InputTag muoLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag jetPFLabel_;
  //  edm::InputTag jetPF2Label_;
  edm::InputTag genmetLabel_;
  edm::InputTag metPFLabel_;
  edm::InputTag metPFTypeILabel_;
  edm::InputTag metTCLabel_;
  edm::InputTag metChargedNeutralPFNoPULabel_;
  edm::InputTag metChargedPFNoPULabel_;
  edm::InputTag PVLabel_;
  edm::InputTag trigLabel_;
  edm::InputTag elecLabel_;     
  //  edm::InputTag pfelecLabel_;
  //  edm::InputTag pfmuoLabel_;
  //  edm::InputTag pftauLabel_;
  edm::InputTag tauLabel_;
  TriggerBitChecker singleMuData_;
  TriggerBitChecker singleElData_;
  TriggerBitChecker doubleMuData_;
  TriggerBitChecker doubleElData_;
  TriggerBitChecker muEGData_;
  TriggerBitChecker singleMuMC_;
  TriggerBitChecker singleElMC_;
  TriggerBitChecker doubleMuMC_;
  TriggerBitChecker doubleElMC_;
  TriggerBitChecker muEGMC_;

 
  edm::ParameterSet jetIdLoose_; 
  
  typedef math::XYZTLorentzVector LorentzVector; 

  //Events
  float T_Event_Rho;
  float T_Event_RhoIso;
  //  float T_Event_RhoNoPu;
  //  float T_Event_RhoIsoNoPu;
  //  float T_Event_RhoCentralNeutral;
  //  float T_Event_RhoCentralNeutralTight;

  bool T_EventF_EcalDeadCell,T_EventF_logErrorTooManyClusters,T_EventF_trackingFailure,T_EventF_hcalLaser;
  bool T_EventF_ecalLaserCorr,T_EventF_toomanystripclus,T_EventF_manystripclus,T_EventF_eeBadSc;

  int T_Event_RunNumber;
  int T_Event_EventNumber;
  int T_Event_LuminosityBlock;
  int T_Event_processID;
  //  float T_Event_KFactorHNLO;
  //float T_Event_PtHat;
  //float T_Event_HiggsPt;
  int T_Event_nPU;
  float T_Event_nTruePU;
  int T_Event_nPUm;
  int T_Event_nPUp;
  float T_Event_AveNTruePU;

  //HLT
  /*  bool T_HLT_Mu8_v1;
      bool T_HLT_Mu8_v2;
      bool T_HLT_Mu8_v3;
      bool T_HLT_Mu8_v4;
      bool T_HLT_Mu8_v5;
      bool T_HLT_Mu8_v6;
      bool T_HLT_Mu8_v7;
      bool T_HLT_Mu8_v8;
      bool T_HLT_Mu8_v9;
      bool T_HLT_Mu8_v10;
      bool T_HLT_Mu8_v11;
      bool T_HLT_Mu8_v12;
  */
  bool T_HLT_Mu8_v16;
  /*
    bool T_HLT_Mu12_v1;
    bool T_HLT_Mu12_v2;
    bool T_HLT_Mu12_v3;
    bool T_HLT_Mu12_v4;
    bool T_HLT_Mu12_v16;

    bool T_HLT_IsoMu24_eta2p1_v11;
    bool T_HLT_IsoMu24_eta2p1_v12;

    bool T_HLT_Mu15_v1;
    bool T_HLT_Mu15_v2;
    bool T_HLT_Mu15_v3;
    bool T_HLT_Mu15_v4;
    bool T_HLT_Mu15_v5;
    bool T_HLT_Mu15_v6;
    bool T_HLT_Mu15_v7;
    bool T_HLT_Mu15_v8;
    bool T_HLT_Mu15_v9;
    bool T_HLT_Mu15_v10;
    bool T_HLT_Mu15_v11;
    bool T_HLT_Mu15_v12;
    bool T_HLT_Mu15_v13;
  */
  bool T_HLT_Mu17_v3;
  /*
    bool T_HLT_Mu24_v1;
    bool T_HLT_Mu24_v2;
    bool T_HLT_Mu24_v3;
    bool T_HLT_Mu24_v4;
    bool T_HLT_Mu24_v5;
    bool T_HLT_Mu24_v6;
    bool T_HLT_Mu24_v7;
    bool T_HLT_Mu24_v8;
    bool T_HLT_Mu24_v9;
    bool T_HLT_Mu24_v10;
    bool T_HLT_Mu24_v11;
    bool T_HLT_Mu24_v12;

    bool T_HLT_Mu9;

    bool T_HLT_Jet30_v1;
    bool T_HLT_Jet30_v2;
    bool T_HLT_Jet30_v3;
    bool T_HLT_Jet30_v4;
    bool T_HLT_Jet30_v5;

    bool T_HLT_Jet60_v1;
    bool T_HLT_Jet60_v2;
    bool T_HLT_Jet60_v3;
    bool T_HLT_Jet60_v4;
    bool T_HLT_Jet60_v5;
  */

  bool T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12;
  bool T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13;
  bool T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14;

  bool T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3;
  bool T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4;
  bool T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5;


  //Gen

  std::vector<float> *T_Gen_StopMass;
  std::vector<float> *T_Gen_Chi0Mass;
  std::vector<float> *T_Gen_CharginoMass;

  std::vector<double> *T_Gen_polWeights;

  std::vector<int> *T_Gen_Muon_PID;
  std::vector<float> *T_Gen_Muon_Px;
  std::vector<float> *T_Gen_Muon_Py;
  std::vector<float> *T_Gen_Muon_Pz;
  std::vector<float> *T_Gen_Muon_Energy;
  
  std::vector<int> *T_Gen_Elec_PID;
  std::vector<float> *T_Gen_Elec_Px;
  std::vector<float> *T_Gen_Elec_Py;
  std::vector<float> *T_Gen_Elec_Pz;
  std::vector<float> *T_Gen_Elec_Energy;
  
  std::vector<int> *T_Gen_b_PID;
  std::vector<float> *T_Gen_b_Px;
  std::vector<float> *T_Gen_b_Py;
  std::vector<float> *T_Gen_b_Pz;
  std::vector<float> *T_Gen_b_Energy;

  std::vector<int> *T_Gen_Muon_MPID;
  std::vector<float> *T_Gen_Muon_MPx;
  std::vector<float> *T_Gen_Muon_MPy;
  std::vector<float> *T_Gen_Muon_MPz;
  std::vector<float> *T_Gen_Muon_MEnergy;
  std::vector<int> *T_Gen_Muon_MSt;

  std::vector<int> *T_Gen_Elec_MPID;
  std::vector<float> *T_Gen_Elec_MPx;
  std::vector<float> *T_Gen_Elec_MPy;
  std::vector<float> *T_Gen_Elec_MPz;
  std::vector<float> *T_Gen_Elec_MEnergy;
  std::vector<int> *T_Gen_Elec_MSt;

  std::vector<int> *T_Gen_b_MPID;
  std::vector<float> *T_Gen_b_MPx;
  std::vector<float> *T_Gen_b_MPy;
  std::vector<float> *T_Gen_b_MPz;
  std::vector<float> *T_Gen_b_MEnergy;
  std::vector<int> *T_Gen_b_MSt;


  /*
  std::vector<SUSYGenParticle> *T_Gen_StopSt3; 
  std::vector<SUSYGenParticle> *T_Gen_Chi0St3; 
  std::vector<SUSYGenParticle> *T_Gen_tSt3;    
  std::vector<SUSYGenParticle> *T_Gen_ChiPMSt3;
  std::vector<SUSYGenParticle> *T_Gen_bSt3;    
  std::vector<SUSYGenParticle> *T_Gen_WSt3;    
  std::vector<SUSYGenParticle> *T_Gen_MuonSt3; 
  std::vector<SUSYGenParticle> *T_Gen_ElecSt3; 
  std::vector<SUSYGenParticle> *T_Gen_TauSt3;  
  */

  std::vector<int>   *T_Gen_StopSt3_pdgId;	
  std::vector<int>   *T_Gen_StopSt3_firstMother;
  std::vector<int>   *T_Gen_StopSt3_i;
  std::vector<float> *T_Gen_StopSt3_energy;	
  std::vector<float> *T_Gen_StopSt3_pt;		
  std::vector<float> *T_Gen_StopSt3_eta;	
  std::vector<float> *T_Gen_StopSt3_phi;        
  
  std::vector<int>   *T_Gen_Chi0St3_pdgId;
  std::vector<int>   *T_Gen_Chi0St3_firstMother;
  std::vector<int>   *T_Gen_Chi0St3_i;
  std::vector<float> *T_Gen_Chi0St3_energy;
  std::vector<float> *T_Gen_Chi0St3_pt;
  std::vector<float> *T_Gen_Chi0St3_eta;
  std::vector<float> *T_Gen_Chi0St3_phi;
  
  std::vector<int>   *T_Gen_tSt3_pdgId;
  std::vector<int>   *T_Gen_tSt3_firstMother;
  std::vector<int>   *T_Gen_tSt3_i;
  std::vector<float> *T_Gen_tSt3_energy;
  std::vector<float> *T_Gen_tSt3_pt;
  std::vector<float> *T_Gen_tSt3_eta;
  std::vector<float> *T_Gen_tSt3_phi;

  std::vector<int>   *T_Gen_ChiPMSt3_pdgId;
  std::vector<int>   *T_Gen_ChiPMSt3_firstMother;
  std::vector<int>   *T_Gen_ChiPMSt3_i;
  std::vector<float> *T_Gen_ChiPMSt3_energy;
  std::vector<float> *T_Gen_ChiPMSt3_pt;
  std::vector<float> *T_Gen_ChiPMSt3_eta;
  std::vector<float> *T_Gen_ChiPMSt3_phi;

  std::vector<int>   *T_Gen_bSt3_pdgId;
  std::vector<int>   *T_Gen_bSt3_firstMother;
  std::vector<int>   *T_Gen_bSt3_i;
  std::vector<float> *T_Gen_bSt3_energy;
  std::vector<float> *T_Gen_bSt3_pt;
  std::vector<float> *T_Gen_bSt3_eta;
  std::vector<float> *T_Gen_bSt3_phi;

  std::vector<int>   *T_Gen_WSt3_pdgId;
  std::vector<int>   *T_Gen_WSt3_firstMother;
  std::vector<int>   *T_Gen_WSt3_i;
  std::vector<float> *T_Gen_WSt3_energy;
  std::vector<float> *T_Gen_WSt3_pt;
  std::vector<float> *T_Gen_WSt3_eta;
  std::vector<float> *T_Gen_WSt3_phi;

  std::vector<int>   *T_Gen_MuonSt3_pdgId;
  std::vector<int>   *T_Gen_MuonSt3_firstMother;
  std::vector<int>   *T_Gen_MuonSt3_i;
  std::vector<float> *T_Gen_MuonSt3_energy;
  std::vector<float> *T_Gen_MuonSt3_pt;
  std::vector<float> *T_Gen_MuonSt3_eta;
  std::vector<float> *T_Gen_MuonSt3_phi;

  std::vector<int>   *T_Gen_ElecSt3_pdgId;
  std::vector<int>   *T_Gen_ElecSt3_firstMother;
  std::vector<int>   *T_Gen_ElecSt3_i;
  std::vector<float> *T_Gen_ElecSt3_energy;
  std::vector<float> *T_Gen_ElecSt3_pt;
  std::vector<float> *T_Gen_ElecSt3_eta;
  std::vector<float> *T_Gen_ElecSt3_phi;

  std::vector<int>   *T_Gen_TauSt3_pdgId;
  std::vector<int>   *T_Gen_TauSt3_firstMother;
  std::vector<int>   *T_Gen_TauSt3_i;
  std::vector<float> *T_Gen_TauSt3_energy;
  std::vector<float> *T_Gen_TauSt3_pt;
  std::vector<float> *T_Gen_TauSt3_eta;
  std::vector<float> *T_Gen_TauSt3_phi;
  




  /*
  std::vector<int> *T_Gen_ElecSt3_PID;
  std::vector<float> *T_Gen_ElecSt3_Px;
  std::vector<float> *T_Gen_ElecSt3_Py;
  std::vector<float> *T_Gen_ElecSt3_Pz;
  std::vector<float> *T_Gen_ElecSt3_Energy;
  std::vector<int> *T_Gen_bSt3_PID;
  std::vector<float> *T_Gen_bSt3_Px;
  std::vector<float> *T_Gen_bSt3_Py;
  std::vector<float> *T_Gen_bSt3_Pz;
  std::vector<float> *T_Gen_bSt3_Energy;
  std::vector<int> *T_Gen_tSt3_PID;
  std::vector<float> *T_Gen_tSt3_Px;
  std::vector<float> *T_Gen_tSt3_Py;
  std::vector<float> *T_Gen_tSt3_Pz;
  std::vector<float> *T_Gen_tSt3_Energy;
  
  std::vector<int> *T_Gen_TauSt3_PID;
  std::vector<float> *T_Gen_TauSt3_Px;
  std::vector<float> *T_Gen_TauSt3_Py;
  std::vector<float> *T_Gen_TauSt3_Pz;
  std::vector<float> *T_Gen_TauSt3_Energy;
  */
  std::vector<bool> *T_Gen_TauSt3_IsLepDec;
  std::vector<int> *T_Gen_TauSt3_LepDec_PID;
  std::vector<float> *T_Gen_TauSt3_LepDec_Px;
  std::vector<float> *T_Gen_TauSt3_LepDec_Py;
  std::vector<float> *T_Gen_TauSt3_LepDec_Pz;
  std::vector<float> *T_Gen_TauSt3_LepDec_Energy;
 
  //PFcandidates
  /*  std::vector<float> *T_PFreducedCand_Px;
      std::vector<float> *T_PFreducedCand_Py;
      std::vector<float> *T_PFreducedCand_Pz;
      std::vector<float> *T_PFreducedCand_E;
      std::vector<int>   *T_PFreducedCand_ID;
      std::vector<float> *T_PFreducedCand_vz;
  */  
  //Muons
  std::vector<float>*T_Muon_Eta;
  std::vector<bool> *T_Muon_IsGlobalMuon;
  //  std::vector<bool> *T_Muon_IsTMLSOLPT;
  //  std::vector<bool> *T_Muon_IsTMLSOLPL;
  std::vector<bool> *T_Muon_IsGMPTMuons;
  //std::vector<bool> *T_Muon_IsAllStandAloneMuons;       // checks isStandAloneMuon flag
  std::vector<bool> *T_Muon_IsAllTrackerMuons;          // checks isTrackerMuon flag
  std::vector<bool> *T_Muon_IsTrackerMuonArbitrated;    // resolve ambiguity of sharing segments
  std::vector<bool> *T_Muon_IsAllArbitrated;            // all muons with the tracker muon arbitrated
  /*std::vector<bool> *T_Muon_IsTMLastStationLoose;       // penetration depth loose selector
    std::vector<bool> *T_Muon_IsTMLastStationTight;       // penetration depth Tight selector
    std::vector<bool> *T_Muon_IsTM2DCompatibilityLoose;   // likelihood based loose selector
    std::vector<bool> *T_Muon_IsTM2DCompatibilityTight;   // likelihood based tight selector
    std::vector<bool> *T_Muon_IsTMOneStationLoose;        // require one well matched segment
    std::vector<bool> *T_Muon_IsTMOneStationTight;        // require one well matched segment
    std::vector<bool> *T_Muon_IsTMLSOPL; // combination of TMLastStation and TMOneStation
    std::vector<bool> *T_Muon_IsGMTkChiCompatibility;  // require tk stub have good chi2 relative to glb track
    std::vector<bool> *T_Muon_IsGMStaChiCompatibility;    // require sta stub have good chi2 compatibility relative to glb track
    std::vector<bool> *T_Muon_IsGMTkKinkTight;
    std::vector<bool> *T_Muon_IsTMLastStationAngLoose;
    std::vector<bool> *T_Muon_IsTMLastStationAngTight;
    std::vector<bool> *T_Muon_IsTMOneStationAngLoose;
    std::vector<bool> *T_Muon_IsTMOneStationAngTight;*/
  std::vector<float> *T_Muon_SegmentCompatibility; 
  std::vector<float> *T_Muon_trkKink;
  std::vector<float> *T_Muon_Px;
  std::vector<float> *T_Muon_Py;
  std::vector<float> *T_Muon_Pz;
  std::vector<float> *T_Muon_Pt;
  std::vector<float> *T_Muon_deltaPt;
  std::vector<float> *T_Muon_Energy;
  std::vector<int> *T_Muon_Charge;
  std::vector<float> *T_Muon_NormChi2GTrk;
  std::vector<int> *T_Muon_NValidHitsInTrk;
  std::vector<int> *T_Muon_NValidPixelHitsInTrk;
  std::vector<int> *T_Muon_NValidHitsSATrk;
  std::vector<int> *T_Muon_NValidHitsGTrk;
  std::vector<int> *T_Muon_NumOfMatchedStations;
  std::vector<int> *T_Muon_InnerTrackFound;
  std::vector<float> *T_Muon_Chi2InTrk;
  std::vector<float> *T_Muon_dofInTrk;
  /*std::vector<float> *T_Muon_SumIsoCalo;
    std::vector<float> *T_Muon_SumIsoTrack;*/
  std::vector<float> *T_Muon_IPAbsGTrack;
  //  std::vector<float> *T_Muon_IPSigGTrack;
  std::vector<float> *T_Muon_IPAbsInTrack;
  //  std::vector<float> *T_Muon_IPSigInTrack;
  //  std::vector<float> *T_Muon_IPwrtBSInTrack;
  std::vector<float> *T_Muon_IPwrtAveBSInTrack;
  //  std::vector<float> *T_Muon_IPwrtBSGTrack;
  /*  std::vector<float> *T_Muon_IP2DBiasedPV;
      std::vector<float> *T_Muon_IP3DBiasedPV;
      std::vector<float> *T_Muon_IP2DUnBiasedPV;
      std::vector<float> *T_Muon_IP3DUnBiasedPV;
      std::vector<float> *T_Muon_dxyPVBiasedPV;
      std::vector<float> *T_Muon_dzPVBiasedPV;
      std::vector<float> *T_Muon_dxyPVUnBiasedPV;
      std::vector<float> *T_Muon_dzPVUnBiasedPV;
      std::vector<float> *T_Muon_IP2DUnBiasedPVnoBS;
      std::vector<float> *T_Muon_IP3DUnBiasedPVnoBS;
      std::vector<float> *T_Muon_dxyPVUnBiasedPVnoBS;
      std::vector<float> *T_Muon_dzPVUnBiasedPVnoBS;
      std::vector<float> *T_Muon_pfCharged;
      std::vector<float> *T_Muon_pfNeutral;
      std::vector<float> *T_Muon_pfPhoton;*/
  /*std::vector<float> *T_Muon_smurfCharged;
    std::vector<float> *T_Muon_smurfNeutral;
    std::vector<float> *T_Muon_smurfPhoton;
    std::vector<float> *T_Muon_smurfNoOverCharged;
    std::vector<float> *T_Muon_smurfNoOverNeutral;
    std::vector<float> *T_Muon_smurfNoOverPhoton;
    std::vector<float> *T_Muon_muSmurfPF;*/
  //  std::vector<float> *T_Muon_chargedParticleIsoR04;
  std::vector<float> *T_Muon_chargedHadronIsoR04;
  std::vector<float> *T_Muon_neutralHadronIsoR04;
  std::vector<float> *T_Muon_photonIsoR04;
  std::vector<float> *T_Muon_sumPUPtR04;
  std::vector<float> *T_Muon_chargedParticleIsoR03;
  std::vector<float> *T_Muon_chargedHadronIsoR03;
  std::vector<float> *T_Muon_neutralHadronIsoR03;
  std::vector<float> *T_Muon_photonIsoR03;
  std::vector<float> *T_Muon_sumPUPtR03;
  std::vector<float> *T_Muon_vz;
  std::vector<float> *T_Muon_vy;
  std::vector<float> *T_Muon_vx;
  std::vector<bool> *T_Muon_isPFMuon;
  std::vector<float> *T_Muon_PFMuonPt; 
  std::vector<float> *T_Muon_PFMuonPx;
  std::vector<float> *T_Muon_PFMuonPy;
  std::vector<float> *T_Muon_PFMuonPz;
  std::vector<float> *T_Muon_PFMuonE;
  std::vector<int> *T_Muon_NLayers;

  /*  std::vector<bool>  *T_Muon_passTriggerSingleMu;
      std::vector<bool>  *T_Muon_passTriggerDoubleMu;
      std::vector<bool>  *T_Muon_passTriggerElMu;
  */  
  // Tau
  /*  std::vector<float> *T_Tau_Px;
      std::vector<float> *T_Tau_Py;
      std::vector<float> *T_Tau_Pz;
      std::vector<float> *T_Tau_Energy;
      std::vector<int> *T_Tau_Charge;
  */
  /*
  // pfTau
  std::vector<float> *T_pfTau_Px;
  std::vector<float> *T_pfTau_Py;
  std::vector<float> *T_pfTau_Pz;
  std::vector<float> *T_pfTau_Energy;
  std::vector<int> *T_pfTau_Charge;
  */

  //pfMuons
  /*
    std::vector<bool> *T_pfMuon_IsGlobalMuon;
    std::vector<bool> *T_pfMuon_IsGMPTMuons;
    std::vector<bool> *T_pfMuon_IsAllStandAloneMuons;       // checks isStandAloneMuon flag
    std::vector<bool> *T_pfMuon_IsAllTrackerMuons;          // checks isTrackerMuon flag
    std::vector<bool> *T_pfMuon_IsTMLastStationAngTight;
    std::vector<float> *T_pfMuon_SegmentCompatibility; 
    std::vector<float> *T_pfMuon_Px;
    std::vector<float> *T_pfMuon_Py;
    std::vector<float> *T_pfMuon_Pz;
    std::vector<float> *T_pfMuon_Pt;
    std::vector<float> *T_pfMuon_Energy;
    std::vector<int> *T_pfMuon_Charge;
    std::vector<float> *T_pfMuon_NormChi2GTrk;
    std::vector<int> *T_pfMuon_NValidHitsInTrk;
    std::vector<int> *T_pfMuon_NValidHitsSATrk;
    std::vector<float> *T_pfMuon_Chi2InTrk;
    std::vector<float> *T_pfMuon_dofInTrk;
    std::vector<float> *T_pfMuon_SumIsoCalo;
    std::vector<float> *T_pfMuon_SumIsoTrack;
    std::vector<float> *T_pfMuon_IPAbsGTrack;
    std::vector<float> *T_pfMuon_IPSigGTrack;
    std::vector<float> *T_pfMuon_IPAbsInTrack;
    std::vector<float> *T_pfMuon_IPSigInTrack;
    //  std::vector<float> *T_pfMuon_IPwrtBSInTrack;
    std::vector<float> *T_pfMuon_IP2DBiasedPV;
    std::vector<float> *T_pfMuon_IP3DBiasedPV;
    std::vector<float> *T_pfMuon_IP2DUnBiasedPV;
    std::vector<float> *T_pfMuon_IP3DUnBiasedPV;
    std::vector<float> *T_pfMuon_vz;
    std::vector<float> *T_pfMuon_vy;
    std::vector<float> *T_pfMuon_vx;
    std::vector<int> *T_pfMuon_NValidHits;
    std::vector<int> *T_pfMuon_NValidPixelHitsInTrk;
    std::vector<float> *T_pfMuon_photonIso;
    std::vector<float> *T_pfMuon_deltaPt;
    std::vector<int> *T_pfMuon_NumOfMatchedStations;
    std::vector<float> *T_pfMuon_pfCharged;
    std::vector<float> *T_pfMuon_pfNeutral;
    std::vector<float> *T_pfMuon_pfPhoton;
    std::vector<float> *T_pfMuon_smurfCharged;
    std::vector<float> *T_pfMuon_smurfNeutral;
    std::vector<float> *T_pfMuon_smurfPhoton;
    std::vector<float> *T_pfMuon_smurfNoOverCharged;
    std::vector<float> *T_pfMuon_smurfNoOverNeutral;
    std::vector<float> *T_pfMuon_smurfNoOverPhoton;
    std::vector<float> *T_pfMuon_muSmurfPF;
    std::vector<float> *T_pfMuon_sumPUPt; 
    std::vector<float> *T_pfMuon_sumPUPtR03;
    std::vector<float> *T_pfMuon_neutralHadronIso;
    std::vector<float> *T_pfMuon_chargedHadronIso;
    std::vector<float> *T_pfMuon_particleIso;
    std::vector<int> * T_pfMuon_NLayers; 
  */

  //Vertex
  std::vector<float> *T_Vertex_z;
  std::vector<float> *T_Vertex_y;
  std::vector<float> *T_Vertex_x;
  std::vector<float> *T_Vertex_Chi2Prob;
  std::vector<float> *T_Vertex_ndof;
  std::vector<float> *T_Vertex_rho;
  std::vector<bool> *T_Vertex_isFake;
  std::vector<int>  *T_Vertex_tracksSize;
 
  
  //Electrons
  std::vector<float> *T_Elec_Eta;
  std::vector<float> *T_Elec_IPwrtAveBS;
  std::vector<float> *T_Elec_IPwrtPV;
  std::vector<float> *T_Elec_dzwrtPV;
  std::vector<float> *T_Elec_Px;
  std::vector<float> *T_Elec_Py;
  std::vector<float> *T_Elec_Pz;
  std::vector<float> *T_Elec_Pt;
  std::vector<float> *T_Elec_Energy;
  std::vector<int> *T_Elec_Charge;
  /*  std::vector<float> *T_Elec_SumIsoCalo;
      std::vector<float> *T_Elec_SumIsoTrack;
      std::vector<float> *T_Elec_IP2DBiasedPV;
      std::vector<float> *T_Elec_IP3DBiasedPV;
      std::vector<float> *T_Elec_IP2DUnBiasedPV;
      std::vector<float> *T_Elec_IP3DUnBiasedPV;
      std::vector<float> *T_Elec_dxyPVBiasedPV;
      std::vector<float> *T_Elec_dzPVBiasedPV;
      std::vector<float> *T_Elec_dxyPVUnBiasedPV;
      std::vector<float> *T_Elec_dzPVUnBiasedPV;
      std::vector<float> *T_Elec_IP2DUnBiasedPVnoBS;
      std::vector<float> *T_Elec_IP3DUnBiasedPVnoBS;
      std::vector<float> *T_Elec_dxyPVUnBiasedPVnoBS;
      std::vector<float> *T_Elec_dzPVUnBiasedPVnoBS;*/
  std::vector<float> *T_Elec_puChargedHadronIso;
  std::vector<float> *T_Elec_chargedHadronIso;
  std::vector<float> *T_Elec_neutralHadronIso;
  std::vector<float> *T_Elec_photonIso;
  std::vector<float> *T_Elec_pfIsoEA03;

  /* std::vector<float> *T_Elec_smurfCharged;
     std::vector<float> *T_Elec_smurfNeutral;
     std::vector<float> *T_Elec_smurfPhoton;
     std::vector<float> *T_Elec_smurfNoOverCharged;
     std::vector<float> *T_Elec_smurfNoOverNeutral;
     std::vector<float> *T_Elec_smurfNoOverPhoton;
     std::vector<float> *T_Elec_eleSmurfPF;*/
  std::vector<bool> *T_Elec_passConversionVeto;
  std::vector<float> *T_Elec_vz;
  std::vector<float> *T_Elec_vy;
  std::vector<float> *T_Elec_vx;
  std::vector<int> *T_Elec_nLost; 
  std::vector<int> *T_Elec_nHits;
  std::vector<float> *T_Elec_SC_Et;
  std::vector<float> *T_Elec_SC_Eta;
  std::vector<int> *T_Elec_nBrems;
  std::vector<float> *T_Elec_fBrem;
  std::vector<float> *T_Elec_eSuperClusterOverP;
  std::vector<float> *T_Elec_ecalEnergy;
  std::vector<float> *T_Elec_dr03TkSumPt;
  std::vector<float> *T_Elec_dr03EcalSumEt;
  std::vector<float> *T_Elec_dr03HcalSumEt;
  std::vector<float> *T_Elec_ConvInfoDist;
  std::vector<float> *T_Elec_ConvInfoDCot;
  std::vector<bool> *T_Elec_isEB;
  std::vector<bool> *T_Elec_isEE;
  std::vector<float> *T_Elec_MVA;
  std::vector<bool> *T_Elec_isPF;
  std::vector<float> *T_Elec_PFElecPt;
  std::vector<float> *T_Elec_PFElecPx;
  std::vector<float> *T_Elec_PFElecPy;
  std::vector<float> *T_Elec_PFElecPz;
  std::vector<float> *T_Elec_PFElecE;
  /*  std::vector<float> *T_Elec_simpleEleId95;
      std::vector<float> *T_Elec_simpleEleId90;
      std::vector<float> *T_Elec_simpleEleId85;*/
  std::vector<float> *T_Elec_simpleEleId80;
  /*  std::vector<float> *T_Elec_simpleEleId70;
      std::vector<float> *T_Elec_simpleEleId60;*/
  //  std::vector<float> *T_Elec_simpleEleId90cIso;
  /*  std::vector<float> *T_Elec_cicVeryLooseMC;
      std::vector<float> *T_Elec_cicLooseMC;
      std::vector<float> *T_Elec_cicMediumMC;
      std::vector<float> *T_Elec_cicSuperTightMC;
      std::vector<float> *T_Elec_cicHyperTight1MC;
      std::vector<float> *T_Elec_cicHyperTight2MC;
      std::vector<float> *T_Elec_cicHyperTight3MC;
      std::vector<float> *T_Elec_cicVeryLooseHWW;
      std::vector<float> *T_Elec_cicLooseHWW;
      std::vector<float> *T_Elec_cicMediumHWW;
      std::vector<float> *T_Elec_cicSuperTightHWW;
      std::vector<float> *T_Elec_cicHyperTight1HWW;
      std::vector<float> *T_Elec_cicHyperTight2HWW;
      std::vector<float> *T_Elec_cicHyperTight3HWW;
      std::vector<float> *T_Elec_cicVeryLoose;
      std::vector<float> *T_Elec_cicLoose;
      std::vector<float> *T_Elec_cicMedium;
      std::vector<float> *T_Elec_cicSuperTight;
      std::vector<float> *T_Elec_cicHyperTight1;
      std::vector<float> *T_Elec_cicHyperTight2;
      std::vector<float> *T_Elec_cicHyperTight3;*/
  /*std::vector<float> *T_Elec_RobustLoose;
    std::vector<float> *T_Elec_Loose;
    std::vector<float> *T_Elec_Tight;
    std::vector<float> *T_Elec_RobustTight;
    std::vector<float> *T_Elec_RobustHighEnergy;
    std::vector<float> *T_Elec_egammaIDLikelihood;*/
  std::vector<float> *T_Elec_deltaPhiIn;
  std::vector<float> *T_Elec_deltaEtaIn;
  std::vector<float> *T_Elec_sigmaIetaIeta;
  std::vector<bool>   *T_Elec_isEcalDriven; 
  std::vector<float> *T_Elec_HtoE;
  /*  std::vector<bool>  *T_Elec_passTriggerDoubleEl;
      std::vector<bool>  *T_Elec_passTriggerElMu;
  */
  /*
  //pfElectrons
  std::vector<float> *T_pfElec_Px;
  std::vector<float> *T_pfElec_Py;
  std::vector<float> *T_pfElec_Pz;
  std::vector<float> *T_pfElec_Pt;
  std::vector<float> *T_pfElec_Energy;
  std::vector<int> *T_pfElec_Charge;
  std::vector<float> *T_pfElec_SumIsoCalo;
  std::vector<float> *T_pfElec_SumIsoTrack;
  std::vector<float> *T_pfElec_IP2DBiasedPV;
  std::vector<float> *T_pfElec_IP3DBiasedPV;
  std::vector<float> *T_pfElec_IP2DUnBiasedPV;
  std::vector<float> *T_pfElec_IP3DUnBiasedPV;
  std::vector<float> *T_pfElec_vz;
  std::vector<float> *T_pfElec_vy;
  std::vector<float> *T_pfElec_vx;
  std::vector<int> *T_pfElec_nBrems;
  std::vector<float> *T_pfElec_dr03TkSumPt;
  std::vector<float> *T_pfElec_dr03EcalSumEt;
  std::vector<float> *T_pfElec_dr03HcalSumEt;
  std::vector<float> *T_pfElec_SC_Et;
  std::vector<float> *T_pfElec_SC_Eta;
  std::vector<bool> *T_pfElec_isEcalDriven;
  std::vector<float> *T_pfElec_photonIso;
  std::vector<float> *T_pfElec_neutralHadronIso;
  std::vector<float> *T_pfElec_chargedHadronIso;
  std::vector<float> *T_pfElec_particleIso;
  std::vector<int> *T_pfElec_nHits;
  std::vector<float> *T_pfElec_ConvInfoDCot;
  std::vector<float> *T_pfElec_ConvInfoDist;
  std::vector<float> *T_pfElec_HtoE;
  */
  //Jets vector
  std::vector<float> *T_Jet_Px[NumberOfJetCollections];
  std::vector<float> *T_Jet_Py[NumberOfJetCollections];
  std::vector<float> *T_Jet_Pz[NumberOfJetCollections];
  std::vector<float> *T_Jet_Et[NumberOfJetCollections];
  /*  std::vector<float> *T_Jet_Px11[NumberOfJetCollections];
      std::vector<float> *T_Jet_Py11[NumberOfJetCollections];
      std::vector<float> *T_Jet_Pz11[NumberOfJetCollections];
      std::vector<float> *T_Jet_Et11[NumberOfJetCollections];*/
  //  std::vector<float> *T_Jet_EtOffset[NumberOfJetCollections];
  //  std::vector<float> *T_Jet_PtOffset[NumberOfJetCollections];
  std::vector<float> *T_Jet_Eta[NumberOfJetCollections];
  std::vector<float> *T_Jet_Energy[NumberOfJetCollections];
  //  std::vector<float> *T_Jet_Corr[NumberOfJetCollections];
  //  std::vector<float> *T_Jet_Corr11[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_HighEffTC[NumberOfJetCollections];  
  std::vector<float> *T_Jet_Tag_CombSVtx[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_CombSVtxMVA[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_TauJet[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_ImpParMVA[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_JetBProb[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_JetProb[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_HighPurSimpSVtx[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_HighEffSimpSVtx[NumberOfJetCollections];
  std::vector<float> *T_Jet_Tag_HighPurTC[NumberOfJetCollections];
  std::vector<float> *T_Jet_Parton_Px[NumberOfJetCollections];
  std::vector<float> *T_Jet_Parton_Py[NumberOfJetCollections];
  std::vector<float> *T_Jet_Parton_Pz[NumberOfJetCollections];
  std::vector<float> *T_Jet_Parton_Energy[NumberOfJetCollections];
  std::vector<int> *T_Jet_Parton_Flavour[NumberOfJetCollections];  
  
  
  std::vector<float> *T_Jet_CharHadEnergyFrac[NumberOfJetCollections];
  std::vector<float> *T_Jet_NeutHadEnergyFrac[NumberOfJetCollections];
  std::vector<float> *T_Jet_CharEmEnergyFrac[NumberOfJetCollections]; 
  std::vector<float> *T_Jet_NeutEmEnergyFrac[NumberOfJetCollections];
  std::vector<float> *T_Jet_CharHadEnergy[NumberOfJetCollections];
  std::vector<float> *T_Jet_NeutHadEnergy[NumberOfJetCollections];
  std::vector<float> *T_Jet_CharEmEnergy[NumberOfJetCollections]; 
  std::vector<float> *T_Jet_NeutEmEnergy[NumberOfJetCollections];
  std::vector<int> *T_Jet_MuonMultiplicity[NumberOfJetCollections];
  std::vector<int> *T_Jet_NeutralMultiplicity[NumberOfJetCollections];
  std::vector<int> *T_Jet_ChargedMultiplicity[NumberOfJetCollections];
  std::vector<bool> *T_Jet_IDLoose[NumberOfJetCollections];
  std::vector<int> *T_Jet_nDaughters[NumberOfJetCollections];
 
  std::vector<float> *T_Jet_GenJet_InvisibleE[NumberOfJetCollections];
  std::vector<float> *T_Jet_GenJet_Px[NumberOfJetCollections];
  std::vector<float> *T_Jet_GenJet_Py[NumberOfJetCollections];
  std::vector<float> *T_Jet_GenJet_Pz[NumberOfJetCollections];
  std::vector<float> *T_Jet_GenJet_Eta[NumberOfJetCollections];
  std::vector<float> *T_Jet_GenJet_Et[NumberOfJetCollections];
  std::vector<float> *T_Jet_GenJet_Energy[NumberOfJetCollections];
  std::vector<bool> *T_Jet_IsGenJet[NumberOfJetCollections];      

  
  //MET: MET 
  /*  float T_MET_ET;
      float T_MET_Phi;
      float T_MET_Sig;
  */  
  float T_METPF_ET;
  float T_METPF_Phi;
  float T_METPF_Sig;
  float T_METPFTypeI_ET;
  float T_METPFTypeI_Phi;
  float T_METPFTypeI_Sig;
  /*
    float T_METChargedNeutralPFNoPU_ET;
    float T_METChargedNeutralPFNoPU_Phi;

    float T_METChargedPFNoPU_ET;
    float T_METChargedPFNoPU_Phi;
  */  
  /*  float T_METtc_ET;
      float T_METtc_Phi;
      float T_METtc_Sig;
  */
  /*
    float T_MET_assocPfMet_ET;
    float T_MET_assocPfMet_Phi;

    float T_MET_trkPfMet_ET;
    float T_MET_trkPfMet_Phi;

    float T_MET_assocOtherVtxPfMet_ET;
    float T_MET_assocOtherVtxPfMet_Phi;

    float T_MET_centralPfMet_ET;
    float T_MET_centralPfMet_Phi;

    float T_MET_cleanPfMet_ET;
    float T_MET_cleanPfMet_Phi;
  */
  float T_METgen_ET;
  float T_METgen_Phi;
  
  bool T_passTriggerDoubleEl;
  bool T_passTriggerDoubleMu;
  //  bool T_passTriggerSingleEl; 
  //  bool T_passTriggerSingleMu; 
  bool T_passTriggerElMu;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SUSYSkimToTree::SUSYSkimToTree(const edm::ParameterSet& iConfig):
  theHistosFileName(iConfig.getUntrackedParameter<string>("histosFileName")),
  
  muoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  jetPFLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetPFTag")),
  //  jetPF2Label_(iConfig.getUntrackedParameter<edm::InputTag>("jetPF2Tag")),
  genmetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("genmetTag")),
  metPFLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metPFTag")),
  metPFTypeILabel_(iConfig.getUntrackedParameter<edm::InputTag>("metPFTypeITag")),
  metTCLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTCTag")),
  metChargedNeutralPFNoPULabel_(iConfig.getUntrackedParameter<edm::InputTag>("metChargedNeutralPFNoPUTag")),
  metChargedPFNoPULabel_(iConfig.getUntrackedParameter<edm::InputTag>("metChargedPFNoPUTag")),
  PVLabel_(iConfig.getUntrackedParameter<edm::InputTag>("PVTag")),
  trigLabel_(iConfig.getUntrackedParameter<edm::InputTag>("trigTag")),
  elecLabel_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
  //  pfelecLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfelectronTag")),
  //   pfmuoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pfmuonTag")),
  //  pftauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("pftauTag")),
  tauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tauTag")),
  singleMuData_  ( iConfig.getParameter<std::vector<std::string> >("singleMuDataPaths") ),
  singleElData_  ( iConfig.getParameter<std::vector<std::string> >("singleElDataPaths") ),
  doubleMuData_  ( iConfig.getParameter<std::vector<std::string> >("doubleMuDataPaths") ),
  doubleElData_  ( iConfig.getParameter<std::vector<std::string> >("doubleElDataPaths") ),
  muEGData_      ( iConfig.getParameter<std::vector<std::string> >("muEGDataPaths") ),
  singleMuMC_    ( iConfig.getParameter<std::vector<std::string> >("singleMuMCPaths") ),
  singleElMC_    ( iConfig.getParameter<std::vector<std::string> >("singleElMCPaths") ),
  doubleMuMC_    ( iConfig.getParameter<std::vector<std::string> >("doubleMuMCPaths") ),
  doubleElMC_    ( iConfig.getParameter<std::vector<std::string> >("doubleElMCPaths") ),
  muEGMC_        ( iConfig.getParameter<std::vector<std::string> >("muEGMCPaths") ) 
{}


SUSYSkimToTree::~SUSYSkimToTree()
{}



void
SUSYSkimToTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  // first: get all objects from the event.
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  
  //for the moment decide if is real data if collection of gen particles is found
  IsRealData = false;
  edm::Handle <reco::GenParticleCollection> genParticles;
  try {
    iEvent.getByLabel("prunedGen", genParticles);
    //I need to call genParticles size to forze the exception
    int aux = genParticles->size();
    //to avoid warnings I add the line aux = 0
    aux = 0+aux;
  }
  catch(...) {IsRealData = true;} 

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muoLabel_,muonHandle);

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  iEvent.getByLabel(elecLabel_, electronHandle);

  //  edm::Handle<edm::View<pat::Tau> > pftauHandle;
  //  iEvent.getByLabel(pftauLabel_,pftauHandle);
 
  /*  edm::Handle<edm::View<pat::Tau> > tauHandle;
      iEvent.getByLabel(tauLabel_,tauHandle);
  */ 
  //  edm::Handle<edm::View<pat::Muon> > pfmuonHandle;
  //  iEvent.getByLabel(pfmuoLabel_,pfmuonHandle);
  
  //  edm::Handle<edm::View<pat::Electron> > pfelectronHandle;
  //  iEvent.getByLabel(pfelecLabel_, pfelectronHandle);
      
  /*  edm::Handle<edm::View<pat::Jet> > jetHandle;
      iEvent.getByLabel(jetLabel_,jetHandle);
      edm::View<pat::Jet> jets = *jetHandle;
  */
  edm::Handle<edm::View<pat::Jet> > jetPFHandle;
  iEvent.getByLabel(jetPFLabel_,jetPFHandle);
  edm::View<pat::Jet> jetsPF = *jetPFHandle;
  /*  edm::Handle<edm::View<pat::Jet> > jetPF2Handle;
      iEvent.getByLabel(jetPF2Label_,jetPF2Handle);
      edm::View<pat::Jet> jetsPF2 = *jetPF2Handle;
  */
 
  /*  edm::Handle<reco::CaloMETCollection> metHandle;
      iEvent.getByLabel(metLabel_,metHandle);
      const CaloMET *mets=	&((metHandle.product())->front());
  */    
  //  edm::Handle<edm::View<pat::Met>> metsPF;
  //  iEvent.getByLabel(metPFLabel_,metsPF);
  //  const PFMET * metsPF= &((metPF.product())->front());

  /*  edm::Handle< std::vector<pat::MET>> metPFTypeI;
      iEvent.getByLabel(metPFLabel_,metPFTypeI);
      pat::MET  metsPFTypeI= *(metPFTypeI->begin());*/

  edm::Handle<std::vector<reco::PFMET>> metsPFTypeI ;
  iEvent.getByLabel(metPFTypeILabel_,metsPFTypeI);
  reco::PFMET metPFTypeI = *(metsPFTypeI->begin());

  edm::Handle<std::vector<reco::PFMET>>  metsPF;
  iEvent.getByLabel(metPFLabel_,metsPF);
  reco::PFMET metPF = *(metsPF->begin());
  
  /*
    edm::Handle<reco::METCollection> metTC;
    iEvent.getByLabel(metTCLabel_,metTC);
    const MET * metsTC = &((metTC.product())->front());
  */  
  
  //  edm::Handle<edm::ValueMap<reco::PFMET> > vmH;
  //  edm::Handle<reco::VertexCollection> vtxH;
  //  iEvent.getByLabel(metChargedNeutralPFNoPULabel_,vmH);
  //  iEvent.getByLabel(edm::InputTag(PVLabel_),vtxH);

  //  edm::Handle<edm::ValueMap<reco::PFMET> > vmH2;
  //  iEvent.getByLabel(metChargedPFNoPULabel_,vmH2);
  edm::Handle<VertexCollection> vertex;
  iEvent.getByLabel(PVLabel_,vertex);
  const reco::VertexCollection& vtxs = *(vertex.product());
 
  // Get beamspot for d0 determination
  //  BeamSpot beamSpot;
  //  Handle<BeamSpot> beamSpotHandle;
  //  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  //  beamSpot = *beamSpotHandle; 

  /*  //PF candidate
  // edm::Handle<reco::CandidateCollection> 
  edm::Handle<reco::CandidateView> candsH;
  iEvent.getByLabel("reducedPFCandsPfNoPU",candsH);
  T_PFreducedCand_Px=new std::vector<float>;
  T_PFreducedCand_Py=new std::vector<float>;
  T_PFreducedCand_Pz=new std::vector<float>;
  T_PFreducedCand_E =new std::vector<float>;
  T_PFreducedCand_ID =new std::vector<int>;
  T_PFreducedCand_vz=new std::vector<float>;

  for(reco::CandidateView::const_iterator it= candsH->begin(); it!=candsH->end(); ++it){
  T_PFreducedCand_Px->push_back(it->p4().px());
  T_PFreducedCand_Py->push_back(it->p4().py());
  T_PFreducedCand_Pz->push_back(it->p4().pz());
  T_PFreducedCand_E->push_back(it->p4().E());
  T_PFreducedCand_ID->push_back(it->pdgId());
  */
  /*
    double dz(99999.);
    if(it->trackRef().isNonnull())
    dz = fabs(it->trackRef()->dz(vtx.position()));
    else if(it->gsfTrackRef().isNonnull())
    dz = fabs(it->gsfTrackRef()->dz(vtx.position()));
    else if(it->muonRef().isNonnull() && it->muonRef()->innerTrack().isNonnull())
    dz = fabs(it->muonRef()->innerTrack()->dz(vtx.position()));
    else if(it->charge()!=0) {
    cout << "WARNING: found charged PF candidate without any track ref" << endl;
    continue;
    }
  */
  /* 
     T_PFreducedCand_vz->push_back(it->vz()); 
     }
  */
    
  
  // Handle<TriggerResults> trh;
  // try {iEvent.getByLabel(trigLabel_,trh);} catch(...) {;}   
  // const edm::TriggerNames & triggerNames_=iEvent.triggerNames(*trh);
 
  
  /*  Handle<TriggerResults> trh;
      try {iEvent.getByLabel(trigLabel_,trh);
      unsigned int aux = trh.product()->size();
      aux = 0;
      } catch(...) {;}
      //    iEvent.getByLabel(trigLabel_,trh);
    
      //    ;}//Chapuza provisional para que no casquen los datos reales
      */  
  
  
  //Events
  EventID eventID = iEvent.id();
  T_Event_EventNumber = eventID.event();
  T_Event_RunNumber = eventID.run();
  T_Event_LuminosityBlock = iEvent.luminosityBlock();  

  T_Event_nPU =-1;
  T_Event_nPUp=-1;
  T_Event_nPUm=-1;
  T_Event_AveNTruePU=-1.;
  float truePu=0.;
  if(!IsRealData){
    Handle<std::vector< PileupSummaryInfo > > puInfo;
    try {
      iEvent.getByLabel("addPileupInfo",puInfo);
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      //The in-time crossing is getBunchCrossing = 0; negative ones are early, positive ones are late.
      for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {

	//    std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
	if(PVI->getBunchCrossing()==0){
	  T_Event_nPU =PVI->getPU_NumInteractions();
	  T_Event_nTruePU=PVI->getTrueNumInteractions();

	}

	else if(PVI->getBunchCrossing()==-1){
	  T_Event_nPUm=PVI->getPU_NumInteractions();
	}
	else if(PVI->getBunchCrossing()==1){
	  T_Event_nPUp=PVI->getPU_NumInteractions();
	}
	truePu += PVI->getTrueNumInteractions();
      }
    } catch (...) {}
	
  }
  T_Event_AveNTruePU=truePu/3.;

  /*
    T_Event_KFactorHNLO = 1;
    Handle<double> weightHandle;
    try {
    iEvent.getByLabel("KFactorProducer",weightHandle);
    T_Event_KFactorHNLO = T_Event_KFactorHNLO*(*weightHandle);
    } catch (...) {}
  
  */  

     
  //Rho 
  edm::Handle<double> rhoH;
  iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"),rhoH);
  T_Event_Rho=*rhoH; 
  
  iEvent.getByLabel(edm::InputTag("kt6PFJetsForIso","rho"),rhoH);
  T_Event_RhoIso=*rhoH; 
  /*
    iEvent.getByLabel(edm::InputTag("kt6PFJetsNoPU","rho"),rhoH);
    T_Event_RhoNoPu=*rhoH; 
  
    iEvent.getByLabel(edm::InputTag("kt6PFJetsForIsoNoPU","rho"),rhoH);
    T_Event_RhoIsoNoPu=*rhoH; 
  */
  //  iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralNeutral","rho"),rhoH);
  //  T_Event_RhoCentralNeutral=*rhoH;
  /*
    iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralNeutralTight","rho"),rhoH);
    T_Event_RhoCentralNeutralTight=*rhoH;
  */

  //MEt filters result

  edm::Handle<bool> result;
  
  try {iEvent.getByLabel(edm::InputTag("EcalDeadCellTriggerPrimitiveFilter",""),result);
    T_EventF_EcalDeadCell=*result;
    iEvent.getByLabel(edm::InputTag("logErrorTooManyClusters",""),result);
    T_EventF_logErrorTooManyClusters=*result;
    iEvent.getByLabel(edm::InputTag("trackingFailureFilter",""),result);
    T_EventF_trackingFailure=*result;
    iEvent.getByLabel(edm::InputTag("hcalLaserEventFilter",""),result);
    T_EventF_hcalLaser=*result;
    iEvent.getByLabel(edm::InputTag("ecalLaserCorrFilter",""),result);
    T_EventF_ecalLaserCorr=*result;
    iEvent.getByLabel(edm::InputTag("toomanystripclus53X",""),result);
    T_EventF_toomanystripclus=*result;
    iEvent.getByLabel(edm::InputTag("manystripclus53X",""),result);
    T_EventF_manystripclus=*result;
    iEvent.getByLabel(edm::InputTag("eeBadScFilter",""),result);
    T_EventF_eeBadSc=*result;
  } catch(...) {;}
 
 
  //HLT

  Handle<TriggerResults> trh;
  try {iEvent.getByLabel(trigLabel_,trh);
    unsigned int aux = trh.product()->size();
    aux = 0 + aux;
  } catch(...) {;}
  const edm::TriggerNames &triggerNames_=iEvent.triggerNames(*trh);
  /*
    T_HLT_Mu8_v1=false;
    T_HLT_Mu8_v2=false;
    T_HLT_Mu8_v3=false;
    T_HLT_Mu8_v4=false;
    T_HLT_Mu8_v5=false;
    T_HLT_Mu8_v6=false;
    T_HLT_Mu8_v7=false;
    T_HLT_Mu8_v8=false;
    T_HLT_Mu8_v9=false;
    T_HLT_Mu8_v10=false;
    T_HLT_Mu8_v11=false;
    T_HLT_Mu8_v12=false;
  */
  T_HLT_Mu8_v16=false;
  /*
    T_HLT_Mu12_v1=false;
    T_HLT_Mu12_v2=false;
    T_HLT_Mu12_v3=false;
    T_HLT_Mu12_v4=false;
    T_HLT_Mu12_v16=false;

    T_HLT_Mu15_v1=false;
    T_HLT_Mu15_v2=false;
    T_HLT_Mu15_v3=false;
    T_HLT_Mu15_v4=false;
    T_HLT_Mu15_v5=false;
    T_HLT_Mu15_v6=false;
    T_HLT_Mu15_v7=false;
    T_HLT_Mu15_v8=false;
    T_HLT_Mu15_v9=false;
    T_HLT_Mu15_v10=false;
    T_HLT_Mu15_v11=false;
    T_HLT_Mu15_v12=false;
    T_HLT_Mu15_v13=false;

    T_HLT_IsoMu24_eta2p1_v11=false;
    T_HLT_IsoMu24_eta2p1_v12=false;
  */
  T_HLT_Mu17_v3=false;



  T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12=false;
  T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13=false;
  T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14=false;
  						   
  T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3=false;
  T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4=false;
  T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5=false;
  


  /*
    T_HLT_Mu24_v1=false;
    T_HLT_Mu24_v2=false;
    T_HLT_Mu24_v3=false;
    T_HLT_Mu24_v4=false;
    T_HLT_Mu24_v5=false;
    T_HLT_Mu24_v6=false;
    T_HLT_Mu24_v7=false;
    T_HLT_Mu24_v8=false;
    T_HLT_Mu24_v9=false;
    T_HLT_Mu24_v10=false;
    T_HLT_Mu24_v11=false;
    T_HLT_Mu24_v12=false;

    T_HLT_Mu9=false;
    T_HLT_Jet30_v1=false;
    T_HLT_Jet30_v2=false;
    T_HLT_Jet30_v3=false;
    T_HLT_Jet30_v4=false;
    T_HLT_Jet30_v5=false;
    T_HLT_Jet60_v1=false;
    T_HLT_Jet60_v2=false;
    T_HLT_Jet60_v3=false;
    T_HLT_Jet60_v4=false;
    T_HLT_Jet60_v5=false;
  */
  for (unsigned int itrig= 0; itrig < trh.product()->size(); ++itrig) { 
    /*
      if (triggerNames_.triggerName(itrig) == "HLT_Mu9") {
      if (trh.product()->accept(itrig)) T_HLT_Mu9 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v1") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v1 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v2") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v2 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v3 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v4") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v4 = true;
      }
      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v5") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v5 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v6") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v6 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v7") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v7 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v8") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v8 = true;
      }
      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v9") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v9 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v10") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v10 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v11") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v11 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v12") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v12 = true;
      }
    */
    if (triggerNames_.triggerName(itrig) == "HLT_Mu8_v16") {
      if (trh.product()->accept(itrig)) T_HLT_Mu8_v16 = true;
    }
    /*
      if (triggerNames_.triggerName(itrig) == "HLT_Mu12_v1") {
      if (trh.product()->accept(itrig)) T_HLT_Mu12_v1 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu12_v2") {
      if (trh.product()->accept(itrig)) T_HLT_Mu12_v2 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu12_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Mu12_v3 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu12_v4") {
      if (trh.product()->accept(itrig)) T_HLT_Mu12_v4 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu12_v16") {
      if (trh.product()->accept(itrig)) T_HLT_Mu12_v16 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v1") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v1 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v2") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v2 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v3 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v4") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v4 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v5") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v5 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v6") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v6 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v7") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v7 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v8") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v8 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v9") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v9 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v10") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v10 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v11") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v11 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v12") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v12 = true;
      }
      if (triggerNames_.triggerName(itrig) == "HLT_Mu15_v13") {
      if (trh.product()->accept(itrig)) T_HLT_Mu15_v13 = true;
      }
    */
    if (triggerNames_.triggerName(itrig) == "HLT_Mu17_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Mu17_v3 = true;
    }



    if (triggerNames_.triggerName(itrig) == "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12") {
      if (trh.product()->accept(itrig)) T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12 = true;
    }

    if (triggerNames_.triggerName(itrig) == "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13") {
      if (trh.product()->accept(itrig)) T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13 = true;
    }

    if (triggerNames_.triggerName(itrig) == "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14") {
      if (trh.product()->accept(itrig)) T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14 = true;
    }


    if (triggerNames_.triggerName(itrig) == "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3 = true;
    }

    if (triggerNames_.triggerName(itrig) == "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4") {
      if (trh.product()->accept(itrig)) T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4 = true;
    }

    if (triggerNames_.triggerName(itrig) == "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5") {
      if (trh.product()->accept(itrig)) T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5 = true;
    }


    /*
      if (triggerNames_.triggerName(itrig) == "HLT_IsoMu24_eta2p1_v11") {
      if (trh.product()->accept(itrig)) T_HLT_IsoMu24_eta2p1_v11 = true;
      }
      if (triggerNames_.triggerName(itrig) == "HLT_IsoMu24_eta2p1_v12") {
      if (trh.product()->accept(itrig)) T_HLT_IsoMu24_eta2p1_v12 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v1") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v1 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v2") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v2 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v3 = true;
      }
      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v4") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v4 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v5") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v5 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v6") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v6 = true;
      }
      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v7") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v7 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v8") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v8 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v9") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v9 = true;
      }
      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v10") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v10 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v11") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v11 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Mu24_v12") {
      if (trh.product()->accept(itrig)) T_HLT_Mu24_v12 = true;
      }


      if (triggerNames_.triggerName(itrig) == "HLT_Jet30_v1") {
      if (trh.product()->accept(itrig)) T_HLT_Jet30_v1 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet30_v2") {
      if (trh.product()->accept(itrig)) T_HLT_Jet30_v2 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet30_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Jet30_v3 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet30_v4") {
      if (trh.product()->accept(itrig)) T_HLT_Jet30_v4 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet30_v5") {
      if (trh.product()->accept(itrig)) T_HLT_Jet30_v5 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet60_v1") {
      if (trh.product()->accept(itrig)) T_HLT_Jet60_v1 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet60_v2") {
      if (trh.product()->accept(itrig)) T_HLT_Jet60_v2 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet60_v3") {
      if (trh.product()->accept(itrig)) T_HLT_Jet60_v3 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet60_v4") {
      if (trh.product()->accept(itrig)) T_HLT_Jet60_v4 = true;
      }

      if (triggerNames_.triggerName(itrig) == "HLT_Jet60_v5") {
      if (trh.product()->accept(itrig)) T_HLT_Jet60_v5 = true;
      }*/
  } // loop in itrig

  //  if(!IsRealData){
  //Gen
  T_Gen_StopMass = new std::vector<float>;
  T_Gen_Chi0Mass = new std::vector<float>;
  T_Gen_CharginoMass = new std::vector<float>;

  T_Gen_polWeights = new std::vector<double>;

  T_Gen_Muon_PID = new std::vector<int>;
  T_Gen_Muon_Px = new std::vector<float>;
  T_Gen_Muon_Py = new std::vector<float>;
  T_Gen_Muon_Pz = new std::vector<float>;
  T_Gen_Muon_Energy = new std::vector<float>;
  
  T_Gen_Muon_MPID = new std::vector<int>;
  T_Gen_Muon_MPx = new std::vector<float>;
  T_Gen_Muon_MPy = new std::vector<float>;
  T_Gen_Muon_MPz = new std::vector<float>;
  T_Gen_Muon_MEnergy = new std::vector<float>;
  T_Gen_Muon_MSt = new std::vector<int>;

  T_Gen_Elec_PID = new std::vector<int>;
  T_Gen_Elec_Px = new std::vector<float>;
  T_Gen_Elec_Py = new std::vector<float>;
  T_Gen_Elec_Pz = new std::vector<float>;
  T_Gen_Elec_Energy = new std::vector<float>;

  T_Gen_Elec_MPID = new std::vector<int>;
  T_Gen_Elec_MPx = new std::vector<float>;
  T_Gen_Elec_MPy = new std::vector<float>;
  T_Gen_Elec_MPz = new std::vector<float>;
  T_Gen_Elec_MEnergy = new std::vector<float>;
  T_Gen_Elec_MSt = new std::vector<int>;

  T_Gen_b_PID = new std::vector<int>;
  T_Gen_b_Px = new std::vector<float>;
  T_Gen_b_Py = new std::vector<float>;
  T_Gen_b_Pz = new std::vector<float>;
  T_Gen_b_Energy = new std::vector<float>;

  T_Gen_b_MPID = new std::vector<int>;
  T_Gen_b_MPx = new std::vector<float>;
  T_Gen_b_MPy = new std::vector<float>;
  T_Gen_b_MPz = new std::vector<float>;
  T_Gen_b_MEnergy = new std::vector<float>;
  T_Gen_b_MSt = new std::vector<int>;

  /*
  T_Gen_StopSt3  = new std::vector<SUSYGenParticle>;
  T_Gen_Chi0St3  = new std::vector<SUSYGenParticle>;
  T_Gen_tSt3     = new std::vector<SUSYGenParticle>;
  T_Gen_ChiPMSt3 = new std::vector<SUSYGenParticle>;
  T_Gen_bSt3     = new std::vector<SUSYGenParticle>;
  T_Gen_WSt3     = new std::vector<SUSYGenParticle>;
  T_Gen_MuonSt3  = new std::vector<SUSYGenParticle>;
  T_Gen_ElecSt3  = new std::vector<SUSYGenParticle>;
  T_Gen_TauSt3   = new std::vector<SUSYGenParticle>;
  */

  T_Gen_StopSt3_pdgId = new std::vector<int>;	   
  T_Gen_StopSt3_firstMother = new std::vector<int>;
  T_Gen_StopSt3_i = new std::vector<int>;
  T_Gen_StopSt3_energy = new std::vector<float>;   
  T_Gen_StopSt3_pt = new std::vector<float>;	   
  T_Gen_StopSt3_eta = new std::vector<float>;	   
  T_Gen_StopSt3_phi = new std::vector<float>;      
  
  T_Gen_Chi0St3_pdgId = new std::vector<int>;	   
  T_Gen_Chi0St3_firstMother = new std::vector<int>;
  T_Gen_Chi0St3_i = new std::vector<int>;
  T_Gen_Chi0St3_energy = new std::vector<float>;   
  T_Gen_Chi0St3_pt = new std::vector<float>;	   
  T_Gen_Chi0St3_eta = new std::vector<float>;	   
  T_Gen_Chi0St3_phi = new std::vector<float>;      
  
  T_Gen_tSt3_pdgId = new std::vector<int>;	   
  T_Gen_tSt3_firstMother = new std::vector<int>;
  T_Gen_tSt3_i = new std::vector<int>;
  T_Gen_tSt3_energy = new std::vector<float>;   
  T_Gen_tSt3_pt = new std::vector<float>;	   
  T_Gen_tSt3_eta = new std::vector<float>;	   
  T_Gen_tSt3_phi = new std::vector<float>;      

  T_Gen_ChiPMSt3_pdgId = new std::vector<int>;	   
  T_Gen_ChiPMSt3_firstMother = new std::vector<int>;
  T_Gen_ChiPMSt3_i = new std::vector<int>;
  T_Gen_ChiPMSt3_energy = new std::vector<float>;   
  T_Gen_ChiPMSt3_pt = new std::vector<float>;	   
  T_Gen_ChiPMSt3_eta = new std::vector<float>;	   
  T_Gen_ChiPMSt3_phi = new std::vector<float>;      

  T_Gen_bSt3_pdgId = new std::vector<int>;	   
  T_Gen_bSt3_firstMother = new std::vector<int>;
  T_Gen_bSt3_i = new std::vector<int>;
  T_Gen_bSt3_energy = new std::vector<float>;   
  T_Gen_bSt3_pt = new std::vector<float>;	   
  T_Gen_bSt3_eta = new std::vector<float>;	   
  T_Gen_bSt3_phi = new std::vector<float>;      

  T_Gen_WSt3_pdgId = new std::vector<int>;	
  T_Gen_WSt3_firstMother = new std::vector<int>;
  T_Gen_WSt3_i = new std::vector<int>;
  T_Gen_WSt3_energy = new std::vector<float>;   
  T_Gen_WSt3_pt = new std::vector<float>;	
  T_Gen_WSt3_eta = new std::vector<float>;	
  T_Gen_WSt3_phi = new std::vector<float>;      

  T_Gen_MuonSt3_pdgId = new std::vector<int>;	
  T_Gen_MuonSt3_firstMother = new std::vector<int>;
  T_Gen_MuonSt3_i = new std::vector<int>;
  T_Gen_MuonSt3_energy = new std::vector<float>;   
  T_Gen_MuonSt3_pt = new std::vector<float>;	
  T_Gen_MuonSt3_eta = new std::vector<float>;	
  T_Gen_MuonSt3_phi = new std::vector<float>;      

  T_Gen_ElecSt3_pdgId = new std::vector<int>;	
  T_Gen_ElecSt3_firstMother = new std::vector<int>;
  T_Gen_ElecSt3_i = new std::vector<int>;
  T_Gen_ElecSt3_energy = new std::vector<float>;   
  T_Gen_ElecSt3_pt = new std::vector<float>;	
  T_Gen_ElecSt3_eta = new std::vector<float>;	
  T_Gen_ElecSt3_phi = new std::vector<float>;      

  T_Gen_TauSt3_pdgId = new std::vector<int>;	
  T_Gen_TauSt3_firstMother = new std::vector<int>;
  T_Gen_TauSt3_i = new std::vector<int>;
  T_Gen_TauSt3_energy = new std::vector<float>;   
  T_Gen_TauSt3_pt = new std::vector<float>;	
  T_Gen_TauSt3_eta = new std::vector<float>;	
  T_Gen_TauSt3_phi = new std::vector<float>;      





  /*
  T_Gen_MuonSt3_PID = new std::vector<int>;
  T_Gen_MuonSt3_Px = new std::vector<float>;
  T_Gen_MuonSt3_Py = new std::vector<float>;
  T_Gen_MuonSt3_Pz = new std::vector<float>;
  T_Gen_MuonSt3_Energy = new std::vector<float>;
  T_Gen_ElecSt3_PID = new std::vector<int>;
  T_Gen_ElecSt3_Px = new std::vector<float>;
  T_Gen_ElecSt3_Py = new std::vector<float>;
  T_Gen_ElecSt3_Pz = new std::vector<float>;
  T_Gen_ElecSt3_Energy = new std::vector<float>;
  T_Gen_bSt3_PID = new std::vector<int>;
  T_Gen_bSt3_Px = new std::vector<float>;
  T_Gen_bSt3_Py = new std::vector<float>;
  T_Gen_bSt3_Pz = new std::vector<float>;
  T_Gen_bSt3_Energy = new std::vector<float>;
  T_Gen_tSt3_PID = new std::vector<int>;
  T_Gen_tSt3_Px = new std::vector<float>;
  T_Gen_tSt3_Py = new std::vector<float>;
  T_Gen_tSt3_Pz = new std::vector<float>;
  T_Gen_tSt3_Energy = new std::vector<float>;
  
  T_Gen_TauSt3_PID = new std::vector<int>;
  T_Gen_TauSt3_Px = new std::vector<float>;
  T_Gen_TauSt3_Py = new std::vector<float>;
  T_Gen_TauSt3_Pz = new std::vector<float>;
  T_Gen_TauSt3_Energy = new std::vector<float>;
  */

  T_Gen_TauSt3_IsLepDec = new std::vector<bool>;
  T_Gen_TauSt3_LepDec_PID = new std::vector<int>;
  T_Gen_TauSt3_LepDec_Px = new std::vector<float>;
  T_Gen_TauSt3_LepDec_Py = new std::vector<float>;
  T_Gen_TauSt3_LepDec_Pz = new std::vector<float>;
  T_Gen_TauSt3_LepDec_Energy = new std::vector<float>;  
  
  if(!IsRealData){
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator", genEvtInfo);
    //T_Event_PtHat =  genEvtInfo->hasBinningValues() ? (genEvtInfo->binningValues())[0] : 0.0;
    T_Event_processID= genEvtInfo->signalProcessID();

    std::vector<SUSYGenParticle> genParVec;
    genParVec.clear();

    for (size_t i = 0; i < genParticles->size(); ++i){
      const Candidate & p = (*genParticles)[i];
      int id = p.pdgId();
      int st = p.status();

      if (abs(id) == 1000006) T_Gen_StopMass->push_back(p.mass());      // enable doSusy in latinosYieldSkim.py
      if (abs(id) == 1000022) T_Gen_Chi0Mass->push_back(p.mass());
      if (abs(id) == 1000024) T_Gen_CharginoMass->push_back(p.mass());

      //if (st ==3 && (abs(id) == 11 || abs(id) == 13 || abs(id) == 5)){
      if( st == 3 && (abs(id) == 11 || abs(id) == 13 || abs(id) == 15 || abs(id)==24 || abs(id)==5 || abs(id)==6 || abs(id)==1000006 || abs(id)==1000022 || abs(id)==1000024) ){
	// 11 == e, 13 == mu, 15 == tau, 24 == W, 5 == b, 6 == t, 1000006 == stop_1, 1000022 == chi0_1, 1000024 == chi+_1
	
	SUSYGenParticle genPar;
	genPar.pdgId  = p.pdgId();
	genPar.firstMother = -1; // placeholder, set later on. It's the index of the first mother in genParVec
	genPar.energy = p.energy();
	genPar.pt     = p.pt();
	genPar.eta    = p.eta();
	genPar.phi    = p.phi();

	genParVec.push_back(genPar);      
	/*
	if (abs(id) == 11) {
	  T_Gen_ElecSt3_PID->push_back(id);
	  T_Gen_ElecSt3_Px->push_back(p.px());
	  T_Gen_ElecSt3_Py->push_back(p.py());
	  T_Gen_ElecSt3_Pz->push_back(p.pz());
	  T_Gen_ElecSt3_Energy->push_back(p.energy());
	}      
	else if (abs(id) == 13) {
	  T_Gen_MuonSt3_PID->push_back(id);
	  T_Gen_MuonSt3_Px->push_back(p.px());
	  T_Gen_MuonSt3_Py->push_back(p.py());
	  T_Gen_MuonSt3_Pz->push_back(p.pz());
	  T_Gen_MuonSt3_Energy->push_back(p.energy());
	}      
	else if (abs(id) == 5) {
	  T_Gen_bSt3_PID->push_back(id);
	  T_Gen_bSt3_Px->push_back(p.px());
	  T_Gen_bSt3_Py->push_back(p.py());
	  T_Gen_bSt3_Pz->push_back(p.pz());
	  T_Gen_bSt3_Energy->push_back(p.energy());
	}
	else if (abs(id) == 6) {
	  T_Gen_tSt3_PID->push_back(id);
	  T_Gen_tSt3_Px->push_back(p.px());
	  T_Gen_tSt3_Py->push_back(p.py());
	  T_Gen_tSt3_Pz->push_back(p.pz());
	  T_Gen_tSt3_Energy->push_back(p.energy());
	}	
	*/
      } // status ==3

      if ((st ==1 && (abs(id) == 11 || abs(id) == 13)) || (st ==2 && abs(id) == 5)){
	// get mother of gen_lept
	const GenParticle* gen_mom = static_cast<const GenParticle*> (p.mother());
	int m_id=id;
	if(gen_mom!=0) m_id = gen_mom -> pdgId();
	else m_id=0;
      
	if(m_id != id);
	else{
	  int id2= m_id;
	  while(id2 == id){
	    gen_mom = static_cast<const GenParticle*> (gen_mom->mother());
	    if(gen_mom!=0) id2=gen_mom->pdgId();
	    else id2=0;
	  }
	}
	if(gen_mom!=0) m_id = gen_mom->pdgId();
	else m_id=0;
            
      
	if (abs(id) == 11) {
	  T_Gen_Elec_PID->push_back(id);
	  T_Gen_Elec_Px->push_back(p.px());
	  T_Gen_Elec_Py->push_back(p.py());
	  T_Gen_Elec_Pz->push_back(p.pz());
	  T_Gen_Elec_Energy->push_back(p.energy());
	
	  T_Gen_Elec_MPID->push_back(m_id);
	  if(gen_mom!=0){
	    T_Gen_Elec_MPx->push_back(gen_mom->px());
	    T_Gen_Elec_MPy->push_back(gen_mom->py());
	    T_Gen_Elec_MPz->push_back(gen_mom->pz());
	    T_Gen_Elec_MEnergy->push_back(gen_mom->energy());
	    T_Gen_Elec_MSt->push_back(gen_mom->status());
	  }
	  else{
	    T_Gen_Elec_MPx->push_back(0);
	    T_Gen_Elec_MPy->push_back(0);
	    T_Gen_Elec_MPz->push_back(0);
	    T_Gen_Elec_MEnergy->push_back(0);
	    T_Gen_Elec_MSt->push_back(0);
	  }
	}
      
	else if (abs(id) == 13) {
	  T_Gen_Muon_PID->push_back(id);
	  T_Gen_Muon_Px->push_back(p.px());
	  T_Gen_Muon_Py->push_back(p.py());
	  T_Gen_Muon_Pz->push_back(p.pz());
	  T_Gen_Muon_Energy->push_back(p.energy());
	
	  T_Gen_Muon_MPID->push_back(m_id);
	  if(gen_mom!=0){
	    T_Gen_Muon_MPx->push_back(gen_mom->px());
	    T_Gen_Muon_MPy->push_back(gen_mom->py());
	    T_Gen_Muon_MPz->push_back(gen_mom->pz());
	    T_Gen_Muon_MEnergy->push_back(gen_mom->energy());
	    T_Gen_Muon_MSt->push_back(gen_mom->status());
	  }
	  else{
	    T_Gen_Muon_MPx->push_back(0);
	    T_Gen_Muon_MPy->push_back(0);
	    T_Gen_Muon_MPz->push_back(0);
	    T_Gen_Muon_MEnergy->push_back(0);
	    T_Gen_Muon_MSt->push_back(0);
	  }
	}
      
	else if (abs(id) == 5) {
	  T_Gen_b_PID->push_back(id);
	  T_Gen_b_Px->push_back(p.px());
	  T_Gen_b_Py->push_back(p.py());
	  T_Gen_b_Pz->push_back(p.pz());
	  T_Gen_b_Energy->push_back(p.energy());
	
	  T_Gen_b_MPID->push_back(m_id);
	  if(gen_mom!=0){
	    T_Gen_b_MPx->push_back(gen_mom->px());
	    T_Gen_b_MPy->push_back(gen_mom->py());
	    T_Gen_b_MPz->push_back(gen_mom->pz());
	    T_Gen_b_MEnergy->push_back(gen_mom->energy());
	    T_Gen_b_MSt->push_back(gen_mom->status());
	  }
	  else{
	    T_Gen_b_MPx->push_back(0);
	    T_Gen_b_MPy->push_back(0);
	    T_Gen_b_MPz->push_back(0);
	    T_Gen_b_MEnergy->push_back(0);
	    T_Gen_b_MSt->push_back(0);
	  }
	}      
      } //status 1
    
      if (st ==3 && abs(id) == 15) {
	/*
	T_Gen_TauSt3_PID->push_back(id);
	T_Gen_TauSt3_Px->push_back(p.px());
	T_Gen_TauSt3_Py->push_back(p.py());
	T_Gen_TauSt3_Pz->push_back(p.pz());
	T_Gen_TauSt3_Energy->push_back(p.energy());
	*/
	bool elecdec = false, muondec = false;
	int pid = 0;
	float px = 0, py = 0, pz = 0, energy = 0;
	LeptonicTauDecay(p, elecdec, muondec, pid, px, py, pz, energy);
	T_Gen_TauSt3_IsLepDec->push_back(elecdec || muondec);
	T_Gen_TauSt3_LepDec_PID->push_back(pid);
	T_Gen_TauSt3_LepDec_Px->push_back(px);
	T_Gen_TauSt3_LepDec_Py->push_back(py);
	T_Gen_TauSt3_LepDec_Pz->push_back(pz);
	T_Gen_TauSt3_LepDec_Energy->push_back(energy);  
      }
    
    } // loop over GenParticles
    

    // another loop to update firstMother for each SUSYGenParticle in genParVec
    for (size_t k = 0; k < genParticles->size(); ++k) {
      const Candidate & pp = (*genParticles)[k];
      int gpid = pp.pdgId();
      int gpst = pp.status();
      
      //make sure that the same status 3 particles are selected as when filling genParVec
      if( gpst != 3 || !( abs(gpid) == 11 || abs(gpid) == 13 || abs(gpid) == 15 || abs(gpid)==24 || abs(gpid)==5 || abs(gpid)==6 || abs(gpid)==1000006 || abs(gpid)==1000022 || abs(gpid)==1000024) ) continue;
      
      int firstMom_pdgId=-1;
      float firstMom_pt = -999., firstMom_eta = -999., firstMom_phi = -999.;
      
      if(pp.numberOfMothers() > 0){

	firstMom_pdgId = pp.mother()->pdgId();
	firstMom_pt    = pp.mother()->pt();
	firstMom_eta   = pp.mother()->eta();
	firstMom_phi   = pp.mother()->phi();
      }
      
      for(unsigned int jm=0; jm<genParVec.size(); jm++){  // looking for the first mother's index in genParVec
	
	if(firstMom_pdgId == genParVec.at(jm).pdgId &&
	   firstMom_pt    == genParVec.at(jm).pt    &&
	   firstMom_eta   == genParVec.at(jm).eta   &&
	   firstMom_phi   == genParVec.at(jm).phi 
	   ){
	  
	  for(unsigned int jd=0; jd<genParVec.size(); jd++){  // looking for the daughter's index in genParVec

	    if(pp.pdgId() == genParVec.at(jd).pdgId &&
	       pp.pt()    == genParVec.at(jd).pt    &&
	       pp.eta()   == genParVec.at(jd).eta   &&
	       pp.phi()   == genParVec.at(jd).phi 
	       ){
	    
	      genParVec.at(jd).firstMother = jm;   // updating firstMother of a SUSYGenParticle with its mother's index
	      
	    } // if found daughter
	  } // loop over possible daughters
	} // if found mother
      } // loop over possible mothers
      
    }    // 2nd loop over genParticles


    for(unsigned int kk=0; kk<genParVec.size(); kk++){ 
      /*
      if( abs(genParVec.at(kk).pdgId ) == 1000006 ) T_Gen_StopSt3  ->push_back( genParVec.at(kk) );
      if( abs(genParVec.at(kk).pdgId ) == 1000022 ) T_Gen_Chi0St3  ->push_back( genParVec.at(kk) ); 
      if( abs(genParVec.at(kk).pdgId ) == 6       ) T_Gen_tSt3     ->push_back( genParVec.at(kk) );    
      if( abs(genParVec.at(kk).pdgId ) == 1000024 ) T_Gen_ChiPMSt3 ->push_back( genParVec.at(kk) );
      if( abs(genParVec.at(kk).pdgId ) == 5       ) T_Gen_bSt3     ->push_back( genParVec.at(kk) );
      if( abs(genParVec.at(kk).pdgId ) == 24      ) T_Gen_WSt3     ->push_back( genParVec.at(kk) );
      if( abs(genParVec.at(kk).pdgId ) == 13      ) T_Gen_MuonSt3  ->push_back( genParVec.at(kk) );
      if( abs(genParVec.at(kk).pdgId ) == 11      ) T_Gen_ElecSt3  ->push_back( genParVec.at(kk) );
      if( abs(genParVec.at(kk).pdgId ) == 15      ) T_Gen_TauSt3   ->push_back( genParVec.at(kk) );
      */      


      if( abs(genParVec.at(kk).pdgId ) == 1000006 ){				   
									   
	T_Gen_StopSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );	   
	T_Gen_StopSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_StopSt3_i           ->push_back( kk );
	T_Gen_StopSt3_energy      ->push_back( genParVec.at(kk).energy );	   
	T_Gen_StopSt3_pt          ->push_back( genParVec.at(kk).pt );	   
	T_Gen_StopSt3_eta         ->push_back( genParVec.at(kk).eta );	   
	T_Gen_StopSt3_phi         ->push_back( genParVec.at(kk).phi );	   
      }                                                                    
      
      if( abs(genParVec.at(kk).pdgId ) == 1000022 ){			      	   
									      
	T_Gen_Chi0St3_pdgId       ->push_back( genParVec.at(kk).pdgId );      	   
	T_Gen_Chi0St3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_Chi0St3_i           ->push_back( kk );
	T_Gen_Chi0St3_energy      ->push_back( genParVec.at(kk).energy );     	   
	T_Gen_Chi0St3_pt          ->push_back( genParVec.at(kk).pt );	      
	T_Gen_Chi0St3_eta         ->push_back( genParVec.at(kk).eta );	      
	T_Gen_Chi0St3_phi         ->push_back( genParVec.at(kk).phi );	      
      }                                                                       
      
      if( abs(genParVec.at(kk).pdgId ) == 6 ){				   
									   
	T_Gen_tSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );	   
	T_Gen_tSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_tSt3_i           ->push_back( kk );
	T_Gen_tSt3_energy      ->push_back( genParVec.at(kk).energy );	   
	T_Gen_tSt3_pt          ->push_back( genParVec.at(kk).pt );	   
	T_Gen_tSt3_eta         ->push_back( genParVec.at(kk).eta );	   
	T_Gen_tSt3_phi         ->push_back( genParVec.at(kk).phi );	   
      }                                                                    

      if( abs(genParVec.at(kk).pdgId ) == 1000024 ){			      
									      
	T_Gen_ChiPMSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );      
	T_Gen_ChiPMSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_ChiPMSt3_i           ->push_back( kk );
	T_Gen_ChiPMSt3_energy      ->push_back( genParVec.at(kk).energy );     
	T_Gen_ChiPMSt3_pt          ->push_back( genParVec.at(kk).pt );	      
	T_Gen_ChiPMSt3_eta         ->push_back( genParVec.at(kk).eta );	      
	T_Gen_ChiPMSt3_phi         ->push_back( genParVec.at(kk).phi );	      
      }                                                                       

      if( abs(genParVec.at(kk).pdgId ) == 5 ){				   
									   
	T_Gen_bSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );	   
	T_Gen_bSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_bSt3_i           ->push_back( kk );
	T_Gen_bSt3_energy      ->push_back( genParVec.at(kk).energy );	   
	T_Gen_bSt3_pt          ->push_back( genParVec.at(kk).pt );	   
	T_Gen_bSt3_eta         ->push_back( genParVec.at(kk).eta );	   
	T_Gen_bSt3_phi         ->push_back( genParVec.at(kk).phi );	   
      }                                                                    

      if( abs(genParVec.at(kk).pdgId ) == 24 ){				   
									   
	T_Gen_WSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );	   
	T_Gen_WSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_WSt3_i           ->push_back( kk );
	T_Gen_WSt3_energy      ->push_back( genParVec.at(kk).energy );	   
	T_Gen_WSt3_pt          ->push_back( genParVec.at(kk).pt );	   
	T_Gen_WSt3_eta         ->push_back( genParVec.at(kk).eta );	   
	T_Gen_WSt3_phi         ->push_back( genParVec.at(kk).phi );	   
      }                                                                    

      if( abs(genParVec.at(kk).pdgId ) == 13 ){				   
									   
	T_Gen_MuonSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );	   
	T_Gen_MuonSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_MuonSt3_i           ->push_back( kk );
	T_Gen_MuonSt3_energy      ->push_back( genParVec.at(kk).energy );	   
	T_Gen_MuonSt3_pt          ->push_back( genParVec.at(kk).pt );	   
	T_Gen_MuonSt3_eta         ->push_back( genParVec.at(kk).eta );	   
	T_Gen_MuonSt3_phi         ->push_back( genParVec.at(kk).phi );	   
      }                                                                    

      if( abs(genParVec.at(kk).pdgId ) == 11 ){				   
									   
	T_Gen_ElecSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );	   
	T_Gen_ElecSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_ElecSt3_i           ->push_back( kk );
	T_Gen_ElecSt3_energy      ->push_back( genParVec.at(kk).energy );	   
	T_Gen_ElecSt3_pt          ->push_back( genParVec.at(kk).pt );	   
	T_Gen_ElecSt3_eta         ->push_back( genParVec.at(kk).eta );	   
	T_Gen_ElecSt3_phi         ->push_back( genParVec.at(kk).phi );	   
      }                                                                    

      if( abs(genParVec.at(kk).pdgId ) == 15 ){				   
									   
	T_Gen_TauSt3_pdgId       ->push_back( genParVec.at(kk).pdgId );	   
	T_Gen_TauSt3_firstMother ->push_back( genParVec.at(kk).firstMother );
	T_Gen_TauSt3_i           ->push_back( kk );
	T_Gen_TauSt3_energy      ->push_back( genParVec.at(kk).energy );	   
	T_Gen_TauSt3_pt          ->push_back( genParVec.at(kk).pt );	   
	T_Gen_TauSt3_eta         ->push_back( genParVec.at(kk).eta );	   
	T_Gen_TauSt3_phi         ->push_back( genParVec.at(kk).phi );	   
      }                                                                    

    }
    
    
    double deltaSinThetaMix = 0.1, sinThetaMix = 0.;
    while(sinThetaMix <= 1.){
      
      double polWeight = Reweight_Stop_to_TopChi0_with_SUSYmodel( genParVec, TMath::ASin(sinThetaMix) );

      T_Gen_polWeights->push_back(polWeight);
      
      sinThetaMix += deltaSinThetaMix;
    }

  }// !IsRealData
    
  //Vertex
  T_Vertex_z = new std::vector<float>;   
  T_Vertex_y = new std::vector<float>;
  T_Vertex_x = new std::vector<float>; 
  T_Vertex_Chi2Prob= new std::vector<float>;
  T_Vertex_rho = new std::vector<float>;
  T_Vertex_ndof = new std::vector<float>;
  T_Vertex_isFake = new std::vector<bool>;
  T_Vertex_tracksSize = new std::vector<int>; 
  
  //Vertex
  
  
  if (vtxs.size() != 0){
    for (size_t i=0; i < vtxs.size(); i++){
      T_Vertex_z->push_back(vtxs[i].z());
      T_Vertex_y->push_back(vtxs[i].y());
      T_Vertex_x->push_back(vtxs[i].x());
      T_Vertex_Chi2Prob->push_back(ChiSquaredProbability(vtxs[i].chi2(),vtxs[i].ndof()));
      T_Vertex_rho->push_back( vtxs[i].position().Rho());
      T_Vertex_ndof->push_back(vtxs[i].ndof());
      T_Vertex_isFake->push_back(vtxs[i].isFake());
      T_Vertex_tracksSize->push_back(vtxs[i].tracksSize());      
    }
  } 
  
  
  
  //Muons        
  T_Muon_Eta = new std::vector<float>;
  T_Muon_IsGlobalMuon = new std::vector<bool>;
  //  T_Muon_IsTMLSOLPT = new std::vector<bool>;
  //  T_Muon_IsTMLSOLPL = new std::vector<bool>;
  T_Muon_IsGMPTMuons = new std::vector<bool>;
  
  //  T_Muon_IsAllStandAloneMuons = new std::vector<bool>;
  T_Muon_IsAllTrackerMuons = new std::vector<bool>;
  T_Muon_IsTrackerMuonArbitrated = new std::vector<bool>;
  T_Muon_IsAllArbitrated = new std::vector<bool>;
  /*  T_Muon_IsTMLastStationLoose = new std::vector<bool>;
      T_Muon_IsTMLastStationTight = new std::vector<bool>;
      T_Muon_IsTM2DCompatibilityLoose = new std::vector<bool>;
      T_Muon_IsTM2DCompatibilityTight = new std::vector<bool>;
      T_Muon_IsTMOneStationLoose = new std::vector<bool>;
      T_Muon_IsTMOneStationTight = new std::vector<bool>;
      T_Muon_IsTMLSOPL = new std::vector<bool>;
      T_Muon_IsGMTkChiCompatibility = new std::vector<bool>;
      T_Muon_IsGMStaChiCompatibility = new std::vector<bool>;
      T_Muon_IsGMTkKinkTight = new std::vector<bool>;
      T_Muon_IsTMLastStationAngLoose= new std::vector<bool>;
      T_Muon_IsTMLastStationAngTight= new std::vector<bool>;
      T_Muon_IsTMOneStationAngLoose= new std::vector<bool>;
      T_Muon_IsTMOneStationAngTight= new std::vector<bool>;
  */
  T_Muon_SegmentCompatibility= new std::vector<float>;
  T_Muon_trkKink  = new std::vector<float>;
  T_Muon_Px = new std::vector<float>;
  T_Muon_Py = new std::vector<float>;
  T_Muon_Pz = new std::vector<float>;
  T_Muon_Pt = new std::vector<float>;
  T_Muon_deltaPt = new std::vector<float>;
  T_Muon_Energy = new std::vector<float>;
  T_Muon_Charge = new std::vector<int>;
  T_Muon_NormChi2GTrk = new std::vector<float>;
  T_Muon_NValidHitsInTrk = new std::vector<int>;
  T_Muon_NValidPixelHitsInTrk = new std::vector<int>;
  T_Muon_NValidHitsSATrk = new std::vector<int>;
  T_Muon_NValidHitsGTrk = new std::vector<int>;
  T_Muon_NumOfMatchedStations = new std::vector<int>;
  T_Muon_Chi2InTrk = new std::vector<float>;
  T_Muon_dofInTrk = new std::vector<float>;
  //T_Muon_SumIsoCalo = new std::vector<float>;
  //T_Muon_SumIsoTrack = new std::vector<float>;
  T_Muon_IPAbsGTrack = new std::vector<float>;
  //  T_Muon_IPSigGTrack = new std::vector<float>;
  T_Muon_IPAbsInTrack = new std::vector<float>;
  //  T_Muon_IPSigInTrack = new std::vector<float>;
  //  T_Muon_IPwrtBSInTrack =  new std::vector<float>;
  T_Muon_IPwrtAveBSInTrack =  new std::vector<float>;
  //  T_Muon_IPwrtBSGTrack =  new std::vector<float>;
  /*  T_Muon_IP2DBiasedPV =  new std::vector<float>;
      T_Muon_IP3DBiasedPV =  new std::vector<float>;
      T_Muon_IP2DUnBiasedPV =  new std::vector<float>;
      T_Muon_IP3DUnBiasedPV =  new std::vector<float>;
      T_Muon_dxyPVBiasedPV =  new std::vector<float>;
      T_Muon_dzPVBiasedPV =  new std::vector<float>;
      T_Muon_dxyPVUnBiasedPV =  new std::vector<float>;
      T_Muon_dzPVUnBiasedPV =  new std::vector<float>;
      T_Muon_IP2DUnBiasedPVnoBS =  new std::vector<float>;
      T_Muon_IP3DUnBiasedPVnoBS =  new std::vector<float>;
      T_Muon_dxyPVUnBiasedPVnoBS =  new std::vector<float>;
      T_Muon_dzPVUnBiasedPVnoBS =  new std::vector<float>;*/
  T_Muon_InnerTrackFound=new std::vector<int>;
  /*T_Muon_smurfCharged= new std::vector<float>;
    T_Muon_smurfNeutral= new std::vector<float>;
    T_Muon_smurfPhoton= new std::vector<float>;
    T_Muon_smurfNoOverCharged= new std::vector<float>;
    T_Muon_smurfNoOverNeutral= new std::vector<float>;
    T_Muon_smurfNoOverPhoton= new std::vector<float>;
    T_Muon_muSmurfPF= new std::vector<float>;*/
  //  T_Muon_chargedParticleIsoR04 = new std::vector<float>;
  T_Muon_chargedHadronIsoR04 = new std::vector<float>;
  T_Muon_neutralHadronIsoR04 = new std::vector<float>;
  T_Muon_photonIsoR04 = new std::vector<float>;
  T_Muon_sumPUPtR04 = new std::vector<float>;
  T_Muon_chargedParticleIsoR03 = new std::vector<float>;
  T_Muon_chargedHadronIsoR03 = new std::vector<float>;
  T_Muon_neutralHadronIsoR03 = new std::vector<float>;
  T_Muon_photonIsoR03 = new std::vector<float>;
  T_Muon_sumPUPtR03 = new std::vector<float>;
  T_Muon_vz = new std::vector<float>;
  T_Muon_vy = new std::vector<float>;  
  T_Muon_vx = new std::vector<float>;
  T_Muon_PFMuonPt = new std::vector<float>;
  T_Muon_PFMuonPx = new std::vector<float>;
  T_Muon_PFMuonPy = new std::vector<float>;
  T_Muon_PFMuonPz = new std::vector<float>;
  T_Muon_PFMuonE = new std::vector<float>;
  T_Muon_isPFMuon = new std::vector<bool>;
  T_Muon_NLayers =  new std::vector<int>;

  /*  T_Muon_passTriggerSingleMu = new std::vector<bool>;
      T_Muon_passTriggerDoubleMu = new std::vector<bool>;
      T_Muon_passTriggerElMu = new std::vector<bool>;
  */ 
  
  //Muons
  //no estoy seguro de que haga falta ordenar (se puede usar para filtrar en el futuro)

  std::map<float,pat::Muon> muonMap;
  for (size_t i = 0; i< muonHandle->size(); ++i) {
    muonMap[(*muonHandle)[i].pt()]=(*muonHandle)[i];  
  }
  std::vector<pat::Muon> selected_muons;
  for( std::map<float,pat::Muon>::reverse_iterator rit=muonMap.rbegin(); rit!=muonMap.rend(); ++rit){
    selected_muons.push_back( (*rit).second );
  }

  
  for (size_t k = 0; k < selected_muons.size(); ++k) {
    /*   leps_.push_back(selected_muons[k]);
	 T_Muon_passTriggerSingleMu->push_back(passTriggerSingleMu(k,IsRealData));
	 T_Muon_passTriggerDoubleMu->push_back(passTriggerDoubleMu(k,IsRealData));
	 T_Muon_passTriggerElMu->push_back(passTriggerElMu(k,IsRealData));
    */
    float IP      = 9999;
    //float IPError = -9999;
    //float IPSignificance = 9999;
    float normchi2 = 9999;
    //quantities wrt global track    
    reco::TrackRef tr_globaltrack = selected_muons[k].globalTrack();       
    if (!tr_globaltrack.isNull() && selected_muons[k].isGlobalMuon()) {
    
      normchi2 = selected_muons[k].globalTrack()->normalizedChi2();

      if (vtxs.size() > 0) {

	//       try{IPTools::ImpactParameterComputer IPComp(vtxs[0]); // use the vtx with the highest pt sum in this example 

	//	Measurement1D mess1D = IPComp.computeIP(iSetup, *tr_globaltrack);
	// to calculate the ImpactParameter, its Error and the Significance use
	IP             =  selected_muons[k].globalTrack()->dxy(vtxs[0].position());

	//	IPSignificance =  mess1D.significance();
	//} catch(...) {
	//	IP=999.;
	//	IPSignificance=-999999.;}
      }
      
    
    }


    //cuantities wrt inner track
    reco::TrackRef tr_innertrack = selected_muons[k].innerTrack(); 
    int nhitsinnertracker = -1,pixelHits=-1, found=-1, nLayers=-1;
    float chi2innertracker=9999;
    float dofinnertracker=9999;
    float IPIn      = 9999;
    //    float IPSignificanceIn = 9999;
    //    float IPwrtBSIn=9999;
    float deltaPt =9999.;
    if (!tr_innertrack.isNull()) {
      if (vtxs.size() > 0) {

	
        //IPTools::ImpactParameterComputer IPComp(vtxs[0]); // use the vtx with the highest pt sum in this example 
	//Measurement1D mess1D = IPComp.computeIP(iSetup, *tr_innertrack);
	

	// to calculate the ImpactParameter, its Error and the Significance use
	IPIn             =   selected_muons[k].innerTrack()->dxy(vtxs[0].position());
	//	IPSignificanceIn =  mess1D.significance();
      }
      nhitsinnertracker = tr_innertrack->hitPattern().numberOfValidTrackerHits();
      pixelHits = tr_innertrack->hitPattern().numberOfValidPixelHits(); 
      nLayers = tr_innertrack->hitPattern().trackerLayersWithMeasurement();
      chi2innertracker=tr_innertrack->chi2();
      dofinnertracker=tr_innertrack->ndof();
      //      IPwrtBSIn=tr_innertrack->dxy(beamSpot.position());
      deltaPt=tr_innertrack->ptError();
      found=tr_innertrack->found();
    }

    
    reco::TrackRef tr_outtrack = selected_muons[k].standAloneMuon(); 

    float nhitsouttrack=9999;
    
    if (!tr_outtrack.isNull()) {
    
      nhitsouttrack=selected_muons[k].standAloneMuon()->hitPattern().numberOfValidMuonHits();
    
    }

    int numOfValidHitsGTrk=0;
    //    float IPwrtBSG = 9999;

    if (selected_muons[k].isGlobalMuon()){

      numOfValidHitsGTrk=selected_muons[k].globalTrack()-> hitPattern().numberOfValidMuonHits();
      //   IPwrtBSG=selected_muons[k].globalTrack()->dxy(beamSpot.position());
    }

    T_Muon_Eta->push_back(selected_muons[k].eta()); 
    
    T_Muon_IsGlobalMuon->push_back(selected_muons[k].isGlobalMuon());
    //    T_Muon_IsTMLSOLPT->push_back(selected_muons[k].muonID("TMLastStationOptimizedLowPtTight"));
    //    T_Muon_IsTMLSOLPL->push_back(selected_muons[k].muonID("TMLastStationOptimizedLowPtLoose"));
    T_Muon_IsGMPTMuons->push_back(selected_muons[k].muonID("GlobalMuonPromptTight"));
    //    T_Muon_IsAllStandAloneMuons->push_back(selected_muons[k].muonID("AllStandAloneMuons"));
    T_Muon_IsAllTrackerMuons->push_back(selected_muons[k].muonID("AllTrackerMuons"));
    T_Muon_IsTrackerMuonArbitrated->push_back(selected_muons[k].muonID("TrackerMuonArbitrated"));
    T_Muon_IsAllArbitrated->push_back(selected_muons[k].muonID("AllArbitrated"));
    /*    T_Muon_IsTMLastStationLoose->push_back(selected_muons[k].muonID("TMLastStationLoose"));
	  T_Muon_IsTMLastStationTight->push_back(selected_muons[k].muonID("TMLastStationTight"));
	  T_Muon_IsTM2DCompatibilityLoose->push_back(selected_muons[k].muonID("TM2DCompatibilityLoose"));
	  T_Muon_IsTM2DCompatibilityTight->push_back(selected_muons[k].muonID("TM2DCompatibilityTight"));
	  T_Muon_IsTMOneStationLoose->push_back(selected_muons[k].muonID("TMOneStationLoose"));
	  T_Muon_IsTMOneStationTight->push_back(selected_muons[k].muonID("TMOneStationTight"));
	  T_Muon_IsTMLSOPL->push_back(selected_muons[k].muonID("TMLastStationOptimizedLowPtLoose"));
	  T_Muon_IsGMTkChiCompatibility->push_back(selected_muons[k].muonID("GMTkChiCompatibility"));
	  T_Muon_IsGMStaChiCompatibility->push_back(selected_muons[k].muonID("GMStaChiCompatibility"));
	  T_Muon_IsGMTkKinkTight->push_back(selected_muons[k].muonID("GMTkKinkTight"));
	  T_Muon_IsTMLastStationAngLoose->push_back(selected_muons[k].muonID("TMLastStationAngLoose"));
	  T_Muon_IsTMLastStationAngTight->push_back(selected_muons[k].muonID("TMLastStationAngTight"));
	  T_Muon_IsTMOneStationAngLoose->push_back(selected_muons[k].muonID("TMOneStationAngLoose"));
	  T_Muon_IsTMOneStationAngTight->push_back(selected_muons[k].muonID("TMOneStationAngTight"));
    */
    T_Muon_Px->push_back(selected_muons[k].px());
    T_Muon_Py->push_back(selected_muons[k].py());
    T_Muon_Pz->push_back(selected_muons[k].pz());
    T_Muon_Pt->push_back(selected_muons[k].pt());
    T_Muon_InnerTrackFound->push_back(found);
    T_Muon_deltaPt->push_back(deltaPt);
    T_Muon_Energy->push_back(selected_muons[k].energy());
    T_Muon_Charge->push_back(selected_muons[k].charge());
    T_Muon_NormChi2GTrk->push_back(normchi2);
    T_Muon_NValidHitsInTrk->push_back(nhitsinnertracker);
    T_Muon_NValidPixelHitsInTrk->push_back(pixelHits);
    T_Muon_Chi2InTrk->push_back(chi2innertracker);
    T_Muon_dofInTrk->push_back(dofinnertracker);
    //T_Muon_SumIsoCalo->push_back(selected_muons[k].caloIso());
    //T_Muon_SumIsoTrack->push_back(selected_muons[k].trackIso());
    T_Muon_IPAbsGTrack->push_back(fabs(IP));
    //    T_Muon_IPSigGTrack->push_back(IPSignificance);
    T_Muon_IPAbsInTrack->push_back(fabs(IPIn));
    //    T_Muon_IPSigInTrack->push_back(IPSignificanceIn);
    //    T_Muon_IPwrtBSInTrack->push_back(IPwrtBSIn);
    T_Muon_IPwrtAveBSInTrack->push_back(selected_muons[k].dB());
    //    T_Muon_IPwrtBSGTrack->push_back(IPwrtBSG);
    /*    T_Muon_IP2DBiasedPV->push_back(selected_muons[k].userFloat("tip"));
	  T_Muon_IP3DBiasedPV->push_back(selected_muons[k].userFloat("ip"));
	  T_Muon_IP2DUnBiasedPV->push_back(selected_muons[k].userFloat("tip2"));
	  T_Muon_IP3DUnBiasedPV->push_back(selected_muons[k].userFloat("ip2"));
	  T_Muon_dxyPVBiasedPV->push_back(selected_muons[k].userFloat("dxyPV"));
	  T_Muon_dxyPVUnBiasedPV->push_back(selected_muons[k].userFloat("dxyPV2"));
	  T_Muon_dzPVBiasedPV->push_back(selected_muons[k].userFloat("dzPV"));
	  T_Muon_dzPVUnBiasedPV->push_back(selected_muons[k].userFloat("dzPV2"));
	  T_Muon_IP2DUnBiasedPVnoBS->push_back(selected_muons[k].userFloat("tip3"));
	  T_Muon_IP3DUnBiasedPVnoBS->push_back(selected_muons[k].userFloat("ip3"));
	  T_Muon_dxyPVUnBiasedPVnoBS->push_back(selected_muons[k].userFloat("dxyPV3"));
	  T_Muon_dzPVUnBiasedPVnoBS->push_back(selected_muons[k].userFloat("dzPV3"));*/
    /*    T_Muon_pfCharged->push_back(selected_muons[k].userFloat("pfCharged"));
	  T_Muon_pfNeutral->push_back(selected_muons[k].userFloat("pfNeutral"));
	  T_Muon_pfPhoton->push_back(selected_muons[k].userFloat("pfPhoton"));
    */
    //    T_Muon_chargedParticleIsoR04->push_back(selected_muons[k].pfIsolationR04().sumChargedParticlePt);
    T_Muon_chargedHadronIsoR04->push_back(selected_muons[k].pfIsolationR04().sumChargedHadronPt);
    T_Muon_neutralHadronIsoR04->push_back(selected_muons[k].pfIsolationR04().sumNeutralHadronEt);
    T_Muon_photonIsoR04->push_back(selected_muons[k].pfIsolationR04().sumPhotonEt);
    /*    T_Muon_smurfCharged->push_back(selected_muons[k].userFloat("smurfCharged"));
	  T_Muon_smurfNeutral->push_back(selected_muons[k].userFloat("smurfNeutral"));
	  T_Muon_smurfPhoton->push_back(selected_muons[k].userFloat("smurfPhoton"));
	  T_Muon_smurfNoOverCharged->push_back(selected_muons[k].userFloat("smurfNoOverCharged"));
	  T_Muon_smurfNoOverNeutral->push_back(selected_muons[k].userFloat("smurfNoOverNeutral"));
	  T_Muon_smurfNoOverPhoton->push_back(selected_muons[k].userFloat("smurfNoOverPhoton"));
	  T_Muon_muSmurfPF->push_back(selected_muons[k].userFloat("muSmurfPF"));
    */
    T_Muon_chargedParticleIsoR03->push_back(selected_muons[k].pfIsolationR03().sumChargedParticlePt);
    T_Muon_chargedHadronIsoR03->push_back(selected_muons[k].pfIsolationR03().sumChargedHadronPt);
    T_Muon_neutralHadronIsoR03->push_back(selected_muons[k].pfIsolationR03().sumNeutralHadronEt);
    T_Muon_photonIsoR03->push_back(selected_muons[k].pfIsolationR03().sumPhotonEt);
    T_Muon_sumPUPtR03->push_back(selected_muons[k].pfIsolationR03().sumPUPt);
    T_Muon_sumPUPtR04->push_back(selected_muons[k].pfIsolationR04().sumPUPt);
    T_Muon_vz->push_back(selected_muons[k].vz());
    T_Muon_vy->push_back(selected_muons[k].vy());
    T_Muon_vx->push_back(selected_muons[k].vx());
    T_Muon_NValidHitsGTrk->push_back(numOfValidHitsGTrk);
    T_Muon_SegmentCompatibility->push_back(muon::segmentCompatibility(selected_muons[k]));
    T_Muon_trkKink ->push_back(selected_muons[k].combinedQuality().trkKink );
    T_Muon_NValidHitsSATrk->push_back(nhitsouttrack);
    //T_Muon_EcalVeto->push_back(selected_muons[k].isolationR03().emVetoEt);
    //T_Muon_HcalVeto->push_back(selected_muons[k].isolationR03().hadVetoEt);
    T_Muon_NumOfMatchedStations->push_back(selected_muons[k].numberOfMatchedStations());
    T_Muon_isPFMuon->push_back(selected_muons[k].isPFMuon ());
    if(selected_muons[k].isPFMuon ())   {
      T_Muon_PFMuonPt->push_back(selected_muons[k].pfP4 ().pt());
      T_Muon_PFMuonPx->push_back(selected_muons[k].pfP4 ().px());
      T_Muon_PFMuonPy->push_back(selected_muons[k].pfP4 ().py());
      T_Muon_PFMuonPz->push_back(selected_muons[k].pfP4 ().pz());
      T_Muon_PFMuonE->push_back(selected_muons[k].pfP4 ().E());
    }
	
    else {
      T_Muon_PFMuonPt->push_back(-99999.);
      T_Muon_PFMuonPx->push_back(-99999.);
      T_Muon_PFMuonPy->push_back(-99999.);
      T_Muon_PFMuonPz->push_back(-99999.);
      T_Muon_PFMuonE->push_back(-99999.);
    }

    T_Muon_NLayers->push_back(nLayers);	
  }
  
  //************ PFTAUS ************************
  /*  T_pfTau_Px = new std::vector<float>;
      T_pfTau_Py = new std::vector<float>;
      T_pfTau_Pz = new std::vector<float>;
      T_pfTau_Energy = new std::vector<float>;
      T_pfTau_Charge = new std::vector<int>;
      //no estoy seguro de que haga falta ordenar (se puede usar para filtrar en el futuro)
      std::map<float,pat::Tau> pftauMap;
      for (size_t i = 0; i< pftauHandle->size(); ++i) {
      pftauMap[(*pftauHandle)[i].pt()]=(*pftauHandle)[i];
      }
      std::vector<pat::Tau> selected_pfTaus;
      for( std::map<float,pat::Tau>::reverse_iterator rit=pftauMap.rbegin(); rit!=pftauMap.rend(); ++rit){
      selected_pfTaus.push_back( (*rit).second );
      }
      for (size_t k = 0; k < selected_pfTaus.size(); ++k) {
      T_pfTau_Px->push_back(selected_pfTaus[k].px());
      T_pfTau_Py->push_back(selected_pfTaus[k].py());
      T_pfTau_Pz->push_back(selected_pfTaus[k].pz());
      T_pfTau_Energy->push_back(selected_pfTaus[k].energy());
      T_pfTau_Charge->push_back(selected_pfTaus[k].charge());
      }
  */  
  //************ TAUS ************************

  /*
    T_Tau_Px = new std::vector<float>;
    T_Tau_Py = new std::vector<float>;
    T_Tau_Pz = new std::vector<float>;
    T_Tau_Energy = new std::vector<float>;
    T_Tau_Charge = new std::vector<int>;
    //no estoy seguro de que haga falta ordenar (se puede usar para filtrar en el futuro)
    std::map<float,pat::Tau> tauMap;
    for (size_t i = 0; i< tauHandle->size(); ++i) {
    tauMap[(*tauHandle)[i].pt()]=(*tauHandle)[i];
    }
    std::vector<pat::Tau> selected_Taus;
    for( std::map<float,pat::Tau>::reverse_iterator rit=tauMap.rbegin(); rit!=tauMap.rend(); ++rit){
    selected_Taus.push_back( (*rit).second );
    }
    for (size_t k = 0; k < selected_Taus.size(); ++k) {
    T_Tau_Px->push_back(selected_Taus[k].px());
    T_Tau_Py->push_back(selected_Taus[k].py());
    T_Tau_Pz->push_back(selected_Taus[k].pz());
    T_Tau_Energy->push_back(selected_Taus[k].energy());
    T_Tau_Charge->push_back(selected_Taus[k].charge());
    }
  */  
  //************ PFMUONS ************************  
  
  //PfMuons
  
  /*       
	   T_pfMuon_IsGlobalMuon = new std::vector<bool>;
	   T_pfMuon_IsGMPTMuons = new std::vector<bool>;
  
	   T_pfMuon_IsAllStandAloneMuons = new std::vector<bool>;
	   T_pfMuon_IsAllTrackerMuons = new std::vector<bool>;
	   T_pfMuon_IsTMLastStationAngTight= new std::vector<bool>;
	   T_pfMuon_SegmentCompatibility= new std::vector<float>;
	   T_pfMuon_Px = new std::vector<float>;
	   T_pfMuon_Py = new std::vector<float>;
	   T_pfMuon_Pz = new std::vector<float>;
	   T_pfMuon_Pt = new std::vector<float>;
	   T_pfMuon_Energy = new std::vector<float>;
	   T_pfMuon_Charge = new std::vector<int>;
	   T_pfMuon_NormChi2GTrk = new std::vector<float>;
	   T_pfMuon_NValidHitsInTrk = new std::vector<int>;
	   T_pfMuon_NValidHitsSATrk = new std::vector<int>;
	   T_pfMuon_NValidPixelHitsInTrk = new std::vector<int>;
	   T_pfMuon_Chi2InTrk = new std::vector<float>;
	   T_pfMuon_dofInTrk = new std::vector<float>;
	   T_pfMuon_SumIsoCalo = new std::vector<float>;
	   T_pfMuon_SumIsoTrack = new std::vector<float>;
	   T_pfMuon_IPAbsGTrack = new std::vector<float>;
	   T_pfMuon_IPSigGTrack = new std::vector<float>;
	   T_pfMuon_IPAbsInTrack = new std::vector<float>;
	   T_pfMuon_IPSigInTrack = new std::vector<float>;
	   //  T_pfMuon_IPwrtBSInTrack =  new std::vector<float>;
	   T_pfMuon_IP2DBiasedPV  =  new std::vector<float>;
	   T_pfMuon_IP3DBiasedPV  =  new std::vector<float>;
	   T_pfMuon_IP2DUnBiasedPV  =  new std::vector<float>;
	   T_pfMuon_IP3DUnBiasedPV  =  new std::vector<float>;

	   T_pfMuon_vz = new std::vector<float>;
	   T_pfMuon_vy = new std::vector<float>;  
	   T_pfMuon_vx = new std::vector<float>;
	   T_pfMuon_NValidHits = new std::vector<int>;
	   T_pfMuon_particleIso = new std::vector<float>;
	   T_pfMuon_chargedHadronIso = new std::vector<float>;
	   T_pfMuon_neutralHadronIso = new std::vector<float>;
	   T_pfMuon_photonIso = new std::vector<float>;
	   T_pfMuon_deltaPt = new std::vector<float>;
	   T_pfMuon_NumOfMatchedStations = new std::vector<int>;
	   T_pfMuon_pfCharged  = new std::vector<float>;
	   T_pfMuon_pfNeutral  = new std::vector<float>;
	   T_pfMuon_pfPhoton = new std::vector<float>;
	   T_pfMuon_smurfCharged = new std::vector<float>;
	   T_pfMuon_smurfNeutral = new std::vector<float>;
	   T_pfMuon_smurfPhoton = new std::vector<float>; 
	   T_pfMuon_smurfNoOverCharged = new std::vector<float>;
	   T_pfMuon_smurfNoOverNeutral = new std::vector<float>;
	   T_pfMuon_smurfNoOverPhoton = new std::vector<float>;
	   T_pfMuon_muSmurfPF = new std::vector<float>;
	   T_pfMuon_sumPUPt = new std::vector<float>;  
	   T_pfMuon_sumPUPtR03 = new std::vector<float>;
	   T_pfMuon_NLayers = new std::vector<int>;

	   //no estoy seguro de que haga falta ordenar (se puede usar para filtrar en el futuro)
	   std::map<float,pat::Muon> pfmuonMap;
	   for (size_t i = 0; i< pfmuonHandle->size(); ++i) {
	   pfmuonMap[(*pfmuonHandle)[i].pt()]=(*pfmuonHandle)[i];  
	   }
	   std::vector<pat::Muon> selected_pfMuons;
	   for( std::map<float,pat::Muon>::reverse_iterator rit=pfmuonMap.rbegin(); rit!=pfmuonMap.rend(); ++rit){
	   selected_pfMuons.push_back( (*rit).second );
	   }  
	   for (size_t k = 0; k < selected_pfMuons.size(); ++k) {
    
	   float IP      = 9999;
	   float IPSignificance = 9999;
	   float normchi2 = 9999;
	   //quantities wrt global track    
	   reco::TrackRef tr_globaltrack = selected_pfMuons[k].globalTrack();       
	   if (!tr_globaltrack.isNull() && selected_pfMuons[k].isGlobalMuon()) {
    
	   normchi2 = selected_pfMuons[k].globalTrack()->normalizedChi2();
	   if (vtxs.size() > 0) {
	
	   IPTools::ImpactParameterComputer IPComp(vtxs[0]); // use the vtx with the highest pt sum in this example 
	   Measurement1D mess1D = IPComp.computeIP(iSetup, *tr_globaltrack);
	
	   // to calculate the ImpactParameter, its Error and the Significance use
	   IP             =  mess1D.value();
	
	   //IPError        =  mess1D.error();
	   IPSignificance =  mess1D.significance();
	
	   }
      
    
	   }
	   //cuantities wrt inner track
	   reco::TrackRef tr_innertrack = selected_pfMuons[k].innerTrack(); 
	   int nhitsinnertracker = -1, pixelHits = -1, nLayers=-1;
	   float chi2innertracker=9999;
	   float dofinnertracker=9999;
	   float IPIn      = 9999;
	   float IPSignificanceIn = 9999, deltaPt=9999;
	   //    float IPwrtBSIn = 9999.;
	   if (!tr_innertrack.isNull()) {
	   if (vtxs.size() > 0) {
	
	   IPTools::ImpactParameterComputer IPComp(vtxs[0]); // use the vtx with the highest pt sum in this example 
	   Measurement1D mess1D = IPComp.computeIP(iSetup, *tr_innertrack);
	
	   // to calculate the ImpactParameter, its Error and the Significance use
	   IPIn             =  mess1D.value();
	   IPSignificanceIn =  mess1D.significance();
	
	   }
	   nhitsinnertracker = tr_innertrack->hitPattern().numberOfValidTrackerHits();
	   pixelHits = tr_innertrack->hitPattern().numberOfValidPixelHits();
	   nLayers=tr_innertrack->hitPattern().trackerLayersWithMeasurement() ;
	   //      nhitsinnertracker = tr_innertrack->numberOfValidHits();
	   chi2innertracker=tr_innertrack->chi2();
	   dofinnertracker=tr_innertrack->ndof();
	   //IPwrtBSIn=tr_innertrack->dxy(beamSpot.position()); 
	   deltaPt=tr_innertrack->ptError();
	   }
    
	   int numOfValidHits=0;
    
	   if (selected_pfMuons[k].isGlobalMuon()){
    
	   numOfValidHits=selected_pfMuons[k].numberOfValidHits();
    
	   }
    
    
	   reco::TrackRef tr_outtrack = selected_pfMuons[k].standAloneMuon(); 
    
	   float nhitsouttrack=9999;
    
	   if (!tr_outtrack.isNull()) {
    
	   nhitsouttrack=selected_pfMuons[k].standAloneMuon()->numberOfValidHits();
    
	   }
    
	   T_pfMuon_IsGlobalMuon->push_back(selected_pfMuons[k].isGlobalMuon());
	   T_pfMuon_IsGMPTMuons->push_back(selected_pfMuons[k].muonID("GlobalMuonPromptTight"));
	   T_pfMuon_IsAllStandAloneMuons->push_back(selected_pfMuons[k].muonID("AllStandAloneMuons"));

	   T_pfMuon_IsAllTrackerMuons->push_back(selected_pfMuons[k].muonID("AllTrackerMuons"));
	   T_pfMuon_IsTMLastStationAngTight->push_back(selected_pfMuons[k].muonID("TMLastStationAngTight"));

	   T_pfMuon_Px->push_back(selected_pfMuons[k].px());
	   T_pfMuon_Py->push_back(selected_pfMuons[k].py());
	   T_pfMuon_Pz->push_back(selected_pfMuons[k].pz());
	   T_pfMuon_Pt->push_back(selected_pfMuons[k].pt());
	   T_pfMuon_Energy->push_back(selected_pfMuons[k].energy());
	   T_pfMuon_Charge->push_back(selected_pfMuons[k].charge());
	   T_pfMuon_NormChi2GTrk->push_back(normchi2);
	   T_pfMuon_NValidHitsInTrk->push_back(nhitsinnertracker);
	   T_pfMuon_Chi2InTrk->push_back(chi2innertracker);
	   T_pfMuon_dofInTrk->push_back(dofinnertracker);
	   T_pfMuon_SumIsoCalo->push_back(selected_pfMuons[k].caloIso());
	   T_pfMuon_SumIsoTrack->push_back(selected_pfMuons[k].trackIso());
	   T_pfMuon_IPAbsGTrack->push_back(IP);
	   T_pfMuon_IPSigGTrack->push_back(IPSignificance);
	   T_pfMuon_IPAbsInTrack->push_back(IPIn);
	   T_pfMuon_IPSigInTrack->push_back(IPSignificanceIn);
	   //    T_pfMuon_IPwrtBSInTrack->push_back(IPwrtBSIn);
	   T_pfMuon_IP2DBiasedPV->push_back(selected_pfMuons[k].userFloat("tip"));
	   T_pfMuon_IP3DBiasedPV->push_back(selected_pfMuons[k].userFloat("ip"));
	   T_pfMuon_IP2DUnBiasedPV->push_back(selected_pfMuons[k].userFloat("tip2"));
	   T_pfMuon_IP3DUnBiasedPV->push_back(selected_pfMuons[k].userFloat("ip2"));

	   T_pfMuon_vz->push_back(selected_pfMuons[k].vz());
	   T_pfMuon_vy->push_back(selected_pfMuons[k].vy());
	   T_pfMuon_vx->push_back(selected_pfMuons[k].vx());
	   T_pfMuon_NValidHits->push_back(numOfValidHits);
	   T_pfMuon_NValidPixelHitsInTrk->push_back(pixelHits);
	   T_pfMuon_SegmentCompatibility->push_back(muon::segmentCompatibility(selected_pfMuons[k]));
	   T_pfMuon_NValidHitsSATrk->push_back(nhitsouttrack);
	   T_pfMuon_pfCharged->push_back(selected_pfMuons[k].userFloat("pfCharged"));
	   T_pfMuon_pfNeutral->push_back(selected_pfMuons[k].userFloat("pfNeutral"));
	   T_pfMuon_pfPhoton->push_back(selected_pfMuons[k].userFloat("pfPhoton"));
	   T_pfMuon_smurfCharged->push_back(selected_pfMuons[k].userFloat("smurfCharged"));
	   T_pfMuon_smurfNeutral->push_back(selected_pfMuons[k].userFloat("smurfNeutral"));
	   T_pfMuon_smurfPhoton->push_back(selected_pfMuons[k].userFloat("smurfPhoton"));
	   T_pfMuon_smurfNoOverCharged->push_back(selected_pfMuons[k].userFloat("smurfNoOverCharged"));
	   T_pfMuon_smurfNoOverNeutral->push_back(selected_pfMuons[k].userFloat("smurfNoOverNeutral"));
	   T_pfMuon_smurfNoOverPhoton->push_back(selected_pfMuons[k].userFloat("smurfNoOverPhoton"));
	   T_pfMuon_muSmurfPF->push_back(selected_pfMuons[k].userFloat("muSmurfPF"));
	   T_pfMuon_sumPUPt->push_back(selected_pfMuons[k].pfIsolationR04().sumPUPt);
	   T_pfMuon_sumPUPtR03->push_back(selected_pfMuons[k].pfIsolationR03().sumPUPt);

	   T_pfMuon_particleIso->push_back(selected_pfMuons[k].particleIso());
	   T_pfMuon_chargedHadronIso->push_back(selected_pfMuons[k].chargedHadronIso());
	   T_pfMuon_neutralHadronIso->push_back(selected_pfMuons[k].neutralHadronIso());
	   T_pfMuon_photonIso->push_back(selected_pfMuons[k].photonIso());

	   T_pfMuon_deltaPt->push_back(deltaPt);
	   T_pfMuon_NumOfMatchedStations->push_back(selected_pfMuons[k].numberOfMatchedStations());
	   T_pfMuon_NLayers->push_back(nLayers);
	   }
  */
  //**** END OF PFMUONS************
  //Electrons
  T_Elec_Eta  = new std::vector<float>;
  T_Elec_IPwrtAveBS = new std::vector<float>;
  T_Elec_IPwrtPV = new std::vector<float>;
  T_Elec_dzwrtPV = new std::vector<float>;
  //  T_Elec_IPwrtBS = new std::vector<float>;
  T_Elec_Px = new std::vector<float>;
  T_Elec_Py = new std::vector<float>;
  T_Elec_Pz = new std::vector<float>;
  T_Elec_Pt = new std::vector<float>;
  T_Elec_Energy = new std::vector<float>;
  T_Elec_Charge = new std::vector<int>;
  /*  T_Elec_SumIsoCalo = new std::vector<float>;
      T_Elec_SumIsoTrack = new std::vector<float>;*/
  /*  T_Elec_IP2DBiasedPV = new std::vector<float>;
      T_Elec_IP3DBiasedPV = new std::vector<float>;
      T_Elec_IP2DUnBiasedPV = new std::vector<float>;
      T_Elec_IP3DUnBiasedPV = new std::vector<float>;
      T_Elec_dxyPVBiasedPV = new std::vector<float>;
      T_Elec_dzPVBiasedPV = new std::vector<float>;
      T_Elec_dxyPVUnBiasedPV = new std::vector<float>;
      T_Elec_dzPVUnBiasedPV = new std::vector<float>;
      T_Elec_IP2DUnBiasedPVnoBS = new std::vector<float>;
      T_Elec_IP3DUnBiasedPVnoBS = new std::vector<float>;
      T_Elec_dxyPVUnBiasedPVnoBS = new std::vector<float>;
      T_Elec_dzPVUnBiasedPVnoBS = new std::vector<float>;*/
  T_Elec_vz = new std::vector<float>;
  T_Elec_vy = new std::vector<float>;  
  T_Elec_vx = new std::vector<float>;
  T_Elec_nLost =  new std::vector<int>; 
  T_Elec_nHits  =  new std::vector<int>;
  T_Elec_SC_Et = new std::vector<float>;
  T_Elec_SC_Eta = new std::vector<float>; 
  T_Elec_nBrems  = new std::vector<int>;
  T_Elec_fBrem = new std::vector<float>;
  T_Elec_eSuperClusterOverP = new std::vector<float>;
  T_Elec_ecalEnergy = new std::vector<float>;
  /*T_Elec_dr03TkSumPt = new std::vector<float>; 
    T_Elec_dr03EcalSumEt = new std::vector<float>; 
    T_Elec_dr03HcalSumEt = new std::vector<float>; 
    T_Elec_ConvInfoDist = new std::vector<float>; 
    T_Elec_ConvInfoDCot = new std::vector<float>; */
  T_Elec_passConversionVeto=new std::vector<bool>;
  T_Elec_isEB = new std::vector<bool>;
  T_Elec_isEE = new std::vector<bool>;
  T_Elec_MVA = new std::vector<float>;
  /*  T_Elec_simpleEleId95 = new std::vector<float>;
      T_Elec_simpleEleId90 = new std::vector<float>;
      T_Elec_simpleEleId85 = new std::vector<float>;
  */  T_Elec_simpleEleId80 = new std::vector<float>;
  /*  T_Elec_simpleEleId70 = new std::vector<float>;
      T_Elec_simpleEleId60 = new std::vector<float>;*/
  //  T_Elec_simpleEleId90cIso = new std::vector<float>;
  /* T_Elec_cicVeryLooseHWW = new std::vector<float>;
     T_Elec_cicLooseHWW = new std::vector<float>;
     T_Elec_cicMediumHWW = new std::vector<float>;
     T_Elec_cicSuperTightHWW  = new std::vector<float>;
     T_Elec_cicHyperTight1HWW  = new std::vector<float>;
     T_Elec_cicHyperTight2HWW  = new std::vector<float>;
     T_Elec_cicHyperTight3HWW  = new std::vector<float>;
     T_Elec_cicVeryLooseMC = new std::vector<float>;
     T_Elec_cicLooseMC = new std::vector<float>;
     T_Elec_cicMediumMC  = new std::vector<float>;
     T_Elec_cicSuperTightMC  = new std::vector<float>;
     T_Elec_cicHyperTight1MC  = new std::vector<float>;
     T_Elec_cicHyperTight2MC  = new std::vector<float>;
     T_Elec_cicHyperTight3MC  = new std::vector<float>;  
     T_Elec_cicVeryLoose = new std::vector<float>;
     T_Elec_cicLoose = new std::vector<float>;
     T_Elec_cicMedium  = new std::vector<float>;
     T_Elec_cicSuperTight  = new std::vector<float>;
     T_Elec_cicHyperTight1  = new std::vector<float>;
     T_Elec_cicHyperTight2  = new std::vector<float>;
     T_Elec_cicHyperTight3  = new std::vector<float>;*/
  /*T_Elec_Loose = new std::vector<float>;
    T_Elec_RobustLoose  = new std::vector<float>;
    T_Elec_Tight  = new std::vector<float>;
    T_Elec_RobustTight  = new std::vector<float>;
    T_Elec_RobustHighEnergy  = new std::vector<float>;
    T_Elec_egammaIDLikelihood = new std::vector<float>;*/
  T_Elec_sigmaIetaIeta = new std::vector<float>;
  T_Elec_deltaPhiIn = new std::vector<float>;
  T_Elec_deltaEtaIn =new  std::vector<float>;
  T_Elec_isEcalDriven = new std::vector<bool>;
  T_Elec_HtoE = new std::vector<float>;
  //  T_Elec_passTriggerDoubleEl = new std::vector<bool>;
  //  T_Elec_passTriggerElMu  = new std::vector<bool>;
  T_Elec_chargedHadronIso= new std::vector<float>;
  T_Elec_neutralHadronIso= new std::vector<float>;
  T_Elec_photonIso= new std::vector<float>;
  T_Elec_puChargedHadronIso = new std::vector<float>;
  T_Elec_isPF = new std::vector<bool>;
  T_Elec_PFElecPt = new std::vector<float>;
  T_Elec_PFElecPx = new std::vector<float>;
  T_Elec_PFElecPy = new std::vector<float>;
  T_Elec_PFElecPz = new std::vector<float>;
  T_Elec_PFElecE  = new std::vector<float>;
  /*   T_Elec_smurfCharged= new std::vector<float>;
       T_Elec_smurfNeutral= new std::vector<float>;
       T_Elec_smurfPhoton= new std::vector<float>;
       T_Elec_smurfNoOverCharged= new std::vector<float>;
       T_Elec_smurfNoOverNeutral= new std::vector<float>;
       T_Elec_smurfNoOverPhoton= new std::vector<float>;
       T_Elec_eleSmurfPF= new std::vector<float>;*/

  std::map<float,pat::Electron> electronMap;
  for (size_t j = 0; j < electronHandle->size(); ++j) {
    electronMap[(*electronHandle)[j].pt()]=(*electronHandle)[j];   
  } 

  std::vector<pat::Electron> selected_electrons;
  for (std::map<float,pat::Electron>::reverse_iterator rit=electronMap.rbegin(); rit!=electronMap.rend(); ++rit) {
    selected_electrons.push_back( (*rit).second );
  }
 
  for (size_t k = 0; k < selected_electrons.size(); ++k) {
    /*    leps_.push_back(selected_electrons[k]);      
	  T_Elec_passTriggerDoubleEl->push_back(passTriggerDoubleEl(selected_muons.size()+k,IsRealData));
	  T_Elec_passTriggerElMu->push_back(passTriggerElMu(selected_muons.size()+k,IsRealData));
    */
    float IP = 9999.;
    float dz = 9999.;
    //float IPError = -9999;
    //    float IPSignificance = 9999;
    int nLost = 9999, nHits=-9999;
    //    float D0= 9999; 
    reco::GsfTrackRef trRef_elec = (*electronHandle)[k].gsfTrack();
    if (!trRef_elec.isNull()) {
      if (vtxs.size() > 0) {
	//        IPTools::ImpactParameterComputer IPComp(vtxs[0]); // use the vtx with the highest pt sum in this example 
	//	Measurement1D mess1D = IPComp.computeIP(iSetup, *trRef_elec);
	
	// to calculate the ImpactParameter, its Error and the Significance use
	IP             = fabs(trRef_elec->dxy(vtxs[0].position())) ;
	dz             = fabs(trRef_elec->dz(vtxs[0].position())) ;
	//	IPSignificance =  mess1D.significance();
      }
     	
      
      nHits = trRef_elec->trackerExpectedHitsInner().numberOfHits();
      nLost= trRef_elec->trackerExpectedHitsInner().numberOfLostHits(); 
      //     D0 = selected_electrons[k].gsfTrack()->dxy(beamSpot.position());
    }

    T_Elec_Eta->push_back(selected_electrons[k].eta());
    T_Elec_IPwrtAveBS->push_back(selected_electrons[k].dB());
    T_Elec_IPwrtPV->push_back(IP); 
    T_Elec_dzwrtPV->push_back(dz); 
    T_Elec_Px->push_back(selected_electrons[k].px());
    T_Elec_Py->push_back(selected_electrons[k].py());
    T_Elec_Pz->push_back(selected_electrons[k].pz());
    T_Elec_Pt->push_back(selected_electrons[k].pt());
    T_Elec_Energy->push_back(selected_electrons[k].energy());
    T_Elec_Charge->push_back(selected_electrons[k].charge());
    /*    T_Elec_SumIsoCalo->push_back(selected_electrons[k].caloIso());
	  T_Elec_SumIsoTrack->push_back(selected_electrons[k].trackIso());*/
    T_Elec_nBrems->push_back(selected_electrons[k].numberOfBrems());
    T_Elec_fBrem->push_back(selected_electrons[k].fbrem());
    T_Elec_eSuperClusterOverP->push_back(selected_electrons[k].eSuperClusterOverP());
    T_Elec_ecalEnergy->push_back(selected_electrons[k].ecalEnergy());
    /*    T_Elec_dr03TkSumPt->push_back(selected_electrons[k].dr03TkSumPt());
	  T_Elec_dr03EcalSumEt->push_back(selected_electrons[k].dr03EcalRecHitSumEt());
	  T_Elec_dr03HcalSumEt->push_back(selected_electrons[k].dr03HcalTowerSumEt());
	  T_Elec_IP2DBiasedPV->push_back(selected_electrons[k].userFloat("tip"));
	  T_Elec_IP3DBiasedPV->push_back(selected_electrons[k].userFloat("ip"));
	  T_Elec_IP2DUnBiasedPV->push_back(selected_electrons[k].userFloat("tip2"));
	  T_Elec_IP3DUnBiasedPV->push_back(selected_electrons[k].userFloat("ip2"));
	  T_Elec_dxyPVBiasedPV->push_back(selected_electrons[k].userFloat("dxyPV"));
	  T_Elec_dzPVBiasedPV->push_back(selected_electrons[k].userFloat("dzPV"));
	  T_Elec_dxyPVUnBiasedPV->push_back(selected_electrons[k].userFloat("dxyPV2"));
	  T_Elec_dzPVUnBiasedPV->push_back(selected_electrons[k].userFloat("dzPV2"));
	  T_Elec_IP2DUnBiasedPVnoBS->push_back(selected_electrons[k].userFloat("tip3"));
	  T_Elec_IP3DUnBiasedPVnoBS->push_back(selected_electrons[k].userFloat("ip3"));
	  T_Elec_dxyPVUnBiasedPVnoBS->push_back(selected_electrons[k].userFloat("dxyPV3"));
	  T_Elec_dzPVUnBiasedPVnoBS->push_back(selected_electrons[k].userFloat("dzPV3"));
    */    T_Elec_vz->push_back(selected_electrons[k].vz());
    T_Elec_vy->push_back(selected_electrons[k].vy());
    T_Elec_vx->push_back(selected_electrons[k].vx());	
    T_Elec_nLost->push_back(nLost); 
    T_Elec_nHits->push_back(nHits);
    T_Elec_SC_Et->push_back( selected_electrons[k].superCluster()->energy()/TMath::CosH(selected_electrons[k].superCluster()->eta()));
    T_Elec_SC_Eta->push_back( selected_electrons[k].superCluster()->eta());
    T_Elec_chargedHadronIso->push_back(selected_electrons[k].chargedHadronIso());
    T_Elec_neutralHadronIso->push_back(selected_electrons[k].neutralHadronIso());
    T_Elec_photonIso->push_back(selected_electrons[k].photonIso());
    T_Elec_puChargedHadronIso->push_back(selected_electrons[k].puChargedHadronIso());
    T_Elec_passConversionVeto->push_back(selected_electrons[k].passConversionVeto());
    //    T_Elec_simpleEleId95 ->push_back(selected_electrons[k].electronID("vbtf11WP95"));
    //    T_Elec_simpleEleId90 ->push_back(selected_electrons[k].electronID("vbtf11WP90"));
    //    T_Elec_simpleEleId85 ->push_back(selected_electrons[k].electronID("vbtf11WP85"));
    T_Elec_simpleEleId80 ->push_back(selected_electrons[k].electronID("vbtf11WP80"));
    //    T_Elec_simpleEleId70 ->push_back(selected_electrons[k].electronID("vbtf11WP70"));
    //    T_Elec_simpleEleId60 ->push_back(selected_electrons[k].electronID("vbtf11WP60"));

    //    T_Elec_simpleEleId90cIso ->push_back(selected_electrons[k].electronID("simpleEleId90cIso"));

    /*    T_Elec_cicVeryLooseHWW ->push_back(selected_electrons[k].electronID("cicVeryLooseHWW"));
	  T_Elec_cicLooseHWW ->push_back(selected_electrons[k].electronID("cicLooseHWW"));
	  T_Elec_cicMediumHWW ->push_back(selected_electrons[k].electronID("cicMediumHWW"));
	  T_Elec_cicSuperTightHWW  ->push_back(selected_electrons[k].electronID("cicSuperTightHWW"));
	  T_Elec_cicHyperTight1HWW  ->push_back(selected_electrons[k].electronID("cicHyperTight1HWW"));
	  T_Elec_cicHyperTight2HWW  ->push_back(selected_electrons[k].electronID("cicHyperTight2HWW"));
	  T_Elec_cicHyperTight3HWW  ->push_back(selected_electrons[k].electronID("cicHyperTight3HWW"));
	  T_Elec_cicVeryLooseMC ->push_back(selected_electrons[k].electronID("cicVeryLooseMC"));
	  T_Elec_cicLooseMC ->push_back(selected_electrons[k].electronID("cicLooseMC"));
	  T_Elec_cicMediumMC  ->push_back(selected_electrons[k].electronID("cicMediumMC"));
	  T_Elec_cicSuperTightMC  ->push_back(selected_electrons[k].electronID("cicSuperTightMC"));
	  T_Elec_cicHyperTight1MC  ->push_back(selected_electrons[k].electronID("cicHyperTight1MC"));
	  T_Elec_cicHyperTight2MC  ->push_back(selected_electrons[k].electronID("cicHyperTight2MC"));
	  T_Elec_cicHyperTight3MC  ->push_back(selected_electrons[k].electronID("cicHyperTight3MC"));
	  T_Elec_cicVeryLoose ->push_back(selected_electrons[k].electronID("cicVeryLoose"));
	  T_Elec_cicLoose ->push_back(selected_electrons[k].electronID("cicLoose"));
	  T_Elec_cicMedium  ->push_back(selected_electrons[k].electronID("cicMedium"));
	  T_Elec_cicSuperTight  ->push_back(selected_electrons[k].electronID("cicSuperTight"));
	  T_Elec_cicHyperTight1  ->push_back(selected_electrons[k].electronID("cicHyperTight1"));
	  T_Elec_cicHyperTight2  ->push_back(selected_electrons[k].electronID("cicHyperTight2"));
	  T_Elec_cicHyperTight3  ->push_back(selected_electrons[k].electronID("cicHyperTight3"));*/
    /*T_Elec_Loose ->push_back(selected_electrons[k].electronID("eidLoose"));
      T_Elec_RobustLoose  ->push_back(selected_electrons[k].electronID("eidRobustLoose"));
      T_Elec_Tight  ->push_back(selected_electrons[k].electronID("eidTight"));
      T_Elec_RobustTight  ->push_back(selected_electrons[k].electronID("eidRobustTight"));
      T_Elec_RobustHighEnergy  ->push_back(selected_electrons[k].electronID("eidRobustHighEnergy"));
      T_Elec_egammaIDLikelihood ->push_back(selected_electrons[k].electronID("egammaIDLikelihood"));*/

    T_Elec_sigmaIetaIeta ->push_back( selected_electrons[k].sigmaIetaIeta());
    T_Elec_deltaPhiIn->push_back( selected_electrons[k].deltaPhiSuperClusterTrackAtVtx());
    T_Elec_deltaEtaIn->push_back( selected_electrons[k].deltaEtaSuperClusterTrackAtVtx());
    T_Elec_isEcalDriven -> push_back(selected_electrons[k].ecalDrivenSeed());
    T_Elec_HtoE ->push_back(selected_electrons[k].hadronicOverEm());
    //T_Elec_ConvInfoDist->push_back(selected_electrons[k].userFloat("convValueMapProd:dist"));
   
    //T_Elec_ConvInfoDCot->push_back(selected_electrons[k].userFloat("convValueMapProd:dcot"));
   
    T_Elec_MVA->push_back(selected_electrons[k].electronID("mvaTrigV0")); 
    T_Elec_isEB->push_back(selected_electrons[k].isEB());
    T_Elec_isEE->push_back(selected_electrons[k].isEE());
    T_Elec_isPF->push_back(selected_electrons[k].isPF());
    reco::GsfElectron::P4Kind pf = reco::GsfElectron::P4_PFLOW_COMBINATION;
    if(selected_electrons[k].isPF ())   {
      T_Elec_PFElecPt->push_back(selected_electrons[k].p4 (pf).pt());
      T_Elec_PFElecPx->push_back(selected_electrons[k].p4 (pf).px());
      T_Elec_PFElecPy->push_back(selected_electrons[k].p4 (pf).py());
      T_Elec_PFElecPz->push_back(selected_electrons[k].p4 (pf).pz());
      T_Elec_PFElecE->push_back(selected_electrons[k].p4 (pf).E());
    }

    else {
      T_Elec_PFElecPt->push_back(-99999.);
      T_Elec_PFElecPx->push_back(-99999.);
      T_Elec_PFElecPy->push_back(-99999.);
      T_Elec_PFElecPz->push_back(-99999.);
      T_Elec_PFElecE->push_back(-99999.);
    }

  }
  
  //**************** PFELECTRONS  **********************
  /*  
  //pfElectrons
  T_pfElec_Px = new std::vector<float>;
  T_pfElec_Py = new std::vector<float>;
  T_pfElec_Pz = new std::vector<float>;
  T_pfElec_Pt = new std::vector<float>;
  T_pfElec_Energy = new std::vector<float>;
  T_pfElec_Charge = new std::vector<int>;
  T_pfElec_SumIsoCalo = new std::vector<float>;
  T_pfElec_SumIsoTrack = new std::vector<float>;
  T_pfElec_IP2DBiasedPV = new std::vector<float>;
  T_pfElec_IP3DBiasedPV = new std::vector<float>;
  T_pfElec_IP2DUnBiasedPV = new std::vector<float>;
  T_pfElec_IP3DUnBiasedPV = new std::vector<float>;

  T_pfElec_vz = new std::vector<float>;
  T_pfElec_vy = new std::vector<float>;  
  T_pfElec_vx = new std::vector<float>;
  T_pfElec_nBrems  = new std::vector<int>;
  T_pfElec_dr03TkSumPt = new std::vector<float>;
  T_pfElec_dr03EcalSumEt = new std::vector<float>;
  T_pfElec_dr03HcalSumEt = new std::vector<float>;
  T_pfElec_SC_Et = new std::vector<float>;
  T_pfElec_SC_Eta = new std::vector<float>;
  T_pfElec_isEcalDriven = new std::vector<bool>;
  T_pfElec_HtoE = new std::vector<float>;
  T_pfElec_particleIso = new std::vector<float>;
  T_pfElec_chargedHadronIso = new std::vector<float>;
  T_pfElec_neutralHadronIso = new std::vector<float>;
  T_pfElec_photonIso = new std::vector<float>;
  T_pfElec_ConvInfoDCot=new std::vector<float>;
  T_pfElec_ConvInfoDist=new std::vector<float>;
  T_pfElec_nHits = new std::vector<int>;
  
  std::map<float,pat::Electron> pfelectronMap;
  for (size_t j = 0; j < pfelectronHandle->size(); ++j) {
  pfelectronMap[(*pfelectronHandle)[j].pt()]=(*pfelectronHandle)[j];
  } 
  std::vector<pat::Electron> selected_pfelectrons;
  for (std::map<float,pat::Electron>::reverse_iterator rit=pfelectronMap.rbegin(); rit!=pfelectronMap.rend(); ++rit) {
  selected_pfelectrons.push_back( (*rit).second );
  }
  
  
  for (size_t k = 0; k < selected_pfelectrons.size(); ++k) {
    
  //    float IP      = 9999;
  //float IPError = -9999;
  //    float IPSignificance = 9999;
  int nHits=-999;
  reco::GsfTrackRef trRef_pfelec = (*pfelectronHandle)[k].gsfTrack();
  if (!trRef_pfelec.isNull()) {

  nHits = trRef_pfelec->trackerExpectedHitsInner().numberOfHits();    
  }
    
  T_pfElec_Px->push_back(selected_pfelectrons[k].px());
  T_pfElec_Py->push_back(selected_pfelectrons[k].py());
  T_pfElec_Pz->push_back(selected_pfelectrons[k].pz());
  T_pfElec_Pt->push_back(selected_pfelectrons[k].pt());
  T_pfElec_Energy->push_back(selected_pfelectrons[k].energy());
  T_pfElec_Charge->push_back(selected_pfelectrons[k].charge());
  T_pfElec_SumIsoCalo->push_back(selected_pfelectrons[k].caloIso());
  T_pfElec_SumIsoTrack->push_back(selected_pfelectrons[k].trackIso());
  T_pfElec_IP2DBiasedPV->push_back(selected_pfelectrons[k].userFloat("tip"));
  T_pfElec_IP3DBiasedPV->push_back(selected_pfelectrons[k].userFloat("ip"));
  T_pfElec_IP2DUnBiasedPV->push_back(selected_pfelectrons[k].userFloat("tip2"));
  T_pfElec_IP3DUnBiasedPV->push_back(selected_pfelectrons[k].userFloat("ip2"));
  T_pfElec_vz->push_back(selected_pfelectrons[k].vz());
  T_pfElec_vy->push_back(selected_pfelectrons[k].vy());
  T_pfElec_vx->push_back(selected_pfelectrons[k].vx());	
  T_pfElec_nBrems->push_back(selected_pfelectrons[k].numberOfBrems());
  T_pfElec_dr03TkSumPt->push_back(selected_pfelectrons[k].dr03TkSumPt());
  T_pfElec_dr03EcalSumEt->push_back(selected_pfelectrons[k].dr03EcalRecHitSumEt());
  T_pfElec_dr03HcalSumEt->push_back(selected_pfelectrons[k].dr03HcalTowerSumEt());        
  T_pfElec_SC_Et->push_back( selected_pfelectrons[k].superCluster()->energy()/TMath::CosH(selected_pfelectrons[k].superCluster()->eta()));
  T_pfElec_SC_Eta->push_back( selected_pfelectrons[k].superCluster()->eta());
  T_pfElec_nHits->push_back(nHits );
  T_pfElec_isEcalDriven -> push_back(selected_pfelectrons[k].ecalDrivenSeed());
  T_pfElec_HtoE->push_back(selected_electrons[k].hadronicOverEm());
  T_pfElec_particleIso->push_back(selected_pfelectrons[k].particleIso());
  T_pfElec_chargedHadronIso->push_back(selected_pfelectrons[k].chargedHadronIso());
  T_pfElec_neutralHadronIso->push_back(selected_pfelectrons[k].neutralHadronIso());
  T_pfElec_photonIso->push_back(selected_pfelectrons[k].photonIso());

  T_pfElec_ConvInfoDist->push_back(selected_pfelectrons[k].userFloat("convValueMapProd:dist"));
   
  T_pfElec_ConvInfoDCot->push_back(selected_pfelectrons[k].userFloat("convValueMapProd:dcot"));


  }
  */
  //************ END OF PFELECTRONS  **************
  
  //  SetJetInfo(0, jets, vtxs, false);
  SetJetInfo(0, jetsPF, vtxs, false);
  //  SetJetInfo(2, jetsPF2, vtxs, false);
  
  
  //MET
  /*  T_MET_ET = mets[0].pt();
      T_MET_Phi = mets[0].phi();
      T_MET_Sig = mets[0].significance();
  */    
  T_METPF_ET = metPF.pt();
  T_METPF_Phi = metPF.phi();
  T_METPF_Sig = metPF.significance();

  T_METPFTypeI_ET = metPFTypeI.pt();
  T_METPFTypeI_Phi = metPFTypeI.phi();
  T_METPFTypeI_Sig = metPFTypeI.significance();

  /*  T_METChargedNeutralPFNoPU_ET = (*vmH)[reco::VertexRef(vtxH,0)].pt();
      T_METChargedNeutralPFNoPU_Phi = (*vmH)[reco::VertexRef(vtxH,0)].phi();

      T_METChargedPFNoPU_ET = (*vmH2)[reco::VertexRef(vtxH,0)].pt();
      T_METChargedPFNoPU_Phi = (*vmH2)[reco::VertexRef(vtxH,0)].phi();
  */
  /* 
     T_METtc_ET = metsTC[0].pt();
     T_METtc_Phi = metsTC[0].phi();
     T_METtc_Sig = metsTC[0].significance();
  */
  /*
    edm::Handle< reco::PFMET > assocPfMetHandle;
    iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd","assocPfMet"),assocPfMetHandle);
    //   assocPfMet - built from the sum of the candidates associated to the primary vertex
    T_MET_assocPfMet_ET = assocPfMetHandle->pt();
    T_MET_assocPfMet_Phi = assocPfMetHandle->phi();

    edm::Handle< reco::PFMET > trkPfMetHandle;
    iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd","trkPfMet"),trkPfMetHandle);
    //    trkPfMet - built from the sum of all charged candidates associated to the primary vertex
    T_MET_trkPfMet_ET= trkPfMetHandle->pt();
    T_MET_trkPfMet_Phi = trkPfMetHandle->phi();

    edm::Handle< reco::PFMET > assocOtherVtxPfMetHandle;
    iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd","assocOtherVtxPfMet"),assocOtherVtxPfMetHandle);
    //    assocOtherVtxPfMet - built from the sum of the candidates associated to all the secondary vertices
    T_MET_assocOtherVtxPfMet_ET = assocOtherVtxPfMetHandle->pt();
    T_MET_assocOtherVtxPfMet_Phi = assocOtherVtxPfMetHandle->phi();

    //edm::Handle< edm::View<reco::MET> > globalPfMetHandle;
    //iEvent.getByLabel("globalPfMetHandle",globalPfMetHandle);
    //    globalPfMet - built from all candidates

    edm::Handle< reco::PFMET > centralPfMetHandle;
    iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd","centralPfMet"),centralPfMetHandle);
    //   centralPfMet - built from candidates in the central region of the detector (i.e. $|\eta|<2.4$)
    T_MET_centralPfMet_ET = centralPfMetHandle->pt();
    T_MET_centralPfMet_Phi = centralPfMetHandle->phi();

    edm::Handle< reco::PFMET > cleanPfMetHandle;
    iEvent.getByLabel(edm::InputTag("ClusteredPFMetProducerStd","cleanPfMet"),cleanPfMetHandle);
    //    cleanPfMet - built from all candidates except the ones associated to secondary vertices
    T_MET_cleanPfMet_ET = cleanPfMetHandle->pt();
    T_MET_cleanPfMet_Phi = cleanPfMetHandle->phi();
  */
  //All details can be found in the draft of CMS-AN 2011/374
 
  if(!IsRealData){
    edm::Handle<reco::GenMETCollection> metGen;
    iEvent.getByLabel(genmetLabel_,metGen);
    const GenMET * metsGen= &((metGen.product())->front());

    T_METgen_ET = metsGen[0].pt();
    T_METgen_Phi = metsGen[0].phi();
  }
 
  /*   Handle<TriggerResults> trh;
       try {iEvent.getByLabel(trigLabel_,trh);
       unsigned int aux = trh.product()->size();
       aux = 0 + aux;
       } catch(...) {;}
  */
  //HLT bits and Trigger Matching
  T_passTriggerDoubleMu=false;
  T_passTriggerDoubleEl=false;
  //T_passTriggerSingleMu=false;
  //T_passTriggerSingleEl=false;
  T_passTriggerElMu=false;

  if(IsRealData){
    //    T_passTriggerSingleMu=singleMuData_.check(iEvent,*trh);
    //    T_passTriggerSingleEl=singleElData_.check(iEvent,*trh);
    //    passBits.push_back( singleElData_.check(iEvent,*trh);
    T_passTriggerDoubleMu= doubleMuData_.check(iEvent,*trh);
    T_passTriggerDoubleEl= doubleElData_.check(iEvent,*trh );
    T_passTriggerElMu = muEGData_.check(    iEvent,*trh) ;
  }
  else{
    //    T_passTriggerSingleMu= singleMuMC_.check(iEvent,*trh) ;
    //    T_passTriggerSingleEl= singleElMC_.check(iEvent,*trh) ;
    //singleElMC_.check(iEvent,*trh) );
    T_passTriggerDoubleMu= doubleMuMC_.check(iEvent,*trh) ;
    T_passTriggerDoubleEl= doubleElMC_.check(iEvent,*trh );
    T_passTriggerElMu = muEGMC_.check(    iEvent,*trh);
  }
  /*
    for (unsigned int i=0;i<leps_.size();++i){

    if(!T_passTriggerDoubleMu) T_passTriggerDoubleMu= passTriggerDoubleMu(i, IsRealData);
    if(!T_passTriggerDoubleEl) T_passTriggerDoubleEl=passTriggerDoubleEl(i,IsRealData);
    if(!T_passTriggerSingleMu) T_passTriggerSingleMu=passTriggerSingleMu(i,IsRealData);
    if(!T_passTriggerElMu)  T_passTriggerElMu=passTriggerElMu(i,IsRealData);
    }

  */
  Tree->Fill();
  //  leps_.clear();
  //Gen

  delete T_Gen_StopMass;
  delete T_Gen_Chi0Mass;
  delete T_Gen_CharginoMass;

  delete T_Gen_polWeights;

  /*
  delete T_Gen_StopSt3 ;
  delete T_Gen_Chi0St3 ;
  delete T_Gen_tSt3    ;
  delete T_Gen_ChiPMSt3;
  delete T_Gen_bSt3    ;
  delete T_Gen_WSt3    ;
  delete T_Gen_MuonSt3 ;
  delete T_Gen_ElecSt3 ;
  delete T_Gen_TauSt3  ;
  */

  delete T_Gen_StopSt3_pdgId;	   
  delete T_Gen_StopSt3_firstMother;
  delete T_Gen_StopSt3_i;
  delete T_Gen_StopSt3_energy;	   
  delete T_Gen_StopSt3_pt;	   
  delete T_Gen_StopSt3_eta;	   
  delete T_Gen_StopSt3_phi;        
  
  delete T_Gen_Chi0St3_pdgId;	   
  delete T_Gen_Chi0St3_firstMother;
  delete T_Gen_Chi0St3_i;
  delete T_Gen_Chi0St3_energy;	   
  delete T_Gen_Chi0St3_pt;	   
  delete T_Gen_Chi0St3_eta;	   
  delete T_Gen_Chi0St3_phi;        
  
  delete T_Gen_tSt3_pdgId;	   
  delete T_Gen_tSt3_firstMother;
  delete T_Gen_tSt3_i;
  delete T_Gen_tSt3_energy;	   
  delete T_Gen_tSt3_pt;	   	
  delete T_Gen_tSt3_eta;	   
  delete T_Gen_tSt3_phi;        

  delete T_Gen_ChiPMSt3_pdgId;	   
  delete T_Gen_ChiPMSt3_firstMother;
  delete T_Gen_ChiPMSt3_i;
  delete T_Gen_ChiPMSt3_energy;	   
  delete T_Gen_ChiPMSt3_pt;	   
  delete T_Gen_ChiPMSt3_eta;	   
  delete T_Gen_ChiPMSt3_phi;        

  delete T_Gen_bSt3_pdgId;	
  delete T_Gen_bSt3_firstMother;
  delete T_Gen_bSt3_i;
  delete T_Gen_bSt3_energy;	
  delete T_Gen_bSt3_pt;	   	
  delete T_Gen_bSt3_eta;	
  delete T_Gen_bSt3_phi;        

  delete T_Gen_WSt3_pdgId;	
  delete T_Gen_WSt3_firstMother;
  delete T_Gen_WSt3_i;
  delete T_Gen_WSt3_energy;	
  delete T_Gen_WSt3_pt;   	
  delete T_Gen_WSt3_eta;	
  delete T_Gen_WSt3_phi;        

  delete T_Gen_MuonSt3_pdgId;	
  delete T_Gen_MuonSt3_firstMother;  
  delete T_Gen_MuonSt3_i;  
  delete T_Gen_MuonSt3_energy;	
  delete T_Gen_MuonSt3_pt;   	
  delete T_Gen_MuonSt3_eta;	
  delete T_Gen_MuonSt3_phi;        

  delete T_Gen_ElecSt3_pdgId;	
  delete T_Gen_ElecSt3_firstMother;
  delete T_Gen_ElecSt3_i;
  delete T_Gen_ElecSt3_energy;	
  delete T_Gen_ElecSt3_pt;   	
  delete T_Gen_ElecSt3_eta;	
  delete T_Gen_ElecSt3_phi;        
	 
  delete T_Gen_TauSt3_pdgId;	  
  delete T_Gen_TauSt3_firstMother;
  delete T_Gen_TauSt3_i;
  delete T_Gen_TauSt3_energy;	  
  delete T_Gen_TauSt3_pt;	  
  delete T_Gen_TauSt3_eta;	  
  delete T_Gen_TauSt3_phi;        

  /*
  delete T_Gen_MuonSt3_PID;
  delete T_Gen_MuonSt3_Px;
  delete T_Gen_MuonSt3_Py;
  delete T_Gen_MuonSt3_Pz;
  delete T_Gen_MuonSt3_Energy;
  
  delete T_Gen_ElecSt3_PID;
  delete T_Gen_ElecSt3_Px;
  delete T_Gen_ElecSt3_Py;
  delete T_Gen_ElecSt3_Pz;
  delete T_Gen_ElecSt3_Energy;

  delete T_Gen_bSt3_PID;
  delete T_Gen_bSt3_Px;
  delete T_Gen_bSt3_Py;
  delete T_Gen_bSt3_Pz;
  delete T_Gen_bSt3_Energy;

  delete T_Gen_tSt3_PID;
  delete T_Gen_tSt3_Px;
  delete T_Gen_tSt3_Py;
  delete T_Gen_tSt3_Pz;
  delete T_Gen_tSt3_Energy;
  */
  

  delete T_Gen_Muon_PID;
  delete T_Gen_Muon_Px;
  delete T_Gen_Muon_Py;
  delete T_Gen_Muon_Pz;
  delete T_Gen_Muon_Energy;

  delete T_Gen_Elec_PID;
  delete T_Gen_Elec_Px;
  delete T_Gen_Elec_Py;
  delete T_Gen_Elec_Pz;
  delete T_Gen_Elec_Energy;

  delete T_Gen_b_PID;
  delete T_Gen_b_Px;
  delete T_Gen_b_Py;
  delete T_Gen_b_Pz;
  delete T_Gen_b_Energy;

  delete T_Gen_Muon_MPID;
  delete T_Gen_Muon_MPx;
  delete T_Gen_Muon_MPy;
  delete T_Gen_Muon_MPz;
  delete T_Gen_Muon_MEnergy;
  delete T_Gen_Muon_MSt;

  delete T_Gen_Elec_MPID;
  delete T_Gen_Elec_MPx;
  delete T_Gen_Elec_MPy;
  delete T_Gen_Elec_MPz;
  delete T_Gen_Elec_MEnergy;
  delete T_Gen_Elec_MSt;

  delete T_Gen_b_MPID;
  delete T_Gen_b_MPx;
  delete T_Gen_b_MPy;
  delete T_Gen_b_MPz;
  delete T_Gen_b_MEnergy;
  delete T_Gen_b_MSt;

  /*
  delete T_Gen_TauSt3_PID;
  delete T_Gen_TauSt3_Px;
  delete T_Gen_TauSt3_Py;
  delete T_Gen_TauSt3_Pz;
  delete T_Gen_TauSt3_Energy;
  */

  delete T_Gen_TauSt3_IsLepDec;
  delete T_Gen_TauSt3_LepDec_PID;
  delete T_Gen_TauSt3_LepDec_Px;
  delete T_Gen_TauSt3_LepDec_Py;
  delete T_Gen_TauSt3_LepDec_Pz;
  delete T_Gen_TauSt3_LepDec_Energy;
  /*
  //PFcand
  delete T_PFreducedCand_Px;
  delete T_PFreducedCand_Py;
  delete T_PFreducedCand_Pz;
  delete T_PFreducedCand_E;
  delete T_PFreducedCand_ID;
  delete T_PFreducedCand_vz;
  */     
  //Vertex 
  delete T_Vertex_z;
  delete T_Vertex_y;
  delete T_Vertex_x;
  delete T_Vertex_Chi2Prob;
  delete T_Vertex_rho;
  delete T_Vertex_ndof;
  delete T_Vertex_isFake; 
  delete T_Vertex_tracksSize;
  //Muons
  delete T_Muon_Eta;
  delete T_Muon_IsGlobalMuon;
  //  delete T_Muon_IsTMLSOLPT;
  //  delete T_Muon_IsTMLSOLPL;
  delete T_Muon_IsGMPTMuons;
  //delete T_Muon_IsAllStandAloneMuons;
  delete T_Muon_IsAllTrackerMuons;
  delete T_Muon_IsTrackerMuonArbitrated;
  delete T_Muon_IsAllArbitrated;
  /*delete T_Muon_IsTMLastStationLoose;
    delete T_Muon_IsTMLastStationTight;
    delete T_Muon_IsTM2DCompatibilityLoose;
    delete T_Muon_IsTM2DCompatibilityTight;
    delete T_Muon_IsTMOneStationLoose;
    delete T_Muon_IsTMOneStationTight;
    delete T_Muon_IsTMLSOPL;
    delete T_Muon_IsGMTkChiCompatibility;
    delete T_Muon_IsGMStaChiCompatibility;
    delete T_Muon_IsGMTkKinkTight;*/
  delete T_Muon_Px;
  delete T_Muon_Py;
  delete T_Muon_Pz;
  delete T_Muon_Pt;
  delete T_Muon_deltaPt;
  delete T_Muon_Energy;
  delete T_Muon_Charge;
  delete T_Muon_NormChi2GTrk;
  delete T_Muon_NValidHitsInTrk;
  delete T_Muon_NValidHitsSATrk;
  delete T_Muon_NValidHitsGTrk;
  delete T_Muon_NValidPixelHitsInTrk;
  delete T_Muon_Chi2InTrk;
  delete T_Muon_dofInTrk;
  delete T_Muon_IPAbsGTrack;
  //  delete T_Muon_IPSigGTrack;
  delete T_Muon_IPAbsInTrack;
  //  delete T_Muon_IPSigInTrack;
  //  delete T_Muon_IPwrtBSInTrack;
  delete T_Muon_IPwrtAveBSInTrack;
  //  delete T_Muon_IPwrtBSGTrack;
  /*  delete T_Muon_IP2DBiasedPV;
      delete T_Muon_IP3DBiasedPV;
      delete T_Muon_IP2DUnBiasedPV;
      delete T_Muon_IP3DUnBiasedPV;
      delete T_Muon_dxyPVBiasedPV;
      delete T_Muon_dxyPVUnBiasedPV;
      delete T_Muon_dzPVBiasedPV;
      delete T_Muon_dzPVUnBiasedPV;
      delete T_Muon_IP2DUnBiasedPVnoBS;
      delete T_Muon_IP3DUnBiasedPVnoBS;
      delete T_Muon_dxyPVUnBiasedPVnoBS;
      delete T_Muon_dzPVUnBiasedPVnoBS;*/
  /*  delete T_Muon_smurfCharged;
      delete T_Muon_smurfNeutral;
      delete T_Muon_smurfPhoton;
      delete T_Muon_smurfNoOverCharged;
      delete T_Muon_smurfNoOverNeutral;
      delete T_Muon_smurfNoOverPhoton;
      delete T_Muon_muSmurfPF;*/
  delete T_Muon_NumOfMatchedStations;
  delete T_Muon_InnerTrackFound;
  delete T_Muon_vz;
  delete T_Muon_vy;
  delete T_Muon_vx;
  /*  delete T_Muon_IsTMLastStationAngLoose;
      delete T_Muon_IsTMLastStationAngTight;
      delete T_Muon_IsTMOneStationAngLoose;
      delete T_Muon_IsTMOneStationAngTight;*/
  delete T_Muon_trkKink ;
  delete T_Muon_SegmentCompatibility;
  delete T_Muon_chargedParticleIsoR03;
  delete T_Muon_chargedHadronIsoR03;
  delete T_Muon_neutralHadronIsoR03;
  delete T_Muon_photonIsoR03;
  //  delete T_Muon_chargedParticleIsoR04;
  delete T_Muon_chargedHadronIsoR04;
  delete T_Muon_neutralHadronIsoR04;
  delete T_Muon_photonIsoR04;
  delete T_Muon_sumPUPtR04;
  delete T_Muon_sumPUPtR03;
  delete T_Muon_PFMuonPt;
  delete T_Muon_PFMuonPx;
  delete T_Muon_PFMuonPy;
  delete T_Muon_PFMuonPz;
  delete T_Muon_PFMuonE;
  delete T_Muon_isPFMuon;
  delete T_Muon_NLayers;
  /*  delete T_Muon_passTriggerDoubleMu;
      delete T_Muon_passTriggerSingleMu;
      delete T_Muon_passTriggerElMu;
  */

  //Electrons
  delete T_Elec_Eta;
  delete T_Elec_IPwrtAveBS;
  delete T_Elec_IPwrtPV;
  delete T_Elec_dzwrtPV;
  delete T_Elec_Px;
  delete T_Elec_Py;
  delete T_Elec_Pz;
  delete T_Elec_Pt;
  delete T_Elec_Energy;
  delete T_Elec_Charge;
  /*  delete T_Elec_SumIsoCalo;
      delete T_Elec_SumIsoTrack;
      delete T_Elec_IP2DBiasedPV;
      delete T_Elec_IP3DBiasedPV;
      delete T_Elec_IP2DUnBiasedPV;
      delete T_Elec_IP3DUnBiasedPV;
      delete T_Elec_dxyPVBiasedPV;
      delete T_Elec_dzPVBiasedPV;
      delete T_Elec_dxyPVUnBiasedPV;
      delete T_Elec_dzPVUnBiasedPV;
      delete T_Elec_IP2DUnBiasedPVnoBS;
      delete T_Elec_IP3DUnBiasedPVnoBS;
      delete T_Elec_dxyPVUnBiasedPVnoBS;
      delete T_Elec_dzPVUnBiasedPVnoBS;
  */  delete T_Elec_vz;
  delete T_Elec_vy;
  delete T_Elec_vx;
  delete T_Elec_nLost;
  delete T_Elec_nHits;
  delete T_Elec_SC_Et;
  delete T_Elec_SC_Eta;
  delete T_Elec_nBrems;
  delete T_Elec_fBrem;
  delete T_Elec_eSuperClusterOverP;
  delete T_Elec_ecalEnergy;
  /*  delete T_Elec_dr03TkSumPt;
      delete T_Elec_dr03EcalSumEt;
      delete T_Elec_dr03HcalSumEt;
      delete T_Elec_ConvInfoDist;
      delete T_Elec_ConvInfoDCot;
  */  delete T_Elec_isEB;
  delete T_Elec_isEE;
  delete T_Elec_MVA;
  //  delete T_Elec_simpleEleId95;
  //  delete T_Elec_simpleEleId90;
  //  delete T_Elec_simpleEleId90cIso;
  //  delete T_Elec_simpleEleId85;
  delete T_Elec_simpleEleId80;
  //  delete T_Elec_simpleEleId70;
  //  delete T_Elec_simpleEleId60;
  /*delete T_Elec_cicVeryLooseMC;
    delete T_Elec_cicLooseMC;
    delete T_Elec_cicMediumMC;
    delete T_Elec_cicSuperTightMC;
    delete T_Elec_cicHyperTight1MC;
    delete T_Elec_cicHyperTight2MC;
    delete T_Elec_cicHyperTight3MC;
    delete T_Elec_cicVeryLooseHWW;
    delete T_Elec_cicLooseHWW;
    delete T_Elec_cicMediumHWW;
    delete T_Elec_cicSuperTightHWW;
    delete T_Elec_cicHyperTight1HWW;
    delete T_Elec_cicHyperTight2HWW;
    delete T_Elec_cicHyperTight3HWW;
    delete T_Elec_cicVeryLoose;
    delete T_Elec_cicLoose;
    delete T_Elec_cicMedium;
    delete T_Elec_cicSuperTight;
    delete T_Elec_cicHyperTight1;
    delete T_Elec_cicHyperTight2;
    delete T_Elec_cicHyperTight3;*/
  //delete T_Elec_egammaIDLikelihood;  
  delete T_Elec_sigmaIetaIeta;
  delete T_Elec_deltaPhiIn;
  delete T_Elec_deltaEtaIn;
  delete T_Elec_isEcalDriven;
  delete T_Elec_HtoE;
  /*  delete T_Elec_passTriggerElMu;
      delete T_Elec_passTriggerDoubleEl;*/
  delete T_Elec_chargedHadronIso;
  delete T_Elec_neutralHadronIso;
  delete T_Elec_photonIso;
  delete T_Elec_puChargedHadronIso;
  /*   delete T_Elec_smurfCharged;
       delete T_Elec_smurfNeutral;
       delete T_Elec_smurfPhoton;
       delete T_Elec_smurfNoOverCharged;
       delete T_Elec_smurfNoOverNeutral;
       delete T_Elec_smurfNoOverPhoton;
       delete T_Elec_eleSmurfPF;*/
  delete T_Elec_passConversionVeto;
  delete T_Elec_PFElecPt;
  delete T_Elec_PFElecPx;
  delete T_Elec_PFElecPy;
  delete T_Elec_PFElecPz;
  delete T_Elec_PFElecE;
  delete T_Elec_isPF;


  //delete T_passTriggerDoubleMu;
  //delete T_passTriggerDoubleEl;
  //delete T_passTriggerSingleMu;
  //delete T_passTriggerSingleEl;
  //delete T_passTriggerElMu;

  /*
    delete T_pfTau_Px;
    delete T_pfTau_Py;
    delete T_pfTau_Pz;
    delete T_pfTau_Energy ;
    delete T_pfTau_Charge ;
  */
  /*  delete T_Tau_Px;
      delete T_Tau_Py;
      delete T_Tau_Pz;
      delete T_Tau_Energy ;
      delete T_Tau_Charge ;
  */

  /*
    delete T_pfMuon_IsGlobalMuon;
    delete T_pfMuon_IsGMPTMuons ;
    delete T_pfMuon_IsAllTrackerMuons ;
    delete T_pfMuon_IsTMLastStationAngTight;
    delete T_pfMuon_IsAllStandAloneMuons;
    delete T_pfMuon_SegmentCompatibility;
    delete T_pfMuon_Px;
    delete T_pfMuon_Py;
    delete T_pfMuon_Pz;
    delete T_pfMuon_Pt;
    delete T_pfMuon_Energy ;
    delete T_pfMuon_Charge ;
    delete T_pfMuon_NormChi2GTrk ;
    delete T_pfMuon_NValidHitsInTrk ;
    delete T_pfMuon_NValidHitsSATrk ;
    delete T_pfMuon_Chi2InTrk ;
    delete T_pfMuon_dofInTrk ;
    delete T_pfMuon_SumIsoCalo ;
    delete T_pfMuon_SumIsoTrack ;
    delete T_pfMuon_IPAbsGTrack ;
    delete T_pfMuon_IPSigGTrack ;
    delete T_pfMuon_IPAbsInTrack ;
    delete T_pfMuon_IPSigInTrack ;
    //  delete T_pfMuon_IPwrtBSInTrack ;
    delete T_pfMuon_IP2DBiasedPV;
    delete T_pfMuon_IP3DBiasedPV;
    delete T_pfMuon_IP2DUnBiasedPV;
    delete T_pfMuon_IP3DUnBiasedPV;
    delete T_pfMuon_vz ;
    delete T_pfMuon_vy ;
    delete T_pfMuon_vx ;
    delete T_pfMuon_NValidHits ;
    delete T_pfMuon_NValidPixelHitsInTrk;
    delete T_pfMuon_particleIso ;
    delete T_pfMuon_chargedHadronIso ;
    delete T_pfMuon_neutralHadronIso ;
    delete T_pfMuon_photonIso ;
    delete T_pfMuon_deltaPt;
    delete T_pfMuon_NumOfMatchedStations;
    delete T_pfMuon_pfCharged;
    delete T_pfMuon_pfNeutral;
    delete T_pfMuon_pfPhoton;
    delete T_pfMuon_smurfCharged;
    delete T_pfMuon_smurfNeutral;
    delete T_pfMuon_smurfPhoton;
    delete T_pfMuon_smurfNoOverCharged;
    delete T_pfMuon_smurfNoOverNeutral;
    delete T_pfMuon_smurfNoOverPhoton;
    delete T_pfMuon_muSmurfPF;
    delete T_pfMuon_sumPUPt;
    delete T_pfMuon_sumPUPtR03;
    delete T_pfMuon_NLayers;
  */
  /*
    delete T_pfElec_Px;
    delete T_pfElec_Py;
    delete T_pfElec_Pz;
    delete T_pfElec_Pt;
    delete T_pfElec_Energy;
    delete T_pfElec_Charge;
    delete T_pfElec_SumIsoCalo;
    delete T_pfElec_SumIsoTrack;
    delete T_pfElec_IP2DBiasedPV;
    delete T_pfElec_IP3DBiasedPV;
    delete T_pfElec_IP2DUnBiasedPV;
    delete T_pfElec_IP3DUnBiasedPV;
    delete T_pfElec_vz;
    delete T_pfElec_vy;
    delete T_pfElec_vx;
    delete T_pfElec_nBrems;
    delete T_pfElec_dr03TkSumPt;
    delete T_pfElec_dr03EcalSumEt;
    delete T_pfElec_dr03HcalSumEt;
    delete T_pfElec_SC_Et;
    delete T_pfElec_SC_Eta;
    delete T_pfElec_nHits;
    delete T_pfElec_ConvInfoDCot;
    delete T_pfElec_ConvInfoDist;
    delete T_pfElec_isEcalDriven;
    delete T_pfElec_HtoE;
    delete T_pfElec_photonIso;
    delete T_pfElec_chargedHadronIso;
    delete T_pfElec_neutralHadronIso;
    delete T_pfElec_particleIso;
  */
  //Jets
  for (int i = 0; i < NumberOfJetCollections; i++) {
    delete T_Jet_Px[i];
    delete T_Jet_Py[i];
    delete T_Jet_Pz[i];
    delete T_Jet_Et[i];
    /*    delete T_Jet_Px11[i];
	  delete T_Jet_Py11[i];
	  delete T_Jet_Pz11[i];
	  delete T_Jet_Et11[i];*/
    //    delete T_Jet_EtOffset[i];
    //    delete T_Jet_PtOffset[i];
    delete T_Jet_Eta[i];
    delete T_Jet_Energy[i];
    //    delete T_Jet_Corr[i];
    //    delete T_Jet_Corr11[i];
    delete T_Jet_Tag_HighEffTC[i];  
    delete T_Jet_Tag_CombSVtx[i];
    delete T_Jet_Tag_CombSVtxMVA[i];
    delete T_Jet_Tag_TauJet[i];
    delete T_Jet_Tag_ImpParMVA[i];
    delete T_Jet_Tag_JetBProb[i];
    delete T_Jet_Tag_JetProb[i];
    delete T_Jet_Tag_HighEffSimpSVtx[i];
    delete T_Jet_Tag_HighPurSimpSVtx[i];
    delete T_Jet_Tag_HighPurTC[i];
    delete T_Jet_Parton_Px[i];
    delete T_Jet_Parton_Py[i];
    delete T_Jet_Parton_Pz[i];
    delete T_Jet_Parton_Energy[i];
    delete T_Jet_Parton_Flavour[i];
    delete T_Jet_CharHadEnergyFrac[i];
    delete T_Jet_NeutHadEnergyFrac[i];
    delete T_Jet_CharEmEnergyFrac[i];
    delete T_Jet_NeutEmEnergyFrac[i];
    delete T_Jet_CharHadEnergy[i];
    delete T_Jet_NeutHadEnergy[i];
    delete T_Jet_CharEmEnergy[i];
    delete T_Jet_NeutEmEnergy[i];
    delete T_Jet_MuonMultiplicity[i];
    delete T_Jet_NeutralMultiplicity[i];
    delete T_Jet_ChargedMultiplicity[i];
    delete T_Jet_IDLoose[i];
    delete T_Jet_nDaughters[i];

    //	if(!IsRealData){
    delete T_Jet_GenJet_InvisibleE[i];
    delete T_Jet_GenJet_Px[i];
    delete T_Jet_GenJet_Py[i];
    delete T_Jet_GenJet_Pz[i];
    delete T_Jet_GenJet_Et[i];
    delete T_Jet_GenJet_Eta[i];
    delete T_Jet_GenJet_Energy[i];
    delete T_Jet_IsGenJet[i];
    //    }
  }
  
  
}


void SUSYSkimToTree::SetJetInfo(int idx, edm::View<pat::Jet> JET, const reco::VertexCollection& vtxs, bool calojet) {
  T_Jet_Px[idx] = new std::vector<float>;
  T_Jet_Py[idx] = new std::vector<float>;
  T_Jet_Pz[idx] = new std::vector<float>;
  T_Jet_Et[idx] = new std::vector<float>;
  /*  T_Jet_Px11[idx] = new std::vector<float>;
      T_Jet_Py11[idx] = new std::vector<float>;
      T_Jet_Pz11[idx] = new std::vector<float>;
      T_Jet_Et11[idx] = new std::vector<float>;*/
  //  T_Jet_EtOffset[idx] = new std::vector<float>;
  //  T_Jet_PtOffset[idx] = new std::vector<float>;
  T_Jet_Eta[idx] = new std::vector<float>;
  T_Jet_Energy[idx] = new std::vector<float>;
  //  T_Jet_Corr[idx] = new std::vector<float>;
  //  T_Jet_Corr11[idx] = new std::vector<float>;
  T_Jet_Parton_Px[idx] = new std::vector<float>;
  T_Jet_Parton_Py[idx] = new std::vector<float>;
  T_Jet_Parton_Pz[idx] = new std::vector<float>;
  T_Jet_Parton_Energy[idx] = new std::vector<float>;
  T_Jet_Parton_Flavour[idx] = new std::vector<int>;
  T_Jet_Tag_HighEffTC[idx] = new std::vector<float>;
  T_Jet_Tag_CombSVtx[idx] = new std::vector<float>;
  T_Jet_Tag_CombSVtxMVA[idx] = new std::vector<float>;
  T_Jet_Tag_TauJet[idx] = new std::vector<float>;
  T_Jet_Tag_ImpParMVA[idx] = new std::vector<float>;
  T_Jet_Tag_JetBProb[idx] = new std::vector<float>;
  T_Jet_Tag_JetProb[idx] = new std::vector<float>;
  T_Jet_Tag_HighPurSimpSVtx[idx] = new std::vector<float>;
  T_Jet_Tag_HighEffSimpSVtx[idx] = new std::vector<float>;
  T_Jet_Tag_HighPurTC[idx] = new std::vector<float>;
  
  T_Jet_CharHadEnergyFrac[idx] = new std::vector<float>;
  T_Jet_NeutHadEnergyFrac[idx] = new std::vector<float>;
  T_Jet_CharEmEnergyFrac[idx] = new std::vector<float>;
  T_Jet_NeutEmEnergyFrac[idx] = new std::vector<float>;
  T_Jet_CharHadEnergy[idx] = new std::vector<float>;
  T_Jet_NeutHadEnergy[idx] = new std::vector<float>;
  T_Jet_CharEmEnergy[idx] = new std::vector<float>;
  T_Jet_NeutEmEnergy[idx] = new std::vector<float>;
  T_Jet_MuonMultiplicity[idx] = new std::vector<int>;
  T_Jet_NeutralMultiplicity[idx] = new std::vector<int>; 
  T_Jet_ChargedMultiplicity[idx] = new std::vector<int>;
  T_Jet_IDLoose[idx] = new std::vector<bool>;
  T_Jet_nDaughters[idx] = new std::vector<int>;
  
  T_Jet_GenJet_InvisibleE[idx] = new std::vector<float>;
  T_Jet_GenJet_Px[idx] = new std::vector<float>;
  T_Jet_GenJet_Py[idx] = new std::vector<float>;
  T_Jet_GenJet_Pz[idx] = new std::vector<float>;
  T_Jet_GenJet_Et[idx] = new std::vector<float>;
  T_Jet_GenJet_Eta[idx] = new std::vector<float>;
  T_Jet_GenJet_Energy[idx] = new std::vector<float>;
  T_Jet_IsGenJet[idx] = new std::vector<bool>;
  
  
  PFJetIDSelectionFunctor PFjetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA,
					PFJetIDSelectionFunctor::LOOSE );

  // New 2012 JEC
  // Create the JetCorrectorParameter objects, the order does not matter.
  /*	JetCorrectorParameters *ResJetPar,*L3JetPar,*L2JetPar,*L1JetPar;
    std::vector<JetCorrectorParameters> vPar;
    if (IsRealData){
    ResJetPar = new JetCorrectorParameters("/gpfs/csic_users/jfernan/CMSSW_5_3_3_patch2/src/GTs/GR_P_V42_AN3_L2L3Residual_AK5PFchs.txt"); 
    L3JetPar  = new JetCorrectorParameters("/gpfs/csic_users/jfernan/CMSSW_5_3_3_patch2/src/GTs/GR_P_V42_AN3_L3Absolute_AK5PFchs.txt");
    L2JetPar  = new JetCorrectorParameters("/gpfs/csic_users/jfernan/CMSSW_5_3_3_patch2/src/GTs/GR_P_V42_AN3_L2Relative_AK5PFchs.txt");
    L1JetPar  = new JetCorrectorParameters("/gpfs/csic_users/jfernan/CMSSW_5_3_3_patch2/src/GTs/GR_P_V42_AN3_L1FastJet_AK5PFchs.txt");
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);
    vPar.push_back(*ResJetPar);
    }
    else{
    L3JetPar  = new JetCorrectorParameters("/gpfs/csic_users/jfernan/CMSSW_5_3_3_patch2/src/GTs/START53_V15_L3Absolute_AK5PFchs.txt");
    L2JetPar  = new JetCorrectorParameters("/gpfs/csic_users/jfernan/CMSSW_5_3_3_patch2/src/GTs/START53_V15_L2Relative_AK5PFchs.txt");
    L1JetPar  = new JetCorrectorParameters("/gpfs/csic_users/jfernan/CMSSW_5_3_3_patch2/src/GTs/START53_V15_L1FastJet_AK5PFchs.txt");
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);
    }

    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);

    float rho = T_Event_Rho;
  */ 
  for (edm::View<pat::Jet>::const_iterator jet_iter=JET.begin(); jet_iter!= JET.end(); jet_iter++) { 
    /*   double correction(1.0);
	 Jet mijetRAW = jet_iter->correctedJet("Uncorrected");
	 JetCorrector->setJetEta(jet_iter->eta());
	 JetCorrector->setJetPt(mijetRAW.pt());
	 JetCorrector->setJetA(jet_iter->jetArea());
	 JetCorrector->setRho(rho);
	 //   cout <<"et: "<<mijetRAW.et()<<" pt: "<<mijetRAW.pt()<<endl;

	 correction = JetCorrector->getCorrection();

	 if((mijetRAW.et()*correction)<10) continue;
	 //	if(jet_iter->et()<10) continue;	
	 */
    T_Jet_Px[idx]->push_back(jet_iter->px());
    T_Jet_Py[idx]->push_back(jet_iter->py());
    T_Jet_Pz[idx]->push_back(jet_iter->pz());
    T_Jet_Et[idx]->push_back(jet_iter->et());
    //    T_Jet_Corr[idx]->push_back(correction);

    T_Jet_Eta[idx]->push_back(jet_iter->eta());
    T_Jet_Energy[idx]->push_back(jet_iter->energy());

    bool jetIdLoose_ = false;

    if (jet_iter->isPFJet()) {
      T_Jet_CharHadEnergyFrac[idx]->push_back(jet_iter->chargedHadronEnergyFraction());
      T_Jet_NeutHadEnergyFrac[idx]->push_back(jet_iter->neutralHadronEnergyFraction()); 
      T_Jet_CharEmEnergyFrac[idx]->push_back(jet_iter->chargedEmEnergyFraction()); 
      T_Jet_NeutEmEnergyFrac[idx]->push_back(jet_iter->neutralEmEnergyFraction());
      T_Jet_CharHadEnergy[idx]->push_back(jet_iter->chargedHadronEnergy());
      T_Jet_NeutHadEnergy[idx]->push_back(jet_iter->neutralHadronEnergy()); 
      T_Jet_CharEmEnergy[idx]->push_back(jet_iter->chargedEmEnergy()); 
      T_Jet_NeutEmEnergy[idx]->push_back(jet_iter->neutralEmEnergy());
      T_Jet_MuonMultiplicity[idx]->push_back(jet_iter->muonMultiplicity());
      T_Jet_NeutralMultiplicity[idx]->push_back(jet_iter->neutralMultiplicity());
      T_Jet_ChargedMultiplicity[idx]->push_back(jet_iter->chargedMultiplicity());
      jetIdLoose_ = PFjetIDLoose(*jet_iter);

      T_Jet_IDLoose[idx]->push_back(jetIdLoose_); 
      T_Jet_nDaughters[idx]->push_back(jet_iter->numberOfDaughters());
    }
    
    /*    try {
	  Jet mijetraw = jet_iter->correctedJet("Uncorrected");
	  float corrJES = jet_iter->et()/mijetraw.et();
	  T_Jet_Corr11[idx]->push_back(corrJES);      
	  } catch(...) {
	  T_Jet_Corr11[idx]->push_back(-9999);
	  }   
    */


    //OffsetJet JES
    /*
      if(idx==0){
      try{
      T_Jet_EtOffset[idx]->push_back(jet_iter->correctedJet("L3Absolute","none","patJetCorrFactors").et());
      T_Jet_PtOffset[idx]->push_back(jet_iter->correctedJet("L3Absolute","none","patJetCorrFactors").pt());
      }catch(...){
      T_Jet_EtOffset[idx]->push_back(-9999);
      T_Jet_PtOffset[idx]->push_back(-9999);
      }}*/
    /*if(idx==0){
      try{
      T_Jet_PtOffset[idx]->push_back(jet_iter->correctedJet("L3Absolute","none","patJetCorrFactorsNoPU").pt());
      T_Jet_EtOffset[idx]->push_back(jet_iter->correctedJet("L3Absolute","none","patJetCorrFactorsNoPU").et());
      }catch(...){
      T_Jet_EtOffset[idx]->push_back(-9999);
      T_Jet_PtOffset[idx]->push_back(-9999);
      }}
    */
    /*if(idx==2){
      try{
      T_Jet_PtOffset[idx]->push_back(jet_iter->correctedJet("L3Absolute","none","patJetCorrFactorsPFlow").pt());
      T_Jet_EtOffset[idx]->push_back(jet_iter->correctedJet("L3Absolute","none","patJetCorrFactorsPFlow").et());
      }catch(...){
      T_Jet_EtOffset[idx]->push_back(-9999);
      T_Jet_PtOffset[idx]->push_back(-9999);
      }}
    */
    //
    //std::cout<<"partonFlavour"<<endl;
    T_Jet_Parton_Flavour[idx]->push_back(jet_iter->partonFlavour());
    // std::cout<<"genParton"<<endl;
    if (jet_iter->genParton()) {
      T_Jet_Parton_Px[idx]->push_back(jet_iter->genParton()->px());
      T_Jet_Parton_Py[idx]->push_back(jet_iter->genParton()->py());
      T_Jet_Parton_Pz[idx]->push_back(jet_iter->genParton()->pz());
      T_Jet_Parton_Energy[idx]->push_back(jet_iter->genParton()->energy());
    }
    else {
      T_Jet_Parton_Px[idx]->push_back(0);
      T_Jet_Parton_Py[idx]->push_back(0);
      T_Jet_Parton_Pz[idx]->push_back(0);
      T_Jet_Parton_Energy[idx]->push_back(0);
    }
   
    if(!IsRealData){
      //   std::cout<<"GenJet"<<endl;
      const reco::GenJet * mygenJet=jet_iter->genJet();
    
      if (mygenJet != 0) {
	T_Jet_IsGenJet[idx]->push_back(true);
	T_Jet_GenJet_Px[idx]->push_back(mygenJet->px());
	T_Jet_GenJet_Py[idx]->push_back(mygenJet->py());
	T_Jet_GenJet_Pz[idx]->push_back(mygenJet->pz());
	T_Jet_GenJet_Et[idx]->push_back(mygenJet->et());
	T_Jet_GenJet_Eta[idx]->push_back(mygenJet->eta());
	T_Jet_GenJet_Energy[idx]->push_back(mygenJet->energy());
	T_Jet_GenJet_InvisibleE[idx]->push_back(mygenJet->invisibleEnergy());
      }
      else {
	T_Jet_IsGenJet[idx]->push_back(false);
	T_Jet_GenJet_Px[idx]->push_back(0);
	T_Jet_GenJet_Py[idx]->push_back(0);
	T_Jet_GenJet_Pz[idx]->push_back(0);
	T_Jet_GenJet_Et[idx]->push_back(0);
	T_Jet_GenJet_Eta[idx]->push_back(0);
	T_Jet_GenJet_Energy[idx]->push_back(0);
	T_Jet_GenJet_InvisibleE[idx]->push_back(0);
      }
    }
   
    
    T_Jet_Tag_HighEffTC[idx]->push_back(jet_iter->bDiscriminator("trackCountingHighEffBJetTags"));
    T_Jet_Tag_CombSVtx[idx]->push_back(jet_iter->bDiscriminator("combinedSecondaryVertexBJetTags"));
    T_Jet_Tag_CombSVtxMVA[idx]->push_back(jet_iter->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
    T_Jet_Tag_TauJet[idx]->push_back(jet_iter->bDiscriminator("coneIsolationTauJetTags"));
    T_Jet_Tag_ImpParMVA[idx]->push_back(jet_iter->bDiscriminator("impactParameterMVABJetTags"));
    T_Jet_Tag_JetBProb[idx]->push_back(jet_iter->bDiscriminator("jetBProbabilityBJetTags"));
    T_Jet_Tag_JetProb[idx]->push_back(jet_iter->bDiscriminator("jetProbabilityBJetTags"));
    T_Jet_Tag_HighEffSimpSVtx[idx]->push_back(jet_iter->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
    T_Jet_Tag_HighPurSimpSVtx[idx]->push_back(jet_iter->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
    T_Jet_Tag_HighPurTC[idx]->push_back(jet_iter->bDiscriminator("trackCountingHighPurBJetTags"));
  }

  /* delete JetCorrector;
     delete L3JetPar;
     delete L2JetPar;
     delete L1JetPar;
     if(IsRealData) delete ResJetPar;*/
}

void SUSYSkimToTree::SetJetBranchAddress(int idx, TString namecol, bool calojet) {
  
  Tree->Branch(TString(namecol + "_Px") , "std::vector<float>", &T_Jet_Px[idx]);
  Tree->Branch(TString(namecol + "_Py"), "std::vector<float>", &T_Jet_Py[idx]);
  Tree->Branch(TString(namecol + "_Pz"), "std::vector<float>", &T_Jet_Pz[idx]);
  Tree->Branch(TString(namecol + "_Et"), "std::vector<float>", &T_Jet_Et[idx]);
  /*  Tree->Branch(TString(namecol + "_Px11") , "std::vector<float>", &T_Jet_Px11[idx]);
      Tree->Branch(TString(namecol + "_Py11"), "std::vector<float>", &T_Jet_Py11[idx]);
      Tree->Branch(TString(namecol + "_Pz11"), "std::vector<float>", &T_Jet_Pz11[idx]);
      Tree->Branch(TString(namecol + "_Et11"), "std::vector<float>", &T_Jet_Et11[idx]);*/
  //  Tree->Branch(TString(namecol + "_EtOffset"), "std::vector<float>", &T_Jet_EtOffset[idx]);
  //  Tree->Branch(TString(namecol + "_PtOffset"), "std::vector<float>", &T_Jet_PtOffset[idx]);
  Tree->Branch(TString(namecol + "_Eta"), "std::vector<float>", &T_Jet_Eta[idx]);
  Tree->Branch(TString(namecol + "_Energy"), "std::vector<float>", &T_Jet_Energy[idx]);
  //  Tree->Branch(TString(namecol + "_Corr"), "std::vector<float>", &T_Jet_Corr[idx]);
  //  Tree->Branch(TString(namecol + "_Corr11"), "std::vector<float>", &T_Jet_Corr11[idx]);
  Tree->Branch(TString(namecol + "_Tag_HighEffTC"), "std::vector<float>", &T_Jet_Tag_HighEffTC[idx]);
  Tree->Branch(TString(namecol + "_Tag_CombSVtx"), "std::vector<float>", &T_Jet_Tag_CombSVtx[idx]);
  Tree->Branch(TString(namecol + "_Tag_CombSVtxMVA"), "std::vector<float>", &T_Jet_Tag_CombSVtxMVA[idx]);
  Tree->Branch(TString(namecol + "_Tag_TauJet"), "std::vector<float>", &T_Jet_Tag_TauJet[idx]);
  Tree->Branch(TString(namecol + "_Tag_ImpParMVA"), "std::vector<float>", &T_Jet_Tag_ImpParMVA[idx]);
  Tree->Branch(TString(namecol + "_Tag_JetBProb"), "std::vector<float>", &T_Jet_Tag_JetBProb[idx]);
  Tree->Branch(TString(namecol + "_Tag_JetProb"), "std::vector<float>", &T_Jet_Tag_JetProb[idx]);
  Tree->Branch(TString(namecol + "_Tag_HighEffSimpSVtx"), "std::vector<float>", &T_Jet_Tag_HighEffSimpSVtx[idx]);
  Tree->Branch(TString(namecol + "_Tag_HighPurSimpSVtx"), "std::vector<float>", &T_Jet_Tag_HighPurSimpSVtx[idx]);
  Tree->Branch(TString(namecol + "_Tag_HighPurTC"), "std::vector<float>", &T_Jet_Tag_HighPurTC[idx]);
  Tree->Branch(TString(namecol + "_Parton_Px"), "std::vector<float>", &T_Jet_Parton_Px[idx]);
  Tree->Branch(TString(namecol + "_Parton_Py"), "std::vector<float>", &T_Jet_Parton_Py[idx]);
  Tree->Branch(TString(namecol + "_Parton_Pz"), "std::vector<float>", &T_Jet_Parton_Pz[idx]);
  Tree->Branch(TString(namecol + "_Parton_Energy"), "std::vector<float>", &T_Jet_Parton_Energy[idx]);
  Tree->Branch(TString(namecol + "_Parton_Flavour"), "std::vector<int>", &T_Jet_Parton_Flavour[idx]);
  
  Tree->Branch(TString(namecol + "_CharHadEnergyFrac"), "std::vector<float>", &T_Jet_CharHadEnergyFrac[idx]);
  Tree->Branch(TString(namecol + "_NeutHadEnergyFrac"), "std::vector<float>", &T_Jet_NeutHadEnergyFrac[idx]);
  Tree->Branch(TString(namecol + "_CharEmEnergyFrac"), "std::vector<float>", &T_Jet_CharEmEnergyFrac[idx]);
  Tree->Branch(TString(namecol + "_NeutEmEnergyFrac"), "std::vector<float>", &T_Jet_NeutEmEnergyFrac[idx]);
  Tree->Branch(TString(namecol + "_CharHadEnergy"), "std::vector<float>", &T_Jet_CharHadEnergy[idx]);
  Tree->Branch(TString(namecol + "_NeutHadEnergy"), "std::vector<float>", &T_Jet_NeutHadEnergy[idx]);
  Tree->Branch(TString(namecol + "_CharEmEnergy"), "std::vector<float>", &T_Jet_CharEmEnergy[idx]);
  Tree->Branch(TString(namecol + "_NeutEmEnergy"), "std::vector<float>", &T_Jet_NeutEmEnergy[idx]);
  Tree->Branch(TString(namecol + "_MuonMultiplicity"), "std::vector<int>", &T_Jet_MuonMultiplicity[idx]);
  Tree->Branch(TString(namecol + "_NeutralMultiplicity"), "std::vector<int>", &T_Jet_NeutralMultiplicity[idx]);
  Tree->Branch(TString(namecol + "_ChargedMultiplicity"), "std::vector<int>", &T_Jet_ChargedMultiplicity[idx]);
  Tree->Branch(TString(namecol + "_IDLoose"), "std::vector<bool>", &T_Jet_IDLoose[idx]);
  Tree->Branch(TString(namecol + "_nDaughters"), "std::vector<int>", &T_Jet_nDaughters[idx]);
 
 
  //   if(!IsRealData){
  Tree->Branch(TString(namecol + "_GenJet_InvisibleE"), "std::vector<float>", &T_Jet_GenJet_InvisibleE[idx]);
  Tree->Branch(TString(namecol + "_GenJet_Px"), "std::vector<float>", &T_Jet_GenJet_Px[idx]);
  Tree->Branch(TString(namecol + "_GenJet_Py"), "std::vector<float>", &T_Jet_GenJet_Py[idx]);
  Tree->Branch(TString(namecol + "_GenJet_Pz"), "std::vector<float>", &T_Jet_GenJet_Pz[idx]);
  Tree->Branch(TString(namecol + "_GenJet_Et"), "std::vector<float>", &T_Jet_GenJet_Et[idx]);
  Tree->Branch(TString(namecol + "_GenJet_Eta"), "std::vector<float>", &T_Jet_GenJet_Eta[idx]);
  Tree->Branch(TString(namecol + "_GenJet_Energy"), "std::vector<float>", &T_Jet_GenJet_Energy[idx]);
  Tree->Branch(TString(namecol + "_IsGenJet"), "std::vector<bool>", &T_Jet_IsGenJet[idx]);
  //  }
  
  
}


//-------

//It checks if there is a muon or electron in the tau decay chain
//if that is the case it stores the info of the first lepton found in the decay chain
void SUSYSkimToTree::LeptonicTauDecay(const reco::Candidate& tau, bool& elecdec, bool& muondec,
				      int &pid, float &px, float &py, float &pz, float &energy)
  
{
  
  // loop on tau decays, check for an electron or muon
  for(reco::Candidate::const_iterator daughter=tau.begin();daughter!=tau.end(); ++daughter){
    // if the tau daughter is a tau, it means the particle has still to be propagated.
    // In that case, return the result of the same method on that daughter.
    if(daughter->pdgId()==tau.pdgId()) return LeptonicTauDecay(*daughter, elecdec, muondec, pid, px, py, pz, energy);
    // check for leptons
    elecdec |= abs(daughter->pdgId()) == 11;
    muondec |= abs(daughter->pdgId()) == 13;
    
    if (abs(daughter->pdgId()) == 11 || abs(daughter->pdgId()) == 13) {
      px = daughter->px();
      py = daughter->py();
      pz = daughter->pz();
      energy = daughter->energy();
      pid = daughter->pdgId();
      return;
    }
  }
  
}

//////


//------



// ------------ method called once each job just before starting event loop  ------------
void 
SUSYSkimToTree::beginJob()
{
  
  theHistosFile = new TFile(theHistosFileName.c_str(), "RECREATE");
  theHistosFile->cd();
  
  Tree = new TTree("Tree","");
  //Events
  Tree->Branch("T_Event_Rho", &T_Event_Rho, "T_Event_Rho/F");
  Tree->Branch("T_Event_RhoIso", &T_Event_RhoIso, "T_Event_RhoIso/F");
  /*  Tree->Branch("T_Event_RhoNoPu", &T_Event_RhoNoPu, "T_Event_RhoNoPu/F");
      Tree->Branch("T_Event_RhoIsoNoPu", &T_Event_RhoIsoNoPu, "T_Event_RhoIsoNoPu/F");*/

  //  Tree->Branch("T_Event_RhoCentralNeutral", &T_Event_RhoCentralNeutral, "T_Event_RhoCentralNeutral/F");
  //  Tree->Branch("T_Event_RhoCentralNeutralTight", &T_Event_RhoCentralNeutralTight, "T_Event_RhoCentralNeutralTight/F");

  Tree->Branch("T_EventF_EcalDeadCell", &T_EventF_EcalDeadCell, "T_EventF_EcalDeadCell/O");
  Tree->Branch("T_EventF_logErrorTooManyClusters", &T_EventF_logErrorTooManyClusters,"T_EventF_logErrorTooManyClusters/O");
  Tree->Branch("T_EventF_trackingFailure", &T_EventF_trackingFailure, "T_EventF_trackingFailure/O");
  Tree->Branch("T_EventF_hcalLaser", &T_EventF_hcalLaser, "T_EventF_hcalLaser/O");
  Tree->Branch("T_EventF_ecalLaserCorr", &T_EventF_ecalLaserCorr, "T_EventF_ecalLaserCorr/O");
  Tree->Branch("T_EventF_toomanystripclus", &T_EventF_toomanystripclus, "T_EventF_toomanystripclus/O");
  Tree->Branch("T_EventF_manystripclus",&T_EventF_manystripclus, "T_EventF_manystripclus/O");
  Tree->Branch("T_EventF_eeBadSc", &T_EventF_eeBadSc, "T_EventF_eeBadSc/O");

  Tree->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
  Tree->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
  Tree->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");
  //  Tree->Branch("T_Event_KFactorHNLO", &T_Event_KFactorHNLO, "T_Event_KFactorHNLO/F");
  //  if(!IsRealData){
  //  Tree->Branch("T_Event_PtHat", &T_Event_PtHat, "T_Event_PtHat/F");
  //  Tree->Branch("T_Event_HiggsPt", &T_Event_HiggsPt, "T_Event_HiggsPt/F");
  Tree->Branch("T_Event_processID", &T_Event_processID, "T_Event_processID/I");
  //  }
  
  Tree->Branch("T_Event_nPU", &T_Event_nPU, "T_Event_nPU/I");
  Tree->Branch("T_Event_nTruePU", &T_Event_nTruePU, "T_Event_nTruePU/F");
  Tree->Branch("T_Event_nPUm", &T_Event_nPUm, "T_Event_nPUm/I");
  Tree->Branch("T_Event_nPUp", &T_Event_nPUp, "T_Event_nPUp/I");
  Tree->Branch("T_Event_AveNTruePU", &T_Event_AveNTruePU, "T_Event_AveNTruePU/F"); 
  //HLT
  /*&
    Tree->Branch("T_HLT_Mu9", &T_HLT_Mu9, "T_HLT_Mu9/O");
    Tree->Branch("T_HLT_Mu8_v1", &T_HLT_Mu8_v1, "T_HLT_Mu8_v1/O");
    Tree->Branch("T_HLT_Mu8_v2", &T_HLT_Mu8_v2, "T_HLT_Mu8_v2/O");
    Tree->Branch("T_HLT_Mu8_v3", &T_HLT_Mu8_v3, "T_HLT_Mu8_v3/O");
    Tree->Branch("T_HLT_Mu8_v4", &T_HLT_Mu8_v4, "T_HLT_Mu8_v4/O");
    Tree->Branch("T_HLT_Mu8_v5", &T_HLT_Mu8_v5, "T_HLT_Mu8_v5/O");
    Tree->Branch("T_HLT_Mu8_v6", &T_HLT_Mu8_v6, "T_HLT_Mu8_v6/O");
    Tree->Branch("T_HLT_Mu8_v7", &T_HLT_Mu8_v7, "T_HLT_Mu8_v7/O");
    Tree->Branch("T_HLT_Mu8_v8", &T_HLT_Mu8_v8, "T_HLT_Mu8_v8/O");
    Tree->Branch("T_HLT_Mu8_v9", &T_HLT_Mu8_v9, "T_HLT_Mu8_v9/O");
    Tree->Branch("T_HLT_Mu8_v10", &T_HLT_Mu8_v10, "T_HLT_Mu8_v10/O");
    Tree->Branch("T_HLT_Mu8_v11", &T_HLT_Mu8_v11, "T_HLT_Mu8_v11/O");
    Tree->Branch("T_HLT_Mu8_v12", &T_HLT_Mu8_v12, "T_HLT_Mu8_v12/O");
  */
  Tree->Branch("T_HLT_Mu8_v16", &T_HLT_Mu8_v16, "T_HLT_Mu8_v16/O");
  /*
    Tree->Branch("T_HLT_Mu12_v1", &T_HLT_Mu12_v1, "T_HLT_Mu12_v1/O");
    Tree->Branch("T_HLT_Mu12_v2", &T_HLT_Mu12_v2, "T_HLT_Mu12_v2/O");
    Tree->Branch("T_HLT_Mu12_v3", &T_HLT_Mu12_v3, "T_HLT_Mu12_v3/O");
    Tree->Branch("T_HLT_Mu12_v4", &T_HLT_Mu12_v4, "T_HLT_Mu12_v4/O");
    Tree->Branch("T_HLT_Mu12_v16", &T_HLT_Mu12_v16, "T_HLT_Mu12_v16/O");

    Tree->Branch("T_HLT_Mu15_v1", &T_HLT_Mu15_v1, "T_HLT_Mu15_v1/O");
    Tree->Branch("T_HLT_Mu15_v2", &T_HLT_Mu15_v2, "T_HLT_Mu15_v2/O");
    Tree->Branch("T_HLT_Mu15_v3", &T_HLT_Mu15_v3, "T_HLT_Mu15_v3/O");
    Tree->Branch("T_HLT_Mu15_v4", &T_HLT_Mu15_v4, "T_HLT_Mu15_v4/O");
    Tree->Branch("T_HLT_Mu15_v5", &T_HLT_Mu15_v5, "T_HLT_Mu15_v5/O");
    Tree->Branch("T_HLT_Mu15_v6", &T_HLT_Mu15_v6, "T_HLT_Mu15_v6/O");
    Tree->Branch("T_HLT_Mu15_v7", &T_HLT_Mu15_v7, "T_HLT_Mu15_v7/O");
    Tree->Branch("T_HLT_Mu15_v8", &T_HLT_Mu15_v8, "T_HLT_Mu15_v8/O");
    Tree->Branch("T_HLT_Mu15_v9", &T_HLT_Mu15_v9, "T_HLT_Mu15_v9/O");
    Tree->Branch("T_HLT_Mu15_v10", &T_HLT_Mu15_v10, "T_HLT_Mu15_v10/O");
    Tree->Branch("T_HLT_Mu15_v11", &T_HLT_Mu15_v11, "T_HLT_Mu15_v11/O");
    Tree->Branch("T_HLT_Mu15_v12", &T_HLT_Mu15_v12, "T_HLT_Mu15_v12/O");
    Tree->Branch("T_HLT_Mu15_v13", &T_HLT_Mu15_v13, "T_HLT_Mu15_v13/O");
  */
  Tree->Branch("T_HLT_Mu17_v3", &T_HLT_Mu17_v3, "T_HLT_Mu17_v3/O");


  Tree->Branch("T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12", &T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12, "T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12/O");
  Tree->Branch("T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13", &T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13, "T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13/O");
  Tree->Branch("T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14", &T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14, "T_HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14/O");

  Tree->Branch("T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3", &T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3, "T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3/O");
  Tree->Branch("T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4", &T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4, "T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4/O");
  Tree->Branch("T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5", &T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5, "T_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5/O");



  /*
    Tree->Branch("T_HLT_IsoMu24_eta2p1_v11", &T_HLT_IsoMu24_eta2p1_v11, "T_HLT_IsoMu24_eta2p1_v11/O");
    Tree->Branch("T_HLT_IsoMu24_eta2p1_v12", &T_HLT_IsoMu24_eta2p1_v12, "T_HLT_IsoMu24_eta2p1_v12/O");

    Tree->Branch("T_HLT_Mu24_v1", &T_HLT_Mu24_v1, "T_HLT_Mu24_v1/O");
    Tree->Branch("T_HLT_Mu24_v2", &T_HLT_Mu24_v2, "T_HLT_Mu24_v2/O");
    Tree->Branch("T_HLT_Mu24_v3", &T_HLT_Mu24_v3, "T_HLT_Mu24_v3/O");
    Tree->Branch("T_HLT_Mu24_v4", &T_HLT_Mu24_v4, "T_HLT_Mu24_v4/O");
    Tree->Branch("T_HLT_Mu24_v5", &T_HLT_Mu24_v5, "T_HLT_Mu24_v5/O");
    Tree->Branch("T_HLT_Mu24_v6", &T_HLT_Mu24_v6, "T_HLT_Mu24_v6/O");
    Tree->Branch("T_HLT_Mu24_v7", &T_HLT_Mu24_v7, "T_HLT_Mu24_v7/O");
    Tree->Branch("T_HLT_Mu24_v8", &T_HLT_Mu24_v8, "T_HLT_Mu24_v8/O");
    Tree->Branch("T_HLT_Mu24_v9", &T_HLT_Mu24_v9, "T_HLT_Mu24_v9/O");
    Tree->Branch("T_HLT_Mu24_v10", &T_HLT_Mu24_v10, "T_HLT_Mu24_v10/O");
    Tree->Branch("T_HLT_Mu24_v11", &T_HLT_Mu24_v11, "T_HLT_Mu24_v11/O");
    Tree->Branch("T_HLT_Mu24_v12", &T_HLT_Mu24_v12, "T_HLT_Mu24_v12/O");

    Tree->Branch("T_HLT_Jet30_v1", &T_HLT_Jet30_v1, "T_HLT_Jet30_v1/O");
    Tree->Branch("T_HLT_Jet30_v2", &T_HLT_Jet30_v2, "T_HLT_Jet30_v2/O");
    Tree->Branch("T_HLT_Jet30_v3", &T_HLT_Jet30_v3, "T_HLT_Jet30_v3/O");
    Tree->Branch("T_HLT_Jet30_v4", &T_HLT_Jet30_v4, "T_HLT_Jet30_v4/O");
    Tree->Branch("T_HLT_Jet30_v5", &T_HLT_Jet30_v5, "T_HLT_Jet30_v5/O");
    Tree->Branch("T_HLT_Jet60_v1", &T_HLT_Jet60_v1, "T_HLT_Jet60_v1/O");
    Tree->Branch("T_HLT_Jet60_v2", &T_HLT_Jet60_v2, "T_HLT_Jet60_v2/O");
    Tree->Branch("T_HLT_Jet60_v3", &T_HLT_Jet60_v3, "T_HLT_Jet60_v3/O");
    Tree->Branch("T_HLT_Jet60_v4", &T_HLT_Jet60_v4, "T_HLT_Jet60_v4/O");
    Tree->Branch("T_HLT_Jet60_v5", &T_HLT_Jet60_v5, "T_HLT_Jet60_v5/O");
  */
  //  if(!IsRealData){
  //Gen 

  Tree->Branch("T_Gen_StopMass", "std::vector<float>", &T_Gen_StopMass);
  Tree->Branch("T_Gen_Chi0Mass", "std::vector<float>", &T_Gen_Chi0Mass);
  Tree->Branch("T_Gen_CharginoMass", "std::vector<float>", &T_Gen_CharginoMass);

  Tree->Branch("T_Gen_polWeights", "std::vector<double>", &T_Gen_polWeights);

  Tree->Branch("T_Gen_Muon_PID", "std::vector<int>", &T_Gen_Muon_PID);
  Tree->Branch("T_Gen_Muon_Px", "std::vector<float>", &T_Gen_Muon_Px);
  Tree->Branch("T_Gen_Muon_Py", "std::vector<float>", &T_Gen_Muon_Py);
  Tree->Branch("T_Gen_Muon_Pz", "std::vector<float>", &T_Gen_Muon_Pz);
  Tree->Branch("T_Gen_Muon_Energy", "std::vector<float>", &T_Gen_Muon_Energy);
  Tree->Branch("T_Gen_Muon_MPID", "std::vector<int>", &T_Gen_Muon_MPID);
  Tree->Branch("T_Gen_Muon_MPx", "std::vector<float>", &T_Gen_Muon_MPx);
  Tree->Branch("T_Gen_Muon_MPy", "std::vector<float>", &T_Gen_Muon_MPy);
  Tree->Branch("T_Gen_Muon_MPz", "std::vector<float>", &T_Gen_Muon_MPz);
  Tree->Branch("T_Gen_Muon_MEnergy", "std::vector<float>", &T_Gen_Muon_MEnergy);
  Tree->Branch("T_Gen_Muon_MSt", "std::vector<int>", &T_Gen_Muon_MSt);

  Tree->Branch("T_Gen_Elec_PID", "std::vector<int>", &T_Gen_Elec_PID);
  Tree->Branch("T_Gen_Elec_Px", "std::vector<float>", &T_Gen_Elec_Px);
  Tree->Branch("T_Gen_Elec_Py", "std::vector<float>", &T_Gen_Elec_Py);
  Tree->Branch("T_Gen_Elec_Pz", "std::vector<float>", &T_Gen_Elec_Pz);
  Tree->Branch("T_Gen_Elec_Energy", "std::vector<float>", &T_Gen_Elec_Energy);
  Tree->Branch("T_Gen_Elec_MPID", "std::vector<int>", &T_Gen_Elec_MPID);
  Tree->Branch("T_Gen_Elec_MPx", "std::vector<float>", &T_Gen_Elec_MPx);
  Tree->Branch("T_Gen_Elec_MPy", "std::vector<float>", &T_Gen_Elec_MPy);
  Tree->Branch("T_Gen_Elec_MPz", "std::vector<float>", &T_Gen_Elec_MPz);
  Tree->Branch("T_Gen_Elec_MEnergy", "std::vector<float>", &T_Gen_Elec_MEnergy);
  Tree->Branch("T_Gen_Elec_MSt", "std::vector<int>", &T_Gen_Elec_MSt);

  Tree->Branch("T_Gen_b_PID", "std::vector<int>", &T_Gen_b_PID);
  Tree->Branch("T_Gen_b_Px", "std::vector<float>", &T_Gen_b_Px);
  Tree->Branch("T_Gen_b_Py", "std::vector<float>", &T_Gen_b_Py);
  Tree->Branch("T_Gen_b_Pz", "std::vector<float>", &T_Gen_b_Pz);
  Tree->Branch("T_Gen_b_Energy", "std::vector<float>", &T_Gen_b_Energy);
  Tree->Branch("T_Gen_b_MPID", "std::vector<int>", &T_Gen_b_MPID);
  Tree->Branch("T_Gen_b_MPx", "std::vector<float>", &T_Gen_b_MPx);
  Tree->Branch("T_Gen_b_MPy", "std::vector<float>", &T_Gen_b_MPy);
  Tree->Branch("T_Gen_b_MPz", "std::vector<float>", &T_Gen_b_MPz);
  Tree->Branch("T_Gen_b_MEnergy", "std::vector<float>", &T_Gen_b_MEnergy);
  Tree->Branch("T_Gen_b_MSt", "std::vector<int>", &T_Gen_b_MSt);

  
  
     
 
     
     
  
  

  /*
  Tree->Branch("T_Gen_StopSt3  ", "std::vector<SUSYGenParticle>", &T_Gen_StopSt3  );
  Tree->Branch("T_Gen_Chi0St3  ", "std::vector<SUSYGenParticle>", &T_Gen_Chi0St3  );
  Tree->Branch("T_Gen_tSt3     ", "std::vector<SUSYGenParticle>", &T_Gen_tSt3     );
  Tree->Branch("T_Gen_ChiPMSt3 ", "std::vector<SUSYGenParticle>", &T_Gen_ChiPMSt3 );
  Tree->Branch("T_Gen_bSt3     ", "std::vector<SUSYGenParticle>", &T_Gen_bSt3     );
  Tree->Branch("T_Gen_WSt3     ", "std::vector<SUSYGenParticle>", &T_Gen_WSt3     );
  Tree->Branch("T_Gen_MuonSt3  ", "std::vector<SUSYGenParticle>", &T_Gen_MuonSt3  );
  Tree->Branch("T_Gen_ElecSt3  ", "std::vector<SUSYGenParticle>", &T_Gen_ElecSt3  );
  Tree->Branch("T_Gen_TauSt3   ", "std::vector<SUSYGenParticle>", &T_Gen_TauSt3   );
  */
  Tree->Branch("T_Gen_MuonSt3_pdgId", "std::vector<int>",         &T_Gen_MuonSt3_pdgId       );
  Tree->Branch("T_Gen_MuonSt3_firstMother", "std::vector<int>",   &T_Gen_MuonSt3_firstMother );	   
  Tree->Branch("T_Gen_MuonSt3_i", "std::vector<int>",             &T_Gen_MuonSt3_i           );	   
  Tree->Branch("T_Gen_MuonSt3_energy", "std::vector<float>",      &T_Gen_MuonSt3_energy      );	    
  Tree->Branch("T_Gen_MuonSt3_pt", "std::vector<float>",          &T_Gen_MuonSt3_pt          );	    
  Tree->Branch("T_Gen_MuonSt3_eta", "std::vector<float>",         &T_Gen_MuonSt3_eta         );	    
  Tree->Branch("T_Gen_MuonSt3_phi", "std::vector<float>",         &T_Gen_MuonSt3_phi         );

  Tree->Branch("T_Gen_ElecSt3_pdgId", "std::vector<int>",         &T_Gen_ElecSt3_pdgId       );
  Tree->Branch("T_Gen_ElecSt3_firstMother", "std::vector<int>",   &T_Gen_ElecSt3_firstMother );
  Tree->Branch("T_Gen_ElecSt3_i", "std::vector<int>",             &T_Gen_ElecSt3_i           );	   
  Tree->Branch("T_Gen_ElecSt3_energy", "std::vector<float>",      &T_Gen_ElecSt3_energy      );
  Tree->Branch("T_Gen_ElecSt3_pt", "std::vector<float>",          &T_Gen_ElecSt3_pt          );
  Tree->Branch("T_Gen_ElecSt3_eta", "std::vector<float>",         &T_Gen_ElecSt3_eta         );
  Tree->Branch("T_Gen_ElecSt3_phi", "std::vector<float>",         &T_Gen_ElecSt3_phi         );

  Tree->Branch("T_Gen_TauSt3_pdgId", "std::vector<int>",          &T_Gen_TauSt3_pdgId       );
  Tree->Branch("T_Gen_TauSt3_firstMother", "std::vector<int>",    &T_Gen_TauSt3_firstMother );
  Tree->Branch("T_Gen_TauSt3_i", "std::vector<int>",              &T_Gen_TauSt3_i           );	   
  Tree->Branch("T_Gen_TauSt3_energy", "std::vector<float>",       &T_Gen_TauSt3_energy      );
  Tree->Branch("T_Gen_TauSt3_pt", "std::vector<float>",           &T_Gen_TauSt3_pt          );
  Tree->Branch("T_Gen_TauSt3_eta", "std::vector<float>",          &T_Gen_TauSt3_eta         );
  Tree->Branch("T_Gen_TauSt3_phi", "std::vector<float>",          &T_Gen_TauSt3_phi         );

  Tree->Branch("T_Gen_StopSt3_pdgId", "std::vector<int>",         &T_Gen_StopSt3_pdgId       );
  Tree->Branch("T_Gen_StopSt3_firstMother", "std::vector<int>",   &T_Gen_StopSt3_firstMother );
  Tree->Branch("T_Gen_StopSt3_i", "std::vector<int>",             &T_Gen_StopSt3_i           );
  Tree->Branch("T_Gen_StopSt3_energy", "std::vector<float>",      &T_Gen_StopSt3_energy      );
  Tree->Branch("T_Gen_StopSt3_pt", "std::vector<float>",          &T_Gen_StopSt3_pt          );
  Tree->Branch("T_Gen_StopSt3_eta", "std::vector<float>",         &T_Gen_StopSt3_eta         );
  Tree->Branch("T_Gen_StopSt3_phi", "std::vector<float>",         &T_Gen_StopSt3_phi         );
  
  Tree->Branch("T_Gen_Chi0St3_pdgId", "std::vector<int>",         &T_Gen_Chi0St3_pdgId       );
  Tree->Branch("T_Gen_Chi0St3_firstMother", "std::vector<int>",   &T_Gen_Chi0St3_firstMother );
  Tree->Branch("T_Gen_Chi0St3_i", "std::vector<int>",             &T_Gen_Chi0St3_i           );
  Tree->Branch("T_Gen_Chi0St3_energy", "std::vector<float>",      &T_Gen_Chi0St3_energy      );
  Tree->Branch("T_Gen_Chi0St3_pt", "std::vector<float>",          &T_Gen_Chi0St3_pt          );
  Tree->Branch("T_Gen_Chi0St3_eta", "std::vector<float>",         &T_Gen_Chi0St3_eta         );
  Tree->Branch("T_Gen_Chi0St3_phi", "std::vector<float>",         &T_Gen_Chi0St3_phi         );
  
  Tree->Branch("T_Gen_ChiPMSt3_pdgId", "std::vector<int>",         &T_Gen_ChiPMSt3_pdgId       );
  Tree->Branch("T_Gen_ChiPMSt3_firstMother", "std::vector<int>",   &T_Gen_ChiPMSt3_firstMother );
  Tree->Branch("T_Gen_ChiPMSt3_i", "std::vector<int>",             &T_Gen_ChiPMSt3_i           );
  Tree->Branch("T_Gen_ChiPMSt3_energy", "std::vector<float>",      &T_Gen_ChiPMSt3_energy      );
  Tree->Branch("T_Gen_ChiPMSt3_pt", "std::vector<float>",          &T_Gen_ChiPMSt3_pt          );
  Tree->Branch("T_Gen_ChiPMSt3_eta", "std::vector<float>",         &T_Gen_ChiPMSt3_eta         );
  Tree->Branch("T_Gen_ChiPMSt3_phi", "std::vector<float>",         &T_Gen_ChiPMSt3_phi         );

  Tree->Branch("T_Gen_bSt3_pdgId", "std::vector<int>",         &T_Gen_bSt3_pdgId       );
  Tree->Branch("T_Gen_bSt3_firstMother", "std::vector<int>",   &T_Gen_bSt3_firstMother );
  Tree->Branch("T_Gen_bSt3_i", "std::vector<int>",             &T_Gen_bSt3_i           );
  Tree->Branch("T_Gen_bSt3_energy", "std::vector<float>",      &T_Gen_bSt3_energy      );
  Tree->Branch("T_Gen_bSt3_pt", "std::vector<float>",          &T_Gen_bSt3_pt          );
  Tree->Branch("T_Gen_bSt3_eta", "std::vector<float>",         &T_Gen_bSt3_eta         );
  Tree->Branch("T_Gen_bSt3_phi", "std::vector<float>",         &T_Gen_bSt3_phi         );

  Tree->Branch("T_Gen_tSt3_pdgId", "std::vector<int>",         &T_Gen_tSt3_pdgId       );
  Tree->Branch("T_Gen_tSt3_firstMother", "std::vector<int>",   &T_Gen_tSt3_firstMother );
  Tree->Branch("T_Gen_tSt3_i", "std::vector<int>",             &T_Gen_tSt3_i           );
  Tree->Branch("T_Gen_tSt3_energy", "std::vector<float>",      &T_Gen_tSt3_energy      );
  Tree->Branch("T_Gen_tSt3_pt", "std::vector<float>",          &T_Gen_tSt3_pt          );
  Tree->Branch("T_Gen_tSt3_eta", "std::vector<float>",         &T_Gen_tSt3_eta         );
  Tree->Branch("T_Gen_tSt3_phi", "std::vector<float>",         &T_Gen_tSt3_phi         );

  Tree->Branch("T_Gen_WSt3_pdgId", "std::vector<int>",         &T_Gen_WSt3_pdgId       );
  Tree->Branch("T_Gen_WSt3_firstMother", "std::vector<int>",   &T_Gen_WSt3_firstMother );
  Tree->Branch("T_Gen_WSt3_i", "std::vector<int>",             &T_Gen_WSt3_i           );
  Tree->Branch("T_Gen_WSt3_energy", "std::vector<float>",      &T_Gen_WSt3_energy      );
  Tree->Branch("T_Gen_WSt3_pt", "std::vector<float>",          &T_Gen_WSt3_pt          );
  Tree->Branch("T_Gen_WSt3_eta", "std::vector<float>",         &T_Gen_WSt3_eta         );
  Tree->Branch("T_Gen_WSt3_phi", "std::vector<float>",         &T_Gen_WSt3_phi         );






  /*
  Tree->Branch("T_Gen_MuonSt3_PID", "std::vector<int>", &T_Gen_MuonSt3_PID);	    
  Tree->Branch("T_Gen_MuonSt3_Px", "std::vector<float>", &T_Gen_MuonSt3_Px);	    
  Tree->Branch("T_Gen_MuonSt3_Py", "std::vector<float>", &T_Gen_MuonSt3_Py);	    
  Tree->Branch("T_Gen_MuonSt3_Pz", "std::vector<float>", &T_Gen_MuonSt3_Pz);	    
  Tree->Branch("T_Gen_MuonSt3_Energy", "std::vector<float>", &T_Gen_MuonSt3_Energy);
  Tree->Branch("T_Gen_ElecSt3_PID", "std::vector<int>", &T_Gen_ElecSt3_PID);
  Tree->Branch("T_Gen_ElecSt3_Px", "std::vector<float>", &T_Gen_ElecSt3_Px);
  Tree->Branch("T_Gen_ElecSt3_Py", "std::vector<float>", &T_Gen_ElecSt3_Py);
  Tree->Branch("T_Gen_ElecSt3_Pz", "std::vector<float>", &T_Gen_ElecSt3_Pz);
  Tree->Branch("T_Gen_ElecSt3_Energy", "std::vector<float>", &T_Gen_ElecSt3_Energy);
  Tree->Branch("T_Gen_bSt3_PID", "std::vector<int>", &T_Gen_bSt3_PID);
  Tree->Branch("T_Gen_bSt3_Px", "std::vector<float>", &T_Gen_bSt3_Px);
  Tree->Branch("T_Gen_bSt3_Py", "std::vector<float>", &T_Gen_bSt3_Py);
  Tree->Branch("T_Gen_bSt3_Pz", "std::vector<float>", &T_Gen_bSt3_Pz);
  Tree->Branch("T_Gen_bSt3_Energy", "std::vector<float>", &T_Gen_bSt3_Energy);
  Tree->Branch("T_Gen_tSt3_PID", "std::vector<int>", &T_Gen_tSt3_PID);
  Tree->Branch("T_Gen_tSt3_Px", "std::vector<float>", &T_Gen_tSt3_Px);
  Tree->Branch("T_Gen_tSt3_Py", "std::vector<float>", &T_Gen_tSt3_Py);
  Tree->Branch("T_Gen_tSt3_Pz", "std::vector<float>", &T_Gen_tSt3_Pz);
  Tree->Branch("T_Gen_tSt3_Energy", "std::vector<float>", &T_Gen_tSt3_Energy);

  Tree->Branch("T_Gen_TauSt3_PID", "std::vector<int>", &T_Gen_TauSt3_PID);
  Tree->Branch("T_Gen_TauSt3_Px", "std::vector<float>", &T_Gen_TauSt3_Px);
  Tree->Branch("T_Gen_TauSt3_Py", "std::vector<float>", &T_Gen_TauSt3_Py);
  Tree->Branch("T_Gen_TauSt3_Pz", "std::vector<float>", &T_Gen_TauSt3_Pz);
  Tree->Branch("T_Gen_TauSt3_Energy", "std::vector<float>", &T_Gen_TauSt3_Energy);
  */

  Tree->Branch("T_Gen_TauSt3_IsLepDec", "std::vector<bool>", &T_Gen_TauSt3_IsLepDec);
  Tree->Branch("T_Gen_TauSt3_LepDec_PID", "std::vector<int>", &T_Gen_TauSt3_LepDec_PID);
  Tree->Branch("T_Gen_TauSt3_LepDec_Px", "std::vector<float>", &T_Gen_TauSt3_LepDec_Px);
  Tree->Branch("T_Gen_TauSt3_LepDec_Py", "std::vector<float>", &T_Gen_TauSt3_LepDec_Py);
  Tree->Branch("T_Gen_TauSt3_LepDec_Pz", "std::vector<float>", &T_Gen_TauSt3_LepDec_Pz);
  Tree->Branch("T_Gen_TauSt3_LepDec_Energy", "std::vector<float>", &T_Gen_TauSt3_LepDec_Energy);
  //}
  //PFcands
  /*  Tree->Branch("T_PFreducedCand_Px", "std::vector<float>", &T_PFreducedCand_Px);
      Tree->Branch("T_PFreducedCand_Py", "std::vector<float>", &T_PFreducedCand_Py);
      Tree->Branch("T_PFreducedCand_Pz", "std::vector<float>", &T_PFreducedCand_Pz);
      Tree->Branch("T_PFreducedCand_E", "std::vector<float>", &T_PFreducedCand_E);
      Tree->Branch("T_PFreducedCand_ID", "std::vector<int>", &T_PFreducedCand_ID);
      Tree->Branch("T_PFreducedCand_vz", "std::vector<float>", &T_PFreducedCand_vz);
  */
   
  //Muons
  Tree->Branch("T_Muon_Eta", "std::vector<float>", &T_Muon_Eta);
  Tree->Branch("T_Muon_IsGlobalMuon", "std::vector<bool>", &T_Muon_IsGlobalMuon);
  //  Tree->Branch("T_Muon_IsTMLSOLPT", "std::vector<bool>", &T_Muon_IsTMLSOLPT);
  //  Tree->Branch("T_Muon_IsTMLSOLPL", "std::vector<bool>", &T_Muon_IsTMLSOLPL);
   
  //  Tree-> Branch("T_Muon_IsAllStandAloneMuons", "std::vector<bool>", &T_Muon_IsAllStandAloneMuons);
  Tree-> Branch("T_Muon_IsAllTrackerMuons", "std::vector<bool>", &T_Muon_IsAllTrackerMuons);
  Tree-> Branch("T_Muon_IsTrackerMuonArbitrated", "std::vector<bool>", &T_Muon_IsTrackerMuonArbitrated);
  /*  Tree-> Branch("T_Muon_IsAllArbitrated", "std::vector<bool>", &T_Muon_IsAllArbitrated);
      Tree-> Branch("T_Muon_IsTMLastStationLoose", "std::vector<bool>", &T_Muon_IsTMLastStationLoose);
      Tree-> Branch("T_Muon_IsTMLastStationTight", "std::vector<bool>", &T_Muon_IsTMLastStationTight);
      Tree-> Branch("T_Muon_IsTM2DCompatibilityLoose", "std::vector<bool>", &T_Muon_IsTM2DCompatibilityLoose);
      Tree-> Branch("T_Muon_IsTM2DCompatibilityTight", "std::vector<bool>", &T_Muon_IsTM2DCompatibilityTight);
      Tree-> Branch("T_Muon_IsTMOneStationLoose", "std::vector<bool>", &T_Muon_IsTMOneStationLoose);
      Tree-> Branch("T_Muon_IsTMOneStationTight", "std::vector<bool>", &T_Muon_IsTMOneStationTight);
      Tree-> Branch("T_Muon_IsTMLSOPL", "std::vector<bool>", &T_Muon_IsTMLSOPL);
      Tree-> Branch("T_Muon_IsGMTkChiCompatibility", "std::vector<bool>", &T_Muon_IsGMTkChiCompatibility);
      Tree-> Branch("T_Muon_IsGMStaChiCompatibility", "std::vector<bool>", &T_Muon_IsGMStaChiCompatibility);
      Tree-> Branch("T_Muon_IsGMTkKinkTight", "std::vector<bool>", &T_Muon_IsGMTkKinkTight);
      Tree-> Branch("T_Muon_IsTMLastStationAngLoose","std::vector<bool>", &T_Muon_IsTMLastStationAngLoose);
      Tree-> Branch("T_Muon_IsTMLastStationAngTight","std::vector<bool>", &T_Muon_IsTMLastStationAngTight);
      Tree-> Branch("T_Muon_IsTMOneStationAngLoose","std::vector<bool>", &T_Muon_IsTMOneStationAngLoose);
      Tree-> Branch("T_Muon_IsTMOneStationAngTight","std::vector<bool>", &T_Muon_IsTMOneStationAngTight);
  */
  Tree-> Branch("T_Muon_IsGMPTMuons", "std::vector<bool>", &T_Muon_IsGMPTMuons);
 
  Tree->Branch("T_Muon_SegmentCompatibility","std::vector<float>", &T_Muon_SegmentCompatibility);
  Tree->Branch("T_Muon_trkKink","std::vector<float>", &T_Muon_trkKink);
  Tree->Branch("T_Muon_Px", "std::vector<float>", &T_Muon_Px);
  Tree->Branch("T_Muon_Py", "std::vector<float>", &T_Muon_Py);
  Tree->Branch("T_Muon_Pz", "std::vector<float>", &T_Muon_Pz);
  Tree->Branch("T_Muon_Pt", "std::vector<float>", &T_Muon_Pt);
  Tree->Branch("T_Muon_deltaPt", "std::vector<float>", &T_Muon_deltaPt);
  Tree->Branch("T_Muon_Energy", "std::vector<float>", &T_Muon_Energy);
  Tree->Branch("T_Muon_Charge", "std::vector<int>", &T_Muon_Charge);
  Tree->Branch("T_Muon_NormChi2GTrk", "std::vector<float>", &T_Muon_NormChi2GTrk);
  Tree->Branch("T_Muon_NValidHitsInTrk", "std::vector<int>", &T_Muon_NValidHitsInTrk);
  Tree->Branch("T_Muon_NValidPixelHitsInTrk", "std::vector<int>", &T_Muon_NValidPixelHitsInTrk);
  Tree->Branch("T_Muon_InnerTrackFound", "std::vector<int>", &T_Muon_InnerTrackFound);
  Tree->Branch("T_Muon_NValidHitsSATrk", "std::vector<int>", &T_Muon_NValidHitsSATrk);
  Tree->Branch("T_Muon_NValidHitsGTrk", "std::vector<int>", &T_Muon_NValidHitsGTrk);
  Tree->Branch("T_Muon_Chi2InTrk", "std::vector<float>", &T_Muon_Chi2InTrk);
  Tree->Branch("T_Muon_dofInTrk", "std::vector<float>", &T_Muon_dofInTrk);
  Tree->Branch("T_Muon_IPAbsGTrack", "std::vector<float>", &T_Muon_IPAbsGTrack);
  //  Tree->Branch("T_Muon_IPSigGTrack", "std::vector<float>", &T_Muon_IPSigGTrack);
  Tree->Branch("T_Muon_IPAbsInTrack", "std::vector<float>", &T_Muon_IPAbsInTrack);
  //  Tree->Branch("T_Muon_IPSigInTrack", "std::vector<float>", &T_Muon_IPSigInTrack);
  //  Tree->Branch("T_Muon_IPwrtBSInTrack", "std::vector<float>", &T_Muon_IPwrtBSInTrack);
  Tree->Branch("T_Muon_IPwrtAveBSInTrack", "std::vector<float>", &T_Muon_IPwrtAveBSInTrack);
  //  Tree->Branch("T_Muon_IPwrtBSGTrack", "std::vector<float>", &T_Muon_IPwrtBSGTrack);
  /*
    Tree->Branch("T_Muon_IP2DBiasedPV", "std::vector<float>", &T_Muon_IP2DBiasedPV);
    Tree->Branch("T_Muon_IP3DBiasedPV", "std::vector<float>", &T_Muon_IP3DBiasedPV);
    Tree->Branch("T_Muon_IP2DUnBiasedPV", "std::vector<float>", &T_Muon_IP2DUnBiasedPV);
    Tree->Branch("T_Muon_IP3DUnBiasedPV", "std::vector<float>", &T_Muon_IP3DUnBiasedPV);
    Tree->Branch("T_Muon_dxyPVBiasedPV", "std::vector<float>", &T_Muon_dxyPVBiasedPV);
    Tree->Branch("T_Muon_dxyPVUnBiasedPV", "std::vector<float>", &T_Muon_dxyPVUnBiasedPV);
    Tree->Branch("T_Muon_dzPVBiasedPV", "std::vector<float>", &T_Muon_dzPVBiasedPV);
    Tree->Branch("T_Muon_dzPVUnBiasedPV", "std::vector<float>", &T_Muon_dzPVUnBiasedPV);
    Tree->Branch("T_Muon_IP2DUnBiasedPVnoBS", "std::vector<float>", &T_Muon_IP2DUnBiasedPVnoBS);
    Tree->Branch("T_Muon_IP3DUnBiasedPVnoBS", "std::vector<float>", &T_Muon_IP3DUnBiasedPVnoBS);
    Tree->Branch("T_Muon_dxyPVUnBiasedPVnoBS", "std::vector<float>", &T_Muon_dxyPVUnBiasedPVnoBS);
    Tree->Branch("T_Muon_dzPVUnBiasedPVnoBS", "std::vector<float>", &T_Muon_dzPVUnBiasedPVnoBS);*/
  /*  Tree->Branch("T_Muon_smurfCharged", "std::vector<float>", &T_Muon_smurfCharged);
      Tree->Branch("T_Muon_smurfNeutral", "std::vector<float>", &T_Muon_smurfNeutral);
      Tree->Branch("T_Muon_smurfPhoton", "std::vector<float>", &T_Muon_smurfPhoton);
      Tree->Branch("T_Muon_smurfNoOverCharged", "std::vector<float>", &T_Muon_smurfNoOverCharged);
      Tree->Branch("T_Muon_smurfNoOverNeutral", "std::vector<float>", &T_Muon_smurfNoOverNeutral);
      Tree->Branch("T_Muon_smurfNoOverPhoton", "std::vector<float>", &T_Muon_smurfNoOverPhoton);
      Tree->Branch("T_Muon_muSmurfPF", "std::vector<float>", &T_Muon_muSmurfPF);
  */
  //  Tree->Branch("T_Muon_chargedParticleIsoR04", "std::vector<float>", &T_Muon_chargedParticleIsoR04);
  Tree->Branch("T_Muon_chargedHadronIsoR04", "std::vector<float>", &T_Muon_chargedHadronIsoR04);
  Tree->Branch("T_Muon_neutralHadronIsoR04", "std::vector<float>", &T_Muon_neutralHadronIsoR04);
  Tree->Branch("T_Muon_photonIsoR04", "std::vector<float>", &T_Muon_photonIsoR04);
  Tree->Branch("T_Muon_chargedParticleIsoR03", "std::vector<float>", &T_Muon_chargedParticleIsoR03);
  Tree->Branch("T_Muon_chargedHadronIsoR03", "std::vector<float>", &T_Muon_chargedHadronIsoR03);
  Tree->Branch("T_Muon_neutralHadronIsoR03", "std::vector<float>", &T_Muon_neutralHadronIsoR03);
  Tree->Branch("T_Muon_photonIsoR03", "std::vector<float>", &T_Muon_photonIsoR03);
  Tree->Branch("T_Muon_sumPUPtR04", "std::vector<float>", &T_Muon_sumPUPtR04);
  Tree->Branch("T_Muon_sumPUPtR03", "std::vector<float>", &T_Muon_sumPUPtR03);
  Tree->Branch("T_Muon_vz", "std::vector<float>", &T_Muon_vz);
  Tree->Branch("T_Muon_vy", "std::vector<float>", &T_Muon_vy);
  Tree->Branch("T_Muon_vx", "std::vector<float>", &T_Muon_vx);
  Tree->Branch("T_Muon_NumOfMatchedStations", "std::vector<int>", &T_Muon_NumOfMatchedStations);
  Tree->Branch("T_Muon_PFMuonPt","std::vector<float>", &T_Muon_PFMuonPt);
  Tree->Branch("T_Muon_PFMuonPx","std::vector<float>", &T_Muon_PFMuonPx);
  Tree->Branch("T_Muon_PFMuonPy","std::vector<float>", &T_Muon_PFMuonPy);
  Tree->Branch("T_Muon_PFMuonPz","std::vector<float>", &T_Muon_PFMuonPx);
  Tree->Branch("T_Muon_PFMuonE","std::vector<float>", &T_Muon_PFMuonE);
  Tree->Branch("T_Muon_isPFMuon","std::vector<bool>", &T_Muon_isPFMuon);
  Tree->Branch("T_Muon_NLayers","std::vector<int>", &T_Muon_NLayers);

  /*  Tree->Branch("T_Muon_passTriggerElMu","std::vector<bool>",&T_Muon_passTriggerElMu);
      Tree->Branch("T_Muon_passTriggerDoubleMu","std::vector<bool>",&T_Muon_passTriggerDoubleMu);
      Tree->Branch("T_Muon_passTriggerSingleMu","std::vector<bool>",&T_Muon_passTriggerSingleMu);
  */

  //pfTaus
  /*  Tree->Branch("T_pfTau_Px", "std::vector<float>", &T_pfTau_Px);
      Tree->Branch("T_pfTau_Py", "std::vector<float>", &T_pfTau_Py);
      Tree->Branch("T_pfTau_Pz", "std::vector<float>", &T_pfTau_Pz);
      Tree->Branch("T_pfTau_Energy", "std::vector<float>", &T_pfTau_Energy);
      Tree->Branch("T_pfTau_Charge", "std::vector<int>", &T_pfTau_Charge);
  */  
  //Taus
  /*  Tree->Branch("T_Tau_Px", "std::vector<float>", &T_Tau_Px);
      Tree->Branch("T_Tau_Py", "std::vector<float>", &T_Tau_Py);
      Tree->Branch("T_Tau_Pz", "std::vector<float>", &T_Tau_Pz);
      Tree->Branch("T_Tau_Energy", "std::vector<float>", &T_Tau_Energy);
      Tree->Branch("T_Tau_Charge", "std::vector<int>", &T_Tau_Charge);
  */  
  
  //pfMuons
  /*
    Tree->Branch("T_pfMuon_IsGlobalMuon", "std::vector<bool>", &T_pfMuon_IsGlobalMuon);
    Tree->Branch("T_pfMuon_IsAllStandAloneMuons", "std::vector<bool>", &T_pfMuon_IsAllStandAloneMuons);
    Tree->Branch("T_pfMuon_IsAllTrackerMuons", "std::vector<bool>", &T_pfMuon_IsAllTrackerMuons);
    Tree->Branch("T_pfMuon_IsTMLastStationAngTight","std::vector<bool>", &T_pfMuon_IsTMLastStationAngTight);
    Tree->Branch("T_pfMuon_IsGMPTMuons", "std::vector<bool>", &T_pfMuon_IsGMPTMuons);
    Tree->Branch("T_pfMuon_SegmentCompatibility","std::vector<float>", &T_pfMuon_SegmentCompatibility);
    Tree->Branch("T_pfMuon_Px", "std::vector<float>", &T_pfMuon_Px);
    Tree->Branch("T_pfMuon_Py", "std::vector<float>", &T_pfMuon_Py);
    Tree->Branch("T_pfMuon_Pz", "std::vector<float>", &T_pfMuon_Pz);
    Tree->Branch("T_pfMuon_Pt", "std::vector<float>", &T_pfMuon_Pt);
    Tree->Branch("T_pfMuon_Energy", "std::vector<float>", &T_pfMuon_Energy);
    Tree->Branch("T_pfMuon_Charge", "std::vector<int>", &T_pfMuon_Charge);
    Tree->Branch("T_pfMuon_NormChi2GTrk", "std::vector<float>", &T_pfMuon_NormChi2GTrk);
    Tree->Branch("T_pfMuon_NValidHitsInTrk", "std::vector<int>", &T_pfMuon_NValidHitsInTrk);
    Tree->Branch("T_pfMuon_NValidHitsSATrk", "std::vector<int>", &T_pfMuon_NValidHitsSATrk);
    Tree->Branch("T_pfMuon_Chi2InTrk", "std::vector<float>", &T_pfMuon_Chi2InTrk);
    Tree->Branch("T_pfMuon_dofInTrk", "std::vector<float>", &T_pfMuon_dofInTrk);
    Tree->Branch("T_pfMuon_SumIsoCalo", "std::vector<float>", &T_pfMuon_SumIsoCalo);
    Tree->Branch("T_pfMuon_SumIsoTrack", "std::vector<float>", &T_pfMuon_SumIsoTrack);
    Tree->Branch("T_pfMuon_IPAbsGTrack", "std::vector<float>", &T_pfMuon_IPAbsGTrack);
    Tree->Branch("T_pfMuon_IPSigGTrack", "std::vector<float>", &T_pfMuon_IPSigGTrack);
    Tree->Branch("T_pfMuon_IPAbsInTrack", "std::vector<float>", &T_pfMuon_IPAbsInTrack);
    Tree->Branch("T_pfMuon_IPSigInTrack", "std::vector<float>", &T_pfMuon_IPSigInTrack);
    //  Tree->Branch("T_pfMuon_IPwrtBSInTrack", "std::vector<float>", &T_pfMuon_IPwrtBSInTrack);
    Tree->Branch("T_pfMuon_IP2DBiasedPV", "std::vector<float>", &T_pfMuon_IP2DBiasedPV);
    Tree->Branch("T_pfMuon_IP3DBiasedPV", "std::vector<float>", &T_pfMuon_IP3DBiasedPV);
    Tree->Branch("T_pfMuon_IP2DUnBiasedPV", "std::vector<float>", &T_pfMuon_IP2DUnBiasedPV);
    Tree->Branch("T_pfMuon_IP3DUnBiasedPV", "std::vector<float>", &T_pfMuon_IP3DUnBiasedPV);

    Tree->Branch("T_pfMuon_vz", "std::vector<float>", &T_pfMuon_vz);
    Tree->Branch("T_pfMuon_vy", "std::vector<float>", &T_pfMuon_vy);
    Tree->Branch("T_pfMuon_vx", "std::vector<float>", &T_pfMuon_vx);
    Tree->Branch("T_pfMuon_NValidHits", "std::vector<int>", &T_pfMuon_NValidHits);
    Tree->Branch("T_pfMuon_NValidPixelHitsInTrk", "std::vector<int>", &T_pfMuon_NValidPixelHitsInTrk);
    Tree->Branch("T_pfMuon_particleIso", "std::vector<float>", &T_pfMuon_particleIso);
    Tree->Branch("T_pfMuon_chargedHadronIso", "std::vector<float>", &T_pfMuon_chargedHadronIso);
    Tree->Branch("T_pfMuon_neutralHadronIso", "std::vector<float>", &T_pfMuon_neutralHadronIso);
    Tree->Branch("T_pfMuon_photonIso", "std::vector<float>", &T_pfMuon_photonIso);
    Tree->Branch("T_pfMuon_deltaPt", "std::vector<float>", &T_pfMuon_deltaPt);
    Tree->Branch("T_pfMuon_NumOfMatchedStations", "std::vector<int>", &T_pfMuon_NumOfMatchedStations);
    Tree->Branch("T_pfMuon_pfCharged", "std::vector<float>", &T_pfMuon_pfCharged);
    Tree->Branch("T_pfMuon_pfNeutral", "std::vector<float>", &T_pfMuon_pfNeutral);  
    Tree->Branch("T_pfMuon_pfPhoton", "std::vector<float>", &T_pfMuon_pfPhoton);
    Tree->Branch("T_pfMuon_smurfCharged", "std::vector<float>", &T_pfMuon_smurfCharged);
    Tree->Branch("T_pfMuon_smurfNeutral", "std::vector<float>", &T_pfMuon_smurfNeutral);
    Tree->Branch("T_pfMuon_smurfPhoton", "std::vector<float>", &T_pfMuon_smurfPhoton);
    Tree->Branch("T_pfMuon_smurfNoOverCharged", "std::vector<float>", &T_pfMuon_smurfNoOverCharged);
    Tree->Branch("T_pfMuon_smurfNoOverNeutral", "std::vector<float>", &T_pfMuon_smurfNoOverNeutral);
    Tree->Branch("T_pfMuon_smurfNoOverPhoton", "std::vector<float>", &T_pfMuon_smurfNoOverPhoton);
    Tree->Branch("T_pfMuon_muSmurfPF", "std::vector<float>", &T_pfMuon_muSmurfPF);
    Tree->Branch("T_pfMuon_sumPUPt", "std::vector<float>", &T_pfMuon_sumPUPt);
    Tree->Branch("T_pfMuon_sumPUPtR03", "std::vector<float>", &T_pfMuon_sumPUPtR03);
    Tree->Branch("T_pfMuon_NLayers", "std::vector<int>", &T_pfMuon_NLayers);
  */
  //Vertex
  Tree->Branch("T_Vertex_z", "std::vector<float>", &T_Vertex_z);
  Tree->Branch("T_Vertex_y", "std::vector<float>", &T_Vertex_y);
  Tree->Branch("T_Vertex_x", "std::vector<float>", &T_Vertex_x);
  Tree->Branch("T_Vertex_Chi2Prob", "std::vector<float>", &T_Vertex_Chi2Prob);
  Tree->Branch("T_Vertex_rho","std::vector<float>", &T_Vertex_rho);
  Tree->Branch("T_Vertex_ndof","std::vector<float>", &T_Vertex_ndof);
  Tree->Branch("T_Vertex_isFake","std::vector<bool>", &T_Vertex_isFake);
  Tree->Branch("T_Vertex_tracksSize","std::vector<int>", &T_Vertex_tracksSize);

  //Electrons  
  Tree->Branch("T_Elec_Eta", "std::vector<float>", &T_Elec_Eta);
  Tree->Branch("T_Elec_IPwrtAveBS", "std::vector<float>", &T_Elec_IPwrtAveBS);
  Tree->Branch("T_Elec_IPwrtPV", "std::vector<float>", &T_Elec_IPwrtPV);
  Tree->Branch("T_Elec_dzwrtPV", "std::vector<float>", &T_Elec_dzwrtPV);
  Tree->Branch("T_Elec_Px", "std::vector<float>", &T_Elec_Px);
  Tree->Branch("T_Elec_Py", "std::vector<float>", &T_Elec_Py);
  Tree->Branch("T_Elec_Pz", "std::vector<float>", &T_Elec_Pz);
  Tree->Branch("T_Elec_Pt", "std::vector<float>", &T_Elec_Pt);
  Tree->Branch("T_Elec_Energy", "std::vector<float>", &T_Elec_Energy);
  Tree->Branch("T_Elec_Charge", "std::vector<int>", &T_Elec_Charge);

  /*  Tree->Branch("T_Elec_SumIsoCalo", "std::vector<float>", &T_Elec_SumIsoCalo);
      Tree->Branch("T_Elec_SumIsoTrack", "std::vector<float>", &T_Elec_SumIsoTrack);*/
  Tree->Branch("T_Elec_nBrems", "std::vector<int>", &T_Elec_nBrems); 
  Tree->Branch("T_Elec_fBrem", "std::vector<float>", &T_Elec_fBrem);
  Tree->Branch("T_Elec_eSuperClusterOverP", "std::vector<float>", &T_Elec_eSuperClusterOverP);
  Tree->Branch("T_Elec_ecalEnergy", "std::vector<float>", &T_Elec_ecalEnergy);
  /*  Tree->Branch("T_Elec_dr03TkSumPt", "std::vector<float>", &T_Elec_dr03TkSumPt);
      Tree->Branch("T_Elec_dr03EcalSumEt", "std::vector<float>", &T_Elec_dr03EcalSumEt);
      Tree->Branch("T_Elec_dr03HcalSumEt", "std::vector<float>", &T_Elec_dr03HcalSumEt);
  
      Tree->Branch("T_Elec_ConvInfoDist", "std::vector<float>", &T_Elec_ConvInfoDist);
      Tree->Branch("T_Elec_ConvInfoDCot", "std::vector<float>", &T_Elec_ConvInfoDCot);*/
  Tree->Branch("T_Elec_isEB", "std::vector<bool>", &T_Elec_isEB);
  Tree->Branch("T_Elec_isEE", "std::vector<bool>", &T_Elec_isEE);
  Tree->Branch("T_Elec_MVA", "std::vector<float>", &T_Elec_MVA);
  //  Tree->Branch("T_Elec_simpleEleId95", "std::vector<float>", &T_Elec_simpleEleId95);
  //  Tree->Branch("T_Elec_simpleEleId90", "std::vector<float>", &T_Elec_simpleEleId90);
  //  Tree->Branch("T_Elec_simpleEleId85", "std::vector<float>", &T_Elec_simpleEleId85);
  Tree->Branch("T_Elec_simpleEleId80", "std::vector<float>", &T_Elec_simpleEleId80);
  //  Tree->Branch("T_Elec_simpleEleId70", "std::vector<float>", &T_Elec_simpleEleId70);
  //  Tree->Branch("T_Elec_simpleEleId60", "std::vector<float>", &T_Elec_simpleEleId60); 
  //  Tree->Branch("T_Elec_simpleEleId90cIso", "std::vector<float>", &T_Elec_simpleEleId90cIso);
  /*
    Tree->Branch("T_Elec_cicVeryLooseMC", "std::vector<float>", &T_Elec_cicVeryLooseMC);
    Tree->Branch("T_Elec_cicLooseMC", "std::vector<float>", &T_Elec_cicLooseMC);
    Tree->Branch("T_Elec_cicMediumMC", "std::vector<float>", &T_Elec_cicMediumMC);
    Tree->Branch("T_Elec_cicSuperTightMC", "std::vector<float>", &T_Elec_cicSuperTightMC);
    Tree->Branch("T_Elec_cicHyperTight1MC", "std::vector<float>", &T_Elec_cicHyperTight1MC);
    Tree->Branch("T_Elec_cicHyperTight2MC", "std::vector<float>", &T_Elec_cicHyperTight2MC);
    Tree->Branch("T_Elec_cicHyperTight3MC", "std::vector<float>", &T_Elec_cicHyperTight3MC);
    Tree->Branch("T_Elec_cicVeryLooseHWW", "std::vector<float>", &T_Elec_cicVeryLooseHWW);
    Tree->Branch("T_Elec_cicLooseHWW", "std::vector<float>", &T_Elec_cicLooseHWW);
    Tree->Branch("T_Elec_cicMediumHWW", "std::vector<float>", &T_Elec_cicMediumHWW);
    Tree->Branch("T_Elec_cicSuperTightHWW", "std::vector<float>", &T_Elec_cicSuperTightHWW);
    Tree->Branch("T_Elec_cicHyperTight1HWW", "std::vector<float>", &T_Elec_cicHyperTight1HWW);
    Tree->Branch("T_Elec_cicHyperTight2HWW", "std::vector<float>", &T_Elec_cicHyperTight2HWW);
    Tree->Branch("T_Elec_cicHyperTight3HWW", "std::vector<float>", &T_Elec_cicHyperTight3HWW);
    Tree->Branch("T_Elec_cicVeryLoose", "std::vector<float>", &T_Elec_cicVeryLoose);
    Tree->Branch("T_Elec_cicLoose", "std::vector<float>", &T_Elec_cicLoose);
    Tree->Branch("T_Elec_cicMedium", "std::vector<float>", &T_Elec_cicMedium);
    Tree->Branch("T_Elec_cicSuperTight", "std::vector<float>", &T_Elec_cicSuperTight);
    Tree->Branch("T_Elec_cicHyperTight1", "std::vector<float>", &T_Elec_cicHyperTight1);
    Tree->Branch("T_Elec_cicHyperTight2", "std::vector<float>", &T_Elec_cicHyperTight2);
    Tree->Branch("T_Elec_cicHyperTight3", "std::vector<float>", &T_Elec_cicHyperTight3);
  */
  /*  Tree->Branch("T_Elec_Loose", "std::vector<float>", &T_Elec_Loose);
      Tree->Branch("T_Elec_RobustLoose", "std::vector<float>", &T_Elec_RobustLoose);
      Tree->Branch("T_Elec_Tight", "std::vector<float>", &T_Elec_Tight);
      Tree->Branch("T_Elec_RobustTight", "std::vector<float>", &T_Elec_RobustTight);
      Tree->Branch("T_Elec_RobustHighEnergy", "std::vector<float>", &T_Elec_RobustHighEnergy);
      Tree->Branch("T_Elec_egammaIDLikelihood", "std::vector<float>", &T_Elec_egammaIDLikelihood ) ;*/
  Tree->Branch("T_Elec_chargedHadronIso", "std::vector<float>", &T_Elec_chargedHadronIso);
  Tree->Branch("T_Elec_neutralHadronIso", "std::vector<float>", &T_Elec_neutralHadronIso);
  Tree->Branch("T_Elec_photonIso", "std::vector<float>", &T_Elec_photonIso);
  Tree->Branch("T_Elec_puChargedHadronIso", "std::vector<float>", &T_Elec_puChargedHadronIso);
  /*   Tree->Branch("T_Elec_smurfCharged", "std::vector<float>", &T_Elec_smurfCharged);
       Tree->Branch("T_Elec_smurfNeutral", "std::vector<float>", &T_Elec_smurfNeutral);
       Tree->Branch("T_Elec_smurfPhoton", "std::vector<float>", &T_Elec_smurfPhoton);
       Tree->Branch("T_Elec_smurfNoOverCharged", "std::vector<float>", &T_Elec_smurfNoOverCharged);
       Tree->Branch("T_Elec_smurfNoOverNeutral", "std::vector<float>", &T_Elec_smurfNoOverNeutral);
       Tree->Branch("T_Elec_smurfNoOverPhoton", "std::vector<float>", &T_Elec_smurfNoOverPhoton);
       Tree->Branch("T_Elec_eleSmurfPF","std::vector<float>", &T_Elec_eleSmurfPF);
  */
  Tree->Branch("T_Elec_passConversionVeto","std::vector<bool>",&T_Elec_passConversionVeto);
    
  Tree->Branch("T_Elec_sigmaIetaIeta", "std::vector<float>", &T_Elec_sigmaIetaIeta);
  Tree->Branch("T_Elec_deltaPhiIn", "std::vector<float>", &T_Elec_deltaPhiIn);
  Tree->Branch("T_Elec_deltaEtaIn", "std::vector<float>", &T_Elec_deltaEtaIn);
  Tree->Branch("T_Elec_isEcalDriven", "std::vector<bool>", &T_Elec_isEcalDriven);
  Tree->Branch("T_Elec_HtoE", "std::vector<float>", &T_Elec_HtoE);
 
  /*  Tree->Branch("T_Elec_IP2DBiasedPV", "std::vector<float>", &T_Elec_IP2DBiasedPV);
      Tree->Branch("T_Elec_IP3DBiasedPV", "std::vector<float>", &T_Elec_IP3DBiasedPV);
      Tree->Branch("T_Elec_IP2DUnBiasedPV", "std::vector<float>", &T_Elec_IP2DUnBiasedPV);
      Tree->Branch("T_Elec_IP3DUnBiasedPV", "std::vector<float>", &T_Elec_IP3DUnBiasedPV);
      Tree->Branch("T_Elec_dxyPVBiasedPV", "std::vector<float>", &T_Elec_dxyPVBiasedPV);
      Tree->Branch("T_Elec_dzPVBiasedPV", "std::vector<float>", &T_Elec_dzPVBiasedPV);
      Tree->Branch("T_Elec_dxyPVUnBiasedPV", "std::vector<float>", &T_Elec_dxyPVUnBiasedPV);
      Tree->Branch("T_Elec_dzPVUnBiasedPV", "std::vector<float>", &T_Elec_dzPVUnBiasedPV);
      Tree->Branch("T_Elec_IP2DUnBiasedPVnoBS", "std::vector<float>", &T_Elec_IP2DUnBiasedPVnoBS);
      Tree->Branch("T_Elec_IP3DUnBiasedPVnoBS", "std::vector<float>", &T_Elec_IP3DUnBiasedPVnoBS);
      Tree->Branch("T_Elec_dxyPVUnBiasedPVnoBS", "std::vector<float>", &T_Elec_dxyPVUnBiasedPVnoBS);
      Tree->Branch("T_Elec_dzPVUnBiasedPVnoBS", "std::vector<float>", &T_Elec_dzPVUnBiasedPVnoBS);
  */
  Tree->Branch("T_Elec_vz", "std::vector<float>", &T_Elec_vz);
  Tree->Branch("T_Elec_vy", "std::vector<float>", &T_Elec_vy);
  Tree->Branch("T_Elec_vx", "std::vector<float>", &T_Elec_vx);
  Tree->Branch("T_Elec_nLost", "std::vector<int>", &T_Elec_nLost);
  Tree->Branch("T_Elec_nHits", "std::vector<int>", &T_Elec_nHits);
  Tree->Branch("T_Elec_SC_Et","std::vector<float>", &T_Elec_SC_Et);
  Tree->Branch("T_Elec_SC_Eta","std::vector<float>", &T_Elec_SC_Eta); 

  Tree->Branch("T_Elec_PFElecPt","std::vector<float>", &T_Elec_PFElecPt);
  Tree->Branch("T_Elec_PFElecPx","std::vector<float>", &T_Elec_PFElecPx);
  Tree->Branch("T_Elec_PFElecPy","std::vector<float>", &T_Elec_PFElecPy);
  Tree->Branch("T_Elec_PFElecPz","std::vector<float>", &T_Elec_PFElecPx);
  Tree->Branch("T_Elec_PFElecE","std::vector<float>", &T_Elec_PFElecE);
  Tree->Branch("T_Elec_isPF","std::vector<bool>", &T_Elec_isPF);

  /*  Tree->Branch("T_Elec_passTriggerElMu","std::vector<bool>",&T_Elec_passTriggerElMu);
      Tree->Branch("T_Elec_passTriggerDoubleEl","std::vector<bool>",&T_Elec_passTriggerDoubleEl);
  */
  /* 
  //pfElectrons
  Tree->Branch("T_pfElec_Px", "std::vector<float>", &T_pfElec_Px);
  Tree->Branch("T_pfElec_Py", "std::vector<float>", &T_pfElec_Py);
  Tree->Branch("T_pfElec_Pz", "std::vector<float>", &T_pfElec_Pz);
  Tree->Branch("T_pfElec_Pt", "std::vector<float>", &T_pfElec_Pt);
  Tree->Branch("T_pfElec_Energy", "std::vector<float>", &T_pfElec_Energy);
  Tree->Branch("T_pfElec_Charge", "std::vector<int>", &T_pfElec_Charge);
  Tree->Branch("T_pfElec_SumIsoCalo", "std::vector<float>", &T_pfElec_SumIsoCalo);
  Tree->Branch("T_pfElec_SumIsoTrack", "std::vector<float>", &T_pfElec_SumIsoTrack);
  Tree->Branch("T_pfElec_IP2DBiasedPV", "std::vector<float>", &T_pfElec_IP2DBiasedPV);
  Tree->Branch("T_pfElec_IP3DBiasedPV", "std::vector<float>", &T_pfElec_IP3DBiasedPV);
  Tree->Branch("T_pfElec_vz", "std::vector<float>", &T_pfElec_vz);
  Tree->Branch("T_pfElec_vy", "std::vector<float>", &T_pfElec_vy);
  Tree->Branch("T_pfElec_vx", "std::vector<float>", &T_pfElec_vx);
  Tree->Branch("T_pfElec_nBrems", "std::vector<int>", &T_pfElec_nBrems);
  Tree->Branch("T_pfElec_dr03TkSumPt", "std::vector<float>", &T_pfElec_dr03TkSumPt);
  Tree->Branch("T_pfElec_dr03EcalSumEt", "std::vector<float>", &T_pfElec_dr03EcalSumEt);
  Tree->Branch("T_pfElec_dr03HcalSumEt", "std::vector<float>", &T_pfElec_dr03HcalSumEt);
  Tree->Branch("T_pfElec_SC_Et", "std::vector<float>", &T_pfElec_SC_Et);
  Tree->Branch("T_pfElec_SC_Eta", "std::vector<float>", &T_pfElec_SC_Eta);
  Tree->Branch("T_pfElec_nHits", "std::vector<int>", &T_pfElec_nHits);
  Tree->Branch("T_pfElec_ConvInfoDCot", "std::vector<float>", &T_pfElec_ConvInfoDCot);
  Tree->Branch("T_pfElec_ConvInfoDist", "std::vector<float>", &T_pfElec_ConvInfoDist);
  Tree->Branch("T_pfElec_isEcalDriven", "std::vector<bool>", &T_pfElec_isEcalDriven);
  Tree->Branch("T_pfElec_HtoE", "std::vector<float>", &T_pfElec_HtoE);
  Tree->Branch("T_pfElec_particleIso", "std::vector<float>", &T_pfElec_particleIso);
  Tree->Branch("T_pfElec_chargedHadronIso", "std::vector<float>", &T_pfElec_chargedHadronIso);
  Tree->Branch("T_pfElec_neutralHadronIso", "std::vector<float>", &T_pfElec_neutralHadronIso);
  Tree->Branch("T_pfElec_photonIso", "std::vector<float>", &T_pfElec_photonIso);
  */
  
  //Jets
  //  SetJetBranchAddress(0, "T_JetAK", true);
  SetJetBranchAddress(0, "T_JetAKCHS", true);
  //  SetJetBranchAddress(2, "T_JetAKPF2PAT", false);

  
  
  //MET 
  /*  Tree->Branch("T_MET_ET", &T_MET_ET, "T_MET_ET/F");
      Tree->Branch("T_MET_Phi", &T_MET_Phi, "T_MET_Phi/F");	
      Tree->Branch("T_MET_Sig", &T_MET_Sig, "T_MET_Sig/F");
  */   
  Tree->Branch("T_METPF_ET", &T_METPF_ET, "T_METPF_ET/F");
  Tree->Branch("T_METPF_Phi", &T_METPF_Phi, "T_METPF_Phi/F");	
  Tree->Branch("T_METPF_Sig", &T_METPF_Sig, "T_METPF_Sig/F");

  Tree->Branch("T_METPFTypeI_ET", &T_METPFTypeI_ET, "T_METPFTypeI_ET/F");
  Tree->Branch("T_METPFTypeI_Phi", &T_METPFTypeI_Phi, "T_METPFTypeI_Phi/F");
  Tree->Branch("T_METPFTypeI_Sig", &T_METPFTypeI_Sig, "T_METPFTypeI_Sig/F");

  /*  Tree->Branch("T_METChargedNeutralPFNoPU_ET", &T_METChargedNeutralPFNoPU_ET, "T_METChargedNeutralPFNoPU_ET/F"); 
      Tree->Branch("T_METChargedNeutralPFNoPU_Phi", &T_METChargedNeutralPFNoPU_Phi, "T_METChargedNeutralPFNoPU_Phi/F");

      Tree->Branch("T_METChargedPFNoPU_ET", &T_METChargedPFNoPU_ET, "T_METChargedPFNoPU_ET/F");
      Tree->Branch("T_METChargedPFNoPU_Phi", &T_METChargedPFNoPU_Phi, "T_METChargedPFNoPU_Phi/F");
  */
  /*  Tree->Branch("T_METtc_ET", &T_METtc_ET, "T_METtc_ET/F");
      Tree->Branch("T_METtc_Phi", &T_METtc_Phi, "T_METtc_Phi/F");	
      Tree->Branch("T_METtc_Sig", &T_METtc_Sig, "T_METtc_Sig/F");
  */
  /*  Tree->Branch("T_MET_assocPfMet_ET", &T_MET_assocPfMet_ET, "T_MET_assocPfMet_ET/F");
      Tree->Branch("T_MET_assocPfMet_Phi",&T_MET_assocPfMet_Phi, "T_MET_assocPfMet_Phi/F");

      Tree->Branch("T_MET_trkPfMet_ET", &T_MET_trkPfMet_ET, "T_MET_trkPfMet_ET/F");
      Tree->Branch("T_MET_trkPfMet_Phi", &T_MET_trkPfMet_Phi, "T_MET_trkPfMet_Phi/F");

      Tree->Branch("T_MET_assocOtherVtxPfMet_ET", &T_MET_assocOtherVtxPfMet_ET, "T_MET_assocOtherVtxPfMet_ET/F");
      Tree->Branch("T_MET_assocOtherVtxPfMet_Phi", &T_MET_assocOtherVtxPfMet_Phi, "T_MET_assocOtherVtxPfMet_Phi/F");

      Tree->Branch("T_MET_centralPfMet_ET", &T_MET_centralPfMet_ET, "T_MET_centralPfMet_ET/F");
      Tree->Branch("T_MET_centralPfMet_Phi",&T_MET_centralPfMet_Phi, "T_MET_centralPfMet_Phi/F");

      Tree->Branch("T_MET_cleanPfMet_ET", &T_MET_cleanPfMet_ET, "T_MET_cleanPfMet_ET/F");
      Tree->Branch("T_MET_cleanPfMet_Phi", &T_MET_cleanPfMet_Phi, "T_MET_cleanPfMet_Phi/F");
  */
  //  if(!IsRealData){  
  Tree->Branch("T_METgen_ET", &T_METgen_ET, "T_METgen_ET/F");
  Tree->Branch("T_METgen_Phi", &T_METgen_Phi, "T_METgen_Phi/F");             
  //  }
  Tree->Branch("T_passTriggerDoubleMu",&T_passTriggerDoubleMu,"T_passTriggerDoubleMu/O");
  Tree->Branch("T_passTriggerDoubleEl",&T_passTriggerDoubleEl,"T_passTriggerDoubleEl/O"); 
  //  Tree->Branch("T_passTriggerSingleMu",&T_passTriggerSingleMu,"T_passTriggerSingleMu/O");
  //  Tree->Branch("T_passTriggerSingleEl",&T_passTriggerSingleEl,"T_passTriggerSingleEl/O");
  Tree->Branch("T_passTriggerElMu"    ,&T_passTriggerElMu    ,"T_passTriggerElMu/O");
}
/*
  bool SUSYSkimToTree::passTriggerSingleMu(size_t i, bool isData) const{ 
  bool result(false);

  if( fabs(leps_[i].pdgId()) != 13 ) return false;

  const pat::Muon& mu = static_cast<const pat::Muon&>(leps_[i]);
  const pat::TriggerObjectStandAlone * match = mu.triggerObjectMatchByCollection("hltL3MuonCandidates");
  if(isData){
  if(match) result=match->hasPathName("HLT_Mu24_v*",false);}
  else{
  if(match) result=(match->hasPathName("HLT_Mu21_v*",false) && match->pt()>24.0);
  }

  return result;
  }

  bool SUSYSkimToTree::passTriggerDoubleMu(size_t i, bool isData) const{
  using namespace std;
  bool result(false);
  
  if( fabs(leps_[i].pdgId()) != 13 ) return false;
  
  const pat::Muon& mu = static_cast<const pat::Muon&>(leps_[i]);
  const pat::TriggerObjectStandAlone * match = mu.triggerObjectMatchByCollection("hltL3MuonCandidates");
  if(isData){
  if(match) result=match->hasPathName("HLT_DoubleMu7_v*",false);  }
  else{
  if(match) result=(match->hasPathName("HLT_DoubleMu5_v*",false) && match->pt()>7.0);  
  }
  return result;
  }

  bool SUSYSkimToTree::passTriggerDoubleEl(size_t i, bool isData) const{ 
  bool result(false);
  
  if( fabs(leps_[i].pdgId()) != 11 ) return false;
  const pat::Electron& el = static_cast<const pat::Electron&>(leps_[i]);
    
  if(isData){  
  if(el.triggerObjectMatchesByPath("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*").size() ||
  el.triggerObjectMatchesByPath("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*").size() ) 
  result=true;
  }else{
  if(el.triggerObjectMatchesByPath("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*").size() )
  result=true;
  }
  
  return result;

  }
  bool SUSYSkimToTree::passTriggerElMu(size_t i, bool isData) const{ 
*/
/*
  vector< std::string > pathNames = match->pathNames(false);
  for(unsigned int i=0; i<pathNames.size();++i){
  cout << "match path name: " << pathNames[i] << endl;
  }
*/

/*  bool result(false);
    using namespace std;
    if( fabs(leps_[i].pdgId()) == 13 ) {
    const pat::Muon& mu = static_cast<const pat::Muon&>(leps_[i]);    
    const pat::TriggerObjectStandAlone * match = mu.triggerObjectMatchByCollection("hltL3MuonCandidates");
    
    if(match){
    if(isData){
    bool res1 = match->hasPathName("HLT_Mu8_Ele17_CaloIdL_v*",false);
    bool res2 = match->hasPathName("HLT_Mu17_Ele8_CaloIdL_v*",false);
    result=( res1 || res2);}
    else{
    bool res1 = (match->hasPathName("HLT_Mu5_Ele17_v*",false) && match->pt()>8.0);
    bool res2 = (match->hasPathName("HLT_Mu11_Ele8_v*",false) && match->pt()>17.0);
    result=( res1 || res2);
    }    
    }
    return result;
    }

    if( fabs(leps_[i].pdgId()) == 11 ) {
    const pat::Electron& el = static_cast<const pat::Electron&>(leps_[i]);
    if(isData){
    if(el.triggerObjectMatchesByPath("HLT_Mu8_Ele17_CaloIdL*").size() ||
    el.triggerObjectMatchesByPath("HLT_Mu17_Ele8_CaloIdL*").size() ) 
    result=true;}
    else{
    const pat::TriggerObjectStandAlone * match1=
    el.triggerObjectMatchByPath("HLT_Mu5_Ele17_v*",true);
    const pat::TriggerObjectStandAlone * match2=
    el.triggerObjectMatchByPath("HLT_Mu11_Ele8_v*",true);     
    result=( match1 || match2 );
    }        
    }

    return result;

    }
*/


// ------------ method called once each job just after ending the event loop  ------------
void 
SUSYSkimToTree::endJob() {
  
  theHistosFile->cd();
  Tree->Write();
  theHistosFile->Close();
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(SUSYSkimToTree);

