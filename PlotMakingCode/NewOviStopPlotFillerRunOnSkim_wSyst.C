
#include "TChain.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TTree.h"
#include "TString.h"
#include "TRegexp.h"
#include "TProfile.h"
#include <fstream>
#include <string>
#include <iostream>
//#include <exception>                                                                                                      
#include <sys/stat.h>

#include "TCut.h"
//#include "StopDict_ver2.h"
#include "../HeaderFiles/StopFillFunctions.h"
#include "../HeaderFiles/StopFunctionDefinitions_v2.h"

//#include "PileUpMC.h"
#include <vector>
#include <cmath>
#include <sstream>
#include <map>
//#ifdef __MAKECINT__
//#pragma link C++ class vector >float< ;
//#endif

typedef struct {
    float stopProdXsec;
    float stopProdXsecUncert;
} StopXSec;
/*
 typedef pair<HistogramT_1D, SampleT> histKey_1D;
 typedef pair<HistogramT_2D, SampleT> histKey_2D;
 */
//typedef pair<HistogramT, SampleT> histKey;
typedef std::pair<HistogramT, SampleT> histKey;
typedef std::map<histKey, TH1 *>      HMap_1D;
typedef std::map<histKey, TH2 *>      HMap_2D;
typedef std::map<histKey, TH3 *>      HMap_3D;
typedef std::map<SampleT, bool>       passCutMap;


StopXSec getCrossSectionStop(float);
float    getMT2(TLorentzVector, TLorentzVector, float, float);

using namespace std;
//different nVtx regions for the plots
int main( int argc, const char* argv[] ) {
    /////////Variable initializations/////////////////
    /////Organization Variables//////////////   
    TFile * inputWeightFile = new TFile("StopABCPUWeight_vers2.root");
    TH3F * weights = (TH3F*) inputWeightFile->Get("WHist");
    TH1F * OneDPUDist = (TH1F*) weights->ProjectionX("OneDWeight");
    TFile * inputPUFileMCOvi = new TFile(TString("OneDMCPURWNewOvi.root"));
    TFile * inputPUFileMCDESY = new TFile(TString("OneDMCPURWNewDESY.root"));
    TFile * inputPUFileMCOviToDESY = new TFile(TString("OneDMCPURW_OviToDESY.root"));

    //    TH1F * truePUDistMC = (TH1F*) inputPUFileMC->Get("MCPU");
    //    TFile * inputPUFileData = new TFile("RunABCPUDist.root");
    //    TH1F * truePUDistData = (TH1F*) inputPUFileData->Get("pileup");
    //    float intMCPU = truePUDistMC->
    TH1F * nVtxSFHist = (TH1F*) inputPUFileMCOvi->Get("nVtxSF_preRW");
    TH1F * nVtxSFHistOviToDESY = (TH1F*) inputPUFileMCOviToDESY->Get("normFrac");
    TFile * inputFileMT2TwoDeeHists = new TFile(TString("MT2UncEnTwoDeeHists.root"));
    TH2F * MT2llUncEnUpDelta2D = (TH2F *) inputFileMT2TwoDeeHists->Get("h_MT2llUncEnUp_Delta2D");
    TH2F * MT2llUncEnDownDelta2D = (TH2F *) inputFileMT2TwoDeeHists->Get("h_MT2llUncEnDown_Delta2D");
    TH2F * MT2lbUncEnUpDelta2D = (TH2F *) inputFileMT2TwoDeeHists->Get("h_MT2lbUncEnUp_Delta2D");
    TH2F * MT2lbUncEnDownDelta2D = (TH2F *) inputFileMT2TwoDeeHists->Get("h_MT2lbUncEnDown_Delta2D");
    vector<TH1F *> * vecOneDeeMT2llUncEnUp = OneDProjectionReturnVec(MT2llUncEnUpDelta2D, 2, 2, 1, 1, 2, "MT2llUncEnUp");
    vector<TH1F *> * vecOneDeeMT2llUncEnDown = OneDProjectionReturnVec(MT2llUncEnDownDelta2D, 2, 2, 1, 1, 2, "MT2llUncEnDown");
    vector<TH1F *> * vecOneDeeMT2lbUncEnUp = OneDProjectionReturnVec(MT2lbUncEnUpDelta2D, 2, 2, 1, 1, 2, "MT2lbUncEnUp");
    vector<TH1F *> * vecOneDeeMT2lbUncEnDown = OneDProjectionReturnVec(MT2lbUncEnDownDelta2D, 2, 2, 1, 1, 2, "MT2lbUncEnDown");
    TString fileTreeName;    
    TString fInName;
    TString fOutName;    
    TH1F * h_eventCount, * h_CutFlow, * h_ElecCharIso, * h_ElecNeutIso, * h_ElecPhotIso, * h_ElecRelIso;
    h_ElecRelIso = new TH1F("h_ElecRelIso", "; Electron Relative Iso.; N_{evts} / bin", 32, 0., 0.16);
    
    /////Event Variables/////////////////////    
    int   Event, Type, nVtx, nVtxTrue;
    bool  passTrigDoubleMu, passTrigDoubleEl, passTrigElMu;
    bool  doEvent;
    float PUWeight;
    float weight, preNVtxRWweight;
    float fillWeight;
    float MET,METPhi,METSig, METPhi_preCorr;
    float METX, METY, METX_preCorr, METY_preCorr;
    float MT2ll, MT2lb; //, MT2lbPair1, MT2lbPair2;
    int caseMT2lb;
    vector<TLorentzVector> vecLepMT2lb(2), vecJetMT2lb(2);
    vector<TLorentzVector> vecBLepsMT2lb(2);
    float DeltaPhiMT2lb_JetsUsed, DeltaPhiMT2lb_BLepsUsed;
    float METdivMeff;//, METdivMeff_PassMT2llCut80, METdivMeff_PassMT2llCut90, METdivMeff_PassMT2llCut100, METdivMeff_PassMT2llCut110, METdivMeff_PassMT2llCut120;
    float MT2llCut = 80;
    float MT2lbCut = 172;
    //    float MT2llIndCut[5] = {80., 90., 100., 110., 120.};
    
    float Lep0Px,Lep0Py,Lep0Pz,Lep0E,Lep1Px,Lep1Py,Lep1Pz,Lep1E;
    int   Lep0PdgId, Lep1PdgId;
    float Lep0RelPFIso, Lep1RelPFIso;
    
    float HT;
    float Jet0Px,Jet0Py,Jet0Pz,Jet0E,Jet1Px,Jet1Py,Jet1Pz,Jet1E;
    float BtagJet0Px,BtagJet0Py,BtagJet0Pz,BtagJet0E,BtagJet1Px,BtagJet1Py,BtagJet1Pz,BtagJet1E;
    int   BtagJet0Index, BtagJet1Index;
    
    
    int TGenStopMass0,TGenStopMass1,TGenChi0Mass0,TGenChi0Mass1, TGenCharginoMass0, TGenCharginoMass1;
    int   grabStopMass, grabChi0Mass, grabCharginoMass;
    int   massDiffThresh = 5;
    float stopWeight    = 0.;
    float stopWeightErr = 0.;
    float stopWeightPlusErr = 0.;
    float stopWeightMinusErr = 0.;
    
    
    const double PI = 3.14159265;
    
    int NJets, NBtagJets;
    
    TLorentzVector Lep0Vec, Lep1Vec, DiLepVec;
    TLorentzVector Jet0Vec, Jet1Vec, DiJetVec, BtagJet0Vec, BtagJet1Vec, DiBJetVec;
    TLorentzVector Lep0Jet0Vec, Lep0Jet1Vec, Lep1Jet0Vec, Lep1Jet1Vec;
    TLorentzVector BLep0Vec, BLep1Vec;
    float diLepInvMass, diLepPt, diLepEta, diLepPhi;
    float diJetInvMass, diJetPt, diJetEta, diJetPhi;
    float diBJetInvMass, diBJetPt, diBJetEta, diBJetPhi;
    
//    LV * Lep0Vec_DESY, * Lep1Vec_DESY;
//    LV * Jet0Vec_DESY, * Jet1Vec_DESY, * BJet0Vec_DESY, * BJet1Vec_DESY;
    
    /// Generator Level Information
//    LV * genTopVec_DESY, * genAntiTopVec_DESY;
    
    bool  hasTopInfo = 0;
    float genTop0Pt, genTop0En, genTop0Eta, genTop0Phi;
    float genTop1Pt, genTop1En, genTop1Eta, genTop1Phi;
    int genTop0PdgId, genTop1PdgId; //genTop0FirstMom, genTop1FirstMom;
    float genTopPt; //, genTopEn, genTopEta, genTopPhi;
    float genAntiTopPt; //, genAntiTopEn, genAntiTopEta, genAntiTopPhi;
    float genMET_Pt, genMET_Phi;
    genTopPt = -1;
    genAntiTopPt = -1;
    
    TH1F * h_genTopPt = new TH1F("h_genTopPt", ";gen-level Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genTopPt->Sumw2();
    TH1F * h_genTopPtRW = new TH1F("h_genTopPtRW", ";gen-level Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genTopPtRW->Sumw2();
    TH1F * h_genAntiTopPt = new TH1F("h_genAntiTopPt", ";gen-level Anti-Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genAntiTopPt->Sumw2();
    TH1F * h_genAntiTopPtRW = new TH1F("h_genAntiTopPtRW", ";gen-level Anti-Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genAntiTopPtRW->Sumw2();
    
    /********************************************************************************************************/
    // Systematics Stuff
    int   runNumber, lumiBlock, eventNumber;
    TString eventCountName;
    
    int   Type_LepESUp, Type_LepESDown;
    bool  doEvent_LepESUp, doEvent_LepESDown;
    float weight_LepESUp, weight_LepESDown;
    float weight_LepEffSFUp, weight_LepEffSFDown;
    //    float weight_LepEffSFUp_LepESUp, weight_LepEffSFDown_LepESUp;
    //    float weight_LepEffSFUp_LepESDown, weight_LepEffSFDown_LepESDown;
    float preNVtxRWweight_LepESUp, preNVtxRWweight_LepESDown;
    float preNVtxRWweight_LepEffSFUp, preNVtxRWweight_LepEffSFDown;
    //    float preNVtxRWweight_LepEffSFUp_LepESUp, preNVtxRWweight_LepEffSFDown_LepESUp;
    //    float preNVtxRWweight_LepEffSFUp_LepESDown, preNVtxRWweight_LepEffSFDown_LepESDown;
    float weight_GenTopReweight, preNVtxRWweight_GenTopReweight;
    float weightSwitch;
    float weight_genStopXSecUp, weight_genStopXSecDown;
    float preNVtxRWweight_genStopXSecUp, preNVtxRWweight_genStopXSecDown;
    
    //    vector<int> * eventJetParams_JetESUp = new vector<int>;
    //    int jet0Index_JetESUp, jet1Index_JetESUp;
    
    //    vector<int> * eventJetParams_JetESDown = new vector<int>;
    //    int jet0Index_JetESDown, jet1Index_JetESDown;
    
    
    float MT2ll_ShiftUp; //, MT2ll_ShiftDown; // TTBar MT2ll smearing
    float MT2ll_UncESUp, MT2ll_UncESDown;
    float MT2ll_LepESUp, MT2ll_LepESDown;
    float MT2ll_JetESUp, MT2ll_JetESDown;
    float MET_LepESUp, MET_LepESDown, MET_JetESUp, MET_JetESDown;//, MET_JetERUp, MET_JetERDown;
    float METPhi_LepESUp, METPhi_LepESDown, METPhi_JetESUp, METPhi_JetESDown;//, METPhi_JetERUp, METPhi_JetERDown;
    float METPhi_preCorr_LepESUp, METPhi_preCorr_LepESDown, METPhi_preCorr_JetESUp, METPhi_preCorr_JetESDown;//, METPhi_preCorr_JetERUp, METPhi_preCorr_JetERDown;
    float METX_LepESUp, METY_LepESUp, METX_LepESDown, METY_LepESDown, METX_JetESUp, METY_JetESUp, METX_JetESDown, METY_JetESDown;//, METX_JetERUp, METY_JetERUp, METX_JetERDown, METY_JetERDown;
    float METX_preCorr_LepESUp, METY_preCorr_LepESUp, METX_preCorr_LepESDown, METY_preCorr_LepESDown, METX_preCorr_JetESUp, METY_preCorr_JetESUp, METX_preCorr_JetESDown, METY_preCorr_JetESDown;//, METX_preCorr_JetERUp, METY_preCorr_JetERUp, METX_preCorr_JetERDown, METY_preCorr_JetERDown;
    
    vector<TLorentzVector> vecLepMT2lb_LepESUp(2), vecLepMT2lb_LepESDown(2);
    vector<TLorentzVector> vecJetMT2lb_JetESUp(2), vecJetMT2lb_JetESDown(2);
    vector<TLorentzVector> vecBLepsMT2lb_LepESUp(2), vecBLepsMT2lb_LepESDown(2);
    vector<TLorentzVector> vecBLepsMT2lb_JetESUp(2), vecBLepsMT2lb_JetESDown(2);
    int caseMT2lb_JetESUp, caseMT2lb_JetESDown;
//    float MT2lbPair1_LepESUp, MT2lbPair2_LepESUp, MT2lbPair1_LepESDown, MT2lbPair2_LepESDown;
//    float MT2lbPair1_JetESUp, MT2lbPair2_JetESUp, MT2lbPair1_JetESDown, MT2lbPair2_JetESDown;
    float MT2lb_LepESUp, MT2lb_LepESDown, MT2lb_JetESUp, MT2lb_JetESDown, MT2lb_UncESUp, MT2lb_UncESDown;//, MT2lb_JetERUp, MT2lb_JetERDown;
    float DeltaPhiMT2lb_BLepsUsed_LepESUp, DeltaPhiMT2lb_BLepsUsed_LepESDown;
    float DeltaPhiMT2lb_JetsUsed_JetESUp, DeltaPhiMT2lb_JetsUsed_JetESDown, DeltaPhiMT2lb_BLepsUsed_JetESUp, DeltaPhiMT2lb_BLepsUsed_JetESDown;
    TLorentzVector BLep0Vec_LepESUp, BLep1Vec_LepESDown, BLep0Vec_JetESUp, BLep1Vec_JetESDown;
    TFile * MT2llSmearFile = new TFile("MT2llSmear.root");
    TH1F * MT2llMeanSmear = (TH1F *) MT2llSmearFile->Get("MT2llSmear");
    float MT2llSmearFactor;
    int MT2llSmearFactorBin;
//    float MT2llSystConst = 1.5;
//    float MT2llSystSlope = 0.0;
    
    float METdivMeff_LepESUp, METdivMeff_LepESDown;
    //    float METdivMeff_PassMT2llCut80_LepESUp, METdivMeff_PassMT2llCut90_LepESUp, METdivMeff_PassMT2llCut100_LepESUp, METdivMeff_PassMT2llCut110_LepESUp, METdivMeff_PassMT2llCut120_LepESUp;
    //    float METdivMeff_PassMT2llCut80_LepESDown, METdivMeff_PassMT2llCut90_LepESDown, METdivMeff_PassMT2llCut100_LepESDown, METdivMeff_PassMT2llCut110_LepESDown, METdivMeff_PassMT2llCut120_LepESDown;
    
    float METdivMeff_JetESUp, METdivMeff_JetESDown;
    //    float METdivMeff_PassMT2llCut80_JetESUp, METdivMeff_PassMT2llCut90_JetESUp, METdivMeff_PassMT2llCut100_JetESUp, METdivMeff_PassMT2llCut110_JetESUp, METdivMeff_PassMT2llCut120_JetESUp;
    //    float METdivMeff_PassMT2llCut80_JetESDown, METdivMeff_PassMT2llCut90_JetESDown, METdivMeff_PassMT2llCut100_JetESDown, METdivMeff_PassMT2llCut110_JetESDown, METdivMeff_PassMT2llCut120_JetESDown;
    
    float Lep0Px_LepESUp, Lep0Py_LepESUp, Lep0Pz_LepESUp, Lep0E_LepESUp, Lep1Px_LepESUp, Lep1Py_LepESUp, Lep1Pz_LepESUp, Lep1E_LepESUp;
    float Lep0Px_LepESDown, Lep0Py_LepESDown, Lep0Pz_LepESDown, Lep0E_LepESDown, Lep1Px_LepESDown, Lep1Py_LepESDown, Lep1Pz_LepESDown, Lep1E_LepESDown;
    int   Lep0PdgId_LepESUp, Lep1PdgId_LepESUp;
    int   Lep0PdgId_LepESDown, Lep1PdgId_LepESDown;
    float Lep0RelPFIso_LepESUp, Lep1RelPFIso_LepESUp;
    float Lep0RelPFIso_LepESDown, Lep1RelPFIso_LepESDown;
    
    
    int NJets_JetESUp, NBtagJets_JetESUp;
    int NJets_JetESDown, NBtagJets_JetESDown;
    float HT_JetESUp;
    float HT_JetESDown;
    float Jet0Px_JetESUp, Jet0Py_JetESUp, Jet0Pz_JetESUp, Jet0E_JetESUp, Jet1Px_JetESUp, Jet1Py_JetESUp, Jet1Pz_JetESUp, Jet1E_JetESUp;
    float Jet0Px_JetESDown, Jet0Py_JetESDown, Jet0Pz_JetESDown, Jet0E_JetESDown, Jet1Px_JetESDown, Jet1Py_JetESDown, Jet1Pz_JetESDown, Jet1E_JetESDown;
    float BtagJet0Px_JetESUp, BtagJet0Py_JetESUp, BtagJet0Pz_JetESUp, BtagJet0E_JetESUp, BtagJet1Px_JetESUp, BtagJet1Py_JetESUp, BtagJet1Pz_JetESUp, BtagJet1E_JetESUp;
    float BtagJet0Px_JetESDown, BtagJet0Py_JetESDown, BtagJet0Pz_JetESDown, BtagJet0E_JetESDown, BtagJet1Px_JetESDown, BtagJet1Py_JetESDown, BtagJet1Pz_JetESDown, BtagJet1E_JetESDown;
    int   BtagJet0Index_JetESUp,  BtagJet1Index_JetESUp;
    int   BtagJet0Index_JetESDown,  BtagJet1Index_JetESDown;
    
    TLorentzVector Lep0Vec_LepESUp, Lep0Vec_LepESDown, Lep0Vec_LepERUp, Lep0Vec_LepERDown;
    TLorentzVector Lep1Vec_LepESUp, Lep1Vec_LepESDown, Lep1Vec_LepERUp, Lep1Vec_LepERDown;
    TLorentzVector DiLepVec_LepESUp, DiLepVec_LepESDown;//, DiLepVec_LepERUp, DiLepVec_LepERDown;
    float diLepInvMass_LepESUp, diLepInvMass_LepESDown;//, diLepInvMass_LepERUp, diLepInvMass_LepERDown;
    float diLepPt_LepESUp, diLepPt_LepESDown;//, diLepPt_LepERUp, diLepPt_LepERDown;
    float diLepEta_LepESUp, diLepEta_LepESDown;//, diLepEta_LepERUp, diLepEta_LepERDown;
    float diLepPhi_LepESUp, diLepPhi_LepESDown;//, diLepPhi_LepERUp, diLepPhi_LepERDown;
    
    TLorentzVector Jet0Vec_JetESUp, Jet0Vec_JetESDown;//, Jet0Vec_JetERUp, Jet0Vec_JetERDown;
    TLorentzVector Jet1Vec_JetESUp, Jet1Vec_JetESDown;//, Jet1Vec_JetERUp, Jet1Vec_JetERDown;
    TLorentzVector DiJetVec_JetESUp, DiJetVec_JetESDown;//, DiJetVec_JetERUp, DiJetVec_JetERDown;
    float diJetInvMass_JetESUp, diJetInvMass_JetESDown;//, diJetInvMass_JetERUp, diJetInvMass_JetERDown;
    float diJetPt_JetESUp, diJetPt_JetESDown;
    float diJetEta_JetESUp, diJetEta_JetESDown;
    float diJetPhi_JetESUp, diJetPhi_JetESDown;
    
    TLorentzVector BtagJet0Vec_JetESUp, BtagJet0Vec_JetESDown;//, BtagJet0Vec_JetERUp, BtagJet0Vec_JetERDown;
    TLorentzVector BtagJet1Vec_JetESUp, BtagJet1Vec_JetESDown;//, BtagJet1Vec_JetERUp, BtagJet1Vec_JetERDown;
    TLorentzVector DiBJetVec_JetESUp, DiBJetVec_JetESDown;//, DiBJetVec_JetERUp, DiBJetVec_JetERDown;
    float diBJetInvMass_JetESUp, diBJetInvMass_JetESDown;//, diBJetInvMass_JetERUp, diBJetInvMass_JetERDown;
    float diBJetPt_JetESUp, diBJetPt_JetESDown;
    float diBJetEta_JetESUp, diBJetEta_JetESDown;
    float diBJetPhi_JetESUp, diBJetPhi_JetESDown;
    
    /********************************************************************************************************/
    
    
//    float lumi = 19300; //5296.3; // ipb                                                                                     
    const float genStopMassMin = 295, genStopMassMax = 355, genDeltaM_stopChi0_Min = 195, genDeltaM_stopChi0_Max = 255; 
    // add 5 GeV safety margin (deltaM = 10 GeV in the FineBin sample)  
//    float Nevt_stop_oneMassPoint = 50000 * ( (genStopMassMax-genStopMassMin)/10. ) * ( (genDeltaM_stopChi0_Max-genDeltaM_stopChi0_Min)/10. );  
    // 50k evts per point x Npoints
    
    ////input cuts/commands    
    
    bool grabOutDir      = 0;      // whether or not to use the file: "outputSavePath.txt" for where to save output
    bool doPhiCorr       = 1;      // whether to do the MetPhi asymmetry correction -- 6/25/13 as of right now parameters need to be updated
    bool doData          = 0;      // Whether you're running on data or not
    bool doVerbosity     = 0;      // prints a lot of debug info if turned on
    int  whichNTupleType = 0;      // 0 IFCA Oviedo; 1 DESY -- 6/25/13 as of right now, doesn't work with Oviedo
    bool doPURW          = 0;      // run pile up reweighting
    bool doHackPURW      = 0;      // called a "hack" because it was bin by bin reweighting to get the nVtx distribution in the "inclusive" channel to exactly match between data and MC
    bool doPURWOviToDESY = 0;      // exactly like doHackPURW but for Oviedo and DESY to enable comparisons across nTuples
    bool doBookSyst      = 0;      // used for deciding whether or not to book systematics
    bool doMETSmear      = false;  // early attempt to mimic jet smearing using gaussians -- didn't work, leave off
    double METSF         = 0.0;
    
    bool blindData       = 1;      // "Until further notice, leave ON -- cuts MT2ll > 80 is data in full object/event selection
    int  nEvents         = -1;     // limits total number of events one is running on
    int subLepPtCut      = 10;    // Sets the pT cut used for subLepPtCut
    
    bool doParallel      = 0;      // Whether or not to run parallel jobs for a given input file. If done, you have to do an additional hadding step
    int  startPointNum   = 1;
    int  numBreakPoints  = 0;
    bool isSignal        = 0;      // Whether or not one is running on signal -- if this is turned on the next two arguments have to be the stop mass to grab and the chi0 mass to grab respectively
    int  startPoint, endPoint;
    ////input cuts/commands    
    /////loop over inputs    
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"-i",2) == 0) {
            fInName = TString(argv[k+1]);   
        }
        else if (strncmp (argv[k],"noPhiCorr",9) == 0) {
            doPhiCorr = 0;   
        }
        else if (strncmp (argv[k],"METSmear",9) == 0) {
            doMETSmear = true;
            METSF = strtod(argv[k+1], NULL);
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
        else if (strncmp (argv[k],"doBookSyst", 10) == 0) {
            doBookSyst = 1;
        }
        else if (strncmp (argv[k],"gOutDir", 7) == 0) {
            grabOutDir = 1;
        }
        else if (strncmp (argv[k],"isSig", 5) == 0) {
            isSignal = 1;
            grabStopMass = strtol(argv[k+1], NULL, 10);
            grabChi0Mass = strtol(argv[k+2], NULL, 10);
            grabCharginoMass = strtol(argv[k+3], NULL, 10);
            cout << "Looking for StopMass: " << grabStopMass << ", Chi0Mass: " << grabChi0Mass << ", and CharginoMass: " << grabCharginoMass << endl;
        }
        else if (strncmp (argv[k],"ReleaseTheKraken", 16) == 0) {
            blindData = 0;
            cout << "RELEASING THE KRAKEN!!! " << endl;
            cout << "http://www.youtube.com/watch?v=gb2zIR2rvRQ " << endl;
        }
        else if (strncmp (argv[k],"limStats",8) == 0) {
            nEvents = strtol(argv[k+1], NULL, 10);   
        }        
        else if (strncmp (argv[k],"subLepPt",8) == 0) {
            subLepPtCut = strtol(argv[k+1], NULL, 10);   
        }
        else if (strncmp (argv[k],"doParallel",10) == 0) {
            doParallel = 1;
            cout << "Make sure you input two numbers after 'doParallel': First is number of break points -- second is which to start on" << endl;
            numBreakPoints = strtol(argv[k+1], NULL, 10 );
            startPointNum = strtol(argv[k+2], NULL, 10 );
        }
    }
    ////input cuts/commands    
    
    char Buffer[500];
    ifstream * outDirFile;
    TRegexp fCutSlash("[^/]+$");
    fOutName = "";
    TString outputPathName = (isSignal) ? "signalOutputSavePath.txt" : "outputSavePath.txt";
    if (grabOutDir) {
        outDirFile = new ifstream(TString(outputPathName));
        if (!(outDirFile->eof())) {
            outDirFile->getline(Buffer,500);
            fOutName += TString(string(Buffer));
            fOutName += "/"; //in case user forgot a slash
        }
    }
    fOutName += fInName(fCutSlash);
    if (fInName.Contains("MuEG") || fInName.Contains("DoubleMu") || fInName.Contains("DoubleEl") || fInName.Contains("run2012")) {
        cout << "Running on Data" << endl;
        doData = 1;
    }
    
    if (doData && blindData) {
        fOutName += "MT2Leq80";   
    }
    else if (doData && !blindData) {
        fOutName += "_NOTBLIND";
    }
    if (doData) doMETSmear = 0;
    if (doMETSmear) {
        fOutName += "METSmear";
        fOutName += METSF * 100.0;
    }
    if (whichNTupleType == 0) {
        fOutName += "_Oviedo";        
    }
    else {
        fOutName += "_DESY";
    }
    if (doData && doBookSyst) doBookSyst = 0;
    if (doPURW && !doData) fOutName += "_PURW";
    if (doPURWOviToDESY && !doData) fOutName += "OviToDESY";
    if (doBookSyst) fOutName += "_wSyst";
    if (isSignal) {
        fOutName += "_SignalStop";
        fOutName += grabStopMass;
        fOutName += "_Chi0";
        fOutName += grabChi0Mass;
        if (grabCharginoMass >= 0) {                        
            fOutName += "_Chargino";
            fOutName += grabCharginoMass;
        }
    }
    if (abs(subLepPtCut - 10.) > 1E-3) {
        fOutName += "_subLepPtCut";
        fOutName += subLepPtCut;
    }
    if (doParallel) {
        fOutName += "_Parallel_";
        fOutName += startPointNum;
        fOutName += ".root";
    }
    else {
        fOutName += "_Output.root";        
    }
    
    cout << "saving to " << fOutName << endl;
    TFile * outputFile;
    outputFile = new TFile(fOutName,"RECREATE");
    if (whichNTupleType == 1) {
        fileTreeName = "DESYSkimTree";
        nVtxSFHist = (TH1F*) inputPUFileMCDESY->Get("nVtxSF_preRW");
    }
    else {
        fileTreeName = "OviSkimTree";           
    }
    TChain fileTree(fileTreeName);
    TFile inputFile(fInName + TString(".root"));
    h_CutFlow = (TH1F *) inputFile.Get("h_CutFlow");
    bool  ZVeto           = true;
    bool  ZVeto_LepESUp   = true;
    bool  ZVeto_LepESDown   = true;
    float ZWindowLB      = 76;
    float ZWindowUB      = 106;   
    float barrelEtaEnd = 1.4442; float endcapEtaStart = 1.566;
    /////Set up the tree////////
    fileTree.Add(fInName + TString(".root"));
    
    if (whichNTupleType == 0) {
        fileTree.SetBranchAddress( "TPassDoubleMu",   &passTrigDoubleMu );
        fileTree.SetBranchAddress( "TPassDoubleEl",   &passTrigDoubleEl );
        fileTree.SetBranchAddress( "TPassElMu",       &passTrigElMu );
    }    
    fileTree.SetBranchAddress( "TWeight",   &PUWeight );
    fileTree.SetBranchAddress( "TChannel",  &Type );
    fileTree.SetBranchAddress( "TNPV",      &nVtx );
    fileTree.SetBranchAddress( "TNPV_True", &nVtxTrue );
    fileTree.SetBranchAddress( "TDoEvent",  &doEvent);
    
    fileTree.SetBranchAddress( "TMET",     &MET );
    fileTree.SetBranchAddress( "TMET_Phi", &METPhi );
    fileTree.SetBranchAddress( "TMETSig",  &METSig );
    
    fileTree.SetBranchAddress( "TLep0Px", &Lep0Px );
    fileTree.SetBranchAddress( "TLep0Py", &Lep0Py );
    fileTree.SetBranchAddress( "TLep0Pz", &Lep0Pz );
    fileTree.SetBranchAddress( "TLep0E",  &Lep0E );
    fileTree.SetBranchAddress( "TLep0PdgId", &Lep0PdgId );
    fileTree.SetBranchAddress( "TLep0RelPFIso",&Lep0RelPFIso);
    
    fileTree.SetBranchAddress( "TLep1Px", &Lep1Px );
    fileTree.SetBranchAddress( "TLep1Py", &Lep1Py );
    fileTree.SetBranchAddress( "TLep1Pz", &Lep1Pz );
    fileTree.SetBranchAddress( "TLep1E",  &Lep1E );
    fileTree.SetBranchAddress( "TLep1PdgId", &Lep1PdgId );
    fileTree.SetBranchAddress( "TLep1RelPFIso",&Lep1RelPFIso);
    
    fileTree.SetBranchAddress( "TNJets",    &NJets );
    //        fileTree.SetBranchAddress( "THT",       &HT );
    fileTree.SetBranchAddress( "TNJetsBtag",&NBtagJets );
    
    fileTree.SetBranchAddress( "HT", &HT );
    
    fileTree.SetBranchAddress( "TJet0Px", &Jet0Px );
    fileTree.SetBranchAddress( "TJet0Py", &Jet0Py );
    fileTree.SetBranchAddress( "TJet0Pz", &Jet0Pz );
    fileTree.SetBranchAddress( "TJet0E", &Jet0E );
    
    fileTree.SetBranchAddress( "TJet1Px", &Jet1Px );
    fileTree.SetBranchAddress( "TJet1Py", &Jet1Py );
    fileTree.SetBranchAddress( "TJet1Pz", &Jet1Pz );
    fileTree.SetBranchAddress( "TJet1E", &Jet1E );                                                                         
    
    fileTree.SetBranchAddress( "TBtagJet0Px", &BtagJet0Px );
    fileTree.SetBranchAddress( "TBtagJet0Py", &BtagJet0Py );
    fileTree.SetBranchAddress( "TBtagJet0Pz", &BtagJet0Pz );
    fileTree.SetBranchAddress( "TBtagJet0E", &BtagJet0E ); 
    fileTree.SetBranchAddress( "TBtagJet0Index", &BtagJet0Index );
    
    fileTree.SetBranchAddress( "TBtagJet1Px", &BtagJet1Px );
    fileTree.SetBranchAddress( "TBtagJet1Py", &BtagJet1Py );
    fileTree.SetBranchAddress( "TBtagJet1Pz", &BtagJet1Pz );
    fileTree.SetBranchAddress( "TBtagJet1E", &BtagJet1E );
    fileTree.SetBranchAddress( "TBtagJet1Index", &BtagJet1Index );
    
    fileTree.SetBranchAddress( "TRunNum", &runNumber );
    fileTree.SetBranchAddress( "TEventNum", &eventNumber );
    fileTree.SetBranchAddress( "TLumiBlock", &lumiBlock );
    if (!doData) {
        fileTree.SetBranchAddress( "TChannel_LepESUp", &Type_LepESUp );
        fileTree.SetBranchAddress( "TChannel_LepESDown", &Type_LepESDown );
        fileTree.SetBranchAddress( "TDoEvent_LepESUp", &doEvent_LepESUp );
        fileTree.SetBranchAddress( "TDoEvent_LepESDown", &doEvent_LepESDown );
        
        fileTree.SetBranchAddress( "TMET_LepESUp",     &MET_LepESUp );
        fileTree.SetBranchAddress( "TMET_LepESDown",     &MET_LepESDown );
        fileTree.SetBranchAddress( "TMET_Phi_LepESUp", &METPhi_LepESUp );
        fileTree.SetBranchAddress( "TMET_Phi_LepESDown", &METPhi_LepESDown );
        fileTree.SetBranchAddress( "TMET_JetESUp",     &MET_JetESUp );
        fileTree.SetBranchAddress( "TMET_JetESDown",     &MET_JetESDown );
        fileTree.SetBranchAddress( "TMET_Phi_JetESUp", &METPhi_JetESUp );
        fileTree.SetBranchAddress( "TMET_Phi_JetESDown", &METPhi_JetESDown );
        
        fileTree.SetBranchAddress( "TLep0Px_LepESUp", &Lep0Px_LepESUp );
        fileTree.SetBranchAddress( "TLep0Px_LepESDown", &Lep0Px_LepESDown );
        fileTree.SetBranchAddress( "TLep0Py_LepESUp", &Lep0Py_LepESUp );
        fileTree.SetBranchAddress( "TLep0Py_LepESDown", &Lep0Py_LepESDown );
        fileTree.SetBranchAddress( "TLep0Pz_LepESUp", &Lep0Pz_LepESUp );
        fileTree.SetBranchAddress( "TLep0Pz_LepESDown", &Lep0Pz_LepESDown );
        fileTree.SetBranchAddress( "TLep0E_LepESUp",  &Lep0E_LepESUp );
        fileTree.SetBranchAddress( "TLep0E_LepESDown",  &Lep0E_LepESDown );
        fileTree.SetBranchAddress( "TLep0PdgId_LepESUp", &Lep0PdgId_LepESUp );
        fileTree.SetBranchAddress( "TLep0PdgId_LepESDown", &Lep0PdgId_LepESDown );
        fileTree.SetBranchAddress( "TLep0RelPFIso_LepESUp",&Lep0RelPFIso_LepESUp );
        fileTree.SetBranchAddress( "TLep0RelPFIso_LepESDown",&Lep0RelPFIso_LepESDown );
        
        fileTree.SetBranchAddress( "TLep1Px_LepESUp", &Lep1Px_LepESUp );
        fileTree.SetBranchAddress( "TLep1Px_LepESDown", &Lep1Px_LepESDown );
        fileTree.SetBranchAddress( "TLep1Py_LepESUp", &Lep1Py_LepESUp );
        fileTree.SetBranchAddress( "TLep1Py_LepESDown", &Lep1Py_LepESDown );
        fileTree.SetBranchAddress( "TLep1Pz_LepESUp", &Lep1Pz_LepESUp );
        fileTree.SetBranchAddress( "TLep1Pz_LepESDown", &Lep1Pz_LepESDown );
        fileTree.SetBranchAddress( "TLep1E_LepESUp",  &Lep1E_LepESUp );
        fileTree.SetBranchAddress( "TLep1E_LepESDown",  &Lep1E_LepESDown );
        fileTree.SetBranchAddress( "TLep1PdgId_LepESUp", &Lep1PdgId_LepESUp );
        fileTree.SetBranchAddress( "TLep1PdgId_LepESDown", &Lep1PdgId_LepESDown );
        fileTree.SetBranchAddress( "TLep1RelPFIso_LepESUp",&Lep1RelPFIso_LepESUp );
        fileTree.SetBranchAddress( "TLep1RelPFIso_LepESDown",&Lep1RelPFIso_LepESDown );
        
        fileTree.SetBranchAddress( "TNJets_JetESUp",    &NJets_JetESUp );
        fileTree.SetBranchAddress( "TNJets_JetESDown",    &NJets_JetESDown );
        fileTree.SetBranchAddress( "TNJetsBtag_JetESUp",&NBtagJets_JetESUp );
        fileTree.SetBranchAddress( "TNJetsBtag_JetESDown",&NBtagJets_JetESDown );
        fileTree.SetBranchAddress( "HT_JetESUp", &HT_JetESUp );
        fileTree.SetBranchAddress( "HT_JetESDown", &HT_JetESDown );
        
        fileTree.SetBranchAddress( "TJet0Px_JetESUp", &Jet0Px_JetESUp );
        fileTree.SetBranchAddress( "TJet0Px_JetESDown", &Jet0Px_JetESDown );
        fileTree.SetBranchAddress( "TJet0Py_JetESUp", &Jet0Py_JetESUp );
        fileTree.SetBranchAddress( "TJet0Py_JetESDown", &Jet0Py_JetESDown );
        fileTree.SetBranchAddress( "TJet0Pz_JetESUp", &Jet0Pz_JetESUp );
        fileTree.SetBranchAddress( "TJet0Pz_JetESDown", &Jet0Pz_JetESDown );
        fileTree.SetBranchAddress( "TJet0E_JetESUp", &Jet0E_JetESUp );
        fileTree.SetBranchAddress( "TJet0E_JetESDown", &Jet0E_JetESDown );
        
        fileTree.SetBranchAddress( "TJet1Px_JetESUp", &Jet1Px_JetESUp );
        fileTree.SetBranchAddress( "TJet1Px_JetESDown", &Jet1Px_JetESDown );
        fileTree.SetBranchAddress( "TJet1Py_JetESUp", &Jet1Py_JetESUp );
        fileTree.SetBranchAddress( "TJet1Py_JetESDown", &Jet1Py_JetESDown );
        fileTree.SetBranchAddress( "TJet1Pz_JetESUp", &Jet1Pz_JetESUp );
        fileTree.SetBranchAddress( "TJet1Pz_JetESDown", &Jet1Pz_JetESDown );
        fileTree.SetBranchAddress( "TJet1E_JetESUp", &Jet1E_JetESUp );
        fileTree.SetBranchAddress( "TJet1E_JetESDown", &Jet1E_JetESDown );
        
        fileTree.SetBranchAddress( "TBtagJet0Px_JetESUp", &BtagJet0Px_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet0Px_JetESDown", &BtagJet0Px_JetESDown );
        fileTree.SetBranchAddress( "TBtagJet0Py_JetESUp", &BtagJet0Py_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet0Py_JetESDown", &BtagJet0Py_JetESDown );
        fileTree.SetBranchAddress( "TBtagJet0Pz_JetESUp", &BtagJet0Pz_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet0Pz_JetESDown", &BtagJet0Pz_JetESDown );
        fileTree.SetBranchAddress( "TBtagJet0E_JetESUp", &BtagJet0E_JetESUp ); 
        fileTree.SetBranchAddress( "TBtagJet0E_JetESDown", &BtagJet0E_JetESDown ); 
        fileTree.SetBranchAddress( "TBtagJet0Index_JetESUp", &BtagJet0Index_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet0Index_JetESDown", &BtagJet0Index_JetESDown );
        
        fileTree.SetBranchAddress( "TBtagJet1Px_JetESUp", &BtagJet1Px_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet1Px_JetESDown", &BtagJet1Px_JetESDown );
        fileTree.SetBranchAddress( "TBtagJet1Py_JetESUp", &BtagJet1Py_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet1Py_JetESDown", &BtagJet1Py_JetESDown );
        fileTree.SetBranchAddress( "TBtagJet1Pz_JetESUp", &BtagJet1Pz_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet1Pz_JetESDown", &BtagJet1Pz_JetESDown );
        fileTree.SetBranchAddress( "TBtagJet1E_JetESUp", &BtagJet1E_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet1E_JetESDown", &BtagJet1E_JetESDown );
        fileTree.SetBranchAddress( "TBtagJet1Index_JetESUp", &BtagJet1Index_JetESUp );
        fileTree.SetBranchAddress( "TBtagJet1Index_JetESDown", &BtagJet1Index_JetESDown );
    }        
    if (fileTree.GetBranch("TGenMET")) {
        fileTree.SetBranchAddress("TGenMET", &genMET_Pt );
        fileTree.SetBranchAddress("TGenMETPhi", &genMET_Phi );
    }
    else {
        genMET_Pt = -99.;
        genMET_Phi = -99.;
    }
    if (fileTree.GetBranch("TGenTopSt3_0_Pt")) {
        hasTopInfo = 1;
        fileTree.SetBranchAddress("TGenTopSt3_0_Energy", &genTop0En);
        fileTree.SetBranchAddress("TGenTopSt3_0_Pt",     &genTop0Pt);
        fileTree.SetBranchAddress("TGenTopSt3_0_Phi",    &genTop0Phi);
        fileTree.SetBranchAddress("TGenTopSt3_0_Eta",    &genTop0Eta);
        fileTree.SetBranchAddress("TGenTopSt3_0_PDGID",  &genTop0PdgId);
        fileTree.SetBranchAddress("TGenTopSt3_1_Energy", &genTop1En);
        fileTree.SetBranchAddress("TGenTopSt3_1_Pt",     &genTop1Pt);
        fileTree.SetBranchAddress("TGenTopSt3_1_Phi",    &genTop1Phi);
        fileTree.SetBranchAddress("TGenTopSt3_1_Eta",    &genTop1Eta);
        fileTree.SetBranchAddress("TGenTopSt3_1_PDGID",  &genTop1PdgId);
    }
    else {
        genTop0Pt = -99.;
        genTop1Pt = -99.;
        genTop0PdgId = -99;
        genTop1PdgId = -99;
    }
    if (fileTree.GetBranch("TGenStopMass0")){
        fileTree.SetBranchAddress( "TGenStopMass0", &TGenStopMass0 );
        fileTree.SetBranchAddress( "TGenStopMass1", &TGenStopMass1 );
        fileTree.SetBranchAddress( "TGenChi0Mass0", &TGenChi0Mass0 );
        fileTree.SetBranchAddress( "TGenChi0Mass1", &TGenChi0Mass1 );            
        fileTree.SetBranchAddress( "TGenCharginoMass0", &TGenCharginoMass0 );
        fileTree.SetBranchAddress( "TGenCharginoMass1", &TGenCharginoMass1 );
    }
    else{
        TGenStopMass0 = -99;
        TGenStopMass1 = -99;
        TGenChi0Mass0 = -99;
        TGenChi0Mass1 = -99;
        TGenCharginoMass0 = -99;
        TGenCharginoMass1 = -99;
    }
    if (fileTree.GetBranch("TGenStopMass0")){
        float stopMassToUseForXSec;
        if (!fInName.Contains("FineBin")) {
            float roundedStopMass = floor((TGenStopMass0+12.5)/25.)*25.; 
            stopMassToUseForXSec = roundedStopMass;
            //                cout << "roundedStopMass " << roundedStopMass << endl;
        }
        else {
            stopMassToUseForXSec = grabStopMass;
        }
        cout << "TEMP FIX !!! " << endl;
        stopMassToUseForXSec = grabStopMass;
        StopXSec theStopXSec = getCrossSectionStop(stopMassToUseForXSec);
        stopWeight = theStopXSec.stopProdXsec;            
        stopWeightErr = theStopXSec.stopProdXsecUncert * stopWeight;
        stopWeightPlusErr = stopWeight + stopWeightErr;
        stopWeightMinusErr = stopWeight - stopWeightErr;            
        cout << "stopWeight " << stopWeight << endl;
        cout << "stopWeightErr " << stopWeightErr << endl;
        cout << "stopWeightErr " << stopWeightPlusErr << endl;
        cout << "stopWeightErr " << stopWeightMinusErr << endl;
    }
    ////Book histograms and histogram names
    outputFile->cd();
    vector<HistogramT> * histVec_1D = OneDeeHistTVec();
    vector<HistogramT> * histVec_2D = TwoDeeHistTVec();
    vector<HistogramT> * histVec_3D = ThreeDeeHistTVec();
    vector<SampleT> * subSampVec    = SubSampVec();    ///Define things necessary for booking histograms
    vector<SystT> * systVec         = SystVec();
    vector<HistogramT> * histVec_1D_Syst = new vector<HistogramT>;
    vector<HistogramT> * histVec_2D_Syst = new vector<HistogramT>;
    vector<HistogramT> * histVec_3D_Syst = new vector<HistogramT>;
    //    vector<HistogramT> * histVec_2D_Syst;
    //    vector<HistogramT> * histVec_3D_Syst;
    HistogramT H_Current; 
    TH1F * h_1DCurr; TH2F * h_2DCurr; TH3F * h_3DCurr;
    HMap_1D histMap_1D; HMap_2D histMap_2D; HMap_3D histMap_3D; passCutMap subSampBool;
    passCutMap subSampBool_LepESUp, subSampBool_LepESDown;
    passCutMap subSampBool_JetESUp, subSampBool_JetESDown;
    SampleT S_Current;
    TString histTitle;
    TString axesTitle;
    float nXBins, nYBins, nZBins;
    float xBinMin, xBinMax;
    float yBinMin, yBinMax;
    float zBinMin, zBinMax;
    if (doBookSyst) {
        histVec_1D_Syst = AddSystHists(histVec_1D, systVec, fInName, isSignal);
        
         histVec_2D_Syst = AddSystHists(histVec_2D, systVec, fInName, isSignal);
        /*
         histVec_3D_Syst = AddSystHists(histVec_3D, systVec, fInName, isSignal);
         */
        //        addSystHists(histVec_2D, systVec);
        //        addSystHists(histVec_3D, systVec);
    }
    for (unsigned int i = 0; i < subSampVec->size(); ++i) {
        S_Current = subSampVec->at(i);
        for (int j = 0; j < (int) histVec_1D->size(); ++j) {            
            H_Current = histVec_1D->at(j);
            histTitle = H_Current.name + S_Current.histNameSuffix;
            if (doVerbosity) cout << "j = " << j << " and the hist title is " << histTitle << endl;
            axesTitle = ";"; axesTitle += H_Current.xLabel;// axesTitle += S_Current.histXaxisSuffix;
            axesTitle += ";"; axesTitle += H_Current.yLabel;
            if (S_Current.doZVeto == 0 && H_Current.xLabel.Contains("M_{ll}")) {
                nXBins = 100;
                xBinMin = ZWindowLB - 10;
                xBinMax = ZWindowUB + 10;
            }
            else {
                nXBins = H_Current.xBinN;
                xBinMin = H_Current.xMin;
                xBinMax = H_Current.xMax;
            }
            if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel && !H_Current.name.Contains("h_PassMT2ll") && !(H_Current.name.Contains("h_MT2llControl"))) {
                xBinMax = xBinMax / 2;
            }
            //            if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.
            h_1DCurr = new TH1F(histTitle, axesTitle, nXBins, xBinMin, xBinMax); h_1DCurr->Sumw2();
            
            if (H_Current.name.Contains("h_ChannelCutFlow")) {
                for (int i = 1; i < h_CutFlow->GetNbinsX()+1; ++i) {
                    h_1DCurr->SetBinContent(i, h_CutFlow->GetBinContent(i));
                }
            }
            
            histMap_1D[histKey(H_Current, S_Current)] = h_1DCurr;
            //will this cause memory leaks?
        }
        for (unsigned int k = 0; k < histVec_2D->size(); ++k) {
            H_Current = histVec_2D->at(k);
            histTitle = H_Current.name + S_Current.histNameSuffix;
            axesTitle = ";"; axesTitle += H_Current.xLabel;// axesTitle += S_Current.histXaxisSuffix;
            axesTitle += ";"; axesTitle += H_Current.yLabel;// axesTitle += S_Current.histYaxisSuffix;
            axesTitle += ";"; axesTitle += H_Current.zLabel;
            nXBins = H_Current.xBinN;
            xBinMin = H_Current.xMin;
            xBinMax = H_Current.xMax;
            nYBins = H_Current.yBinN;
            yBinMin = H_Current.yMin;
            yBinMax = H_Current.yMax;
            if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel && !H_Current.name.Contains("h_PassMT2ll") && !(H_Current.name.Contains("h_MT2llControl"))) {
                xBinMax = xBinMax / 2;
            }
            h_2DCurr = new TH2F(histTitle, axesTitle, nXBins, xBinMin, xBinMax, nYBins, yBinMin, yBinMax); h_2DCurr->Sumw2();
            histMap_2D[histKey(H_Current, S_Current)] = h_2DCurr;
        }        
        for (unsigned int l = 0; l < histVec_3D->size(); ++l) {
            H_Current = histVec_3D->at(l);
            histTitle = H_Current.name + S_Current.histNameSuffix;
            axesTitle = ";"; axesTitle += H_Current.xLabel;// axesTitle += S_Current.histXaxisSuffix;
            axesTitle += ";"; axesTitle += H_Current.yLabel;// axesTitle += S_Current.histYaxisSuffix;
            axesTitle += ";"; axesTitle += H_Current.zLabel;
            nXBins = H_Current.xBinN;
            xBinMin = H_Current.xMin;
            xBinMax = H_Current.xMax;
            nYBins = H_Current.yBinN;
            yBinMin = H_Current.yMin;
            yBinMax = H_Current.yMax;
            nZBins = H_Current.zBinN;
            zBinMin = H_Current.zMin;
            zBinMax = H_Current.zMax;
            if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel && !H_Current.name.Contains("h_PassMT2ll") && !(H_Current.name.Contains("h_MT2llControl"))) {
                xBinMax = xBinMax / 2;
            }
            h_3DCurr = new TH3F(histTitle, axesTitle, nXBins, xBinMin, xBinMax, nYBins, yBinMin, yBinMax, nZBins, zBinMin, zBinMax); h_3DCurr->Sumw2(); 
            histMap_3D[histKey(H_Current, S_Current)] = h_3DCurr;
        }
        if (doBookSyst) {
            for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {                
                H_Current = histVec_1D_Syst->at(js);
                histTitle = H_Current.name + S_Current.histNameSuffix;
                axesTitle = ";"; axesTitle += H_Current.xLabel;// axesTitle += S_Current.histXaxisSuffix;
                axesTitle += ";"; axesTitle += H_Current.yLabel;  
                if (S_Current.doZVeto == 0 && H_Current.xLabel.Contains("M_{ll}")) {
                    nXBins = 100;
                    xBinMin = ZWindowLB - 10;
                    xBinMax = ZWindowUB + 10;
                }
                else {
                    nXBins = H_Current.xBinN;
                    xBinMin = H_Current.xMin;
                    xBinMax = H_Current.xMax;
                }
                if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel && !H_Current.name.Contains("h_PassMT2ll") && !(H_Current.name.Contains("h_MT2llControl"))) {
                    xBinMax = xBinMax / 2;
                }
                h_1DCurr = new TH1F(histTitle, axesTitle, nXBins, xBinMin, xBinMax); h_1DCurr->Sumw2();
                histMap_1D[histKey(H_Current, S_Current)] = h_1DCurr;
                //will this cause memory leaks?
            }
             for (unsigned int ks = 0; ks < histVec_2D_Syst->size(); ++ks) {
             H_Current = histVec_2D_Syst->at(ks);
             histTitle = H_Current.name + S_Current.histNameSuffix;
             axesTitle = ";"; axesTitle += H_Current.xLabel;// axesTitle += S_Current.histXaxisSuffix;
             axesTitle += ";"; axesTitle += H_Current.yLabel;// axesTitle += S_Current.histYaxisSuffix;
             axesTitle += ";"; axesTitle += H_Current.zLabel;
             nXBins = H_Current.xBinN;
             xBinMin = H_Current.xMin;
             xBinMax = H_Current.xMax;
             nYBins = H_Current.yBinN;
             yBinMin = H_Current.yMin;
             yBinMax = H_Current.yMax;
             if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel && !H_Current.name.Contains("h_PassMT2ll") && !(H_Current.name.Contains("h_MT2llControl"))) {
             xBinMax = xBinMax / 2;
             }
             h_2DCurr = new TH2F(histTitle, axesTitle, nXBins, xBinMin, xBinMax, nYBins, yBinMin, yBinMax); h_2DCurr->Sumw2();
             histMap_2D[histKey(H_Current, S_Current)] = h_2DCurr;
             }
            /*
             for (unsigned int ls = 0; ls < histVec_3D_Syst->size(); ++ls) {
             H_Current = histVec_3D_Syst->at(ls);
             histTitle = H_Current.name + S_Current.histNameSuffix;
             axesTitle = ";"; axesTitle += H_Current.xLabel;// axesTitle += S_Current.histXaxisSuffix;
             axesTitle += ";"; axesTitle += H_Current.yLabel;// axesTitle += S_Current.histYaxisSuffix;
             axesTitle += ";"; axesTitle += H_Current.zLabel;
             nXBins = H_Current.xBinN;
             xBinMin = H_Current.xMin;
             xBinMax = H_Current.xMax;
             nYBins = H_Current.yBinN;
             yBinMin = H_Current.yMin;
             yBinMax = H_Current.yMax;
             nZBins = H_Current.zBinN;
             zBinMin = H_Current.zMin;
             zBinMax = H_Current.zMax;
             if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel && !H_Current.name.Contains("h_PassMT2ll") && !(H_Current.name.Contains("h_MT2llControl"))) {
             xBinMax = xBinMax / 2;
             }
             h_3DCurr = new TH3F(histTitle, axesTitle, nXBins, xBinMin, xBinMax, nYBins, yBinMin, yBinMax, nZBins, zBinMin, zBinMax); h_3DCurr->Sumw2(); 
             histMap_3D[histKey(H_Current, S_Current)] = h_3DCurr;
             }
             */
        }
    }   
    /////
    TRandom rand;
    map<string, float>::iterator xIter;
    map<string, float>::iterator yIter;
    map<string, float>::iterator zIter;    
    /////Iterate over events  
    cout << "--- Total Events in file: " << fileTree.GetEntries() << " events" << endl;
    if (nEvents < 0) {
        cout << "running over all events " << endl;  
        nEvents = fileTree.GetEntries();
    }
    else {
        if (nEvents > fileTree.GetEntries()) nEvents = fileTree.GetEntries();
        cout << "running on just " << nEvents << " events " << endl;
    }
    vector<int> * breakPoints = new vector<int>;
    int currBreakPoint;
    int endBreakPoint = -1;
    breakPoints->push_back(0);
    if (numBreakPoints > 0) {
        endBreakPoint = (numBreakPoints + 1 > nEvents) ? nEvents : numBreakPoints + 1;
        for (int i = 1; i < endBreakPoint; ++i) {
            currBreakPoint = i * nEvents / endBreakPoint;
            breakPoints->push_back(currBreakPoint);
        }
    }
    if (startPointNum > (int) breakPoints->size()) return 0;
    startPoint = breakPoints->at(startPointNum - 1);
    cout << " startpoint (event number 'startpoint'): " << startPoint << endl;
    endPoint = ((int) breakPoints->size() <= startPointNum) ? nEvents : breakPoints->at(startPointNum);
    cout << " endpoint (the 'endpoint + 1' th event won't be run over -- i.e. event number 'endpoint'): " << endPoint << endl;    
    if (whichNTupleType == 0) {
        eventCountName = "h_eventCount";
    }
    else {
        eventCountName = "weightedEvents";
    }
    h_eventCount = (TH1F *) inputFile.Get(eventCountName);
    h_ElecCharIso = (TH1F *) inputFile.Get("h_ElecCharIso");
    h_ElecNeutIso = (TH1F *) inputFile.Get("h_ElecNeutIso");
    h_ElecPhotIso = (TH1F *) inputFile.Get("h_ElecPhotIso");
    for (Long64_t ievt = startPoint; ievt < nEvents;ievt++) {
        //    for (Long64_t ievt=0; ievt<1000;ievt++)
        if (startPoint > 0 ) {
            if ((ievt / startPoint == 1) && (ievt % startPoint == 0)) cout << "ievt at start point: " << ievt << endl;
        }
        if (ievt%10000 == 0) cout << ievt << endl;
        if (ievt == endPoint) {
            cout << "ievt at end point (note, not running on this event!):" << ievt << endl;
            break;
        }
        diBJetPhi = 0; diBJetPt = 0; diBJetEta = 0; diBJetInvMass = 0;
        diJetPhi = 0; diJetPt = 0; diJetEta = 0; diJetInvMass = 0;
        MT2lb = 0;
        //        doEvent = true;
        map<string, float> stringKeyToVar;
        fileTree.GetEntry(ievt);
        
        if (whichNTupleType == 0) {
            if (doData) {
                if (Type == 0 && !passTrigDoubleMu) continue;
                else if (Type == 1 && !passTrigDoubleEl) continue;
                else if (Type == 2 && !passTrigElMu) continue;
            }
            else {
                if (!passTrigDoubleMu) {
                    if (Type == 0) doEvent = false;
                    if (Type_LepESUp == 0) doEvent_LepESUp = false;
                    if (Type_LepESDown == 0) doEvent_LepESDown = false;
                }
                if (!passTrigDoubleEl) {
                    if (Type == 1) doEvent = false;
                    if (Type_LepESUp == 1) doEvent_LepESUp = false;
                    if (Type_LepESDown == 1) doEvent_LepESDown = false;
                }
                if (!passTrigElMu) {
                    if (Type == 2) doEvent = false;
                    if (Type_LepESUp == 2) doEvent_LepESUp = false;
                    if (Type_LepESDown == 2) doEvent_LepESDown = false;
                }
            }
        }
        /*
         cout << "TGenStopMass0 " << TGenStopMass0 << endl;
         cout << "TGenStopMass1 " << TGenStopMass1 << endl;
         cout << "TGenChi0Mass0 " << TGenChi0Mass0 << endl;
         cout << "TGenChi0Mass1 " << TGenChi0Mass1 << endl;
         cout << "TGenCharginoMass0 " << TGenCharginoMass0 << endl;
         cout << "TGenCharginoMass1 " << TGenCharginoMass1 << endl;
         */        
        if (doData) doEvent = true;
        if (isSignal) {
            if (abs(TGenStopMass0 - grabStopMass) > massDiffThresh) continue;
            if (abs(TGenStopMass1 - grabStopMass) > massDiffThresh) continue;
            if (abs(TGenChi0Mass0 - grabChi0Mass) > massDiffThresh) continue;
            if (abs(TGenChi0Mass1 - grabChi0Mass) > massDiffThresh) continue;
            if (abs(TGenCharginoMass0 - grabCharginoMass) > massDiffThresh) continue;
            if (abs(TGenCharginoMass1 - grabCharginoMass) > massDiffThresh) continue;
        }
        //        cout << "Lep0 Px " << Lep0Px << endl;
        ZVeto = true;
        ZVeto_LepESUp = true;
        ZVeto_LepESDown = true;
        //*****************************************                                                                                                                   
        // SELECT THE SIGNAL MASS POINT (if signal)                                                                                                                   
        //*****************************************                                                                                                                   
        weight = 1.;
        weight = PUWeight;
        preNVtxRWweight = weight;
        ////Calculate all the event variables needed
        Lep0Vec.SetPxPyPzE(Lep0Px,Lep0Py,Lep0Pz,Lep0E);
        Lep1Vec.SetPxPyPzE(Lep1Px,Lep1Py,Lep1Pz,Lep1E);
        if (Lep1Vec.Pt() < subLepPtCut) {
            doEvent = false;
        }
        if (NJets > 0) {
            Jet0Vec.SetPxPyPzE(Jet0Px,Jet0Py,Jet0Pz,Jet0E);
            if (NJets > 1) {
                Jet1Vec.SetPxPyPzE(Jet1Px,Jet1Py,Jet1Pz,Jet1E);
                DiJetVec = Jet0Vec + Jet1Vec;
                diJetInvMass = DiJetVec.M();
                diJetPt = DiJetVec.Pt();
                diJetEta = DiJetVec.Eta();
                diJetPhi = DiJetVec.Phi();
            }
            if (NBtagJets > 0) {
                BtagJet0Vec.SetPxPyPzE(BtagJet0Px,BtagJet0Py,BtagJet0Pz,BtagJet0E);
                if (NBtagJets > 1) {
                    BtagJet1Vec.SetPxPyPzE(BtagJet1Px,BtagJet1Py,BtagJet1Pz,BtagJet1E);
                    DiBJetVec = BtagJet0Vec + BtagJet1Vec;
                    diBJetInvMass = DiBJetVec.M();
                    diBJetPt = DiBJetVec.Pt();
                    diBJetEta = DiBJetVec.Eta();
                    diBJetPhi = DiBJetVec.Phi();
                }
            }
        }
        if (!doData) {
            Lep0Vec_LepESUp.SetPxPyPzE(Lep0Px_LepESUp, Lep0Py_LepESUp, Lep0Pz_LepESUp, Lep0E_LepESUp);
            Lep0Vec_LepESDown.SetPxPyPzE(Lep0Px_LepESDown, Lep0Py_LepESDown, Lep0Pz_LepESDown, Lep0E_LepESDown);
            Lep1Vec_LepESUp.SetPxPyPzE(Lep1Px_LepESUp, Lep1Py_LepESUp, Lep1Pz_LepESUp, Lep1E_LepESUp);
            Lep1Vec_LepESDown.SetPxPyPzE(Lep1Px_LepESDown, Lep1Py_LepESDown, Lep1Pz_LepESDown, Lep1E_LepESDown);
            if (Lep1Vec_LepESUp.Pt() < subLepPtCut) doEvent_LepESUp = false;
            if (Lep1Vec_LepESDown.Pt() < subLepPtCut) doEvent_LepESDown = false;                
            //                cout << "doEvent? " << doEvent << endl;
            //                cout << "doEvent_LepESUp? " << doEvent_LepESUp << endl;
            //                cout << "doEvent_LepESDown? " << doEvent_LepESDown << endl;
            if (NJets_JetESUp > 0) {
                Jet0Vec_JetESUp.SetPxPyPzE(Jet0Px_JetESUp, Jet0Py_JetESUp, Jet0Pz_JetESUp, Jet0E_JetESUp);
                if (NJets_JetESUp > 1) {
                    Jet1Vec_JetESUp.SetPxPyPzE(Jet1Px_JetESUp, Jet1Py_JetESUp, Jet1Pz_JetESUp, Jet1E_JetESUp);
                    DiJetVec_JetESUp = Jet0Vec_JetESUp + Jet1Vec_JetESUp;
                    diJetInvMass_JetESUp = DiJetVec_JetESUp.M();
                    //                cout << "diJetInvMass " << diJetInvMass << endl;
                    diJetPt_JetESUp = DiJetVec_JetESUp.Pt();
                    diJetEta_JetESUp = DiJetVec_JetESUp.Eta();
                    diJetPhi_JetESUp = DiJetVec_JetESUp.Phi();
                }
                if (NBtagJets_JetESUp > 0) {
                    BtagJet0Vec_JetESUp.SetPxPyPzE(BtagJet0Px_JetESUp, BtagJet0Py_JetESUp, BtagJet0Pz_JetESUp, BtagJet0E_JetESUp);
                    if (NBtagJets > 1) {
                        BtagJet1Vec_JetESUp.SetPxPyPzE(BtagJet1Px_JetESUp, BtagJet1Py_JetESUp, BtagJet1Pz_JetESUp, BtagJet1E_JetESUp);
                        DiBJetVec_JetESUp = BtagJet0Vec_JetESUp + BtagJet1Vec_JetESUp;
                        diBJetInvMass_JetESUp = DiBJetVec_JetESUp.M();
                        diBJetPt_JetESUp = DiBJetVec_JetESUp.Pt();
                        diBJetEta_JetESUp = DiBJetVec_JetESUp.Eta();
                        diBJetPhi_JetESUp = DiBJetVec_JetESUp.Phi();
                    }
                }
            }
            if (NJets_JetESDown > 0) {
                Jet0Vec_JetESDown.SetPxPyPzE(Jet0Px_JetESDown, Jet0Py_JetESDown, Jet0Pz_JetESDown, Jet0E_JetESDown);
                if (NJets_JetESDown > 1) {
                    Jet1Vec_JetESDown.SetPxPyPzE(Jet1Px_JetESDown, Jet1Py_JetESDown, Jet1Pz_JetESDown, Jet1E_JetESDown);
                    DiJetVec_JetESDown = Jet0Vec_JetESDown + Jet1Vec_JetESDown;
                    diJetInvMass_JetESDown = DiJetVec_JetESDown.M();
                    //                cout << "diJetInvMass " << diJetInvMass << endl;
                    diJetPt_JetESDown = DiJetVec_JetESDown.Pt();
                    diJetEta_JetESDown = DiJetVec_JetESDown.Eta();
                    diJetPhi_JetESDown = DiJetVec_JetESDown.Phi();
                }
                if (NBtagJets_JetESDown > 0) {
                    BtagJet0Vec_JetESDown.SetPxPyPzE(BtagJet0Px_JetESDown, BtagJet0Py_JetESDown, BtagJet0Pz_JetESDown, BtagJet0E_JetESDown);
                    if (NBtagJets > 1) {
                        BtagJet1Vec_JetESDown.SetPxPyPzE(BtagJet1Px_JetESDown, BtagJet1Py_JetESDown, BtagJet1Pz_JetESDown, BtagJet1E_JetESDown);
                        DiBJetVec_JetESDown = BtagJet0Vec_JetESDown + BtagJet1Vec_JetESDown;
                        diBJetInvMass_JetESDown = DiBJetVec_JetESDown.M();
                        diBJetPt_JetESDown = DiBJetVec_JetESDown.Pt();
                        diBJetEta_JetESDown = DiBJetVec_JetESDown.Eta();
                        diBJetPhi_JetESDown = DiBJetVec_JetESDown.Phi();
                    }
                }
            }
            //        cout << "nVtx " << nVtx << endl;
            //        cout << "Lep0Pt " << Lep0Vec.Pt() << endl;
            nVtxTrue = nVtx;
            if (doPURWOviToDESY && !doData)  weight *= PileupRW(nVtxSFHistOviToDESY, nVtx);
            //        cout << "weight pre doHack " << weight << endl;
            if (doHackPURW && !doData) weight = PileupRW(nVtxSFHist, nVtx);            
            if (!doData) {
                if (doVerbosity) {
                    cout << "Type " << Type << endl;
                    cout << "weight pre scale " << weight << endl;
                }
                if (Type_LepESUp == -2) Type_LepESUp = 2;
                if (Type_LepESDown == -2) Type_LepESDown = 2;
                //                cout << "weight pre scale factor MC " << weight << endl;
                //                cout << "weight scalefactor for Type = " << Type << " is " << ScaleFactorMC(Type) << endl;
                //                cout << "ScaleFactorMC_0 for Type " << Type << " is " << ScaleFactorMC(Type, 0) << endl;
                //                cout << "ScaleFactorMC_+ for Type " << Type << " is " << ScaleFactorMC(Type, +1) << endl;
                //                cout << "ScaleFactorMC_- for Type " << Type << " is " << ScaleFactorMC(Type, -1) << endl;
                //                cout << " weight " << weight << endl;
//                weight_LepESUp = weight;
//                weight_LepESDown = weight;
//                preNVtxRWweight_LepESUp = preNVtxRWweight;
//                preNVtxRWweight_LepESDown = preNVtxRWweight;
                weight *= ScaleFactorMC(Type, 0);
                /*
                weight_LepESUp *= ScaleFactorMC(Type_LepESUp, 0);
                weight_LepESDown *= ScaleFactorMC(Type_LepESDown, 0);
                */
                preNVtxRWweight *= ScaleFactorMC(Type, 0);
                /*
                preNVtxRWweight_LepESUp *= ScaleFactorMC(Type_LepESUp, 0);
                preNVtxRWweight_LepESDown *= ScaleFactorMC(Type_LepESDown, 0);
                */
                //                cout << "weight post scale factor MC " << weight << endl;
                
            }
            
            if (hasTopInfo && (fInName.Contains("TT") || fInName.Contains("ttbar"))) {
                if (genTop0PdgId * genTop1PdgId > 0) {
                    cout << "something is funky with the genTop PDGIDs" << endl;
                    continue;
                }
                else {
                    if (genTop0PdgId > 0) {
                        genTopPt = genTop0Pt;
                        genAntiTopPt = genTop1Pt;
                    }
                    else {
                        genTopPt = genTop1Pt;
                        genAntiTopPt = genTop0Pt;   
                    }                
                }
                weight_GenTopReweight = GenLevelTopPtWeight(genTopPt, genAntiTopPt);
                h_genTopPt->Fill(genTopPt, 1.);
                h_genTopPtRW->Fill(genTopPt, weight_GenTopReweight);
                h_genAntiTopPt->Fill(genAntiTopPt, 1.);
                h_genAntiTopPtRW->Fill(genAntiTopPt, weight_GenTopReweight);
            }
            else {
                weight_GenTopReweight = 1.;
            }
            preNVtxRWweight_GenTopReweight = weight_GenTopReweight;
            h_ElecRelIso->Fill(Lep0RelPFIso, weight);
            h_ElecRelIso->Fill(Lep1RelPFIso, weight);
            //        cout << "weight post doHack " << weight << endl;
            
        }
        // basic condition, if running on data, only select type of events that are relevant to prevent double counting
        if (doData) {
            if (Type == 0) {
                if (!(fInName.Contains("DoubleMu") || fInName.Contains("mumu_run2012"))) continue;
            }
            else if (Type == 1) {
                if (!(fInName.Contains("DoubleEl") || fInName.Contains("ee_run2012"))) continue;
            }
            else if (Type == 2) {
                if (!(fInName.Contains("MuEG") || fInName.Contains("emu_run2012"))) continue;
            }
        }    
        //^^^^ Maybe check the above
        DiLepVec = Lep0Vec + Lep1Vec;
        diLepInvMass = DiLepVec.M();
        if (diLepInvMass > ZWindowLB && diLepInvMass < ZWindowUB) ZVeto = false; //Mass is in the ZVeto region
        diLepPt = DiLepVec.Pt();
        diLepEta = DiLepVec.Eta();
        diLepPhi = DiLepVec.Phi();
        METPhi_preCorr = METPhi;
        METX_preCorr = MET*TMath::Cos(METPhi_preCorr);
        METY_preCorr = MET*TMath::Sin(METPhi_preCorr);
        METX = METX_preCorr;
        METY = METY_preCorr;
        if (doPhiCorr) MetPhiCorrect(doData, METX, METY, nVtx);        
        if (doMETSmear) {
            METX *= rand.Gaus(1, METSF * METX);   
            METY *= rand.Gaus(1, METSF * METY);   
        }
        METPhi = TMath::ATan2(METY, METX);
        MET = TMath::Sqrt(METX * METX + METY * METY);
        MT2ll=getMT2(Lep0Vec, Lep1Vec, MET, METPhi);
        switch (NJets) {
            case 0:
                METdivMeff = MET / (MET + Lep0Vec.Pt() + Lep1Vec.Pt());
                break;                
            case 1:
                METdivMeff = MET / (MET + Lep0Vec.Pt() + Lep1Vec.Pt() + Jet0Vec.Pt());
                break;
            default:
                METdivMeff = MET / (MET + Lep0Vec.Pt() + Lep1Vec.Pt() + Jet0Vec.Pt() + Jet1Vec.Pt());
                break;
        }
        
        /******************************************************/
        //        MT2lb=getMT2(Lep0Vec, BtagJet0Vec, MET, METPhi);  
        //        cout << "NJets " << NJets << endl;
        if (NJets > 1) {
            vecLepMT2lb[0] = Lep0Vec;
            vecLepMT2lb[1] = Lep1Vec;
            if (NBtagJets > 1) {
                vecJetMT2lb[0] = BtagJet0Vec;
                vecJetMT2lb[1] = BtagJet1Vec;
                caseMT2lb = 0;
                DeltaPhiMT2lb_JetsUsed = dPhi(BtagJet0Vec.Phi(), BtagJet1Vec.Phi());
            }
            else if (NBtagJets == 1) {
                vecJetMT2lb[0] = BtagJet0Vec;
                if (BtagJet0Index == 0) {
                    vecJetMT2lb[1] = Jet1Vec;
                    caseMT2lb = 1;
                    DeltaPhiMT2lb_JetsUsed = dPhi(BtagJet0Vec.Phi(), Jet1Vec.Phi());
                }
                else {
                    vecJetMT2lb[1] = Jet0Vec;
                    caseMT2lb = 2;
                    DeltaPhiMT2lb_JetsUsed = dPhi(BtagJet0Vec.Phi(), Jet0Vec.Phi());
                }
            }
            else {
                vecJetMT2lb[0] = Jet0Vec;
                vecJetMT2lb[1] = Jet1Vec;
                caseMT2lb = 3;
                DeltaPhiMT2lb_JetsUsed = dPhi(Jet0Vec.Phi(), Jet1Vec.Phi());
            }
            MT2lb = MT2lbCalculator(vecLepMT2lb, vecJetMT2lb, MET, METPhi, vecBLepsMT2lb);
            DeltaPhiMT2lb_BLepsUsed = dPhi(vecBLepsMT2lb[0].Phi(), vecBLepsMT2lb[1].Phi());
        }
        else {
            MT2lb = -99.;
            DeltaPhiMT2lb_JetsUsed = -99.;
            DeltaPhiMT2lb_BLepsUsed = -99.;
        }
        
        
        /******************************************************/
        //Systematics//
        if (!doData) {
            DiLepVec_LepESUp = Lep0Vec_LepESUp + Lep1Vec_LepESUp;
            diLepInvMass_LepESUp = DiLepVec_LepESUp.M();
            if (diLepInvMass_LepESUp > ZWindowLB && diLepInvMass_LepESUp < ZWindowUB) ZVeto_LepESUp = false; //Mass is in the ZVeto region
            diLepPt_LepESUp = DiLepVec_LepESUp.Pt();
            diLepEta_LepESUp = DiLepVec_LepESUp.Eta();
            diLepPhi_LepESUp = DiLepVec_LepESUp.Phi();
            DiLepVec_LepESDown = Lep0Vec_LepESDown + Lep1Vec_LepESDown;
            diLepInvMass_LepESDown = DiLepVec_LepESDown.M();
            if (diLepInvMass_LepESDown > ZWindowLB && diLepInvMass_LepESDown < ZWindowUB) ZVeto_LepESDown = false; //Mass is in the ZVeto region
            diLepPt_LepESDown = DiLepVec_LepESDown.Pt();
            diLepEta_LepESDown = DiLepVec_LepESDown.Eta();
            diLepPhi_LepESDown = DiLepVec_LepESDown.Phi();
            
            switch (NJets) {
                case 0:
                    METdivMeff_LepESUp = MET_LepESUp / (MET_LepESUp + Lep0Vec_LepESUp.Pt() + Lep1Vec_LepESUp.Pt());
                    METdivMeff_LepESDown = MET_LepESDown / (MET_LepESDown + Lep0Vec_LepESDown.Pt() + Lep1Vec_LepESDown.Pt());
                    break;                
                case 1:
                    METdivMeff_LepESUp = MET_LepESUp / (MET_LepESUp + Lep0Vec_LepESUp.Pt() + Lep1Vec_LepESUp.Pt() + Jet0Vec.Pt());
                    METdivMeff_LepESDown = MET_LepESDown / (MET_LepESDown + Lep0Vec_LepESDown.Pt() + Lep1Vec_LepESDown.Pt() + Jet0Vec.Pt());
                    break;
                default:
                    METdivMeff_LepESUp = MET_LepESUp / (MET_LepESUp + Lep0Vec_LepESUp.Pt() + Lep1Vec_LepESUp.Pt() + Jet0Vec.Pt() + Jet1Vec.Pt());
                    METdivMeff_LepESDown = MET_LepESDown / (MET_LepESDown + Lep0Vec_LepESDown.Pt() + Lep1Vec_LepESDown.Pt() + Jet0Vec.Pt() + Jet1Vec.Pt());
                    break;
            }
            switch (NJets_JetESUp) {
                case 0:
                    METdivMeff_JetESUp = MET_JetESUp / (MET_JetESUp + Lep0Vec.Pt() + Lep1Vec.Pt());
                    break;                
                case 1:
                    METdivMeff_JetESUp = MET_JetESUp / (MET_JetESUp + Lep0Vec.Pt() + Lep1Vec.Pt() + Jet0Vec_JetESUp.Pt());
                    break;
                default:
                    METdivMeff_JetESUp = MET_JetESUp / (MET_JetESUp + Lep0Vec.Pt() + Lep1Vec.Pt() + Jet0Vec_JetESUp.Pt() + Jet1Vec_JetESUp.Pt());
                    break;
            }
            switch (NJets_JetESDown) {
                case 0:
                    METdivMeff_JetESDown = MET_JetESDown / (MET_JetESDown + Lep0Vec.Pt() + Lep1Vec.Pt());
                    break;                
                case 1:
                    METdivMeff_JetESDown = MET_JetESDown / (MET_JetESDown + Lep0Vec.Pt() + Lep1Vec.Pt() + Jet0Vec_JetESDown.Pt());
                    break;
                default:
                    METdivMeff_JetESDown = MET_JetESDown / (MET_JetESDown + Lep0Vec.Pt() + Lep1Vec.Pt() + Jet0Vec_JetESDown.Pt() + Jet1Vec_JetESDown.Pt());
                    break;
            }
            //correct the systematic shifted METs
            METPhi_preCorr_LepESUp = METPhi_LepESUp;
            METX_preCorr_LepESUp = MET_LepESUp*TMath::Cos(METPhi_preCorr_LepESUp);
            METY_preCorr_LepESUp = MET_LepESUp*TMath::Sin(METPhi_preCorr_LepESUp);
            METX_LepESUp = METX_preCorr_LepESUp;
            METY_LepESUp = METY_preCorr_LepESUp;
            if (doPhiCorr) MetPhiCorrect(doData, METX_LepESUp, METY_LepESUp, nVtx);        
            if (doMETSmear) {
                METX_LepESUp *= rand.Gaus(1, METSF * METX_LepESUp);   
                METY_LepESUp *= rand.Gaus(1, METSF * METY_LepESUp);   
            }
            METPhi_LepESUp = TMath::ATan2(METY_LepESUp, METX_LepESUp);
            MET_LepESUp = TMath::Sqrt(METX_LepESUp * METX_LepESUp + METY_LepESUp * METY_LepESUp);
            
            METPhi_preCorr_LepESDown = METPhi_LepESDown;
            METX_preCorr_LepESDown = MET_LepESDown*TMath::Cos(METPhi_preCorr_LepESDown);
            METY_preCorr_LepESDown = MET_LepESDown*TMath::Sin(METPhi_preCorr_LepESDown);
            METX_LepESDown = METX_preCorr_LepESDown;
            METY_LepESDown = METY_preCorr_LepESDown;
            if (doPhiCorr) MetPhiCorrect(doData, METX_LepESDown, METY_LepESDown, nVtx);        
            if (doMETSmear) {
                METX_LepESDown *= rand.Gaus(1, METSF * METX_LepESDown);   
                METY_LepESDown *= rand.Gaus(1, METSF * METY_LepESDown);   
            }
            METPhi_LepESDown = TMath::ATan2(METY_LepESDown, METX_LepESDown);
            MET_LepESDown = TMath::Sqrt(METX_LepESDown * METX_LepESDown + METY_LepESDown * METY_LepESDown);
            
            METPhi_preCorr_JetESUp = METPhi_JetESUp;
            METX_preCorr_JetESUp = MET_JetESUp*TMath::Cos(METPhi_preCorr_JetESUp);
            METY_preCorr_JetESUp = MET_JetESUp*TMath::Sin(METPhi_preCorr_JetESUp);
            METX_JetESUp = METX_preCorr_JetESUp;
            METY_JetESUp = METY_preCorr_JetESUp;
            if (doPhiCorr) MetPhiCorrect(doData, METX_JetESUp, METY_JetESUp, nVtx);        
            if (doMETSmear) {
                METX_JetESUp *= rand.Gaus(1, METSF * METX_JetESUp);   
                METY_JetESUp *= rand.Gaus(1, METSF * METY_JetESUp);   
            }
            METPhi_JetESUp = TMath::ATan2(METY_JetESUp, METX_JetESUp);
            MET_JetESUp = TMath::Sqrt(METX_JetESUp * METX_JetESUp + METY_JetESUp * METY_JetESUp);
            
            METPhi_preCorr_JetESDown = METPhi_JetESDown;
            METX_preCorr_JetESDown = MET_JetESDown*TMath::Cos(METPhi_preCorr_JetESDown);
            METY_preCorr_JetESDown = MET_JetESDown*TMath::Sin(METPhi_preCorr_JetESDown);
            METX_JetESDown = METX_preCorr_JetESDown;
            METY_JetESDown = METY_preCorr_JetESDown;
            if (doPhiCorr) MetPhiCorrect(doData, METX_JetESDown, METY_JetESDown, nVtx);        
            if (doMETSmear) {
                METX_JetESDown *= rand.Gaus(1, METSF * METX_JetESDown);   
                METY_JetESDown *= rand.Gaus(1, METSF * METY_JetESDown);   
            }
            METPhi_JetESDown = TMath::ATan2(METY_JetESDown, METX_JetESDown);
            MET_JetESDown = TMath::Sqrt(METX_JetESDown * METX_JetESDown + METY_JetESDown * METY_JetESDown);
            MT2llSmearFactorBin = MT2llMeanSmear->FindBin(MT2ll);
            if (MT2llSmearFactorBin > 20) MT2llSmearFactorBin = 21;
            MT2llSmearFactor = MT2llMeanSmear->GetBinContent(MT2llSmearFactorBin);
            //        cout << "MT2ll " << MT2ll << endl;
            //        cout << "MT2ll Smear " << MT2llSmearFactor << endl;
            MT2ll_ShiftUp = MT2ll + rand.Gaus(0, MT2llSmearFactor);
            if (MT2ll_ShiftUp < 0) MT2ll_ShiftUp = 0;
            MT2ll_LepESUp = getMT2(Lep0Vec_LepESUp, Lep1Vec_LepESUp, MET_LepESUp, METPhi_LepESUp);
            MT2ll_LepESDown = getMT2(Lep0Vec_LepESDown, Lep1Vec_LepESDown, MET_LepESDown, METPhi_LepESDown);
            MT2ll_JetESUp = getMT2(Lep0Vec, Lep1Vec, MET_JetESUp, METPhi_JetESUp);
            MT2ll_JetESDown = getMT2(Lep0Vec, Lep1Vec, MET_JetESDown, METPhi_JetESDown);
            MT2ll_UncESUp = MT2ll - DeltaMT2UncEn(vecOneDeeMT2llUncEnUp, MT2llUncEnUpDelta2D, MT2ll);
            MT2ll_UncESDown = MT2ll - DeltaMT2UncEn(vecOneDeeMT2llUncEnDown, MT2llUncEnDownDelta2D, MT2ll);
            if (MT2ll_UncESUp < 0) MT2ll_UncESUp = 0;
            if (MT2ll_UncESDown < 0) MT2ll_UncESDown = 0;
            if (NJets > 1) {
                vecLepMT2lb_LepESUp[0] = Lep0Vec_LepESUp;
                vecLepMT2lb_LepESUp[1] = Lep1Vec_LepESUp;
                MT2lb_LepESUp = MT2lbCalculator(vecLepMT2lb_LepESUp, vecJetMT2lb, MET_LepESUp, METPhi_LepESUp, vecBLepsMT2lb_LepESUp);
                DeltaPhiMT2lb_BLepsUsed_LepESUp = dPhi(vecBLepsMT2lb_LepESUp[0].Phi(), vecBLepsMT2lb_LepESUp[1].Phi());
                vecLepMT2lb_LepESDown[0] = Lep0Vec_LepESDown;
                vecLepMT2lb_LepESDown[1] = Lep1Vec_LepESDown;
                MT2lb_LepESDown = MT2lbCalculator(vecLepMT2lb_LepESDown, vecJetMT2lb, MET_LepESDown, METPhi_LepESDown, vecBLepsMT2lb_LepESDown);
                DeltaPhiMT2lb_BLepsUsed_LepESDown = dPhi(vecBLepsMT2lb_LepESDown[0].Phi(), vecBLepsMT2lb_LepESDown[1].Phi());                
                MT2lb_UncESUp = MT2lb - DeltaMT2UncEn(vecOneDeeMT2lbUncEnUp, MT2lbUncEnUpDelta2D, MT2lb);                
                MT2lb_UncESDown = MT2lb - DeltaMT2UncEn(vecOneDeeMT2lbUncEnDown, MT2lbUncEnDownDelta2D, MT2lb);
                if (MT2lb_UncESUp < 0) MT2lb_UncESUp = 0;
                if (MT2lb_UncESDown < 0) MT2lb_UncESDown = 0;
            }
            else {
                MT2lb_LepESUp = -99.;
                MT2lb_LepESDown = -99.;
                MT2lb_UncESUp = -99.;
                MT2lb_UncESDown = -99.;
                DeltaPhiMT2lb_BLepsUsed_LepESUp = -99.;
                DeltaPhiMT2lb_BLepsUsed_LepESDown = -99.;
            }
            if (NJets_JetESUp > 1) {
                if (NBtagJets_JetESUp > 1) {
                    vecJetMT2lb_JetESUp[0] = BtagJet0Vec_JetESUp;
                    vecJetMT2lb_JetESUp[1] = BtagJet1Vec_JetESUp;
                    caseMT2lb_JetESUp = 0;
                    DeltaPhiMT2lb_JetsUsed_JetESUp = dPhi(BtagJet0Vec_JetESUp.Phi(), BtagJet1Vec_JetESUp.Phi());
                }
                else if (NBtagJets_JetESUp == 1) {                    
                    vecJetMT2lb_JetESUp[0] = BtagJet0Vec_JetESUp;
                    if (BtagJet0Index_JetESUp == 0) {
                        vecJetMT2lb_JetESUp[1] = Jet1Vec_JetESUp;
                        caseMT2lb_JetESUp = 1;
                        DeltaPhiMT2lb_JetsUsed_JetESUp = dPhi(BtagJet0Vec_JetESUp.Phi(), Jet1Vec_JetESUp.Phi());                        
                    }
                    else {
                        vecJetMT2lb_JetESUp[1] = Jet0Vec_JetESUp;
                        caseMT2lb_JetESUp = 2;
                        DeltaPhiMT2lb_JetsUsed_JetESUp = dPhi(BtagJet0Vec_JetESUp.Phi(), Jet0Vec_JetESUp.Phi()); 
                    }
                }
                else {
                    vecJetMT2lb_JetESUp[0] = Jet0Vec_JetESUp;
                    vecJetMT2lb_JetESUp[1] = Jet1Vec_JetESUp;
                    caseMT2lb_JetESUp = 3;
                    DeltaPhiMT2lb_JetsUsed_JetESUp = dPhi(Jet0Vec_JetESUp.Phi(), Jet1Vec_JetESUp.Phi());
                }
                MT2lb_JetESUp = MT2lbCalculator(vecLepMT2lb, vecJetMT2lb_JetESUp, MET_JetESUp, METPhi_JetESUp, vecBLepsMT2lb_JetESUp);
                DeltaPhiMT2lb_BLepsUsed_JetESUp = dPhi(vecBLepsMT2lb_JetESUp[0].Phi(), vecBLepsMT2lb_JetESUp[1].Phi());
            }
            else {
                MT2lb_JetESUp = -99.;
                DeltaPhiMT2lb_JetsUsed_JetESUp = -99.;
                DeltaPhiMT2lb_BLepsUsed_JetESUp = -99.;
            }
            if (NJets_JetESDown > 1) {
                if (NBtagJets_JetESDown > 1) {
                    vecJetMT2lb_JetESDown[0] = BtagJet0Vec_JetESDown;
                    vecJetMT2lb_JetESDown[1] = BtagJet1Vec_JetESDown;
                    caseMT2lb_JetESDown = 0;
                    DeltaPhiMT2lb_JetsUsed_JetESDown = dPhi(BtagJet0Vec_JetESDown.Phi(), BtagJet1Vec_JetESDown.Phi());
                }
                else if (NBtagJets_JetESDown == 1) {                    
                    vecJetMT2lb_JetESDown[0] = BtagJet0Vec_JetESDown;
                    if (BtagJet0Index_JetESDown == 0) {
                        vecJetMT2lb_JetESDown[1] = Jet1Vec_JetESDown;
                        caseMT2lb_JetESDown = 1;
                        DeltaPhiMT2lb_JetsUsed_JetESDown = dPhi(BtagJet0Vec_JetESDown.Phi(), Jet1Vec_JetESDown.Phi());                        
                    }
                    else {
                        vecJetMT2lb_JetESDown[1] = Jet0Vec_JetESDown;
                        caseMT2lb_JetESDown = 2;
                        DeltaPhiMT2lb_JetsUsed_JetESDown = dPhi(BtagJet0Vec_JetESDown.Phi(), Jet0Vec_JetESDown.Phi()); 
                    }
                }
                else {
                    vecJetMT2lb_JetESDown[0] = Jet0Vec_JetESDown;
                    vecJetMT2lb_JetESDown[1] = Jet1Vec_JetESDown;
                    caseMT2lb_JetESDown = 3;
                    DeltaPhiMT2lb_JetsUsed_JetESDown = dPhi(Jet0Vec_JetESDown.Phi(), Jet1Vec_JetESDown.Phi());
                }
                MT2lb_JetESDown = MT2lbCalculator(vecLepMT2lb, vecJetMT2lb_JetESDown, MET_JetESDown, METPhi_JetESDown, vecBLepsMT2lb_JetESDown);
                DeltaPhiMT2lb_BLepsUsed_JetESDown = dPhi(vecBLepsMT2lb_JetESDown[0].Phi(), vecBLepsMT2lb_JetESDown[1].Phi());
            }
            else {
                MT2lb_JetESDown = -99.;
                DeltaPhiMT2lb_JetsUsed_JetESDown = -99.;
                DeltaPhiMT2lb_BLepsUsed_JetESDown = -99.;
            }            
        }        
        //        cout << "Mt2 " << MT2ll << endl;
        //        cout << "weight " << weight << endl;
        ///Set up the mapping of string keys to the appropriate event variables
        /*
         list of keys needed taken from command:
         cat StopFunctionDefinitions.h | grep VarKey
         _______________________________
         H_diLepPt.xVarKey = "diLepPt";
         H_diLepInvMass.xVarKey = "diLepInvMass";
         H_diLepEta.xVarKey = "diLepEta";
         H_diLepPhi.xVarKey = "diLepPhi";
         H_MT2ll.xVarKey = "MT2ll";
         H_MT2lb.xVarKey = "MT2lb";
         H_MET.xVarKey = "MET";
         H_METPhi.xVarKey = "METPhi";
         H_METPhi_noCorr.xVarKey = "METPhi_noPhiCorr";
         H_NJets.xVarKey = "NJets";
         H_NJetswBTag.xVarKey = "NBJets";
         H_diJetPt.xVarKey = "diJetPt";
         H_diJetInvMass.xVarKey = "diJetInvMass";
         H_diJetEta.xVarKey = "diJetEta";
         H_diJetPhi.xVarKey = "diJetPhi";
         H_diBJetPt.xVarKey = "diBJetPt";
         H_diBJetInvMass.xVarKey = "diBJetInvMass";
         H_diBJetEta.xVarKey = "diBJetEta";
         H_diBJetPhi.xVarKey = "diBJetPhi";
         H_DeltaPhiLep0Lep1.xVarKey = "DPhiLep0Lep1";
         H_DeltaPhiLep0Jet0.xVarKey = "DPhiLep0Jet0";
         H_DeltaPhiLep0Jet1.xVarKey = "DPhiLep0Jet1";
         H_DeltaPhiLep0BJet0.xVarKey = "DPhiLep0BJet0";
         H_DeltaPhiLep0BJet1.xVarKey = "DPhiLep0BJet1";
         H_DeltaPhiJet0BJet0.xVarKey = "DPhiJet0BJet0";
         H_DeltaPhiJet1BJet1.xVarKey = "DPhiJet1BJet1";
         H_DeltaPhiJet1BJet0.xVarKey = "DPhiJet1BJet0";
         H_nVtx.xVarKey = "nVtx";
         H_METX_vs_nVtx.xVarKey = "nVtx";
         H_METX_vs_nVtx.yVarKey = "METX";
         H_METX_vs_nVtx.xVarKey = "nVtx";
         H_METX_vs_nVtx.yVarKey = "METY";
         H_METX_vs_nVtx.xVarKey = "nVtx";
         H_METX_vs_nVtx.yVarKey = "METX_noPhiCorr";
         H_METX_vs_nVtx.xVarKey = "nVtx";
         H_METX_vs_nVtx.yVarKey = "METY_noPhiCorr";
         -----------------------------
         */
        stringKeyToVar["leadLepPt"] = Lep0Vec.Pt();
        stringKeyToVar["leadLepEta"] = Lep0Vec.Eta();
        stringKeyToVar["subLepPt"] = Lep1Vec.Pt();
        stringKeyToVar["subLepEta"] = Lep1Vec.Eta();
        stringKeyToVar["lep0RelPFIso"] = Lep0RelPFIso;
        stringKeyToVar["lep1RelPFIso"] = Lep1RelPFIso;
        stringKeyToVar["CutFlowEntry"] = 4;
        stringKeyToVar["diLepPt"] = diLepPt;
        stringKeyToVar["diLepInvMass"] = diLepInvMass;
        stringKeyToVar["diLepEta"] = diLepEta;
        stringKeyToVar["diLepPhi"] = diLepPhi;
        stringKeyToVar["MT2ll"] = MT2ll;
        //        cout << "MT2ll " << MT2ll << endl;
        //        cout << "PassMT2llCut80 " << (MT2ll > 80.) << endl;
        stringKeyToVar["PassMT2llCut80"] = (MT2ll > 80.);
        stringKeyToVar["PassMT2llCut90"] = (MT2ll > 90.);
        stringKeyToVar["PassMT2llCut100"] = (MT2ll > 100.);
        stringKeyToVar["PassMT2llCut110"] = (MT2ll > 110.);
        stringKeyToVar["PassMT2llCut120"] = (MT2ll > 120.);
        stringKeyToVar["MT2lb"] = MT2lb;
        stringKeyToVar["MET"] = MET;
        stringKeyToVar["METPhi"] = METPhi;
        stringKeyToVar["METPhi_noPhiCorr"] = METPhi_preCorr;
        stringKeyToVar["NJets"] = NJets;
        stringKeyToVar["NBJets"] = NBtagJets;
        stringKeyToVar["DPhiLep0Lep1"] = dPhi(Lep0Vec.Phi(), Lep1Vec.Phi());
        stringKeyToVar["DPhiLep0MET"] = dPhi((float) Lep0Vec.Phi(), METPhi);
        stringKeyToVar["DPhiLep1MET"] = dPhi((float) Lep0Vec.Phi(), METPhi);
        stringKeyToVar["DPhiLep0MET_PreCorr"] = dPhi((float) Lep0Vec.Phi(), METPhi_preCorr);
        stringKeyToVar["DPhiLep1MET_PreCorr"] = dPhi((float) Lep0Vec.Phi(), METPhi_preCorr);
        stringKeyToVar["DPhiZMET"] = dPhi(diLepPhi, METPhi);
        stringKeyToVar["DPhiZMET_PreCorr"] = dPhi(diLepPhi, METPhi_preCorr);
        stringKeyToVar["nVtx"] = (float) nVtx;
        stringKeyToVar["nVtxTrue"] = (float) nVtxTrue;
        stringKeyToVar["METX"] = METX;
        stringKeyToVar["METY"] = METY;
        stringKeyToVar["METX_noPhiCorr"] = METX_preCorr;
        stringKeyToVar["METY_noPhiCorr"] = METY_preCorr;
        stringKeyToVar["METdivMeff"] = METdivMeff;
        stringKeyToVar["METdivMeff_PassMT2llCut80"] = MT2ll > 80 ? METdivMeff : -99;
        stringKeyToVar["METdivMeff_PassMT2llCut90"] = MT2ll > 90 ? METdivMeff : -99;
        stringKeyToVar["METdivMeff_PassMT2llCut190"] = MT2ll > 100 ? METdivMeff : -99;
        stringKeyToVar["METdivMeff_PassMT2llCut110"] = MT2ll > 110 ? METdivMeff : -99;
        stringKeyToVar["METdivMeff_PassMT2llCut120"] = MT2ll > 120 ? METdivMeff : -99;
        stringKeyToVar["HT"] = HT;                    
        ///Conditional setting of maps for jet info -- depends upon number of jets in event
        if (NJets > 0) {
            stringKeyToVar["leadJetPt"] = Jet0Vec.Pt();
            stringKeyToVar["leadJetEta"] = Jet0Vec.Eta();
            stringKeyToVar["DPhiLep0Jet0"] = dPhi(Lep0Vec.Phi(), Jet0Vec.Phi());
            if (NJets > 1) {
                stringKeyToVar["subJetPt"] = Jet1Vec.Pt();
                stringKeyToVar["subJetEta"] = Jet1Vec.Eta();                
                stringKeyToVar["diJetPt"] = diJetPt;
                stringKeyToVar["diJetInvMass"] = diJetInvMass;
                stringKeyToVar["diJetEta"] = diJetEta;
                stringKeyToVar["diJetPhi"] = diJetPhi;
                stringKeyToVar["DPhiLep0Jet1"] = dPhi(Lep0Vec.Phi(), Jet1Vec.Phi());
                stringKeyToVar["ELepEJet"] = Lep0Vec.E() + Lep1Vec.E() - Jet0Vec.E() - Jet1Vec.E();
                stringKeyToVar["DPhiLepB0LepB1"] = dPhi(vecBLepsMT2lb[0].Phi(), vecBLepsMT2lb[1].Phi());
            }
            if (NBtagJets > 0) {
                stringKeyToVar["leadBJetPt"] = BtagJet0Vec.Pt();
                stringKeyToVar["leadBJetEta"] = BtagJet0Vec.Eta();
                stringKeyToVar["leadBJetEn"] = BtagJet0Vec.E();
                //                cout << "leadBJet En " << BtagJet0Vec.E() << endl;
                stringKeyToVar["DPhiLep0BJet0"] = dPhi(Lep0Vec.Phi(), BtagJet0Vec.Phi());
                stringKeyToVar["DPhiJet0BJet0"] = dPhi(Jet0Vec.Phi(), BtagJet0Vec.Phi());                
                if (NBtagJets > 1) {
                    stringKeyToVar["subBJetPt"] = BtagJet1Vec.Pt();
                    stringKeyToVar["subBJetEta"] = BtagJet1Vec.Eta(); 
                    stringKeyToVar["subBJetEn"] = BtagJet1Vec.E();
                    stringKeyToVar["diBJetPt"] = diBJetPt;
                    stringKeyToVar["diBJetInvMass"] = diBJetInvMass;
                    stringKeyToVar["diBJetEta"] = diBJetEta;
                    stringKeyToVar["diBJetPhi"] = diBJetPhi;                    
                    stringKeyToVar["DPhiLep0BJet1"] = dPhi(Lep0Vec.Phi(), BtagJet1Vec.Phi());
                    stringKeyToVar["DPhiJet1BJet0"] = dPhi(Jet1Vec.Phi(), BtagJet0Vec.Phi());
                    stringKeyToVar["DPhiJet1BJet1"] = dPhi(Jet1Vec.Phi(), BtagJet1Vec.Phi());
                }
            }   
        }
        if (doBookSyst) {
            stringKeyToVar["leadLepPt_LepESShiftUp"] = Lep0Vec_LepESUp.Pt();
            stringKeyToVar["leadLepPt_LepESShiftDown"] = Lep0Vec_LepESDown.Pt();
            stringKeyToVar["subLepPt_LepESShiftUp"] = Lep1Vec_LepESUp.Pt();
            stringKeyToVar["subLepPt_LepESShiftDown"] = Lep1Vec_LepESDown.Pt();
            stringKeyToVar["diLepPt_LepESShiftUp"] = diLepPt_LepESUp;
            stringKeyToVar["diLepPt_LepESShiftDown"] = diLepPt_LepESDown;
            stringKeyToVar["diLepInvMass_LepESShiftUp"] = diLepInvMass_LepESUp;         
            stringKeyToVar["diLepInvMass_LepESShiftDown"] = diLepInvMass_LepESDown;   
            stringKeyToVar["diLepEta_LepESShiftUp"] = diLepEta_LepESUp;
            stringKeyToVar["diLepEta_LepESShiftDown"] = diLepEta_LepESDown;
            stringKeyToVar["diLepPhi_LepESShiftUp"] = diLepPhi_LepESUp;
            stringKeyToVar["diLepPhi_LepESShiftDown"] = diLepPhi_LepESDown;
            stringKeyToVar["HT_JetESShiftUp"] = HT_JetESUp;
            stringKeyToVar["HT_JetESShiftDown"] = HT_JetESDown;
            /*
            stringKeyToVar["leadLepPt_JetESShiftUp"] = Lep0Vec.Pt();
            stringKeyToVar["leadLepPt_JetESShiftDown"] = Lep0Vec.Pt();
            stringKeyToVar["subLepPt_JetESShiftUp"] = Lep1Vec.Pt();
            stringKeyToVar["subLepPt_JetESShiftDown"] = Lep1Vec.Pt();
            stringKeyToVar["diLepPt_JetESShiftUp"] = diLepPt;
            stringKeyToVar["diLepPt_JetESShiftDown"] = diLepPt;
            stringKeyToVar["diLepInvMass_JetESShiftUp"] = diLepInvMass;         
            stringKeyToVar["diLepInvMass_JetESShiftDown"] = diLepInvMass;   
            stringKeyToVar["diLepEta_JetESShiftUp"] = diLepEta;
            stringKeyToVar["diLepEta_JetESShiftDown"] = diLepEta;
            stringKeyToVar["diLepPhi_JetESShiftUp"] = diLepPhi;
            stringKeyToVar["diLepPhi_JetESShiftDown"] = diLepPhi;
            
            stringKeyToVar["NJets_LepESShiftUp"] = NJets;
            stringKeyToVar["NJets_LepESShiftDown"] = NJets;
            stringKeyToVar["NBJets_LepESShiftUp"] = NBtagJets;
            stringKeyToVar["NBJets_LepESShiftDown"] = NBtagJets;
             */
            stringKeyToVar["NJets_JetESShiftUp"] = NJets_JetESUp;
            stringKeyToVar["NJets_JetESShiftDown"] = NJets_JetESDown;
            stringKeyToVar["NBJets_JetESShiftUp"] = NBtagJets_JetESUp;
            stringKeyToVar["NBJets_JetESShiftDown"] = NBtagJets_JetESDown;
            
            if (NJets > 0) {
                stringKeyToVar["DPhiLep0Jet0_LepESShiftUp"] = dPhi(Lep0Vec_LepESUp.Phi(), Jet0Vec.Phi());
                stringKeyToVar["DPhiLep0Jet0_LepESShiftDown"] = dPhi(Lep0Vec_LepESDown.Phi(), Jet0Vec.Phi());
                stringKeyToVar["DPhiLep0Jet1_LepESShiftUp"] = dPhi(Lep0Vec_LepESUp.Phi(), Jet1Vec.Phi());
                stringKeyToVar["DPhiLep0Jet1_LepESShiftDown"] = dPhi(Lep0Vec_LepESDown.Phi(), Jet1Vec.Phi());
                stringKeyToVar["DPhiLepB0LepB1_LepESShiftUp"] = DeltaPhiMT2lb_BLepsUsed_LepESUp;
                stringKeyToVar["DPhiLepB0LepB1_LepESShiftDown"] = DeltaPhiMT2lb_BLepsUsed_LepESDown;
                stringKeyToVar["ELepEJet_LepESShiftUp"] = Lep0Vec_LepESUp.E() + Lep1Vec_LepESUp.E() - Jet0Vec.E() - Jet1Vec.E();
                stringKeyToVar["ELepEJet_LepESShiftUp"] = Lep0Vec_LepESDown.E() + Lep1Vec_LepESDown.E() - Jet0Vec.E() - Jet1Vec.E();
                if (NBtagJets > 0) {
                    stringKeyToVar["DPhiLep0BJet0_LepESShiftUp"] = dPhi(Lep0Vec_LepESUp.Phi(), BtagJet0Vec.Phi());
                    stringKeyToVar["DPhiLep0BJet0_LepESShiftDown"] = dPhi(Lep0Vec_LepESDown.Phi(), BtagJet0Vec.Phi());
                    if (NBtagJets > 1) {
                        stringKeyToVar["DPhiLep0BJet1_LepESShiftUp"] = dPhi(Lep0Vec_LepESUp.Phi(), BtagJet1Vec.Phi());
                        stringKeyToVar["DPhiLep0BJet1_LepESShiftDown"] = dPhi(Lep0Vec_LepESDown.Phi(), BtagJet1Vec.Phi());
                    }
                }
            }
            if (NJets_JetESUp > 0) {
                stringKeyToVar["leadJetPt_JetESShiftUp"] = Jet0Vec_JetESUp.Pt();
                stringKeyToVar["leadJetEta_JetESShiftUp"] = Jet0Vec_JetESUp.Eta();
                stringKeyToVar["DPhiLep0Jet0_JetESShiftUp"] = dPhi(Lep0Vec.Phi(), Jet0Vec_JetESUp.Phi());
                if (NJets_JetESUp > 1) {
                    stringKeyToVar["subJetPt_JetESShiftUp"] = Jet1Vec_JetESUp.Pt();
                    stringKeyToVar["subJetEta_JetESShiftUp"] = Jet1Vec_JetESUp.Eta();                
                    stringKeyToVar["diJetPt_JetESShiftUp"] = diJetPt_JetESUp;
                    stringKeyToVar["diJetInvMass_JetESShiftUp"] = diJetInvMass_JetESUp;
                    stringKeyToVar["diJetEta_JetESShiftUp"] = diJetEta_JetESUp;
                    stringKeyToVar["diJetPhi_JetESShiftUp"] = diJetPhi_JetESUp;
                    stringKeyToVar["DPhiLep0Jet1_JetESShiftUp"] = dPhi(Lep0Vec.Phi(), Jet1Vec_JetESUp.Phi());
                    stringKeyToVar["ELepEJet_JetESShiftUp"] = Lep0Vec.E() + Lep1Vec.E() - Jet0Vec_JetESUp.E() - Jet1Vec_JetESUp.E();
                    stringKeyToVar["DPhiLepB0LepB1_JetESShiftUp"] = DeltaPhiMT2lb_BLepsUsed_JetESUp;
                }
                if (NBtagJets_JetESUp > 0) {
                    stringKeyToVar["leadBJetPt_JetESShiftUp"] = BtagJet0Vec_JetESUp.Pt();
                    stringKeyToVar["leadBJetEta_JetESShiftUp"] = BtagJet0Vec_JetESUp.Eta();
                    stringKeyToVar["leadBJetEn_JetESShiftUp"] = BtagJet0Vec_JetESUp.E();
                    //                cout << "leadBJet En " << BtagJet0Vec.E() << endl;
                    stringKeyToVar["DPhiLep0BJet0_JetESShiftUp"] = dPhi(Lep0Vec.Phi(), BtagJet0Vec_JetESUp.Phi());
                    stringKeyToVar["DPhiJet0BJet0_JetESShiftUp"] = dPhi(Jet0Vec.Phi(), BtagJet0Vec_JetESUp.Phi());                
                    if (NBtagJets_JetESUp > 1) {
                        stringKeyToVar["subBJetPt_JetESShiftUp"] = BtagJet1Vec_JetESUp.Pt();
                        stringKeyToVar["subBJetEta_JetESShiftUp"] = BtagJet1Vec_JetESUp.Eta(); 
                        stringKeyToVar["subBJetEn_JetESShiftUp"] = BtagJet1Vec_JetESUp.E();
                        stringKeyToVar["diBJetPt_JetESShiftUp"] = diBJetPt_JetESUp;
                        stringKeyToVar["diBJetInvMass_JetESShiftUp"] = diBJetInvMass_JetESUp;
                        stringKeyToVar["diBJetEta_JetESShiftUp"] = diBJetEta_JetESUp;
                        stringKeyToVar["diBJetPhi_JetESShiftUp"] = diBJetPhi_JetESUp;
                        stringKeyToVar["DPhiLep0BJet1_JetESShiftUp"] = dPhi(Lep0Vec.Phi(), BtagJet1Vec_JetESUp.Phi());
                        stringKeyToVar["DPhiJet1BJet0_JetESShiftUp"] = dPhi(Jet1Vec_JetESUp.Phi(), BtagJet0Vec_JetESUp.Phi());
                        stringKeyToVar["DPhiJet1BJet1_JetESShiftUp"] = dPhi(Jet1Vec_JetESUp.Phi(), BtagJet1Vec_JetESUp.Phi());
                    }
                }   
            }
            if (NJets_JetESDown > 0) {
                stringKeyToVar["leadJetPt_JetESShiftDown"] = Jet0Vec_JetESDown.Pt();
                stringKeyToVar["leadJetEta_JetESShiftDown"] = Jet0Vec_JetESDown.Eta();
                stringKeyToVar["DPhiLep0Jet0_JetESShiftDown"] = dPhi(Lep0Vec.Phi(), Jet0Vec_JetESDown.Phi());
                if (NJets_JetESDown > 1) {
                    stringKeyToVar["subJetPt_JetESShiftDown"] = Jet1Vec_JetESDown.Pt();
                    stringKeyToVar["subJetEta_JetESShiftDown"] = Jet1Vec_JetESDown.Eta();                
                    stringKeyToVar["diJetPt_JetESShiftDown"] = diJetPt_JetESDown;
                    stringKeyToVar["diJetInvMass_JetESShiftDown"] = diJetInvMass_JetESDown;
                    stringKeyToVar["diJetEta_JetESShiftDown"] = diJetEta_JetESDown;
                    stringKeyToVar["diJetPhi_JetESShiftDown"] = diJetPhi_JetESDown;
                    stringKeyToVar["DPhiLep0Jet1_JetESShiftDown"] = dPhi(Lep0Vec.Phi(), Jet1Vec_JetESDown.Phi());
                    stringKeyToVar["ELepEJet_JetESShiftDown"] = Lep0Vec.E() + Lep1Vec.E() - Jet0Vec_JetESDown.E() - Jet1Vec_JetESDown.E();
                    stringKeyToVar["DPhiLepB0LepB1_JetESShiftDown"] = DeltaPhiMT2lb_BLepsUsed_JetESDown;
                }
                if (NBtagJets_JetESDown > 0) {
                    stringKeyToVar["leadBJetPt_JetESShiftDown"] = BtagJet0Vec_JetESDown.Pt();
                    stringKeyToVar["leadBJetEta_JetESShiftDown"] = BtagJet0Vec_JetESDown.Eta();
                    stringKeyToVar["leadBJetEn_JetESShiftDown"] = BtagJet0Vec_JetESDown.E();
                    //                cout << "leadBJet En " << BtagJet0Vec.E() << endl;
                    stringKeyToVar["DPhiLep0BJet0_JetESShiftDown"] = dPhi(Lep0Vec.Phi(), BtagJet0Vec_JetESDown.Phi());
                    stringKeyToVar["DPhiJet0BJet0_JetESShiftDown"] = dPhi(Jet0Vec.Phi(), BtagJet0Vec_JetESDown.Phi());                
                    if (NBtagJets_JetESDown > 1) {
                        stringKeyToVar["subBJetPt_JetESShiftDown"] = BtagJet1Vec_JetESDown.Pt();
                        stringKeyToVar["subBJetEta_JetESShiftDown"] = BtagJet1Vec_JetESDown.Eta(); 
                        stringKeyToVar["subBJetEn_JetESShiftDown"] = BtagJet1Vec_JetESDown.E();
                        stringKeyToVar["diBJetPt_JetESShiftDown"] = diBJetPt_JetESDown;
                        stringKeyToVar["diBJetInvMass_JetESShiftDown"] = diBJetInvMass_JetESDown;
                        stringKeyToVar["diBJetEta_JetESShiftDown"] = diBJetEta_JetESDown;
                        stringKeyToVar["diBJetPhi_JetESShiftDown"] = diBJetPhi_JetESDown;
                        stringKeyToVar["DPhiLep0BJet1_JetESShiftDown"] = dPhi(Lep0Vec.Phi(), BtagJet1Vec_JetESDown.Phi());
                        stringKeyToVar["DPhiJet1BJet0_JetESShiftDown"] = dPhi(Jet1Vec_JetESDown.Phi(), BtagJet0Vec_JetESDown.Phi());
                        stringKeyToVar["DPhiJet1BJet1_JetESShiftDown"] = dPhi(Jet1Vec_JetESDown.Phi(), BtagJet1Vec_JetESDown.Phi());
                    }
                }   
            }
            
            
            if (NJets > 1) {
                stringKeyToVar["ELepEJet_LepESShiftUp"] = Lep0Vec_LepESUp.E() + Lep1Vec_LepESUp.E() - Jet0Vec.E() - Jet1Vec.E();
                stringKeyToVar["ELepEJet_LepESShiftDown"] = Lep0Vec_LepESDown.E() + Lep1Vec_LepESDown.E() - Jet0Vec.E() - Jet1Vec.E();
                stringKeyToVar["DPhiLepB0LepB1_LepESShiftUp"] = DeltaPhiMT2lb_BLepsUsed_LepESUp;
                stringKeyToVar["DPhiLepB0LepB1_LepESShiftDown"] = DeltaPhiMT2lb_BLepsUsed_LepESDown;
            }
            if (NJets_JetESUp > 1) {
                stringKeyToVar["ELepEJet_JetESShiftUp"] = Lep0Vec.E() + Lep1Vec.E() - Jet0Vec_JetESUp.E() - Jet1Vec_JetESUp.E();
            }
            if (NJets_JetESDown > 1) {
                stringKeyToVar["ELepEJet_JetESShiftUp"] = Lep0Vec.E() + Lep1Vec.E() - Jet0Vec_JetESDown.E() - Jet1Vec_JetESDown.E();
            }
            stringKeyToVar["METdivMeff_LepESShiftUp"] = METdivMeff_LepESUp;
            stringKeyToVar["METdivMeff_PassMT2llCut80_LepESShiftUp"] = MT2ll_LepESUp > 80 ? METdivMeff_LepESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut90_LepESShiftUp"] = MT2ll_LepESUp > 90 ? METdivMeff_LepESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut190_LepESShiftUp"] = MT2ll_LepESUp > 100 ? METdivMeff_LepESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut110_LepESShiftUp"] = MT2ll_LepESUp > 110 ? METdivMeff_LepESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut120_LepESShiftUp"] = MT2ll_LepESUp > 120 ? METdivMeff_LepESUp : -99;                        
            stringKeyToVar["METdivMeff_LepESShiftDown"] = METdivMeff_LepESDown;
            stringKeyToVar["METdivMeff_PassMT2llCut80_LepESShiftDown"] = MT2ll_LepESDown > 80 ? METdivMeff_LepESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut90_LepESShiftDown"] = MT2ll_LepESDown > 90 ? METdivMeff_LepESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut190_LepESShiftDown"] = MT2ll_LepESDown > 100 ? METdivMeff_LepESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut110_LepESShiftDown"] = MT2ll_LepESDown > 110 ? METdivMeff_LepESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut120_LepESShiftDown"] = MT2ll_LepESDown > 120 ? METdivMeff_LepESDown : -99;
            
            stringKeyToVar["METdivMeff_JetESShiftUp"] = METdivMeff_JetESUp;
            stringKeyToVar["METdivMeff_PassMT2llCut80_JetESShiftUp"] = MT2ll_JetESUp > 80 ? METdivMeff_JetESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut90_JetESShiftUp"] = MT2ll_JetESUp > 90 ? METdivMeff_JetESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut190_JetESShiftUp"] = MT2ll_JetESUp > 100 ? METdivMeff_JetESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut110_JetESShiftUp"] = MT2ll_JetESUp > 110 ? METdivMeff_JetESUp : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut120_JetESShiftUp"] = MT2ll_JetESUp > 120 ? METdivMeff_JetESUp : -99;                        
            stringKeyToVar["METdivMeff_JetESShiftDown"] = METdivMeff_JetESDown;
            stringKeyToVar["METdivMeff_PassMT2llCut80_JetESShiftDown"] = MT2ll_JetESDown > 80 ? METdivMeff_JetESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut90_JetESShiftDown"] = MT2ll_JetESDown > 90 ? METdivMeff_JetESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut190_JetESShiftDown"] = MT2ll_JetESDown > 100 ? METdivMeff_JetESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut110_JetESShiftDown"] = MT2ll_JetESDown > 110 ? METdivMeff_JetESDown : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut120_JetESShiftDown"] = MT2ll_JetESDown > 120 ? METdivMeff_JetESDown : -99;
            
            stringKeyToVar["METdivMeff_PassMT2llCut80_MT2UncESShiftUp"] = MT2ll_UncESUp > 80 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut90_MT2UncESShiftUp"] = MT2ll_UncESUp > 90 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut190_MT2UncESShiftUp"] = MT2ll_UncESUp > 100 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut110_MT2UncESShiftUp"] = MT2ll_UncESUp > 110 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut120_MT2UncESShiftUp"] = MT2ll_UncESUp > 120 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut80_MT2UncESShiftDown"] = MT2ll_UncESDown > 80 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut90_MT2UncESShiftDown"] = MT2ll_UncESDown > 90 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut190_MT2UncESShiftDown"] = MT2ll_UncESDown > 100 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut110_MT2UncESShiftDown"] = MT2ll_UncESDown > 110 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut120_MT2UncESShiftDown"] = MT2ll_UncESDown > 120 ? METdivMeff : -99;
            
            
            stringKeyToVar["METdivMeff_PassMT2llCut80_MT2llShiftUp"] = MT2ll_ShiftUp > 80 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut90_MT2llShiftUp"] = MT2ll_ShiftUp > 90 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut100_MT2llShiftUp"] = MT2ll_ShiftUp > 100 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut110_MT2llShiftUp"] = MT2ll_ShiftUp > 110 ? METdivMeff : -99;
            stringKeyToVar["METdivMeff_PassMT2llCut120_MT2llShiftUp"] = MT2ll_ShiftUp > 120 ? METdivMeff : -99;
            
            stringKeyToVar["MT2ll_MT2llShiftUp"] = MT2ll_ShiftUp;
            stringKeyToVar["PassMT2llCut80_MT2llShiftUp"] = (MT2ll_ShiftUp > 80.);
            stringKeyToVar["PassMT2llCut90_MT2llShiftUp"] = (MT2ll_ShiftUp > 90.);
            stringKeyToVar["PassMT2llCut100_MT2llShiftUp"] = (MT2ll_ShiftUp > 100.);
            stringKeyToVar["PassMT2llCut110_MT2llShiftUp"] = (MT2ll_ShiftUp > 110.);
            stringKeyToVar["PassMT2llCut120_MT2llShiftUp"] = (MT2ll_ShiftUp > 120.);
            
            stringKeyToVar["MT2ll_LepESShiftUp"] = MT2ll_LepESUp;
            stringKeyToVar["PassMT2llCut80_LepESShiftUp"] = (MT2ll_LepESUp > 80.);
            stringKeyToVar["PassMT2llCut90_LepESShiftUp"] = (MT2ll_LepESUp > 90.);
            stringKeyToVar["PassMT2llCut100_LepESShiftUp"] = (MT2ll_LepESUp > 100.);
            stringKeyToVar["PassMT2llCut110_LepESShiftUp"] = (MT2ll_LepESUp > 110.);
            stringKeyToVar["PassMT2llCut120_LepESShiftUp"] = (MT2ll_LepESUp > 120.);
            
            stringKeyToVar["MT2ll_LepESShiftDown"] = MT2ll_LepESDown;
            stringKeyToVar["PassMT2llCut80_LepESShiftDown"] = (MT2ll_LepESDown > 80.);
            stringKeyToVar["PassMT2llCut90_LepESShiftDown"] = (MT2ll_LepESDown > 90.);
            stringKeyToVar["PassMT2llCut100_LepESShiftDown"] = (MT2ll_LepESDown > 100.);
            stringKeyToVar["PassMT2llCut110_LepESShiftDown"] = (MT2ll_LepESDown > 110.);
            stringKeyToVar["PassMT2llCut120_LepESShiftDown"] = (MT2ll_LepESDown > 120.);
            
            stringKeyToVar["MT2ll_JetESShiftUp"] = MT2ll_JetESUp;
            stringKeyToVar["PassMT2llCut80_JetESShiftUp"] = (MT2ll_JetESUp > 80.);
            stringKeyToVar["PassMT2llCut90_JetESShiftUp"] = (MT2ll_JetESUp > 90.);
            stringKeyToVar["PassMT2llCut100_JetESShiftUp"] = (MT2ll_JetESUp > 100.);
            stringKeyToVar["PassMT2llCut110_JetESShiftUp"] = (MT2ll_JetESUp > 110.);
            stringKeyToVar["PassMT2llCut120_JetESShiftUp"] = (MT2ll_JetESUp > 120.);
            
            stringKeyToVar["MT2ll_JetESShiftDown"] = MT2ll_JetESDown;
            stringKeyToVar["PassMT2llCut80_JetESShiftDown"] = (MT2ll_JetESDown > 80.);
            stringKeyToVar["PassMT2llCut90_JetESShiftDown"] = (MT2ll_JetESDown > 90.);
            stringKeyToVar["PassMT2llCut100_JetESShiftDown"] = (MT2ll_JetESDown > 100.);
            stringKeyToVar["PassMT2llCut110_JetESShiftDown"] = (MT2ll_JetESDown > 110.);
            stringKeyToVar["PassMT2llCut120_JetESShiftDown"] = (MT2ll_JetESDown > 120.);
            
            stringKeyToVar["MT2ll_MT2UncESShiftUp"] = MT2ll_UncESUp;
            stringKeyToVar["PassMT2llCut80_MT2UncESShiftUp"] = (MT2ll_UncESUp > 80.);
            stringKeyToVar["PassMT2llCut90_MT2UncESShiftUp"] = (MT2ll_UncESUp > 90.);
            stringKeyToVar["PassMT2llCut100_MT2UncESShiftUp"] = (MT2ll_UncESUp > 100.);
            stringKeyToVar["PassMT2llCut110_MT2UncESShiftUp"] = (MT2ll_UncESUp > 110.);
            stringKeyToVar["PassMT2llCut120_MT2UncESShiftUp"] = (MT2ll_UncESUp > 120.);
            
            stringKeyToVar["MT2ll_MT2UncESShiftDown"] = MT2ll_UncESDown;
            stringKeyToVar["PassMT2llCut80_MT2UncESShiftDown"] = (MT2ll_UncESDown > 80.);
            stringKeyToVar["PassMT2llCut90_MT2UncESShiftDown"] = (MT2ll_UncESDown > 90.);
            stringKeyToVar["PassMT2llCut100_MT2UncESShiftDown"] = (MT2ll_UncESDown > 100.);
            stringKeyToVar["PassMT2llCut110_MT2UncESShiftDown"] = (MT2ll_UncESDown > 110.);
            stringKeyToVar["PassMT2llCut120_MT2UncESShiftDown"] = (MT2ll_UncESDown > 120.);
            
            stringKeyToVar["MT2lb_LepESShiftUp"] = MT2lb_LepESUp;
            stringKeyToVar["MT2lb_LepESShiftDown"] = MT2lb_LepESDown;
            stringKeyToVar["MT2lb_JetESShiftUp"] = MT2lb_JetESUp;
            stringKeyToVar["MT2lb_JetESShiftDown"] = MT2lb_JetESDown;
            stringKeyToVar["MT2lb_MT2UncESShiftUp"] = MT2lb_UncESUp;
            stringKeyToVar["MT2lb_MT2UncESShiftDown"] = MT2lb_UncESDown;
            
            stringKeyToVar["MET_LepESShiftUp"] = MET_LepESUp;
            stringKeyToVar["MET_LepESShiftDown"] = MET_LepESDown;
            stringKeyToVar["METPhi_LepESShiftUp"] = METPhi_LepESUp;
            stringKeyToVar["METPhi_LepESShiftDown"] = METPhi_LepESDown;
            
            stringKeyToVar["MET_JetESShiftUp"] = MET_JetESUp;
            stringKeyToVar["MET_JetESShiftDown"] = MET_JetESDown;
            stringKeyToVar["METPhi_JetESShiftUp"] = METPhi_JetESUp;
            stringKeyToVar["METPhi_JetESShiftDown"] = METPhi_JetESDown;
        }
        /*#######################
         MAKE SELECTION CUTS
         #####################*/
        for (unsigned int iCentVal = 0; iCentVal < subSampVec->size(); ++iCentVal) {
            S_Current = subSampVec->at(iCentVal);
            subSampBool[S_Current] = false;
            if (!doEvent) continue;
            if (!(S_Current.whichdiLepType < 0 || Type == S_Current.whichdiLepType)) continue;
            if (!(NJets >= S_Current.cutNJets)) continue;
            if (!(NBtagJets >= S_Current.cutNBJets)) continue;
            if (!(Type == 2 && S_Current.histNameSuffix.Contains("FullCut"))) {
                if (MET < S_Current.cutMET) continue;
                if (!(S_Current.doZVeto < 0 || ZVeto == S_Current.doZVeto)) continue;
            }
            if (S_Current.histNameSuffix.Contains("BothinBarrel")) {
                if (!(TMath::Abs(Lep0Vec.Eta()) < barrelEtaEnd && TMath::Abs(Lep1Vec.Eta()) < barrelEtaEnd)) continue;
            }
            if (S_Current.histNameSuffix.Contains("OneinBarrel")) {            
                if (!(TMath::Abs(Lep0Vec.Eta()) > barrelEtaEnd || TMath::Abs(Lep1Vec.Eta()) > barrelEtaEnd)) continue;
                if (!(TMath::Abs(Lep0Vec.Eta()) < barrelEtaEnd || TMath::Abs(Lep1Vec.Eta()) < barrelEtaEnd)) continue;                                        
            }   
            if (S_Current.histNameSuffix.Contains("BothinEndcap")) {
                if (!(TMath::Abs(Lep0Vec.Eta()) > endcapEtaStart && TMath::Abs(Lep1Vec.Eta()) > endcapEtaStart)) continue;
            }
            if (S_Current.histNameSuffix.Contains("0BJets")) {
                if (NBtagJets > 0) {
                    if (!(S_Current.histNameSuffix.Contains("inZMass") && !ZVeto)) continue;
                }
            }
            if (S_Current.histNameSuffix.Contains("_0Jets") && NJets != 0) continue;
            if (S_Current.histNameSuffix.Contains("_1Jet") && NJets != 1) continue;
            if (S_Current.histNameSuffix.Contains("FullCutBlind")) {
                if (MT2ll > 80) continue;
            }
            subSampBool[S_Current] = true;
            if (doVerbosity) {
                if (S_Current.histNameSuffix == "_NJetsGeq2") {
                    cout << "" << endl;
                    if (subSampBool[S_Current] == true) cout << "true!" << endl;
                    if (subSampBool[S_Current] == false) cout << "false!" << endl;
                    cout << NJets << endl;
                    cout << "" << endl;
                }
            }
        }
        
        for (unsigned int iLepESUp = 0; iLepESUp < subSampVec->size(); ++iLepESUp) {
            S_Current = subSampVec->at(iLepESUp);
            subSampBool_LepESUp[S_Current] = false;
            if (!doEvent_LepESUp) continue;
            if (!(S_Current.whichdiLepType < 0 || Type_LepESUp == S_Current.whichdiLepType)) continue;
            if (!(NJets >= S_Current.cutNJets)) continue;
            if (!(NBtagJets >= S_Current.cutNBJets)) continue;
            if (!(Type_LepESUp == 2 && S_Current.histNameSuffix.Contains("FullCut"))) {
                if (MET_LepESUp < S_Current.cutMET) continue;
                if (!(S_Current.doZVeto < 0 || ZVeto_LepESUp == S_Current.doZVeto)) continue;
            }
            if (S_Current.histNameSuffix.Contains("BothinBarrel")) {
                if (!(TMath::Abs(Lep0Vec_LepESUp.Eta()) < barrelEtaEnd && TMath::Abs(Lep1Vec_LepESUp.Eta()) < barrelEtaEnd)) continue;
            }
            if (S_Current.histNameSuffix.Contains("OneinBarrel")) {            
                if (!(TMath::Abs(Lep0Vec_LepESUp.Eta()) > barrelEtaEnd || TMath::Abs(Lep1Vec_LepESUp.Eta()) > barrelEtaEnd)) continue;
                if (!(TMath::Abs(Lep0Vec_LepESUp.Eta()) < barrelEtaEnd || TMath::Abs(Lep1Vec_LepESUp.Eta()) < barrelEtaEnd)) continue;                                        
            }   
            if (S_Current.histNameSuffix.Contains("BothinEndcap")) {
                if (!(TMath::Abs(Lep0Vec_LepESUp.Eta()) > endcapEtaStart && TMath::Abs(Lep1Vec_LepESUp.Eta()) > endcapEtaStart)) continue;
            }
            if (S_Current.histNameSuffix.Contains("0BJets")) {
                if (NBtagJets > 0) {
                    if (!(S_Current.histNameSuffix.Contains("inZMass") && !ZVeto_LepESUp)) continue;
                }
            }
            if (S_Current.histNameSuffix.Contains("_0Jets") && NJets != 0) continue;
            if (S_Current.histNameSuffix.Contains("_1Jet") && NJets != 1) continue;
            if (S_Current.histNameSuffix.Contains("FullCutBlind")) {
                if (MT2ll > 80) continue;
            }
            subSampBool_LepESUp[S_Current] = true;
        }
        
        for (unsigned int iLepESDown = 0; iLepESDown < subSampVec->size(); ++iLepESDown) {
            S_Current = subSampVec->at(iLepESDown);
            subSampBool_LepESDown[S_Current] = false;
            if (!doEvent_LepESDown) continue;
            if (!(S_Current.whichdiLepType < 0 || Type_LepESDown == S_Current.whichdiLepType)) continue;
            if (!(NJets >= S_Current.cutNJets)) continue;
            if (!(NBtagJets >= S_Current.cutNBJets)) continue;
            if (!(Type_LepESDown == 2 && S_Current.histNameSuffix.Contains("FullCut"))) {
                if (MET_LepESDown < S_Current.cutMET) continue;
                if (!(S_Current.doZVeto < 0 || ZVeto_LepESDown == S_Current.doZVeto)) continue;
            }
            if (S_Current.histNameSuffix.Contains("BothinBarrel")) {
                if (!(TMath::Abs(Lep0Vec_LepESDown.Eta()) < barrelEtaEnd && TMath::Abs(Lep1Vec_LepESDown.Eta()) < barrelEtaEnd)) continue;
            }
            if (S_Current.histNameSuffix.Contains("OneinBarrel")) {            
                if (!(TMath::Abs(Lep0Vec_LepESDown.Eta()) > barrelEtaEnd || TMath::Abs(Lep1Vec_LepESDown.Eta()) > barrelEtaEnd)) continue;
                if (!(TMath::Abs(Lep0Vec_LepESDown.Eta()) < barrelEtaEnd || TMath::Abs(Lep1Vec_LepESDown.Eta()) < barrelEtaEnd)) continue;                                        
            }   
            if (S_Current.histNameSuffix.Contains("BothinEndcap")) {
                if (!(TMath::Abs(Lep0Vec_LepESDown.Eta()) > endcapEtaStart && TMath::Abs(Lep1Vec_LepESDown.Eta()) > endcapEtaStart)) continue;
            }
            if (S_Current.histNameSuffix.Contains("0BJets")) {
                if (NBtagJets > 0) {
                    if (!(S_Current.histNameSuffix.Contains("inZMass") && !ZVeto_LepESDown)) continue;
                }
            }
            if (S_Current.histNameSuffix.Contains("_0Jets") && NJets != 0) continue;
            if (S_Current.histNameSuffix.Contains("_1Jet") && NJets != 1) continue;
            if (S_Current.histNameSuffix.Contains("FullCutBlind")) {
                if (MT2ll > 80) continue;
            }
            subSampBool_LepESDown[S_Current] = true;
        }
        
        
        for (unsigned int iJetESUp = 0; iJetESUp < subSampVec->size(); ++iJetESUp) {
            S_Current = subSampVec->at(iJetESUp);
            subSampBool_JetESUp[S_Current] = false;
            if (!doEvent) continue;
            if (!(S_Current.whichdiLepType < 0 || Type == S_Current.whichdiLepType)) continue;
            if (!(NJets_JetESUp >= S_Current.cutNJets)) continue;
            if (!(NBtagJets_JetESUp >= S_Current.cutNBJets)) continue;
            if (!(Type == 2 && S_Current.histNameSuffix.Contains("FullCut"))) {
                if (MET_JetESUp < S_Current.cutMET) continue;
                if (!(S_Current.doZVeto < 0 || ZVeto == S_Current.doZVeto)) continue;
            }
            if (S_Current.histNameSuffix.Contains("BothinBarrel")) {
                if (!(TMath::Abs(Lep0Vec.Eta()) < barrelEtaEnd && TMath::Abs(Lep1Vec.Eta()) < barrelEtaEnd)) continue;
            }
            if (S_Current.histNameSuffix.Contains("OneinBarrel")) {            
                if (!(TMath::Abs(Lep0Vec.Eta()) > barrelEtaEnd || TMath::Abs(Lep1Vec.Eta()) > barrelEtaEnd)) continue;
                if (!(TMath::Abs(Lep0Vec.Eta()) < barrelEtaEnd || TMath::Abs(Lep1Vec.Eta()) < barrelEtaEnd)) continue;                                        
            }   
            if (S_Current.histNameSuffix.Contains("BothinEndcap")) {
                if (!(TMath::Abs(Lep0Vec.Eta()) > endcapEtaStart && TMath::Abs(Lep1Vec.Eta()) > endcapEtaStart)) continue;
            }
            if (S_Current.histNameSuffix.Contains("0BJets")) {
                if (NBtagJets_JetESUp > 0) {
                    if (!(S_Current.histNameSuffix.Contains("inZMass") && !ZVeto)) continue;
                }
            }
            if (S_Current.histNameSuffix.Contains("_0Jets") && NJets_JetESUp != 0) continue;
            if (S_Current.histNameSuffix.Contains("_1Jet") && NJets_JetESUp != 1) continue;
            if (S_Current.histNameSuffix.Contains("FullCutBlind")) {
                if (MT2ll > 80) continue;
            }
            subSampBool_JetESUp[S_Current] = true;
        }
        
        for (unsigned int iJetESDown = 0; iJetESDown < subSampVec->size(); ++iJetESDown) {
            S_Current = subSampVec->at(iJetESDown);
            subSampBool_JetESDown[S_Current] = false;
            if (!doEvent) continue;
            if (!(S_Current.whichdiLepType < 0 || Type == S_Current.whichdiLepType)) continue;
            if (!(NJets_JetESDown >= S_Current.cutNJets)) continue;
            if (!(NBtagJets_JetESDown >= S_Current.cutNBJets)) continue;
            if (!(Type == 2 && S_Current.histNameSuffix.Contains("FullCut"))) {
                if (MET_JetESDown < S_Current.cutMET) continue;
                if (!(S_Current.doZVeto < 0 || ZVeto == S_Current.doZVeto)) continue;
            }
            if (S_Current.histNameSuffix.Contains("BothinBarrel")) {
                if (!(TMath::Abs(Lep0Vec.Eta()) < barrelEtaEnd && TMath::Abs(Lep1Vec.Eta()) < barrelEtaEnd)) continue;
            }
            if (S_Current.histNameSuffix.Contains("OneinBarrel")) {            
                if (!(TMath::Abs(Lep0Vec.Eta()) > barrelEtaEnd || TMath::Abs(Lep1Vec.Eta()) > barrelEtaEnd)) continue;
                if (!(TMath::Abs(Lep0Vec.Eta()) < barrelEtaEnd || TMath::Abs(Lep1Vec.Eta()) < barrelEtaEnd)) continue;                                        
            }   
            if (S_Current.histNameSuffix.Contains("BothinEndcap")) {
                if (!(TMath::Abs(Lep0Vec.Eta()) > endcapEtaStart && TMath::Abs(Lep1Vec.Eta()) > endcapEtaStart)) continue;
            }
            if (S_Current.histNameSuffix.Contains("0BJets")) {
                if (NBtagJets_JetESDown > 0) {
                    if (!(S_Current.histNameSuffix.Contains("inZMass") && !ZVeto)) continue;
                }
            }
            if (S_Current.histNameSuffix.Contains("_0Jets") && NJets_JetESDown != 0) continue;
            if (S_Current.histNameSuffix.Contains("_1Jet") && NJets_JetESDown != 1) continue;
            if (S_Current.histNameSuffix.Contains("FullCutBlind")) {
                if (MT2ll > 80) continue;
            }
            subSampBool_JetESDown[S_Current] = true;
        }
        
        
        /*#######################
         FILL PLOTS
         #####################*/   
        
        /// Set genTop weights to be appropriate values
        preNVtxRWweight_GenTopReweight *= preNVtxRWweight;
        weight_GenTopReweight *= weight;
        //        cout << "weight " << weight << endl;
        //        cout << "weight_GenTopReweight  " << weight_GenTopReweight << endl;
        weightSwitch = weight;
        weight = weight_GenTopReweight;
        weight_GenTopReweight = weightSwitch;
        //        cout << "weight post switch " << weight << endl;
        //        cout << "weight_GenTopReweight post switch " << weight_GenTopReweight << endl;
        weightSwitch = preNVtxRWweight;
        preNVtxRWweight = preNVtxRWweight_GenTopReweight;
        preNVtxRWweight_GenTopReweight = weightSwitch;
        
        if (!doData) {
            weight_LepEffSFUp = weight;
            weight_LepESUp = weight;
            //            weight_LepEffSFUp_LepESUp = weight_LepESUp;
            //            weight_LepEffSFUp_LepESDown = weight_LepESDown;
            weight_LepEffSFDown = weight;
            weight_LepESDown = weight;
            //            weight_LepEffSFDown_LepESUp = weight_LepESUp;
            //            weight_LepEffSFDown_LepESDown = weight_LepESDown;
            preNVtxRWweight_LepEffSFUp = preNVtxRWweight;
            preNVtxRWweight_LepESUp = preNVtxRWweight;
            //            preNVtxRWweight_LepEffSFUp_LepESUp = preNVtxRWweight_LepESUp;
            //            preNVtxRWweight_LepEffSFUp_LepESDown = preNVtxRWweight_LepESDown;
            preNVtxRWweight_LepEffSFDown = preNVtxRWweight;
            preNVtxRWweight_LepESDown = preNVtxRWweight;
            //            preNVtxRWweight_LepEffSFDown_LepESUp = preNVtxRWweight_LepESUp;
            //            preNVtxRWweight_LepEffSFDown_LepESDown = preNVtxRWweight_LepESDown;
            
            weight_LepEffSFUp /= ScaleFactorMC(Type, 0);
            weight_LepEffSFDown /= ScaleFactorMC(Type, 0);
            weight_LepEffSFUp *= ScaleFactorMC(Type, 1);
            weight_LepEffSFDown *= ScaleFactorMC(Type, -1);
            if (!doEvent_LepESUp) {
                weight_LepESUp = 0;            
                preNVtxRWweight_LepESUp = 0;
            }
            else {
                weight_LepESUp /= ScaleFactorMC(Type, 0);    
                weight_LepESUp *= ScaleFactorMC(Type_LepESUp, 0);
                preNVtxRWweight_LepESUp /= ScaleFactorMC(Type, 0);
                preNVtxRWweight_LepESUp *= ScaleFactorMC(Type_LepESUp, 0);
            }
            
            if (!doEvent_LepESDown) {
                weight_LepESDown = 0;
                preNVtxRWweight_LepESDown = 0;
            }
            else {
                weight_LepESDown /= ScaleFactorMC(Type, 0);    
                weight_LepESDown *= ScaleFactorMC(Type_LepESDown, 0);
                preNVtxRWweight_LepESDown /= ScaleFactorMC(Type, 0);
                preNVtxRWweight_LepESDown *= ScaleFactorMC(Type_LepESDown, 0);    
            }
            //            weight_LepEffSFUp_LepESUp /= ScaleFactorMC(Type_LepESUp, 0);
            //            weight_LepEffSFUp_LepESUp *= ScaleFactorMC(Type_LepESUp, 1);
            //            weight_LepEffSFUp_LepESDown /= ScaleFactorMC(Type_LepESDown, 0);
            //            weight_LepEffSFDown_LepESUp /= ScaleFactorMC(Type_LepESUp, 0);
            //            weight_LepEffSFDown_LepESDown /= ScaleFactorMC(Type_LepESDown, 0);
            //stopped with this halfway
            
            if (doVerbosity) cout << "weight post scale " << weight << endl;
            preNVtxRWweight_LepEffSFUp /= ScaleFactorMC(Type, 0);
            preNVtxRWweight_LepEffSFDown /= ScaleFactorMC(Type, 0);
            preNVtxRWweight_LepEffSFUp *= ScaleFactorMC(Type, 1);
            preNVtxRWweight_LepEffSFDown *= ScaleFactorMC(Type, -1);
        }
        
        
        
        /// Set StopWeights to be appropriate values
        if (isSignal) {
            weight_genStopXSecUp = weight;
            weight_genStopXSecDown = weight;
            weight *= stopWeight;
            weight_LepEffSFUp *= stopWeight;
            weight_LepEffSFDown *= stopWeight;
            weight_LepESUp *= stopWeight;
            weight_LepESDown *= stopWeight;
            weight_genStopXSecUp *= stopWeightPlusErr;
            weight_genStopXSecDown *= stopWeightMinusErr;
            /*
             cout << "weight " << weight << endl;
             cout << "weightup " << weight_genStopXSecUp << endl;
             cout << "weightdown " << weight_genStopXSecDown << endl;
             */
            preNVtxRWweight_genStopXSecUp = preNVtxRWweight;
            preNVtxRWweight_genStopXSecDown = preNVtxRWweight;
            preNVtxRWweight *= stopWeight;
            preNVtxRWweight_LepEffSFUp *= stopWeight;
            preNVtxRWweight_LepEffSFDown *= stopWeight;
            preNVtxRWweight_LepESUp *= stopWeight;
            preNVtxRWweight_LepESDown *= stopWeight;
            preNVtxRWweight_genStopXSecUp *= stopWeightPlusErr;
            preNVtxRWweight_genStopXSecDown *= stopWeightMinusErr;
        }
        //require MT2 < 80 for data (for now, I'm allowing it for MC)           
        for (unsigned int i = 0; i < subSampVec->size(); ++i) {
            S_Current = subSampVec->at(i);
            if (subSampBool[S_Current] == true) {
                for (unsigned int j = 0; j < histVec_1D->size(); ++j) {
                    H_Current = histVec_1D->at(j);
                    xIter = stringKeyToVar.find(H_Current.xVarKey);
                    if (doVerbosity) {
                        cout << "" << endl;
                        cout << "ievt " << ievt << endl;
                        cout << "i " << i << endl;
                        cout << "j " << j << endl;
                        cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                        cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                    }
                    if (xIter != stringKeyToVar.end()) {
                        if (doVerbosity) {
                            cout << "xIter first " << xIter->first << endl;
                            cout << "xIter second " << xIter->second << endl;
                        }
                        ///Some necessary continue checks
                        //                        if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                        if (S_Current.blindDataChannel && TString(H_Current.xVarKey).Contains("MT2ll")) {
                            if (blindData && doData && MT2ll > MT2llCut) continue;
                        }
                        if (S_Current.blindDataChannel && TString(H_Current.xVarKey).Contains("MT2lb")) {
                            if (blindData && doData && MT2lb > MT2lbCut) continue;
                        }
                        if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                            if (dPhi(diLepPhi, METPhi) > 1./3. * PI) continue;
                        }
                        else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                            if (dPhi(diLepPhi, METPhi) < 1./3. * PI || dPhi(diLepPhi, METPhi) > 2./3. * PI) continue;                         
                        }                        
                        else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                            if (dPhi(diLepPhi, METPhi) < 2./3. * PI) continue;
                        }                        
                        if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                            if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 1./3. * PI) continue;
                        }
                        else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                            if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 1./3. * PI || dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 2./3. * PI) continue;                         
                        }                        
                        else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                            if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 2./3. * PI) continue;
                        }
                        
                        if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                            if (DeltaPhiMT2lb_JetsUsed > 1./3. * PI) continue;
                        }
                        else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                            if (DeltaPhiMT2lb_JetsUsed < 1./3. * PI || DeltaPhiMT2lb_JetsUsed > 2./3. * PI) continue;                         
                        }                        
                        else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                            if (DeltaPhiMT2lb_JetsUsed < 2./3. * PI) continue;
                        }
                        
                        
                        if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                            if (DeltaPhiMT2lb_BLepsUsed > 1./3. * PI) continue;
                        }
                        else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                            if (DeltaPhiMT2lb_BLepsUsed < 1./3. * PI || DeltaPhiMT2lb_BLepsUsed > 2./3. * PI) continue;                         
                        }                        
                        else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                            if (DeltaPhiMT2lb_BLepsUsed < 2./3. * PI) continue;
                        }
                        fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                        if (H_Current.name.Contains("h_ChannelCutFlow")) fillWeight = 1.;
                        histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                    }
                    if (doVerbosity) cout << "" << endl;
                }            
                for (unsigned int k = 0; k < histVec_2D->size(); ++k) {
                    H_Current = histVec_2D->at(k);
                    xIter = stringKeyToVar.find(H_Current.xVarKey);
                    yIter = stringKeyToVar.find(H_Current.yVarKey);
                    if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end()) { 
                        //                        if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                        if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2ll") || TString(H_Current.yVarKey).Contains("MT2ll"))) {
                            if (blindData && doData && MT2ll > MT2llCut) continue;
                        }
                        if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2lb") || TString(H_Current.yVarKey).Contains("MT2lb"))) {
                            if (blindData && doData && MT2lb > MT2lbCut) continue;
                        }                    
                        histMap_2D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, weight); //there could be shenanigans with this one
                    }
                }
                
                for (unsigned int m = 0; m < histVec_3D->size(); ++m) {
                    H_Current = histVec_3D->at(m);
                    xIter = stringKeyToVar.find(H_Current.xVarKey);
                    yIter = stringKeyToVar.find(H_Current.yVarKey); 
                    zIter = stringKeyToVar.find(H_Current.zVarKey);
                    if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end() && zIter != stringKeyToVar.end()) {
                        //                        if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                        if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2ll") || TString(H_Current.yVarKey).Contains("MT2ll") || TString(H_Current.zVarKey).Contains("MT2ll"))) {
                            if (blindData && doData && MT2ll > MT2llCut) continue;
                        }
                        if (S_Current.blindDataChannel && (TString(H_Current.xVarKey).Contains("MT2lb") || TString(H_Current.yVarKey).Contains("MT2lb") || TString(H_Current.zVarKey).Contains("MT2lb"))) {
                            if (blindData && doData && MT2lb > MT2lbCut) continue;
                        }
                        if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                            if (H_Current.name.Contains("21to30")) {
                                if ((nVtx > 30 || nVtx < 21)) continue;
                            }
                            else if (H_Current.name.Contains("11to20")) {
                                if ((nVtx > 20 || nVtx < 11)) continue;
                            }
                            else if (H_Current.name.Contains("1to10")) {
                                if ((nVtx > 10 || nVtx < 1)) continue;
                            }
                        }
                        histMap_3D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, weight); //there could be shenanigans with this one   
                    }
                }
                if (!doData) {
                    for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {
                        H_Current = histVec_1D_Syst->at(js);
                        if (H_Current.name.Contains("LepESShift") || H_Current.name.Contains("JetESShift")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "js " << js << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " << xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                                if (dPhi(diLepPhi, METPhi) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                                if (dPhi(diLepPhi, METPhi) < 1./3. * PI || dPhi(diLepPhi, METPhi) > 2./3. * PI) continue;
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                                if (dPhi(diLepPhi, METPhi) < 2./3. * PI) continue;
                            }
                            if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 1./3. * PI || dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 2./3. * PI) continue;
                            }
                            if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                                if (DeltaPhiMT2lb_JetsUsed > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                                if (DeltaPhiMT2lb_JetsUsed < 1./3. * PI || DeltaPhiMT2lb_JetsUsed > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                                if (DeltaPhiMT2lb_JetsUsed < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                                if (DeltaPhiMT2lb_BLepsUsed > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                                if (DeltaPhiMT2lb_BLepsUsed < 1./3. * PI || DeltaPhiMT2lb_BLepsUsed > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                                if (DeltaPhiMT2lb_BLepsUsed < 2./3. * PI) continue;
                            }
                            if (H_Current.name.Contains("LepEffSFShiftUp")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepEffSFUp : weight_LepEffSFUp;
                            }
                            else if (H_Current.name.Contains("LepEffSFShiftDown")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepEffSFDown : weight_LepEffSFDown;
                            }
                            else if (H_Current.name.Contains("genTopRW")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_GenTopReweight : weight_GenTopReweight;
                            }
                            else if (H_Current.name.Contains("genStopXSecShiftUp")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_genStopXSecUp : weight_genStopXSecUp;
                            }
                            else if (H_Current.name.Contains("genStopXSecShiftDown")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_genStopXSecDown : weight_genStopXSecDown;
                            }
                            else {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            }
                            histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    
                    for (unsigned int ks = 0; ks < histVec_2D_Syst->size(); ++ks) {                  
                        H_Current = histVec_2D_Syst->at(ks);
                        if (H_Current.name.Contains("LepESShift") || H_Current.name.Contains("JetESShift")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ks " << ks << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            ///condition for which stuff to run on
                            if (H_Current.name.Contains("LepEffSFShiftUp")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepEffSFUp : weight_LepEffSFUp;
                            }
                            else if (H_Current.name.Contains("LepEffSFShiftDown")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepEffSFDown : weight_LepEffSFDown;
                            }
                            else if (H_Current.name.Contains("genTopRW")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_GenTopReweight : weight_GenTopReweight;
                            }
                            else if (H_Current.name.Contains("genStopXSecShiftUp")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_genStopXSecUp : weight_genStopXSecUp;
                            }
                            else if (H_Current.name.Contains("genStopXSecShiftDown")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_genStopXSecDown : weight_genStopXSecDown;
                            }
                            else {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            }
                            histMap_2D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ls = 0; ls < histVec_3D_Syst->size(); ++ls) {                  
                        H_Current = histVec_3D_Syst->at(ls);
                        if (H_Current.name.Contains("LepESShift") || H_Current.name.Contains("JetESShift")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        zIter = stringKeyToVar.find(H_Current.zVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ls " << ls << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                            cout << "H_Current.zVarKey " << H_Current.zVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end() && zIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;                            
                                cout << "zIter first " <<  zIter->first << endl;
                                cout << "zIter second " << zIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                                if (H_Current.name.Contains("21to30")) {
                                    if ((nVtx > 30 || nVtx < 21)) continue;
                                }
                                else if (H_Current.name.Contains("11to20")) {
                                    if ((nVtx > 20 || nVtx < 11)) continue;
                                }
                                else if (H_Current.name.Contains("1to10")) {
                                    if ((nVtx > 10 || nVtx < 1)) continue;
                                }
                            }
                            if (H_Current.name.Contains("LepEffSFShiftUp")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepEffSFUp : weight_LepEffSFUp;
                            }
                            else if (H_Current.name.Contains("LepEffSFShiftDown")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepEffSFDown : weight_LepEffSFDown;
                            }
                            else if (H_Current.name.Contains("genTopRW")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_GenTopReweight : weight_GenTopReweight;
                            }
                            else if (H_Current.name.Contains("genStopXSecShiftUp")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_genStopXSecUp : weight_genStopXSecUp;
                            }
                            else if (H_Current.name.Contains("genStopXSecShiftDown")) {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_genStopXSecDown : weight_genStopXSecDown;
                            }
                            else {
                                fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            }
                            histMap_3D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                }
            }
            if (!doData) {                
                if (subSampBool_LepESUp[S_Current] == true) {
                    for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {
                        H_Current = histVec_1D_Syst->at(js);
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "js " << js << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                        }
                        if (!H_Current.name.Contains("LepESShiftUp")) continue;
                        if (xIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " << xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                                if (dPhi(diLepPhi_LepESUp, METPhi_LepESUp) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                                if (dPhi(diLepPhi_LepESUp, METPhi_LepESUp) < 1./3. * PI || dPhi(diLepPhi_LepESUp, METPhi_LepESUp) > 2./3. * PI) continue;
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                                if (dPhi(diLepPhi_LepESUp, METPhi_LepESUp) < 2./3. * PI) continue;
                            }
                            if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                                if (dPhi(Lep0Vec_LepESUp.Phi(), Lep1Vec_LepESUp.Phi()) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                                if (dPhi(Lep0Vec_LepESUp.Phi(), Lep1Vec_LepESUp.Phi()) < 1./3. * PI || dPhi(Lep0Vec_LepESUp.Phi(), Lep1Vec_LepESUp.Phi()) > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                                if (dPhi(Lep0Vec_LepESUp.Phi(), Lep1Vec_LepESUp.Phi()) < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                                if (DeltaPhiMT2lb_JetsUsed > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                                if (DeltaPhiMT2lb_JetsUsed < 1./3. * PI || DeltaPhiMT2lb_JetsUsed > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                                if (DeltaPhiMT2lb_JetsUsed < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                                if (DeltaPhiMT2lb_BLepsUsed_LepESUp > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                                if (DeltaPhiMT2lb_BLepsUsed_LepESUp < 1./3. * PI || DeltaPhiMT2lb_BLepsUsed_LepESUp > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                                if (DeltaPhiMT2lb_BLepsUsed_LepESUp < 2./3. * PI) continue;
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepESUp : weight_LepESUp;
                            histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ks = 0; ks < histVec_2D_Syst->size(); ++ks) {                  
                        H_Current = histVec_2D_Syst->at(ks);
                        if (!H_Current.name.Contains("LepESShiftUp")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ks " << ks << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepESUp : weight_LepESUp;
                            histMap_2D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ls = 0; ls < histVec_3D_Syst->size(); ++ls) {                  
                        H_Current = histVec_3D_Syst->at(ls);
                        if (!H_Current.name.Contains("LepESShiftUp")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        zIter = stringKeyToVar.find(H_Current.zVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ls " << ls << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                            cout << "H_Current.zVarKey " << H_Current.zVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end() && zIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;                            
                                cout << "zIter first " <<  zIter->first << endl;
                                cout << "zIter second " << zIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                                if (H_Current.name.Contains("21to30")) {
                                    if ((nVtx > 30 || nVtx < 21)) continue;
                                }
                                else if (H_Current.name.Contains("11to20")) {
                                    if ((nVtx > 20 || nVtx < 11)) continue;
                                }
                                else if (H_Current.name.Contains("1to10")) {
                                    if ((nVtx > 10 || nVtx < 1)) continue;
                                }
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepESUp : weight_LepESUp;
                            histMap_3D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                } 
                if (subSampBool_LepESDown[S_Current] == true) {
                    for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {
                        H_Current = histVec_1D_Syst->at(js);
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "js " << js << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "histname contains LepESShiftDown? " << H_Current.name.Contains("LepESShiftDown") << endl;
                        }
                        if (!H_Current.name.Contains("LepESShiftDown")) continue;
                        if (xIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " << xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                                if (dPhi(diLepPhi_LepESDown, METPhi_LepESDown) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                                if (dPhi(diLepPhi_LepESDown, METPhi_LepESDown) < 1./3. * PI || dPhi(diLepPhi_LepESDown, METPhi_LepESDown) > 2./3. * PI) continue;
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                                if (dPhi(diLepPhi_LepESDown, METPhi_LepESDown) < 2./3. * PI) continue;
                            }
                            if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                                if (dPhi(Lep0Vec_LepESDown.Phi(), Lep1Vec_LepESDown.Phi()) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                                if (dPhi(Lep0Vec_LepESDown.Phi(), Lep1Vec_LepESDown.Phi()) < 1./3. * PI || dPhi(Lep0Vec_LepESDown.Phi(), Lep1Vec_LepESDown.Phi()) > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                                if (dPhi(Lep0Vec_LepESDown.Phi(), Lep1Vec_LepESDown.Phi()) < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                                if (DeltaPhiMT2lb_JetsUsed > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                                if (DeltaPhiMT2lb_JetsUsed < 1./3. * PI || DeltaPhiMT2lb_JetsUsed > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                                if (DeltaPhiMT2lb_JetsUsed < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                                if (DeltaPhiMT2lb_BLepsUsed_LepESDown > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                                if (DeltaPhiMT2lb_BLepsUsed_LepESDown < 1./3. * PI || DeltaPhiMT2lb_BLepsUsed_LepESDown > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                                if (DeltaPhiMT2lb_BLepsUsed_LepESDown < 2./3. * PI) continue;
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepESDown : weight_LepESDown;
                            histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ks = 0; ks < histVec_2D_Syst->size(); ++ks) {                  
                        H_Current = histVec_2D_Syst->at(ks);
                        if (!H_Current.name.Contains("LepESShiftDown")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ks " << ks << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepESDown : weight_LepESDown;
                            histMap_2D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ls = 0; ls < histVec_3D_Syst->size(); ++ls) {                  
                        H_Current = histVec_3D_Syst->at(ls);
                        if (!H_Current.name.Contains("LepESShiftDown")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        zIter = stringKeyToVar.find(H_Current.zVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ls " << ls << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                            cout << "H_Current.zVarKey " << H_Current.zVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end() && zIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;                            
                                cout << "zIter first " <<  zIter->first << endl;
                                cout << "zIter second " << zIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                                if (H_Current.name.Contains("21to30")) {
                                    if ((nVtx > 30 || nVtx < 21)) continue;
                                }
                                else if (H_Current.name.Contains("11to20")) {
                                    if ((nVtx > 20 || nVtx < 11)) continue;
                                }
                                else if (H_Current.name.Contains("1to10")) {
                                    if ((nVtx > 10 || nVtx < 1)) continue;
                                }
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight_LepESDown : weight_LepESDown;
                            histMap_3D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                }
                
                if (subSampBool_JetESUp[S_Current] == true) {
                    for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {
                        H_Current = histVec_1D_Syst->at(js);
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        if (!H_Current.name.Contains("JetESShiftUp")) continue;
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "js " << js << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " << xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                                if (dPhi(diLepPhi, METPhi_JetESUp) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                                if (dPhi(diLepPhi, METPhi_JetESUp) < 1./3. * PI || dPhi(diLepPhi, METPhi_JetESUp) > 2./3. * PI) continue;
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                                if (dPhi(diLepPhi, METPhi_JetESUp) < 2./3. * PI) continue;
                            }
                            if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 1./3. * PI || dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                                if (DeltaPhiMT2lb_JetsUsed_JetESUp > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                                if (DeltaPhiMT2lb_JetsUsed_JetESUp < 1./3. * PI || DeltaPhiMT2lb_JetsUsed_JetESUp > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                                if (DeltaPhiMT2lb_JetsUsed_JetESUp < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                                if (DeltaPhiMT2lb_BLepsUsed_JetESUp > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                                if (DeltaPhiMT2lb_BLepsUsed_JetESUp < 1./3. * PI || DeltaPhiMT2lb_BLepsUsed_JetESUp > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                                if (DeltaPhiMT2lb_BLepsUsed_JetESUp < 2./3. * PI) continue;
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ks = 0; ks < histVec_2D_Syst->size(); ++ks) {                  
                        H_Current = histVec_2D_Syst->at(ks);
                        if (!H_Current.name.Contains("JetESShiftUp")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ks " << ks << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            histMap_2D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ls = 0; ls < histVec_3D_Syst->size(); ++ls) {                  
                        H_Current = histVec_3D_Syst->at(ls);
                        if (!H_Current.name.Contains("JetESShiftUp")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        zIter = stringKeyToVar.find(H_Current.zVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ls " << ls << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                            cout << "H_Current.zVarKey " << H_Current.zVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end() && zIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;                            
                                cout << "zIter first " <<  zIter->first << endl;
                                cout << "zIter second " << zIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                                if (H_Current.name.Contains("21to30")) {
                                    if ((nVtx > 30 || nVtx < 21)) continue;
                                }
                                else if (H_Current.name.Contains("11to20")) {
                                    if ((nVtx > 20 || nVtx < 11)) continue;
                                }
                                else if (H_Current.name.Contains("1to10")) {
                                    if ((nVtx > 10 || nVtx < 1)) continue;
                                }
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            histMap_3D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                }
                
                
                if (subSampBool_JetESDown[S_Current] == true) {
                    for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {
                        H_Current = histVec_1D_Syst->at(js);
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        if (!H_Current.name.Contains("JetESShiftDown")) continue;
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "js " << js << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " << xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (TString(H_Current.name).Contains("MT2ll_DPhiZMETClose")) {
                                if (dPhi(diLepPhi, METPhi_JetESDown) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETMid")) {
                                if (dPhi(diLepPhi, METPhi_JetESDown) < 1./3. * PI || dPhi(diLepPhi, METPhi_JetESDown) > 2./3. * PI) continue;
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiZMETFar")) {
                                if (dPhi(diLepPhi, METPhi_JetESDown) < 2./3. * PI) continue;
                            }
                            if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Close")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 1./3. * PI || dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                                if (dPhi(Lep0Vec.Phi(), Lep1Vec.Phi()) < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Close")) {
                                if (DeltaPhiMT2lb_JetsUsed_JetESDown > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Mid")) {
                                if (DeltaPhiMT2lb_JetsUsed_JetESDown < 1./3. * PI || DeltaPhiMT2lb_JetsUsed_JetESDown > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiJet0Jet1Far")) {
                                if (DeltaPhiMT2lb_JetsUsed_JetESDown < 2./3. * PI) continue;
                            }                        
                            if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Close")) {
                                if (DeltaPhiMT2lb_BLepsUsed_JetESDown > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Mid")) {
                                if (DeltaPhiMT2lb_BLepsUsed_JetESDown < 1./3. * PI || DeltaPhiMT2lb_BLepsUsed_JetESDown > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2lb_DPhiBLep0BLep1Far")) {
                                if (DeltaPhiMT2lb_BLepsUsed_JetESDown < 2./3. * PI) continue;
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ks = 0; ks < histVec_2D_Syst->size(); ++ks) {                  
                        H_Current = histVec_2D_Syst->at(ks);
                        if (!H_Current.name.Contains("JetESShiftDown")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ks " << ks << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            histMap_2D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                    for (unsigned int ls = 0; ls < histVec_3D_Syst->size(); ++ls) {                  
                        H_Current = histVec_3D_Syst->at(ls);
                        if (!H_Current.name.Contains("JetESShiftDown")) continue;
                        xIter = stringKeyToVar.find(H_Current.xVarKey);
                        yIter = stringKeyToVar.find(H_Current.yVarKey);
                        zIter = stringKeyToVar.find(H_Current.zVarKey);
                        if (doVerbosity) {
                            cout << "" << endl;
                            cout << "ievt " << ievt << endl;
                            cout << "i " << i << endl;
                            cout << "ls " << ls << endl;
                            cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                            cout << "H_Current.name " << H_Current.name << endl;
                            cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                            cout << "H_Current.yVarKey " << H_Current.yVarKey << endl;
                            cout << "H_Current.zVarKey " << H_Current.zVarKey << endl;
                        }
                        if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end() && zIter != stringKeyToVar.end()) {
                            if (doVerbosity) {
                                cout << "xIter first " <<  xIter->first << endl;
                                cout << "xIter second " << xIter->second << endl;
                                cout << "yIter first " <<  yIter->first << endl;
                                cout << "yIter second " << yIter->second << endl;                            
                                cout << "zIter first " <<  zIter->first << endl;
                                cout << "zIter second " << zIter->second << endl;
                            }
                            ///Some necessary continue checks
                            //                            if (TString(H_Current.xVarKey).Contains("MT2lb") && NJets < 2) continue;
                            if (H_Current.name.Contains("h_MT2ll_vs_DeltaPhiZMET_vs_NJets_nVtx")) {
                                if (H_Current.name.Contains("21to30")) {
                                    if ((nVtx > 30 || nVtx < 21)) continue;
                                }
                                else if (H_Current.name.Contains("11to20")) {
                                    if ((nVtx > 20 || nVtx < 11)) continue;
                                }
                                else if (H_Current.name.Contains("1to10")) {
                                    if ((nVtx > 10 || nVtx < 1)) continue;
                                }
                            }
                            fillWeight = H_Current.name.Contains("preRW") ? preNVtxRWweight : weight;
                            histMap_3D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, fillWeight); //there could be shenanigans with this one
                        }
                        if (doVerbosity) cout << "" << endl;                    
                    }
                }
            }
        }        
    }
    cout << "All events done" << endl;
    outputFile->cd();
    cout << "cd-ing to output directory" << endl;
    TH1F * h_numParFiles = new TH1F("h_numParFiles", "", 2, -0.5, 1.5);
    h_numParFiles->Fill(1);
    h_numParFiles->SetEntries(1);
    /*
     if (endBreakPoint > 0) {
     h_numParFiles->SetEntries(1);
     }
     else {
     h_numParFiles->SetEntries(0);
     }
     */
    h_numParFiles->Write();
    h_eventCount->Write();
    h_CutFlow->Write();
    h_ElecCharIso->Write();
    h_ElecNeutIso->Write();
    h_ElecPhotIso->Write();
    h_ElecRelIso->Write();
    if (fInName.Contains("TT") || fInName.Contains("ttbarsignal")) {
        h_genTopPt->Write();
        h_genTopPtRW->Write();
        h_genAntiTopPt->Write();
        h_genAntiTopPtRW->Write();
    }
    outputFile->Write();
    cout << "Writing of output file done" << endl;
    outputFile->Close();
    cout << "end of code" << endl;
}
StopXSec getCrossSectionStop(float stopMass){
    
    // taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVstopsbottom                                                                  
    // values are in pb                                                                                                                                             
    
    float xsecStopPairProd=0., xsecUncert=0.;
    
    if(stopMass==100) {xsecStopPairProd = 559.757;       xsecUncert = 16.1085/100.;}
    if(stopMass==105) {xsecStopPairProd = 448.456;       xsecUncert = 15.9732/100.;}
    if(stopMass==110) {xsecStopPairProd = 361.917;       xsecUncert = 16.1134/100.;}
    if(stopMass==115) {xsecStopPairProd = 293.281;       xsecUncert = 15.9763/100.;}
    if(stopMass==120) {xsecStopPairProd = 240.077;       xsecUncert = 15.9212/100.;}
    if(stopMass==125) {xsecStopPairProd = 197.122;       xsecUncert = 15.7303/100.;}
    if(stopMass==130) {xsecStopPairProd = 163.376;       xsecUncert = 15.8101/100.;}
    if(stopMass==135) {xsecStopPairProd = 135.791;       xsecUncert = 15.8086/100.;}
    if(stopMass==140) {xsecStopPairProd = 113.319;       xsecUncert = 15.7234/100.;}
    if(stopMass==145) {xsecStopPairProd =  95.0292;      xsecUncert = 15.649/100.;}
    if(stopMass==150) {xsecStopPairProd =  80.268;       xsecUncert = 15.5946/100.;}
    if(stopMass==155) {xsecStopPairProd =  68.0456;      xsecUncert = 15.5232/100.;}
    if(stopMass==160) {xsecStopPairProd =  58.01;        xsecUncert = 15.3899/100.;}
    if(stopMass==165) {xsecStopPairProd =  49.6639;      xsecUncert = 15.3711/100.;}
    if(stopMass==170) {xsecStopPairProd =  42.6441;      xsecUncert = 15.3017/100.;}
    if(stopMass==175) {xsecStopPairProd =  36.7994;      xsecUncert = 15.1749/100.;}
    if(stopMass==180) {xsecStopPairProd =  31.8695;      xsecUncert = 15.2449/100.;}
    if(stopMass==185) {xsecStopPairProd =  27.7028;      xsecUncert = 15.063/100.;}
    if(stopMass==190) {xsecStopPairProd =  24.1585;      xsecUncert = 15.16/100.;}
    if(stopMass==195) {xsecStopPairProd =  21.1597;      xsecUncert = 14.9422/100.;}
    if(stopMass==200) {xsecStopPairProd =  18.5245;      xsecUncert = 14.9147/100.;}
    if(stopMass==205) {xsecStopPairProd =  16.2439;      xsecUncert = 15.117/100.;}
    if(stopMass==210) {xsecStopPairProd =  14.3201;      xsecUncert = 14.8495/100.;}
    if(stopMass==215) {xsecStopPairProd =  12.6497;      xsecUncert = 14.8689/100.;}
    if(stopMass==220) {xsecStopPairProd =  11.1808;      xsecUncert = 14.9108/100.;}
    if(stopMass==225) {xsecStopPairProd =   9.90959;     xsecUncert = 14.9662/100.;}
    if(stopMass==230) {xsecStopPairProd =   8.78125;     xsecUncert = 14.796/100.;}
    if(stopMass==235) {xsecStopPairProd =   7.81646;     xsecUncert = 14.7983/100.;}
    if(stopMass==240) {xsecStopPairProd =   6.96892;     xsecUncert = 14.7878/100.;}
    if(stopMass==245) {xsecStopPairProd =   6.22701;     xsecUncert = 14.7897/100.;}
    if(stopMass==250) {xsecStopPairProd =   5.57596;     xsecUncert = 14.7529/100.;}
    if(stopMass==255) {xsecStopPairProd =   5.00108;     xsecUncert = 14.729/100.;}
    if(stopMass==260) {xsecStopPairProd =   4.48773;     xsecUncert = 14.6782/100.;}
    if(stopMass==265) {xsecStopPairProd =   4.03416;     xsecUncert = 14.7964/100.;}
    if(stopMass==270) {xsecStopPairProd =   3.63085;     xsecUncert = 14.6565/100.;}
    if(stopMass==275) {xsecStopPairProd =   3.2781;      xsecUncert = 14.7341/100.;}
    if(stopMass==280) {xsecStopPairProd =   2.95613;     xsecUncert = 14.7816/100.;}
    if(stopMass==285) {xsecStopPairProd =   2.67442;     xsecUncert = 14.7661/100.;}
    if(stopMass==290) {xsecStopPairProd =   2.42299;     xsecUncert = 14.6805/100.;}
    if(stopMass==295) {xsecStopPairProd =   2.19684;     xsecUncert = 14.8465/100.;}
    if(stopMass==300) {xsecStopPairProd =   1.99608;     xsecUncert = 14.6905/100.;}
    if(stopMass==305) {xsecStopPairProd =   1.81486;     xsecUncert = 14.4434/100.;}
    if(stopMass==310) {xsecStopPairProd =   1.64956;     xsecUncert = 14.4769/100.;}
    if(stopMass==315) {xsecStopPairProd =   1.50385;     xsecUncert = 14.4549/100.;}
    if(stopMass==320) {xsecStopPairProd =   1.3733;      xsecUncert = 14.7503/100.;}
    if(stopMass==325) {xsecStopPairProd =   1.25277;     xsecUncert = 14.2875/100.;}
    if(stopMass==330) {xsecStopPairProd =   1.14277;     xsecUncert = 14.578/100.;}
    if(stopMass==335) {xsecStopPairProd =   1.04713;     xsecUncert = 14.3659/100.;}
    if(stopMass==340) {xsecStopPairProd =   0.959617;    xsecUncert = 14.3896/100.;}
    if(stopMass==345) {xsecStopPairProd =   0.879793;    xsecUncert = 14.3881/100.;}
    if(stopMass==350) {xsecStopPairProd =   0.807323;    xsecUncert = 14.3597/100.;}
    if(stopMass==355) {xsecStopPairProd =   0.74141;     xsecUncert = 14.368/100.;}
    if(stopMass==360) {xsecStopPairProd =   0.681346;    xsecUncert = 14.3357/100.;}
    if(stopMass==365) {xsecStopPairProd =   0.626913;    xsecUncert = 14.3627/100.;}
    if(stopMass==370) {xsecStopPairProd =   0.576882;    xsecUncert = 14.2712/100.;}
    if(stopMass==375) {xsecStopPairProd =   0.531443;    xsecUncert = 14.266/100.;}
    if(stopMass==380) {xsecStopPairProd =   0.489973;    xsecUncert = 14.3962/100.;}
    if(stopMass==385) {xsecStopPairProd =   0.452072;    xsecUncert = 14.2234/100.;}
    if(stopMass==390) {xsecStopPairProd =   0.4176;      xsecUncert = 14.3166/100.;}
    if(stopMass==395) {xsecStopPairProd =   0.385775;    xsecUncert = 14.3112/100.;}
    if(stopMass==400) {xsecStopPairProd =   0.35683;     xsecUncert = 14.2848/100.;}
    if(stopMass==405) {xsecStopPairProd =   0.329881;    xsecUncert = 14.2072/100.;}
    if(stopMass==410) {xsecStopPairProd =   0.305512;    xsecUncert = 14.2648/100.;}
    if(stopMass==415) {xsecStopPairProd =   0.283519;    xsecUncert = 14.102/100.;}
    if(stopMass==420) {xsecStopPairProd =   0.262683;    xsecUncert = 14.3075/100.;}
    if(stopMass==425) {xsecStopPairProd =   0.243755;    xsecUncert = 14.0504/100.;}
    if(stopMass==430) {xsecStopPairProd =   0.226367;    xsecUncert = 14.0494/100.;}
    if(stopMass==435) {xsecStopPairProd =   0.209966;    xsecUncert = 14.0334/100.;}
    if(stopMass==440) {xsecStopPairProd =   0.195812;    xsecUncert = 14.0772/100.;}
    if(stopMass==445) {xsecStopPairProd =   0.181783;    xsecUncert = 14.1771/100.;}
    if(stopMass==450) {xsecStopPairProd =   0.169668;    xsecUncert = 14.2368/100.;}
    if(stopMass==455) {xsecStopPairProd =   0.158567;    xsecUncert = 14.2609/100.;}
    if(stopMass==460) {xsecStopPairProd =   0.147492;    xsecUncert = 14.4105/100.;}
    if(stopMass==465) {xsecStopPairProd =   0.137392;    xsecUncert = 14.4772/100.;}
    if(stopMass==470) {xsecStopPairProd =   0.128326;    xsecUncert = 14.5144/100.;}
    if(stopMass==475) {xsecStopPairProd =   0.119275;    xsecUncert = 14.6664/100.;}
    if(stopMass==480) {xsecStopPairProd =   0.112241;    xsecUncert = 14.6307/100.;}
    if(stopMass==485) {xsecStopPairProd =   0.104155;    xsecUncert = 14.7581/100.;}
    if(stopMass==490) {xsecStopPairProd =   0.0977878;   xsecUncert = 14.7977/100.;}
    if(stopMass==495) {xsecStopPairProd =   0.091451;    xsecUncert = 14.8963/100.;}
    if(stopMass==500) {xsecStopPairProd =   0.0855847;   xsecUncert = 14.9611/100.;}
    if(stopMass==505) {xsecStopPairProd =   0.0801322;   xsecUncert = 15.0389/100.;}
    if(stopMass==510) {xsecStopPairProd =   0.0751004;   xsecUncert = 15.1402/100.;}
    if(stopMass==515) {xsecStopPairProd =   0.0703432;   xsecUncert = 15.2139/100.;}
    if(stopMass==520) {xsecStopPairProd =   0.0660189;   xsecUncert = 15.3368/100.;}
    if(stopMass==525) {xsecStopPairProd =   0.0618641;   xsecUncert = 15.4135/100.;}
    if(stopMass==530) {xsecStopPairProd =   0.0580348;   xsecUncert = 15.4422/100.;}
    if(stopMass==535) {xsecStopPairProd =   0.0545113;   xsecUncert = 15.5446/100.;}
    if(stopMass==540) {xsecStopPairProd =   0.0511747;   xsecUncert = 15.6283/100.;}
    if(stopMass==545) {xsecStopPairProd =   0.0481537;   xsecUncert = 15.726/100.;}
    if(stopMass==550) {xsecStopPairProd =   0.0452067;   xsecUncert = 15.8177/100.;}
    if(stopMass==555) {xsecStopPairProd =   0.0424781;   xsecUncert = 15.9022/100.;}
    if(stopMass==560) {xsecStopPairProd =   0.0399591;   xsecUncert = 16.0067/100.;}
    if(stopMass==565) {xsecStopPairProd =   0.0376398;   xsecUncert = 16.0367/100.;}
    if(stopMass==570) {xsecStopPairProd =   0.0354242;   xsecUncert = 16.137/100.;}
    if(stopMass==575) {xsecStopPairProd =   0.0333988;   xsecUncert = 16.2132/100.;}
    if(stopMass==580) {xsecStopPairProd =   0.0313654;   xsecUncert = 16.3135/100.;}
    if(stopMass==585) {xsecStopPairProd =   0.0295471;   xsecUncert = 16.4264/100.;}
    if(stopMass==590) {xsecStopPairProd =   0.0279395;   xsecUncert = 16.4546/100.;}
    if(stopMass==595) {xsecStopPairProd =   0.0263263;   xsecUncert = 16.567/100.;}
    if(stopMass==600) {xsecStopPairProd =   0.0248009;   xsecUncert = 16.6406/100.;}
    if(stopMass==605) {xsecStopPairProd =   0.0233806;   xsecUncert = 16.7295/100.;}
    if(stopMass==610) {xsecStopPairProd =   0.0220672;   xsecUncert = 16.8447/100.;}
    if(stopMass==615) {xsecStopPairProd =   0.0208461;   xsecUncert = 16.9276/100.;}
    if(stopMass==620) {xsecStopPairProd =   0.0196331;   xsecUncert = 17.0459/100.;}
    if(stopMass==625) {xsecStopPairProd =   0.0185257;   xsecUncert = 17.0835/100.;}
    if(stopMass==630) {xsecStopPairProd =   0.0175075;   xsecUncert = 17.1478/100.;}
    if(stopMass==635) {xsecStopPairProd =   0.0164955;   xsecUncert = 17.2753/100.;}
    if(stopMass==640) {xsecStopPairProd =   0.0155809;   xsecUncert = 17.3814/100.;}
    if(stopMass==645) {xsecStopPairProd =   0.0147721;   xsecUncert = 17.4885/100.;}
    if(stopMass==650) {xsecStopPairProd =   0.0139566;   xsecUncert = 17.56/100.;}
    if(stopMass==655) {xsecStopPairProd =   0.0132456;   xsecUncert = 17.6129/100.;}
    if(stopMass==660) {xsecStopPairProd =   0.0125393;   xsecUncert = 17.7363/100.;}
    if(stopMass==665) {xsecStopPairProd =   0.0118287;   xsecUncert = 17.7959/100.;}
    if(stopMass==670) {xsecStopPairProd =   0.0112223;   xsecUncert = 17.8974/100.;}
    if(stopMass==675) {xsecStopPairProd =   0.0106123;   xsecUncert = 17.9891/100.;}
    if(stopMass==680) {xsecStopPairProd =   0.0100516;   xsecUncert = 18.0618/100.;}
    if(stopMass==685) {xsecStopPairProd =   0.0095256;   xsecUncert = 18.1714/100.;}
    if(stopMass==690) {xsecStopPairProd =   0.0090306;   xsecUncert = 18.2108/100.;}
    if(stopMass==695) {xsecStopPairProd =   0.00856339;  xsecUncert = 18.3365/100.;}
    if(stopMass==700) {xsecStopPairProd =   0.0081141;   xsecUncert = 18.4146/100.;}
    if(stopMass==705) {xsecStopPairProd =   0.00769525;  xsecUncert = 18.4937/100.;}
    if(stopMass==710) {xsecStopPairProd =   0.00730084;  xsecUncert = 18.6195/100.;}
    if(stopMass==715) {xsecStopPairProd =   0.00692243;  xsecUncert = 18.7005/100.;}
    if(stopMass==720) {xsecStopPairProd =   0.00656729;  xsecUncert = 18.819/100.;}
    if(stopMass==725) {xsecStopPairProd =   0.00623244;  xsecUncert = 18.8796/100.;}
    if(stopMass==730) {xsecStopPairProd =   0.00591771;  xsecUncert = 18.996/100.;}
    if(stopMass==735) {xsecStopPairProd =   0.00561049;  xsecUncert = 19.0787/100.;}
    if(stopMass==740) {xsecStopPairProd =   0.00532605;  xsecUncert = 19.1995/100.;}
    if(stopMass==745) {xsecStopPairProd =   0.00506044;  xsecUncert = 19.2916/100.;}
    if(stopMass==750) {xsecStopPairProd =   0.00480639;  xsecUncert = 19.4088/100.;}
    if(stopMass==755) {xsecStopPairProd =   0.00455979;  xsecUncert = 19.508/100.;}
    if(stopMass==760) {xsecStopPairProd =   0.00433688;  xsecUncert = 19.632/100.;}
    if(stopMass==765) {xsecStopPairProd =   0.00412174;  xsecUncert = 19.7141/100.;}
    if(stopMass==770) {xsecStopPairProd =   0.00391839;  xsecUncert = 19.8299/100.;}
    if(stopMass==775) {xsecStopPairProd =   0.00372717;  xsecUncert = 19.9097/100.;}
    if(stopMass==780) {xsecStopPairProd =   0.00354211;  xsecUncert = 20.0016/100.;}
    if(stopMass==785) {xsecStopPairProd =   0.00336904;  xsecUncert = 20.123/100.;}
    if(stopMass==790) {xsecStopPairProd =   0.00320476;  xsecUncert = 20.2271/100.;}
    if(stopMass==795) {xsecStopPairProd =   0.00304935;  xsecUncert = 20.4479/100.;}
    if(stopMass==800) {xsecStopPairProd =   0.00289588;  xsecUncert = 20.516/100.;}
    if(stopMass==805) {xsecStopPairProd =   0.00275424;  xsecUncert = 20.5444/100.;}
    if(stopMass==810) {xsecStopPairProd =   0.0026184;   xsecUncert = 20.8204/100.;}
    if(stopMass==815) {xsecStopPairProd =   0.00249291;  xsecUncert = 21.0063/100.;}
    if(stopMass==820) {xsecStopPairProd =   0.00237168;  xsecUncert = 21.0865/100.;}
    if(stopMass==825) {xsecStopPairProd =   0.00226163;  xsecUncert = 21.0511/100.;}
    if(stopMass==830) {xsecStopPairProd =   0.00214607;  xsecUncert = 21.3705/100.;}
    if(stopMass==835) {xsecStopPairProd =   0.00204589;  xsecUncert = 21.3026/100.;}
    if(stopMass==840) {xsecStopPairProd =   0.00195172;  xsecUncert = 21.6053/100.;}
    if(stopMass==845) {xsecStopPairProd =   0.0018573;   xsecUncert = 21.8177/100.;}
    if(stopMass==850) {xsecStopPairProd =   0.00176742;  xsecUncert = 21.836/100.;}
    if(stopMass==855) {xsecStopPairProd =   0.00168383;  xsecUncert = 22.1411/100.;}
    if(stopMass==860) {xsecStopPairProd =   0.00160403;  xsecUncert = 22.0506/100.;}
    if(stopMass==865) {xsecStopPairProd =   0.00153063;  xsecUncert = 22.3461/100.;}
    if(stopMass==870) {xsecStopPairProd =   0.00145772;  xsecUncert = 22.6206/100.;}
    if(stopMass==875) {xsecStopPairProd =   0.0013878;   xsecUncert = 22.5422/100.;}
    if(stopMass==880) {xsecStopPairProd =   0.00132077;  xsecUncert = 23.2161/100.;}
    if(stopMass==885) {xsecStopPairProd =   0.00126234;  xsecUncert = 23.1283/100.;}
    if(stopMass==890) {xsecStopPairProd =   0.00120568;  xsecUncert = 23.8404/100.;}
    if(stopMass==895) {xsecStopPairProd =   0.00114627;  xsecUncert = 23.7327/100.;}
    if(stopMass==900) {xsecStopPairProd =   0.00109501;  xsecUncert = 23.9439/100.;}
    if(stopMass==905) {xsecStopPairProd =   0.001044;    xsecUncert = 24.1049/100.;}
    if(stopMass==910) {xsecStopPairProd =   0.000996193; xsecUncert = 24.2789/100.;}
    if(stopMass==915) {xsecStopPairProd =   0.00095071;  xsecUncert = 24.5443/100.;}
    if(stopMass==920) {xsecStopPairProd =   0.000907494; xsecUncert = 24.7597/100.;}
    if(stopMass==925) {xsecStopPairProd =   0.000866391; xsecUncert = 24.877/100.;}
    if(stopMass==930) {xsecStopPairProd =   0.000826533; xsecUncert = 25.0849/100.;}
    if(stopMass==935) {xsecStopPairProd =   0.000789573; xsecUncert = 25.2885/100.;}
    if(stopMass==940) {xsecStopPairProd =   0.000753768; xsecUncert = 25.4768/100.;}
    if(stopMass==945) {xsecStopPairProd =   0.000719675; xsecUncert = 25.6582/100.;}
    if(stopMass==950) {xsecStopPairProd =   0.000687022; xsecUncert = 25.8341/100.;}
    if(stopMass==955) {xsecStopPairProd =   0.000656279; xsecUncert = 26.0372/100.;}
    if(stopMass==960) {xsecStopPairProd =   0.000626876; xsecUncert = 26.2059/100.;}
    if(stopMass==965) {xsecStopPairProd =   0.000598955; xsecUncert = 26.3653/100.;}
    if(stopMass==970) {xsecStopPairProd =   0.000571551; xsecUncert = 26.5169/100.;}
    if(stopMass==975) {xsecStopPairProd =   0.000546728; xsecUncert = 26.7985/100.;}
    if(stopMass==980) {xsecStopPairProd =   0.000522495; xsecUncert = 26.9218/100.;}
    if(stopMass==985) {xsecStopPairProd =   0.000499017; xsecUncert = 27.1036/100.;}
    if(stopMass==990) {xsecStopPairProd =   0.000476255; xsecUncert = 27.3032/100.;}
    if(stopMass==995) {xsecStopPairProd =   0.000455959; xsecUncert = 27.4544/100.;}
    if(stopMass==1000) {xsecStopPairProd =  0.000435488; xsecUncert = 27.6595/100.;}
    StopXSec stopCrossSection;
    stopCrossSection.stopProdXsec = xsecStopPairProd;
    stopCrossSection.stopProdXsecUncert = xsecUncert;
    
    if(xsecStopPairProd==0. || xsecUncert==0.)
        std::cout << "no xsec available for this mass point - choose another one!" << std::endl;
    
    return stopCrossSection;
}