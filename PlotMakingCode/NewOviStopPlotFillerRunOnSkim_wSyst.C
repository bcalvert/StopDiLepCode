
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
#include <fstream>
#include <string>
#include <iostream>
//#include <exception>                                                                                                      
#include <sys/stat.h>

#include "TCut.h"
#include "mt2bisect.h"
//#include "StopDict_ver2.h"
#include "../HeaderFiles/StopFunctionDefinitions_v2.h"

//#include "PileUpMC.h"
#include <iostream>
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
    TFile * inputPUFileMCOvi = new TFile(TString("OneDMCPURWNewOvi.root"));
    TFile * inputPUFileMCDESY = new TFile(TString("OneDMCPURWNewDESY.root"));
    TFile * inputPUFileMCOviToDESY = new TFile(TString("OneDMCPURW_OviToDESY.root"));
    //    TH1F * truePUDistMC = (TH1F*) inputPUFileMC->Get("MCPU");
    //    TFile * inputPUFileData = new TFile("RunABCPUDist.root");
    //    TH1F * truePUDistData = (TH1F*) inputPUFileData->Get("pileup");
    //    float intMCPU = truePUDistMC->
    TH1F * nVtxSFHist = (TH1F*) inputPUFileMCOvi->Get("nVtxSF_preRW");
    TH1F * nVtxSFHistOviToDESY = (TH1F*) inputPUFileMCOviToDESY->Get("normFrac");
    TString fileTreeName;    
    TString fInName;
    TString fOutName;    
    TH1F * eventCount;
    
    /////Event Variables/////////////////////    
    int   Event, Type, nVtx, nVtxTrue;
    float PUWeight;
    float weight, preNVtxRWweight;
    float fillWeight;
    float MET,MET_Phi,METSig, MET_Phi_preCorr;
    float METX, METY, METX_preCorr, METY_preCorr;
    float MT2ll, MT2lb, MT2lbPair1, MT2lbPair2;
    float MT2llCut = 80;
    float MT2lbCut = 172;
    float Lep0Px,Lep0Py,Lep0Pz,Lep0E,Lep1Px,Lep1Py,Lep1Pz,Lep1E;
    int   Lep0PdgId, Lep1PdgId;
    float HT,Btag_j0,Btag_j1;
    float Jet0Px,Jet0Py,Jet0Pz,Jet0E,Jet0Et,Jet1Px,Jet1Py,Jet1Pz,Jet1E,Jet1Et;
    float BtagJet0Px,BtagJet0Py,BtagJet0Pz,BtagJet0E,BtagJet0Et,BtagJet1Px,BtagJet1Py,BtagJet1Pz,BtagJet1E,BtagJet1Et;
    float BDT,BDTDis;
    
    float TGenStopMass0,TGenStopMass1,TGenChi0Mass0,TGenChi0Mass1;
    int NJets, NBtagJets;
    float SysVar;
    
    TLorentzVector Lep0Vec, Lep1Vec, DiLepVec;
    TLorentzVector Jet0Vec, Jet1Vec, DiJetVec, BtagJet0Vec, BtagJet1Vec, DiBJetVec;
    TLorentzVector Lep0Jet0Vec, Lep0Jet1Vec, Lep1Jet0Vec, Lep1Jet1Vec;
    TLorentzVector BLep0Vec, BLep1Vec;
    float diLepInvMass, diLepPt, diLepEta, diLepPhi;
    float diJetInvMass, diJetPt, diJetEta, diJetPhi;
    float diBJetInvMass, diBJetPt, diBJetEta, diBJetPhi;
    
    LV * Lep0Vec_DESY, * Lep1Vec_DESY;
    LV * Jet0Vec_DESY, * Jet1Vec_DESY, * BJet0Vec_DESY, * BJet1Vec_DESY;
    
    /********************************************************************************************************/
    // Systematics Stuff
    float weight_LepEffSFUp, weight_LepEffSFDown;
    float preNVtxRWweight_LepEffSFUp, preNVtxRWweight_LepEffSFDown;
    float MT2ll_ShiftUp, MT2ll_ShiftDown;
    float MET_LepESUp, MET_LepESDown, MET_JetESUp, MET_JetESDown, MET_JetERUp, MET_JetERDown;
    float MET_Phi_LepESUp, MET_Phi_LepESDown, MET_Phi_JetESUp, MET_Phi_JetESDown, MET_Phi_JetERUp, MET_Phi_JetERDown;
    float MT2lb_LepESUp, MT2lb_LepESDown, MT2lb_JetESUp, MT2lb_JetESDown, MT2lb_JetERUp, MT2lb_JetERDown;
    TFile * MT2llSmearFile = new TFile("MT2llSmear.root");
    TH1D * MT2llMeanSmear = (TH1D*) MT2llSmearFile->Get("MT2llSmear");
    float MT2llSmearFactor;
    int MT2llSmearFactorBin;
    float MT2llSystConst = 1.5;
    float MT2llSystSlope = 0.0;
    
    TLorentzVector Lep0Vec_LepESUp, Lep0Vec_LepESDown, Lep0Vec_LepERUp, Lep0Vec_LepERDown;
    TLorentzVector Lep1Vec_LepESUp, Lep1Vec_LepESDown, Lep1Vec_LepERUp, Lep1Vec_LepERDown;
    TLorentzVector DiLepVec_LepESUp, DiLepVec_LepESDown, DiLepVec_LepERUp, DiLepVec_LepERDown;
    float diLepInvMass_LepESUp, diLepInvMass_LepESDown, diLepInvMass_LepERUp, diLepInvMass_LepERDown;
    float diLepPt_LepESUp, diLepPt_LepESDown, diLepPt_LepERUp, diLepPt_LepERDown;
    float diLepEta_LepESUp, diLepEta_LepESDown, diLepEta_LepERUp, diLepEta_LepERDown;
    float diLepPhi_LepESUp, diLepPhi_LepESDown, diLepPhi_LepERUp, diLepPhi_LepERDown;

    TLorentzVector Jet0Vec_JetESUp, Jet0Vec_JetESDown, Jet0Vec_JetERUp, Jet0Vec_JetERDown;
    TLorentzVector Jet1Vec_JetESUp, Jet1Vec_JetESDown, Jet1Vec_JetERUp, Jet1Vec_JetERDown;
    TLorentzVector DiJetVec_JetESUp, DiJetVec_JetESDown, DiJetVec_JetERUp, DiJetVec_JetERDown;
    float diJetInvMass_JetESUp, diJetInvMass_JetESDown, diJetInvMass_JetERUp, diJetInvMass_JetERDown;

    TLorentzVector BtagJet0Vec_JetESUp, BtagJet0Vec_JetESDown, BtagJet0Vec_JetERUp, BtagJet0Vec_JetERDown;
    TLorentzVector BtagJet1Vec_JetESUp, BtagJet1Vec_JetESDown, BtagJet1Vec_JetERUp, BtagJet1Vec_JetERDown;
    TLorentzVector DiBJetVec_JetESUp, DiBJetVec_JetESDown, DiBJetVec_JetERUp, DiBJetVec_JetERDown;
    float diBJetInvMass_JetESUp, diBJetInvMass_JetESDown, diBJetInvMass_JetERUp, diBJetInvMass_JetERDown;

    /********************************************************************************************************/
    
    
    float lumi = 19300; //5296.3; // ipb                                                                                     
    const float genStopMassMin = 295, genStopMassMax = 355, genDeltaM_stopChi0_Min = 195, genDeltaM_stopChi0_Max = 255; 
    // add 5 GeV safety margin (deltaM = 10 GeV in the FineBin sample)  
    float Nevt_stop_oneMassPoint = 50000 * ( (genStopMassMax-genStopMassMin)/10. ) * ( (genDeltaM_stopChi0_Max-genDeltaM_stopChi0_Min)/10. );  
    // 50k evts per point x Npoints
    
    ////input cuts/commands    
    //    const double PI = 3.14159265;
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
        else if (strncmp (argv[k],"ReleaseTheKraken", 16) == 0) {
            blindData = 0;
            cout << "RELEASING THE KRAKEN!!! " << endl;
            cout << "http://www.youtube.com/watch?v=gb2zIR2rvRQ " << endl;
        }
        else if (strncmp (argv[k],"limStats",8) == 0) {
            nEvents = strtol(argv[k+1], NULL, 10);   
        }
    }
        ////input cuts/commands    
    
    char Buffer[500];
    char MyRootFile[2000];
    ifstream * outDirFile;
    TRegexp fCutSlash("[^/]+$");
    if (!grabOutDir) {
        fOutName = "";   
    }
    else {
        outDirFile = new ifstream(TString("outputSavePath.txt"));
        if (!(outDirFile->eof())) {
            outDirFile->getline(Buffer,500);
            fOutName = TString(string(Buffer));        
        }
    }    
            /*
    if (fInName.Contains("tkolberg")) {
        fOutName = fInName;
        fOutName.Replace(12, 8, "bcalvert");   
    }
    else {
        fOutName = fInName;
    }
    */
    fOutName += fInName(fCutSlash);
    if (fInName.Contains("MuEG") || fInName.Contains("DoubleMu") || fInName.Contains("DoubleEl") || fInName.Contains("run2012")) {
        cout << "Running on Data" << endl;
        doData = 1;
    }
    if (doData) fOutName += "MT2Leq80";
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
    if (!blindData && doData) fOutName += "_NOTBLIND";
    fOutName += "_Output.root";
    cout << "saving to " << fOutName << endl;
    TFile * outputFile;
    outputFile = new TFile(fOutName,"RECREATE");
    if (whichNTupleType == 1) {
        fileTreeName = "writeNTuple/NTuple";
        nVtxSFHist = (TH1F*) inputPUFileMCDESY->Get("nVtxSF_preRW");
    }
    else {
        fileTreeName = "OviSkimTree";           
    }
    TChain fileTree(fileTreeName);
    TFile inputFile(fInName + TString(".root"));
    float   BTagCut = 0.679;  //CSV Middle working point, see (remove underscore in address): h_ttps://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    float   JetPtCutforCount = 30;
    bool ZVeto           = true;
    float ZWindowLB      = 76;
    float ZWindowUB      = 106;   
    float barrelEtaEnd = 1.442; float endcapEtaStart = 1.566; 
    //////////DESY NTupleStuff
    
    ///Branch definitions
    /// basic things
    TBranch        *b_runNumber;   //!
    TBranch        *b_lumiBlock;   //!
    TBranch        *b_eventNumber;   //!
    TBranch        *b_vertMulti;   //!
    TBranch        *b_vertMultiTrue;   //!
    
    UInt_t          runNumber;
    UInt_t          lumiBlock;
    UInt_t          eventNumber;
    
    Int_t           vertMulti;
    Int_t           vertMultiTrue;
    /// Lepton branches
    TBranch        *b_lepton;
    TBranch        *b_leptons_;   //!
    TBranch        *b_leptons_fCoordinates_fPt;   //!
    TBranch        *b_leptons_fCoordinates_fEta;   //!
    TBranch        *b_leptons_fCoordinates_fPhi;   //!
    TBranch        *b_leptons_fCoordinates_fM;   //!
    TBranch        *b_lepPdgId;   //!
    TBranch        *b_lepPfIso;
    VLV            *leptons = 0;
    vector<int>    *lepPdgId = 0;
    vector<double> *lepPfIso = 0;
    /// Jet branches
    TBranch        *b_jet;
    TBranch        *b_jets_;   //!
    TBranch        *b_jets_fCoordinates_fPt;   //!
    TBranch        *b_jets_fCoordinates_fEta;   //!
    TBranch        *b_jets_fCoordinates_fPhi;   //!
    TBranch        *b_jets_fCoordinates_fM;   //!
    TBranch        *b_jetBTagCSV;   //!
    VLV            *jets = 0;
    VLV            *jetsForMET = 0;
    vector<double>  *jetBTagCSV = 0;
    //    float   BTagCut = 0.679;  //CSV Middle working point, see (remove underscore in address): h_ttps://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    //    float   JetPtCutforCount = 30;
    ///sub branch for gen jets for met smearing
    /*
     TBranch        *b_allGenJets_;   //!
     TBranch        *b_allGenJets_fCoordinates_fPt;   //!
     TBranch        *b_allGenJets_fCoordinates_fEta;   //!
     TBranch        *b_allGenJets_fCoordinates_fPhi;   //!
     TBranch        *b_allGenJets_fCoordinates_fM;   //!
     TBranch        *b_associatedGenJet_;   //!
     TBranch        *b_associatedGenJet_fCoordinates_fPt;   //!
     TBranch        *b_associatedGenJet_fCoordinates_fEta;   //!
     TBranch        *b_associatedGenJet_fCoordinates_fPhi;   //!
     TBranch        *b_associatedGenJet_fCoordinates_fM;   //!
     TBranch        *b_associatedGenJetForMET_;   //!
     TBranch        *b_associatedGenJetForMET_fCoordinates_fPt;   //!
     TBranch        *b_associatedGenJetForMET_fCoordinates_fEta;   //!
     TBranch        *b_associatedGenJetForMET_fCoordinates_fPhi;   //!
     TBranch        *b_associatedGenJetForMET_fCoordinates_fM;   //!
     VLV            *associatedGenJet;
     VLV            *associatedGenJetForMET;
     */
    /// MET branches
    TBranch        *b_met;
    TBranch        *b_met_fCoordinates_fPt;   //!
    TBranch        *b_met_fCoordinates_fEta;   //!
    TBranch        *b_met_fCoordinates_fPhi;   //!
    TBranch        *b_met_fCoordinates_fM;   //!
    LV             *met = 0;
    
    TBranch        *b_GenWeight;
    double         GenWeight;
    
    //////////////////////////
    
    /////Set up the tree////////    
    //    fileTree.Add(fInName + TString(".root"));
    fileTree.Add(fInName + TString(".root"));
    if (whichNTupleType == 0) {
        fileTree.SetBranchAddress( "TWeight",  &PUWeight );
        fileTree.SetBranchAddress( "TChannel", &Type );
        fileTree.SetBranchAddress( "TNPV",     &nVtx );
        
        fileTree.SetBranchAddress( "TMET",     &MET );
        fileTree.SetBranchAddress( "TMET_Phi", &MET_Phi );
        fileTree.SetBranchAddress( "TMETSig",  &METSig );
        
        fileTree.SetBranchAddress( "TLep0Px", &Lep0Px );
        fileTree.SetBranchAddress( "TLep0Py", &Lep0Py );
        fileTree.SetBranchAddress( "TLep0Pz", &Lep0Pz );
        fileTree.SetBranchAddress( "TLep0E",  &Lep0E );
        fileTree.SetBranchAddress( "TLep0PdgId", &Lep0PdgId );
        
        fileTree.SetBranchAddress( "TLep1Px", &Lep1Px );
        fileTree.SetBranchAddress( "TLep1Py", &Lep1Py );
        fileTree.SetBranchAddress( "TLep1Pz", &Lep1Pz );
        fileTree.SetBranchAddress( "TLep1E",  &Lep1E );
        fileTree.SetBranchAddress( "TLep1PdgId", &Lep1PdgId );
        
        fileTree.SetBranchAddress( "TNJets",    &NJets );
        //        fileTree.SetBranchAddress( "THT",       &HT );
        fileTree.SetBranchAddress( "TNJetsBtag",&NBtagJets );
        //        fileTree.SetBranchAddress( "TBtagJet0", &Btag_j0 );
        //        fileTree.SetBranchAddress( "TBtagJet1", &Btag_j1 );
        
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
        
        fileTree.SetBranchAddress( "TBtagJet1Px", &BtagJet1Px );
        fileTree.SetBranchAddress( "TBtagJet1Py", &BtagJet1Py );
        fileTree.SetBranchAddress( "TBtagJet1Pz", &BtagJet1Pz );
        fileTree.SetBranchAddress( "TBtagJet1E", &BtagJet1E );
        
        if (fileTree.GetBranch("TGenStopMass0")){
            fileTree.SetBranchAddress( "TGenStopMass0", &TGenStopMass0 );
            fileTree.SetBranchAddress( "TGenStopMass1", &TGenStopMass1 );
            fileTree.SetBranchAddress( "TGenChi0Mass0", &TGenChi0Mass0 );
            fileTree.SetBranchAddress( "TGenChi0Mass1", &TGenChi0Mass1 );
        }
        else{
            TGenStopMass0 = -99.;
            TGenStopMass1 = -99.;
            TGenChi0Mass0 = -99.;
            TGenChi0Mass1 = -99.;
        }
        
        if(fInName.Contains("TMVA")) fileTree.SetBranchAddress( "BDT",&BDT);
        if(fInName.Contains("Dis")) fileTree.SetBranchAddress( "BDTDis",&BDTDis);
        
        if(fInName.Contains("Up") || fInName.Contains("Down")) fileTree.SetBranchAddress( "TSysVar", &SysVar );
        else SysVar=1.0;
    }
    else if (whichNTupleType == 1) {
        fileTree.SetBranchAddress("runNumber", &runNumber, &b_runNumber);
        fileTree.SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
        fileTree.SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
        fileTree.SetBranchAddress("vertMulti", &vertMulti, &b_vertMulti);
        fileTree.SetBranchAddress("vertMultiTrue", &vertMultiTrue, &b_vertMultiTrue);
        
        fileTree.SetBranchAddress("lepPdgId", &lepPdgId, &b_lepPdgId );	
        fileTree.SetBranchAddress("lepPfIso", &lepPfIso, &b_lepPfIso );	
        fileTree.SetBranchAddress("jets", &jets, &b_jet );
        fileTree.SetBranchAddress("leptons", &leptons, &b_lepton );
        fileTree.SetBranchAddress("jetBTagCSV", &jetBTagCSV, &b_jetBTagCSV );
        fileTree.SetBranchAddress("met", &met, &b_met );
        fileTree.SetBranchAddress("weightGenerator", &GenWeight, &b_GenWeight);
        SysVar = 1.0;
    }
    ////Book histograms and histogram names
    outputFile->cd();
    vector<HistogramT> * histVec_1D = OneDeeHistTVec();
    vector<HistogramT> * histVec_2D = TwoDeeHistTVec();
    //    vector<HistogramT> * histVec_3D;
    vector<SampleT> * subSampVec    = SubSampVec();    ///Define things necessary for booking histograms
    vector<SystT> * systVec         = SystVec();
    vector<HistogramT> * histVec_1D_Syst = new vector<HistogramT>;
//    vector<HistogramT> * histVec_2D_Syst;
//    vector<HistogramT> * histVec_3D_Syst;
    HistogramT H_Current; 
    TH1D * h_1DCurr; TH2D * h_2DCurr; TH3D * h_3DCurr;
    HMap_1D histMap_1D; HMap_2D histMap_2D; HMap_3D histMap_3D; passCutMap subSampBool;
    SampleT S_Current;
    TString histTitle;
    TString axesTitle;
    float nXBins, nYBins, nZBins;
    float xBinMin, xBinMax;
    float yBinMin, yBinMax;
    float zBinMin, zBinMax;
    if (doBookSyst) {
        histVec_1D_Syst = AddSystHists(histVec_1D, systVec);
//        addSystHists(histVec_2D, systVec);
//        addSystHists(histVec_3D, systVec);
    }
    // cout << "test " << endl;
    for (unsigned int i = 0; i < subSampVec->size(); ++i) {
        S_Current = subSampVec->at(i);
        for (int j = 0; j < (int) histVec_1D->size(); ++j) {            
            H_Current = histVec_1D->at(j);
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
            if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel) {
                xBinMax = xBinMax / 2;
            }
//            if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.
            h_1DCurr = new TH1D(histTitle, axesTitle, nXBins, xBinMin, xBinMax); h_1DCurr->Sumw2();
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
            h_2DCurr = new TH2D(histTitle, axesTitle, nXBins, xBinMin, xBinMax, nYBins, yBinMin, yBinMax); h_2DCurr->Sumw2();
            histMap_2D[histKey(H_Current, S_Current)] = h_2DCurr;
        }
        /*
         for (int l = 0; l < histVec_3D->size(); ++l) {
         H_Current = histVec_3D->at(l);
         }
         */
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
                if (H_Current.xLabel.Contains("MT2_{ll}") && S_Current.blindDataChannel) {
                    xBinMax = xBinMax / 2;
                }
                h_1DCurr = new TH1D(histTitle, axesTitle, nXBins, xBinMin, xBinMax); h_1DCurr->Sumw2();
                histMap_1D[histKey(H_Current, S_Current)] = h_1DCurr;
                //will this cause memory leaks?
            }
        }
    }   
    /////
    TRandom rand;
    map<string, float>::iterator xIter;
    map<string, float>::iterator yIter;
    map<string, float>::iterator zIter;
    
    vector<int> * eventJetParams = new vector<int>;
    int lep0Index, lep1Index;
    int jet0Index, jet1Index;
    bool doEvent;
    // cout << "test2 " << endl;
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
    for (Long64_t ievt=0; ievt<nEvents;ievt++) {
        //    for (Long64_t ievt=0; ievt<1000;ievt++) 
        diBJetPhi = 0; diBJetPt = 0; diBJetEta = 0; diBJetInvMass = 0;
        diJetPhi = 0; diJetPt = 0; diJetEta = 0; diJetInvMass = 0;
        MT2lb = 0; MT2lbPair1 = 0; MT2lbPair2 = 0;
        lep0Index = 0; lep1Index = 1;
        doEvent = true;
        map<string, float> stringKeyToVar;
        fileTree.GetEntry(ievt);
        //        cout << "Lep0 Px " << Lep0Px << endl;
        ZVeto = true;
        //*****************************************                                                                                                                   
        // SELECT THE SIGNAL MASS POINT (if signal)                                                                                                                   
        //*****************************************                                                                                                                   
        weight = 1.;
        // cout << "test 5" << endl;
        if (whichNTupleType == 0) {
            if( (TGenStopMass0 > genStopMassMin && TGenStopMass0 < genStopMassMax &&
                 TGenStopMass1 > genStopMassMin && TGenStopMass1 < genStopMassMax && // Stop mass limits                                                                  
                 (TGenStopMass0-TGenChi0Mass0) > genDeltaM_stopChi0_Min && (TGenStopMass0-TGenChi0Mass0) < genDeltaM_stopChi0_Max &&
                 (TGenStopMass1-TGenChi0Mass1) > genDeltaM_stopChi0_Min && (TGenStopMass1-TGenChi0Mass1) < genDeltaM_stopChi0_Max) || // Neutralino mass limits           
               (TGenStopMass0 == -99. && TGenStopMass1 == -99. &&
                TGenChi0Mass0 == -99. && TGenChi0Mass1 == -99.) ) { // Background or data                                                                                
                   
                   //*****************************************                                                                                                                 
                   // NORMALIZE STOP SIGNAL TO INTEGRATED LUMI                                                                                                                 
                   //*****************************************                                                                                                                 
                   
                   float stopWeight = 0.;
                   if (fileTree.GetBranch("TGenStopMass0")){
                       float roundedStopMass = floor((TGenStopMass0+12.5)/25.)*25.; 
                       // assuming delta_Mstop = 25 GeV, this rounds the stop mass to the closest value divisible by 25                                                                                                                                                              
                       StopXSec theStopXSec = getCrossSectionStop(roundedStopMass);
                       stopWeight = theStopXSec.stopProdXsec * lumi / Nevt_stop_oneMassPoint; 
                       //equivalent to G_Event_Weight= dm->GetCrossSection()* G_Event_Lumi/ dm->GetEventsInTheSample();                                                                                                                                                    
                   }
                   
                   if(fileTree.GetBranch("TGenStopMass0"))
                       weight = PUWeight * stopWeight;
                   else
                       weight = PUWeight;
                   
               }        // cout << "test 6" << endl;
        }
        else if (whichNTupleType == 1 ) {
            weight = (float) GenWeight;
        }
        preNVtxRWweight = weight;
        weight_LepEffSFUp = weight;
        weight_LepEffSFDown = weight;
        preNVtxRWweight_LepEffSFUp = weight;
        preNVtxRWweight_LepEffSFDown = weight;
        //        cout << "weight from GenWeight " << weight << endl;
        ////Calculate all the event variables needed
        if (whichNTupleType == 0) {
            Lep0Vec.SetPxPyPzE(Lep0Px,Lep0Py,Lep0Pz,Lep0E);
            Lep1Vec.SetPxPyPzE(Lep1Px,Lep1Py,Lep1Pz,Lep1E);
            //        cout << "NJets " << NJets << endl;
            if (NJets > 0) {
                Jet0Vec.SetPxPyPzE(Jet0Px,Jet0Py,Jet0Pz,Jet0E);
                if (NJets > 1) {
                    Jet1Vec.SetPxPyPzE(Jet1Px,Jet1Py,Jet1Pz,Jet1E);
                    DiJetVec = Jet0Vec + Jet1Vec;
                    diJetInvMass = DiJetVec.M();
                    //                cout << "diJetInvMass " << diJetInvMass << endl;
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
                cout << "ScaleFactorMC_0 for Type " << Type << " is " << ScaleFactorMC(Type, 0);
                cout << "ScaleFactorMC_+ for Type " << Type << " is " << ScaleFactorMC(Type, +1);
                cout << "ScaleFactorMC_- for Type " << Type << " is " << ScaleFactorMC(Type, -1);
                weight *= ScaleFactorMC(Type, 0);
                weight_LepEffSFUp *= ScaleFactorMC(Type, 1);
                weight_LepEffSFDown *= ScaleFactorMC(Type, -1);
                if (doVerbosity) cout << "weight post scale " << weight << endl;
                preNVtxRWweight *= ScaleFactorMC(Type, 0);
                preNVtxRWweight_LepEffSFUp *= ScaleFactorMC(Type, 1);
                preNVtxRWweight_LepEffSFDown *= ScaleFactorMC(Type, -1);
            }
            //        cout << "weight post doHack " << weight << endl;
            
        }
        else if (whichNTupleType == 1) {
            //            cout << "sub test " << endl;
            LeptonInfo(leptons, lepPdgId, lepPfIso, lep0Index, lep1Index, doEvent, Type, Lep0PdgId, Lep1PdgId);
            if (!doEvent) continue;
            //            cout << "sub test 1" << endl;
            eventJetParams = NumJets(jets, JetPtCutforCount, leptons, lep0Index, lep1Index, jetBTagCSV, HT, BTagCut);
            //            // cout << "test 6a" << endl;
            //            // cout << "test 6b" << endl;
            NJets = eventJetParams->at(0);
            NBtagJets = eventJetParams->at(1);
            MET = met->Pt();
            MET_Phi = met->Phi();
            nVtx = vertMulti;
            nVtxTrue = vertMultiTrue;
            if (doPURW && !doData) {
                weight = PileupRW(OneDPUDist, nVtxTrue);
                if (doHackPURW) weight = PileupRW(nVtxSFHist, nVtx);
            }
            //            cout << "weight from hackPURW " << weight << endl;
            //            cout << "test " << endl;
            Lep0Vec_DESY = &leptons->at(lep0Index);
            Lep1Vec_DESY = &leptons->at(lep1Index);
            //            cout << "test 1 " << endl;
            Lep0Vec.SetPxPyPzE(Lep0Vec_DESY->Px(), Lep0Vec_DESY->Py(), Lep0Vec_DESY->Pz(), Lep0Vec_DESY->E());
            Lep1Vec.SetPxPyPzE(Lep1Vec_DESY->Px(), Lep1Vec_DESY->Py(), Lep1Vec_DESY->Pz(), Lep1Vec_DESY->E());
            if (NJets > 0) {
                //                cout << "test 2 " << endl;
                Jet0Vec_DESY = &jets->at(eventJetParams->at(2));
                Jet0Vec.SetPxPyPzE(Jet0Vec_DESY->Px(), Jet0Vec_DESY->Py(), Jet0Vec_DESY->Pz(), Jet0Vec_DESY->E());
                if (NJets > 1) {
                    //                    cout << "test 3" << endl;
                    Jet1Vec_DESY = &jets->at(eventJetParams->at(3));   
                    Jet1Vec.SetPxPyPzE(Jet1Vec_DESY->Px(), Jet1Vec_DESY->Py(), Jet1Vec_DESY->Pz(), Jet1Vec_DESY->E());
                    DiJetVec = Jet0Vec + Jet1Vec;
                    diJetInvMass = DiJetVec.M();
                    diJetPt = DiJetVec.Pt();
                    diJetEta = DiJetVec.Eta();
                    diJetPhi = DiJetVec.Phi();
                }
                if (NBtagJets > 0) {
                    //                    cout << "test 4" << endl;
                    //                    cout << "eventJetParams->at(4) " << endl;
                    BJet0Vec_DESY = &jets->at(eventJetParams->at(4));
                    BtagJet0Vec.SetPxPyPzE(BJet0Vec_DESY->Px(), BJet0Vec_DESY->Py(), BJet0Vec_DESY->Pz(), BJet0Vec_DESY->E());
                    if (NBtagJets > 1) {                   
                        //                        cout << "test 5" << endl;
                        BJet1Vec_DESY = &jets->at(eventJetParams->at(5));
                        BtagJet1Vec.SetPxPyPzE(BJet1Vec_DESY->Px(), BJet1Vec_DESY->Py(), BJet1Vec_DESY->Pz(), BJet1Vec_DESY->E());
                        DiBJetVec = BtagJet0Vec + BtagJet1Vec;
                        diBJetInvMass = DiBJetVec.M();
                        diBJetPt = DiBJetVec.Pt();
                        diBJetEta = DiBJetVec.Eta();
                        diBJetPhi = DiBJetVec.Phi();
                    }
                }
            }
            if (!doData) {
                if (doVerbosity) {
                    cout << "Type " << Type << endl;
                    cout << "weight pre scale " << weight << endl;
                }
                //                cout << "weight pre scale factor MC " << weight << endl;
                //                cout << "weight scalefactor for Type = " << Type << " is " << ScaleFactorMC(Type) << endl;
//                cout << "ScaleFactorMC_0 for Type " << Type << " is " << ScaleFactorMC(Type, 0) << endl;
//                cout << "ScaleFactorMC_+ for Type " << Type << " is " << ScaleFactorMC(Type, +1) << endl;
//                cout << "ScaleFactorMC_- for Type " << Type << " is " << ScaleFactorMC(Type, -1) << endl;
//                cout << " weight " << weight << endl;
                weight *= ScaleFactorMC(Type, 0);
//                cout << " weight post " << weight << endl;
//                cout << " weight_LepEffSFUp " << weight_LepEffSFUp << endl;
                weight_LepEffSFUp *= ScaleFactorMC(Type, 1);
//                cout << " weight_LepEffSFUp post " << weight_LepEffSFUp << endl;
//                cout << " weight_LepEffSFDown " << weight_LepEffSFDown << endl;
                weight_LepEffSFDown *= ScaleFactorMC(Type, -1);
//                cout << " weight_LepEffSFDown post " << weight_LepEffSFDown << endl;
                if (doVerbosity) cout << "weight post scale " << weight << endl;
                preNVtxRWweight *= ScaleFactorMC(Type, 0);
                preNVtxRWweight_LepEffSFUp *= ScaleFactorMC(Type, 1);
                preNVtxRWweight_LepEffSFDown *= ScaleFactorMC(Type, -1);
                //                cout << "weight post scale factor MC " << weight << endl;
            }
        }            
        // cout << "test 7" << endl;
        if(Type==-2) Type=2;// Fill Histograms MuonElectron=ElectronMuon 
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
        MET_Phi_preCorr = MET_Phi;
        METX_preCorr = MET*TMath::Cos(MET_Phi_preCorr);
        METY_preCorr = MET*TMath::Sin(MET_Phi_preCorr);
        METX = METX_preCorr;
        METY = METY_preCorr;
        if (doPhiCorr) MetPhiCorrect(doData, METX, METY, nVtx);        
        if (doMETSmear) {
            METX *= rand.Gaus(1, METSF * METX);   
            METY *= rand.Gaus(1, METSF * METY);   
        }
        MET_Phi = TMath::ATan2(METY, METX);
        MET = TMath::Sqrt(METX * METX + METY * METY);
        MT2ll=getMT2(Lep0Vec, Lep1Vec, MET, MET_Phi);
        /******************************************************/
        //Systematics//
        Lep0Vec_LepESUp = LeptonScaleSystShift(Lep0Vec, Lep0PdgId, 1.0);
        Lep0Vec_LepESDown = LeptonScaleSystShift(Lep0Vec, Lep0PdgId, -1.0);
        Lep1Vec_LepESUp = LeptonScaleSystShift(Lep1Vec, Lep1PdgId, 1.0);
        Lep1Vec_LepESDown = LeptonScaleSystShift(Lep1Vec, Lep1PdgId, -1.0);
        DiLepVec_LepESUp = Lep0Vec_LepESUp + Lep1Vec_LepESUp;
        diLepInvMass_LepESUp = DiLepVec_LepESUp.M();
        diLepPt_LepESUp = DiLepVec_LepESUp.Pt();
        diLepEta_LepESUp = DiLepVec_LepESUp.Eta();
        diLepPhi_LepESUp = DiLepVec_LepESUp.Phi();
        DiLepVec_LepESDown = Lep0Vec_LepESDown + Lep1Vec_LepESDown;
        diLepInvMass_LepESDown = DiLepVec_LepESDown.M();
        diLepPt_LepESDown = DiLepVec_LepESDown.Pt();
        diLepEta_LepESDown = DiLepVec_LepESDown.Eta();
        diLepPhi_LepESDown = DiLepVec_LepESDown.Phi();
        
        vector<TLorentzVector> * leptonVec = new vector<TLorentzVector>;
        vector<TLorentzVector> * leptonVecEnUp = new vector<TLorentzVector>;
        vector<TLorentzVector> * leptonVecEnDown = new vector<TLorentzVector>;
        leptonVec->push_back(Lep0Vec); leptonVec->push_back(Lep1Vec);
        leptonVecEnUp->push_back(Lep0Vec_LepESUp); leptonVecEnUp->push_back(Lep1Vec_LepESUp);
        leptonVecEnDown->push_back(Lep0Vec_LepESDown); leptonVecEnDown->push_back(Lep1Vec_LepESDown);
        MET_LepESUp = MET; MET_Phi_LepESUp = MET_Phi;
        MET_LepESDown = MET; MET_Phi_LepESDown = MET_Phi;
        METSystShift(leptonVec, leptonVecEnUp, MET_LepESUp, MET_Phi_LepESUp, MET, MET_Phi);
        METSystShift(leptonVec, leptonVecEnDown, MET_LepESDown, MET_Phi_LepESDown, MET, MET_Phi);
        MT2llSmearFactorBin = MT2llMeanSmear->FindBin(MT2ll);
        if (MT2llSmearFactorBin > 20) MT2llSmearFactorBin = 21;
        MT2llSmearFactor = MT2llMeanSmear->GetBinContent(MT2llSmearFactorBin);
//        cout << "MT2ll " << MT2ll << endl;
//        cout << "MT2ll Smear " << MT2llSmearFactor << endl;
        MT2ll_ShiftUp = MT2ll + rand.Gaus(0, MT2llSmearFactor);
        /******************************************************/
        //        MT2lb=getMT2(Lep0Vec, BtagJet0Vec, MET, MET_Phi);  
        //        cout << "NJets " << NJets << endl;
        if (NJets > 1) {
            if (NBtagJets > 1) {
                MT2lbPair1 = getMT2(Lep0Vec + BtagJet0Vec, Lep1Vec + BtagJet1Vec, MET, MET_Phi);
                MT2lbPair2 = getMT2(Lep0Vec + BtagJet1Vec, Lep1Vec + BtagJet0Vec, MET, MET_Phi);
                if (MT2lbPair1 > MT2lbPair2) {
                    BLep0Vec = Lep0Vec + BtagJet1Vec;
                    BLep1Vec = Lep1Vec + BtagJet0Vec;
                }
                else {
                    BLep0Vec = Lep0Vec + BtagJet0Vec;
                    BLep1Vec = Lep1Vec + BtagJet1Vec;                    
                }
                //                cout << "BLep0Vec Energy " << BLep0Vec.E() << endl;
                //                cout << "BLep1Vec Energy " << BLep1Vec.E() << endl;
            }
            else if (NBtagJets == 1) {
                if (eventJetParams->at(2) == 0) { //Btag jet is lead jet
                    MT2lbPair1 = getMT2(Lep0Vec + BtagJet0Vec, Lep1Vec + Jet1Vec, MET, MET_Phi);
                    MT2lbPair2 = getMT2(Lep0Vec + Jet1Vec, Lep1Vec + BtagJet0Vec, MET, MET_Phi);
                    if (MT2lbPair1 > MT2lbPair2) {
                        BLep0Vec = Lep0Vec + Jet1Vec;
                        BLep1Vec = Lep1Vec + BtagJet0Vec;
                    }
                    else {
                        BLep0Vec = Lep0Vec + BtagJet0Vec;
                        BLep1Vec = Lep1Vec + Jet1Vec;                        
                    }
                    //                    cout << "BLep0Vec Energy " << BLep0Vec.E() << endl;
                    //                    cout << "BLep1Vec Energy " << BLep1Vec.E() << endl;
                }
                else if (eventJetParams->at(2) != 0) { //Btag jet is not lead jet
                    MT2lbPair1 = getMT2(Lep0Vec + BtagJet0Vec, Lep1Vec + Jet0Vec, MET, MET_Phi);
                    MT2lbPair2 = getMT2(Lep0Vec + Jet0Vec, Lep1Vec + BtagJet0Vec, MET, MET_Phi);
                    if (MT2lbPair1 > MT2lbPair2) {
                        BLep0Vec = Lep0Vec + Jet0Vec;
                        BLep1Vec = Lep1Vec + BtagJet0Vec;
                    }
                    else {
                        BLep0Vec = Lep0Vec + BtagJet0Vec;
                        BLep1Vec = Lep1Vec + Jet0Vec;                        
                    }
                }
            }
            else {
                MT2lbPair1 = getMT2(Lep0Vec + Jet0Vec, Lep1Vec + Jet1Vec, MET, MET_Phi);
                MT2lbPair2 = getMT2(Lep0Vec + Jet1Vec, Lep1Vec + Jet0Vec, MET, MET_Phi);
                if (MT2lbPair1 > MT2lbPair2) {
                    BLep0Vec = Lep0Vec + Jet1Vec;
                    BLep1Vec = Lep1Vec + Jet0Vec;
                }
                else {
                    BLep0Vec = Lep0Vec + Jet0Vec;
                    BLep1Vec = Lep1Vec + Jet1Vec;                        
                }
            }
            /*
             cout << "BLep0Vec Energy " << BLep0Vec.E() << endl;
             cout << "BLep1Vec Energy " << BLep1Vec.E() << endl;
             cout << "BLep1Vec dpHI " << dPhi(BLep1Vec.Phi(), BLep0Vec.Phi()) << endl;
             */
            //            cout << "MT2lbPair1 " << MT2lbPair1 << endl;
            //            cout << "MT2lbPair2 " << MT2lbPair2 << endl;
            MT2lb = TMath::Min(MT2lbPair1, MT2lbPair2);
            
            //            cout << "MT2lb " << MT2lb << endl;
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
        stringKeyToVar["diLepPt"] = diLepPt;
        stringKeyToVar["diLepInvMass"] = diLepInvMass;
        stringKeyToVar["diLepEta"] = diLepEta;
        stringKeyToVar["diLepPhi"] = diLepPhi;
        stringKeyToVar["MT2ll"] = MT2ll;
        stringKeyToVar["MT2lb"] = MT2lb;
        stringKeyToVar["MET"] = MET;
        stringKeyToVar["METPhi"] = MET_Phi;
        stringKeyToVar["METPhi_noPhiCorr"] = MET_Phi_preCorr;
        stringKeyToVar["NJets"] = NJets;
        stringKeyToVar["NBJets"] = NBtagJets;
        stringKeyToVar["DPhiLep0Lep1"] = dPhi(Lep0Vec.Phi(), Lep1Vec.Phi());
        stringKeyToVar["DPhiLep0MET"] = dPhi((float) Lep0Vec.Phi(), MET_Phi);
        stringKeyToVar["DPhiLep1MET"] = dPhi((float) Lep0Vec.Phi(), MET_Phi);
        stringKeyToVar["DPhiLep0MET_PreCorr"] = dPhi((float) Lep0Vec.Phi(), MET_Phi_preCorr);
        stringKeyToVar["DPhiLep1MET_PreCorr"] = dPhi((float) Lep0Vec.Phi(), MET_Phi_preCorr);
        stringKeyToVar["nVtx"] = (float) nVtx;
        stringKeyToVar["nVtxTrue"] = (float) nVtxTrue;
        stringKeyToVar["METX"] = METX;
        stringKeyToVar["METY"] = METY;
        stringKeyToVar["METX_noPhiCorr"] = METX_preCorr;
        stringKeyToVar["METY_noPhiCorr"] = METY_preCorr;
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
                stringKeyToVar["DPhiLepB0LepB1"] = dPhi(BLep0Vec.Phi(), BLep1Vec.Phi());
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
            stringKeyToVar["ELepEJet_LepESShiftUp"] = Lep0Vec_LepESUp.E() + Lep1Vec_LepESUp.E() - Jet0Vec.E() - Jet1Vec.E();
            stringKeyToVar["ELepEJet_LepESShiftDown"] = Lep0Vec_LepESDown.E() + Lep1Vec_LepESDown.E() - Jet0Vec.E() - Jet1Vec.E();
            stringKeyToVar["MT2ll_MT2llShiftUp"] = MT2ll_ShiftUp;            
            stringKeyToVar["MET_LepESShiftUp"] = MET_LepESUp;
            stringKeyToVar["MET_LepESShiftDown"] = MET_LepESDown;
            stringKeyToVar["METPhi_LepESShiftUp"] = MET_Phi_LepESUp;
            stringKeyToVar["METPhi_LepESShiftDown"] = MET_Phi_LepESDown;
        }
        /*#######################
         MAKE SELECTION CUTS
         #####################*/
        for (unsigned int i = 0; i < subSampVec->size(); ++i) {
            S_Current = subSampVec->at(i);
            subSampBool[S_Current] = false;            
//            if (S_Current.histNameSuffix == "_inclusive") subSampBool[S_Current] = true;
            if (!(S_Current.whichdiLepType < 0 || Type == S_Current.whichdiLepType)) continue;
            if (!(NJets >= S_Current.cutNJets)) continue;
            if (!(NBtagJets >= S_Current.cutNBJets)) continue;
            if (!(S_Current.doZVeto < 0 || ZVeto == S_Current.doZVeto)) continue;
            if ((Type == 0 || Type == 1) && MET < S_Current.cutMET) continue;
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
            subSampBool[S_Current] = true;
            /*
             switch (Type) {
             case 0:
             if (ZVeto) {
                        if (S_Current.histNameSuffix == "_mumu_ZVeto") { 
                            subSampBool[S_Current] = true;   
                        }
                        else if (S_Current.histNameSuffix == "_mumu_ZVeto_Jet2BJet1" && NJets > 1 && NBtagJets > 0) {
                            subSampBool[S_Current] = true; 
                        }
                        else if (S_Current.histNameSuffix == "_mumu_ZVeto_Jet2" && NJets > 1) {
                            subSampBool[S_Current] = true;
                        }
                        else if ((S_Current.histNameSuffix == "_mumu_ZVeto_METGeq40_Jet2BJet1" || S_Current.histNameSuffix == "_FullCut") && NJets > 1 && NBtagJets > 0 && MET > 40) {
                            subSampBool[S_Current] = true; 
                        }
                    }
                    else {                        
                        if (S_Current.histNameSuffix == "_mumu") { 
                            subSampBool[S_Current] = true;   
                        }
                        else if (S_Current.histNameSuffix == "_mumu_Jet2BJet1" && NJets > 1 && NBtagJets > 0) {
                            subSampBool[S_Current] = true; 
                        }
                        else if (S_Current.histNameSuffix == "_mumu_Jet2" && NJets > 1) {
                            subSampBool[S_Current] = true;
                        }
                        else if (S_Current.histNameSuffix == "_mumu_METGeq40_Jet2BJet1" && NJets > 1 && NBtagJets > 0 && MET > 40) {
                            subSampBool[S_Current] = true; 
                        }
                    }
                    break;
                case 1:
                    if (ZVeto) {
                        if (S_Current.histNameSuffix == "_ee_ZVeto") { 
                            subSampBool[S_Current] = true;   
                        }
                        else if (S_Current.histNameSuffix == "_ee_ZVeto_Jet2BJet1" && NJets > 1 && NBtagJets > 0) {
                            subSampBool[S_Current] = true; 
                        }
                        else if (S_Current.histNameSuffix == "_ee_ZVeto_Jet2" && NJets > 1) {
                            subSampBool[S_Current] = true;
                        }
                        else if ((S_Current.histNameSuffix == "_ee_ZVeto_METGeq40_Jet2BJet1" || S_Current.histNameSuffix == "_FullCut") && NJets > 1 && NBtagJets > 0 && MET > 40) {
                            subSampBool[S_Current] = true; 
                        }
                    }
                    else {
                        if (S_Current.histNameSuffix == "_ee") { 
                            subSampBool[S_Current] = true;   
                        }
                        else if (S_Current.histNameSuffix == "_ee_Jet2BJet1" && NJets > 1 && NBtagJets > 0) {
                            subSampBool[S_Current] = true; 
                        }
                        else if (S_Current.histNameSuffix == "_ee_Jet2" && NJets > 1) {
                            subSampBool[S_Current] = true;
                        }
                        else if (S_Current.histNameSuffix == "_ee_METGeq40_Jet2BJet1" && NJets > 1 && NBtagJets > 0 && MET > 40) {
                            subSampBool[S_Current] = true; 
                        }
                    }
                    break;
                case 2:
                    if (S_Current.histNameSuffix == "_emu") { 
                        subSampBool[S_Current] = true;   
                    }
                    else if (S_Current.histNameSuffix == "_emu_Jet2" && NJets > 1) {
                        subSampBool[S_Current] = true;
                    }
                    else if ((S_Current.histNameSuffix == "_emu_Jet2BJet1" || S_Current.histNameSuffix == "_FullCut") && NJets > 1 && NBtagJets > 0) {
                        subSampBool[S_Current] = true;   
                    }
                    break;    
            }
            if (NJets > 1 && S_Current.histNameSuffix == "_NJetsGeq2") subSampBool[S_Current] = true;
            if (NBtagJets > 0 && S_Current.histNameSuffix == "_NBJetsGeq1") subSampBool[S_Current] = true;
            if (NBtagJets > 1 && S_Current.histNameSuffix == "_NBJetsGeq2") subSampBool[S_Current] = true;
            if (!ZVeto) {
                if (NJets > 1 && NBtagJets > 0) { 
                    if (S_Current.histNameSuffix == "_inZ_w_JetCuts") subSampBool[S_Current] = true;   
                    if ((Type == 0 || Type == 1) && MET > 40 && S_Current.histNameSuffix == "_inZ_w_JetCuts_w_METCut") subSampBool[S_Current] = true;
                }
                else {
                    if (S_Current.histNameSuffix == "_inZ_w_o_JetCuts") subSampBool[S_Current] = true;
                    if ((Type == 0 || Type == 1) && MET > 40 && S_Current.histNameSuffix == "_inZ_w_o_JetCuts_w_METCut") subSampBool[S_Current] = true;
                }
            }
            */
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
        /*#######################
         FILL PLOTS
         #####################*/          
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
                    //                    cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                    /*
                     if (H_Current.xVarKey == "leadBJetEn" && NBtagJets > 0) {
                     cout << "variable energy " << xIter->second << endl;
                     }
                     */
                    if (xIter != stringKeyToVar.end()) {
                        if (doVerbosity) {
                            cout << "xIter first " << xIter->first << endl;
                            cout << "xIter second " << xIter->second << endl;
                        }
                        ///Some necessary continue checks
                        if (H_Current.xVarKey == "MT2lb" && NJets < 2) continue;
                        if (S_Current.blindDataChannel && H_Current.xVarKey == "MT2ll") {
                            if (blindData && doData && MT2ll > MT2llCut) continue;
                        }
                        if (S_Current.blindDataChannel && H_Current.xVarKey == "MT2lb") {
                            if (blindData && doData && MT2lb > MT2lbCut) continue;
                        }
                        /*
                        if ((S_Current.histNameSuffix == "_FullCut" || S_Current.histNameSuffix == "_emu_Jet2BJet1" || S_Current.histNameSuffix == "_ee_ZVeto_METGeq40_Jet2BJet1" || S_Current.histNameSuffix == "_mumu_ZVeto_METGeq40_Jet2BJet1") && H_Current.xVarKey == "MT2ll") {
                            if (doData && MT2ll > MT2llCut) continue;                            
                        }
                        if ((S_Current.histNameSuffix == "_FullCut" || S_Current.histNameSuffix == "_emu_Jet2BJet1" || S_Current.histNameSuffix == "_ee_ZVeto_METGeq40_Jet2BJet1" || S_Current.histNameSuffix == "_mumu_ZVeto_METGeq40_Jet2BJet1") && H_Current.xVarKey == "MT2lb") {
                            if (doData && MT2lb > MT2lbCut) continue;                            
                        }
                        */
                        ///condition for which stuff to run on
//                        cout << "H_Current Name " << H_Current.name << endl;
                        if (H_Current.name.Contains("LepEffSFShiftUp")) {                            
                            fillWeight = ((H_Current.name.Contains("preRW")) ? preNVtxRWweight_LepEffSFUp : weight_LepEffSFUp);
//                            cout << "weight_LepEffSFUp " << weight_LepEffSFUp << endl;
//                            cout << "fillWeight " << fillWeight;
                        }
                        else if (H_Current.name.Contains("LepEffSFShiftDown")) {
                            fillWeight = ((H_Current.name.Contains("preRW")) ? preNVtxRWweight_LepEffSFDown : weight_LepEffSFDown);
//                            cout << "weight_LepEffSFDown " << weight_LepEffSFDown << endl;
//                            cout << "fillWeight " << fillWeight;
                        }
                        else {
                            fillWeight = ((H_Current.name.Contains("preRW")) ? preNVtxRWweight : weight);
                        }
                        histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                    }
                    if (doVerbosity) cout << "" << endl;
                    
                }            
                for (unsigned int k = 0; k < histVec_2D->size(); ++k) {
                    H_Current = histVec_2D->at(k);
                    xIter = stringKeyToVar.find(H_Current.xVarKey);
                    yIter = stringKeyToVar.find(H_Current.yVarKey);
                    if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end()) { 
                        if (H_Current.xVarKey == "MT2lb" && NJets < 2) continue;
                        /*
                         if (NJets > 1 && H_Current.xVarKey == "MT2lb") {
                         cout << "HCurrent xvar key " << H_Current.xVarKey << endl;
                         cout << "HCurrent yvar key " << H_Current.yVarKey << endl;
                         cout << "xiter second " << xIter->second << endl;
                         cout << "yiter second " << yIter->second << endl;
                         }
                         */
                        if (S_Current.blindDataChannel && H_Current.xVarKey == "MT2ll") {
                            if (blindData && doData && MT2ll > MT2llCut) continue;
                        }
                        if (S_Current.blindDataChannel && H_Current.xVarKey == "MT2lb") {
                            if (blindData && doData && MT2lb > MT2lbCut) continue;
                        }
                        /*
                        if ((S_Current.histNameSuffix == "_FullCut" || S_Current.histNameSuffix == "_emu_Jet2BJet1" || S_Current.histNameSuffix == "_ee_ZVeto_METGeq40_Jet2BJet1" || S_Current.histNameSuffix == "_mumu_ZVeto_METGeq40_Jet2BJet1") && H_Current.xVarKey == "MT2ll") {
                            if (doData && MT2ll > MT2llCut) continue;                            
                        }
                        if ((S_Current.histNameSuffix == "_FullCut" || S_Current.histNameSuffix == "_emu_Jet2BJet1" || S_Current.histNameSuffix == "_ee_ZVeto_METGeq40_Jet2BJet1" || S_Current.histNameSuffix == "_mumu_ZVeto_METGeq40_Jet2BJet1") && H_Current.xVarKey == "MT2lb") {
                            if (doData && MT2lb > MT2lbCut) continue;                            
                        }
                         */
                        histMap_2D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, weight); //there could be shenanigans with this one
                    }
                }
                /*
                 for (unsigned int m = 0; m < histVec_3D->size(); ++m) {
                 H_Current = histVec_3D->at(m);
                 xIter = stringKeyToVar.find(H_Current.xVarKey);
                 yIter = stringKeyToVar.find(H_Current.yVarKey); 
                 zIter = stringKeyToVar.find(H_Current.zVarKey);
                 if (xIter != stringKeyToVar.end() && yIter != stringKeyToVar.end() && zIter != stringKeyToVar.end()) {
                 histMap_3D[histKey(H_Current, S_Current)]->Fill(xIter->second, yIter->second, zIter->second, weight); //there could be shenanigans with this one   
                 }
                 }
                 */
                for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {
                    H_Current = histVec_1D_Syst->at(js);
                    xIter = stringKeyToVar.find(H_Current.xVarKey);
                    if (doVerbosity) {
                        cout << "" << endl;
                        cout << "ievt " << ievt << endl;
                        cout << "i " << i << endl;
                        cout << "js " << js << endl;
                        cout << "S_Current.histNamesuffix " << S_Current.histNameSuffix << endl;
                        cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                    }
                    //                    cout << "H_Current.xVarKey " << H_Current.xVarKey << endl;
                    /*
                     if (H_Current.xVarKey == "leadBJetEn" && NBtagJets > 0) {
                     cout << "variable energy " << xIter->second << endl;
                     }
                     */
                    if (xIter != stringKeyToVar.end()) {
                        if (doVerbosity) {
                            cout << "xIter first " << xIter->first << endl;
                            cout << "xIter second " << xIter->second << endl;
                        }
                        ///Some necessary continue checks
                        if (H_Current.xVarKey == "MT2lb" && NJets < 2) continue;
                        if (S_Current.blindDataChannel && H_Current.xVarKey == "MT2ll") {
                            if (blindData && doData && MT2ll > MT2llCut) continue;
                        }
                        if (S_Current.blindDataChannel && H_Current.xVarKey == "MT2lb") {
                            if (blindData && doData && MT2lb > MT2lbCut) continue;
                        }
                        /*
                        if ((S_Current.histNameSuffix == "_FullCut" || S_Current.histNameSuffix == "_emu_Jet2_BJet1" || S_Current.histNameSuffix == "_ee_ZVeto_Jet2_BJet1_METGeq40" || S_Current.histNameSuffix == "_mumu_ZVeto_Jet2_BJet1_METGeq40") && H_Current.xVarKey == "MT2ll") {
                            if (doData && MT2ll > MT2llCut) continue;                            
                        }
                        if ((S_Current.histNameSuffix == "_FullCut" || S_Current.histNameSuffix == "_emu_Jet2_BJet1" || S_Current.histNameSuffix == "_ee_ZVeto_Jet2_BJet1_METGeq40" || S_Current.histNameSuffix == "_mumu_ZVeto_Jet2_BJet1_METGeq40") && H_Current.xVarKey == "MT2lb") {
                            if (doData && MT2lb > MT2lbCut) continue;                            
                        }
                        */
                        ///condition for which stuff to run on
                        if (H_Current.name.Contains("LepEffSFShiftUp")) {
                            fillWeight = ((H_Current.name.Contains("preRW")) ? preNVtxRWweight_LepEffSFUp : weight_LepEffSFUp);
                        }
                        else if (H_Current.name.Contains("LepEffSFShiftDown")) {
                            fillWeight = ((H_Current.name.Contains("preRW")) ? preNVtxRWweight_LepEffSFDown : weight_LepEffSFDown);
                        }
                        else {
                            fillWeight = ((H_Current.name.Contains("preRW")) ? preNVtxRWweight : weight);
                        }
                        histMap_1D[histKey(H_Current, S_Current)]->Fill(xIter->second, fillWeight);
                    }
                    if (doVerbosity) cout << "" << endl;                    
                }
            }
        }        
    }
    cout << "All events done" << endl;
    outputFile->cd();
    cout << "cd-ing to output directory" << endl;
    
    if (whichNTupleType == 1 && !doData) {
        eventCount = (TH1F*) inputFile.Get("EventsBeforeSelection/weightedEvents");
        eventCount->Write();
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
float getMT2(TLorentzVector lept1, TLorentzVector lept2, float theMET, float theMETphi){
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
