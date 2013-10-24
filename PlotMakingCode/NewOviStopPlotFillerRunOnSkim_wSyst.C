
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
#include <map>

#include "TCut.h"
//#include "StopDict_ver2.h"
#include "../HeaderFiles/StopStructDefinitions.h"
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

typedef std::pair<HistogramT, SampleT> histKey;
typedef std::map<histKey, TH1 *>      HMap_1D;
typedef std::map<histKey, TH2 *>      HMap_2D;
typedef std::map<histKey, TH3 *>      HMap_3D;
typedef std::map<SampleT, bool>       passCutMap;

StopXSec getCrossSectionStop(float);

using namespace std;
//different nVtx regions for the plots
int main( int argc, const char* argv[] ) {
    /////////Variable initializations/////////////////
    /////Organization Variables//////////////   
    TFile * inputWeightFile = new TFile("StopABCPUWeight_vers2.root");
    TH3F * weights = (TH3F*) inputWeightFile->Get("WHist");
    TH1F * OneDPUDist = (TH1F*) weights->ProjectionX("OneDWeight");
    TFile * inputPUFileMCOvi_v2 = new TFile(TString("OneDMCPURWNewOvi.root"));
    TFile * inputPUFileMCOvi = new TFile(TString("OneDMCPURWNewOvi_v2.root"));
    TFile * inputPUFileMCDESY = new TFile(TString("OneDMCPURWNewDESY.root"));
    TFile * inputPUFileMCOviToDESY = new TFile(TString("OneDMCPURW_OviToDESY.root"));
    TFile * inputPUS7toS10 = new TFile(TString("S7toS10RW.root"));

    //    TH1F * truePUDistMC = (TH1F*) inputPUFileMC->Get("MCPU");
    //    TFile * inputPUFileData = new TFile("RunABCPUDist.root");
    //    TH1F * truePUDistData = (TH1F*) inputPUFileData->Get("pileup");
    //    float intMCPU = truePUDistMC->
    TH1F * nVtxSFHist = (TH1F*) inputPUFileMCOvi->Get("DataDivZDY");
    TH1F * nVtxSFHist_v2 = (TH1F*) inputPUFileMCOvi_v2->Get("nVtxSF_preRW");
    TH1F * nVtxSFHistOviToDESY = (TH1F *) inputPUFileMCOviToDESY->Get("normFrac");
    TH1F * h_S7toS10RWHist = (TH1F *) inputPUS7toS10->Get("S7toS10RWHist");
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
    TH1F * h_eventCount, * h_CutFlow;    
    /////Event Variables/////////////////////
    bool  keepEvent;

    BasicEventInfo       BEI;
    EventLepInfo         ELI; ELI.ELIDefaultVarVals();
    EventJetInfo         EJI; EJI.EJIDefaultVarVals();
    EventMETInfo         EMI; EMI.EMIDefaultVarVals();
    EventDiStructureInfo EDSI; EDSI.EDSIDefaultVarVals();
    EventGenParticleInfo EGPI; EGPI.EGPIDefaultVarVals();
    
    BasicEventInfo       BEI_LepESUp, BEI_LepESDown;
    BasicEventInfo       BEI_GenTopRW;
    BasicEventInfo       BEI_StopXSecUp, BEI_StopXSecDown;
    BasicEventInfo       BEI_LepEffSFUp, BEI_LepEffSFDown;
    
    EventLepInfo         ELI_LepESUp; ELI_LepESUp.ELIDefaultVarVals();
    EventLepInfo         ELI_LepESDown; ELI_LepESDown.ELIDefaultVarVals();
    
    EventJetInfo         EJI_JetESUp; EJI_JetESUp.EJIDefaultVarVals();
    EventJetInfo         EJI_JetESDown; EJI_JetESDown.EJIDefaultVarVals();
    EventJetInfo         EJI_BTagSFUp; EJI_BTagSFUp.EJIDefaultVarVals();
    EventJetInfo         EJI_BTagSFDown; EJI_BTagSFDown.EJIDefaultVarVals();
    
    EventMETInfo         EMI_LepESUp; EMI_LepESUp.EMIDefaultVarVals();
    EventMETInfo         EMI_LepESDown; EMI_LepESDown.EMIDefaultVarVals();
    EventMETInfo         EMI_JetESUp; EMI_JetESUp.EMIDefaultVarVals();
    EventMETInfo         EMI_JetESDown; EMI_JetESDown.EMIDefaultVarVals();
    EventMETInfo         EMI_BTagSFUp; EMI_BTagSFUp.EMIDefaultVarVals();
    EventMETInfo         EMI_BTagSFDown; EMI_BTagSFDown.EMIDefaultVarVals();
    
    EventDiStructureInfo EDSI_LepESUp; EDSI_LepESUp.EDSIDefaultVarVals();
    EventDiStructureInfo EDSI_LepESDown; EDSI_LepESDown.EDSIDefaultVarVals();
    EventDiStructureInfo EDSI_JetESUp; EDSI_JetESUp.EDSIDefaultVarVals();
    EventDiStructureInfo EDSI_JetESDown; EDSI_JetESDown.EDSIDefaultVarVals();
    EventDiStructureInfo EDSI_BTagSFUp; EDSI_BTagSFUp.EDSIDefaultVarVals();
    EventDiStructureInfo EDSI_BTagSFDown; EDSI_BTagSFDown.EDSIDefaultVarVals();
    
    EventSpecialMT2Info  ESMT2I;
    
    int TGenStopMass0,TGenStopMass1,TGenChi0Mass0,TGenChi0Mass1, TGenCharginoMass0, TGenCharginoMass1;
    int   grabStopMass, grabChi0Mass, grabCharginoMass;
    int   massDiffThresh = 5;
    /*
    float stopWeight    = 0.;
    float stopWeightErr = 0.;
    float stopWeightPlusErr = 0.;
    float stopWeightMinusErr = 0.;
    */
    
    
    const double PI = 3.14159265;
    
    TLorentzVector genMETVec, genTop0Vec, genTop1Vec;
    TLorentzVector genStop0Vec, genStop1Vec;
    TLorentzVector genChiZero0Vec, genChiZero1Vec;
    TLorentzVector genChargino0Vec, genChargino1Vec;
    
    float genTop0Pt, genTop0En, genTop0Eta, genTop0Phi;
    float genTop1Pt, genTop1En, genTop1Eta, genTop1Phi;
    int genTop0PdgId, genTop1PdgId; //genTop0FirstMom, genTop1FirstMom;
    float genMET_Pt, genMET_Phi;
    
    TH1F * h_genTopPt = new TH1F("h_genTopPt", ";gen-level Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genTopPt->Sumw2();
    TH1F * h_genTopPtRW = new TH1F("h_genTopPtRW", ";gen-level Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genTopPtRW->Sumw2();
    TH1F * h_genAntiTopPt = new TH1F("h_genAntiTopPt", ";gen-level Anti-Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genAntiTopPt->Sumw2();
    TH1F * h_genAntiTopPtRW = new TH1F("h_genAntiTopPtRW", ";gen-level Anti-Top p_{T} [GeV];Number of Events / NUM GeV", 500, 0, 1000); h_genAntiTopPtRW->Sumw2();
    
    /********************************************************************************************************/
    TString eventCountName;
    TFile * MT2llSmearFile = new TFile("MT2llSmear.root");
    TH1F * MT2llMeanSmear = (TH1F *) MT2llSmearFile->Get("MT2llSmear");        
    /********************************************************************************************************/
    
    
//    float lumi = 19300; //5296.3; // ipb                                                                                 
    // add 5 GeV safety margin (deltaM = 10 GeV in the FineBin sample)  
//    float Nevt_stop_oneMassPoint = 50000 * ( (genStopMassMax-genStopMassMin)/10. ) * ( (genDeltaM_stopChi0_Max-genDeltaM_stopChi0_Min)/10. );  
    // 50k evts per point x Npoints
    
    ////input cuts/commands    
    
    bool grabOutDir      = 0;      // whether or not to use the file: "outputSavePath.txt" for where to save output
    BEI.doPhiCorr       = 1;      // whether to do the MetPhi asymmetry correction -- 6/25/13 as of right now parameters need to be updated
    BEI.doData          = 0;      // Whether you're running on data or not
    bool doVerbosity     = 0;      // prints a lot of debug info if turned on
    int  whichNTupleType = 0;      // 0 IFCA Oviedo; 1 DESY -- 6/25/13 as of right now, doesn't work with Oviedo
    BEI.doPURW          = 0;      // run pile up reweighting
    BEI.doHackPURW      = 0;      // called a "hack" because it was bin by bin reweighting to get the nVtx distribution in the "inclusive" channel to exactly match between data and MC
    BEI.doPURWOviToDESY = 0;      // exactly like BEI.doHackPURW but for Oviedo and DESY to enable comparisons across nTuples
    bool doBookSyst      = 0;      // used for deciding whether or not to book systematics
    bool doMETSmear      = false;  // early attempt to mimic jet smearing using gaussians -- didn't work, leave off
    double METSF         = 0.0;
    
    BEI.blindData       = 1;      // "Until further notice, leave ON -- cuts MT2ll > 80 is data in full object/event selection
    int  nEvents         = -1;     // limits total number of events one is running on
    int subLepPtCut      = 10;    // Sets the pT cut used for subLepPtCut
    
    bool doParallel      = 0;      // Whether or not to run parallel jobs for a given input file. If done, you have to do an additional hadding step
    bool isLimStats      = 0;
    int  startPointNum   = 1;
    int  numBreakPoints  = 0;
    BEI.isSignal         = 0;
    //bool isSignal  = 0;      // Whether or not one is running on signal -- if this is turned on the next two arguments have to be the stop mass to grab and the chi0 mass to grab respectively --- NOTE Now contained in the BasicEventInfo structure
    BEI.doReReco        = 0;
    int  startPoint, endPoint;
    TString outputSavePathString = "outputSavePath";
    ////input cuts/commands    
    /////loop over inputs    
    for (int k = 0; k < argc; ++k) {
        cout << "argv[k] for k = " << k << " is: " << argv[k] << endl;
        if (strncmp (argv[k],"-i",2) == 0) {
            fInName = TString(argv[k+1]);   
        }
        else if (strncmp (argv[k],"noPhiCorr",9) == 0) {
            BEI.doPhiCorr = 0;   
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
            BEI.doPURW = 1;           
        }
        else if (strncmp (argv[k],"doHackPURW",10) == 0) {
            BEI.doHackPURW = 1;
        }
        else if (strncmp (argv[k],"doPURWOviToDESY",15) == 0) {
            BEI.doPURWOviToDESY = 1;   
        }
        else if (strncmp (argv[k],"doBookSyst", 10) == 0) {
            doBookSyst = 1;
        }
        else if (strncmp (argv[k],"gOutDir", 7) == 0) {
            grabOutDir = 1;
        }        
        else if (strncmp (argv[k], "OutDirName", 10) == 0) {
            outputSavePathString = TString(argv[k+1]);
        }
        else if (strncmp (argv[k],"isSig", 5) == 0) {
            BEI.isSignal = 1;
            grabStopMass = strtol(argv[k+1], NULL, 10);
            grabChi0Mass = strtol(argv[k+2], NULL, 10);
            grabCharginoMass = strtol(argv[k+3], NULL, 10);
            cout << "Looking for StopMass: " << grabStopMass << ", Chi0Mass: " << grabChi0Mass << ", and CharginoMass: " << grabCharginoMass << endl;
        }
        else if (strncmp (argv[k],"ReleaseTheKraken", 16) == 0) {
            BEI.blindData = 0;
            cout << "RELEASING THE KRAKEN!!! " << endl;
            cout << "http://www.youtube.com/watch?v=gb2zIR2rvRQ " << endl;
        }
        else if (strncmp (argv[k],"limStats",8) == 0) {
            isLimStats = 1;
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
    TString outputPathName = (BEI.isSignal) ? "signalOutputSavePath" : outputSavePathString;
    if (grabOutDir) {
        outDirFile = new ifstream(outputPathName + TString(".txt"));
        if (!(outDirFile->eof())) {
            outDirFile->getline(Buffer,500);
            fOutName += TString(string(Buffer));
            fOutName += "/"; //in case user forgot a slash
        }
    }
    fOutName += fInName(fCutSlash);
    if (fInName.Contains("MuEG") || fInName.Contains("DoubleMu") || fInName.Contains("DoubleEl") || fInName.Contains("run2012")) {
        cout << "Running on Data" << endl;
        BEI.doData = 1;
    }
    if (fInName.Contains("ReReco")) BEI.doReReco = 1;
    if (BEI.doData && BEI.blindData) {
        fOutName += "MT2Leq80";   
    }
    else if (BEI.doData && !BEI.blindData) {
        fOutName += "_NOTBLIND";
    }
    if (BEI.doData) doMETSmear = 0;
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
    if (BEI.doData && doBookSyst) doBookSyst = 0;
    if (BEI.doPURW && !BEI.doData) fOutName += "_PURW";
    if (BEI.doPURWOviToDESY && !BEI.doData) fOutName += "OviToDESY";
    if (doBookSyst) fOutName += "_wSyst";
    if (BEI.isSignal) {
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
    if (isLimStats) {
        fOutName += "_isLimStats_";
        fOutName += nEvents;
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
    float ZWindowLB      = 76;
    float ZWindowUB      = 106;   
    /////Set up the tree////////
    fileTree.Add(fInName + TString(".root"));
    
    if (whichNTupleType == 0) {
        fileTree.SetBranchAddress( "TPassDoubleMu",     &BEI.passTrigDoubleMu );
        fileTree.SetBranchAddress( "TPassDoubleEl",     &BEI.passTrigDoubleEl );
        fileTree.SetBranchAddress( "TPassElMu",         &BEI.passTrigElMu );
    }    
    fileTree.SetBranchAddress( "TWeight",               &BEI.weight );
    fileTree.SetBranchAddress( "TNPV",                  &BEI.nVtx );
    fileTree.SetBranchAddress( "TNPV_True",             &BEI.nVtxTrue );
    fileTree.SetBranchAddress( "TRunNum",               &BEI.runNumber );
    fileTree.SetBranchAddress( "TEventNum",             &BEI.eventNumber );
    fileTree.SetBranchAddress( "TLumiBlock",            &BEI.lumiBlock );
    
    fileTree.SetBranchAddress( "TMET",     &EMI.EventMET_preCorr);
    fileTree.SetBranchAddress( "TMET_Phi", &EMI.EventMETPhi_preCorr);
    fileTree.SetBranchAddress( "TMETSig",  &EMI.EventMETSig_preCorr);
    
    fileTree.SetBranchAddress( "TDoEvent",      &ELI.doEvent);    
    fileTree.SetBranchAddress( "TChannel",      &ELI.EventDiLepType);
    fileTree.SetBranchAddress( "TDiLepMass",    &ELI.EventDiLepMass);
    
    fileTree.SetBranchAddress( "TLep0Px",       &ELI.EventLep0Px);
    fileTree.SetBranchAddress( "TLep0Py",       &ELI.EventLep0Py);
    fileTree.SetBranchAddress( "TLep0Pz",       &ELI.EventLep0Pz);
    fileTree.SetBranchAddress( "TLep0E",        &ELI.EventLep0E);
    fileTree.SetBranchAddress( "TLep0PdgId",    &ELI.EventLep0PDGID);
    fileTree.SetBranchAddress( "TLep0RelPFIso", &ELI.EventLep0RelPFIso);
    
    fileTree.SetBranchAddress( "TLep1Px",       &ELI.EventLep1Px);
    fileTree.SetBranchAddress( "TLep1Py",       &ELI.EventLep1Py);
    fileTree.SetBranchAddress( "TLep1Pz",       &ELI.EventLep1Pz);
    fileTree.SetBranchAddress( "TLep1E",        &ELI.EventLep1E);
    fileTree.SetBranchAddress( "TLep1PdgId",    &ELI.EventLep1PDGID);
    fileTree.SetBranchAddress( "TLep1RelPFIso", &ELI.EventLep1RelPFIso);
    
    
    fileTree.SetBranchAddress("TNJets",           &EJI.EventNJets); 
    fileTree.SetBranchAddress("TNJetsBtag",       &EJI.EventNBtagJets);     
    fileTree.SetBranchAddress("THT",              &EJI.EventHT);
    
    fileTree.SetBranchAddress("TJet0Px",          &EJI.EventJet0Px); 
    fileTree.SetBranchAddress("TJet0Py",          &EJI.EventJet0Py); 
    fileTree.SetBranchAddress("TJet0Pz",          &EJI.EventJet0Pz); 
    fileTree.SetBranchAddress("TJet0E",           &EJI.EventJet0E);
    
    fileTree.SetBranchAddress("TJet1Px",          &EJI.EventJet1Px); 
    fileTree.SetBranchAddress("TJet1Py",          &EJI.EventJet1Py); 
    fileTree.SetBranchAddress("TJet1Pz",          &EJI.EventJet1Pz); 
    fileTree.SetBranchAddress("TJet1E",           &EJI.EventJet1E);  
    
    fileTree.SetBranchAddress("TBtagJet0Px",      &EJI.EventBtagJet0Px); 
    fileTree.SetBranchAddress("TBtagJet0Py",      &EJI.EventBtagJet0Py); 
    fileTree.SetBranchAddress("TBtagJet0Pz",      &EJI.EventBtagJet0Pz); 
    fileTree.SetBranchAddress("TBtagJet0E",       &EJI.EventBtagJet0E);
    fileTree.SetBranchAddress("TBtagJet0Index",   &EJI.EventBtagJet0Index);
    
    fileTree.SetBranchAddress("TBtagJet1Px",      &EJI.EventBtagJet1Px); 
    fileTree.SetBranchAddress("TBtagJet1Py",      &EJI.EventBtagJet1Py); 
    fileTree.SetBranchAddress("TBtagJet1Pz",      &EJI.EventBtagJet1Pz); 
    fileTree.SetBranchAddress("TBtagJet1E",       &EJI.EventBtagJet1E);
    fileTree.SetBranchAddress("TBtagJet1Index",   &EJI.EventBtagJet1Index);
    
    if (!BEI.doData) {
        
        fileTree.SetBranchAddress("TDoEvent_LepESUp",       &ELI_LepESUp.doEvent);
        fileTree.SetBranchAddress("TChannel_LepESUp",       &ELI_LepESUp.EventDiLepType);
        fileTree.SetBranchAddress("TDiLepMass_LepESUp",     &ELI_LepESUp.EventDiLepMass);
        
        fileTree.SetBranchAddress("TLep0Px_LepESUp",        &ELI_LepESUp.EventLep0Px); 
        fileTree.SetBranchAddress("TLep0Py_LepESUp",        &ELI_LepESUp.EventLep0Py); 
        fileTree.SetBranchAddress("TLep0Pz_LepESUp",        &ELI_LepESUp.EventLep0Pz); 
        fileTree.SetBranchAddress("TLep0E_LepESUp",         &ELI_LepESUp.EventLep0E);
        fileTree.SetBranchAddress("TLep0PdgId_LepESUp",       &ELI_LepESUp.EventLep0PDGID);
        fileTree.SetBranchAddress("TLep0RelPFIso_LepESUp",    &ELI_LepESUp.EventLep0RelPFIso);
        
        fileTree.SetBranchAddress("TLep1Px_LepESUp",          &ELI_LepESUp.EventLep1Px); 
        fileTree.SetBranchAddress("TLep1Py_LepESUp",          &ELI_LepESUp.EventLep1Py); 
        fileTree.SetBranchAddress("TLep1Pz_LepESUp",          &ELI_LepESUp.EventLep1Pz); 
        fileTree.SetBranchAddress("TLep1E_LepESUp",           &ELI_LepESUp.EventLep1E);
        fileTree.SetBranchAddress("TLep1PdgId_LepESUp",       &ELI_LepESUp.EventLep1PDGID);
        fileTree.SetBranchAddress("TLep1RelPFIso_LepESUp",    &ELI_LepESUp.EventLep1RelPFIso);        
        
        fileTree.SetBranchAddress("TDoEvent_LepESDown",       &ELI_LepESDown.doEvent);
        fileTree.SetBranchAddress("TChannel_LepESDown",       &ELI_LepESDown.EventDiLepType);        
        fileTree.SetBranchAddress("TDiLepMass_LepESDown",     &ELI_LepESDown.EventDiLepMass);
        
        fileTree.SetBranchAddress("TLep0Px_LepESDown",        &ELI_LepESDown.EventLep0Px); 
        fileTree.SetBranchAddress("TLep0Py_LepESDown",        &ELI_LepESDown.EventLep0Py); 
        fileTree.SetBranchAddress("TLep0Pz_LepESDown",        &ELI_LepESDown.EventLep0Pz); 
        fileTree.SetBranchAddress("TLep0E_LepESDown",         &ELI_LepESDown.EventLep0E);
        fileTree.SetBranchAddress("TLep0PdgId_LepESDown",     &ELI_LepESDown.EventLep0PDGID);
        fileTree.SetBranchAddress("TLep0RelPFIso_LepESDown",  &ELI_LepESDown.EventLep0RelPFIso);
        
        fileTree.SetBranchAddress("TLep1Px_LepESDown",        &ELI_LepESDown.EventLep1Px); 
        fileTree.SetBranchAddress("TLep1Py_LepESDown",        &ELI_LepESDown.EventLep1Py); 
        fileTree.SetBranchAddress("TLep1Pz_LepESDown",        &ELI_LepESDown.EventLep1Pz); 
        fileTree.SetBranchAddress("TLep1E_LepESDown",         &ELI_LepESDown.EventLep1E);
        fileTree.SetBranchAddress("TLep1PdgId_LepESDown",     &ELI_LepESDown.EventLep1PDGID);
        fileTree.SetBranchAddress("TLep1RelPFIso_LepESDown",  &ELI_LepESDown.EventLep1RelPFIso);
        
        fileTree.SetBranchAddress("TNJets_JetESUp",           &EJI_JetESUp.EventNJets); 
        fileTree.SetBranchAddress("TNJetsBtag_JetESUp",       &EJI_JetESUp.EventNBtagJets);     
        fileTree.SetBranchAddress("THT_JetESUp",              &EJI_JetESUp.EventHT);
        
        fileTree.SetBranchAddress("TJet0Px_JetESUp",          &EJI_JetESUp.EventJet0Px); 
        fileTree.SetBranchAddress("TJet0Py_JetESUp",          &EJI_JetESUp.EventJet0Py); 
        fileTree.SetBranchAddress("TJet0Pz_JetESUp",          &EJI_JetESUp.EventJet0Pz); 
        fileTree.SetBranchAddress("TJet0E_JetESUp",           &EJI_JetESUp.EventJet0E);
        fileTree.SetBranchAddress("TJet1Px_JetESUp",          &EJI_JetESUp.EventJet1Px); 
        fileTree.SetBranchAddress("TJet1Py_JetESUp",          &EJI_JetESUp.EventJet1Py); 
        fileTree.SetBranchAddress("TJet1Pz_JetESUp",          &EJI_JetESUp.EventJet1Pz); 
        fileTree.SetBranchAddress("TJet1E_JetESUp",           &EJI_JetESUp.EventJet1E);
        
        fileTree.SetBranchAddress("TBtagJet0Px_JetESUp",      &EJI_JetESUp.EventBtagJet0Px); 
        fileTree.SetBranchAddress("TBtagJet0Py_JetESUp",      &EJI_JetESUp.EventBtagJet0Py); 
        fileTree.SetBranchAddress("TBtagJet0Pz_JetESUp",      &EJI_JetESUp.EventBtagJet0Pz); 
        fileTree.SetBranchAddress("TBtagJet0E_JetESUp",       &EJI_JetESUp.EventBtagJet0E);
        fileTree.SetBranchAddress("TBtagJet0Index_JetESUp",   &EJI_JetESUp.EventBtagJet0Index);    
        fileTree.SetBranchAddress("TBtagJet1Px_JetESUp",      &EJI_JetESUp.EventBtagJet1Px); 
        fileTree.SetBranchAddress("TBtagJet1Py_JetESUp",      &EJI_JetESUp.EventBtagJet1Py); 
        fileTree.SetBranchAddress("TBtagJet1Pz_JetESUp",      &EJI_JetESUp.EventBtagJet1Pz); 
        fileTree.SetBranchAddress("TBtagJet1E_JetESUp",       &EJI_JetESUp.EventBtagJet1E);
        fileTree.SetBranchAddress("TBtagJet1Index_JetESUp",   &EJI_JetESUp.EventBtagJet1Index);        
        
        fileTree.SetBranchAddress("TNJets_JetESDown",           &EJI_JetESDown.EventNJets); 
        fileTree.SetBranchAddress("TNJetsBtag_JetESDown",       &EJI_JetESDown.EventNBtagJets);     
        fileTree.SetBranchAddress("THT_JetESDown",              &EJI_JetESDown.EventHT);
        
        fileTree.SetBranchAddress("TJet0Px_JetESDown",          &EJI_JetESDown.EventJet0Px); 
        fileTree.SetBranchAddress("TJet0Py_JetESDown",          &EJI_JetESDown.EventJet0Py); 
        fileTree.SetBranchAddress("TJet0Pz_JetESDown",          &EJI_JetESDown.EventJet0Pz); 
        fileTree.SetBranchAddress("TJet0E_JetESDown",           &EJI_JetESDown.EventJet0E);
        fileTree.SetBranchAddress("TJet1Px_JetESDown",          &EJI_JetESDown.EventJet1Px); 
        fileTree.SetBranchAddress("TJet1Py_JetESDown",          &EJI_JetESDown.EventJet1Py); 
        fileTree.SetBranchAddress("TJet1Pz_JetESDown",          &EJI_JetESDown.EventJet1Pz); 
        fileTree.SetBranchAddress("TJet1E_JetESDown",           &EJI_JetESDown.EventJet1E);
        
        fileTree.SetBranchAddress("TBtagJet0Px_JetESDown",      &EJI_JetESDown.EventBtagJet0Px); 
        fileTree.SetBranchAddress("TBtagJet0Py_JetESDown",      &EJI_JetESDown.EventBtagJet0Py); 
        fileTree.SetBranchAddress("TBtagJet0Pz_JetESDown",      &EJI_JetESDown.EventBtagJet0Pz); 
        fileTree.SetBranchAddress("TBtagJet0E_JetESDown",       &EJI_JetESDown.EventBtagJet0E);
        fileTree.SetBranchAddress("TBtagJet0Index_JetESDown",   &EJI_JetESDown.EventBtagJet0Index);    
        fileTree.SetBranchAddress("TBtagJet1Px_JetESDown",      &EJI_JetESDown.EventBtagJet1Px); 
        fileTree.SetBranchAddress("TBtagJet1Py_JetESDown",      &EJI_JetESDown.EventBtagJet1Py); 
        fileTree.SetBranchAddress("TBtagJet1Pz_JetESDown",      &EJI_JetESDown.EventBtagJet1Pz); 
        fileTree.SetBranchAddress("TBtagJet1E_JetESDown",       &EJI_JetESDown.EventBtagJet1E);
        fileTree.SetBranchAddress("TBtagJet1Index_JetESDown",   &EJI_JetESDown.EventBtagJet1Index);
        
        
        fileTree.SetBranchAddress("TNJetsBtag_BTagSFUp",       &EJI_BTagSFUp.EventNBtagJets);
        
        fileTree.SetBranchAddress("TBtagJet0Px_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet0Px); 
        fileTree.SetBranchAddress("TBtagJet0Py_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet0Py); 
        fileTree.SetBranchAddress("TBtagJet0Pz_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet0Pz); 
        fileTree.SetBranchAddress("TBtagJet0E_BTagSFUp",       &EJI_BTagSFUp.EventBtagJet0E);
        fileTree.SetBranchAddress("TBtagJet0Index_BTagSFUp",   &EJI_BTagSFUp.EventBtagJet0Index);    
        fileTree.SetBranchAddress("TBtagJet1Px_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet1Px); 
        fileTree.SetBranchAddress("TBtagJet1Py_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet1Py); 
        fileTree.SetBranchAddress("TBtagJet1Pz_BTagSFUp",      &EJI_BTagSFUp.EventBtagJet1Pz); 
        fileTree.SetBranchAddress("TBtagJet1E_BTagSFUp",       &EJI_BTagSFUp.EventBtagJet1E);
        fileTree.SetBranchAddress("TBtagJet1Index_BTagSFUp",   &EJI_BTagSFUp.EventBtagJet1Index);
        
        fileTree.SetBranchAddress("TNJetsBtag_BTagSFDown",       &EJI_BTagSFDown.EventNBtagJets);
        
        fileTree.SetBranchAddress("TBtagJet0Px_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet0Px); 
        fileTree.SetBranchAddress("TBtagJet0Py_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet0Py); 
        fileTree.SetBranchAddress("TBtagJet0Pz_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet0Pz); 
        fileTree.SetBranchAddress("TBtagJet0E_BTagSFDown",       &EJI_BTagSFDown.EventBtagJet0E);
        fileTree.SetBranchAddress("TBtagJet0Index_BTagSFDown",   &EJI_BTagSFDown.EventBtagJet0Index);    
        fileTree.SetBranchAddress("TBtagJet1Px_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet1Px); 
        fileTree.SetBranchAddress("TBtagJet1Py_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet1Py); 
        fileTree.SetBranchAddress("TBtagJet1Pz_BTagSFDown",      &EJI_BTagSFDown.EventBtagJet1Pz); 
        fileTree.SetBranchAddress("TBtagJet1E_BTagSFDown",       &EJI_BTagSFDown.EventBtagJet1E);
        fileTree.SetBranchAddress("TBtagJet1Index_BTagSFDown",   &EJI_BTagSFDown.EventBtagJet1Index);        
        
        
        fileTree.SetBranchAddress("TMET_LepESUp",             &EMI_LepESUp.EventMET_preCorr);
        fileTree.SetBranchAddress("TMET_LepESDown",           &EMI_LepESDown.EventMET_preCorr);
        fileTree.SetBranchAddress("TMET_JetESUp",             &EMI_JetESUp.EventMET_preCorr);
        fileTree.SetBranchAddress("TMET_JetESDown",           &EMI_JetESDown.EventMET_preCorr);
        fileTree.SetBranchAddress("TMET_Phi_LepESUp",         &EMI_LepESUp.EventMETPhi_preCorr);
        fileTree.SetBranchAddress("TMET_Phi_LepESDown",       &EMI_LepESDown.EventMETPhi_preCorr);
        fileTree.SetBranchAddress("TMET_Phi_LepESUp",         &EMI_JetESUp.EventMETPhi_preCorr);
        fileTree.SetBranchAddress("TMET_Phi_LepESDown",       &EMI_JetESDown.EventMETPhi_preCorr);
    }        
    if (fileTree.GetBranch("TGenMET")) {
        BEI.hasMETInfo = true;
        fileTree.SetBranchAddress("TGenMET", &genMET_Pt );
        fileTree.SetBranchAddress("TGenMETPhi", &genMET_Phi );
    }
    if (fileTree.GetBranch("TGenTopSt3_0_Pt")) {
        BEI.hasTopInfo = true;
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
    if (fileTree.GetBranch("TGenStopMass0")){
        BEI.hasStopInfo = true;
        fileTree.SetBranchAddress( "TGenStopMass0",     &TGenStopMass0 );
        fileTree.SetBranchAddress( "TGenStopMass1",     &TGenStopMass1 );
        fileTree.SetBranchAddress( "TGenChi0Mass0",     &TGenChi0Mass0 );
        fileTree.SetBranchAddress( "TGenChi0Mass1",     &TGenChi0Mass1 );            
        fileTree.SetBranchAddress( "TGenCharginoMass0", &TGenCharginoMass0 );
        fileTree.SetBranchAddress( "TGenCharginoMass1", &TGenCharginoMass1 );
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
        EGPI.StopWeight = theStopXSec.stopProdXsec;
        EGPI.StopWeightErr = theStopXSec.stopProdXsecUncert * EGPI.StopWeight;
        EGPI.StopWeightPlusErr = EGPI.StopWeight + EGPI.StopWeightErr;
        EGPI.StopWeightMinusErr = EGPI.StopWeight - EGPI.StopWeightErr;
        cout << "stopWeight " << EGPI.StopWeight << endl;
        cout << "stopWeightErr " << EGPI.StopWeightErr << endl;
        cout << "stopWeightErr " << EGPI.StopWeightPlusErr << endl;
        cout << "stopWeightErr " << EGPI.StopWeightMinusErr << endl;
    }
    ////Book histograms and histogram names
    outputFile->cd();
    vector<HistogramT> * histVec_1D = OneDeeHistTVec();
    vector<HistogramT> * histVec_2D = TwoDeeHistTVec();
    vector<HistogramT> * histVec_3D = ThreeDeeHistTVec();
    vector<HistogramT> * histVec_1D_Smear, * histVec_1D_Smear_Syst;
    /*
    if (!BEI.doData) {
        histVec_1D_Smear = OneDeeHistTVecSmearJets();
        histVec_1D_Smear_Syst = new vector<HistogramT>;
    }
    */
    vector<SampleT> * subSampVec    = SubSampVec();    ///Define things necessary for booking histograms
    vector<SystT> * systVec         = SystVec();
    vector<HistogramT> * histVec_1D_Syst = new vector<HistogramT>;
    vector<HistogramT> * histVec_2D_Syst = new vector<HistogramT>;
    vector<HistogramT> * histVec_3D_Syst = new vector<HistogramT>;
    HistogramT H_Current; 
    TH1F * h_1DCurr; TH2F * h_2DCurr; TH3F * h_3DCurr;
    HMap_1D histMap_1D; HMap_2D histMap_2D; HMap_3D histMap_3D; passCutMap subSampBool;
    passCutMap subSampBool_LepESUp, subSampBool_LepESDown;
    passCutMap subSampBool_JetESUp, subSampBool_JetESDown;
    passCutMap subSampBool_BTagSFUp, subSampBool_BTagSFDown;
    SampleT S_Current;
    TString histTitle;
    TString axesTitle;
    float nXBins, nYBins, nZBins;
    float xBinMin, xBinMax;
    float yBinMin, yBinMax;
    float zBinMin, zBinMax;
    if (doBookSyst) {
        histVec_1D_Syst = AddSystHists(histVec_1D, systVec, fInName, BEI.isSignal);
        
        histVec_2D_Syst = AddSystHists(histVec_2D, systVec, fInName, BEI.isSignal);        
        /*
         histVec_3D_Syst = AddSystHists(histVec_3D, systVec, fInName, BEI.isSignal);
         */
//        histVec_1D_Smear_Syst = AddSystHists(histVec_1D_Smear, systVec, fInName, BEI.isSignal);
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
        map<string, float> stringKeyToVar;
        
        ELI.ELIDefaultVarVals();
        EJI.EJIDefaultVarVals();
        EMI.EMIDefaultVarVals();
        EDSI.EDSIDefaultVarVals();

        if (!BEI.doData) {
            EGPI.EGPIDefaultVarVals();
            ELI_LepESUp.ELIDefaultVarVals();
            ELI_LepESDown.ELIDefaultVarVals();
            
            EJI_JetESUp.EJIDefaultVarVals();
            EJI_JetESDown.EJIDefaultVarVals();
            EJI_BTagSFUp.EJIDefaultVarVals();
            EJI_BTagSFDown.EJIDefaultVarVals();
            
            EMI_LepESUp.EMIDefaultVarVals();
            EMI_LepESDown.EMIDefaultVarVals();
            EMI_JetESUp.EMIDefaultVarVals();
            EMI_JetESDown.EMIDefaultVarVals();
            EMI_BTagSFUp.EMIDefaultVarVals();
            EMI_BTagSFDown.EMIDefaultVarVals();
            
            EDSI_LepESUp.EDSIDefaultVarVals();
            EDSI_LepESDown.EDSIDefaultVarVals();
            EDSI_JetESUp.EDSIDefaultVarVals();
            EDSI_JetESDown.EDSIDefaultVarVals();
            EDSI_BTagSFUp.EDSIDefaultVarVals();
            EDSI_BTagSFDown.EDSIDefaultVarVals();
        }
        
        
        fileTree.GetEntry(ievt);
//        cout << "test 1 " << endl;
        if (whichNTupleType == 0) {
            ELI.doEvent = EventPassTrigger(&BEI, &ELI);
            if (!BEI.doData) {
                ELI_LepESUp.doEvent = EventPassTrigger(&BEI, &ELI_LepESUp);
                ELI_LepESDown.doEvent = EventPassTrigger(&BEI, &ELI_LepESDown);
                keepEvent = (ELI.doEvent || ELI_LepESUp.doEvent || ELI_LepESDown.doEvent);
            }
            else {
                keepEvent = ELI.doEvent;
            }
        }
        else {
            keepEvent = true;
        }
        if (!keepEvent) continue;
//        cout << "test 2 " << endl;
        if (!BEI.doData) {
            EGPI.weight_GenTopReweight = 1.; // set it to 1 for default -- will switch later if "hasTopInfo"
            if (BEI.isSignal) {
                if (!BEI.hasStopInfo) continue;
                if (fabs(TGenStopMass0 - grabStopMass) > massDiffThresh) continue;
                if (fabs(TGenStopMass1 - grabStopMass) > massDiffThresh) continue;
                if (fabs(TGenChi0Mass0 - grabChi0Mass) > massDiffThresh) continue;
                if (fabs(TGenChi0Mass1 - grabChi0Mass) > massDiffThresh) continue;
                if (fabs(TGenCharginoMass0 - grabCharginoMass) > massDiffThresh) continue;
                if (fabs(TGenCharginoMass1 - grabCharginoMass) > massDiffThresh) continue;
            }
            if (BEI.hasMETInfo) {
                genMETVec.SetPtEtaPhiM(genMET_Pt, 0., genMET_Phi, 0.);
                EGPI.genMET.P4 = genMETVec;
                EGPI.genMET.PDGID = 0;
                EGPI.genMET.PDGStatus = 0;
                EGPI.genMET.Mass = 0.;
            }
            if (BEI.hasTopInfo) {
                genTop0Vec.SetPtEtaPhiE(genTop0Pt, genTop0Eta, genTop0Phi, genTop0En);
                EGPI.genTop0.P4 = genTop0Vec;
                EGPI.genTop0.PDGID = genTop0PdgId;
                EGPI.genTop0.PDGStatus = 3;
                EGPI.genTop0.Mass = genTop0Vec.M();                
                genTop1Vec.SetPtEtaPhiE(genTop1Pt, genTop1Eta, genTop1Phi, genTop1En);
                EGPI.genTop1.P4 = genTop1Vec;
                EGPI.genTop1.PDGID = genTop1PdgId;
                EGPI.genTop1.PDGStatus = 3;
                EGPI.genTop1.Mass = genTop1Vec.M();                
                if (genTop0PdgId * genTop1PdgId > 0) {
                    cout << "something is funky with the genTop PDGIDs" << endl;
                    continue;
                }
                else {
                    if (genTop0PdgId > 0) {
                        EGPI.genTopPt = genTop0Pt;
                        EGPI.genAntiTopPt = genTop1Pt;
                    }
                    else {
                        EGPI.genTopPt = genTop1Pt;
                        EGPI.genAntiTopPt = genTop0Pt;   
                    }                
                }
                EGPI.weight_GenTopReweight = GenLevelTopPtWeight(EGPI.genTopPt, EGPI.genAntiTopPt);
                h_genTopPt->Fill(EGPI.genTopPt, 1.);
                h_genTopPtRW->Fill(EGPI.genTopPt, EGPI.weight_GenTopReweight);
                h_genAntiTopPt->Fill(EGPI.genAntiTopPt, 1.);
                h_genAntiTopPtRW->Fill(EGPI.genAntiTopPt, EGPI.weight_GenTopReweight);
            }
            if (BEI.hasStopInfo) {
                EGPI.genStop0.Mass = TGenStopMass0;
                EGPI.genStop1.Mass = TGenStopMass1;
                EGPI.genChiZero0.Mass = TGenChi0Mass0;
                EGPI.genChiZero1.Mass = TGenChi0Mass1;
                EGPI.genChargino0.Mass = TGenCharginoMass0;
                EGPI.genChargino1.Mass = TGenCharginoMass1;                
            }
        }
        ///****************************************
        // Call the code to set up the event information
        ///****************************************
        
        BEI.preNVtxRWWeight = BEI.weight;
        if (!BEI.doData) {
            BEI.weight = nVtxWeight(&BEI, nVtxSFHistOviToDESY, nVtxSFHist_v2, h_S7toS10RWHist, nVtxSFHist);
            BEI_GenTopRW.SetVars(&BEI); // i.e. no GenTopReweighting systematic.
//            cout << "EGPI.weight_GenTopReweight " << EGPI.weight_GenTopReweight << endl;
            BEI.weight *= EGPI.weight_GenTopReweight;
            BEI.preNVtxRWWeight *= EGPI.weight_GenTopReweight;
            BEI_LepESUp.SetVars(&BEI);
            BEI_LepESDown.SetVars(&BEI);
            BEI_StopXSecUp.SetVars(&BEI);
            BEI_StopXSecDown.SetVars(&BEI);
            BEI_LepEffSFUp.SetVars(&BEI);
            BEI_LepEffSFDown.SetVars(&BEI);
        }
        SetEventInformation(&BEI, &ELI, &EJI, EMI, EDSI);

        ///****************************************
        // Call the code to set up the event information
        ///****************************************        
        if (ELI.Lep1.P4.Pt() < subLepPtCut) {
            ELI.doEvent = false;
        }
        if (ELI.EventDiLepType == -2) ELI.EventDiLepType = 2;
        if (BEI.doData) {
            if (!ELI.doEvent) continue;
//            cout << "test 3 " << endl;
        }
        else {            
            if (doVerbosity) {
                cout << "Type " << ELI.EventDiLepType << endl;
                cout << "weight pre scale " << BEI.weight << endl;
            }
            SetEventInformation(&BEI, &ELI_LepESUp,   &EJI,            EMI_LepESUp,    EDSI_LepESUp);
            SetEventInformation(&BEI, &ELI_LepESDown, &EJI,            EMI_LepESDown,  EDSI_LepESDown);
            SetEventInformation(&BEI, &ELI,           &EJI_JetESUp,    EMI_JetESUp,    EDSI_JetESUp);
            SetEventInformation(&BEI, &ELI,           &EJI_JetESDown,  EMI_JetESDown,  EDSI_JetESDown);
            EJI_BTagSFUp = EJI;
            EJI_BTagSFDown = EJI;
            SetEventInformation(&BEI, &ELI,           &EJI_BTagSFUp,   EMI_BTagSFUp,   EDSI_BTagSFUp);
            SetEventInformation(&BEI, &ELI,           &EJI_BTagSFDown, EMI_BTagSFDown, EDSI_BTagSFDown);
            if (ELI_LepESUp.Lep1.P4.Pt() < subLepPtCut) {
                ELI_LepESUp.doEvent = false;
            }
            if (ELI_LepESDown.Lep1.P4.Pt() < subLepPtCut) {
                ELI_LepESDown.doEvent = false;
            }
            if (ELI_LepESUp.EventDiLepType == -2) ELI_LepESUp.EventDiLepType = 2;
            if (ELI_LepESDown.EventDiLepType == -2) ELI_LepESDown.EventDiLepType = 2;
            if (!ELI.doEvent && !ELI_LepESUp.doEvent && !ELI_LepESDown.doEvent) continue;
//            cout << " for ievt " << ievt << endl;
//            cout << "Type " << ELI.EventDiLepType << endl;
//            cout << "weight pre scale " << BEI.weight << endl;
            BEI.weight *= ScaleFactorMC(ELI.EventDiLepType, 0);
            BEI.preNVtxRWWeight *= ScaleFactorMC(ELI.EventDiLepType, 0);
//            cout << "Type " << ELI.EventDiLepType << endl;
//            cout << "weight post scale " << BEI.weight << endl;
            
            
            BEI_GenTopRW.weight *= ScaleFactorMC(ELI.EventDiLepType, 0);
            BEI_GenTopRW.preNVtxRWWeight *= ScaleFactorMC(ELI.EventDiLepType, 0);
            
            BEI_LepESUp.weight *= ScaleFactorMC(ELI_LepESUp.EventDiLepType, 0);
            BEI_LepESUp.preNVtxRWWeight *= ScaleFactorMC(ELI_LepESUp.EventDiLepType, 0);
            BEI_LepESDown.weight *= ScaleFactorMC(ELI_LepESDown.EventDiLepType, 0);
            BEI_LepESDown.preNVtxRWWeight *= ScaleFactorMC(ELI_LepESDown.EventDiLepType, 0);
            
            BEI_LepEffSFUp.weight *= ScaleFactorMC(ELI.EventDiLepType, 1);
            BEI_LepEffSFUp.preNVtxRWWeight *= ScaleFactorMC(ELI.EventDiLepType, 1);            
            BEI_LepEffSFDown.weight *= ScaleFactorMC(ELI.EventDiLepType, -1);
            BEI_LepEffSFDown.preNVtxRWWeight *= ScaleFactorMC(ELI.EventDiLepType, -1);
            
            BEI_StopXSecUp.weight *= ScaleFactorMC(ELI.EventDiLepType, 0);
            BEI_StopXSecUp.preNVtxRWWeight *= ScaleFactorMC(ELI.EventDiLepType, 0);
            BEI_StopXSecDown.weight *= ScaleFactorMC(ELI.EventDiLepType, 0);
            BEI_StopXSecDown.preNVtxRWWeight *= ScaleFactorMC(ELI.EventDiLepType, 0);

            if (BEI.isSignal) {
                BEI.weight *= EGPI.StopWeight;
                BEI.preNVtxRWWeight *= EGPI.StopWeight;
                
                BEI_LepEffSFUp.weight *= EGPI.StopWeight;
                BEI_LepEffSFUp.preNVtxRWWeight *= EGPI.StopWeight;
                BEI_LepEffSFDown.weight *= EGPI.StopWeight;
                BEI_LepEffSFDown.preNVtxRWWeight *= EGPI.StopWeight;
                
                BEI_StopXSecUp.weight *= EGPI.StopWeightPlusErr;
                BEI_StopXSecUp.preNVtxRWWeight *= EGPI.StopWeightPlusErr;                
                BEI_StopXSecDown.weight *= EGPI.StopWeightMinusErr;
                BEI_StopXSecDown.preNVtxRWWeight *= EGPI.StopWeightMinusErr;
                
                BEI_LepESUp.weight *= EGPI.StopWeight;
                BEI_LepESUp.preNVtxRWWeight *= EGPI.StopWeight;
                BEI_LepESDown.weight *= EGPI.StopWeight;
                BEI_LepESDown.preNVtxRWWeight *= EGPI.StopWeight;
                
                BEI_GenTopRW.weight *= EGPI.StopWeight;
                BEI_GenTopRW.preNVtxRWWeight *= EGPI.StopWeight;
            }         
        }
        // basic condition, if running on data, only select type of events that are relevant to prevent double counting
        if (BEI.doData) {
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
//        cout << "test 4 " << endl;
        ESMT2I.SetVars(&EMI, &EJI, MT2llMeanSmear, MT2llUncEnUpDelta2D, MT2llUncEnDownDelta2D, MT2lbUncEnUpDelta2D, MT2lbUncEnDownDelta2D, vecOneDeeMT2llUncEnUp, vecOneDeeMT2llUncEnDown, vecOneDeeMT2lbUncEnUp, vecOneDeeMT2lbUncEnDown);
        
        /******************************************************/
        ///Set up the mapping of string keys to the appropriate event variables
        /*
         list of keys needed taken from command:
         cat StopFunctionDefinitions.h | grep VarKey
         _______________________________
         -----------------------------
         */
        stringKeyToVar["CutFlowEntry"] = 4;        
        stringKeyToVar["nVtx"]      = (float) BEI.nVtx;
        stringKeyToVar["nVtxTrue"]  = (float) BEI.nVtxTrue;
        
        SetStringKeyMap(stringKeyToVar, &ELI, &EJI, &EMI, &EDSI, systVec);
        
        if (doBookSyst) {
            SetStringKeyMap(stringKeyToVar, &ELI_LepESUp,    &EJI,           &EMI_LepESUp,   &EDSI_LepESUp,   systVec,  1);
            SetStringKeyMap(stringKeyToVar, &ELI_LepESDown,  &EJI,           &EMI_LepESDown, &EDSI_LepESDown, systVec, -1);            
            SetStringKeyMap(stringKeyToVar, &ELI,            &EJI_JetESUp,   &EMI_JetESUp,   &EDSI_JetESUp,   systVec,  2);
            SetStringKeyMap(stringKeyToVar, &ELI,            &EJI_JetESDown, &EMI_JetESDown, &EDSI_JetESDown, systVec, -2);            
            SetStringKeyMap(stringKeyToVar, &ELI,            &EJI_BTagSFUp,  &EMI_BTagSFUp,  &EDSI_BTagSFUp,  systVec,  3);
            SetStringKeyMap(stringKeyToVar, &ELI,            &EJI_BTagSFUp,  &EMI_BTagSFUp,  &EDSI_BTagSFUp,  systVec, -3);
            SetStringKeyMapSpecial(stringKeyToVar, &EMI, &ESMT2I, systVec, 5);
            SetStringKeyMapSpecial(stringKeyToVar, &EMI, &ESMT2I, systVec, -5);
            SetStringKeyMapSpecial(stringKeyToVar, &EMI, &ESMT2I, systVec, 8);
        }

        /*#######################
         MAKE SELECTION CUTS
         #####################*/
        SetPassCutMap(subSampBool, subSampVec, &ELI, &EJI, &EMI);
        /*
        for (unsigned int iSamp = 0; iSamp < subSampVec->size(); ++iSamp) {
            
            cout << "inputCut Map is " << subSampBool[subSampVec->at(iSamp)] << " for S_Current = " << subSampVec->at(iSamp).histNameSuffix << endl;
        }
        */
        if (doBookSyst) {
            SetPassCutMap(subSampBool_LepESUp,    subSampVec, &ELI_LepESUp,   &EJI,            &EMI_LepESUp);
            SetPassCutMap(subSampBool_LepESDown,  subSampVec, &ELI_LepESDown, &EJI,            &EMI_LepESDown);            
            SetPassCutMap(subSampBool_JetESUp,    subSampVec, &ELI,           &EJI_JetESUp,    &EMI_JetESUp);
            SetPassCutMap(subSampBool_JetESDown,  subSampVec, &ELI,           &EJI_JetESDown,  &EMI_JetESDown);            
            SetPassCutMap(subSampBool_BTagSFUp,   subSampVec, &ELI,           &EJI_BTagSFUp,   &EMI_BTagSFUp);
            SetPassCutMap(subSampBool_BTagSFDown, subSampVec, &ELI,           &EJI_BTagSFDown, &EMI_BTagSFDown);
        }
        HistogramFillOneDee(&subSampBool, &stringKeyToVar, subSampVec, histVec_1D, &histMap_1D, &BEI, &ELI, &EJI, &EMI, &EDSI, doVerbosity);
        HistogramFillTwoDee(&subSampBool, &stringKeyToVar, subSampVec, histVec_2D, &histMap_2D, &BEI, &ELI, &EJI, &EMI, &EDSI, doVerbosity);
        HistogramFillThreeDee(&subSampBool, &stringKeyToVar, subSampVec, histVec_3D, &histMap_3D, &BEI, &ELI, &EJI, &EMI, &EDSI, doVerbosity);
        if (!BEI.doData) {
            HistogramFillOneDeeSyst(&subSampBool_LepESUp, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI_LepESUp, &ELI_LepESUp, &EJI, &EMI_LepESUp, &EDSI_LepESUp, "LepESShiftUp", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool_LepESUp, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI_LepESUp, &ELI_LepESUp, &EJI, &EMI_LepESUp, &EDSI_LepESUp, "LepESShiftUp", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool_LepESUp, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI_LepESUp, &ELI_LepESUp, &EJI, &EMI_LepESUp, &EDSI_LepESUp, "LepESShiftUp", doVerbosity);            
            HistogramFillOneDeeSyst(&subSampBool_LepESDown, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI_LepESDown, &ELI_LepESDown, &EJI, &EMI_LepESDown, &EDSI_LepESDown, "LepESShiftDown", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool_LepESDown, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI_LepESDown, &ELI_LepESDown, &EJI, &EMI_LepESDown, &EDSI_LepESDown, "LepESShiftDown", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool_LepESDown, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI_LepESDown, &ELI_LepESDown, &EJI, &EMI_LepESDown, &EDSI_LepESDown, "LepESShiftDown", doVerbosity);
            
            HistogramFillOneDeeSyst(&subSampBool_JetESUp, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI, &ELI, &EJI_JetESUp, &EMI_JetESUp, &EDSI_JetESUp, "JetESShiftUp", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool_LepESDown, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI, &ELI, &EJI_JetESUp, &EMI_JetESUp, &EDSI_JetESUp, "JetESShiftUp", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool_JetESUp, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI, &ELI, &EJI_JetESUp, &EMI_JetESUp, &EDSI_JetESUp, "JetESShiftUp", doVerbosity);            
            HistogramFillOneDeeSyst(&subSampBool_JetESDown, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI, &ELI, &EJI_JetESDown, &EMI_JetESDown, &EDSI_JetESDown, "JetESShiftDown", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool_LepESDown, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI, &ELI, &EJI_JetESDown, &EMI_JetESDown, &EDSI_JetESDown, "JetESShiftDown", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool_JetESDown, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI, &ELI, &EJI_JetESDown, &EMI_JetESDown, &EDSI_JetESDown, "JetESShiftDown", doVerbosity);
            
            HistogramFillOneDeeSyst(&subSampBool_BTagSFUp, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI, &ELI, &EJI_BTagSFUp, &EMI_BTagSFUp, &EDSI_BTagSFUp, "BTagSFShiftUp", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool_BTagSFUp, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI, &ELI, &EJI_BTagSFUp, &EMI_BTagSFUp, &EDSI_BTagSFUp, "BTagSFShiftUp", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool_BTagSFUp, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI, &ELI, &EJI_BTagSFUp, &EMI_BTagSFUp, &EDSI_BTagSFUp, "BTagSFShiftUp", doVerbosity);
            HistogramFillOneDeeSyst(&subSampBool_BTagSFDown, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI, &ELI, &EJI_BTagSFDown, &EMI_BTagSFDown, &EDSI_BTagSFDown, "BTagSFShiftDown", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool_BTagSFDown, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI, &ELI, &EJI_BTagSFDown, &EMI_BTagSFDown, &EDSI_BTagSFDown, "BTagSFShiftDown", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool_BTagSFDown, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI, &ELI, &EJI_BTagSFDown, &EMI_BTagSFDown, &EDSI_BTagSFDown, "BTagSFShiftDown", doVerbosity);
            
            HistogramFillOneDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI_LepEffSFUp, &ELI, &EJI, &EMI, &EDSI, "LepEffSFShiftUp", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI_LepEffSFUp, &ELI, &EJI, &EMI, &EDSI, "LepEffSFShiftUp", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI_LepEffSFUp, &ELI, &EJI, &EMI, &EDSI, "LepEffSFShiftUp", doVerbosity);            
            HistogramFillOneDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI_LepEffSFDown, &ELI, &EJI, &EMI, &EDSI, "LepEffSFShiftDown", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI_LepEffSFDown, &ELI, &EJI, &EMI, &EDSI, "LepEffSFShiftDown", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI_LepEffSFDown, &ELI, &EJI, &EMI, &EDSI, "LepEffSFShiftDown", doVerbosity);
            
            HistogramFillOneDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI_GenTopRW, &ELI, &EJI, &EMI, &EDSI, "genTopRW", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI_GenTopRW, &ELI, &EJI, &EMI, &EDSI, "genTopRW", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI_GenTopRW, &ELI, &EJI, &EMI, &EDSI, "genTopRW", doVerbosity);           
            
            HistogramFillOneDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI_StopXSecUp, &ELI, &EJI, &EMI, &EDSI, "genStopXSecShiftUp", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI_StopXSecUp, &ELI, &EJI, &EMI, &EDSI, "genStopXSecShiftUp", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI_StopXSecUp, &ELI, &EJI, &EMI, &EDSI, "genStopXSecShiftUp", doVerbosity);
            HistogramFillOneDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_1D_Syst, &histMap_1D, &BEI_StopXSecDown, &ELI, &EJI, &EMI, &EDSI, "genStopXSecShiftDown", doVerbosity);
            HistogramFillTwoDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_2D_Syst, &histMap_2D, &BEI_StopXSecDown, &ELI, &EJI, &EMI, &EDSI, "genStopXSecShiftDown", doVerbosity);
            HistogramFillThreeDeeSyst(&subSampBool, &stringKeyToVar, subSampVec, histVec_3D_Syst, &histMap_3D, &BEI_StopXSecDown, &ELI, &EJI, &EMI, &EDSI, "genStopXSecShiftDown", doVerbosity);                                    
        }
        /*
        for (unsigned int i = 0; i < subSampVec->size(); ++i) {
            S_Current = subSampVec->at(i);
            if (subSampBool[S_Current] == true) {                            
                if (!doData) {
                    for (unsigned int js = 0; js < histVec_1D_Syst->size(); ++js) {
                        H_Current = histVec_1D_Syst->at(js);
                        if (H_Current.name.Contains("LepESShift") || H_Current.name.Contains("JetESShift") || H_Current.name.Contains("BTagSFShift")) continue;
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
                                if (dPhi(ELI.Lep0.P4.Phi(), ELI.Lep1.P4.Phi()) > 1./3. * PI) continue;
                            }
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Mid")) {
                                if (dPhi(ELI.Lep0.P4.Phi(), ELI.Lep1.P4.Phi()) < 1./3. * PI || dPhi(ELI.Lep0.P4.Phi(), ELI.Lep1.P4.Phi()) > 2./3. * PI) continue;                         
                            }                        
                            else if (TString(H_Current.name).Contains("MT2ll_DPhiLep0Lep1Far")) {
                                if (dPhi(ELI.Lep0.P4.Phi(), ELI.Lep1.P4.Phi()) < 2./3. * PI) continue;
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
                        if (H_Current.name.Contains("LepESShift") || H_Current.name.Contains("JetESShift") || H_Current.name.Contains("BTagSFShift")) continue;
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
                        if (H_Current.name.Contains("LepESShift") || H_Current.name.Contains("JetESShift")|| H_Current.name.Contains("BTagSFShift")) continue;
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
        }        
        */
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