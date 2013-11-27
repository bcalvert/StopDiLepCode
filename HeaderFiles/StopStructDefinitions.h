
#include "TLorentzVector.h"
#include "TRandom.h"
#include <iostream>
#include <vector>
#include <TROOT.h>
#include <iostream>
#include <fstream>
//#include <vector>
#include <cmath>
#include <sstream>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

using namespace std;
typedef struct {
    TString name;
    float   fVar;
    double  dVar;
    string  systVarKey;
    int     whichSystType;   // 0 = universal systematic, 1 = lepton systematic, 2 = jet systematic, 3 = other systematic
} SystT;

typedef struct HistAxis {
    int AxisBinN;
    int AxisRBN;
    float AxisMin, AxisMax;
    string AxisVarKey;
    bool AxisDoSyst;
    float AxisLB, AxisUB;
    void ClearVars() {
        AxisBinN = 0;
        AxisRBN = 0;
        AxisLB = -99999.;
        AxisLB = 99999.;
    }
} HistAxis;

typedef struct {
    TString name;
    TString stringAxis;
    TString xLabel;
    int xBinN;
    int RBNX;
    float xMin, xMax;
    string xVarKey;
    bool doXSyst;
    float xLB, xUB;
    
    TString yLabel;
    int yBinN;
    int RBNY;
    float yMin, yMax; 
    string yVarKey;
    bool doYSyst;
    float yLB, yUB;
    
    TString zLabel;
    int zBinN;
    float zMin, zMax; 
    int RBNZ;
    string zVarKey;
    bool doZSyst;
    float zLB, zUB;
    
    bool logY1D;
    
    void SetAxisString(TString XString, TString YString, TString ZString, TString titleString = "") {
        stringAxis = titleString;
        stringAxis += ";";
        stringAxis += XString;
        stringAxis += ";";
        stringAxis += YString;
        stringAxis += ";";
        stringAxis += ZString;
        stringAxis += ";";
    }
    void SetName(TString XString, TString YString, TString ZString, TString AppendString = "", int numDims = 1) {
        name = "h_";
        if (numDims > 0) {
            name += XString;
            if (numDims > 1) {
                name += "_vs_";
                name += YString;
                if (numDims > 2) {
                    name += "_vs_";
                    name += ZString;
                }
            }
        }
        name += AppendString;
    }
    void SetIndAxisLabel(TString inputString, map<TString, TString> * mapVartoLabel, int whichDim = 1) {
        map<TString, TString>::iterator xIter;
        xIter = mapVartoLabel->find(inputString);
        if (xIter != mapVartoLabel->end()) {
            if (whichDim == 1) {
                xLabel = xIter->second;
            }
            else if (whichDim == 2) {
                yLabel = xIter->second;
            }
            else if (whichDim == 3) {
                zLabel = xIter->second;
            }
        }
    }
} HistogramT;

typedef struct {
    //    SampleT() : histNameSuffix(""), whichdiLepType(-1), doZVeto(-1), cutNJets(-1), cutNBJets(-1), cutMET(0) {}
    TString histNameSuffix;
    TString histXaxisSuffix;
    TString histYaxisSuffix;
    TString histZaxisSuffix;
    /// variables to store what kinds of
    int     whichdiLepType; // -1: inclusive, 0: MuMu, 1: EE, 2: EMu
    int     doZVeto;   // -1: inclusive, 0: ZMass window, 1: outside ZMass window;
    int     cutNJets;  // -1: inclusive, # > 0: require NJets >= #
    int     cutNBJets;  // -1: inclusive, # > 0: require NBJets >= #
    float   cutMET; // # > 0: require MET >= #
    bool    blindDataChannel;
    
} SampleT;

typedef struct {
    TString HistTXLabel;
    TString HistTYLabel;
    TString HistTZLabel;
    TString SampTLabel;
    
    int newXBinN;
    float newXMin, newXMax;
    
    int newYBinN;
    float newYMin, newYMax;
    
    int newZBinN;
    float newZMin, newZMax;    
    
} SpecHistBinT;

typedef std::pair<HistogramT, SampleT> histKey;
inline bool operator<(const histKey &a, const histKey &b)
{
    return (a.first.name < b.first.name) || (a.first.name == b.first.name && a.second.histNameSuffix < b.second.histNameSuffix);
}
inline bool operator<(const SampleT &a, const SampleT &b)
{
    return (a.histNameSuffix < b.histNameSuffix) || (a.histXaxisSuffix < b.histXaxisSuffix && a.histNameSuffix == b.histNameSuffix);
}
inline SampleT operator+(const SampleT &a, const SampleT &b)
{
    SampleT outSampleT;
    outSampleT.histNameSuffix = a.histNameSuffix + b.histNameSuffix;
    outSampleT.histXaxisSuffix = a.histXaxisSuffix + b.histXaxisSuffix;
    outSampleT.histYaxisSuffix = a.histYaxisSuffix + b.histYaxisSuffix;
    outSampleT.histZaxisSuffix = a.histZaxisSuffix + b.histZaxisSuffix;
    return outSampleT;
}


//typedef pair<HistogramT, SampleT> histKey;


inline float DeltaMT2UncEn(std::vector<TH1F *> * vecOneDeeHists, TH2F * TwoDeeHist, float inputMT2Value) {
    //returns float which is a random number drawn from the distribution of MT2 Central Value minus MT2 Unclustered ES shifted version
    
    int whichOneDeeHist = TwoDeeHist->GetXaxis()->FindBin(inputMT2Value) - 1;
    unsigned int numOneDeeHists = vecOneDeeHists->size();
    if (whichOneDeeHist < 0) {
        std::cout << inputMT2Value << std::endl;
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

inline float dPhi(float phi1, float phi2) {
    float result = phi1-phi2;
    while (result >= TMath::Pi()) result -= 2*TMath::Pi();
    while (result < -1*TMath::Pi()) result += 2*TMath::Pi();
    return fabs(result);
}
inline float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float dphi = dPhi(phi1,phi2);
    float deta = eta1-eta2;
    float result = dphi*dphi+deta*deta;
    result = sqrt(result);
    return result;
}

inline double dPhi(double phi1, double phi2) {
    double result = phi1-phi2;
    while (result >= TMath::Pi()) result -= 2*TMath::Pi();
    while (result < -1*TMath::Pi()) result += 2*TMath::Pi();
    return fabs(result);
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dphi = dPhi(phi1,phi2);
    double deta = eta1-eta2;
    double result = dphi*dphi+deta*deta;
    result = sqrt(result);
    return result;
}
typedef struct PFJetEventPointers {
    std::vector<float> * JetPx, * JetPy, * JetPz, * JetE, * JetNHF, * JetNEF, * JetCHF, * JetCEF, * JetBTag;
    std::vector<int> * JetNDaug, * JetCharMult, * JetPartFlav;
    unsigned int numPFJets;
    //    std::vector<float> * genJetPx, * genJetPy, * genJetPz, * genJetEn, * genJetInvisE; //* genJetEt, * genJetEta;
    //    std::vector<bool> * genJetIsGenJet;
}PFJetEventPointers;

typedef struct Lepton{
    TLorentzVector  P4;
    int             PDGID;
    float           relPFLepIso;
    void ClearVars() {
        P4.SetPxPyPzE(0., 0., 0., 0.);
        PDGID = 0;
        relPFLepIso = -1;
    }
    void isBadLep() {
        P4.SetPxPyPzE(-99999., -99999., -99999., -99999.);
        PDGID = 0;
        relPFLepIso = -99999;
    }
    bool isElec() {
        return (abs(PDGID) == 11);
    }
    bool isMuon() {
        return (abs(PDGID) == 13);
    }
    /*
     Lepton(TLorentzVector input4Vec, int pdgid) : P4(input4Vec), PDGID(pdgid) {}
     bool operator < (const Lepton& otherLep) const
     {
     return (P4.Pt() < Lepton.P4.Pt());
     }
     */
} Lepton;

inline bool operator<(const Lepton &a, const Lepton &b)
{
    return (a.P4.Pt() < b.P4.Pt());
}

inline bool operator>(const Lepton &a, const Lepton &b)
{
    return (a.P4.Pt() > b.P4.Pt());
}

typedef struct MuonEventPointers {
    std::vector<float> * MuonPt, * MuonPx, * MuonPy, * MuonPz, * MuonEn, * MuonPFCharHadIso, * MuonPFNeutHadIso, * MuonPFPhotIso, * MuonSumPUPt, * MuonD0, * MuonVertZ;
    std::vector<float> * PFMuonPt;
    std::vector<bool> * isPFMuon, * isGMPTMuon, * isGlobMuon; //, * isTrackArbitMuon;
    std::vector<int> * MuonNumMatchStations, * MuonNumLayers, * MuonNumValidPixHitsinTrack;
    std::vector<int> * MuonCharge;
    unsigned int numMuons;
        //ambiguity in which Muon Pt...rolling with PFMuonPt...never mind ambiguity resolved see https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleLepton2012#Muon_Selection
}MuonEventPointers;

typedef struct ElectronEventPointers {
    std::vector<float> * ElecPt, * ElecPx, * ElecPy, * ElecPz, * ElecEn, * ElecPFCharHadIso, * ElecPFNeutHadIso, * ElecPFPhotIso; 
    std::vector<float> * PFElecPt;
    std::vector<float> * ElecSCEta, * ElecDeltaPhiIn, * ElecDeltaEtaIn, * ElecSigIetaIeta, * ElecHtoERatio, * ElecIP, * ElecDZ, * ElecECalE, * ElecSCEOverP;
    std::vector<int> * ElecNumMissHits;
    std::vector<bool> * ElecisEB, * ElecisEE;
    std::vector<int> * ElecCharge;
    //ambiguity in which Electron Pt...rolling with PFElecPt...never mind ambiguity resolved    
    std::vector<bool> * isPFElectron, * passConvVeto;
    unsigned int numElectrons;
}ElectronEventPointers;


typedef struct JetCutInfo {
    // Look here for CutID definitions: https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    float cutNHF, cutNEF, cutCHF, cutCEF;
    int cutNConst, cutChMult;
    float cutBTagDisc, cutJetPt, cutJetEta;
    void SetJetCutVars(int cutStrength) {
        // cutStrength 0: Loose, 1: Medium, 2: Tight
        float cutNHFStrength[3] = {0.99, 0.95, 0.90};
        float cutNEFStrength[3] = {0.99, 0.95, 0.90};
        float cutCHFStrength[3] = {0., 0., 0.};
        float cutCEFStrength[3] = {0.99, 0.99, 0.99};
        int cutNConstStrength[3] = {1, 1, 1};
        int cutChMultStrength[3] = {0, 0, 0};
        cutNHF = cutNHFStrength[cutStrength];
        cutNEF = cutNEFStrength[cutStrength];
        cutCHF = cutCHFStrength[cutStrength];
        cutCEF = cutCEFStrength[cutStrength];
        cutNConst = cutNConstStrength[cutStrength];
        cutChMult = cutChMultStrength[cutStrength];
    }
    void DefaultCutVarVals() {
        SetJetCutVars(0);
        cutJetPt = 30.;
        cutJetEta = 2.4;
        cutBTagDisc = 0.679;  //CSV Middle working point, see (remove underscore in address): h_ttps://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    }
    void PrintVals() {
        std::cout << "cutNHF " << cutNHF << std::endl;
        std::cout << "cutNEF " << cutNEF << std::endl;
        std::cout << "cutCHF " << cutCHF << std::endl;
        std::cout << "cutCEF " << cutCEF << std::endl;
        std::cout << "cutNConst " << cutNConst << std::endl;
        std::cout << "cutChMult " << cutChMult << std::endl;
        std::cout << "cutBTagDisc " << cutBTagDisc << std::endl;
        std::cout << "cutJetPt " << cutJetPt << std::endl;
        std::cout << "cutJetEta " << cutJetEta << std::endl;
    }
} JetCutInfo;

typedef struct PFJet{
    TLorentzVector  P4;
    float           valBTagDisc;
    int             partonFlavor;
    bool            isBJet;
    bool            isGenJetMatched;
    float           dEnRecoGen;
    float NeutHadFrac, ChargeHadFrac, NeutEMFrac, ChargeEMFrac;
    int   NumDaughters, ChargeMult;
    bool  passJetID;
    void ClearVars() {
        P4.SetPxPyPzE(1E-8, 1E-8, 1E-4, 1E-8);
        valBTagDisc = 0.0;
        partonFlavor = -999999;
        isBJet = false;
        isGenJetMatched = false;
        dEnRecoGen      = 0.0;
        NeutHadFrac     = 0.0;
        ChargeHadFrac   = 0.0;
        NeutEMFrac      = 0.0;
        ChargeEMFrac    = 0.0;
        NumDaughters    = 0;
        ChargeMult      = 0;
        passJetID       = true;
    }
    void isBadJet() {
        P4.SetPxPyPzE(-99999., -99999., -99999., -99999.);
        valBTagDisc = -1;
        partonFlavor = -999999;
        isBJet = false;
        isGenJetMatched = false;
        dEnRecoGen = 0.0;
    }
    bool PassesID(JetCutInfo * inJCI) {
        bool passJetID = false;
        passJetID = (NeutHadFrac < inJCI->cutNHF && NeutEMFrac < inJCI->cutNEF);
//        std::cout << "passJetIDStage 1? " << passJetID << std::endl;
        if (fabs(P4.Eta()) < inJCI->cutJetEta) {
//            std::cout << "P4.Eta() " << P4.Eta() << std::endl;
            passJetID &= (ChargeEMFrac < inJCI->cutCEF && ChargeHadFrac > inJCI->cutCHF);
            passJetID &= (ChargeMult > inJCI->cutChMult);
//            std::cout << "passJetIDStage 2? " << passJetID << std::endl;
        }
        return passJetID;
    }
    void SetJetVars(PFJetEventPointers inPFJEPs, int iJet) {
        //Set LorentzVector
        P4.SetPxPyPzE(inPFJEPs.JetPx->at(iJet), inPFJEPs.JetPy->at(iJet), inPFJEPs.JetPz->at(iJet), inPFJEPs.JetE->at(iJet));
        // Set BTagging info
        valBTagDisc = inPFJEPs.JetBTag->at(iJet);
        partonFlavor = inPFJEPs.JetPartFlav->at(iJet);
        
        //Set Jet quality cut variables
        NeutHadFrac = inPFJEPs.JetNHF->at(iJet);
        ChargeHadFrac = inPFJEPs.JetCHF->at(iJet);
        NeutEMFrac = inPFJEPs.JetNEF->at(iJet);
        ChargeEMFrac = inPFJEPs.JetCEF->at(iJet);
        NumDaughters = inPFJEPs.JetNDaug->at(iJet);
        ChargeMult      = inPFJEPs.JetCharMult->at(iJet);
    }
    void PrintVals() {
        std::cout << "P4.Pt() " << P4.Pt() << std::endl;
        std::cout << "P4.Eta() " << P4.Eta() << std::endl;
        std::cout << "P4.Phi() " << P4.Phi() << std::endl;
        std::cout << "P4.E() " << P4.E() << std::endl;
        
        std::cout << "valBTagDisc " << valBTagDisc << std::endl;
        std::cout << "partonFlavor " << partonFlavor << std::endl;
        
        std::cout << "NeutHadFrac " << NeutHadFrac << std::endl;        
        std::cout << "ChargeHadFrac " << ChargeHadFrac << std::endl;
        std::cout << "NeutEMFrac " << NeutEMFrac << std::endl;
        std::cout << "ChargeEMFrac " << ChargeEMFrac << std::endl;
        std::cout << "NumDaughters " << NumDaughters << std::endl;
        std::cout << "ChargeMult " << ChargeMult << std::endl;
    }
} PFJet;
inline bool operator<(const PFJet &a, const PFJet &b)
{
    return (a.P4.Pt() < b.P4.Pt());
}
inline bool operator>(const PFJet &a, const PFJet &b)
{
    return (a.P4.Pt() > b.P4.Pt());
}
typedef struct GenJet{
    TLorentzVector  P4;
    int             seedPDGID;
    float           invisE;
    bool            isBJet;
    void ClearVars() {
        P4.SetPxPyPzE(0., 0., 0., 0.);
        seedPDGID = 0;
        invisE = 0.0;
        isBJet = false;
    }
} GenJet;
inline bool operator<(const GenJet &a, const GenJet &b)
{
    return (a.P4.Pt() < b.P4.Pt());
}

inline bool operator>(const GenJet &a, const GenJet &b)
{
    return (a.P4.Pt() > b.P4.Pt());
}

typedef struct GenJetEventPointers {
    std::vector<float> * genJetPx, * genJetPy, * genJetPz, * genJetEn, * genJetInvisE; //* genJetEt, * genJetEta;
    std::vector<bool> * genJetIsGenJet;
    unsigned int numGenJets;
}GenJetEventPointers;

typedef struct GenParticle{
    TLorentzVector P4;
    int            PDGID;
    int            PDGStatus;
    float          Mass;
    void ClearVals() {        
        P4.SetPxPyPzE(-99999., -99999., -99999., -99999.);
        PDGID = 0;
        PDGStatus = -1;
        Mass = 0.;
    }
}GenParticle;

inline bool operator<(const GenParticle &a, const GenParticle &b)
{
    // sort based on PDG status (status 1 takes priority, ratcheting down after that)
    // then among those with same PDGstatus, sort based on PDGID number
    // then among those with same PDGID number -- either anti- or normal particle versions
    // sort based on Pt
    bool LessThan;
    if (a.PDGStatus == b.PDGStatus) {
        if (abs(a.PDGID) == abs(b.PDGID)) {
            LessThan = a.P4.Pt() < b.P4.Pt();
        }
        else {
            LessThan = (abs(a.PDGID) > abs(b.PDGID));
        }
    }
    else {
        LessThan = (a.PDGStatus > b.PDGStatus);
    }
    return LessThan;
}

inline bool operator>(const GenParticle &a, const GenParticle &b)
{
    // sort based on PDG status (status 1 takes priority, ratcheting down after that)
    // then among those with same PDGstatus, sort based on PDGID number
    // then among those with same PDGID number -- either anti- or normal particle versions
    // sort based on Pt
    bool GreaterThan;
    if (a.PDGStatus == b.PDGStatus) {
        if (abs(a.PDGID) == abs(b.PDGID)) {
            GreaterThan = a.P4.Pt() > b.P4.Pt();
        }
        else {
            GreaterThan = (abs(a.PDGID) < abs(b.PDGID));
        }
    }
    else {
        GreaterThan = (a.PDGStatus < b.PDGStatus);
    }
    return GreaterThan;
}

inline int DiLeptonEventType(int Lep0PDGID, int Lep1PDGID) {
    int productPdgId = Lep0PDGID * Lep1PDGID;
    if (productPdgId == -169) {
        return 0;  
    }
    else if (productPdgId == -121) {
        return 1;  
    }
    else if (productPdgId == -143) {
        return 2;
    }
    else {
        std::cout << "funky business in productPdgId " << productPdgId << std::endl;
    } 
    return -1;
}
typedef struct EventLepInfo{
    float EventDiLepMass, EventLepST;
    bool  EventDiLepinZMass; //When true, means the diLep Mass is inside the ZMassWindow (76:106) GeV
    Lepton Lep0, Lep1;
    std::vector<float> DiLeptonTrigSF;
    std::vector<float> DiLeptonIDIsoSF;
    std::vector<float> DiLeptonTotSF;
    /*
    float DiLeptonTrigSF, DiLeptonTrigSFUpErr, DiLeptonTrigSFDownErr;
    float DiLeptonIDIsoSF, DiLeptonIDIsoSFUpErr, DiLeptonIDIsoSFDownErr;
    */
    
    float EventLep0Px, EventLep0Py, EventLep0Pz, EventLep0E;
    float EventLep1Px, EventLep1Py, EventLep1Pz, EventLep1E;
    int   EventLep0PDGID, EventLep1PDGID;
    float EventLep0RelPFIso, EventLep1RelPFIso;
    
    int   EventNIsoElecs_pT20, EventNIsoMuons_pT20, EventNIsoPosits_pT20, EventNIsoMubars_pT20;
    int   EventNIsoElecs_pT10to20, EventNIsoMuons_pT10to20, EventNIsoPosits_pT10to20, EventNIsoMubars_pT10to20;
    int   EventNViableLepPairsPreMassCut;
    int   EventDiLepType;
    bool  doEvent;
    void ELIDefaultVarVals() {        
        EventDiLepMass = 0.; EventDiLepType = -1; EventLepST = 0.;
        Lep0.ClearVars(); Lep1.ClearVars();
        EventLep0Px = Lep0.P4.Px(); EventLep0Py = Lep0.P4.Py(); EventLep0Pz = Lep0.P4.Pz(); EventLep0E = Lep0.P4.E();
        EventLep1Px = Lep1.P4.Px(); EventLep1Py = Lep1.P4.Py(); EventLep1Pz = Lep1.P4.Pz(); EventLep1E = Lep1.P4.E();
        EventLep0PDGID = Lep0.PDGID; EventLep0RelPFIso = Lep0.relPFLepIso; EventLep1PDGID = Lep1.PDGID; EventLep1RelPFIso = Lep1.relPFLepIso;
        
        EventNIsoElecs_pT20 = 0; EventNIsoMuons_pT20 = 0; EventNIsoPosits_pT20 = 0; EventNIsoMubars_pT20 = 0;
        EventNIsoElecs_pT10to20 = 0; EventNIsoMuons_pT10to20 = 0; EventNIsoPosits_pT10to20 = 0; EventNIsoMubars_pT10to20 = 0;
        EventNViableLepPairsPreMassCut = 0;
        doEvent = false;
        EventDiLepinZMass = false;
        DiLeptonTrigSF.clear();
        DiLeptonIDIsoSF.clear();
        DiLeptonTotSF.clear();
    }
    void EventFails() {
        EventDiLepMass = -99999.; EventDiLepType = -1; EventLepST = -1.;
        Lep0.isBadLep(); Lep1.isBadLep();  
        EventLep0Px = Lep0.P4.Px(); EventLep0Py = Lep0.P4.Py(); EventLep0Pz = Lep0.P4.Pz(); EventLep0E = Lep0.P4.E();
        EventLep1Px = Lep1.P4.Px(); EventLep1Py = Lep1.P4.Py(); EventLep1Pz = Lep1.P4.Pz(); EventLep1E = Lep1.P4.E();
        EventLep0PDGID = Lep0.PDGID; EventLep0RelPFIso = Lep0.relPFLepIso; EventLep1PDGID = Lep1.PDGID; EventLep1RelPFIso = Lep1.relPFLepIso;
        EventNIsoElecs_pT20 = -1; EventNIsoMuons_pT20 = -1; EventNIsoPosits_pT20 = -1; EventNIsoMubars_pT20 = -1;
        EventNIsoElecs_pT10to20 = -1; EventNIsoMuons_pT10to20 = -1; EventNIsoPosits_pT10to20 = -1; EventNIsoMubars_pT10to20 = -1;
        EventNViableLepPairsPreMassCut = -1;
        doEvent = false;
        EventDiLepinZMass = false;
        DiLeptonTrigSF.clear();
        DiLeptonIDIsoSF.clear();
        DiLeptonTotSF.clear();
    }
    void EventPasses(int indexLep0, int indexLep1, std::vector<Lepton> * inputLepVec) {
        Lep0 = inputLepVec->at(indexLep0); Lep1 = inputLepVec->at(indexLep1);
        
        EventLep0Px = Lep0.P4.Px(); EventLep0Py = Lep0.P4.Py(); EventLep0Pz = Lep0.P4.Pz(); EventLep0E = Lep0.P4.E();
        EventLep1Px = Lep1.P4.Px(); EventLep1Py = Lep1.P4.Py(); EventLep1Pz = Lep1.P4.Pz(); EventLep1E = Lep1.P4.E();
        EventLep0PDGID = Lep0.PDGID; EventLep0RelPFIso = Lep0.relPFLepIso; EventLep1PDGID = Lep1.PDGID; EventLep1RelPFIso = Lep1.relPFLepIso;        
        EventDiLepType = DiLeptonEventType(Lep0.PDGID, Lep1.PDGID);
        EventDiLepMass = (Lep0.P4 + Lep1.P4).M();
        EventDiLepinZMass = EventDiLepMass > 76 && EventDiLepMass < 106;
        EventLepST = Lep0.P4.Pt() + Lep1.P4.Pt();
        doEvent = true;
    }
    void GrabbedFromTree() {
        Lep0.P4.SetPxPyPzE(EventLep0Px, EventLep0Py, EventLep0Pz, EventLep0E);
        Lep0.PDGID = EventLep0PDGID; Lep0.relPFLepIso = EventLep0RelPFIso;
        Lep1.P4.SetPxPyPzE(EventLep1Px, EventLep1Py, EventLep1Pz, EventLep1E);
        Lep1.PDGID = EventLep1PDGID; Lep1.relPFLepIso = EventLep1RelPFIso;
        EventLepST = Lep0.P4.Pt() + Lep1.P4.Pt();
        EventDiLepinZMass = EventDiLepMass > 76 && EventDiLepMass < 106;
    }
    void ScaleFactorTrigMC(vector<TH2F *> * vecTrigSF) {

        //    float DiLeptonTrigSFCentVal, DiLeptonTrigSFErr;
        int XBin, YBin;
        float LepEtaX = (abs(Lep0.PDGID) < abs(Lep1.PDGID)) ? Lep0.P4.Eta() : Lep1.P4.Eta();
        float LepEtaY = (abs(Lep0.PDGID) > abs(Lep1.PDGID)) ? Lep1.P4.Eta() : Lep0.P4.Eta();
        XBin = vecTrigSF->at(EventDiLepType)->GetXaxis()->FindBin(fabs(LepEtaX));
        YBin = vecTrigSF->at(EventDiLepType)->GetXaxis()->FindBin(fabs(LepEtaY));
        if (XBin > vecTrigSF->at(EventDiLepType)->GetNbinsX()) XBin = vecTrigSF->at(EventDiLepType)->GetNbinsX();
        if (YBin > vecTrigSF->at(EventDiLepType)->GetNbinsY()) YBin = vecTrigSF->at(EventDiLepType)->GetNbinsY();
        
        //    DiLeptonTrigSFCentVal = vecTrigSF->at(EventDiLepType)->GetBinContent(XBin, YBin);
        //    DiLeptonTrigSFErr = vecTrigSF->at(EventDiLepType)->GetBinError(XBin, YBin);
        //    return DiLeptonTrigSFCentVal + Syst * DiLeptonTrigSFErr;
        DiLeptonTrigSF.resize(3);
        DiLeptonTrigSF[0] = vecTrigSF->at(EventDiLepType)->GetBinContent(XBin, YBin);
        DiLeptonTrigSF[1] = vecTrigSF->at(EventDiLepType)->GetBinError(XBin, YBin);
        DiLeptonTrigSF[2] = vecTrigSF->at(EventDiLepType)->GetBinError(XBin, YBin);
    }
    void ScaleFactorIDIsoMC(bool doSignal, TH2F * histElecFastSimSF, TH2F * histMuonFastSimSF) {
        float Lep0FastSimSF, Lep1FastSimSF;
        float Lep0FastSimSFErr = 0, Lep1FastSimSFErr = 0;
        float TotLepFastSimSF, TotLepFastSimSFUpErr = 0, TotLepFastSimSFDownErr = 0;
        float FastSimSystErrBase = 0.02;
        TH2F * Lep0FastSimHist, * Lep1FastSimHist;
        float SF_IDISO[3]           = {0.997, 0.975, 0.986};
        float SF_IDISOSystUp[3]     = {0.00085, 0.006, 0.007};
        float SF_IDISOSystDown[3]   = {0.00085, 0.006, 0.007};        
        
        float baseSF = SF_IDISO[EventDiLepType];
        float baseSFUpErr = 0, baseSFDownErr = 0;
        int XBin, YBin;
        baseSFUpErr = SF_IDISOSystUp[EventDiLepType];
        baseSFDownErr = SF_IDISOSystDown[EventDiLepType];
        
        DiLeptonIDIsoSF.resize(3);
        if (doSignal) {
            Lep0FastSimHist = Lep0.isElec() ? histElecFastSimSF : histMuonFastSimSF;
            XBin = Lep0FastSimHist->GetXaxis()->FindBin(Lep0.P4.Pt());
            YBin = Lep0FastSimHist->GetYaxis()->FindBin(fabs(Lep0.P4.Eta()));
            if (XBin > Lep0FastSimHist->GetNbinsX()) XBin = Lep0FastSimHist->GetNbinsX();
            if (YBin > Lep0FastSimHist->GetNbinsY()) YBin = Lep0FastSimHist->GetNbinsY();
            Lep0FastSimSF = Lep0FastSimHist->GetBinContent(XBin, YBin);
            
            Lep1FastSimHist = Lep1.isElec() ? histElecFastSimSF : histMuonFastSimSF;
            XBin = Lep1FastSimHist->GetXaxis()->FindBin(Lep1.P4.Pt());
            YBin = Lep1FastSimHist->GetYaxis()->FindBin(fabs(Lep1.P4.Eta()));
            if (XBin > Lep1FastSimHist->GetNbinsX()) XBin = Lep1FastSimHist->GetNbinsX();
            if (YBin > Lep1FastSimHist->GetNbinsY()) YBin = Lep1FastSimHist->GetNbinsY();
            Lep1FastSimSF = Lep1FastSimHist->GetBinContent(XBin, YBin);
            
            TotLepFastSimSF = Lep0FastSimSF * Lep1FastSimSF;
            
            Lep0FastSimSFErr = FastSimSystErrBase * Lep0FastSimSF;
            Lep1FastSimSFErr = FastSimSystErrBase * Lep1FastSimSF;
            TotLepFastSimSFUpErr = Lep0FastSimSFErr * Lep1FastSimSFErr;
            TotLepFastSimSFUpErr += Lep0FastSimSF * Lep1FastSimSFErr + Lep1FastSimSF * Lep0FastSimSFErr;
            TotLepFastSimSFDownErr = Lep0FastSimSFErr * Lep1FastSimSFErr;
            TotLepFastSimSFDownErr -= Lep0FastSimSF * Lep1FastSimSFErr + Lep1FastSimSF * Lep0FastSimSFErr;
            TotLepFastSimSFDownErr = fabs(TotLepFastSimSFDownErr);
            
            DiLeptonIDIsoSF[0] = baseSF * TotLepFastSimSF;
            DiLeptonIDIsoSF[1] = TotLepFastSimSFUpErr * baseSFUpErr; 
            DiLeptonIDIsoSF[1] += TotLepFastSimSF * baseSFUpErr + TotLepFastSimSFUpErr * baseSF;
            DiLeptonIDIsoSF[2] = TotLepFastSimSFDownErr * baseSFDownErr;
            DiLeptonIDIsoSF[2] -= TotLepFastSimSF * baseSFDownErr + TotLepFastSimSFDownErr * baseSF;
            DiLeptonIDIsoSF[2] = fabs(DiLeptonIDIsoSF[2]);
        }
        else {        
            DiLeptonIDIsoSF[0] = baseSF;
            DiLeptonIDIsoSF[1] = baseSFUpErr;
            DiLeptonIDIsoSF[2] = baseSFDownErr;
        }
    }
    void SetScaleFactors(bool doSignal, std::vector<TH2F *> * vecTrigSFBkg, std::vector<TH2F *> * vecTrigSFSig, TH2F * histElecFastSimSF, TH2F * histMuonFastSimSF) {
        if (doSignal) {
            ScaleFactorTrigMC(vecTrigSFSig);
        }
        else {
            ScaleFactorTrigMC(vecTrigSFBkg);
        }
        ScaleFactorIDIsoMC(doSignal, histElecFastSimSF, histMuonFastSimSF);
        DiLeptonTotSF.resize(3);
        if (DiLeptonIDIsoSF.size() < 3 || DiLeptonTrigSF.size() < 3) {
            cout << "issue: sizes are less than 3! check that they got set " << endl;
            cout << "DiLeptonIDIsoSF.size()  " << DiLeptonIDIsoSF.size() << endl;
            cout << "DiLeptonTrigSF.size()  " << DiLeptonTrigSF.size() << endl;
        }
        else {
            /*
            cout << "DiLeptonTrigSF[0] " << DiLeptonTrigSF[0] << endl;
            cout << "DiLeptonTrigSF[1] " << DiLeptonTrigSF[1] << endl;
            cout << "DiLeptonTrigSF[2] " << DiLeptonTrigSF[2] << endl;
            cout << "DiLeptonIDIsoSF[0] " << DiLeptonIDIsoSF[0] << endl;
            cout << "DiLeptonIDIsoSF[1] " << DiLeptonIDIsoSF[1] << endl;
            cout << "DiLeptonIDIsoSF[2] " << DiLeptonIDIsoSF[2] << endl;
             */
            DiLeptonTotSF[0] = DiLeptonTrigSF[0] * DiLeptonIDIsoSF[0];
            DiLeptonTotSF[1] = TMath::Sqrt(DiLeptonTrigSF[1] * DiLeptonTrigSF[1] + DiLeptonIDIsoSF[1] * DiLeptonIDIsoSF[1]);
            DiLeptonTotSF[2] = TMath::Sqrt(DiLeptonTrigSF[2] * DiLeptonTrigSF[2] + DiLeptonIDIsoSF[2] * DiLeptonIDIsoSF[2]);
        }
    }
    void SetScaleFactorsFailedEvent() {
        DiLeptonTrigSF.resize(3);
        DiLeptonTrigSF[0] = 0.0;
        DiLeptonTrigSF[1] = 0.0;
        DiLeptonTrigSF[1] = 0.0;
        DiLeptonIDIsoSF.resize(3);
        DiLeptonIDIsoSF[0] = 0.0;
        DiLeptonIDIsoSF[1] = 0.0;
        DiLeptonIDIsoSF[2] = 0.0;
        DiLeptonTotSF.resize(3);
        DiLeptonTotSF[0] = 0.0;
        DiLeptonTotSF[1] = 0.0;
        DiLeptonTotSF[2] = 0.0;
    }
    float GetSF(int Syst = 0) {
        if (DiLeptonTotSF.size() < 3) cout << "issue with DiLeptonTotSF.size(): " << DiLeptonTotSF.size() << endl;
        float outSF = DiLeptonTotSF[0];
//        cout << "DiLeptonTotSF[0] " << DiLeptonTotSF[0] << endl;
        if (Syst > 0) {
            outSF += DiLeptonTotSF[1];
//            cout << "DiLeptonTotSF[1] " << DiLeptonTotSF[1] << endl;
        }
        else if (Syst < 0) {
            outSF -= DiLeptonTotSF[2];
//            cout << "DiLeptonTotSF[2] " << DiLeptonTotSF[1] << endl;
        }
        return outSF;
    }
    /*
    ELISetValsInput(int indexLep0, int indexLep1, std::vector<Lepton> * inputLepVec) {
        Lep0 = inputLepVec->at(indexLep0); Lep1 = inputLepVec->at(indexLep1);
        EventLep0Px = Lep0.P4.Px(); EventLep0Py = Lep0.P4.Py(); EventLep0Pz = Lep0.P4.Pz(); EventLep0E = Lep0.P4.E();
        EventLep1Px = Lep1.P4.Px(); EventLep1Py = Lep1.P4.Py(); EventLep1Pz = Lep1.P4.Pz(); EventLep1E = Lep1.P4.E();
        EventLep0PDGID = Lep0.PDGID; EventLep0RelPFIso = Lep0.relPFLepIso; EventLep1PDGID = Lep1.PDGID; EventLep1RelPFIso = Lep1.relPFLepIso;        
        EventDiLepType = DiLeptonEventType(Lep0.PDGID, Lep1.PDGID);
        EventDiLepMass = (Lep0.P4 + Lep1.P4).M();
        doEvent = true;
    }
    */
} EventLepInfo;

typedef struct EventJetInfo{
    float EventHT, EventJetST;
    int   EventNJets, EventNBtagJets;
    PFJet Jet0, Jet1;
    PFJet BtagJet0, BtagJet1;
    int   EventBtagJet0Index, EventBtagJet1Index;
    float EventJet0Px, EventJet0Py, EventJet0Pz, EventJet0E;
    float EventJet1Px, EventJet1Py, EventJet1Pz, EventJet1E;
    float EventBtagJet0Px, EventBtagJet0Py, EventBtagJet0Pz, EventBtagJet0E;
    float EventBtagJet1Px, EventBtagJet1Py, EventBtagJet1Pz, EventBtagJet1E;
    bool  EventJet0isGenMatched, EventJet1isGenMatched;
    float EventJet0DeltaEnRecoGen, EventJet1DeltaEnRecoGen;
    
    void SetNonBTagInfo(EventJetInfo * inEJI) {
        EventHT = inEJI->EventHT;
        EventJetST = inEJI->EventJetST;
        EventNJets = inEJI->EventNJets;
        Jet0 = inEJI->Jet0;
        Jet1 = inEJI->Jet1;
        EventJet0Px = inEJI->EventJet0Px;
        EventJet0Py = inEJI->EventJet0Py;
        EventJet0Pz = inEJI->EventJet0Pz;
        EventJet0E =  inEJI->EventJet0E;
        EventJet1Px = inEJI->EventJet1Px;
        EventJet1Py = inEJI->EventJet1Py;
        EventJet1Pz = inEJI->EventJet1Pz;
        EventJet1E =  inEJI->EventJet1E;
        EventJet0isGenMatched = inEJI->EventJet0isGenMatched;
        EventJet1isGenMatched = inEJI->EventJet0isGenMatched;
        EventJet0DeltaEnRecoGen = inEJI->EventJet0DeltaEnRecoGen;
        EventJet1DeltaEnRecoGen = inEJI->EventJet1DeltaEnRecoGen;
        
        BtagJet0.isBadJet(); BtagJet1.isBadJet();
        if (EventNBtagJets > 0) {
            BtagJet0.P4.SetPxPyPzE(EventBtagJet0Px, EventBtagJet0Py, EventBtagJet0Pz, EventBtagJet0E);
            if (EventNBtagJets > 1) {
                BtagJet1.P4.SetPxPyPzE(EventBtagJet1Px, EventBtagJet1Py, EventBtagJet1Pz, EventBtagJet1E);
            }
        }
    }
    void EJIDefaultVarVals() {
        EventHT = 0.; EventJetST = 0.; EventNJets = 0; EventNBtagJets = 0;
        Jet0.ClearVars(); Jet1.ClearVars();
        BtagJet0.ClearVars(); BtagJet1.ClearVars();        
        EventBtagJet0Index = -1; EventBtagJet1Index = -1;
        EventJet0Px = Jet0.P4.Px(); EventJet0Py = Jet0.P4.Py(); EventJet0Pz = Jet0.P4.Pz(); EventJet0E = Jet0.P4.E();
        EventJet1Px = Jet1.P4.Px(); EventJet1Py = Jet1.P4.Py(); EventJet1Pz = Jet1.P4.Pz(); EventJet1E = Jet1.P4.E();
        EventBtagJet0Px = BtagJet0.P4.Px(); EventBtagJet0Py = BtagJet0.P4.Py(); EventBtagJet0Pz = BtagJet0.P4.Pz(); EventBtagJet0E = BtagJet0.P4.E();
        EventBtagJet1Px = BtagJet1.P4.Px(); EventBtagJet1Py = BtagJet1.P4.Py(); EventBtagJet1Pz = BtagJet1.P4.Pz(); EventBtagJet1E = BtagJet1.P4.E();
        EventJet0isGenMatched = false; EventJet1isGenMatched = false;
        EventJet0DeltaEnRecoGen = 0.0; EventJet1DeltaEnRecoGen = 0.0;
    }
    void EJISetValsInput(float inputHT, int inputNJets, int indexJet0, int indexJet1, int inputNBtagJets, int inputBtagJet0Index, int inputBtagJet0SubIndex, int inputBtagJet1Index, int inputBtagJet1SubIndex, std::vector<PFJet> * inputPFJetVec) {
        EventHT = inputHT; EventJetST = 0.;
        EventNJets = inputNJets; EventNBtagJets = inputNBtagJets;
        EventBtagJet0Index = inputBtagJet0SubIndex; EventBtagJet1Index = inputBtagJet1SubIndex;
        EventJet0isGenMatched = false; EventJet1isGenMatched = false;
        EventJet0DeltaEnRecoGen = 0.0; EventJet1DeltaEnRecoGen = 0.0;
        
        //Default behavior, assume all jets are bad -- i.e. no good jets
        Jet0.isBadJet(); Jet1.isBadJet();
        BtagJet0.isBadJet(); BtagJet1.isBadJet();
        if (inputNJets > 0) {
            Jet0 = inputPFJetVec->at(indexJet0);
            EventJetST += Jet0.P4.Pt();
            if (inputNBtagJets > 0) {
                BtagJet0 = inputPFJetVec->at(inputBtagJet0Index);
            }
            if (inputNJets > 1) {
                Jet1 = inputPFJetVec->at(indexJet1);
                EventJetST += Jet1.P4.Pt();
                if (inputNBtagJets > 1) {
                    BtagJet1 = inputPFJetVec->at(inputBtagJet1Index);
                }
            }
        }
        else if (inputNJets < 0) {
            std::cout << "something weird with number of input jets " << inputNJets << std::endl;
        }
        EventJet0Px = Jet0.P4.Px(); EventJet0Py = Jet0.P4.Py(); EventJet0Pz = Jet0.P4.Pz(); EventJet0E = Jet0.P4.E();
        EventJet1Px = Jet1.P4.Px(); EventJet1Py = Jet1.P4.Py(); EventJet1Pz = Jet1.P4.Pz(); EventJet1E = Jet1.P4.E();
        EventBtagJet0Px = BtagJet0.P4.Px(); EventBtagJet0Py = BtagJet0.P4.Py(); EventBtagJet0Pz = BtagJet0.P4.Pz(); EventBtagJet0E = BtagJet0.P4.E();
        EventBtagJet1Px = BtagJet1.P4.Px(); EventBtagJet1Py = BtagJet1.P4.Py(); EventBtagJet1Pz = BtagJet1.P4.Pz(); EventBtagJet1E = BtagJet1.P4.E();
    }
    void EJISetGenMatchValsInput(int inputNJets, int indexJet0, int indexJet1, std::vector<PFJet> * inputPFJetVec) {
        if (inputNJets > 0) {
            EventJet0isGenMatched = inputPFJetVec->at(indexJet0).isGenJetMatched;
            EventJet0DeltaEnRecoGen = inputPFJetVec->at(indexJet0).dEnRecoGen; 
            if (inputNJets > 1) {
                EventJet1isGenMatched = inputPFJetVec->at(indexJet1).isGenJetMatched;
                EventJet1DeltaEnRecoGen = inputPFJetVec->at(indexJet1).dEnRecoGen;
            }
        }
        else if (inputNJets < 0) {
            std::cout << "something weird with number of input jets " << inputNJets << std::endl;
        }        
    }    
    void GrabbedFromTree() {
            //Default behavior, assume all jets are bad -- i.e. no good jets
            Jet0.isBadJet(); Jet1.isBadJet();
            BtagJet0.isBadJet(); BtagJet1.isBadJet();
            if (EventNJets > 0) {
                Jet0.P4.SetPxPyPzE(EventJet0Px, EventJet0Py, EventJet0Pz, EventJet0E);
                EventJetST += Jet0.P4.Pt();
                if (EventNBtagJets > 0) {
                    BtagJet0.P4.SetPxPyPzE(EventBtagJet0Px, EventBtagJet0Py, EventBtagJet0Pz, EventBtagJet0E);
                }
                if (EventNJets > 1) {                    
                    Jet1.P4.SetPxPyPzE(EventJet1Px, EventJet1Py, EventJet1Pz, EventJet1E);
                    EventJetST += Jet1.P4.Pt();
                    if (EventNBtagJets > 1) {
                        BtagJet1.P4.SetPxPyPzE(EventBtagJet1Px, EventBtagJet1Py, EventBtagJet1Pz, EventBtagJet1E);
                    }
                }
            }
            else if (EventNJets < 0) {
                std::cout << "something weird with number of input jets " << EventNJets << std::endl;
            }
    }
    /*
     NumJetsEqZero() {
     EventHT = 0.; EventNJets = 0; EventNBtagJets = 0;
     Jet0.isBadJet(); Jet1.isBadJet();
     BtagJet0.isBadJet(); BtagJet1.isBadJet();
     EventBtagJet0Index = -1; EventBtagJet1Index = -1;
     }
     NumJetsEqOne(float inputHT, int indexJet0, std::vector<PFJet> * inputPFJetVec) {
     EventHT = inputHT; EventNJets = 1;
     Jet0 = inputPFJetVec->at(indexJet0);
     Jet1.isBadJet(); BtagJet1.isBadJet();
     EventBtagJet1Index = -1;
     }
     NumJetsGEqTwo(float inputHT, int inputNJets, int indexJet0, int indexJet1, std::vector<PFJet> * inputPFJetVec) {
     EventHT = inputHT; EventNJets = inputNJets;
     Jet0 = inputPFJetVec->at(indexJet0); Jet1 = inputPFJetVec->at(indexJet1);
     }
     NumBJetsEqOne(int inputBtagJet0Index, std::vector<PFJet> * inputPFJetVec) {
     EventNBtagJets = 1; 
     EventBtagJet0Index = inputBtagJet0Index;
     BtagJet0 = inputPFJetVec->at(inputBtagJet0Index);
     }
     NumBJetsGEqTwo(int inputNBtagJets. int inputBtagJet0Index, int inputBtagJet1Index, std::vector<PFJet> * inputPFJetVec) {
     EventNBtagJets = inputNBtagJets;
     EventBtagJet0Index = inputBtagJet0Index; EventBtagJet1Index = inputBtagJet1Index;
     BtagJet0 = inputPFJetVec->at(inputBtagJet0Index); BtagJet1 = inputPFJetVec->at(inputBtagJet1Index);
     }
     */
} EventJetInfo;

typedef struct EventMETInfo {
    TVector3 P3;
    float EventMET, EventMETPhi, EventMETX, EventMETY, EventMETSig;
    float EventMET_preCorr, EventMETPhi_preCorr, EventMETX_preCorr, EventMETY_preCorr, EventMETSig_preCorr;
    float EventMT2ll, EventMT2lb;
    float EventDeltaPhiMT2lb_JetsUsed, EventDeltaPhiMT2lb_BLepsUsed;
    //    std::vector<TLorentzVector> EventVecLepMT2lb(2), EventVecJetMT2lb(2), EventVecBLepsMT2lb(2);
    std::vector<TLorentzVector> EventVecBLepsMT2lb;
    float EventMETdivMeff;
    int caseMT2lb;
    void EMIDefaultVarVals() {   
        EventMET = 0.; EventMETPhi = 0.; EventMETX = 0.; EventMETY = 0.; EventMETSig = 0.;
        EventMET_preCorr = 0.; EventMETPhi_preCorr = 0.; EventMETX_preCorr = 0.; EventMETY_preCorr = 0.; EventMETSig_preCorr = 0.;
        EventMT2ll = 0.; EventMT2lb = 0.;
        EventDeltaPhiMT2lb_JetsUsed = 0.; EventDeltaPhiMT2lb_BLepsUsed= 0.;
        EventVecBLepsMT2lb.clear();
        EventVecBLepsMT2lb.resize(2);
        EventMETdivMeff = 0.;
        caseMT2lb = 0;
    }
    
    void PrintVals() {
        std::cout << "EventMET " << EventMET << std::endl;
        std::cout << "EventMETPhi " << EventMETPhi << std::endl;
        std::cout << "EventMETX " << EventMETX << std::endl;
        std::cout << "EventMETY " << EventMETY << std::endl;
        
        std::cout << "EventMET_preCorr " << EventMET_preCorr << std::endl;
        std::cout << "EventMETPhi_preCorr " << EventMETPhi_preCorr << std::endl;
        std::cout << "EventMETX_preCorr " << EventMETX_preCorr << std::endl;
        std::cout << "EventMETY_preCorr " << EventMETY_preCorr << std::endl;
        std::cout << "EventMETdivMeff " << EventMETdivMeff << std::endl;
    }
    void GrabbedFromTree() {
        P3.SetPtEtaPhi(EventMET, 0., EventMETPhi);
        EventMETX = P3.Px(); EventMETY = P3.Py();
    }
    void SetVarsfromVec() {
        EventMET = P3.Pt();
        EventMETPhi = P3.Phi();
        EventMETX = P3.Px();
        EventMETY = P3.Py();
    }
    void SetVarsfromXY() {
        P3.SetXYZ(EventMETX, EventMETY, 0.);
        EventMET = P3.Pt();
        EventMET = P3.Phi();
    }
} EventMETInfo;

typedef struct EventSpecialMT2Info {
    float EventMT2ll_ShiftUp; //, MT2ll_ShiftDown; // TTBar MT2ll smearing
    float EventMT2ll_UncESUp, EventMT2ll_UncESDown;
    float EventMT2lb_UncESUp, EventMT2lb_UncESDown;    
    float EventMT2llSmearFactor;
    int   EventMT2llSmearFactorBin;
    void SetVars(EventMETInfo * inEMI, EventJetInfo * inEJI, TH1F * MT2llMeanSmear, TH2F * MT2llUncEnUpDelta2D, TH2F * MT2llUncEnDownDelta2D, TH2F * MT2lbUncEnUpDelta2D, TH2F * MT2lbUncEnDownDelta2D, std::vector<TH1F *> * vecOneDeeMT2llUncEnUp, std::vector<TH1F *> * vecOneDeeMT2llUncEnDown, std::vector<TH1F *> * vecOneDeeMT2lbUncEnUp, std::vector<TH1F *> * vecOneDeeMT2lbUncEnDown) {  

        TRandom rand;
        EventMT2llSmearFactorBin = MT2llMeanSmear->FindBin(inEMI->EventMT2ll);
        if (EventMT2llSmearFactorBin > 20) EventMT2llSmearFactorBin = 21;        
        EventMT2llSmearFactor = MT2llMeanSmear->GetBinContent(EventMT2llSmearFactorBin);
        EventMT2ll_ShiftUp = inEMI->EventMT2ll + rand.Gaus(0, EventMT2llSmearFactor);
        if (EventMT2ll_ShiftUp < 0) EventMT2ll_ShiftUp = 0;
        
        EventMT2ll_UncESUp = inEMI->EventMT2ll - DeltaMT2UncEn(vecOneDeeMT2llUncEnUp, MT2llUncEnUpDelta2D, inEMI->EventMT2ll);
        EventMT2ll_UncESDown = inEMI->EventMT2ll - DeltaMT2UncEn(vecOneDeeMT2llUncEnDown, MT2llUncEnDownDelta2D, inEMI->EventMT2ll);
        if (EventMT2ll_UncESUp < 0) EventMT2ll_UncESUp = 0;
        if (EventMT2ll_UncESDown < 0) EventMT2ll_UncESDown = 0;
        
        if (inEJI->EventNJets > 1) {
            EventMT2lb_UncESUp = inEMI->EventMT2lb - DeltaMT2UncEn(vecOneDeeMT2lbUncEnUp, MT2lbUncEnUpDelta2D, inEMI->EventMT2lb);
            EventMT2lb_UncESDown = inEMI->EventMT2lb - DeltaMT2UncEn(vecOneDeeMT2lbUncEnDown, MT2lbUncEnDownDelta2D, inEMI->EventMT2lb);
            if (EventMT2lb_UncESUp < 0) EventMT2lb_UncESUp = 0;
            if (EventMT2lb_UncESDown < 0) EventMT2lb_UncESDown = 0;
        }
        else {
            EventMT2lb_UncESUp = -99.;
            EventMT2lb_UncESDown = -99.;
        }
    }
} EventSpecialMT2Info;

typedef struct EventDiStructureInfo{
    float diLepInvMass, diLepPt, diLepEta, diLepPhi;
    float diJetInvMass, diJetPt, diJetEta, diJetPhi;
    float diBJetInvMass, diBJetPt, diBJetEta, diBJetPhi;
    float ELepEJet;
    float DPhiLep0Lep1;
    float DPhiLep0MET, DPhiLep1MET;
    float DPhiLep0MET_PreCorr, DPhiLep1MET_PreCorr;
    float DPhiZMET, DPhiZMET_PreCorr;
    float DPhiLep0Jet0, DPhiLep0Jet1;
    float DPhiBLep0BLep1;
    float DPhiLep0BJet0, DPhiLep0BJet1;
    float DPhiJet0BJet0, DPhiJet0BJet1;
    float DPhiJet1BJet0, DPhiJet1BJet1;
    void EDSIDefaultVarVals() {
        diLepInvMass = 0.; diLepPt = 0.;
        diLepEta = 0.; diLepPhi = 0.;
        diJetInvMass = 0.; diJetPt = 0.;
        diJetEta = 0.; diJetPhi = 0.;
        diBJetInvMass = 0.; diBJetPt = 0.;
        diBJetEta = 0.; diBJetPhi = 0.;
        ELepEJet = 0.;
        DPhiLep0Lep1 = -99.;
        DPhiLep0MET = -99.; DPhiLep1MET = -99.;
        DPhiLep0MET_PreCorr = -99.; DPhiLep1MET_PreCorr = -99.;
        DPhiZMET = -99.; DPhiZMET_PreCorr = -99.;
        DPhiLep0Jet0 = -99.; DPhiLep0Jet1 = -99.;
        DPhiBLep0BLep1 = -99.;
        DPhiLep0BJet0 = -99.; DPhiLep0BJet1 = -99.;
        DPhiJet0BJet0 = -99.; DPhiJet0BJet1 = -99.;
        DPhiJet1BJet0 = -99.; DPhiJet1BJet1 = -99.;
    }
    void SetVars(EventLepInfo * inELI, EventJetInfo * inEJI, EventMETInfo * inEMI) {
//        cout << "test a" << endl;
        diLepInvMass = inELI->EventDiLepMass;
//        cout << "test b" << endl;
        diLepPt  = (inELI->Lep0.P4 + inELI->Lep1.P4).Pt();
//        cout << "test c" << endl;
        diLepEta = (inELI->Lep0.P4 + inELI->Lep1.P4).Eta();
//        cout << "test cd" << endl;
        diLepPhi = (inELI->Lep0.P4 + inELI->Lep1.P4).Phi();
//        cout << "test e" << endl;
        DPhiLep0Lep1        = dPhi(inELI->Lep0.P4.Phi(), inELI->Lep1.P4.Phi());
//        cout << "test e1" << endl;
        DPhiLep0MET         = dPhi((float) inELI->Lep0.P4.Phi(), inEMI->EventMETPhi);
//        cout << "test e2 " << DPhiLep0MET << endl;
//        cout << "additional test " << inELI->Lep0.P4.Phi() << endl;
//        cout << "additional test " << inEMI->EventMETPhi_preCorr << endl;
        DPhiLep0MET_PreCorr = dPhi((float) inELI->Lep0.P4.Phi(), inEMI->EventMETPhi_preCorr);
//        cout << "test e3" << endl;
        DPhiLep1MET         = dPhi((float) inELI->Lep1.P4.Phi(), inEMI->EventMETPhi);
//        cout << "test e34" << endl;
        DPhiLep1MET_PreCorr = dPhi((float) inELI->Lep1.P4.Phi(), inEMI->EventMETPhi_preCorr);        
//        cout << "test e5" << endl;
        DPhiZMET            = dPhi(diLepPhi, inEMI->EventMETPhi);
//        cout << "test e6" << endl;
        DPhiZMET_PreCorr    = dPhi(diLepPhi, inEMI->EventMETPhi_preCorr);
//        cout << "test f" << endl;
        if (inEJI->EventNJets > 0) {
            DPhiLep0Jet0 = dPhi(inELI->Lep0.P4.Phi(), inEJI->Jet0.P4.Phi());
            if (inEJI->EventNJets > 1) {
                DPhiLep0Jet1 = dPhi(inELI->Lep0.P4.Phi(), inEJI->Jet1.P4.Phi());
                ELepEJet = inELI->Lep0.P4.E() + inELI->Lep1.P4.E() - inEJI->Jet0.P4.E() - inEJI->Jet1.P4.E();
                DPhiBLep0BLep1 = dPhi(inEMI->EventVecBLepsMT2lb[0].Phi(), inEMI->EventVecBLepsMT2lb[1].Phi());
                diJetPt      = (inEJI->Jet0.P4 + inEJI->Jet1.P4).Pt();
                diJetEta     = (inEJI->Jet0.P4 + inEJI->Jet1.P4).Eta();
                diJetPhi     = (inEJI->Jet0.P4 + inEJI->Jet1.P4).Phi();
                diJetInvMass = (inEJI->Jet0.P4 + inEJI->Jet1.P4).M();
            }
//            cout << "test g" << endl;
            if (inEJI->EventNBtagJets > 0) {
                DPhiLep0BJet0 = dPhi(inELI->Lep0.P4.Phi(), inEJI->BtagJet0.P4.Phi());
                DPhiJet0BJet0 = dPhi(inEJI->Jet0.P4.Phi(), inEJI->BtagJet0.P4.Phi());
                if (inEJI->EventNBtagJets > 1) {
                    DPhiLep0BJet1 = dPhi(inELI->Lep0.P4.Phi(), inEJI->BtagJet1.P4.Phi());
                    DPhiJet0BJet1 = dPhi(inEJI->Jet0.P4.Phi(), inEJI->BtagJet1.P4.Phi());
                    DPhiJet1BJet1 = dPhi(inEJI->Jet1.P4.Phi(), inEJI->BtagJet1.P4.Phi());
                    diBJetPt      = (inEJI->BtagJet0.P4 + inEJI->BtagJet1.P4).Pt();
                    diBJetEta     = (inEJI->BtagJet0.P4 + inEJI->BtagJet1.P4).Eta();
                    diBJetPhi     = (inEJI->BtagJet0.P4 + inEJI->BtagJet1.P4).Phi();
                    diBJetInvMass = (inEJI->BtagJet0.P4 + inEJI->BtagJet1.P4).M();
                }
            }
        }
    }
} EventDiStructureInfo;


typedef struct BasicEventInfo {
    bool  doData, doReReco, doPhiCorr, isSignal;
    bool  doPURW, doPURWOviToDESY, doHackPURW;
    bool  blindData;
    int   nVtx, nVtxTrue;
    bool  passTrigDoubleMu, passTrigDoubleEl, passTrigElMu;
    float weight, preNVtxRWWeight;
    int   runNumber, lumiBlock, eventNumber;
    
    bool hasMETInfo, hasStopInfo, hasTopInfo;
    
    void  SetVars(BasicEventInfo * inBEI) {
        doData              = inBEI->doData;
        doReReco            = inBEI->doReReco; 
        doPhiCorr           = inBEI->doPhiCorr;
        isSignal            = inBEI->isSignal;
        doPURW              = inBEI->doPURW;
        doPURWOviToDESY     = inBEI->doPURWOviToDESY;
        doHackPURW          = inBEI->doHackPURW;
        nVtx                = inBEI->nVtx;
        nVtxTrue            = inBEI->nVtxTrue;
        passTrigDoubleEl    = inBEI->passTrigDoubleEl;
        passTrigDoubleMu    = inBEI->passTrigDoubleMu;
        passTrigElMu        = inBEI->passTrigElMu;
        weight              = inBEI->weight;
        preNVtxRWWeight     = inBEI->preNVtxRWWeight;
        runNumber           = inBEI->runNumber;
        lumiBlock           = inBEI->lumiBlock;
        eventNumber         = inBEI->eventNumber;
        hasMETInfo          = inBEI->hasMETInfo;
        hasStopInfo         = inBEI->hasStopInfo;
        hasTopInfo          = inBEI->hasTopInfo;
    }
} BasicEventInfo;

typedef struct EventGenParticleInfo{
    GenParticle genMET;
    GenParticle genTop0, genTop1;
    GenParticle genStop0, genStop1;
    GenParticle genChiZero0, genChiZero1;
    GenParticle genChargino0, genChargino1;
    float       genTopPt, genAntiTopPt;
    float       StopWeight, StopWeightErr;
    float       StopWeightPlusErr, StopWeightMinusErr;
    float       weight_GenTopReweight;
    bool hasMETInfo, hasTopInfo, hasStopInfo;
    
    void EGPIDefaultVarVals() {
        genMET.ClearVals();
        genTop0.ClearVals(); genTop1.ClearVals();
        genStop0.ClearVals(); genStop1.ClearVals();
        genChiZero0.ClearVals(); genChiZero1.ClearVals();
        genChargino0.ClearVals(); genChargino1.ClearVals();
        genTopPt = 0.; genAntiTopPt = 0.;
        StopWeight = 0.; StopWeightErr = 0.;
        StopWeightPlusErr = 0.; StopWeightMinusErr = 0.;
        weight_GenTopReweight = 0.;
        hasMETInfo = false; hasTopInfo = false; hasStopInfo = false;

    }
}EventGenParticleInfo;




