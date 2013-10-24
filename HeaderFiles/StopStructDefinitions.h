
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

typedef struct PFJet{
    TLorentzVector  P4;
    float           valBTagDisc;
    int             partonFlavor;
    bool            isBJet;
    bool            isGenJetMatched;
    float           dEnRecoGen;
    void ClearVars() {
        P4.SetPxPyPzE(0., 0., 0., 0.);
        valBTagDisc = 0.0;
        partonFlavor = -999999;
        isBJet = false;
        isGenJetMatched = false;
        dEnRecoGen = 0.0;
    }
    void isBadJet() {
        P4.SetPxPyPzE(-99999., -99999., -99999., -99999.);
        valBTagDisc = -1;
        partonFlavor = -999999;
        isBJet = false;
        isGenJetMatched = false;
        dEnRecoGen = 0.0;
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
typedef struct PFJetEventPointers {
    std::vector<float> * JetPx, * JetPy, * JetPz, * JetE, * JetNHF, * JetNEF, * JetCHF, * JetCEF, * JetBTag;
    std::vector<int> * JetNDaug, * JetCharMult, * JetPartFlav;
    unsigned int numPFJets;
//    std::vector<float> * genJetPx, * genJetPy, * genJetPz, * genJetEn, * genJetInvisE; //* genJetEt, * genJetEta;
//    std::vector<bool> * genJetIsGenJet;
}PFJetEventPointers;


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
        diLepInvMass = inELI->EventDiLepMass;
        diLepPt  = (inELI->Lep0.P4 + inELI->Lep1.P4).Pt();
        diLepEta = (inELI->Lep0.P4 + inELI->Lep1.P4).Eta();
        diLepPhi = (inELI->Lep0.P4 + inELI->Lep1.P4).Phi();
        DPhiLep0Lep1        = dPhi(inELI->Lep0.P4.Phi(), inELI->Lep1.P4.Phi());
        DPhiLep0MET         = dPhi((float) inELI->Lep0.P4.Phi(), inEMI->EventMETPhi);
        DPhiLep0MET_PreCorr = dPhi((float) inELI->Lep0.P4.Phi(), inEMI->EventMETPhi_preCorr);
        DPhiLep1MET         = dPhi((float) inELI->Lep1.P4.Phi(), inEMI->EventMETPhi);
        DPhiLep1MET_PreCorr = dPhi((float) inELI->Lep1.P4.Phi(), inEMI->EventMETPhi_preCorr);        
        DPhiZMET            = dPhi(diLepPhi, inEMI->EventMETPhi);
        DPhiZMET_PreCorr    = dPhi(diLepPhi, inEMI->EventMETPhi_preCorr);
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




