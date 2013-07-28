import FWCore.ParameterSet.Config as cms

process = cms.Process("Tree")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')



process.demo = cms.EDAnalyzer('SUSYSkimToTree',
 
  histosFileName = cms.untracked.string("Tree.root"),
#    histosFileName = cms.untracked.string("#outputfile#"),
 
  trigTag = cms.untracked.InputTag('TriggerResults::HLT'),
  muonTag = cms.untracked.InputTag('boostedMuons'),
  pfmuonTag = cms.untracked.InputTag('selectedPatMuonsPF'),
  jetTag=cms.untracked.InputTag('slimPatJetsTriggerMatch'),
  jetPFTag=cms.untracked.InputTag('slimPatJetsTriggerMatchNoPU'),
#  jetPF2Tag=cms.untracked.InputTag('slimPatJetsTriggerMatchPFlow'),
  genmetTag=cms.untracked.InputTag('genMetTrue'),
  metPFTag=cms.untracked.InputTag('pfMet'),
  metPFTypeITag=cms.untracked.InputTag('pfType1CorrectedMet'),
  metChargedNeutralPFNoPUTag=cms.untracked.InputTag('chargedMetProducer'),
  metChargedPFNoPUTag=cms.untracked.InputTag('trackMetProducer'),
  metTCTag=cms.untracked.InputTag('tcMet'),
  PVTag=cms.untracked.InputTag('goodPrimaryVertices'),
  electronTag = cms.untracked.InputTag('boostedElectrons'),
#  pfelectronTag = cms.untracked.InputTag('selectedPatElectronsPF'),
#  pftauTag = cms.untracked.InputTag('cleanPatTausTriggerMatchPFlow'),
  tauTag = cms.untracked.InputTag('cleanPatTausTriggerMatch'),
    singleMuDataPaths = cms.vstring(
#         "1-163261:HLT_Mu15_v*",
#         "163262-165099:HLT_Mu24_v*",
#         "165102-173235:HLT_Mu30_v*",
#         "173236-175972:HLT_Mu40_v*",
#         "175973-180252:HLT_Mu40_eta2p1_v*",
#         "163262-170901:HLT_IsoMu17_v*",
#         "171050-175910:HLT_IsoMu20_v*",
#         "175911-175921:HLT_IsoMu24_v*",
#         "175922-180252:HLT_IsoMu24_eta2p1_v*",
# #       end of 2011 Data
# 	"190456-999999:HLT_IsoMu24_eta2p1_v*",

        # end of 2011 Data
        "190456-999999:HLT_Mu8_v*",
        "190456-999999:HLT_Mu17_v*"

    ),
    doubleMuDataPaths = cms.vstring(
        "1-165208:HLT_DoubleMu7_v*",
        "165364-178419:HLT_Mu13_Mu8_v*",
        "178420-180252:HLT_Mu17_Mu8_v*",
        "178420-180252:HLT_Mu17_TkMu8_v*",
#       end of 2011 Data
	"190456-999999:HLT_Mu17_Mu8_v*",
	"190456-999999:HLT_Mu17_TkMu8_v*",
    ),
    doubleElDataPaths = cms.vstring(
        "1-170901:HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*",
        "171050-180252:HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
#       end of 2011 Data
	"190456-999999:HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
    ),
    muEGDataPaths     = cms.vstring(
        "1-175972:HLT_Mu17_Ele8_CaloIdL_v*",
        "175973-180252:HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*",
        "1-167913:HLT_Mu8_Ele17_CaloIdL_v*",
        "167914-180252:HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*",
#       end of 2011 Data
        "190456-999999:HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
	"190456-999999:HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
    ),
    singleElDataPaths = cms.vstring(
#         "1-164237:HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
#         "165085-166967:HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
#         "166968-170901:HLT_Ele52_CaloIdVT_TrkIdT_v*",
#         "170902-178419:HLT_Ele65_CaloIdVT_TrkIdT_v*",
#         "178420-180252:HLT_Ele80_CaloIdVT_TrkIdT_v*",
# #       end of 2011 Data
# 	"190456-999999:HLT_Ele27_WP80_v*",

#       end of 2011 Data
        "190456-999999:HLT_Ele8_CaloIdT_TrkIdVL_v*",
        "190456-999999:HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*",
        "190456-999999:HLT_Ele8_CaloIdL_CaloIsoVL_v*",
        "190456-999999:HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
        "190456-999999:HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*",
        "190456-999999:HLT_Ele17_CaloIdL_CaloIsoVL_v*",
        "190456-999999:HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
        "190456-999999:HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*"
    ),

    singleMuMCPaths   = cms.vstring("*"),
    singleElMCPaths   = cms.vstring("*"),
    doubleMuMCPaths   = cms.vstring("1-999999:HLT_Mu17_Mu8_v*",
				    "1-999999:HLT_Mu17_TkMu8_v*"),
    doubleElMCPaths   = cms.vstring("1-999999:HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"),
    muEGMCPaths       = cms.vstring("1-999999:HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*",
				    "1-999999:HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*",
					#51X
				    "1-999999:HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
				    "1-999999:HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"
					#52X
)

   )



process.p = cms.Path(
process.demo
)


process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring("#inputfiles#")
)

#process.source.fileNames = cms.untracked.vstring('file:/gpfs/csic_users/jfernan/latinosYieldSkim.root')
process.source.fileNames = cms.untracked.vstring('file:/gpfs/csic_users/albertog/fromJaviF_STOP_newProducers_18mar2013/latinosYieldSkim.root')
#process.source.fileNames = cms.untracked.vstring('file:/gpfs/res_projects/csic/cms-lhc/jfernan/NoSkim/SUSYPAT.root')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('143833:517-143833:522','143833:537-143833:540','143833:545-143833:546','143833:549-143833:552','143833:557-143833:560')
#import FWCore.PythonUtilities.LumiList as LumiList
#import FWCore.ParameterSet.Types as CfgTypes
#myLumis = LumiList.LumiList(filename = '/gpfs/csic_users/jfernan/CMSSW_5_3_1/src/WWAnalysis/Misc/Jsons/HCP.json').getCMSSWString().split(',')
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)

#process.source.noEventSort = cms.untracked.bool(True)
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

#Global Tag Stuff
process.GlobalTag.globaltag = 'GR_R_52_V7::All'
#process.GlobalTag.globaltag = 'GR_R_311_V2::All'



