import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

#===================== Message Logger =============================
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
        limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(1)
            )
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


# Check for ny duplicates
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

runningOnMC = False
process.GlobalTag.globaltag = "START52_V11C::All"
if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_52_V9D::All"

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
          REPLACEWITHFILENAMES   

    )


  )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#=============== configure cleaning sequence ================================

process.lostLeptonTree = cms.EDAnalyzer("LostLeptonTree",
                                        Debug           = cms.bool    (False),
                                        VertexSource    = cms.InputTag('goodVertices'),
                                        DoPUReweight    = cms.bool    (True),
                                        PUWeigthSource  = cms.InputTag("puWeight"),
                                        PFMetSource     = cms.InputTag("patMETsPF"),
                                        JetAllSource    = cms.InputTag("patJetsPF"),
                                        #bTagName        = cms.string  ("trackCountingHighEffBJetTags"),
                                        bTagName        = cms.string  ("combinedSecondaryVertexBJetTags"),
                                        MHTSource       = cms.InputTag("mhtPFchs"),
                                        HTSource        = cms.InputTag("htPFchs"),
                                        MuonVetoSrc     = cms.InputTag("patMuonsPFIDIso"),
                                        EleVetoSrc      = cms.InputTag("patElectronsIDIso"),
                                        SaveAllMuons    = cms.bool    (True),
                                        MuonSource      = cms.InputTag('patMuonsPF'),
                                        MinMuPt         = cms.double  (5.0),
                                        MaxMuEta        = cms.double  (2.4),
                                        MaxMuD0         = cms.double  (0.2),
                                        MaxMuDz         = cms.double  (0.5),
                                        MaxMuRelIso     = cms.double  (0.20),
                                        PFCandidateSrc  = cms.InputTag("pfNoPileUpIsoPF"),
                                        #ElePFSrc        = cms.InputTag("patElectronsPFIDIso"),
                                        SaveAllElectrons  = cms.bool(True),
                                        ElectronSource    = cms.InputTag('gsfElectrons'),
                                        ConversionsSource = cms.InputTag("allConversions"),
                                        IsoValInputTags   = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                                                          cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                                                          cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
                                        RhoIsoSrcEle      = cms.InputTag("kt6PFJetsForIsolation", "rho"),
                                        BeamSpotSource    = cms.InputTag("offlineBeamSpot"),
                                        MinElePt          = cms.double(5.0)
)






# an example sequence to create skimmed susypat-tuples
process.analysisSeq = cms.Sequence(
    process.lostLeptonTree
)



process.load("RecoJets.JetProducers.HEPTopTagger_cfi")
from RecoJets.JetProducers.HEPTopTagger_cfi import HEPTopTagJets
from RecoJets.JetProducers.HEPTopTagger_cfi import HEPTopTagInfos

process.HEPTopTag15Jets = HEPTopTagJets.clone(
    )
process.HEPTopTag15Infos = HEPTopTagInfos.clone(
src = cms.InputTag("HEPTopTag15Jets")
)


from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca15PFJetsPFlow = ca4PFJets.clone(
rParam = cms.double(1.5),
src = cms.InputTag('particleFlow'),
doAreaFastjet = cms.bool(True),
doRhoFastjet = cms.bool(True),
Rho_EtaMax = cms.double(6.0),
Ghost_EtaMax = cms.double(7.0)
)
 


##============ configure output module configuration ========================
process.TFileService = cms.Service("TFileService",

                                   fileName = cms.string('OUTFILENAME.root')
)

process.load("SandBox.Skims.RA2Selection_cff")



process.stopPFMuonVeto = cms.Sequence(
     ~process.countPFMuonsIDIso
)

process.stopElectronVeto = cms.Sequence(
     ~process.countElectronsIDIso
)

process.load("SandBox.Stop.StopSelection_cff")



process.stopdPhiFilter = cms.Sequence(
    process.dPhiFilter
)


process.load("SandBox.Stop.trackIsolationMaker_cfi")
process.stopTrackIsolation = cms.Sequence(
    process.trackIsolationFilter
)

process.allEvents = cms.EDAnalyzer('EventCount'
)


process.afterMuonVetoEvents = cms.EDAnalyzer('EventCount'
)

process.afterElectronVetoEvents = cms.EDAnalyzer('EventCount'
)

process.afterdPhiFilterEvents = cms.EDAnalyzer('EventCount'
)

process.afterTrackIsolationEvents = cms.EDAnalyzer('EventCount'
)

#process.ppf = cms.Path( process.allEvents * process.stopPFMuonVeto* process.afterMuonVetoEvents * process.stopElectronVeto  * process.afterElectronVetoEvents* process.stopTrackIsolation * process.afterTrackIsolationEvents * process.stopdPhiFilter * process.afterdPhiFilterEvents *process.HEPTopTag15Jets * process.HEPTopTag15Infos *  process.ca15PFJetsPFlow * process.analysisSeq)
process.ppf = cms.Path( process.allEvents * process.stopPFMuonVeto* process.afterMuonVetoEvents * process.stopElectronVeto  * process.afterElectronVetoEvents* process.stopTrackIsolation * process.afterTrackIsolationEvents * process.stopdPhiFilter * process.afterdPhiFilterEvents * process.analysisSeq)

