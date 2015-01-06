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

	  "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888//susypat_419_1_ng8.root",
       "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888//susypat_162_1_FC2.root",
       "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888//susypat_277_1_sBw.root",
       "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888//susypat_339_1_9Tb.root",
       "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888//susypat_426_1_Ho7.root",
       "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888//susypat_367_1_51v.root"


#    "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DoubleMu_Run2012A-13Jul2012-v1_lpc1/vchetlur/DoubleMu/DoubleMu_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3_lpc1/062a2f358645c26cb65ce944f0ba30c4/susypat_129_1_3yW.root",
#     "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DoubleMu_Run2012A-13Jul2012-v1_lpc1/vchetlur/DoubleMu/DoubleMu_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3_lpc1/062a2f358645c26cb65ce944f0ba30c4/susypat_12_1_4i7.root",
#     "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DoubleMu_Run2012A-13Jul2012-v1_lpc1/vchetlur/DoubleMu/DoubleMu_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3_lpc1/062a2f358645c26cb65ce944f0ba30c4/susypat_130_1_mHz.root",
#    "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DoubleMu_Run2012A-13Jul2012-v1_lpc1/vchetlur/DoubleMu/DoubleMu_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3_lpc1/062a2f358645c26cb65ce944f0ba30c4/susypat_131_1_DBi.root"

#	"dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/sudan/DoubleMu/Run2012B-PromptReco-v1/82c5ab7277a61a422bf46e5043813259/susypat_999_1_V55.root",
 #       "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/sudan/DoubleMu/Run2012B-PromptReco-v1/82c5ab7277a61a422bf46e5043813259/susypat_99_3_fsN.root",
 #       "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/sudan/DoubleMu/Run2012B-PromptReco-v1/82c5ab7277a61a422bf46e5043813259/susypat_9_3_u9H.root",

   
#    "dcap:///cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcsusyhad/sudan/DoubleMu/Run2012A-PromptReco-v1/82c5ab7277a61a422bf46e5043813259/susypat_9_2_TrZ.root"

    )


  )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#=============== configure cleaning sequence ================================

process.STopAnalysisTree = cms.EDAnalyzer("DiSTopStudyTree",
                                        Debug           = cms.bool    (False),
                                        isDimuonSample           = cms.bool    (True),
					isMC           = cms.bool    (True),

                                        metSrc = cms.InputTag("patMETsPF"),
					jetSrc = cms.InputTag("patJetsAK5PFPt30"),
					evtWeightInput = cms.InputTag(""),
					vtxSrc = cms.InputTag("goodVertices"),


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
                                        MinMuPt         = cms.double  (20.0),
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
    process.STopAnalysisTree
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

                                   fileName = cms.string('outTest_T2tt_mStop550_mLSP250_NOCUTS_All.root')
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

#process.load("UserCode/TopTagger/

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

process.afterType3TopTaggerEvents = cms.EDAnalyzer('EventCount'
)


#process.load('SandBox.Skims.RA2Jets_cff')

#process.count5JetsAK5PFPt30 = process.countPatJets.clone()
#process.count5JetsAK5PFPt30.src = cms.InputTag('patJetsAK5PFPt30')
#process.count5JetsAK5PFPt30.minNumber = cms.uint32(5)


# Sequence to clean events based on boolean filter results stored on pattuples

#from SandBox.Skims.filterBoolean_cfi import *

process.load("SandBox.Skims.filterBoolean_cfi")

process.RA2_HBHENoiseFilterRA2    = process.booleanFilter.clone()
process.RA2_HBHENoiseFilterRA2.ResultSource = cms.InputTag("HBHENoiseFilterRA2","HBHENoiseFilterResult","PAT")

process.RA2_beamHaloFilter        = process.booleanFilter.clone()
process.RA2_beamHaloFilter.ResultSource = cms.InputTag("beamHaloFilter")

# Don't use this filter for now, needs to understand performance on physics events
process.RA2_eeNoiseFilter         = process.booleanFilter.clone()
process.RA2_eeNoiseFilter.ResultSource = cms.InputTag("eeNoiseFilter")

process.RA2_trackingFailureFilter = process.booleanFilter.clone()
process.RA2_trackingFailureFilter.ResultSource = cms.InputTag("trackingFailureFilter")

process.RA2_inconsistentMuons     = process.booleanFilter.clone()
process.RA2_inconsistentMuons.ResultSource = cms.InputTag("inconsistentMuons")

process.RA2_greedyMuons           = process.booleanFilter.clone()
process.RA2_greedyMuons.ResultSource = cms.InputTag("greedyMuons")

process.RA2_EcalTPFilter          = process.booleanFilter.clone()
process.RA2_EcalTPFilter.ResultSource = cms.InputTag("ra2EcalTPFilter")

process.RA2_EcalBEFilter          = process.booleanFilter.clone()
process.RA2_EcalBEFilter.ResultSource = cms.InputTag("ra2EcalBEFilter")

process.HcalLaserEventFilter      = process.booleanFilter.clone()
process.HcalLaserEventFilter.ResultSource = cms.InputTag("hcalLaserEventFilter")

process.EEBadScFilter             = process.booleanFilter.clone()
process.EEBadScFilter.ResultSource= cms.InputTag("eeBadScFilter")

process.EcalLaserFilter           = process.booleanFilter.clone()
process.EcalLaserFilter.ResultSource= cms.InputTag("ecalLaserCorrFilter")

process.RA2_PBNRFilter            = process.booleanFilter.clone()
process.RA2_PBNRFilter.ResultSource = cms.InputTag("ra2PBNR")

process.cleaningOnFilterResults  = cms.Sequence(
    process.RA2_HBHENoiseFilterRA2
  * process.RA2_beamHaloFilter
  * process.RA2_trackingFailureFilter
  * process.RA2_inconsistentMuons
  * process.RA2_greedyMuons
#  * process.RA2_EcalTPFilter
#  * process.RA2_EcalBEFilter
#  * process.HcalLaserEventFilter
#  * process.EEBadScFilter
  #* process.RA2_eeNoiseFilter
#  * process.EcalLaserFilter
  #* process.RA2_PBNRFilter #(to be applied jet-by-jet)
)



process.load("UserCode.TopTagger.topTagger_cfi")

##process.ppf = cms.Path( process.allEvents * process.stopPFMuonVeto* process.afterMuonVetoEvents * process.stopElectronVeto  * process.afterElectronVetoEvents* process.stopTrackIsolation * process.afterTrackIsolationEvents * process.stopdPhiFilter * process.afterdPhiFilterEvents *process.HEPTopTag15Jets * process.HEPTopTag15Infos *  process.ca15PFJetsPFlow * process.analysisSeq)


if (process.STopAnalysisTree.isMC==True) & (process.STopAnalysisTree.isDimuonSample==True):
   print 'is MC and Dimuon so loading reweiting module'
   process.QCDScaleWeightProducer = cms.EDProducer('ShowerWeightProducer'
   )
   process.ppf = cms.Path( process.allEvents * process.cleaningOnFilterResults *  process.afterMuonVetoEvents *  process.afterElectronVetoEvents* process.afterTrackIsolationEvents * process.afterdPhiFilterEvents *  process.afterType3TopTaggerEvents * process.QCDScaleWeightProducer * process.analysisSeq)
else:
   print 'No rewighting'
   process.ppf = cms.Path( process.allEvents * process.cleaningOnFilterResults *  process.afterMuonVetoEvents *  process.afterElectronVetoEvents* process.afterTrackIsolationEvents * process.afterdPhiFilterEvents *  process.afterType3TopTaggerEvents *process.analysisSeq)

