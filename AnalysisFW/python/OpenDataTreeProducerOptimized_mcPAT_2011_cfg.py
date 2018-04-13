
# Forked from SMPJ Analysis Framework
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMPJAnalysisFW
# https://github.com/cms-smpj/SMPJ/tree/v1.0
# (further optimized to improve performance)


## Skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
import FWCore.Utilities.FileUtils as FileUtils

####
# test flavour
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()

from PhysicsTools.JetMCAlgos.AK5PFJetsMCFlavourInfos_cfi import ak5JetFlavourInfos
process.jetFlavourInfosAK5PFJets = ak5JetFlavourInfos.clone()

# To print Jets info using class printJetFlavourInfo in 5_3_X release cmssw/PhysicsTools/JetExamples/test/
#process.printEventAK5PFJets = cms.EDAnalyzer("printJetFlavourInfo",
#    jetFlavourInfos = cms.InputTag("jetFlavourInfosAK5PFJets")
#)
####
# test SV

# Impact Parameter Tag Collection Info
from RecoBTag.ImpactParameter.impactParameter_cfi import impactParameterTagInfos
process.impactParameterTagInfo_track = impactParameterTagInfos.clone()
#process.impactParameterTagInfosV2 = impactParameterTagInfos.clone()

# Secondary Vertes Tag Collection Info
from RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi import secondaryVertexTagInfos
process.secondaryVertexTagInfosV2 = secondaryVertexTagInfos.clone()
#process.secondaryVertexTagInfosV2.trackSelection.qualityClass = cms.string('any')
#process.secondaryVertexTagInfosV2.trackIPTagInfos = cms.InputTag("impactParameterTagInfosV2")
#process.SecondaryVertexTagInfosV2.trackIPTagInfos = "newImpactParameterTagInfos"

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# True : when running in OpenData virtual machine
# False: when runing in lxplus 
runOnVM = False
#runOnVM = True

# Local input
######################################################
# run with the bash script Full dataset
fileList = FileUtils.loadListFromFile(NAMEOFINPUTFILE)
#####################################################
#run partial dataset
#fileList = FileUtils.loadListFromFile('CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt-80to120_TuneZ2_7TeV_pythia6_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt')
process.source.fileNames = cms.untracked.vstring(*fileList)

if runOnVM:
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db')

# Global tag for Summer11LegDR-PU_S13_START53_LV6-v1
process.GlobalTag.globaltag = cms.string('START53_LV6A1::All')

# Select good vertices
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
    )

# -------- The Tracking failure filter ------#
from RecoMET.METFilters.trackingFailureFilter_cfi import trackingFailureFilter
process.trackingFailureFilter = trackingFailureFilter.clone()
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")

################### EDAnalyzer ##############################3
process.ak5ak7 = cms.EDAnalyzer('OpenDataTreeProducerOptimized',
    ## jet collections ###########################
    pfak7jets       = cms.InputTag('ak7PFJets'),
    pfak5jets       = cms.InputTag('ak5PFJets'),
    ## MET collection ####
    pfmet           = cms.InputTag('pfMET7'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string('AK7PF'),

    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('goodOfflinePrimaryVertices'),
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## rho #######################################
    srcPFRho        = cms.InputTag('kt6PFJets','rho'),
    ## preselection cuts #########################
    maxY            = cms.double(5.0), 
    minPFPt         = cms.double(30),
    minNPFJets      = cms.int32(1),
    minGenPt        = cms.untracked.double(30),
    minJJMass       = cms.double(-1),
    isMCarlo        = cms.untracked.bool(True),
    genjets         = cms.untracked.InputTag('ak7GenJets'),
    useGenInfo      = cms.untracked.bool(True),
    ## trigger ###################################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerNames    = cms.vstring(
                                'HLT_Jet30', 'HLT_Jet60', 'HLT_Jet80', 'HLT_Jet110', 
                                'HLT_Jet150','HLT_Jet190','HLT_Jet240', 'HLT_Jet300 ','HLT_Jet370', # HLT_Jet300 added
                                ),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jet energy correction labels ##############
    jetCorr_ak5      = cms.string('ak5PFL1FastL2L3Residual'),
    jetCorr_ak7      = cms.string('ak7PFL1FastL2L3Residual'),
    
    ##Test flavour
    jetFlavourInfos = cms.InputTag("jetFlavourInfosAK5PFJets"),
    ##Test SV
    #impactParameterTagInfos = cms.InputTag("impactParameterTagInfosV2"),  
    secondaryVertexTagInfos = cms.InputTag("secondaryVertexTagInfosV2"),

    ##track IP 
    ipassociation = cms.InputTag("impactParameterTagInfo_track") 
)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_Jet*', 'HLT_DiJetAve*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)



# Let it run
process.p = cms.Path(
    #process.impactParameterTagInfosV2*    ### Test CSV  
    process.impactParameterTagInfo_track*  ## IP for tracks
    process.secondaryVertexTagInfosV2*  ### Test CSV
    process.goodOfflinePrimaryVertices*
    process.hltFilter *
    process.trackingFailureFilter *
    process.selectedHadronsAndPartons*  ### Test flavour
    process.jetFlavourInfosAK5PFJets*   ### Test flavour
    process.ak5ak7   
   
)

# Approximate processing time on VM (Intel Core i5-5300U 2.3GHz laptop):
# 50000 events per 1 hour (both for DATA and MC)

# Change number of events here:

process.maxEvents.input = -1

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

# Output file
#######################################################################################################################################
process.TFileService = cms.Service("TFileService", fileName = cms.string("/eos/user/b/bchazinq/QCDPt470to600/" + NAMEROOTOFOUTPUTROOT ))
#######################################################################################################################################

# To suppress long output at the end of the job
#process.options.wantSummary = False   

del process.outpath
