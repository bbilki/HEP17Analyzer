# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("P5HFAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

import sys
runNumber = sys.argv[2]
runType = int(sys.argv[3])

process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
        'file:./Data/USC_'+runNumber+'.root'
    )
)

#process.source = cms.Source("HcalTBSource",fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/group/dpg_hcal/comm_hcal/SX5/HF/PMT/SX5_'+runNumber+'.root'))

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False)
        )

process.tbunpack = cms.EDProducer("HcalTBObjectUnpacker",
        IncludeUnmatchedHits = cms.untracked.bool(False),
	ConfigurationFile=cms.untracked.string('HCALCommissioning2017/HEP17Analyzer/test/configQADCTDC.txt'),
        HcalSlowDataFED = cms.untracked.int32(14),
        HcalTriggerFED = cms.untracked.int32(1),
        HcalTDCFED = cms.untracked.int32(8),
        HcalQADCFED = cms.untracked.int32(8),
        fedRawDataCollectionTag = cms.InputTag('source')
)

process.hcalDigis = cms.EDProducer("HcalRawToDigi",
                                   #       UnpackHF = cms.untracked.bool(True),
                                   ### Falg to enable unpacking of TTP channels(default = false)
                                   ### UnpackTTP = cms.untracked.bool(True),
                                   FilterDataQuality = cms.bool(False),
                                   InputLabel = cms.InputTag('source'),
                                   #HcalFirstFED = cms.untracked.int32(928),
                                   HcalFirstFED = cms.untracked.int32(1100),
                                   #ComplainEmptyData = cms.untracked.bool(False),
                                   ComplainEmptyData = cms.untracked.bool(True),
                                   #       UnpackCalib = cms.untracked.bool(True),
                                   #FEDs = cms.untracked.vint32(928,930,932),
                                   #FEDs = cms.untracked.vint32(928,930,932,938),
                                   #FEDs = cms.untracked.vint32(928),
                                   FEDs = cms.untracked.vint32(1100,1102,1104,1106,1108,1110,1112,1114,1116,1118,1119,1120,1121,1122,1123,1134),
                                   firstSample = cms.int32(0),
                                   lastSample = cms.int32(10)
                                   )


process.hcalAnalyzer = cms.EDAnalyzer('HEP17Analyzer',
        OutFileName = cms.untracked.string('N_'+runNumber+'.root'),
        RunType = cms.int32(runType),
	histoFED =  cms.int32(74)
)

process.load('Configuration.Geometry.GeometryIdeal_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
#from CondCore.DBCommon.CondDBSetup_cfi import *

from CondCore.CondDB.CondDB_cfi import *

process.GlobalTag.globaltag = autoCond['startup'] 

#   EMAP Needed for H2 DATA
process.es_ascii = cms.ESSource('HcalTextCalibrations',
        input = cms.VPSet(
               cms.PSet(
                object = cms.string('ElectronicsMap'),
		file = cms.FileInPath('HCALCommissioning2017/HEP17Analyzer/test/emap_HEP17.txt')
               )
        )
)

process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.p = cms.Path(process.hcalDigis*process.hcalAnalyzer)


