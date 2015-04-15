# To run execute
# cmsRun tmttt_tf_analysis_cfg.py someList=/path/to/sample.txt outputFile=/path/to/file.root

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load('Configuration.Geometry.GeometryExtended2023Pixel_cff')
process.load('Configuration.Geometry.GeometryExtended2023PixelReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

options = VarParsing.VarParsing ('analysis')
options.outputFile = 'test.root'
options.register('someList', 'Samples/Muons/PU0.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
options.register('HBins', 25, VarParsing.VarParsing.multiplicity.singleton,  VarParsing.VarParsing.varType.int, "Number of HT Array Bins")

options.parseArguments()

list = FileUtils.loadListFromFile(options.someList)
readFiles = cms.untracked.vstring(*list)

secFiles = cms.untracked.vstring()

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
             )


process.source = cms.Source("PoolSource",
                            fileNames = readFiles,
#                            fileNames = cms.untracked.vstring('file:ttbar_PU140.root'),
                            #fileNames = cms.untracked.vstring('/store/user/dcieri/MuonsPu140_test-newbias/Muons/MuonsPu140_test-newbias/150324_100930/0000/step2_1.root'),
                            secondaryFileNames = secFiles
                            )

process.demo = cms.EDAnalyzer('TMTrackFinder',
        HoughPrint   = cms.bool(False),  # Print the Hough Transform histograms (Need to change some lines in TMTrackFinderAlgorithm.cc and .h)
        Debug        = cms.bool(False),                      # Print some Debug information
        StubMinPt    = cms.double(-2.0),                     # Select just stubs which have a truth pT > StubMinPt (-2 if you want to keep all stubs)
        StubEtaRange = cms.double(2.55),                     # Check that the beta segment is within the CMS Tracker acceptance
        EtaRegions   = cms.vdouble(-3.2,-1.45,-0.61,0.61,1.45,3.2),   #Define the 5 trigger regions
        BeamErrorZ   = cms.double(15),                       #
        BetaRange    = cms.double(0.2),                      # Lenght of the Beta segment
        MinHoughPt   = cms.double(3.0),                      # Min stubs pT to use in the Hough Transform
        HoughBinsX   = cms.int32(options.HBins),             # N. of Hough X bins
        HoughBinsY   = cms.int32(options.HBins),             # N. of Hough Y bins
        SampleType   = cms.int32(0),                         # 0 For Muons, 1 For Electrons. To be used only  with RECO datasets
        PtLow        = cms.double(3),                        # Limits for the pT histograms
        PtHigh       = cms.double(20),
        RECO         = cms.bool(False),                      # True if you are analyzing a RECO sample
        EtaRegister  = cms.bool(False)                       # Enable the new register cut in the transition region
                              )

secFiles.extend([
        ])
process.p = cms.Path(process.demo)
