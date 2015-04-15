import FWCore.ParameterSet.Config as cms

# configuration to model pileup for initial physics phase
from SimGeneral.MixingModule.mixObjects_cfi import theMixObjects
from SimGeneral.MixingModule.mixPoolSource_cfi import *
from SimGeneral.MixingModule.digitizers_cfi import *


readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

mix = cms.EDProducer("MixingModule",
    digitizers = cms.PSet(theDigitizers),
    LabelPlayback = cms.string(''),
    maxBunch = cms.int32(3),
    minBunch = cms.int32(-12), ## in terms of 25 nsec

    bunchspace = cms.int32(25), ##ns
    mixProdStep1 = cms.bool(False),
    mixProdStep2 = cms.bool(False),

    playback = cms.untracked.bool(False),
    useCurrentProcessOnly = cms.bool(False),
    input = cms.SecSource("PoolSource",
        nbPileupEvents = cms.PSet(
            averageNumber = cms.double(60.)
        ),
        type = cms.string('fixed'),
	    sequential = cms.untracked.bool(True),
        fileNames = readFiles,
        secondaryFileNames = secFiles
        ),

    mixObjects = cms.PSet(theMixObjects)
)

readFiles.extend([
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_101_2_5X4.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_117_2_1OJ.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_150_1_P85.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_191_1_EUR.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_24_1_XS2.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_68_1_CEQ.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_76_1_DxE.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_92_2_tPw.root',
'/store/user/dcieri/MINBIAS_SIM/MINBIAS_SIM/3b9061f60047bdbeab21f315154b001d/minbias_9_1_V6E.root'
    ]);


            #            'file:minbias.root',',',',',,
#'file:minbias1.root',',',',',
#'file:file3.root',',',',',,
#'file:fileN.root',',',',',

secFiles.extend([
        ])


