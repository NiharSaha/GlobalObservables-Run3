import FWCore.ParameterSet.VarParsing as VarParsing
#import CondCore.CondDB.CondDB_cfi

ivars = VarParsing.VarParsing('standard')

ivars.register ('outputTag',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="for testing")


ivars.register ('inputFile',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="for testing")

ivars.register ('outputFile',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="for testing")

#ivars.inputFile="./CentralityTable_HFtowers200_DataPbPb_usingMC_2024Run_HYDMC_xSF1.01_pphfCoincFilterPF3Th5_Threshold100_Nominal_Normalisation1000_4000_HIPhysicsRawPrime0_Oct22.root

##### for DATA ##########
ivars.inputFile="./CentralityTable_HFtowers200_DataPbPb_usingMC_2024Run_HYDMC_xSF1.01_pphfCoincFilterPF3Th5_Threshold50_SysDown_Normalisation1000_4000_HIPhysicsRawPrime0_usinghiHF_Oct29.root"
ivars.outputTag="CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal" #Has to be the inputFile tag
ivars.outputFile="./CentralityTable_HFtowers200_DataPbPb2024_periHYDJETshape_run3v140x01_offline_SysDown.db"




##### for MC ##########
#ivars.inputFile="./CentralityTable_hiHF200_HydjetCello_official_MC2024_Oct28.root"
#ivars.outputTag="CentralityTable_HFtowers200_HydjetCello_v1401x0_official_MC2024"
#ivars.outputFile="./CentralityTable_HFtowers200_HydjetCello_v1401x0_official_MC2024.db"


ivars.parseArguments()



import FWCore.ParameterSet.Config as cms

process = cms.Process('DUMMY')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("EmptyIOVSource",
                            timetype = cms.string("runnumber"),
                            firstValue = cms.uint64(1),
                            lastValue = cms.uint64(1),
                            interval = cms.uint64(1)
                        )

#process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.CondDB.CondDB_cfi")

#process.CondDBCommon.connect = "sqlite_file:" + ivars.outputFile
process.CondDB.connect = "sqlite_file:" + ivars.outputFile

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
                                          #process.CondDBCommon,
                                          process.CondDB,
                                          timetype = cms.untracked.string("runnumber"),
                                          toPut = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRcd'),
                                          tag = cms.string(ivars.outputTag)
                                          )
                                          )
)



process.makeCentralityTableDB = cms.EDAnalyzer('CentralityTableProducer',
                                               makeDBFromTFile = cms.untracked.bool(True),
                                               inputTFile = cms.string(ivars.inputFile),
                                               rootTag = cms.string(ivars.outputTag)
)


process.step  = cms.Path(process.makeCentralityTableDB)

