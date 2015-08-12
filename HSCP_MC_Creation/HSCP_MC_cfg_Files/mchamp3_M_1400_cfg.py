# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/ThirteenTeV/HSCPmchamp3_M_1400_TuneZ2star_13TeV_pythia6_cff.py --fileout file:mchamp3_M_1400_AOD.root --mc --eventcontent AODSIM --datatier GEN-SIM-DIGI-AOD --conditions MCRUN2_74_V8 --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:GRun,RAW2DIGI,L1Reco,RECO --python_filename mchamp3_M_1400_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi.customise,SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n 1000
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/ThirteenTeV/HSCPmchamp3_M_1400_TuneZ2star_13TeV_pythia6_cff.py nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-AOD'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('file:mchamp3_M_1400_AOD.root'),
    outputCommands = process.AODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V8', '')

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters'),
        processParameters = cms.vstring('MSEL=0          ! User defined processes', 
            'MSUB(1)=1 !', 
            'MSTP(43)    = 3   ! complete Z0/gamma* interference', 
            'MSTP(1)=4 !fourth generation', 
            'CKIN(1)=1400.000000 !min sqrt(s hat)', 
            'CKIN(2)= -1  ! (no) max sqrt(s hat) (GeV)', 
            'KCHG(17,1)=3 !charge of tauprime', 
            'PMAS(17,1)=1400.000000 !tauprime mass', 
            'MDME(174,1) = 0   !Z decay into d dbar', 
            'MDME(175,1) = 0   !Z decay into u ubar', 
            'MDME(176,1) = 0   !Z decay into s sbar', 
            'MDME(177,1) = 0   !Z decay into c cbar', 
            'MDME(178,1) = 0   !Z decay into b bbar', 
            'MDME(179,1) = 0   !Z decay into t tbar', 
            'MDME(180,1) = 0   !Z decay into bprime bprimebar', 
            'MDME(181,1) = 0   !Z decay into tprime tprimebar', 
            'MDME(182,1) = 0   !Z decay into e- e+', 
            'MDME(183,1) = 0   !Z decay into nu_e nu_ebar', 
            'MDME(184,1) = 0   !Z decay into mu- mu+', 
            'MDME(185,1) = 0   !Z decay into nu_mu nu_mubar', 
            'MDME(186,1) = 0   !Z decay into tau- tau+', 
            'MDME(187,1) = 0   !Z decay into nu_tau nu_taubar', 
            'MDME(188,1) = 1   !Z decay into tauprime tauprimebar', 
            'MDME(189,1) = 0   !Z decay into nu_tauprime nu_tauprimebar', 
            'MDCY(17,1)=0    ! set tauprime stable', 
            'MWID(17)=0      ! set tauprime width 0'),
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.921 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.227 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model')
    ),
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(-1),
    filterEfficiency = cms.untracked.double(1.0),
    hscpFlavor = cms.untracked.string('stau'),
    massPoint = cms.untracked.int32(1400),
    maxEventsToPrint = cms.untracked.int32(0),
    particleFile = cms.untracked.string('Configuration/Generator/data/particles_HIP3_stau_1400_GeV.txt'),
    pdtFile = cms.FileInPath('Configuration/Generator/data/hscppythiapdtHIP3stau1400.tbl'),
    processFile = cms.untracked.string('SimG4Core/CustomPhysics/data/RhadronProcessList.txt'),
    slhaFile = cms.untracked.string('None'),
    useregge = cms.bool(False)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SimG4Core.CustomPhysics.Exotica_HSCP_SIM_cfi
from SimG4Core.CustomPhysics.Exotica_HSCP_SIM_cfi import customise 

#call to customisation function customise imported from SimG4Core.CustomPhysics.Exotica_HSCP_SIM_cfi
process = customise(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforFullSim 

#call to customisation function customizeHLTforFullSim imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforFullSim(process)

# End of customisation functions

