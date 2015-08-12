import sys, os
import FWCore.ParameterSet.Config as cms
#Makes EDM from AOD
isSignal = True
isBckg = False
isData = False
isSkimmedSample = False
GTAG = 'MCRUN2_74_V8'
OUTPUTFILE = '/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_Root_Files/mchamp24_M_300_EDM.root'

InputFileList = cms.untracked.vstring(
'file:/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_Root_Files/mchamp24_M_300_AOD.root'
)

execfile( '/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py' )
