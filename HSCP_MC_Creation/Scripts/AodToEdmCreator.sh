#!/bin/bash

#Populate an array of all of the AOD files
shopt -s nullglob
filearray=( "HSCP_MC_Root_Files"/*AOD* )
shopt -u nullglob

for file in "${filearray[@]}"
do
    parts=(${file//_/ })
    charge=${parts[3]}
    #Extract the number from the charge
    chargeFixed=$(echo $charge | tr -dc '0-9')
    mass=${parts[5]}    
    aod_file="mchamp${chargeFixed}_M_${mass}_AOD.root"
    root_file="mchamp${chargeFixed}_M_${mass}_EDM.root"
    python_file="mchamp${chargeFixed}_M_${mass}_cfg.py"
    
    cat > /afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_AODtoEDM_Python_Files/${python_file} << EOF
import sys, os
import FWCore.ParameterSet.Config as cms
#Makes EDM from AOD
isSignal = True
isBckg = False
isData = False
isSkimmedSample = False
GTAG = 'MCRUN2_74_V8'
OUTPUTFILE = '/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_Root_Files/${root_file}'

InputFileList = cms.untracked.vstring(
'file:/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_Root_Files/${aod_file}'
)

execfile( '${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py' )
EOF
    
done
