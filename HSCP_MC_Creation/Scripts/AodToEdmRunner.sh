#!/bin/bash

shopt -s nullglob
filearray=( "HSCP_MC_Root_Files"/*AOD* )
shopt -u nullglob
#printf "%s\n" "${filearray[@]}"


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
   
#Replace the standard configuration file with the one currently being ran
    cp /afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_AODtoEDM_Python_Files/${python_file} /afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_Signal_cfg.py
    cd /afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/
    cmsRun HSCParticleProducer_Signal_cfg.py
    wait
done