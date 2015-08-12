#!/bin/bash

shopt -s nullglob
filearray=( "HSCP_MC_Root_Files"/*EDM* )
shopt -u nullglob
#printf "%s\n" "${filearray[@]}"
appendTo="/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/Analysis_Samples.txt"
#Create a python file for each config file    
for file in "${filearray[@]}"
do
    parts=(${file//_/ })
    charge=${parts[3]}
    #Extract the number from the charge
    chargeFixed=$(echo $charge | tr -dc '0-9')
    mass=${parts[5]}    
    aod_file="mchamp${chargeFixed}_M_${mass}_AOD.root"
    edm_file="mchamp${chargeFixed}_M_${mass}_EDM"
    gen_file="mchamp${chargeFixed}_M_${mass}"
    #We need to append to the Analysis_Samples.txt file
    cat >> /afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/Analysis_Samples.txt << EOF
"CMSSW_7_4",   2, "$gen_file"    ,"$edm_file"    , "MC: mchamp${chargeFixed} ${mass} GeV/#font[12]{c}^{2}" , "S10", $mass, +9.8480000000E+01, 0, 1.000, 1.000, 1.000
EOF
done
