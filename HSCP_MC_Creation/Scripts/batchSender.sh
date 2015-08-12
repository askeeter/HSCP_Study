#!/bin/bash

#Populate an array of all of the configuration files
shopt -s nullglob
filearray=( "HSCP_MC_sh_Files"/* )
shopt -u nullglob
#printf "%s\n" "${filearray[@]}"

#Create a bash file for each config file
for file in "${filearray[@]}"
do
    #Strip off the chracters that we dont need
    fileFixed=${file:17}
    #Send to the one day queue
    bsub -R "pool>20000" -q 2nd -J $fileFixed < /afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_sh_Files/$fileFixed
    
    #echo ${file:17}
done