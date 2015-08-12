#!/bin/bash

#Populate an array of all of the configuration files
shopt -s nullglob
filearray=( "HSCP_MC_cfg_Files"/* )
shopt -u nullglob

#Create a bash file for each config file
for file in "${filearray[@]}"
do
    parts=(${file//_/ })
    charge=${parts[3]}
    #Extract the number from the charge
    chargeFixed=$(echo $charge | tr -dc '0-9')
    mass=${parts[5]}

    #All of the important data has been stripped from the config filename
    #Now to create the bath scripts
    filename="mchamp${chargeFixed}_M_${mass}.sh"
    cfgfile="mchamp${chargeFixed}_M_${mass}_cfg.py"
    rootfile="mchamp${chargeFixed}_M_${mass}_AOD.root"
    #Create an empty file to be filled
    touch $filename
    #Use echo to populate the file contents. Not the cleanest way, but it works for a file this short.
    echo "#!/bin/sh">$filename
    echo 'CMSSW_PROJECT_SRC="/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/"'>>$filename
    echo """CFG_FILE='/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_cfg_Files/${cfgfile}'""">>$filename
    echo """OUTPUT_FILE='/afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_Root_Files/$rootfile'""">>$filename
    echo "OUT_FILE='$rootfile'">>$filename
    echo 'TOP="$PWD"'>>$filename
    echo 'cd $CMSSW_PROJECT_SRC'>>$filename
    echo 'eval `scramv1 runtime -sh`'>>$filename
    echo 'cd $TOP'>>$filename
    echo 'cmsRun $CFG_FILE'>>$filename
    echo 'rfcp $OUT_FILE $OUTPUT_FILE'>>$filename
    #DO NOT FORGET to change the config file permissions if you are creating these by hand.
    chmod 744 $filename
done
