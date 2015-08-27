
        #!/bin/bash

        #Populate an array of all of the batch scripts 
        shopt -s nullglob
        filearray=( "HSCP_MC_sh_Files"/* )
        shopt -u nullglob

        for file in "${filearray[@]}"
        do
            #Strip off the chracters that we dont need
            fileFixed=${file:17}
            #Send to the two day queue. This can be changed
            bsub -R "pool>20000" -q 2nd -J $fileFixed < /afs/cern.ch/work/a/askeeter/private/CMSSW_7_4_4_patch4/src/HSCP_MC_sh_Files/$fileFixed
        done
