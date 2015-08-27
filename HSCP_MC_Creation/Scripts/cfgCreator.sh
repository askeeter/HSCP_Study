
      #!/bin/bash
      CHARGE=$1
      MASS=$2
      EVENTS=$3
      cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/HSCPmchamp${CHARGE}_M_${MASS}_TuneZ2star_13TeV_pythia6_cff.py --fileout file:mchamp${CHARGE}_M_${MASS}_AOD.root --mc --eventcontent AODSIM --datatier GEN-SIM-DIGI-AOD --conditions MCRUN2_74_V8 --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:GRun,RAW2DIGI,L1Reco,RECO --python_filename mchamp${CHARGE}_M_${MASS}_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi.customise,SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n ${EVENTS}
