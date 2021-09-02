#!/bin/bash
cd /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/condor_logs
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
#cd $TMP
cp -r /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/data .
echo "/storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/HHLooper $1 $2 $3 $4 $5 $6 $7 $8"
/storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/HHLooper $1 $2 $3 $4 $5 $6 $7 $8

#$1 -  ${inputBase}/${year}/${process}
#$2 - outputfilename
#$3 ${TAG}
#$4 isData
#$5 - ${doSyst} 
#$6 - ${doShapeSyst}
#$7 -  ${doTrigSyst} 
#$8 - ${doPNetSFSyst}
