#!/bin/bash
v=$2
proc=${v::-5}
#mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/$3_$proc
cd /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/$3_$proc
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
#cd $TMP
cp -r /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/data .
cp /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/HHLooper .
echo "./HHLooper $1 $2 $3 $4 $5 $6 $7 $8"
./HHLooper $1 $2 $3 $4 $5 $6 $7 $8

#$1 -  ${inputBase}/${year}/${process}
#$2 - outputfilename
#$3 ${TAG}
#$4 isData
#$5 - ${doSyst} 
#$6 - ${doShapeSyst}
#$7 -  ${doTrigSyst} 
#$8 - ${doPNetSFSyst}
