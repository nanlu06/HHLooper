
TAG=$1

mkdir -p hists/${TAG}/
rm -rf hists/${TAG}/*
inputBase=/home/nlu/ntuples/HH4brun3/combined/

doSyst=no
if [ "$#" -gt 1 ]; then
    doSyst=$2
fi

doShapeSyst=no
if [ "$#" -gt 2 ]; then
    doShapeSyst=$3
fi

doTrigSyst=no
if [ "$#" -gt 3 ]; then
    doTrigSyst=$4
fi

doPNetSFSyst=no
if [ "$#" -gt 4 ]; then
    doPNetSFSyst=$5
fi

./HHLooper ${inputBase}/2022/data/JetMET_2022C.root data_20122C.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2022/data/JetMET_2022D-v1.root data_2022D-v1.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2022/data/JetMET_2022D-v2.root data_2022D-v2.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2022/data/JetMET_2022E.root data_20122E.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2022/data/JetMET_2022F.root data_20122F.root ${TAG} 1 0 0 0 0 >&1

for year in 2022

do
 (set -x ;./HHLooper ${inputBase}/${year}/HHc1/ HHc1.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &

done

sleep 1
wait

mkdir -p hists/${TAG}/combine/
hadd -k -f hists/${TAG}/2022/data.root hists/${TAG}/2022/data_*.root 