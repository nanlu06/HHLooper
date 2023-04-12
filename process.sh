
TAG=$1

mkdir -p hists/${TAG}/
rm -rf hists/${TAG}/*
inputBase=/eos/cms/store/user/nlu/HH4b/

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

./HHLooper ${inputBase}/2016/data/JetHT_2016B-ver2_BDTs.root data_2016B.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2016/data/JetHT_2016C_BDTs.root data_2016C.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2016/data/JetHT_2016D_BDTs.root data_2016D.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2016/data/JetHT_2016E_BDTs.root data_2016E.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2016/data/JetHT_2016F_BDTs.root data_2016F.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2016/data/JetHT_2016G_BDTs.root data_2016G.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2016/data/JetHT_2016H_BDTs.root data_2016H.root ${TAG} 1 0 0 0 0 >&1

./HHLooper ${inputBase}/2017/data/JetHT_2017B_BDTs.root data_2017B.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2017/data/JetHT_2017C_BDTs.root data_2017C.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2017/data/JetHT_2017D_BDTs.root data_2017D.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2017/data/JetHT_2017E_BDTs.root data_2017E.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2017/data/JetHT_2017F_BDTs.root data_2017F.root ${TAG} 1 0 0 0 0 >&1

./HHLooper ${inputBase}/2018/data/JetHT_2018A_BDTs.root data_2018A.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2018/data/JetHT_2018B_BDTs.root data_2018B.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2018/data/JetHT_2018C_BDTs.root data_2018C.root ${TAG} 1 0 0 0 0 >&1
./HHLooper ${inputBase}/2018/data/JetHT_2018D_BDTs.root data_2018D.root ${TAG} 1 0 0 0 0 >&1

for year in 2016 2017 2018

do
 (set -x ;./HHLooper ${inputBase}/${year}/qcd/  qcd.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/ttbar/ ttbar.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHSM/ HHSM.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHc0/ HHc0.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHc1/ HHc1.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHc2p45/ HHc2p45.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHc5/ HHc5.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHVBFSM/ HHVBFSM.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHVBF0p511/ HHVBF0p511.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHVBF1p511/ HHVBF1p511.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHVBF110/ HHVBF110.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHVBF112/ HHVBF112.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHVBF121/ HHVBF121.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/HHVBF101/ HHVBF101.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/VH/ VH.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/Higgs/ Higgs.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/ttH/ ttH.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/others/VV/ vv.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/others/VJets/ vjets.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &

done

sleep 1
wait

mkdir -p hists/${TAG}/combine/
hadd -k -f hists/${TAG}/2016/data.root hists/${TAG}/2016/data_*.root 
hadd -k -f hists/${TAG}/2017/data.root hists/${TAG}/2017/data_*.root 
hadd -k -f hists/${TAG}/2018/data.root hists/${TAG}/2018/data_*.root
hadd -k -f hists/${TAG}/combine/data.root hists/${TAG}/2016/data.root hists/${TAG}/2017/data.root hists/${TAG}/2018/data.root
hadd -k -f hists/${TAG}/combine/qcd.root hists/${TAG}/2016/qcd.root hists/${TAG}/2017/qcd.root hists/${TAG}/2018/qcd.root
hadd -k -f hists/${TAG}/combine/ttbar.root hists/${TAG}/2016/ttbar.root hists/${TAG}/2017/ttbar.root hists/${TAG}/2018/ttbar.root
hadd -k -f hists/${TAG}/combine/HHSM.root hists/${TAG}/2016/HHSM.root hists/${TAG}/2017/HHSM.root hists/${TAG}/2018/HHSM.root
hadd -k -f hists/${TAG}/combine/HHc0.root hists/${TAG}/2016/HHc0.root hists/${TAG}/2017/HHc0.root hists/${TAG}/2018/HHc0.root
hadd -k -f hists/${TAG}/combine/HHc1.root hists/${TAG}/2016/HHc1.root hists/${TAG}/2017/HHc1.root hists/${TAG}/2018/HHc1.root
hadd -k -f hists/${TAG}/combine/HHc2p45.root hists/${TAG}/2016/HHc2p45.root hists/${TAG}/2017/HHc2p45.root hists/${TAG}/2018/HHc2p45.root
hadd -k -f hists/${TAG}/combine/HHc5.root hists/${TAG}/2016/HHc5.root hists/${TAG}/2017/HHc5.root hists/${TAG}/2018/HHc5.root
hadd -k -f hists/${TAG}/combine/HHVBFSM.root hists/${TAG}/2016/HHVBFSM.root hists/${TAG}/2017/HHVBFSM.root hists/${TAG}/2018/HHVBFSM.root
hadd -k -f hists/${TAG}/combine/HHVBF0p511.root hists/${TAG}/2016/HHVBF0p511.root hists/${TAG}/2017/HHVBF0p511.root hists/${TAG}/2018/HHVBF0p511.root
hadd -k -f hists/${TAG}/combine/HHVBF1p511.root hists/${TAG}/2016/HHVBF1p511.root hists/${TAG}/2017/HHVBF1p511.root hists/${TAG}/2018/HHVBF1p511.root
hadd -k -f hists/${TAG}/combine/HHVBF110.root hists/${TAG}/2016/HHVBF110.root hists/${TAG}/2017/HHVBF110.root hists/${TAG}/2018/HHVBF110.root
hadd -k -f hists/${TAG}/combine/HHVBF112.root hists/${TAG}/2016/HHVBF112.root hists/${TAG}/2017/HHVBF112.root hists/${TAG}/2018/HHVBF112.root
hadd -k -f hists/${TAG}/combine/HHVBF121.root hists/${TAG}/2016/HHVBF121.root hists/${TAG}/2017/HHVBF121.root hists/${TAG}/2018/HHVBF121.root
hadd -k -f hists/${TAG}/combine/HHVBF101.root hists/${TAG}/2016/HHVBF101.root hists/${TAG}/2017/HHVBF101.root hists/${TAG}/2018/HHVBF101.root
hadd -k -f hists/${TAG}/combine/VH.root hists/${TAG}/2016/VH.root hists/${TAG}/2017/VH.root hists/${TAG}/2018/VH.root
hadd -k -f hists/${TAG}/combine/Higgs.root hists/${TAG}/2016/Higgs.root hists/${TAG}/2017/Higgs.root hists/${TAG}/2018/Higgs.root
hadd -k -f hists/${TAG}/combine/ttH.root hists/${TAG}/2016/ttH.root hists/${TAG}/2017/ttH.root hists/${TAG}/2018/ttH.root
hadd -k -f hists/${TAG}/2016/others.root hists/${TAG}/2016/vv.root hists/${TAG}/2016/vjets.root 
hadd -k -f hists/${TAG}/2017/others.root hists/${TAG}/2017/vv.root hists/${TAG}/2017/vjets.root 
hadd -k -f hists/${TAG}/2018/others.root hists/${TAG}/2018/vv.root hists/${TAG}/2018/vjets.root 
hadd -k -f hists/${TAG}/combine/others.root hists/${TAG}/2016/others.root hists/${TAG}/2017/others.root hists/${TAG}/2018/others.root

hadd -k -f hists/${TAG}/combine/QCDggHVBF.root hists/${TAG}/combine/qcd.root hists/${TAG}/combine/Higgs.root

hadd -k -f hists/${TAG}/combine/bkg.root hists/${TAG}/combine/qcd.root hists/${TAG}/combine/ttbar.root  hists/${TAG}/combine/VH.root hists/${TAG}/combine/Higgs.root hists/${TAG}/combine/ttH.root hists/${TAG}/combine/others.root 

