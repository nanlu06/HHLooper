TAG=$1

mkdir -p hists/${TAG}/
rm -rf hists/${TAG}/*

inputBase=/eos/user/n/nlu/HH4b/

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

hadd -k -f hists/${TAG}/2016/data.root hists/${TAG}/2016/data_*.root 
hadd -k -f hists/${TAG}/2017/data.root hists/${TAG}/2017/data_*.root 
hadd -k -f hists/${TAG}/2018/data.root hists/${TAG}/2018/data_*.root
mkdir hists/${TAG}/combine/
hadd -k -f hists/${TAG}/combine/data.root hists/${TAG}/2016/data.root hists/${TAG}/2017/data.root hists/${TAG}/2018/data.root
