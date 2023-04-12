TAG=$1

mkdir -p hists/${TAG}/
rm -rf hists/${TAG}/*

inputBase=/eos/cms/store/group/phys_susy/razor/Run3Analysis/HH/HHTo4BNtupler/option5/nano/run3/combined/

./HHLooper ${inputBase}/2022/data/JetMET_2022C.root data_2022C.root ${TAG} 1 0 0 0 0 >&1

#hadd -k -f hists/${TAG}/2016/data.root hists/${TAG}/2016/data_*.root 
#hadd -k -f hists/${TAG}/2017/data.root hists/${TAG}/2017/data_*.root 
#hadd -k -f hists/${TAG}/2018/data.root hists/${TAG}/2018/data_*.root
#mkdir hists/${TAG}/combine/
#hadd -k -f hists/${TAG}/combine/data.root hists/${TAG}/2016/data.root hists/${TAG}/2017/data.root hists/${TAG}/2018/data.root
