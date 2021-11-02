
TAG=$1

#mkdir -p hists/${TAG}/
#rm -rf hists/${TAG}/*

#inputBase=/storage/af/user/nlu/work/HH/ntuples/20210310/option10/combined/Tau3Over20p54Skim/BDT/
inputBase=/storage/af/user/idutta/work/HH/ntuple/20210712_regression/option10/combined/Tau32TopSkim/BDT/
#/eos/cms/store/group/phys_susy/razor/Run2Analysis/HHBoost/Tau3Over20p54Skim/
#inputBase=/eos/cms/store/group/phys_susy/razor/Run2Analysis/HHBoost/HHTo4BNtupler_20200902_option5_combined_BDT_Jet2Xbb0p8Skim/ 
#inputBase=data/Jet2Xbb0p8Skim/
#inputBase=data/Option5NonSkim/
#inputBase=data/Tau3Over20p54Skim/

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

for year in 2016 2017 2018

do
 (set -x ;./HHLooper ${inputBase}/${year}/data/ data.root  ${TAG} 1 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/qcd/  qcd.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/ttbar/ ttbar.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &

 #(set -x ;./HHLooper ${inputBase}/${year}/ttbarMtt700to1000/ ttbarMtt700to1000.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/ttbarMtt1000toInf/ ttbarMtt1000toInf.root ${TAG} 0 ${doSyst} >&1) &

 #(set -x ;./HHLooper ${inputBase}/${year}/HHSM/ HHSM.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/HHc0/ HHc0.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/HHc1/ HHc1.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/HHc2p45/ HHc2p45.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/HHc5/ HHc5.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/VH/ VH.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/Higgs/ Higgs.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/ttH/ ttH.root ${TAG} 0 ${doSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/others/ others.root ${TAG} 0 ${doSyst} ${doShapeSyst} ${doTrigSyst} ${doPNetSFSyst} >&1) &

 #(set -x ;./HHLooper ${inputBase}/${year}/tW/ tW.root ${TAG} 0 ${doSyst} >&1) &

done

sleep 1
wait

mkdir -p hists/${TAG}/combine/
hadd -k -f hists/${TAG}/combine/data.root hists/${TAG}/2016/data.root hists/${TAG}/2017/data.root hists/${TAG}/2018/data.root
hadd -k -f hists/${TAG}/combine/qcd.root hists/${TAG}/2016/qcd.root hists/${TAG}/2017/qcd.root hists/${TAG}/2018/qcd.root
hadd -k -f hists/${TAG}/combine/ttbar.root hists/${TAG}/2016/ttbar.root hists/${TAG}/2017/ttbar.root hists/${TAG}/2018/ttbar.root
#hadd -k -f hists/${TAG}/combine/HHSM.root hists/${TAG}/2016/HHSM.root hists/${TAG}/2017/HHSM.root hists/${TAG}/2018/HHSM.root
#hadd -k -f hists/${TAG}/combine/HHc0.root hists/${TAG}/2016/HHc0.root hists/${TAG}/2017/HHc0.root hists/${TAG}/2018/HHc0.root
#hadd -k -f hists/${TAG}/combine/HHc1.root hists/${TAG}/2016/HHc1.root hists/${TAG}/2017/HHc1.root hists/${TAG}/2018/HHc1.root
#hadd -k -f hists/${TAG}/combine/HHc2p45.root hists/${TAG}/2016/HHc2p45.root hists/${TAG}/2017/HHc2p45.root hists/${TAG}/2018/HHc2p45.root
#hadd -k -f hists/${TAG}/combine/HHc5.root hists/${TAG}/2016/HHc5.root hists/${TAG}/2017/HHc5.root hists/${TAG}/2018/HHc5.root
#hadd -k -f hists/${TAG}/combine/VH.root hists/${TAG}/2016/VH.root hists/${TAG}/2017/VH.root hists/${TAG}/2018/VH.root
#hadd -k -f hists/${TAG}/combine/Higgs.root hists/${TAG}/2016/Higgs.root hists/${TAG}/2017/Higgs.root hists/${TAG}/2018/Higgs.root
#hadd -k -f hists/${TAG}/combine/ttH.root hists/${TAG}/2016/ttH.root hists/${TAG}/2017/ttH.root hists/${TAG}/2018/ttH.root
hadd -k -f hists/${TAG}/combine/others.root hists/${TAG}/2016/others.root hists/${TAG}/2017/others.root hists/${TAG}/2018/others.root
#hadd -k -f hists/${TAG}/combine/tW.root hists/${TAG}/2016/tW.root hists/${TAG}/2017/tW.root hists/${TAG}/2018/tW.root

#hadd -k -f hists/${TAG}/combine/bkg.root hists/${TAG}/combine/qcd.root hists/${TAG}/combine/ttbar.root hists/${TAG}/combine/VH.root hists/${TAG}/combine/Higgs.root hists/${TAG}/combine/ttH.root hists/${TAG}/combine/others.root 

