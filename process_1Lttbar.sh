
TAG=$1

mkdir -p hists/${TAG}/
rm -rf hists/${TAG}/*

inputBase=/eos/user/n/nlu/HH4b/20210712_regression/option10/combined/1LTopSkim/

doSyst=no

if [ "$#" -gt 1 ]; then
    doSyst=$2
fi

for year in 2016 2017 2018

do
 (set -x ;./HHLooper ${inputBase}/${year}/qcd/  qcd.root ${TAG} 0 ${doSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/tt1L/ tt1L.root ${TAG} 0 ${doSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/tt2L/ tt2L.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/ttbar/ ttbar.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/HHc1/ HHc1.root ${TAG} 0 ${doSyst} >&1) &
 #(set -x ;./HHLooper ${inputBase}/${year}/WJets/ WJets.root ${TAG} 0 ${doSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/others/ others.root ${TAG} 0 ${doSyst} >&1) &
 (set -x ;./HHLooper ${inputBase}/${year}/data/ data.root ${TAG} 1 ${doSyst} >&1)
#(set -x ;./HHLooper ${inputBase}/${year}/tW/ tW.root ${TAG} 0 ${doSyst} >&1) &

done

sleep 1
wait

mkdir -p hists/${TAG}/combine/
hadd -k -f hists/${TAG}/combine/data.root hists/${TAG}/2016/data.root hists/${TAG}/2017/data.root hists/${TAG}/2018/data.root
hadd -k -f hists/${TAG}/combine/qcd.root hists/${TAG}/2016/qcd.root hists/${TAG}/2017/qcd.root hists/${TAG}/2018/qcd.root
hadd -k -f hists/${TAG}/combine/tt1L.root hists/${TAG}/2016/tt1L.root hists/${TAG}/2017/tt1L.root hists/${TAG}/2018/tt1L.root
hadd -k -f hists/${TAG}/combine/tt2L.root hists/${TAG}/2016/tt2L.root hists/${TAG}/2017/tt2L.root hists/${TAG}/2018/tt2L.root
#hadd -k -f hists/${TAG}/combine/HHc1.root hists/${TAG}/2016/HHc1.root hists/${TAG}/2017/HHc1.root hists/${TAG}/2018/HHc1.root
#hadd -k -f hists/${TAG}/combine/WJets.root hists/${TAG}/2016/WJets.root hists/${TAG}/2017/WJets.root hists/${TAG}/2018/WJets.root
hadd -k -f hists/${TAG}/combine/others.root hists/${TAG}/2016/others.root hists/${TAG}/2017/others.root hists/${TAG}/2018/others.root
#hadd -k -f hists/${TAG}/combine/tW.root hists/${TAG}/2016/tW.root hists/${TAG}/2017/tW.root hists/${TAG}/2018/tW.root


