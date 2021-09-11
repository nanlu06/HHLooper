#weight sys, trigger sys, spectral treatment for the ttbar bkg in BDT bin1
#python prepare_card_SR_final.py yield_AN_sr_sys_0830oldVJetsXS yield_AN_sr_sys_0819_trig v8p2 /storage/af/user/nlu/hh/looper_output/hists/ PNetp9

python prepare_card_SR_final_ntuple0403.py yield_AN_sr_sys_0824_ntuple_20210403_HHLooper_sysTest_pT300_ttbarbin1PNet9 yield_AN_sr_sys_0819_trig v8p2 /storage/af/user/nlu/hh/looper_output/hists/ PNetp9

#final trig SF  in 2017
python prepare_card_SR_final.py yield_AN_sr_sys_0830 yield_AN_sr_sys_0830_trig v8p2 /storage/af/user/nlu/hh/looper_output/hists/result0908/ PNetp9

python prepare_card_SR_final.py yield_AN_sr_sys_0830 yield_AN_sr_sys_0830_trig v8p2 /storage/af/user/nlu/hh/looper_output/hists/ PNetp9
