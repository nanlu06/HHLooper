#weight sys, trigger sys, spectral treatment for the ttbar bkg in BDT bin1, isblinded?

#preapproval final trig SF  in 2017
python prepare_card_SR_final.py yield_AN_sr_sys_0830 yield_AN_sr_sys_0830_trig v8p2 /storage/af/user/nlu/hh/looper_output/hists/result0908/ PNetp9 yes

#new trig SF in 2016
python prepare_card_SR_final.py yield_AN_sr_sys_0830 yield_AN_sr_sys_0830_trig v8p2 /storage/af/user/nlu/hh/looper_output/hists/result1027/ PNetp9 no

#new trig correction and uncertainty, new ttbar recoil correction
python prepare_card_SR_final.py yield_AN_sr_sys_0830 yield_AN_sr_sys_0830_trig v8p2 /storage/af/user/nlu/hh/looper_output/hists/result1105/ PNetp9 no
