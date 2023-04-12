TAG=yield_AN_sr_sys_0830_nominal_nosys
for region in SRv8p2Bin1 SRv8p2Bin2 SRv8p2Bin3
do
    python plot.py -i ../hists/${TAG}/combine/ ${region}__fatJet2MassSD -n 30 -d -bd -O ${region}__fatJet2MassSD_SideBand
done

TAG=yield_AN_1Lttbar
for year in 2016 2017 2018
do
    for var in MET fatJet1Tau3OverTau2 fatJet1MassSD
    do
        python plot_1Lttbar.py -i ../hists/${TAG}/${year}/ -w CutLepEta__${var} -n 25 -d
        python plot_1Lttbar.py -i ../hists/${TAG}/${year}/ -w TTBarLepJetCR__${var} -n 25 -d
    done
    python plot_1Lttbar.py -i ../hists/${TAG}/${year}/ -w TTBarLepJetCR__fatJet1PNetXbb_Bin -d -au
done

#recoil plots, Figure 19
TAG=yield_AN_ttbar_cor
for year in 2016 2017 2018
do
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__hh_pt  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1MassSD  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__MET  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__hh_mass  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetXbb  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetQCDb  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetQCDbb  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetQCDothers  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__ptj1_over_mhh  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__ptj2_over_mhh  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__ptj2_over_ptj1  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1Tau3OverTau2  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet2Tau3OverTau2  -d

done

#figure 25
python prepare_card_ttbar_jet.py yield_AN_ttbar_reg_cor v8p2

#the output is HHTo4BPlots_Run2_ttbarSkim_BDTv8p2.root, now go to the combine-hh directory and follow the instructions in the README.md on how to make the datacards figure 27 in ANv4. How to install combined-hh package see https://github.com/LPC-HH/combine-hh

#then go back to HHLooper director do the following command to produce figure 27:

# python makePostFitPlot_TTCR.py directory of the fitDiagnosticsSBfitonly.root file

# for example: python makePostFitPlot_TTCR.py /storage/af/user/nlu/work/HH/CMSSW_10_2_13/src/combine-hh/cards_shapes_TTBarCR/HHModel/fitDiagnosticsSBfitonly.root

#figure 28
TAG=yield_AN_ttbar_cor
python plot.py -i ../hists/${TAG}/combine/ -s 1 -w TTBarCRTight__EventBDTv8p2 -d
