#output of plot_1Lttbar.py and plot.py will be in python/plots folder by default if no output folder directory is given directly
#figure 16,17,18
#I added the red line in the AN by hand

#figure17
#\includegraphics[width=0.32\textwidth]{figs/yield_AN_1Lttbar/2016/TTBarLepJetCR__fatJet1MassSDLR_linY.pdf}
#\includegraphics[width=0.32\textwidth]{figs/yield_AN_1Lttbar/2017/TTBarLepJetCR__fatJet1MassSDLR_linY.pdf}
#\includegraphics[width=0.32\textwidth]{figs/yield_AN_1Lttbar/2018/TTBarLepJetCR__fatJet1MassSDLR_linY.pdf}

#figure18
#\includegraphics[width=0.32\textwidth]{figs/yield_AN_1Lttbar/2016/TTBarLepJetCR__fatJet1PNetXbb_Bin2_linY.pdf}
#\includegraphics[width=0.32\textwidth]{figs/yield_AN_1Lttbar/2017/TTBarLepJetCR__fatJet1PNetXbb_Bin2_linY.pdf}
#\includegraphics[width=0.32\textwidth]{figs/yield_AN_1Lttbar/2018/TTBarLepJetCR__fatJet1PNetXbb_Bin2_linY.pdf}

TAG=yield_AN_1Lttbar_reg
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
TAG=yield_AN_ttbar_ptcheck_1Nov
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

#Figure 20 is a manual fit from  Figure 19 using script scripts/ttbarscripts/PNetcorrection.C

#recoil plots, Figure 21,22,23
TAG=yield_AN_ttbar_ptcheck_1Nov_cor
for year in 2016 2017 2018
do
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__hh_pt  -d 
    #other BDT vars
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__hh_eta  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__hh_mass  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__met  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1Tau3OverTau2  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet2Tau3OverTau2  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1MassSD  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1Eta  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetXbb  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetQCDb  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetQCDbb  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1PNetQCDothers  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__ptj1_over_mhh  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__ptj2_over_mhh  -d
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__ptj2_over_ptj1  -d

    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet2Pt  -d -n 40
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__fatJet1Pt  -d -n 40
done

#figure 24
python plot.py -i ../hists/${TAG}/combine/ -s 1 TTBarCR__hh_mass_lr  -d

#Figure 25 top mass peak vs BDT
#\includegraphics[width=0.4\textwidth]{figs/yield_AN_ttbar_cor/combine/TTBarCRBDT1v8p2__fatJet2MassSD_linY.pdf}
#\includegraphics[width=0.4\textwidth]{figs/yield_AN_ttbar_cor/combine/TTBarCRBDT2v8p2__fatJet2MassSD_linY.pdf}
#\includegraphics[width=0.4\textwidth]{figs/yield_AN_ttbar_cor/combine/TTBarCRBDT3v8p2__fatJet2MassSD_linY.pdf}
#\includegraphics[width=0.4\textwidth]{figs/yield_AN_ttbar_cor/combine/TTBarCRBDT4v8p2__fatJet2MassSD_linY.pdf}

TAG=yield_AN_ttbar_reg_cor
for BDTv in v8p2
do
    for  BDTbin in  1 2 3 4
    do
        python plot.py -i ../hists/${TAG}/combine/ -s 1 -w TTBarCRBDT${BDTbin}${BDTv}__fatJet2MassSD -n 40 -d
	python plot.py -i ../hists/${TAG}/combine/ -s 1 -w TTBarCRBDT${BDTbin}${BDTv}__fatJet2MassRegressed -n 40 -d
    done
done

#figure 26 scripts/ttbar_scripts/fit_topmass_inBDTbins.C

#####ttbar BDT shapes, Figure 28
#\includegraphics[width=0.45\textwidth]{figs/yield_AN_ttbar_cor/combine/TTBarCRTight__EventBDTv8p2_logY.pdf}
#\includegraphics[width=0.45\textwidth]{figs/yield_AN_ttbar_cor/combine/TTBarCRTight__EventBDTv8p2v3_logY.pdf}

#figure 25
python prepare_card_ttbar_jet.py yield_AN_ttbar_reg_cor v8p2

#the output is HHTo4BPlots_Run2_ttbarSkim_BDTv8p2.root, now go to the combine-hh directory and follow the instructions in the README.md on how to make the datacards figure 27 in ANv4. How to install combined-hh package see https://github.com/LPC-HH/combine-hh

#then go back to HHLooper director do the following command to produce figure 27:

# python makePostFitPlot_TTCR.py directory of the fitDiagnosticsSBfitonly.root file

# for example: python makePostFitPlot_TTCR.py /storage/af/user/nlu/work/HH/CMSSW_10_2_13/src/combine-hh/cards_shapes_TTBarCR/HHModel/fitDiagnosticsSBfitonly.root

#figure 28
TAG=yield_AN_ttbar_reg_cor
python plot.py -i ../hists/${TAG}/combine/ -s 1 -w TTBarCRTight__EventBDTv8p2 -d


