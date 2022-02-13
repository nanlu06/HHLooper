
########This is the instruction to reproduce Figure 7 to  Figure 25 in this version of  AN-20-201: https://zhicaiz.web.cern.ch/zhicaiz/sharebox/notes/AN-20-201-20210115.pdf


#### first, run then tupler in the following order:

cd ..
source process.sh yield_AN_sr
source process.sh yield_AN_sr_sys yes
source process_1Lttbar.sh yield_AN_1Lttbar
source process_ttbar.sh yield_AN_ttbar_cor
source process_ttbar.sh yield_AN_ttbar_cor_sys yes

###now comment Line 298, and uncomment Line  299 in ../app/HHLooper.cc to remove recoil correction;  
make clean; make
source process_ttbar.sh yield_AN_ttbar
###now comment line 299 and uncommend line 298  in ../app/HHLooper.cc
###now make data card:
cd python
python prepare_card.py yield_AN_sr_sys v24
python prepare_card.py yield_AN_sr_sys v8p2
python prepare_card.py yield_AN_ttbar_cor_sys v24
python prepare_card.py yield_AN_ttbar_cor_sys v8p2


###now begin to make plots:

####Figure 7 (J2 mass in SR for BDT v1)
TAG=yield_AN_sr
for region in SRv24Bin1 SRv24Bin2 SRv24Bin3 SRv24Bin4
do
    python plot.py -i ../hists/${TAG}/combine/ ${region}__fatJet2MassSD -n 30 -d -bd -O ${region}__fatJet2MassSD_SideBand
done


####Figure 8 (J2 mass in SR for BDT v2)
TAG=yield_AN_sr
for region in SRv8p2Bin1 SRv8p2Bin2 SRv8p2Bin3
do
    python plot.py -i ../hists/${TAG}/combine/ ${region}__fatJet2MassSD -n 30 -d -bd -O ${region}__fatJet2MassSD_SideBand
done



####1L top plots, Figure 9, 10, 11

TAG=yield_AN_1Lttbar
for year in 2016 2017  2018
do
    for var in MET fatJet1Tau3OverTau2 fatJet1MassSD  
    do
        python plot_1Lttbar.py -i ../hists/${TAG}/${year}/ -w CutLepEta__${var} -n 25 -d
        python plot_1Lttbar.py -i ../hists/${TAG}/${year}/ -w TTBarLepJetCR__${var} -n 25 -d
    done
    python plot_1Lttbar.py -i ../hists/${TAG}/${year}/ -w TTBarLepJetCR__fatJet1PNetXbb_Bin -d -au
done


####ttbar recoil plots, Figure 12
TAG=yield_AN_ttbar
for year in 2016 2017  2018
do
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 TTBarCR__hh_pt  -d 
done
#Figure 13 is just a manual fit from  Figure 12


#####ttbar BDT shapes, Figure 14
TAG=yield_AN_ttbar_cor
python plot.py -i ../hists/${TAG}/combine/ -s 1 -w TTBarCRTight__EventBDTv24 -d 
python plot.py -i ../hists/${TAG}/combine/ -s 1 -w TTBarCRTight__EventBDTv8p2 -d 


#####top mass peak vs BDT, figure 15 and 16
TAG=yield_AN_ttbar_cor
for BDTv in v24 v8p2
do
    for  BDTbin in  1 2 3 4
    do
        python plot.py -i ../hists/${TAG}/combine/ -s 1 -w TTBarCRBDT${BDTbin}${BDTv}__fatJet2MassSD -n 20 -d
    done
done


#### top mass vs BDT fit, Figure 17
python J2Mass_fit_data_vs_MC.py 
#change BDTbins and  x_bounds (Line 15-19) in  the code to choose  v24 or v8p2 BDT 


#### J2 SD mass distributions in ttbar control region  after all corrections, Figure 18
TAG=yield_AN_ttbar_cor
for year in 2016 2017 2018 combine
do
    python plot.py -i ../hists/${TAG}/${year}/ -s 1 -w TTBarCR__fatJet2MassSD -n 20 -d
done



#### J2 mass in pass  and fail regions, figure 19 and 20
TAG=yield_AN_sr_sys
for region in \
SRv24Bin1 \
SRv24Bin2 \
SRv24Bin3 \
SRv24Bin4 \
SRv8p2Bin1 \
SRv8p2Bin2 \
SRv8p2Bin3
do
    python plot.py -i ../hists/${TAG}/combine/ ${region}__fatJet2MassSD -d -bd -O ${region}__fatJet2MassSD_SideBand
done

for region in FailSRv24 FailSRv8p2
do
    python plot.py -i ../hists/${TAG}/combine/ ${region}__fatJet2MassSD -d -O ${region}__fatJet2MassSD
done


#### J2 mass for QCD events in pass and fail regions and the transfer factors, Figure 21 and 22
python plot_qcd_shape_pass_fail.py
#now change vbdt  to "v8p2" in L94 and do the above step again

