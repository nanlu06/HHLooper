
TAG=$1 

python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet1SDMassLeft__FatJet1_msoftdrop -s 100 -b 1.5
python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet1SDMassRight__FatJet1_msoftdrop -s 100 -b 1.5 -R
python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet1Pt__FatJet1_pt -s 100 -b 1.5
python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet1DDB__FatJet1_btagDDBvL -s 100 -b 1.5

python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet2SDMassLeft__FatJet2_msoftdrop -s 100 -b 1.5
python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet2SDMassRight__FatJet2_msoftdrop -s 100 -b 1.5 -R
python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet2Pt__FatJet2_pt -s 100 -b 1.5
python plot.py -OP -i ../hists/${TAG}/ Nm1FatJet2DDB__FatJet2_btagDDBvL -s 100 -b 1.5

python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet1SDMassLeft__FatJet1_msoftdrop -n 40
python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet1SDMassRight__FatJet1_msoftdrop -n 40
python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet1Pt__FatJet1_pt -n 40
python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet1DDB__FatJet1_btagDDBvL -n 40

python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet2SDMassLeft__FatJet2_msoftdrop -n 40
python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet2SDMassRight__FatJet2_msoftdrop -n 40
python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet2Pt__FatJet2_pt -n 40
python plot.py -i ../hists/${TAG}/ -s 100 Nm1FatJet2DDB__FatJet2_btagDDBvL -n 40
