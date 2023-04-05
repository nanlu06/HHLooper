TAG=yield_AN_sr_sys_0830_nominal_nosys
for region in SRv8p2Bin1 SRv8p2Bin2 SRv8p2Bin3
do
    python plot.py -i ../hists/${TAG}/combine/ ${region}__fatJet2MassSD -x 0,250 -n 30 -d -bd -O ${region}__fatJet2MassSD_SideBand
done
