#outputname no/yesforsyst systname no/yestrigsys no/yesPNetSFsys
#note the "nominal" here is essential and this format is used in the command used to make histograms for datacards 
source process.sh yield_AN_sr_sys_0830_nominal_nosys no nominal no no
source process_noData.sh yield_AN_sr_sys_0830_JER_Up no JER_Up no no
source process_noData.sh yield_AN_sr_sys_0830_JER_Down no JER_Down no no
source process_noData.sh yield_AN_sr_sys_0830_JES_Up no JES_Up no no
source process_noData.sh yield_AN_sr_sys_0830_JES_Down no JES_Down no no
source process_noData.sh yield_AN_sr_sys_0830_JMR_Up no JMR_Up no no
source process_noData.sh yield_AN_sr_sys_0830_JMR_Down no JMR_Down no no
source process_noData.sh yield_AN_sr_sys_0830_JMS_Up no JMS_Up no no
source process_noData.sh yield_AN_sr_sys_0830_JMS_Down no JMS_Down no no
source process_noData.sh yield_AN_sr_sys_0830_nominal yes nominal no yes
cp hists/yield_AN_sr_sys_0830_nominal_nosys/2016/data.root hists/yield_AN_sr_sys_0830_nominal/2016/
cp hists/yield_AN_sr_sys_0830_nominal_nosys/2017/data.root hists/yield_AN_sr_sys_0830_nominal/2017/
cp hists/yield_AN_sr_sys_0830_nominal_nosys/2018/data.root hists/yield_AN_sr_sys_0830_nominal/2018/
cp hists/yield_AN_sr_sys_0830_nominal_nosys/combine/data.root hists/yield_AN_sr_sys_0830_nominal/combine/ 
#source process_noData.sh yield_AN_sr_sys_0830_trig_nominal yes nominal yes no
