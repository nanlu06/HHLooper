import os,sys,re,fileinput,string,shutil
from datetime import date
##             Dataset        Name   
datasets = ["qcd",
"ttbar",
"HHSM",
"HHc0",
"HHc1",
"HHc2p45",
"HHc5",
"HHVBFSM",
"HHVBF0p511",
"HHVBF1p511",
"HHVBF110",
"HHVBF112",
"HHVBF121",
"HHVBF101",
"VH",
"Higgs",
"ttH",
"others"
]
Tags =[
    ["yield_AN_sr_sys_0830_nominal_nosys","no","nominal","no","no"],
    ["yield_AN_sr_sys_0830_JER_Up","no","JER_Up","no","no"],
    ["yield_AN_sr_sys_0830_JER_Down","no","JER_Down","no","no"],
    ["yield_AN_sr_sys_0830_JES_Up","no","JES_Up","no","no"],
    ["yield_AN_sr_sys_0830_JES_Down","no","JES_Down","no","no"],
    ["yield_AN_sr_sys_0830_JMR_Up","no","JMR_Up","no","no"],
    ["yield_AN_sr_sys_0830_JMR_Down","no","JMR_Down","no","no"],
    ["yield_AN_sr_sys_0830_JMS_Up","no","JMS_Up","no","no"],
    ["yield_AN_sr_sys_0830_JMS_Down","no","JMS_Down","no","no"],
    ["yield_AN_sr_sys_0830_nominal","yes","nominal","no","yes"],
    #["yield_AN_sr_sys_0830_trig_nominal","yes","nominal","yes","no"],
]
inputBase = "/storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/" 
condordir = inputBase+"condor/condor_output/condor_logs"
os.chdir(condordir)
os.system("mkdir -p hists/yield_AN_sr_sys_0830_nominal_nosys/combine/" )
os.system("hadd -k -f hists/yield_AN_sr_sys_0830_nominal_nosys/combine/data.root hists/yield_AN_sr_sys_0830_nominal_nosys/2016/data.root hists/yield_AN_sr_sys_0830_nominal_nosys/2017/data.root hists/yield_AN_sr_sys_0830_nominal_nosys/2018/data.root")
for tag in Tags:
    mkdircomm= "mkdir -p hists/"+tag[0]+"/combine/"
    os.system(mkdircomm)
    os.system("cp hists/yield_AN_sr_sys_0830_nominal_nosys/combine/data.root hists/"+tag[0]+"/combine/data.root")
    os.system("cp hists/yield_AN_sr_sys_0830_nominal_nosys/2016/data.root hists/"+tag[0]+"/2016/data.root") ## copy the nominal data file here
    os.system("cp hists/yield_AN_sr_sys_0830_nominal_nosys/2017/data.root hists/"+tag[0]+"/2017/data.root") ## copy the nominal data file here
    os.system("cp hists/yield_AN_sr_sys_0830_nominal_nosys/2018/data.root hists/"+tag[0]+"/2018/data.root") ## copy the nominal data file here       
    for data in datasets:
        haddCommand = "hadd -k -f hists/"+tag[0]+"/combine/"+data+".root hists/"+tag[0]+"/2016/"+data+".root hists/"+tag[0]+"/2017/"+data+".root hists/"+tag[0]+"/2018/"+data+".root"
        print haddCommand
        os.system(haddCommand)
    os.system("hadd -k -f hists/"+tag[0]+"/combine/QCDggHVBF.root hists/"+tag[0]+"/combine/qcd.root hists/"+tag[0]+"/combine/Higgs.root")
    os.system("hadd -k -f hists/"+tag[0]+"/combine/bkg.root hists/"+tag[0]+"/combine/qcd.root hists/"+tag[0]+"/combine/ttbar.root  hists/"+tag[0]+"/combine/VH.root hists/"+tag[0]+"/combine/Higgs.root hists/"+tag[0]+"/combine/ttH.root hists/"+tag[0]+"/combine/others.root")

os.system("mv hists "+inputBase)
