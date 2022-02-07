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
#"others"
"vjets",
"vv"
]
Tags =[
    ["yield_AN_sr_sys_010222_nominal_nosys","no","nominal","no","no"],
    ["yield_AN_sr_sys_010222_JER_Up","no","JER_Up","no","no"],
    ["yield_AN_sr_sys_010222_JER_Down","no","JER_Down","no","no"],
    #["yield_AN_sr_sys_010222_JES_Up","no","JES_Up","no","no"],
    #["yield_AN_sr_sys_010222_JES_Down","no","JES_Down","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_Abs","no","JESUp_Abs","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_Abs","no","JESDown_Abs","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_Abs_2016","no","JESUp_Abs_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_Abs_2016","no","JESDown_Abs_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_Abs_2017","no","JESUp_Abs_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_Abs_2017","no","JESDown_Abs_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_Abs_2018","no","JESUp_Abs_2018","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_Abs_2018","no","JESDown_Abs_2018","no","no"],
    
    ["yield_AN_sr_sys_010222_JESUp_BBEC1","no","JESUp_BBEC1","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_BBEC1","no","JESDown_BBEC1","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_BBEC1_2016","no","JESUp_BBEC1_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_BBEC1_2016","no","JESDown_BBEC1_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_BBEC1_2017","no","JESUp_BBEC1_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_BBEC1_2017","no","JESDown_BBEC1_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_BBEC1_2018","no","JESUp_BBEC1_2018","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_BBEC1_2018","no","JESDown_BBEC1_2018","no","no"],

    ["yield_AN_sr_sys_010222_JESUp_EC2","no","JESUp_EC2","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_EC2","no","JESDown_EC2","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_EC2_2016","no","JESUp_EC2_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_EC2_2016","no","JESDown_EC2_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_EC2_2017","no","JESUp_EC2_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_EC2_2017","no","JESDown_EC2_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_EC2_2018","no","JESUp_EC2_2018","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_EC2_2018","no","JESDown_EC2_2018","no","no"],
    
    ["yield_AN_sr_sys_010222_JESUp_FlavQCD","no","JESUp_FlavQCD","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_FlavQCD","no","JESDown_FlavQCD","no","no"],

    ["yield_AN_sr_sys_010222_JESUp_HF","no","JESUp_HF","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_HF","no","JESDown_HF","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_HF_2016","no","JESUp_HF_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_HF_2016","no","JESDown_HF_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_HF_2017","no","JESUp_HF_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_HF_2017","no","JESDown_HF_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_HF_2018","no","JESUp_HF_2018","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_HF_2018","no","JESDown_HF_2018","no","no"],

    ["yield_AN_sr_sys_010222_JESUp_RelBal","no","JESUp_RelBal","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_RelBal","no","JESDown_RelBal","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_RelSample_2016","no","JESUp_RelSample_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_RelSample_2016","no","JESDown_RelSample_2016","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_RelSample_2017","no","JESUp_RelSample_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_RelSample_2017","no","JESDown_RelSample_2017","no","no"],
    ["yield_AN_sr_sys_010222_JESUp_RelSample_2018","no","JESUp_RelSample_2018","no","no"],
    ["yield_AN_sr_sys_010222_JESDown_RelSample_2018","no","JESDown_RelSample_2018","no","no"],

    ["yield_AN_sr_sys_010222_JMR_Up","no","JMR_Up","no","no"],
    ["yield_AN_sr_sys_010222_JMR_Down","no","JMR_Down","no","no"],
    ["yield_AN_sr_sys_010222_JMS_Up","no","JMS_Up","no","no"],
    ["yield_AN_sr_sys_010222_JMS_Down","no","JMS_Down","no","no"],
    ["yield_AN_sr_sys_010222_nominal","yes","nominal","no","yes"],
    ["yield_AN_sr_sys_010222_trig_nominal","yes","nominal","yes","no"],
]

inputBase = "/storage/af/user/idutta/work/HH/ntuple/20211209_regression/option5/combined/BDT/"

for data in datasets:
    for year in ["2016","2017","2018"]:
    #for year in ["2016"]:
        for tag in Tags:
            if "vjets" in data or "vv" in data:
                condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/"+tag[0]+"_others"
                os.system(condorout)
                condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/"+tag[0]+"_others/"+year
                os.system(condorout)
            else:
                condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/"+tag[0]+"_"+data
                os.system(condorout)
                condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/"+tag[0]+"_"+data+"/"+year
                os.system(condorout)
            condorSubmit = "condor_submit/submitCondor_"+data+"_"+tag[0]+"_"+year
            jobName      = str(date.today())+"_"+data+"_"+tag[0]+"_"+year
            outHistFile = data+".root"
            if("data" in data):
                isData       =  "1"
            else:
                isData       =  "0"
            if "vjets" in data or "vv" in data:
                if ("vjets" in data):
                    inputDir = inputBase+year+"/others/VJets"
                else:
                    inputDir = inputBase+year+"/others/VV"
                proc="others/"
            else:
                inputDir = inputBase+year+"/"+data
                proc = data
            print(tag,data)
            print(inputDir) 
            
            shutil.copyfile("proto_condor_submit",condorSubmit)
            for line in fileinput.FileInput(condorSubmit, inplace=1):
                line=line.replace("JOBNAME", jobName)
                line=line.replace("INPUTDIR",inputDir)
                line=line.replace("OUTFILENAME",outHistFile)
                line=line.replace("TAG",tag[0])
                line=line.replace("ISDATA",isData)
                line=line.replace("DOSYST",tag[1])
                line=line.replace("DOSHAPESYST",tag[2])
                line=line.replace("DOTRIGSYST",tag[3])
                line=line.replace("DOPNETSFSYST",tag[4])
                line=line.replace("YEAR",year)
                line=line.replace("PROC",proc)
                print line.rstrip()
                        
            submitCommand = "condor_submit "+condorSubmit
            print submitCommand
            os.system(submitCommand)
inputBase_data = "/storage/af/user/nlu/work/HH/ntuples/20210712_regression_v2/option5/combined/BDT/"
for year in ["2016","2017","2018"]:
    condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/yield_AN_sr_sys_010222_nominal_nosys_data"
    os.system(condorout)
    condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/yield_AN_sr_sys_010222_nominal_nosys_data/"+year
    os.system(condorout)
    condorSubmit = "condor_submit/submitCondor_data_yield_AN_sr_sys_010222_nominal_nosys_"+year
    jobName      = str(date.today())+"_data_yield_AN_sr_sys_010222_nominal_nosys_"+year
    outHistFile = "data.root"
    isData       =  "1"
    inputDir = inputBase_data+year+"/data"
    proc ='data'
    shutil.copyfile("proto_condor_submit",condorSubmit)
    for line in fileinput.FileInput(condorSubmit, inplace=1):
        line=line.replace("JOBNAME", jobName)
        line=line.replace("INPUTDIR",inputDir)
        line=line.replace("OUTFILENAME",outHistFile)
        line=line.replace("TAG","yield_AN_sr_sys_010222_nominal_nosys")
        line=line.replace("ISDATA",isData)
        line=line.replace("DOSYST","no")
        line=line.replace("DOSHAPESYST","nominal")
        line=line.replace("DOTRIGSYST","no")
        line=line.replace("DOPNETSFSYST","no")
        line=line.replace("YEAR",year)
        line=line.replace("PROC",proc)
        print line.rstrip()
                        
    submitCommand = "condor_submit "+condorSubmit
    print submitCommand
    os.system(submitCommand)
