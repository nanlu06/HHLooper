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
    #["yield_AN_sr_sys_1111_nominal_nosys","no","nominal","no","no"],
    #["yield_AN_sr_sys_1111_JER_Up","no","JER_Up","no","no"],
    #["yield_AN_sr_sys_1111_JER_Down","no","JER_Down","no","no"],
    #["yield_AN_sr_sys_1111_JES_Up","no","JES_Up","no","no"],
    #["yield_AN_sr_sys_1111_JES_Down","no","JES_Down","no","no"],
    #["yield_AN_sr_sys_1111_JMR_Up","no","JMR_Up","no","no"],
    #["yield_AN_sr_sys_1111_JMR_Down","no","JMR_Down","no","no"],
    #["yield_AN_sr_sys_1111_JMS_Up","no","JMS_Up","no","no"],
    #["yield_AN_sr_sys_1111_JMS_Down","no","JMS_Down","no","no"],
    #["yield_AN_sr_sys_1111_nominal","yes","nominal","no","yes"],
    ["yield_AN_sr_sys_1111_trig_nominal","yes","nominal","yes","no"],
]
inputBase = "/storage/af/user/nlu/work/HH/ntuples/20210712_regression_v2/option5/combined/BDT/"

for data in datasets:
    for year in ["2016","2017","2018"]:
    #for year in ["2016"]:
        for tag in Tags:
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
            inputDir = inputBase+year+"/"+data
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
                print line.rstrip()
                        
            submitCommand = "condor_submit "+condorSubmit
            print submitCommand
            os.system(submitCommand)

for year in ["2016","2017","2018"]:
    condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/"+tag[0]+"_data"
    os.system(condorout)
    condorout = "mkdir /storage/af/user/idutta/HH/CMSSW_9_4_2/src/V2/HHLooper/condor/condor_output/"+tag[0]+"_data/"+year
    os.system(condorout)
    condorSubmit = "condor_submit/submitCondor_data_yield_AN_sr_sys_1111_nominal_nosys_"+year
    jobName      = str(date.today())+"_data_yield_AN_sr_sys_1111_nominal_nosys_"+year
    outHistFile = "data.root"
    isData       =  "1"
    inputDir = inputBase+year+"/data"
    shutil.copyfile("proto_condor_submit",condorSubmit)
    for line in fileinput.FileInput(condorSubmit, inplace=1):
        line=line.replace("JOBNAME", jobName)
        line=line.replace("INPUTDIR",inputDir)
        line=line.replace("OUTFILENAME",outHistFile)
        line=line.replace("TAG","yield_AN_sr_sys_1111_nominal_nosys")
        line=line.replace("ISDATA",isData)
        line=line.replace("DOSYST","no")
        line=line.replace("DOSHAPESYST","nominal")
        line=line.replace("DOTRIGSYST","no")
        line=line.replace("DOPNETSFSYST","no")
        line=line.replace("YEAR",year)
        print line.rstrip()
                        
    submitCommand = "condor_submit "+condorSubmit
    print submitCommand
    os.system(submitCommand)
    
