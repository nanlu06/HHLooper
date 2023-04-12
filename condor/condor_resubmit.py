import os,sys,re,fileinput,string,shutil
from datetime import date
import subprocess

#com1 = 'grep -lR "Error" condor_output/condor_logs/*.err'
#com2 = 'grep -lR "cannot" condor_output/condor_logs/*.err'
com3 = 'grep -lR "error" condor_output/condor_logs/*.err'
#files = subprocess.check_output(com1,shell=True)
#files_1 = subprocess.check_output(com2,shell=True)
files_2 = subprocess.check_output(com3,shell=True)
#print(files.split(".err\n"))

#files = files.split(".err\n")
#files_1 = files_1.split(".err\n")
files_2 = files_2.split(".err\n")

print(files_2)
'''
for f in files:
    #print(f)
    if (f not in files_1) and ( f not in files_2):
        #print(f)
        f_n = f.replace("condor_output/condor_logs/condor_2021-11-07_","")
        print(f_n)
        os.system("condor_submit condor_submit/submitCondor_"+f_n)

for f in files_1:
    if ( f not in files_2):
        f_n = f.replace("condor_output/condor_logs/condor_2021-11-07_","")
        print(f_n)
        os.system("condor_submit condor_submit/submitCondor_"+f_n)
'''
for f in files_2:
    f_n = f.replace("condor_output/condor_logs/condor_2021-11-07_","")
    print(f_n)
    os.system("condor_submit condor_submit/submitCondor_"+f_n)
