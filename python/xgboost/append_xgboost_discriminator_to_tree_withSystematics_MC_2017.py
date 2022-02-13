#!/usr/bin/env python

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt
from xgboost import XGBClassifier
from xgboost import plot_tree
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_curve
#import root_pandas as rp
import numpy as np
import pandas as pd
import math
import pickle as pickle
import ROOT as root
import uproot3 as uproot

#################
##Preliminaries
#################

#USAGE
#python append_xgboost_discriminator_to_tree.py <folder_name> <input_file_name> <year> <bkg_type> <bdt_type>

root.gROOT.Reset()

print(sys.argv)
input_file  = sys.argv[1]
output_file = sys.argv[2]
_year       = '2017'
_bkg_type   = 'ttbar'##options are <qcd> and <ttbar>
_bdt_type   = 'basic0'

if len(sys.argv) >= 4:
    _year  = sys.argv[3]#<2016> <2017> <2018>
if len(sys.argv) >= 5:
    _bkg_type  = sys.argv[4]##options are <qcd>, <ttbar>, and <qcd_and_ttbar>
if len(sys.argv) >= 6:
    _bdt_type  = sys.argv[5]#<basic0>, <basic1>, <basic2>, <enhanced>, <enhanced_v2>

print(input_file, output_file, _year, _bkg_type, _bdt_type)

FileName = input_file
File = root.TFile(FileName,'READ')
#inputTree = File.Get('tree')
inputTree = File.Get('Events') #for new ntuples from Cristina
Nevents = File.Get('NEvents')



test_name = 'ReadingXgBoostModel'

_model_name = 'models/model_xgboost_training_weights_' + _bkg_type + '_' + _year + '_bdt_' + _bdt_type +'.pkl'

if _bdt_type == 'basic0':
    variables =   [
    ['fatJet1MassSD', 'j1_mass_sd', '$M_{j1}$ (GeV)', 40,  0.,   5000.],
    ['fatJet1Pt', 'j1_pt', '$p_{T}^{j1}$ (GeV)', 40,  0.,   5000.],

    ['fatJet2Pt', 'j2_pt', '$p_{T}^{j2}$ (GeV)', 40,  0.,   5000.],

    ]
elif _bdt_type == 'enhanced_v8p2':
    variables =   [
    ['hh_pt', 'hh_pt', '$p_{T}^{HH}$ (GeV)', 40, 0, 5000],
    ['hh_eta', 'hh_eta', '$\eta^{HH}$', 40, -5.0, 5.0],
    #['hh_phi', 'hh_phi', '$\phi^{HH}$', 40, -3.2, 3.2],
    ['hh_mass', 'hh_mass', '$m_{HH}$ (GeV)', 40, 0, 1500],
    #['MET', 'MET', '$MET$ (GeV)', 60, 0, 600],
    ['met', 'met', '$MET$ (GeV)', 60, 0, 600], #for new ntuples


    #                ['FatJet1_area', 'j1_area', 'fat j1 area', 40, 1.85, 2.15],
    #                ['FatJet2_area', 'j2_area', 'fat j2 area', 40,  1.85, 2.15],
    #['fatJet1MassSD', 'j1_m', 'j1 soft drop mass (GeV)', 40,  0.,   200.],
    #                ['fatJet1Mass', 'j2_m', 'j2 soft drop mass (GeV)', 40,  0.,   200.],

    #                ['fatJet1HasBJetCSVLoose', 'j1_CSVLoose', 'j1 DDB tagger', 3,  0,  3],
    # ['fatJet2HasBJetCSVLoose', 'j2_CSVLoose', 'j2 DDB tagger', 40,  0.78,  1.0],

    ['fatJet1Tau3OverTau2', 'fatJet1Tau3OverTau2', 'fatJet1Tau3OverTau2', 50,  0.0,  1.0],
    ['fatJet2Tau3OverTau2', 'fatJet2Tau3OverTau2', 'fatJet2Tau3OverTau2', 50,  0.0,  1.0],
    ['fatJet1MassSD', 'j1_mass_sd', '$M_{j1}$ (GeV)', 40,  0.,   5000.],
    ['fatJet1Pt', 'j1_pt', '$p_{T}^{j1}$ (GeV)', 40,  0.,   5000.],
    ['fatJet1Eta', 'j1_eta', '$\eta^{j1}$', 40,  -2.5,  2.5],
    #['fatJet1Phi', 'j1_phi', '$\phi^{j1}$', 40,  -3.2,   3.2],
    ['fatJet1PNetXbb', 'fatJet1PNetXbb', 'fatJet1PNetXbb', 40,  -100,   100],
    ['fatJet1PNetQCDb', 'fatJet1PNetQCDb', 'fatJet1PNetQCDb', 40,  -100,   100],
    ['fatJet1PNetQCDbb', 'fatJet1PNetQCDbb', 'fatJet1PNetQCDbb', 40,  -100,   100],
    ['fatJet1PNetQCDothers', 'fatJet1PNetQCDothers', 'fatJet1PNetQCDothers', 40,  -100,   100],
    ['fatJet2Pt', 'j2_pt', '$p_{T}^{j2}$ (GeV)', 40,  0.,   500.],
    #['fatJet2PNetXbb', 'fatJet2PNetXbb', 'fatJet2PNetXbb', 40,  -100,   100],
    #['fatJet2PNetQCDb', 'fatJet2PNetQCDb', 'fatJet2PNetQCDb', 40,  -100,   100],
    #['fatJet2PNetQCDbb', 'fatJet2PNetQCDbb', 'fatJet2PNetQCDbb', 40,  -100,   100],
    #['fatJet2PNetQCDothers', 'fatJet2PNetQCDothers', 'fatJet2PNetQCDothers', 40,  -100,   100],
    #['deltaEta_j1j2', 'dEta_j1j2', '$\Delta\eta(j_{1}, j_{2})$', 40,  0.,   5.],
    #['deltaPhi_j1j2', 'dPhi_j1j2', '$\Delta\phi(j_{1}, j_{2})$', 40,  2.,   4.5],
    #['deltaR_j1j2', 'dR_j1j2', '$\Delta R(j_{1}, j_{2})$', 40,  0.,   5.],
    ['fatJet1PtOverMHH', 'ptj1Omhh', '$p_{T}^{j1}/m_{HH}$', 40,   0.,   1.],
    ['fatJet2PtOverMHH', 'ptj2Omhh', '$p_{T}^{j2}/m_{HH}$', 40,  0.,  0.7],

    #                ['ptj1_over_mj1', 'ptj1Omj1', '$p_{T}^{j1}/m_{j1}$', 40,  0.,   10.],
    #                ['ptj2_over_mj2', 'ptj2Omj2', '$p_{T}^{j2}/m_{j2}$', 40,  0.5,  10.],

    ['ptj2_over_ptj1', 'ptj2Optj1', '$p_{T}^{j2}/p_{T}^{j1}$', 40,  0.5,  1.],
    #['mj2_over_mj1', 'mj2Omj1', '$m^{j2}/m^{j1}$', 40,  0.0,  1.5],
    ]
elif _bdt_type == 'enhanced_v24':
    variables =   [
    ['hh_pt', 'hh_pt', '$p_{T}^{HH}$ (GeV)', 40, 0, 5000],
    ['hh_eta', 'hh_eta', '$\eta^{HH}$', 40, -5.0, 5.0],
    #['hh_phi', 'hh_phi', '$\phi^{HH}$', 40, -3.2, 3.2],
    ['hh_mass', 'hh_mass', '$m_{HH}$ (GeV)', 40, 0, 1500],
    #['MET', 'MET', '$MET$ (GeV)', 60, 0, 600],
    ['met', 'met', '$MET$ (GeV)', 60, 0, 600], #for new ntuples
    #                ['FatJet1_area', 'j1_area', 'fat j1 area', 40, 1.85, 2.15],
    #                ['FatJet2_area', 'j2_area', 'fat j2 area', 40,  1.85, 2.15],
    #['fatJet1MassSD', 'j1_m', 'j1 soft drop mass (GeV)', 40,  0.,   200.],
    #                ['fatJet1Mass', 'j2_m', 'j2 soft drop mass (GeV)', 40,  0.,   200.],

    #                ['fatJet1HasBJetCSVLoose', 'j1_CSVLoose', 'j1 DDB tagger', 3,  0,  3],
    # ['fatJet2HasBJetCSVLoose', 'j2_CSVLoose', 'j2 DDB tagger', 40,  0.78,  1.0],

    #['fatJet1Tau3OverTau2', 'fatJet1Tau3OverTau2', 'fatJet1Tau3OverTau2', 50,  0.0,  1.0],
    #['fatJet2Tau3OverTau2', 'fatJet2Tau3OverTau2', 'fatJet2Tau3OverTau2', 50,  0.0,  1.0],
    ['fatJet1MassSD', 'j1_mass_sd', '$M_{j1}$ (GeV)', 40,  0.,   5000.],
    ['fatJet1Pt', 'j1_pt', '$p_{T}^{j1}$ (GeV)', 40,  0.,   5000.],
    ['fatJet1Eta', 'j1_eta', '$\eta^{j1}$', 40,  -2.5,  2.5],
    #['fatJet1Phi', 'j1_phi', '$\phi^{j1}$', 40,  -3.2,   3.2],
    #['fatJet1PNetXbb', 'fatJet1PNetXbb', 'fatJet1PNetXbb', 40,  -100,   100],
    #['fatJet1PNetQCDb', 'fatJet1PNetQCDb', 'fatJet1PNetQCDb', 40,  -100,   100],
    #['fatJet1PNetQCDbb', 'fatJet1PNetQCDbb', 'fatJet1PNetQCDbb', 40,  -100,   100],
    #['fatJet1PNetQCDothers', 'fatJet1PNetQCDothers', 'fatJet1PNetQCDothers', 40,  -100,   100],
    ['fatJet2Pt', 'j2_pt', '$p_{T}^{j2}$ (GeV)', 40,  0.,   500.],
    #['fatJet2PNetXbb', 'fatJet2PNetXbb', 'fatJet2PNetXbb', 40,  -100,   100],
    #['fatJet2PNetQCDb', 'fatJet2PNetQCDb', 'fatJet2PNetQCDb', 40,  -100,   100],
    #['fatJet2PNetQCDbb', 'fatJet2PNetQCDbb', 'fatJet2PNetQCDbb', 40,  -100,   100],
    #['fatJet2PNetQCDothers', 'fatJet2PNetQCDothers', 'fatJet2PNetQCDothers', 40,  -100,   100],
    #['deltaEta_j1j2', 'dEta_j1j2', '$\Delta\eta(j_{1}, j_{2})$', 40,  0.,   5.],
    #['deltaPhi_j1j2', 'dPhi_j1j2', '$\Delta\phi(j_{1}, j_{2})$', 40,  2.,   4.5],
    #['deltaR_j1j2', 'dR_j1j2', '$\Delta R(j_{1}, j_{2})$', 40,  0.,   5.],
    ['fatJet1PtOverMHH', 'ptj1Omhh', '$p_{T}^{j1}/m_{HH}$', 40,   0.,   1.],
    ['fatJet2PtOverMHH', 'ptj2Omhh', '$p_{T}^{j2}/m_{HH}$', 40,  0.,  0.7],

    #                ['ptj1_over_mj1', 'ptj1Omj1', '$p_{T}^{j1}/m_{j1}$', 40,  0.,   10.],
    #                ['ptj2_over_mj2', 'ptj2Omj2', '$p_{T}^{j2}/m_{j2}$', 40,  0.5,  10.],

    ['ptj2_over_ptj1', 'ptj2Optj1', '$p_{T}^{j2}/p_{T}^{j1}$', 40,  0.5,  1.],
    #['mj2_over_mj1', 'mj2Omj1', '$m^{j2}/m^{j1}$', 40,  0.0,  1.5],
    ]
elif _bdt_type == 'mass_sculpting_control':
    variables =   [
    ['hh_pt', 'hh_pt', '$p_{T}^{HH}$ (GeV)', 40, 0, 5000],
    ['hh_eta', 'hh_eta', '$\eta^{HH}$', 40, -5.0, 5.0],
    ['hh_phi', 'hh_phi', '$\phi^{HH}$', 40, -3.2, 3.2],
    ['hh_mass', 'hh_mass', '$m_{HH}$ (GeV)', 40, 0, 1500],
    #['MET', 'MET', '$MET$ (GeV)', 60, 0, 600],
    ['met', 'met', '$MET$ (GeV)', 60, 0, 600], #for new ntuples


    #                ['FatJet1_area', 'j1_area', 'fat j1 area', 40, 1.85, 2.15],
    #                ['FatJet2_area', 'j2_area', 'fat j2 area', 40,  1.85, 2.15],
    #['fatJet1MassSD', 'j1_m', 'j1 soft drop mass (GeV)', 40,  0.,   200.],
    #                ['fatJet1Mass', 'j2_m', 'j2 soft drop mass (GeV)', 40,  0.,   200.],

    #                ['fatJet1HasBJetCSVLoose', 'j1_CSVLoose', 'j1 DDB tagger', 3,  0,  3],
    # ['fatJet2HasBJetCSVLoose', 'j2_CSVLoose', 'j2 DDB tagger', 40,  0.78,  1.0],

    ['fatJet1Tau3OverTau2', 'fatJet1Tau3OverTau2', 'fatJet1Tau3OverTau2', 50,  0.0,  1.0],
    ['fatJet2Tau3OverTau2', 'fatJet2Tau3OverTau2', 'fatJet2Tau3OverTau2', 50,  0.0,  1.0],
    ['fatJet1MassSD', 'j1_mass_sd', '$M_{j1}$ (GeV)', 40,  0.,   5000.],
    ['fatJet1Pt', 'j1_pt', '$p_{T}^{j1}$ (GeV)', 40,  0.,   5000.],
    ['fatJet1Eta', 'j1_eta', '$\eta^{j1}$', 40,  -2.5,  2.5],
    ['fatJet1Phi', 'j1_phi', '$\phi^{j1}$', 40,  -3.2,   3.2],
    ['fatJet1PNetXbb', 'fatJet1PNetXbb', 'fatJet1PNetXbb', 40,  -100,   100],
    ['fatJet1PNetQCDb', 'fatJet1PNetQCDb', 'fatJet1PNetQCDb', 40,  -100,   100],
    ['fatJet1PNetQCDbb', 'fatJet1PNetQCDbb', 'fatJet1PNetQCDbb', 40,  -100,   100],
    ['fatJet1PNetQCDothers', 'fatJet1PNetQCDothers', 'fatJet1PNetQCDothers', 40,  -100,   100],
    ['fatJet2Pt', 'j2_pt', '$p_{T}^{j2}$ (GeV)', 40,  0.,   500.],
    ['fatJet2MassSD', 'j2_mass_sd', '$M_{j2}$ (GeV)', 40,  0.,   5000.],
    #['fatJet2PNetXbb', 'fatJet2PNetXbb', 'fatJet2PNetXbb', 40,  -100,   100],
    #['fatJet2PNetQCDb', 'fatJet2PNetQCDb', 'fatJet2PNetQCDb', 40,  -100,   100],
    #['fatJet2PNetQCDbb', 'fatJet2PNetQCDbb', 'fatJet2PNetQCDbb', 40,  -100,   100],
    #['fatJet2PNetQCDothers', 'fatJet2PNetQCDothers', 'fatJet2PNetQCDothers', 40,  -100,   100],
    ['deltaEta_j1j2', 'dEta_j1j2', '$\Delta\eta(j_{1}, j_{2})$', 40,  0.,   5.],
    ['deltaPhi_j1j2', 'dPhi_j1j2', '$\Delta\phi(j_{1}, j_{2})$', 40,  2.,   4.5],
    ['deltaR_j1j2', 'dR_j1j2', '$\Delta R(j_{1}, j_{2})$', 40,  0.,   5.],
    ['fatJet1PtOverMHH', 'ptj1Omhh', '$p_{T}^{j1}/m_{HH}$', 40,   0.,   1.],
    ['fatJet2PtOverMHH', 'ptj2Omhh', '$p_{T}^{j2}/m_{HH}$', 40,  0.,  0.7],

    #                ['ptj1_over_mj1', 'ptj1Omj1', '$p_{T}^{j1}/m_{j1}$', 40,  0.,   10.],
    #                ['ptj2_over_mj2', 'ptj2Omj2', '$p_{T}^{j2}/m_{j2}$', 40,  0.5,  10.],

    ['ptj2_over_ptj1', 'ptj2Optj1', '$p_{T}^{j2}/p_{T}^{j1}$', 40,  0.5,  1.],
    #['mj2_over_mj1', 'mj2Omj1', '$m^{j2}/m^{j1}$', 40,  0.0,  1.5],

    #['totalWeight', 'totalWeight', 'totalWeight', 100, -1000, 1000]
    ]



print(len(variables))
my_variables = []
for var in variables:
    print (var[0])
    my_variables.append(var[0])
#df = uproot.open(FileName)['tree'].pandas.df([row[0] for row in variables], flatten=False)
df = uproot.open(FileName)['Events'].pandas.df([row[0] for row in variables], flatten=False) #for new ntuples


#getting a numpy array from two pandas data frames
x_test = df.values
#creating numpy array for target variables
y_test = np.zeros(len(df))


############################
# get model from file
############################
with open(_model_name,'rb') as pkl_file:
    model = pickle.load(pkl_file)

# make predictions for test data
y_pred = model.predict_proba(x_test)[:, 1]
print ("y_pred", y_pred)
predictions = [round(value) for value in y_pred]

#print y_pred
##########################################################
# make histogram of discriminator value for signal and bkg
##########################################################
y_frame = pd.DataFrame({'truth':y_test, 'disc':y_pred})
print ('y_frame',y_frame)
disc    = y_frame[y_frame['truth'] == 0]['disc'].values
#disc    = y_frame[y_frame['truth'] == 1]['disc'].values
plt.figure()
plt.hist(disc, bins=50, alpha=0.3)
plt.savefig('mydiscriminator_' + test_name + '.png')
print ("disc_bkg: ", disc)
#print y_pred


#create output file
outputFilename = output_file
print("output file: " + outputFilename + "\n")
outputFile = root.TFile( outputFilename,  "RECREATE")
outputFile.cd()
outputTree = inputTree.CloneTree(0);

nEntries = inputTree.GetEntries()
print ("nEntries = ", nEntries)
_disc_var_name = 'disc_' + _bkg_type + '_' + _year + '_' + _bdt_type 

root.gROOT.ProcessLine("struct MyStruct{float disc;};")
from ROOT import MyStruct

nSamplings = 100
random = root.TRandom3(1000)

# Create branches in the tree
my_s = MyStruct()
my_s_JESUp = MyStruct()
my_s_JESDown = MyStruct()
my_s_JESUp_Abs = MyStruct()
my_s_JESDown_Abs = MyStruct()
my_s_JESUp_Abs_2017 = MyStruct()
my_s_JESDown_Abs_2017 = MyStruct()
my_s_JESUp_BBEC1 = MyStruct()
my_s_JESDown_BBEC1 = MyStruct()
my_s_JESUp_BBEC1_2017 = MyStruct()
my_s_JESDown_BBEC1_2017 = MyStruct()
my_s_JESUp_EC2 = MyStruct()
my_s_JESDown_EC2 = MyStruct()
my_s_JESUp_EC2_2017 = MyStruct()
my_s_JESDown_EC2_2017 = MyStruct()
my_s_JESUp_HF = MyStruct()
my_s_JESDown_HF = MyStruct()
my_s_JESUp_HF_2017 = MyStruct()
my_s_JESDown_HF_2017 = MyStruct()
my_s_JESUp_FlavQCD = MyStruct()
my_s_JESDown_FlavQCD = MyStruct()
my_s_JESUp_RelBal = MyStruct()
my_s_JESDown_RelBal = MyStruct()
my_s_JESUp_RelSample_2017 = MyStruct()
my_s_JESDown_RelSample_2017 = MyStruct()
my_s_JMSUp = MyStruct()
my_s_JMSDown = MyStruct()
my_s_JMRUp = MyStruct()
my_s_JMRDown = MyStruct()

#print(my_s)
#print (root.addressof(my_s))
#print (root.addressof(my_s,'disc'))
#root.addressof(my_s,'disc')

branch_disc = outputTree.Branch(_disc_var_name,root.addressof(my_s,'disc'), _disc_var_name+'/F');
branch_disc_JESUp = outputTree.Branch(_disc_var_name+"_JESUp",root.addressof(my_s_JESUp,'disc'), _disc_var_name+'_JESUp/F');
branch_disc_JESDown = outputTree.Branch(_disc_var_name+"_JESDown",root.addressof(my_s_JESDown,'disc'), _disc_var_name+'_JESDown/F');
branch_disc_JESUp_Abs = outputTree.Branch(_disc_var_name+"_JESUp_Abs",root.addressof(my_s_JESUp_Abs,'disc'), _disc_var_name+'_JESUp_Abs/F');
branch_disc_JESDown_Abs = outputTree.Branch(_disc_var_name+"_JESDown_Abs",root.addressof(my_s_JESDown_Abs,'disc'), _disc_var_name+'_JESDown_Abs/F');
branch_disc_JESUp_Abs_2017 = outputTree.Branch(_disc_var_name+"_JESUp_Abs_2017",root.addressof(my_s_JESUp_Abs_2017,'disc'), _disc_var_name+'_JESUp_Abs_2017/F');
branch_disc_JESDown_Abs_2017 = outputTree.Branch(_disc_var_name+"_JESDown_Abs_2017",root.addressof(my_s_JESDown_Abs_2017,'disc'), _disc_var_name+'_JESDown_Abs_2017/F');
branch_disc_JESUp_BBEC1 = outputTree.Branch(_disc_var_name+"_JESUp_BBEC1",root.addressof(my_s_JESUp_BBEC1,'disc'), _disc_var_name+'_JESUp_BBEC1/F');
branch_disc_JESDown_BBEC1 = outputTree.Branch(_disc_var_name+"_JESDown_BBEC1",root.addressof(my_s_JESDown_BBEC1,'disc'), _disc_var_name+'_JESDown_BBEC1/F');
branch_disc_JESUp_BBEC1_2017 = outputTree.Branch(_disc_var_name+"_JESUp_BBEC1_2017",root.addressof(my_s_JESUp_BBEC1_2017,'disc'), _disc_var_name+'_JESUp_BBEC1_2017/F');
branch_disc_JESDown_BBEC1_2017 = outputTree.Branch(_disc_var_name+"_JESDown_BBEC1_2017",root.addressof(my_s_JESDown_BBEC1_2017,'disc'), _disc_var_name+'_JESDown_BBEC1_2017/F');
branch_disc_JESUp_EC2 = outputTree.Branch(_disc_var_name+"_JESUp_EC2",root.addressof(my_s_JESUp_EC2,'disc'), _disc_var_name+'_JESUp_EC2/F');
branch_disc_JESDown_EC2 = outputTree.Branch(_disc_var_name+"_JESDown_EC2",root.addressof(my_s_JESDown_EC2,'disc'), _disc_var_name+'_JESDown_EC2/F');
branch_disc_JESUp_EC2_2017 = outputTree.Branch(_disc_var_name+"_JESUp_EC2_2017",root.addressof(my_s_JESUp_EC2_2017,'disc'), _disc_var_name+'_JESUp_EC2_2017/F');
branch_disc_JESDown_EC2_2017 = outputTree.Branch(_disc_var_name+"_JESDown_EC2_2017",root.addressof(my_s_JESDown_EC2_2017,'disc'), _disc_var_name+'_JESDown_EC2_2017/F');
branch_disc_JESUp_HF = outputTree.Branch(_disc_var_name+"_JESUp_HF",root.addressof(my_s_JESUp_HF,'disc'), _disc_var_name+'_JESUp_HF/F');
branch_disc_JESDown_HF = outputTree.Branch(_disc_var_name+"_JESDown_HF",root.addressof(my_s_JESDown_HF,'disc'), _disc_var_name+'_JESDown_HF/F');
branch_disc_JESUp_HF_2017 = outputTree.Branch(_disc_var_name+"_JESUp_HF_2017",root.addressof(my_s_JESUp_HF_2017,'disc'), _disc_var_name+'_JESUp_HF_2017/F');
branch_disc_JESDown_HF_2017 = outputTree.Branch(_disc_var_name+"_JESDown_HF_2017",root.addressof(my_s_JESDown_HF_2017,'disc'), _disc_var_name+'_JESDown_HF_2017/F');
branch_disc_JESUp_FlavQCD = outputTree.Branch(_disc_var_name+"_JESUp_FlavQCD",root.addressof(my_s_JESUp_FlavQCD,'disc'), _disc_var_name+'_JESUp_FlavQCD/F');
branch_disc_JESDown_FlavQCD = outputTree.Branch(_disc_var_name+"_JESDown_FlavQCD",root.addressof(my_s_JESDown_FlavQCD,'disc'), _disc_var_name+'_JESDown_FlavQCD/F');
branch_disc_JESUp_RelBal = outputTree.Branch(_disc_var_name+"_JESUp_RelBal",root.addressof(my_s_JESUp_RelBal,'disc'), _disc_var_name+'_JESUp_RelBal/F');
branch_disc_JESDown_RelBal = outputTree.Branch(_disc_var_name+"_JESDown_RelBal",root.addressof(my_s_JESDown_RelBal,'disc'), _disc_var_name+'_JESDown_RelBal/F');
branch_disc_JESUp_RelSample_2017 = outputTree.Branch(_disc_var_name+"_JESUp_RelSample_2017",root.addressof(my_s_JESUp_RelSample_2017,'disc'), _disc_var_name+'_JESUp_RelSample_2017/F');
branch_disc_JESDown_RelSample_2017 = outputTree.Branch(_disc_var_name+"_JESDown_RelSample_2017",root.addressof(my_s_JESDown_RelSample_2017,'disc'), _disc_var_name+'_JESDown_RelSample_2017/F');
branch_disc_JMSUp = outputTree.Branch(_disc_var_name+"_JMSUp",root.addressof(my_s_JMSUp,'disc'), _disc_var_name+'_JMSUp/F');
branch_disc_JMSDown = outputTree.Branch(_disc_var_name+"_JMSDown",root.addressof(my_s_JMSDown,'disc'), _disc_var_name+'_JMSDown/F');
branch_disc_JMRUp = outputTree.Branch(_disc_var_name+"_JMRUp",root.addressof(my_s_JMRUp,'disc'), _disc_var_name+'_JMRUp/F');
branch_disc_JMRDown = outputTree.Branch(_disc_var_name+"_JMRDown",root.addressof(my_s_JMRDown,'disc'), _disc_var_name+'_JMRDown/F');


for i in range(nEntries):
    inputTree.GetEntry(i)
    if i%100==0:
        print ("Processing event nr. %i of %i" % (i,nEntries))
        
    tmpXList=[]
    ##############################
    #Only For enhanced_v8p2
    ##############################
    tmpXList.append(inputTree.hh_pt)
    tmpXList.append(inputTree.hh_eta)
    tmpXList.append(inputTree.hh_mass)
    #tmpXList.append(inputTree.MET)
    tmpXList.append(inputTree.met) #for new ntuples
    tmpXList.append(inputTree.fatJet1Tau3OverTau2)
    tmpXList.append(inputTree.fatJet2Tau3OverTau2)
    tmpXList.append(inputTree.fatJet1MassSD)
    tmpXList.append(inputTree.fatJet1Pt),
    tmpXList.append(inputTree.fatJet1Eta),
    tmpXList.append(inputTree.fatJet1PNetXbb)
    tmpXList.append(inputTree.fatJet1PNetQCDb)
    tmpXList.append(inputTree.fatJet1PNetQCDbb)
    tmpXList.append(inputTree.fatJet1PNetQCDothers)
    tmpXList.append(inputTree.fatJet2Pt)
    tmpXList.append(inputTree.fatJet1PtOverMHH)
    tmpXList.append(inputTree.fatJet2PtOverMHH)
    tmpXList.append(inputTree.ptj2_over_ptj1)

    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s.disc = tmpY

    ################ JES Up  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp
    tmpXList[2] = inputTree.hh_mass_JESUp
    #tmpXList[0] = inputTree.hh_pt #for new ntuples
    #tmpXList[1] = inputTree.hh_eta #for new ntuples
    #tmpXList[2] = inputTree.hh_mass #for new ntuples
    #tmpXList[7] = inputTree.fatJet1Pt_JES_Up
    #tmpXList[13] = inputTree.fatJet2Pt_JES_Up
    tmpXList[7] = inputTree.fatJet1Pt_JESUp #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp.disc = tmpY

    ################ JES Down  ###################
    tmpXList[0] = inputTree.hh_pt_JESDown
    tmpXList[1] = inputTree.hh_eta_JESDown
    tmpXList[2] = inputTree.hh_mass_JESDown
    #tmpXList[0] = inputTree.hh_pt #for new ntuples
    #tmpXList[1] = inputTree.hh_eta #for new ntuples
    #tmpXList[2] = inputTree.hh_mass #for new ntuples
    #tmpXList[7] = inputTree.fatJet1Pt_JES_Down
    #tmpXList[13] = inputTree.fatJet2Pt_JES_Down
    tmpXList[7] = inputTree.fatJet1Pt_JESDown #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown.disc = tmpY

    ################ JES Up Abs  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_Abs #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_Abs
    tmpXList[2] = inputTree.hh_mass_JESUp_Abs
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_Abs #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_Abs #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_Abs
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_Abs
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_Abs.disc = tmpY

    ################ JES Down Abs ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_Abs
    tmpXList[1] = inputTree.hh_eta_JESDown_Abs
    tmpXList[2] = inputTree.hh_mass_JESDown_Abs
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_Abs #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_Abs #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_Abs
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_Abs
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_Abs.disc = tmpY

    ################ JES Up Abs_2017  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_Abs_2017 #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_Abs_2017
    tmpXList[2] = inputTree.hh_mass_JESUp_Abs_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_Abs_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_Abs_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_Abs_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_Abs_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_Abs_2017.disc = tmpY

    ################ JES Down Abs_2017 ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_Abs_2017
    tmpXList[1] = inputTree.hh_eta_JESDown_Abs_2017
    tmpXList[2] = inputTree.hh_mass_JESDown_Abs_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_Abs_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_Abs_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_Abs_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_Abs_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_Abs_2017.disc = tmpY

    ################ JES Up BBEC1  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_BBEC1 #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_BBEC1
    tmpXList[2] = inputTree.hh_mass_JESUp_BBEC1
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_BBEC1 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_BBEC1 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_BBEC1
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_BBEC1
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_BBEC1.disc = tmpY

    ################ JES Down BBEC1 ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_BBEC1
    tmpXList[1] = inputTree.hh_eta_JESDown_BBEC1
    tmpXList[2] = inputTree.hh_mass_JESDown_BBEC1
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_BBEC1 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_BBEC1 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_BBEC1
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_BBEC1
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_BBEC1.disc = tmpY

    ################ JES Up BBEC1_2017  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_BBEC1_2017 #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_BBEC1_2017
    tmpXList[2] = inputTree.hh_mass_JESUp_BBEC1_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_BBEC1_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_BBEC1_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_BBEC1_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_BBEC1_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_BBEC1_2017.disc = tmpY

    ################ JES Down BBEC1_2017 ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_BBEC1_2017
    tmpXList[1] = inputTree.hh_eta_JESDown_BBEC1_2017
    tmpXList[2] = inputTree.hh_mass_JESDown_BBEC1_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_BBEC1_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_BBEC1_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_BBEC1_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_BBEC1_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_BBEC1_2017.disc = tmpY

    ################ JES Up EC2  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_EC2 #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_EC2
    tmpXList[2] = inputTree.hh_mass_JESUp_EC2
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_EC2 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_EC2 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_EC2
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_EC2
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_EC2.disc = tmpY

    ################ JES Down EC2 ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_EC2
    tmpXList[1] = inputTree.hh_eta_JESDown_EC2
    tmpXList[2] = inputTree.hh_mass_JESDown_EC2
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_EC2 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_EC2 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_EC2
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_EC2
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_EC2.disc = tmpY

    ################ JES Up EC2_2017  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_EC2_2017 #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_EC2_2017
    tmpXList[2] = inputTree.hh_mass_JESUp_EC2_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_EC2_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_EC2_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_EC2_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_EC2_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_EC2_2017.disc = tmpY

    ################ JES Down EC2_2017 ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_EC2_2017
    tmpXList[1] = inputTree.hh_eta_JESDown_EC2_2017
    tmpXList[2] = inputTree.hh_mass_JESDown_EC2_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_EC2_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_EC2_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_EC2_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_EC2_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_EC2_2017.disc = tmpY

    ################ JES Up HF  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_HF #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_HF
    tmpXList[2] = inputTree.hh_mass_JESUp_HF
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_HF #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_HF #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_HF
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_HF
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_HF.disc = tmpY

    ################ JES Down HF ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_HF
    tmpXList[1] = inputTree.hh_eta_JESDown_HF
    tmpXList[2] = inputTree.hh_mass_JESDown_HF
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_HF #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_HF #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_HF
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_HF
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_HF.disc = tmpY

    ################ JES Up HF_2017  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_HF_2017 #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_HF_2017
    tmpXList[2] = inputTree.hh_mass_JESUp_HF_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_HF_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_HF_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_HF_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_HF_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_HF_2017.disc = tmpY

    ################ JES Down HF_2017 ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_HF_2017
    tmpXList[1] = inputTree.hh_eta_JESDown_HF_2017
    tmpXList[2] = inputTree.hh_mass_JESDown_HF_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_HF_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_HF_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_HF_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_HF_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_HF_2017.disc = tmpY

    ################ JES Up FlavQCD  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_FlavQCD #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_FlavQCD
    tmpXList[2] = inputTree.hh_mass_JESUp_FlavQCD
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_FlavQCD #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_FlavQCD #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_FlavQCD
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_FlavQCD
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_FlavQCD.disc = tmpY

    ################ JES Down FlavQCD ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_FlavQCD
    tmpXList[1] = inputTree.hh_eta_JESDown_FlavQCD
    tmpXList[2] = inputTree.hh_mass_JESDown_FlavQCD
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_FlavQCD #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_FlavQCD #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_FlavQCD
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_FlavQCD
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_FlavQCD.disc = tmpY

    ################ JES Up RelBal  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_RelBal #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_RelBal
    tmpXList[2] = inputTree.hh_mass_JESUp_RelBal
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_RelBal #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_RelBal #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_RelBal
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_RelBal
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_RelBal.disc = tmpY

    ################ JES Down RelBal ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_RelBal
    tmpXList[1] = inputTree.hh_eta_JESDown_RelBal
    tmpXList[2] = inputTree.hh_mass_JESDown_RelBal
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_RelBal #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_RelBal #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_RelBal
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_RelBal
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_RelBal.disc = tmpY

    ################ JES Up RelSample_2017  ###################
    tmpXList[0] = inputTree.hh_pt_JESUp_RelSample_2017 #new ntuples didn't have this
    tmpXList[1] = inputTree.hh_eta_JESUp_RelSample_2017
    tmpXList[2] = inputTree.hh_mass_JESUp_RelSample_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESUp_RelSample_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESUp_RelSample_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESUp_RelSample_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESUp_RelSample_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESUp_RelSample_2017.disc = tmpY

    ################ JES Down RelSample_2017 ###################
    tmpXList[0] = inputTree.hh_pt_JESDown_RelSample_2017
    tmpXList[1] = inputTree.hh_eta_JESDown_RelSample_2017
    tmpXList[2] = inputTree.hh_mass_JESDown_RelSample_2017
    tmpXList[7] = inputTree.fatJet1Pt_JESDown_RelSample_2017 #for new ntuples
    tmpXList[13] = inputTree.fatJet2Pt_JESDown_RelSample_2017 #for new ntuples
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JESDown_RelSample_2017
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JESDown_RelSample_2017
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JESDown_RelSample_2017.disc = tmpY

    ################ JMS Up  ###################
    tmpXList[0] = inputTree.hh_pt_JMS_Up
    tmpXList[1] = inputTree.hh_eta_JMS_Up
    tmpXList[2] = inputTree.hh_mass_JMS_Up   
    #tmpXList[0] = inputTree.hh_pt #for new ntuples
    #tmpXList[1] = inputTree.hh_eta #for new ntuples
    #tmpXList[2] = inputTree.hh_mass #for new ntuples
    tmpXList[6] = inputTree.fatJet1MassSD_JMS_Up
    tmpXList[7] = inputTree.fatJet1Pt
    tmpXList[13] = inputTree.fatJet2Pt
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JMS_Up
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JMS_Up
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JMSUp.disc = tmpY

    ################ JMS Down  ###################
    tmpXList[0] = inputTree.hh_pt_JMS_Down
    tmpXList[1] = inputTree.hh_eta_JMS_Down
    tmpXList[2] = inputTree.hh_mass_JMS_Down
    #tmpXList[0] = inputTree.hh_pt #for new ntuples
    #tmpXList[1] = inputTree.hh_eta #for new ntuples
    #tmpXList[2] = inputTree.hh_mass #for new ntuples
    tmpXList[6] = inputTree.fatJet1MassSD_JMS_Down
    tmpXList[7] = inputTree.fatJet1Pt
    tmpXList[13] = inputTree.fatJet2Pt
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JMS_Down
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JMS_Down
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JMSDown.disc = tmpY

    ################ JMR Up  ###################
    tmpXList[0] = inputTree.hh_pt_JMR_Up
    tmpXList[1] = inputTree.hh_eta_JMR_Up
    tmpXList[2] = inputTree.hh_mass_JMR_Up
    #tmpXList[0] = inputTree.hh_pt #for new ntuples
    #tmpXList[1] = inputTree.hh_eta #for new ntuples
    #tmpXList[2] = inputTree.hh_mass #for new ntuples
    tmpXList[6] = inputTree.fatJet1MassSD_JMR_Up
    tmpXList[7] = inputTree.fatJet1Pt
    tmpXList[13] = inputTree.fatJet2Pt
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JMR_Up
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JMR_Up
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JMRUp.disc = tmpY

    ################ JMR Down  ###################
    tmpXList[0] = inputTree.hh_pt_JMR_Down
    tmpXList[1] = inputTree.hh_eta_JMR_Down
    tmpXList[2] = inputTree.hh_mass_JMR_Down
    #tmpXList[0] = inputTree.hh_pt #for new ntuples
    #tmpXList[1] = inputTree.hh_eta #for new ntuples
    #tmpXList[2] = inputTree.hh_mass #for new ntuples
    tmpXList[6] = inputTree.fatJet1MassSD_JMR_Down
    tmpXList[7] = inputTree.fatJet1Pt
    tmpXList[13] = inputTree.fatJet2Pt
    tmpXList[14] = inputTree.fatJet1PtOverMHH_JMR_Down
    tmpXList[15] = inputTree.fatJet2PtOverMHH_JMR_Down
    tmpXNPArray = np.array([ tmpXList ])
    tmpY = model.predict_proba(tmpXNPArray)[:, 1]
    my_s_JMRDown.disc = tmpY

    outputTree.Fill()

    

outputFile.Write()
outputFile.Close()
