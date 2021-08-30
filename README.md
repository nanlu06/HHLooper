# HHLooper
HHLooper

# Setup

```
cmsrel CMSSW_9_4_2
cd CMSSW_9_4_2/src
cmsenv
git clone git@github.com:LPC-HH/HHLooper.git
cd HHLooper
make
```

# About  the looper:

All the selections are defined  in app/HHLooper.cc in the cutflow object

`cutflow.addCut`: start a cut

`cutflow.addCutToLastActiveCut`:  add a cut on top of the cut in the previous line

`cutflow.getCut`: get a cut in a cutflow

To define a cut  in addCut or addCutToLastActiveCut, you need to give the cutname, a lambda function for a bool return, and a weight 

For example, cut `SRBin2` contains the following cutflow: `CutWeight` -> `CutHLT` -> `CutfatJetsPt` -> `CutfatJetsMassSD` -> `SRBin2`


All the output histograms are defined in  app/HHLooper.cc in the histograms object

Each output histogram will have the name of `CutName__HistogramName`, i.e. by default the program saves all histograms for all cuts in the cut flow defined above


# Run looper

```
./process.sh tag0 no nominal no no

```

The first argument is a tag for the folder name of the output histograms, and the second argument tells the program to include all the weight systematic Up/Down histograms, the third argument tells the program to run a specific shape systematic (reading a different branch in the ntuples. the systematics are JES_Up, JES_Down, JMR_Up, JMR_Down, JMS_Up, JMS_Down), then the next two arguments tells the programs to do trigsys, PNetSFsys or not. the trigsys is very time comsuming. 


## Command to run for the signal region

./run_syst_test.sh 

cd python

./prepare_histogram_forcard.sh

## Command to run for the ttbar one lepton CR

./process_1Lttbar.sh  yield_AN_1Lttbar

## Command to run for the ttbar two jet CR

./run_syst_ttbar_test.sh (including the recoil correction)

To remove the recoil correction, now comment Line 382, and uncomment Line  384 in ../app/HHLooper.cc to remove recoil correction;  

make clean; make

source process_ttbar.sh yield_AN_ttbar no nominal no no

# Make root file for datacard
```
cd python

```

The argument of the python script is the tag name you provided in the previous step

# make plots

```
cd python

```
