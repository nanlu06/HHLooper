change the path to the file that you want to use:

e.g.

TFile *f = new TFile("../../hists/yield_AN_ttbar_0329pileupweight/"+year+"/TTBarCR__hh_pt.root");


RampModel.C: produce figure 20, result of the ramp model fit for the recoil correction<br />
BDTshapeunc.C: produce BDT shape uncertainties for ttbar<br />
fit_topmass_inBDTbins.C: Figure 26 fitted top mass as a function of event BDT<br />
PNetcorrection.C : produce PNet ratio comparisons for ttbar after all correction