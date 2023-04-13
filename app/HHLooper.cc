//C++ INCLUDES
#include <sys/stat.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
//ROOT
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
//LOCAL INCLUDES
#include "hhtree.hh"
#include "anautil.h"
#include "bbbb_vs_bkg.h"
#include "scalefactors.h"

using namespace std;

int lumi = 137650.0;
TTJetsScaleFactors ttjets_sf;
TRandom3* r_nominal = new TRandom3(1);
const float pi = 3.141592653;

float DeltaPhi(float phi1, float phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

float DeltaR(float eta1, float phi1, float eta2, float phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}

float FatJetMassCorrection(string year, float mass, int type=0)
{
    //TRandom3* r_nominal = new TRandom3(1);
    //TRandomMT64* r_nominal = new TRandomMT64(1);
    //type=0: nominal
    //type=-1: jms only
    //type=-2: jmr only
    //type=1: jms down
    //type=2: jms up
    //type=3: jmr down
    //type=4: jmr up
    float* jmsValues;//{nominal, down, up}
    if(year == "2016")
      {
        float tmp_jms[] = {1.00, 0.9906, 1.0094};
        jmsValues = tmp_jms;
      }
    else if(year == "2017")
      {
        //float tmp_jms[] = {0.982, 0.978, 0.986};
        float tmp_jms[] = {1.0016, 0.978, 0.986};
        jmsValues = tmp_jms;
      }
    else if(year == "2018")
      {
        float tmp_jms[] = {0.997, 0.993, 1.001};
        jmsValues = tmp_jms;
      }
    else
      {
        std::cout << "year is not acceptable! Use: 2016, 2017, 2018" << std::endl;
        exit(EXIT_FAILURE);
      }

    float* jmrValues;//{nominal, down, up}
    if(year == "2016")
      {
        //float tmp_jmr[] = {1.00, 1.0, 1.2};
        float tmp_jmr[] = {1.00, 1.0, 1.09};
        jmrValues = tmp_jmr;
      }
    else if(year == "2017")
      {
        //float tmp_jmr[] = {1.09, 1.04, 1.14}; // percent of resolution
        //float tmp_jmr[] = {1.043, 1.00, 1.09};
        float tmp_jmr[] = {1.03, 1.00, 1.07};
        jmrValues = tmp_jmr;
      }
    else if(year == "2018")
      {
        //float tmp_jmr[] = {1.24, 1.20, 1.28}; //percent of resolution
        float tmp_jmr[] = {1.065, 1.031, 1.099};
        jmrValues = tmp_jmr;
      }
    else
      {
        std::cout << "year is not acceptable! Use: 2016, 2017, 2018" << std::endl;
        exit(EXIT_FAILURE);
      }

    float result = mass;
    float random = r_nominal->Gaus(0.0, 1); 
    float jmr_cor = (jmrValues[0]-1.0)*random; //r_nominal->Gaus( 0.0, jmrValues[0] -1.0 );
    float jmr_cor_dn = (jmrValues[1]-1.0)*random;
    float jmr_cor_up = (jmrValues[2]-1.0)*random;
    if(type==0) result = jmsValues[0]*mass*( 1.0 + jmr_cor ); //JMS + JMR correction
    if(type==-1) result = jmsValues[0]*mass; //JMS correction only
    if(type==-2) result = mass*( 1.0 + jmr_cor); //JMR correction only
    if(type==1) result = jmsValues[1]*mass*( 1.0 + jmr_cor ); // JMS down
    if(type==2) result = jmsValues[2]*mass*( 1.0 + jmr_cor ); // JMS up
    if(type==3) result = jmsValues[0]*mass*( 1.0 + jmr_cor_dn ); // JMR down 
    if(type==4) result = jmsValues[0]*mass*( 1.0 + jmr_cor_up ); // JMR up
    //delete r_nominal;
    return result;
}

int main ( int argc, char* argv[])
{

if(argc < 5) 
{
 cout<<"ERROR  - not enough arguments: need input(file, dir, or list) outputFileName label isData"<<endl;
 return 0;
}

std::string input = argv[1];
std::string outputFileName = argv[2];
std::string label = argv[3];
std::string isData_ = argv[4];
bool saveSkim = false;
bool doSystematics = false;
bool dotrigsys = false;
bool doPNetSFsys = false;
std::string syst_name = "";

if(argc > 5)
{
  std::string s_syst = argv[5];
  if(s_syst == "1" || s_syst == "syst" || s_syst == "yes") doSystematics = true; 
}

if(argc > 6)
{
  syst_name = argv[6]; //_JMS_Up, _JMS_Down, _JMR_Up, _JMR_Down, _JES_Up, _JES_Down, _JER_Up, _JER_Down for vars (fatJet2MassSD,fatJet1MassSD; fatJet2Pt,fatJet1Pt)
}

if(argc > 7)
{
  std::string s_trigsys = argv[7];
  if(s_trigsys == "1" || s_trigsys == "syst" || s_trigsys == "yes") dotrigsys = true; 
}

if(argc > 8)
{
  std::string s_PNetSFsys = argv[8];
  if(s_PNetSFsys == "1" || s_PNetSFsys == "syst" || s_PNetSFsys == "yes") doPNetSFsys = true;
}
    
if(argc > 9)
{
  std::string s_skim = argv[9];
  if(s_skim == "1" || s_skim == "skim" || s_skim == "yes") saveSkim = true; 
}

cout<< "debug: "<<input<<" outputFileName "<<outputFileName<<" label "<<label<<" isData_ "<<isData_<<" syst_name "<<syst_name<<" doSystematics "<<doSystematics<<" doTrigsys "<<dotrigsys <<" do PNet Hbb SF unc "<< doPNetSFsys <<endl;
    
system("mkdir -p hists");
system(("mkdir -p hists/"+label).c_str());

//luminosity
std::string year_ = "2022";
if(input.find("2022") != std::string::npos) {year_ = "2022"; lumi = 34000.0;}

//make folder to store the output histograms
system(("mkdir -p hists/"+label+"/"+year_).c_str());

//efficiency scale factors
mHH_THunc_ScaleFactors mhh_thunc_sf;
TrigEffScaleFactors trig_sf(year_);
miniIsoEleScaleFactors miniIsoEle_sf(year_);
miniIsoMuScaleFactors miniIsoMu_sf(year_);
MuTrigScaleFactors muTrig_sf(year_);
EleTrigScaleFactors elTrig_sf(year_);
MuIDScaleFactors muID_sf(year_);
EleIDScaleFactors elID_sf(year_);
vector<TString> trig_unc_names_up;
vector<TString> trig_unc_names_dn;
vector<int> trig_index;
trig_unc_names_up.clear(); 
trig_unc_names_dn.clear();
trig_index.clear();

//ParticleNet scale factors for H->bb 
PNetHbbScaleFactors PNet_sf(year_);

bool isData = false;
if(isData_ == "1" || isData_ == "true" || isData_ == "yes" || isData_ == "True" || isData_ == "Yes") isData = true;
if(isData) lumi = 1.0;
//if(input.find("qcd") != std::string::npos) lumi = lumi*0.72;

if(input.find("qcd") != std::string::npos and input.find("1LTopSkim") != std::string::npos) 
{
    if(year_ == "2016") lumi = lumi*0.8;
    if(year_ == "2017") lumi = lumi*2.1;
    if(year_ == "2018") lumi = lumi*1.9;
}

bool isTTJets = false;
if((input.find("ttbar") != std::string::npos) || (input.find("tt1L") != std::string::npos) ||  (input.find("tt2L") != std::string::npos)) isTTJets = true;

bool isHH = false;
if(outputFileName.find("HH") != std::string::npos) isHH = true;

bool isH = false;    
if(outputFileName.find("VH") != std::string::npos || outputFileName.find("ttH") != std::string::npos) isH = true;
 
bool isZJets = false;
if(outputFileName.find("vjets") != std::string::npos) isZJets = true;  

bool isVV = false;
if(outputFileName.find("vv") != std::string::npos) isVV = true;  
   
std::vector<std::string> list_chain;

if(input.find(".root") != std::string::npos) list_chain.push_back(input);//a file is given
else if (input.find(",") != std::string::npos)//a list is given, separated by ","
{
  int a =0;
  int b=input.find(",");
  while(1)
  {
	std::string thisFile = input.substr(a,b-a);
	list_chain.push_back(thisFile);
	a = b+1;
	if(input.find(",", a) != std::string::npos) b = input.find(",", a);
	else
	{
		b = input.size();
		std::string thisFile = input.substr(a,b-a);
		list_chain.push_back(thisFile);
	}
  }
  
}
else //a directory is given
{
  DIR* dirp = opendir(input.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != NULL)
  {
	std::string thisFile = dp->d_name;
	if(thisFile.find(".root") == std::string::npos) continue;
        cout <<"debug: input files: "<<thisFile<<endl;
	list_chain.push_back((input+"/"+thisFile));
  }

}

TChain * chain = new TChain("tree");
for(int idx = 0; idx < list_chain.size(); idx++) chain->Add(list_chain[idx].c_str());
int nEntries = chain->GetEntries();
cout<<"total number of events to process: "<<nEntries<<endl;

TFile *outfile = new TFile(("hists/"+label+"/"+year_+"/"+outputFileName).c_str(), "recreate");
TFile *outfile_skim;
if(saveSkim) outfile_skim = new TFile(("hists/"+label+"/"+year_+"/"+outputFileName.replace(outputFileName.find(".root"), 5, "_skim.root")).c_str(), "recreate");

RooUtil::Cutflow cutflow;
RooUtil::Histograms histograms;

//************define histograms**********//
histograms.addHistogram("yield",               "; yield; Events",                      1,    0.,   1.,    [&]() { return 0; } );

if(doSystematics)
{
    histograms.addHistogram("fatJet2MassSD",   "; j_{2} soft drop mass (GeV); Events", 46,   40.,    500.,  [&]() { return hh.fatJet2MassSD();} );
    if(input.find("Tau32TopSkim") != std::string::npos){
        histograms.addHistogram("fatJet1Pt",          "; p_{T}^{j1} (GeV); Events",           200,   300.,   2300.,  [&]() { return  hh.fatJet1Pt(); });
        histograms.addHistogram("fatJet2Pt",          "; p_{T}^{j2} (GeV); Events",           200,   300.,   2300.,  [&]() { return  hh.fatJet2Pt(); });
    }       
}
else
{
//final fit discriminant 
histograms.addHistogram("fatJet2MassSD",   "; j_{2} soft drop mass (GeV); Events", 46,   40.,    500.,  [&]() {
    return hh.fatJet2MassSD();
});
    
histograms.addHistogram("fatJet1MassSD",   "; j_{1} soft drop mass (GeV); Events", 46,   40.,    500.,  [&]() {
    return hh.fatJet1MassSD();
});

histograms.addHistogram("fatJet1Pt",          "; p_{T}^{j1} (GeV); Events",           200,   300.,   2300.,  [&]() { 
    return  hh.fatJet1Pt(); 
});

histograms.addHistogram("fatJet2Pt",          "; p_{T}^{j2} (GeV); Events",           200,   300.,   2300.,  [&]() { 
    return  hh.fatJet2Pt(); 
});

//other variables does not have correct JMS/JMR variations    
histograms.add2DHistogram("fat Jet2 MassSD vs pT", "mj2", 30,   50.,   200., "ptj2",  25,   250.,   750.,  [&]() { return  hh.fatJet2MassSD(); }, [&]() { return hh.fatJet2Pt();} );
histograms.add2DHistogram("fat Jet1 MassSD vs pT", "mj1", 30,   50.,   200., "ptj1",  25,   250.,   750.,  [&]() { return  hh.fatJet1MassSD(); }, [&]() { return hh.fatJet1Pt();} );
histograms.add2DHistogram("fat Jet1 pT vs Jet2 pT", "ptj1", 72,   200.,   2000., "ptj2",  72,   200.,   2000.,  [&]() { return  hh.fatJet1Pt(); }, [&]() { return hh.fatJet2Pt();} );

histograms.addHistogram("fatJet1MassSDLR",   "; j_{1} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet1MassSD(); } );
//histograms.addHistogram("fatJet2MassSD",   "; j_{2} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet2MassSD(); } );
histograms.addHistogram("fatJet1MassSD_raw",   "; j_{1} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet1MassSD_UnCorrected(); } );
histograms.addHistogram("fatJet2MassSD_raw",   "; j_{2} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet2MassSD_UnCorrected(); } );
histograms.addHistogram("MET",           "; p_{T}^{miss} (GeV); Events",         200,   0.,   500.,  [&]() { return hh.MET(); } );
histograms.addHistogram("hh_pt",               "; p_{T}^{jj} (GeV); Events",           {0.,50., 100., 150., 200., 250., 300., 400., 500., 600., 800., 1000.},  [&]() { return hh.hh_pt(); } );
//histograms.addHistogram("hh_pt",               "; p_{T}^{jj} (GeV); Events",           {0.,50., 100., 150., 300., 1000.},  [&]() { return hh.hh_pt(); } );
histograms.addHistogram("hh_eta",               "; #eta^{jj}; Events",                 200,   -5.0,  5.0,  [&]() { return hh.hh_eta(); } );
histograms.addHistogram("hh_phi",               "; #Phi^{jj}; Events",                 200,   -3.2,  3.2,  [&]() { return hh.hh_phi(); } );
histograms.addHistogram("hh_mass",             "; m_{jj} (GeV); Events",               200,   0.,  1500.,  [&]() { return hh.hh_mass(); } );
histograms.addHistogram("hh_mass_lr",             "; m_{jj} (GeV); Events",               200,   0.,  3000.,  [&]() { return hh.hh_mass(); } );
histograms.addHistogram("fatJet1PNetXbb",   "; j_{1} PNet Xbb tagger; Events",           200,   0.0,  1.0,   [&]() { return  hh.fatJet1PNetXbb(); } );
histograms.addHistogram("fatJet2PNetXbb",   "; j_{2} PNet Xbb tagger; Events",           200,   0.0,  1.0,   [&]() { return  hh.fatJet2PNetXbb(); } );
histograms.addHistogram("fatJet1PNetXbb_Bin1",   "; j_{1} PNet Xbb tagger; Events",    {0.90, 0.95,  0.975, 0.985,  1.00} ,   [&]() { return  hh.fatJet1PNetXbb(); } );
histograms.addHistogram("fatJet2PNetXbb_Bin1",   "; j_{2} PNet Xbb tagger; Events",    {0.90, 0.95,  0.975, 0.985,  1.00} ,   [&]() { return  hh.fatJet2PNetXbb(); } );
histograms.addHistogram("fatJet1PNetXbb_Bin2",   "; j_{1} PNet Xbb tagger; Events",    {0.90, 0.945, 0.955,  0.975, 0.985,  1.00} ,   [&]() { return  hh.fatJet1PNetXbb(); } );
histograms.addHistogram("fatJet2PNetXbb_Bin2",   "; j_{2} PNet Xbb tagger; Events",    {0.90, 0.945, 0.955,  0.975, 0.985,  1.00} ,   [&]() { return  hh.fatJet2PNetXbb(); } );
histograms.addHistogram("fatJet1Eta",          "; #eta^{j1}; Events",                 200,   -2.5,  2.5,  [&]() { return  hh.fatJet1Eta(); } );
histograms.addHistogram("fatJet1Phi",          "; #Phi^{j1}; Events",                 200,  -3.2,   3.2,  [&]() { return  hh.fatJet1Phi(); } );
histograms.addHistogram("fatJet2Eta",          "; #eta^{j2}; Events",                 200,   -2.5,  2.5,  [&]() { return  hh.fatJet2Eta(); } );
histograms.addHistogram("fatJet2Phi",          "; #Phi^{j2}; Events",                 200,  -3.2,   3.2,  [&]() { return  hh.fatJet2Phi(); } );
histograms.addHistogram("abs_dEta_j1j2",       "; #Delta#eta(j_{1}, j_{2}); Events",   200,   0.,   5.,    [&]() { return  fabs(hh.fatJet1Eta() - hh.fatJet2Eta()); } );
histograms.addHistogram("abs_dPhi_j1j2",       "; #Delta#Phi(j_{1}, j_{2}); Events",   200,   2.,   4.5,    [&]() { return  fabs(hh.fatJet1Phi() - hh.fatJet2Phi()); } );
histograms.addHistogram("abs_dR_j1j2",       "; #DeltaR(j_{1}, j_{2}); Events",        200,   0.,   5.0,    [&]() { return  sqrt((hh.fatJet1Eta() - hh.fatJet2Eta())*(hh.fatJet1Eta() - hh.fatJet2Eta())  + pow(fabs(hh.fatJet1Phi()-hh.fatJet2Phi()) > pi ? fabs(hh.fatJet1Phi()-hh.fatJet2Phi()) - 2*pi : hh.fatJet1Phi()-hh.fatJet2Phi(), 2)); } );
histograms.addHistogram("ptj1_over_mhh",       "; p_{T}^{j1}/m_{HH}; Events",         200,   0.,   1.,    [&]() { return  hh.fatJet1PtOverMHH(); } );
histograms.addHistogram("ptj2_over_mhh",       "; p_{T}^{j2}/m_{HH}; Events",         200,   0.,   1.,    [&]() { return  hh.fatJet2PtOverMHH(); } );
histograms.addHistogram("ptj1_over_mj1",       "; p_{T}^{j1}/m_{j1}; Events",         200,   0.,   10.,   [&]() { return  hh.fatJet1PtOverMSD(); } );
histograms.addHistogram("ptj2_over_mj2",       "; p_{T}^{j2}/m_{j2}; Events",         200,   0.,   10.,   [&]() { return  hh.fatJet2PtOverMSD(); } );
//histograms.addHistogram("ptj1_over_mregj1",       "; p_{T}^{j1}/Mreg_{j1}; Events",         200,   0.,   10.,   [&]() { return  hh.fatJet1PtOverMRegressed(); } );
//histograms.addHistogram("ptj2_over_mregj2",       "; p_{T}^{j2}/Mreg_{j2}; Events",         200,   0.,   10.,   [&]() { return  hh.fatJet2PtOverMRegressed(); } );
histograms.addHistogram("ptj2_over_ptj1",      "; p_{T}^{j2}/p_{T}^{j1}; Events",     200,   0.5,  1.,    [&]() { return  hh.fatJet2Pt() / hh.fatJet1Pt(); } );
histograms.addHistogram("mj2_over_mj1",      "; m^{j2}/m^{j1}; Events",               200,   0.0,  1.5,   [&]() { return  hh.fatJet2MassSD() / hh.fatJet1MassSD(); } );
histograms.addHistogram("fatJet1Tau3OverTau2",   "; j_{1} Tau3/2; Events",           200,   0.0,  1.0,   [&]() { return  hh.fatJet1Tau3OverTau2(); } );
histograms.addHistogram("fatJet2Tau3OverTau2",   "; j_{2} Tau3/2; Events",           200,   0.0,  1.0,   [&]() { return  hh.fatJet2Tau3OverTau2(); } );
histograms.addHistogram("lep1Pt",          "; #p_{T}^{l1}; Events",                 200,   0.,  500.,  [&]() { return  hh.lep1Pt(); } );
histograms.addHistogram("lep1Eta",          "; #eta^{j1}; Events",                 200,   -2.5,  2.5,  [&]() { return  hh.lep1Eta(); } );
histograms.addHistogram("lep1Phi",          "; #Phi^{j1}; Events",                 200,  -3.2,   3.2,  [&]() { return  hh.lep1Phi(); } );
 histograms.addHistogram("gen_mHH",          "; gen m_{HH}; Events",                 200, 0.,   1500.,[&]() {  TLorentzVector gh1,gh2;
	  gh1.SetPtEtaPhiM(hh.genHiggs1Pt(), hh.genHiggs1Eta(),hh.genHiggs1Phi(),125.0);
	  gh2.SetPtEtaPhiM(hh.genHiggs2Pt(), hh.genHiggs2Eta(),hh.genHiggs2Phi(),125.0);
	  return (gh1+gh2).M();} );
if(input.find("1LTopSkim") != std::string::npos) histograms.addHistogram("abs_dR_l1j1",       "; #DeltaR(l_{1}, j_{1}); Events",        200,   0.,   5.0,    [&]() { return  sqrt(pow(hh.lep1Eta() - hh.fatJet1Eta(), 2) + pow(fabs(hh.lep1Phi()-hh.fatJet1Phi()) > pi ? fabs(hh.lep1Phi()-hh.fatJet1Phi()) - 2*pi : hh.lep1Phi()-hh.fatJet1Phi(), 2)); } );
 
}
    
//************define cuts**********//

cutflow.setTFile(outfile);

if(input.find("1LTopSkim") != std::string::npos) // this is 1LTopSkim input
{ 
//Pre-selection cuts
cutflow.addCut("CutWeight", [&](){ return 1; },   [&](){ return isData ?  lumi : lumi * hh.weight() *(abs(hh.lep1Id()) == 11 ? miniIsoEle_sf.getminiIsoScaleFactors(hh.lep1Pt(),hh.lep1Eta()) : miniIsoMu_sf.getminiIsoScaleFactors(hh.lep1Pt(),hh.lep1Eta())) * (abs(hh.lep1Id()) == 11 ? elTrig_sf.getScaleFactors(hh.lep1Pt(),hh.lep1Eta()) : muTrig_sf.getTrigScaleFactors(hh.lep1Pt(),hh.lep1Eta(), year_)) * (abs(hh.lep1Id()) == 11 ? elID_sf.getScaleFactors(hh.lep1Pt(),hh.lep1Eta()) : muID_sf.getIDScaleFactors(hh.lep1Pt(),hh.lep1Eta(), year_))*(year_ == "2016" ? 0.813367:1.0)*(year_ == "2017" ?0.955044 :1.0)*(year_ == "2018" ?0.950837 :1.0) ;}); //before correction //  * (isTTJets ?ttjets_sf.getScaleFactors_VBF(year_, hh.fatJet1Pt()) :1.0)
//cutflow.addCut("CutWeight", [&](){ return 1; },   [&](){ return isData ?  lumi : lumi * hh.weight() * hh.pileupWeight() * (isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0) : 1.0); });//after correction
//if(input.find("HHc1") == std::string::npos) cutflow.addCutToLastActiveCut("CutHLT",       [&](){ return abs(hh.lep1Id()) == 11 ? (hh.HLT_Ele27_WPTight_Gsf() || hh.HLT_Ele32_WPTight_Gsf() || hh.HLT_Ele32_WPTight_Gsf_L1DoubleEG()) : (hh.HLT_IsoMu24()  ||  hh.HLT_IsoMu27()  || hh.HLT_Mu50() ); },   UNITY);
cutflow.addCutToLastActiveCut("CutLepJetPt",       [&](){ return hh.fatJet1Pt() > 300.0 && hh.lep1Pt() > 50.0 && hh.lep2Pt() <=0 ; },   UNITY);
cutflow.addCutToLastActiveCut("CutfatJetMassSD",       [&](){ return hh.fatJet1MassSD() > 50.0; },   UNITY);
//ttbar 1L+jet control region
cutflow.addCutToLastActiveCut("CutLepEta",       [&](){ return (abs(hh.lep1Id()) == 11 && fabs(hh.lep1Eta()) <  2.5) || (abs(hh.lep1Id()) == 13 && fabs(hh.lep1Eta()) <  2.4); },   UNITY);
cutflow.addCutToLastActiveCut("CutMET",       [&](){ return hh.MET() > 100.0; },   UNITY);
//cutflow.addCutToLastActiveCut("CutHEM2018",       [&](){ return ((year_ == "2018" && isData && hh.run() >=319077) ? !((abs(hh.lep1Id()) == 11 && hh.lep1Phi() > -1.57 && hh.lep1Phi() < -0.87 && hh.lep1Eta() < -1.3) || (hh.fatJet1Phi() > -1.57 && hh.fatJet1Phi() < -0.87 && hh.fatJet1Eta()< -1.3 ) ): true); },   UNITY);
cutflow.addCutToLastActiveCut("CutTau3Over2",       [&](){ return hh.fatJet1MassSD() > 140.0 && hh.fatJet1Tau3OverTau2() < 0.46 ; },   [&]() {return isTTJets ? TopTagSF("0.46", year_, hh.fatJet1Pt() ) : 1.0;} );
cutflow.addCutToLastActiveCut("TTBarLepJetCR",       [&](){ return sqrt(pow(hh.lep1Eta() - hh.fatJet1Eta(), 2) + pow(fabs(hh.lep1Phi()-hh.fatJet1Phi()) > pi ? fabs(hh.lep1Phi()-hh.fatJet1Phi()) - 2*pi : hh.lep1Phi()-hh.fatJet1Phi(), 2)) > 1.0; },   UNITY); //delta R > 1.0
cutflow.addCutToLastActiveCut("TTBarLepJetCRElectron",       [&](){ return abs(hh.lep1Id()) == 11; },   UNITY);
cutflow.getCut("TTBarLepJetCR");
cutflow.addCutToLastActiveCut("TTBarLepJetCRMuon",       [&](){ return abs(hh.lep1Id()) == 13; },   UNITY);
}
//signal region and ttbar two jet control region
else
{
//Pre-selection cuts
  cutflow.addCut("CutWeight", [&](){ return 1; },  [&](){
      //after ttbar recoil correction
      float ttbar_factor = input.find("option5") == std::string::npos ? 1.0 : 2.0;
     
      float total_weight = isData ?  lumi :lumi * hh.weight() * (isTTJets ? ttbar_factor*ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 0) : 1.0)*((isHH && (outputFileName.find("VBF")== std::string::npos))? mhh_thunc_sf.getmHHTHuncScaleFactors(hh.hh_mass(),0):1.0);
      //If applying correction from VBF analysis
      //float total_weight = isData ?  lumi :lumi * hh.l1PreFiringWeight() * hh.pileupWeight() * hh.xsecWeight() * (isHH? hh.weight() : hh.genWeight());
      //if(isTTJets){
      //	if (hh.fatJet1Pt()>hh.fatJet2Pt())total_weight *= ttjets_sf.getScaleFactors_VBF(year_, hh.fatJet1Pt()) ;
      //	else total_weight *= ttjets_sf.getScaleFactors_VBF(year_, hh.fatJet2Pt()) ;
      //}

    //before ttbar recoil correction
      // float total_weight = isData ?  lumi :lumi * hh.l1PreFiringWeight() * hh.pileupWeight() * hh.xsecWeight() * (isHH? hh.weight() : hh.genWeight()* (isTTJets ? ttbar_factor:1.0));
    
    //apply PNet SF for Hbb
    if(isH){ //genHiggs1Eta
        float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
        float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
        if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) )  total_weight = total_weight * PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
        else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) )  total_weight = total_weight * PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
    }
    else if(isHH){
        total_weight = total_weight * PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
    } 
       
    return total_weight;
}); 
 
cutflow.addCutToLastActiveCut("CutHLT",       [&](){ 
   return isData ? ((year_ == "2022" && (hh.HLT_PFHT780() || hh.HLT_PFHT890() || hh.HLT_PFHT1050() || hh.HLT_AK8PFJet360_TrimMass30() || hh.HLT_AK8PFJet380_TrimMass30()
   ||hh.HLT_AK8PFJet400_TrimMass30() || hh.HLT_AK8PFJet420_TrimMass30()|| hh.HLT_AK8PFHT750_TrimMass50() || hh.HLT_AK8PFHT800_TrimMass50()
   ||hh.HLT_AK8PFHT850_TrimMass50() || hh.HLT_AK8PFHT900_TrimMass50() || hh.HLT_PFJet450() || hh.HLT_PFJet500() || hh.HLT_PFJet550()
   ||hh.HLT_AK8PFJet400() || hh.HLT_AK8PFJet450() || hh.HLT_AK8PFJet500() || hh.HLT_AK8PFJet550() || hh.HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17()
   ||hh.HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1() || hh.HLT_AK8PFJet330_PFAK8BTagCSV_p17() || hh.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02()
   ||hh.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2() ||hh.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4() ||hh.HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20() 
   ||hh.HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087() ||hh.HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087() ||hh.HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20() 
   ||hh.HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20() ||hh.HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20())) ) : 1.0; },   UNITY);
    
cutflow.addCutToLastActiveCut("CutfatJetsPt",       [&](){
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.isVBFtag() < 1 && hh.fatJet1Pt() > 300.0 && hh.fatJet2Pt() > 300.0;
    else return hh.isVBFtag() < 1 && hh.fatJet1Pt() > 300.0 && hh.fatJet2Pt() > 300.0;
},   UNITY);
cutflow.addCutToLastActiveCut("CutfatJetsMassSD",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.fatJet1MassSD() > 50.0 && hh.fatJet2MassSD() > 50.0;  
    else return hh.fatJet1MassSD() > 50.0 && hh.fatJet2MassSD() > 50.0;
},   UNITY);

//Signal regions - pass - BDT v8p2
if(input.find("Tau32TopSkim") == std::string::npos){

cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1",       [&](){    
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.fatJet2PNetXbb() > 0.980; 
    else return hh.fatJet2PNetXbb() > 0.980; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
//cutflow.addCutToLastActiveCut("J2MassSideBandv8p2Bin1",       [&](){ return hh.fatJet2MassSD() <= 110 || hh.fatJet2MassSD() >= 140.0; },   UNITY);

cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin2",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return (!(hh.fatJet2PNetXbb() > 0.980))  && ((hh.fatJet2PNetXbb() > 0.980)||(hh.fatJet2PNetXbb() > 0.950)); 
    else return (!(hh.fatJet2PNetXbb() > 0.980))  && ((hh.fatJet2PNetXbb() > 0.980)||(hh.fatJet2PNetXbb() > 0.950)); 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
//cutflow.addCutToLastActiveCut("J2MassSideBandv8p2Bin2",       [&](){ return hh.fatJet2MassSD() <= 110 || hh.fatJet2MassSD() >= 140.0; },   UNITY);

cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin3",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return (!((hh.fatJet2PNetXbb() > 0.980)||(hh.fatJet2PNetXbb() > 0.950)))  && hh.fatJet2PNetXbb() > 0.950; 
    else return (!(( hh.fatJet2PNetXbb() > 0.980)||( hh.fatJet2PNetXbb() > 0.950)))  && hh.fatJet2PNetXbb() > 0.950; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
//cutflow.addCutToLastActiveCut("J2MassSideBandv8p2Bin3",       [&](){ return hh.fatJet2MassSD() <= 110 || hh.fatJet2MassSD() >= 140.0; },   UNITY);

cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("FailSRv8p2",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return  hh.fatJet2PNetXbb() < 0.950; 
    else return hh.fatJet2PNetXbb() < 0.950; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0) : 1.0; }); 
//cutflow.addCutToLastActiveCut("J2MassSideBandFailSRv8p2",       [&](){ return hh.fatJet2MassSD() <= 110 || hh.fatJet2MassSD() >= 140.0; },   UNITY);

////BDT inverted regions for F-test
//cutflow.getCut("CutfatJetsMassSD");
//cutflow.addCutToLastActiveCut("FitCRv24",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v24() < 0.0024 && hh.fatJet1PNetXbb() > 0.945 && hh.fatJet2PNetXbb() > 0.945; },   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0) * ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)     : 1.0; });
//cutflow.getCut("CutfatJetsMassSD");
//cutflow.addCutToLastActiveCut("FailFitCRv24",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v24() < 0.0024 && hh.fatJet1PNetXbb() > 0.945 && hh.fatJet2PNetXbb() < 0.945; },   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0) * ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)     : 1.0; });

//cutflow.getCut("CutfatJetsMassSD");
//cutflow.addCutToLastActiveCut("FitCRv8p2",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() < 0.03 && hh.fatJet2PNetXbb() > 0.95; },   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)     : 1.0; });
//cutflow.getCut("CutfatJetsMassSD");
//cutflow.addCutToLastActiveCut("FailFitCRv8p2",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() < 0.03 && hh.fatJet2PNetXbb() < 0.95; },   [&](){ return isTTJets  ? //ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)     : 1.0; });

        
}

////Top control region
else{
    cutflow.getCut("CutfatJetsMassSD");
    cutflow.addCutToLastActiveCut("CutfatJetsXbb",       [&](){ return hh.fatJet1PNetXbb() > 0.1 && hh.fatJet2PNetXbb() > 0.1; },   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0) * ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)     : 1.0; });
    cutflow.addCutToLastActiveCut("TTBarCR",       [&](){ 
    return hh.fatJet1Tau3OverTau2() < 0.46 && hh.fatJet2Tau3OverTau2() < 0.46 &&  hh.fatJet1HasBJetCSVLoose() && hh.fatJet2HasBJetCSVLoose(); },   [&]() {
        return isTTJets ? (TopTagSF("0.46", year_, hh.fatJet1Pt()) * TopTagSF("0.46", year_, hh.fatJet2Pt())) : 1.0;
    } );
    cutflow.getCut("CutfatJetsXbb");
    cutflow.addCutToLastActiveCut("TTBarCRTight",       [&](){ return hh.fatJet1Tau3OverTau2() < 0.39 && hh.fatJet2Tau3OverTau2() < 0.39 &&  hh.fatJet1HasBJetCSVLoose() && hh.fatJet2HasBJetCSVLoose(); },   [&]() {return isTTJets ? (TopTagSF("0.40", year_, hh.fatJet1Pt()) * TopTagSF("0.40", year_, hh.fatJet2Pt())) : 1.0;} );  
    }
}


/****Systematics******/
if(doSystematics && (outputFileName.find("qcd") == std::string::npos ) && (outputFileName.find("data") == std::string::npos ) && (syst_name.find("nominal") != std::string::npos) )
{   
    //ttbar recoil correction uncertainty for ttbar
    cutflow.addWgtSyst("ttJetsCorrUp",  [&](){return isTTJets ?  ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 1)/ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 0) : 1.0;});
    cutflow.addWgtSyst("ttJetsCorrDown",  [&](){return isTTJets ? ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), -1)/ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 0) : 1.0;});
    
    //PNet shape uncertainty for ttbar 
    cutflow.addWgtSyst("PNetShapeUp",  [&](){return isTTJets ?  ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 1)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 1)/(ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)) : 1.0;});
    cutflow.addWgtSyst("PNetShapeDown",  [&](){return isTTJets ?  ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), -1)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), -1)/(ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)) : 1.0;});
    
    //pileup reweighting uncertainty
    //cutflow.addWgtSyst("pileupWeightUp",  [&](){return hh.puWeight()!=0 ? hh.puWeightUp()/hh.puWeight() : 1.0;});
    //cutflow.addWgtSyst("pileupWeightDown",  [&](){return hh.puWeight()!=0 ? hh.puWeightDown()/hh.puWeight() : 1.0;});
}


//book histograms for cuts
if(not doSystematics) cutflow.bookHistogramsForCutAndBelow(histograms, "CutWeight");
else
{
    if(input.find("Tau32TopSkim") != std::string::npos) // this is ttbar enriched input
    {
        cutflow.bookHistogramsForCut(histograms, "TTBarCR");
    }
    else // normal  signal enriched input
    { 
        cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1");
        cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin2");
        cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin3");
        cutflow.bookHistogramsForCut(histograms, "FailSRv8p2");

    }
}
cutflow.bookCutflows();

int iEntry = 0;

if(saveSkim) outfile_skim->cd();
TTree *tree_out;
    
int BDTcat_index;
UInt_t run;
UInt_t luminosityBlock;
ULong64_t event;
float weight;
float fatJet1MassSD;
float fatJet2MassSD;
float fatJet1MassRegressed;
float fatJet2MassRegressed;
float fatJet1PNetXbb;
float fatJet2PNetXbb;
float fatJet1Pt;
float fatJet2Pt;
float fatJet1Eta;
float fatJet2Eta;
float fatJet1Phi;
float fatJet2Phi;
float fatJet1PtOverMHH;
float fatJet2PtOverMHH;
float fatJet1PtOverMSD;
float fatJet2PtOverMSD;
float fatJet1PtOverMRegressed;
float fatJet2PtOverMRegressed;
float abs_dEta_j1j2;
float abs_dPhi_j1j2;
float abs_dR_j1j2;
float ptj2_over_ptj1;
float mj2_over_mj1;
float mregj2_over_mregj1;
float hh_pt;
float hh_eta;
float hh_phi;
float hh_mass;
float gen_hh_mass;

if(saveSkim)
{ 
    tree_out = new TTree("hh", "output skim tree");
    tree_out->Branch("BDTcat_index", &BDTcat_index, "BDTcat_index/I");
    tree_out->Branch("run", &run, "run/I");
    tree_out->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/I");
    tree_out->Branch("event", &event, "event/l");
    tree_out->Branch("weight", &weight, "weight/F");
    tree_out->Branch("fatJet1MassSD", &fatJet1MassSD, "fatJet1MassSD/F");
    //tree_out->Branch("fatJet1MassRegressed", &fatJet1MassRegressed, "fatJet1MassRegressed/F");
    tree_out->Branch("fatJet1PNetXbb", &fatJet1PNetXbb, "fatJet1PNetXbb/F");
    tree_out->Branch("fatJet1Pt", &fatJet1Pt, "fatJet1Pt/F");
    tree_out->Branch("fatJet1Eta", &fatJet1Eta, "fatJet1Eta/F");
    tree_out->Branch("fatJet1Phi", &fatJet1Phi, "fatJet1Phi/F");
    tree_out->Branch("fatJet1PtOverMHH", &fatJet1PtOverMHH, "fatJet1PtOverMHH/F");
    tree_out->Branch("fatJet1PtOverMSD", &fatJet1PtOverMSD, "fatJet1PtOverMSD/F");
    //tree_out->Branch("fatJet1PtOverMRegressed", &fatJet1PtOverMRegressed, "fatJet1PtOverMRegressed/F");
    tree_out->Branch("fatJet2MassSD", &fatJet2MassSD, "fatJet2MassSD/F");
    //tree_out->Branch("fatJet2MassRegressed", &fatJet2MassRegressed, "fatJet2MassRegressed/F");
    tree_out->Branch("fatJet2PNetXbb", &fatJet2PNetXbb, "fatJet2PNetXbb/F");
    tree_out->Branch("fatJet2Pt", &fatJet2Pt, "fatJet2Pt/F");
    tree_out->Branch("fatJet2Eta", &fatJet2Eta, "fatJet2Eta/F");
    tree_out->Branch("fatJet2Phi", &fatJet2Phi, "fatJet2Phi/F");
    tree_out->Branch("fatJet2PtOverMHH", &fatJet2PtOverMHH, "fatJet2PtOverMHH/F");
    tree_out->Branch("fatJet2PtOverMSD", &fatJet2PtOverMSD, "fatJet2PtOverMSD/F");
    //tree_out->Branch("fatJet2PtOverMRegressed", &fatJet2PtOverMRegressed, "fatJet2PtOverMRegressed/F");
    tree_out->Branch("abs_dEta_j1j2", &abs_dEta_j1j2, "abs_dEta_j1j2/F");
    tree_out->Branch("abs_dPhi_j1j2", &abs_dPhi_j1j2, "abs_dPhi_j1j2/F");
    tree_out->Branch("abs_dR_j1j2", &abs_dR_j1j2, "abs_dR_j1j2/F");
    tree_out->Branch("ptj2_over_ptj1", &ptj2_over_ptj1, "ptj2_over_ptj1/F");
    tree_out->Branch("mj2_over_mj1", &mj2_over_mj1, "mj2_over_mj1/F");
    tree_out->Branch("mregj2_over_mregj1", &mregj2_over_mregj1, "mj2_over_mregj1/F");
    tree_out->Branch("hh_pt", &hh_pt, "hh_pt/F");
    tree_out->Branch("hh_eta", &hh_eta, "hh_eta/F");
    tree_out->Branch("hh_phi", &hh_phi, "hh_phi/F");
    tree_out->Branch("hh_mass", &hh_mass, "hh_mass/F");
    tree_out->Branch("gen_hh_mass", &gen_hh_mass, "gen_hh_mass/F");
}

for(int idx = 0; idx < list_chain.size(); idx++)
{
  cout<<"[INFO] processing file: "<<list_chain[idx]<<endl;
  TTree * tree_this;
  TFile * file_this = new TFile(list_chain[idx].c_str(), "READ");
  tree_this = (TTree*)file_this->Get("tree");
  hh.Init(tree_this);
  int nEntries_this = tree_this->GetEntries();
  for(int iEntry_this=0; iEntry_this<nEntries_this; iEntry_this++)
  {
	hh.GetEntry(iEntry_this);
        if(saveSkim) outfile->cd();
	cutflow.fill();
	if(saveSkim && (cutflow.getCut("SRv8p2Bin1").pass || cutflow.getCut("SRv8p2Bin2").pass || cutflow.getCut("SRv8p2Bin3").pass))
	{
	  outfile_skim->cd();	
	  BDTcat_index = -1;
	  if(cutflow.getCut("SRv8p2Bin1").pass) BDTcat_index = 1;
	  else if(cutflow.getCut("SRv8p2Bin2").pass) BDTcat_index = 2;
	  else if(cutflow.getCut("SRv8p2Bin3").pass) BDTcat_index = 3;
	  weight = isData ?  1.0 : lumi*hh.weight();
	  fatJet1MassSD = hh.fatJet1MassSD();
	  //fatJet1MassRegressed = hh.fatJet1MassRegressed();
	  fatJet1PNetXbb = hh.fatJet1PNetXbb();
	  fatJet1Pt = hh.fatJet1Pt();
	  fatJet1Eta = hh.fatJet1Eta();
	  fatJet1Phi = hh.fatJet1Phi();
	  fatJet1PtOverMHH = hh.fatJet1PtOverMHH();
	  fatJet1PtOverMSD = hh.fatJet1PtOverMSD();
	  //fatJet1PtOverMRegressed = hh.fatJet1PtOverMRegressed();
	  fatJet2MassSD = hh.fatJet2MassSD();
	  //fatJet2MassRegressed = hh.fatJet2MassRegressed();
	  fatJet2PNetXbb = hh.fatJet2PNetXbb();
	  fatJet2Pt = hh.fatJet2Pt();
	  fatJet2Eta = hh.fatJet2Eta();
	  fatJet2Phi = hh.fatJet2Phi();
	  fatJet2PtOverMHH = hh.fatJet2PtOverMHH();
	  fatJet2PtOverMSD = hh.fatJet2PtOverMSD();
	  //fatJet2PtOverMRegressed = hh.fatJet2PtOverMRegressed();
	  abs_dEta_j1j2 = fabs(hh.fatJet1Eta() - hh.fatJet2Eta());
	  abs_dPhi_j1j2 = fabs(hh.fatJet1Phi() - hh.fatJet2Phi());
	  abs_dR_j1j2 = sqrt((hh.fatJet1Eta() - hh.fatJet2Eta())*(hh.fatJet1Eta() - hh.fatJet2Eta())  + (hh.fatJet1Phi() - hh.fatJet2Phi())*(hh.fatJet1Phi() - hh.fatJet2Phi()));
	  ptj2_over_ptj1 = hh.fatJet2Pt() / hh.fatJet1Pt();
	  mj2_over_mj1 = fatJet2MassSD/fatJet1MassSD;
	  //mregj2_over_mregj1 = fatJet2MassRegressed/fatJet1MassRegressed;
	  hh_pt = hh.hh_pt();
	  hh_eta = hh.hh_eta();
	  hh_phi = hh.hh_phi();
	  hh_mass = hh.hh_mass();
	  TLorentzVector gh1,gh2;
	  gh1.SetPtEtaPhiM(hh.genHiggs1Pt(), hh.genHiggs1Eta(),hh.genHiggs1Phi(),125.0);
	  gh2.SetPtEtaPhiM(hh.genHiggs2Pt(), hh.genHiggs2Eta(),hh.genHiggs2Phi(),125.0);
	  gen_hh_mass = (gh1+gh2).M();
	  run = hh.run();
	  luminosityBlock = hh.lumi();
	  event = hh.event();
	  tree_out->Fill();
	}
	if(iEntry%100000 == 0) cout<<"[INFO] processing event "<<iEntry<<" / "<<nEntries<<endl;
	iEntry ++;
  }
  delete tree_this;
  file_this->Close();
  delete file_this;
}

//save histograms
cutflow.saveOutput();
outfile->Close();
if(saveSkim) 
{
 outfile_skim->cd();
 tree_out->Write();
 outfile_skim->Close();
}
cout<<"[INFO]: all  files successfully processed... "<<endl;

return 0;
}
