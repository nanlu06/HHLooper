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

int lumi = 137000.0;
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

std::string year_ = "2016";
if(input.find("2016") != std::string::npos) {year_ = "2016"; lumi = 35922.0;}
if(input.find("2017") != std::string::npos) {year_ = "2017"; lumi = 41480.0;}
if(input.find("2018") != std::string::npos) {year_ = "2018"; lumi = 59741.0;}
system(("mkdir -p hists/"+label+"/"+year_).c_str());

//trigger efficiency scale factors
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

bool isTTJets =  false;
if((input.find("ttbar") != std::string::npos) || (input.find("tt1L") != std::string::npos) ||  (input.find("tt2L") != std::string::npos)) isTTJets = true;

bool isHH = false;
if(outputFileName.find("HH") != std::string::npos) isHH = true;

bool isH = false;    
if(outputFileName.find("VH") != std::string::npos || outputFileName.find("ttH") != std::string::npos) isH = true;
   
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

TChain * chain = new TChain("Events");
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
    //histograms.addHistogram("fatJet1MassSD",   "; j_{1} soft drop mass (GeV); Events", 46,   40.,    500.,  [&]() { return hh.fatJet1MassSD();} );
    histograms.addHistogram("fatJet2MassRegressed",   "; j_{2} regressed mass (GeV); Events", 46,   40.,    500.,  [&]() { return hh.fatJet2MassRegressed();} );
    //histograms.addHistogram("fatJet1MassRegressed",   "; j_{1} regressed mass (GeV); Events", 46,   40.,    500.,  [&]() { return hh.fatJet1MassRegressed();} );
    //histograms.addHistogram("fatJet1PNetXbb",   "; j_{1} PNet Xbb tagger; Events",           200,   0.95, 1.0,   [&]() { return  hh.fatJet1PNetXbb(); } );
    //histograms.addHistogram("fatJet2PNetXbb",   "; j_{2} PNet Xbb tagger; Events",           200,   0.95, 1.0,   [&]() { return  hh.fatJet2PNetXbb(); } );

    if(input.find("Tau32TopSkim") != std::string::npos){
        histograms.addHistogram("fatJet1Pt",          "; p_{T}^{j1} (GeV); Events",           200,   300.,   2300.,  [&]() { return  hh.fatJet1Pt(); });
        histograms.addHistogram("fatJet2Pt",          "; p_{T}^{j2} (GeV); Events",           200,   300.,   2300.,  [&]() { return  hh.fatJet2Pt(); });
    }   
    
    //histograms.addHistogram("EventBDTv8p2",   "; Event BDT; Events",           200, 0.0, 1.0,   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(); } );
}

else
{
//final fit discriminant 
histograms.addHistogram("fatJet2MassSD",   "; j_{2} soft drop mass (GeV); Events", 46,   40.,    500.,  [&]() {
    if(syst_name.find("JMS_Down") != std::string::npos) return hh.fatJet2MassSD_JMS_Down();
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.fatJet2MassSD_JMS_Up();
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.fatJet2MassSD_JMR_Down();
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.fatJet2MassSD_JMR_Up();
    else return hh.fatJet2MassSD();
});
    
histograms.addHistogram("fatJet2MassRegressed",   "; j_{2} regressed mass (GeV); Events", 46,   40.,    500.,  [&]() {
    if(syst_name.find("JMS_Down") != std::string::npos) return hh.fatJet2MassRegressed_JMS_Down();
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.fatJet2MassRegressed_JMS_Up();
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.fatJet2MassRegressed_JMR_Down();
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.fatJet2MassRegressed_JMR_Up();
    else return hh.fatJet2MassRegressed();
});
    
histograms.addHistogram("fatJet1MassSD",   "; j_{1} soft drop mass (GeV); Events", 46,   40.,    500.,  [&]() {
    if(syst_name.find("JMS_Down") != std::string::npos) return hh.fatJet1MassSD_JMS_Down();
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.fatJet1MassSD_JMS_Up();
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.fatJet1MassSD_JMR_Down();
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.fatJet1MassSD_JMR_Up();
    else return hh.fatJet1MassSD();
});
    
histograms.addHistogram("fatJet1MassRegressed",   "; j_{1} regressed mass (GeV); Events", 46,   40.,    500.,  [&]() {
    if(syst_name.find("JMS_Down") != std::string::npos) return hh.fatJet1MassRegressed_JMS_Down();
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.fatJet1MassRegressed_JMS_Up();
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.fatJet1MassRegressed_JMR_Down();
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.fatJet1MassRegressed_JMR_Up();
    else return hh.fatJet1MassRegressed();
});

histograms.addHistogram("fatJet1Pt",          "; p_{T}^{j1} (GeV); Events",           200,   300.,   2300.,  [&]() { 
    if(syst_name.find("JES_Down") != std::string::npos) return hh.fatJet1Pt_JESDown();
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.fatJet1Pt_JESUp();
    else if(syst_name.find("JER_Down") != std::string::npos) return hh.fatJet1Pt_JERDown();
    else if(syst_name.find("JER_Up") != std::string::npos) return hh.fatJet1Pt_JERUp();
    else return  hh.fatJet1Pt(); 
});

histograms.addHistogram("fatJet2Pt",          "; p_{T}^{j2} (GeV); Events",           200,   300.,   2300.,  [&]() { 
    if(syst_name.find("JES_Down") != std::string::npos) return hh.fatJet2Pt_JESDown();
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.fatJet2Pt_JESUp();
    else if(syst_name.find("JER_Down") != std::string::npos) return hh.fatJet2Pt_JERDown();
    else if(syst_name.find("JER_Up") != std::string::npos) return hh.fatJet2Pt_JERUp();
    else return  hh.fatJet2Pt(); 
});

//other variables does not have correct JMS/JMR variations    
histograms.add2DHistogram("fat Jet2 MassSD vs pT", "mj2", 30,   50.,   200., "ptj2",  25,   250.,   750.,  [&]() { return  hh.fatJet2MassSD(); }, [&]() { return hh.fatJet2Pt();} );
histograms.add2DHistogram("fat Jet1 MassSD vs pT", "mj1", 30,   50.,   200., "ptj1",  25,   250.,   750.,  [&]() { return  hh.fatJet1MassSD(); }, [&]() { return hh.fatJet1Pt();} );
histograms.add2DHistogram("fat Jet2 Mass regressed vs pT", "mregj2", 30,   50.,   200., "ptj2",  25,   250.,   750.,  [&]() { return  hh.fatJet2MassRegressed(); }, [&]() { return hh.fatJet2Pt();} );
histograms.add2DHistogram("fat Jet1 Mass regressed vs pT", "mregj1", 30,   50.,   200., "ptj1",  25,   250.,   750.,  [&]() { return  hh.fatJet1MassRegressed(); }, [&]() { return hh.fatJet1Pt();} );
histograms.add2DHistogram("fat Jet1 pT vs Jet2 pT", "ptj1", 72,   200.,   2000., "ptj2",  72,   200.,   2000.,  [&]() { return  hh.fatJet1Pt(); }, [&]() { return hh.fatJet2Pt();} );

histograms.addHistogram("fatJet1MassSDLR",   "; j_{1} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet1MassSD(); } );
 histograms.addHistogram("fatJet1MassRegressedLR",   "; j_{1} regressed mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet1MassRegressed(); } );
//histograms.addHistogram("fatJet2MassSD",   "; j_{2} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet2MassSD(); } );
histograms.addHistogram("fatJet1MassSD_raw",   "; j_{1} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet1MassSD_UnCorrected(); } );
histograms.addHistogram("fatJet2MassSD_raw",   "; j_{2} soft drop mass (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet2MassSD_UnCorrected(); } );
histograms.addHistogram("fatJet1MassRegressed_raw",   "; j_{1} regressed mass raw (GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet1MassRegressed_UnCorrected(); } );
histograms.addHistogram("fatJet2MassRegressed_raw",   "; j_{2} regressed mass raw(GeV); Events", 300,   0.,   300.,  [&]() { return  hh.fatJet2MassRegressed_UnCorrected(); } );
histograms.addHistogram("MET",           "; p_{T}^{miss} (GeV); Events",         200,   0.,   500.,  [&]() { return hh.met(); } );
 histograms.addHistogram("MET_phi ",           "; p_{T}^{miss} #phi; Events",         200,   -3.2,   3.2,  [&]() { return hh.metphi(); } );
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
histograms.addHistogram("ptj1_over_mregj1",       "; p_{T}^{j1}/Mreg_{j1}; Events",         200,   0.,   10.,   [&]() { return  hh.fatJet1PtOverMRegressed(); } );
histograms.addHistogram("ptj2_over_mregj2",       "; p_{T}^{j2}/Mreg_{j2}; Events",         200,   0.,   10.,   [&]() { return  hh.fatJet2PtOverMRegressed(); } );
histograms.addHistogram("ptj2_over_ptj1",      "; p_{T}^{j2}/p_{T}^{j1}; Events",     200,   0.5,  1.,    [&]() { return  hh.fatJet2Pt() / hh.fatJet1Pt(); } );
histograms.addHistogram("mj2_over_mj1",      "; m^{j2}/m^{j1}; Events",               200,   0.0,  1.5,   [&]() { return  hh.fatJet2MassSD() / hh.fatJet1MassSD(); } );
histograms.addHistogram("mregj2_over_mregj1",      "; Mreg^{j2}/Mreg^{j1}; Events",               200,   0.0,  1.5,   [&]() { return  hh.fatJet2MassRegressed() / hh.fatJet1MassRegressed(); } );
histograms.addHistogram("fatJet1Tau3OverTau2",   "; j_{1} Tau3/2; Events",           200,   0.0,  1.0,   [&]() { return  hh.fatJet1Tau3OverTau2(); } );
histograms.addHistogram("fatJet2Tau3OverTau2",   "; j_{2} Tau3/2; Events",           200,   0.0,  1.0,   [&]() { return  hh.fatJet2Tau3OverTau2(); } );
histograms.addHistogram("lep1Eta",          "; #eta^{j1}; Events",                 200,   -2.5,  2.5,  [&]() { return  hh.lep1Eta(); } );
histograms.addHistogram("lep1Phi",          "; #Phi^{j1}; Events",                 200,  -3.2,   3.2,  [&]() { return  hh.lep1Phi(); } );
 histograms.addHistogram("gen_mHH",          "; gen m_{HH}; Events",                 200, 0.,   1500.,[&]() {  TLorentzVector gh1,gh2;
	  gh1.SetPtEtaPhiM(hh.genHiggs1Pt(), hh.genHiggs1Eta(),hh.genHiggs1Phi(),125.0);
	  gh2.SetPtEtaPhiM(hh.genHiggs2Pt(), hh.genHiggs2Eta(),hh.genHiggs2Phi(),125.0);
	  return (gh1+gh2).M();} );
if(input.find("1LTopSkim") != std::string::npos) histograms.addHistogram("abs_dR_l1j1",       "; #DeltaR(l_{1}, j_{1}); Events",        200,   0.,   5.0,    [&]() { return  sqrt(pow(hh.lep1Eta() - hh.fatJet1Eta(), 2) + pow(fabs(hh.lep1Phi()-hh.fatJet1Phi()) > pi ? fabs(hh.lep1Phi()-hh.fatJet1Phi()) - 2*pi : hh.lep1Phi()-hh.fatJet1Phi(), 2)); } );
else
{
    //histograms.addHistogram("EventBDTNoBinv24",   "; Event BDT; Events",           1000,   0.0,  0.2,   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v24(); } );
    //histograms.addHistogram("EventBDTv24",   "; Event BDT; Events",           {0.000, 0.0024, 0.0054, 0.0118, 0.10},   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v24(); } );
    histograms.addHistogram("EventBDTNoBinv8p2",   "; Event BDT; Events",           1000,   0.0,  0.2,   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(); } );
    histograms.addHistogram("EventBDTv8p2v3",   "; Event BDT; Events",           {0.000, 0.00008, 0.0002, 0.0004, 0.0005, 0.0007, 0.01},   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(); } );
    histograms.addHistogram("EventBDTv8p2v2",   "; Event BDT; Events",           {0.000, 0.00008, 0.0002, 0.0004, 0.0005, 0.01},   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(); } );
    histograms.addHistogram("EventBDTv8p2v1",   "; Event BDT; Events",           {0.000, 0.00008, 0.0002, 0.0004, 1.0},   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(); } );
    histograms.addHistogram("EventBDTv8p2",   "; Event BDT; Events",           {0.000, 0.00008, 0.0002, 0.0004, 0.01},   [&]() { return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(); } );
 }
 
 }
    
//************define cuts**********//

cutflow.setTFile(outfile);

if(input.find("1LTopSkim") != std::string::npos) // this is 1LTopSkim input
{ 
//Pre-selection cuts
  cutflow.addCut("CutWeight", [&](){ return 1; },   [&](){ return isData ?  lumi : lumi * hh.weight() *hh.l1PreFiringWeight()*hh.xsecWeight()* hh.puWeight() * hh.genWeight() * (abs(hh.lep1Id()) == 11 ? miniIsoEle_sf.getminiIsoScaleFactors(hh.lep1Pt(),hh.lep1Eta()) : miniIsoMu_sf.getminiIsoScaleFactors(hh.lep1Pt(),hh.lep1Eta())) * (abs(hh.lep1Id()) == 11 ? elTrig_sf.getScaleFactors(hh.lep1Pt(),hh.lep1Eta()) : muTrig_sf.getTrigScaleFactors(hh.lep1Pt(),hh.lep1Eta(), year_)) * (abs(hh.lep1Id()) == 11 ? elID_sf.getScaleFactors(hh.lep1Pt(),hh.lep1Eta()) : muID_sf.getIDScaleFactors(hh.lep1Pt(),hh.lep1Eta(), year_)) ;}); //before correction
//cutflow.addCut("CutWeight", [&](){ return 1; },   [&](){ return isData ?  lumi : lumi * hh.weight() * hh.puWeight() * (isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0) : 1.0); });//after correction
  if(input.find("HHc1") == std::string::npos) cutflow.addCutToLastActiveCut("CutHLT",       [&](){ return abs(hh.lep1Id()) == 11 ? (hh.HLT_Ele27_WPTight_Gsf() || hh.HLT_Ele32_WPTight_Gsf() || hh.HLT_Ele32_WPTight_Gsf_L1DoubleEG()) : (hh.HLT_IsoMu24()  ||  hh.HLT_IsoMu27()  || hh.HLT_Mu50() ); },   UNITY);
cutflow.addCutToLastActiveCut("CutLepJetPt",       [&](){ return hh.fatJet1Pt() > 300.0 && hh.lep1Pt() > 50.0 && hh.lep2Pt() <=0 ; },   UNITY);
cutflow.addCutToLastActiveCut("CutfatJetMassSD",       [&](){ return hh.fatJet1MassSD() > 50.0; },   UNITY);
//ttbar 1L+jet control region
cutflow.addCutToLastActiveCut("CutLepEta",       [&](){ return (abs(hh.lep1Id()) == 11 && fabs(hh.lep1Eta()) <  2.5) || (abs(hh.lep1Id()) == 13 && fabs(hh.lep1Eta()) <  2.4); },   UNITY);
cutflow.addCutToLastActiveCut("CutMET",       [&](){ return hh.met() > 100.0; },   UNITY);
//cutflow.addCutToLastActiveCut("CutHEM2018",       [&](){ return ((year_ == "2018" && isData && hh.run() >=319077) ? !((abs(hh.lep1Id()) == 11 && hh.lep1Phi() > -1.57 && hh.lep1Phi() < -0.87 && hh.lep1Eta() < -1.3) || (hh.fatJet1Phi() > -1.57 && hh.fatJet1Phi() < -0.87 && hh.fatJet1Eta()< -1.3 ) ): true); },   UNITY);
cutflow.addCutToLastActiveCut("CutTau3Over2",       [&](){ return hh.fatJet1MassSD() > 140.0 && hh.fatJet1Tau3OverTau2() < 0.46 ; },   [&]() {return isTTJets ? TopTagSF("0.46", year_, hh.fatJet1Pt() ) : 1.0;} );
cutflow.addCutToLastActiveCut("TTBarLepJetCR",       [&](){ return sqrt(pow(hh.lep1Eta() - hh.fatJet1Eta(), 2) + pow(fabs(hh.lep1Phi()-hh.fatJet1Phi()) > pi ? fabs(hh.lep1Phi()-hh.fatJet1Phi()) - 2*pi : hh.lep1Phi()-hh.fatJet1Phi(), 2)) > 1.0; },   UNITY); //delta R > 1.0
cutflow.addCutToLastActiveCut("TTBarLepJetCRElectron",       [&](){ return abs(hh.lep1Id()) == 11; },   UNITY);
cutflow.getCut("TTBarLepJetCR");
cutflow.addCutToLastActiveCut("TTBarLepJetCRMuon",       [&](){ return abs(hh.lep1Id()) == 13; },   UNITY);
}

//two jet ttbar CR and SR
else
{

//Pre-selection cuts
//cutflow.addCut("CutWeight", [&](){ return 1; },   [&](){ return isData ?  lumi : lumi*hh.weight()*hh.xsecWeight()*hh.l1PreFiringWeight()*hh.triggerEff3DWeight()*hh.puWeight() * (isTTJets  ? ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 0) : 1.0); }); //after correction
    
//use the trigger weight directly from the ntuple, verified to be the same as running the function getTrigEffEvt()
//cutflow.addCut("CutWeight", [&](){ return 1; },   [&](){ return isData ?  lumi : lumi*hh.weight()*hh.xsecWeight()*hh.l1PreFiringWeight()*hh.triggerEff3DWeight()*hh.puWeight(); }); //before correction

//Pre-selection cuts
  cutflow.addCut("CutWeight", [&](){ return 1; },  [&](){
      //after ttbar recoil correction
      float ttbar_factor = input.find("option5") == std::string::npos ? 1.0 : 2.0;
     
      float total_weight = isData ?  lumi :lumi * hh.l1PreFiringWeight() * hh.puWeight() * hh.xsecWeight() * (isHH? hh.weight() : hh.genWeight()) * (isTTJets ? ttbar_factor*ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 0) : 1.0);

      //If applying correction from VBF analysis
      //float total_weight = isData ?  lumi :lumi * hh.l1PreFiringWeight() * hh.puWeight() * hh.xsecWeight() * (isHH? hh.weight() : hh.genWeight());
      //if(isTTJets){
      //	if (hh.fatJet1Pt()>hh.fatJet2Pt())total_weight *= ttjets_sf.getScaleFactors_VBF(year_, hh.fatJet1Pt()) ;
      //	else total_weight *= ttjets_sf.getScaleFactors_VBF(year_, hh.fatJet2Pt()) ;
      //}

    //before ttbar recoil correction
      //float total_weight = isData ?  lumi :lumi * hh.l1PreFiringWeight() * hh.puWeight() * hh.xsecWeight() * (isHH? hh.weight() : hh.genWeight()* (isTTJets ? ttbar_factor:1.0));
    
      
    if(!isData){

    //apply trigger SF
    if(syst_name.find("JER_Up") != std::string::npos){
      total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt_JERUp(), hh.fatJet1MassSD(), hh.fatJet1PNetXbb(), hh.fatJet2Pt_JERUp(), hh.fatJet2MassSD(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    }
    else if(syst_name.find("JER_Down") != std::string::npos){
      total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt_JERDown(), hh.fatJet1MassSD(), hh.fatJet1PNetXbb(), hh.fatJet2Pt_JERDown(), hh.fatJet2MassSD(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    } 
    else if(syst_name.find("JES_Up") != std::string::npos){
        total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt_JESUp(), hh.fatJet1MassSD(), hh.fatJet1PNetXbb(), hh.fatJet2Pt_JESUp(), hh.fatJet2MassSD(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    }
    else if(syst_name.find("JES_Down") != std::string::npos){
        total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt_JESDown(), hh.fatJet1MassSD(), hh.fatJet1PNetXbb(), hh.fatJet2Pt_JESDown(), hh.fatJet2MassSD(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    } 
    else if(syst_name.find("JMS_Down") != std::string::npos){
        total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt(), hh.fatJet1MassSD_JMS_Down(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2MassSD_JMS_Down(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    }
    else if(syst_name.find("JMS_Up") != std::string::npos){
        total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt(), hh.fatJet1MassSD_JMS_Up(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2MassSD_JMS_Up(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    }
    else if(syst_name.find("JMR_Down") != std::string::npos){
        total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt(), hh.fatJet1MassSD_JMR_Down(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2MassSD_JMR_Down(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    }
    else if(syst_name.find("JMR_Up") != std::string::npos){
        total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt(), hh.fatJet1MassSD_JMR_Up(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2MassSD_JMR_Up(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");
    }
    else{
      total_weight = total_weight * trig_sf.getTrigEffEvt(hh.fatJet1Pt(), hh.fatJet1MassSD(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2MassSD(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");        
    }
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
    
    }
    
    return total_weight;
}); 
 
cutflow.addCutToLastActiveCut("CutHLT",       [&](){ 
   //cout <<"test "<<hh.HLT_AK8PFJet330_PFAK8BTagCSV_p17()<<endl;
   return isData ? ((year_ == "2016" && (hh.HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20() || hh.HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20() || hh.HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20() || hh.HLT_AK8PFJet360_TrimMass30() || hh.HLT_AK8PFJet450() || hh.HLT_PFJet450() )) || (year_ == "2017" && (hh.HLT_PFJet450() || hh.HLT_PFJet500() || hh.HLT_AK8PFJet500() || hh.HLT_PFHT1050() || hh.HLT_AK8PFJet360_TrimMass30() || hh.HLT_AK8PFJet380_TrimMass30() || hh.HLT_AK8PFJet400_TrimMass30() || hh.HLT_AK8PFHT800_TrimMass50() || hh.HLT_AK8PFHT750_TrimMass50() || hh.HLT_AK8PFJet330_PFAK8BTagCSV_p17())) || (year_ == "2018" && (hh.HLT_PFHT1050() || hh.HLT_PFJet500() || hh.HLT_AK8PFJet500() || hh.HLT_AK8PFJet400_TrimMass30() || hh.HLT_AK8PFHT800_TrimMass50() || hh.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4()))) : 1.0; },   UNITY);
    

//cutflow.addCutToLastActiveCut("CutfatJetsPt",       [&](){ return hh.fatJet1Pt() > 300.0 && hh.fatJet2Pt() > 300.0; },   UNITY);
cutflow.addCutToLastActiveCut("CutfatJetsPt",       [&](){
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.isVBFtag() < 1 && hh.fatJet1Pt() > 300.0 && hh.fatJet2Pt() > 300.0;
    else if(syst_name.find("JER_Up") != std::string::npos) return hh.isVBFtag() < 1 && hh.fatJet1Pt_JERUp() > 300.0 && hh.fatJet2Pt_JERUp() > 300.0; 
    else if(syst_name.find("JER_Down") != std::string::npos) return hh.isVBFtag() < 1 && hh.fatJet1Pt_JERDown() > 300.0 && hh.fatJet2Pt_JERDown() > 300.0; 
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.isVBFtag() < 1 && hh.fatJet1Pt_JESUp() > 300.0 && hh.fatJet2Pt_JESUp() > 300.0; 
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.isVBFtag() < 1 && hh.fatJet1Pt_JESDown() > 300.0 && hh.fatJet2Pt_JESDown() > 300.0; 
    else return hh.isVBFtag() < 1 && hh.fatJet1Pt() > 300.0 && hh.fatJet2Pt() > 300.0;
},   UNITY);
cutflow.addCutToLastActiveCut("CutfatJetsMassSD",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.fatJet1MassSD() > 50.0 && hh.fatJet2MassSD() > 50.0; 
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.fatJet1MassSD_JMS_Up() > 50.0 && hh.fatJet2MassSD_JMS_Up() > 50.0; 
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.fatJet1MassSD_JMS_Down() > 50.0 && hh.fatJet2MassSD_JMS_Down() > 50.0; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.fatJet1MassSD_JMR_Up() > 50.0 && hh.fatJet2MassSD_JMR_Up() > 50.0; 
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.fatJet1MassSD_JMR_Down() > 50.0 && hh.fatJet2MassSD_JMR_Down() > 50.0;
    else return hh.fatJet1MassSD() > 50.0 && hh.fatJet2MassSD() > 50.0;
},   UNITY);

  
//Signal regions - pass - BDT v8p2
if(input.find("Tau32TopSkim") == std::string::npos){

cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1",       [&](){    
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.980; 
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.980; 
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.980; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.980; 
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.980;
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.980;
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.980;
    else return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.980; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
//cutflow.addCutToLastActiveCut("J2MassSideBandv8p2Bin1",       [&](){ return hh.fatJet2MassSD() <= 110 || hh.fatJet2MassSD() >= 140.0; },   UNITY);


cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin2",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.950)); 
    else if(syst_name.find("JMS_Up") != std::string::npos) return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.950)); 
    else if(syst_name.find("JMS_Down") != std::string::npos) return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.950));
    else if(syst_name.find("JMR_Up") != std::string::npos) return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.950));
    else if(syst_name.find("JMR_Down") != std::string::npos) return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.950)); 
    else if(syst_name.find("JES_Up") != std::string::npos) return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.950));
    else if(syst_name.find("JES_Down") != std::string::npos) return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.950));
    else return (!(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.980))  && ((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.950)); 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
//cutflow.addCutToLastActiveCut("J2MassSideBandv8p2Bin2",       [&](){ return hh.fatJet2MassSD() <= 110 || hh.fatJet2MassSD() >= 140.0; },   UNITY);

cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin3",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
    else if(syst_name.find("JMS_Up") != std::string::npos) return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
    else if(syst_name.find("JMS_Down") != std::string::npos) return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
    else if(syst_name.find("JMR_Down") != std::string::npos) return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
    else if(syst_name.find("JES_Up") != std::string::npos) return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
    else if(syst_name.find("JES_Down") != std::string::npos) return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
    else return (!((hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.11 && hh.fatJet2PNetXbb() > 0.980)||(hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.950)))  && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.03 && hh.fatJet2PNetXbb() > 0.950; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
//cutflow.addCutToLastActiveCut("J2MassSideBandv8p2Bin3",       [&](){ return hh.fatJet2MassSD() <= 110 || hh.fatJet2MassSD() >= 140.0; },   UNITY);

cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1PNetp95",       [&](){    
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.95; 
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.95; 
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.95; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.95; 
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.95;
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.95;
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.95;
    else return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.95; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
    
cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1PNetp92",       [&](){    
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.92; 
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.92; 
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.92; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.92; 
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.92;
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.92;
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.92;
    else return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.92; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });

    cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1PNetp9",       [&](){
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.9;
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.9;
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.9; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.9;
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.9;
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.9;
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.9;
    else return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.9; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
    
    cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1PNetp5",       [&](){
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.5;
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.5;
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.5; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.5;
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.5;
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.5;
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.5;
    else return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.5; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
    
    cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1PNetp2",       [&](){
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.2;
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.2;
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.2; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.2;
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.2;
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.2;
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.2;
    else return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.2; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
    
    cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("SRv8p2Bin1PNetp0",       [&](){
    if(isData || syst_name.find("nominal") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.0;
    else if(syst_name.find("JMS_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.43 && hh.fatJet2PNetXbb() > 0.0;
    else if(syst_name.find("JMS_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.43 && hh.fatJet2PNetXbb() > 0.0; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.43 && hh.fatJet2PNetXbb() > 0.0;
    else if(syst_name.find("JMR_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.43 && hh.fatJet2PNetXbb() > 0.0;
    else if(syst_name.find("JES_Up") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.43 && hh.fatJet2PNetXbb() > 0.0;
    else if(syst_name.find("JES_Down") != std::string::npos) return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.43 && hh.fatJet2PNetXbb() > 0.0;
    else return hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.43 && hh.fatJet2PNetXbb() > 0.0; 
},   [&](){ return isTTJets  ? ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0): 1.0; });
    
cutflow.getCut("CutfatJetsMassSD");
cutflow.addCutToLastActiveCut("FailSRv8p2",       [&](){ 
    if(isData || syst_name.find("nominal") != std::string::npos) return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
    else if(syst_name.find("JMS_Up") != std::string::npos) return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
    else if(syst_name.find("JMS_Down") != std::string::npos) return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
    else if(syst_name.find("JMR_Up") != std::string::npos) return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
    else if(syst_name.find("JMR_Down") != std::string::npos) return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
    else if(syst_name.find("JES_Up") != std::string::npos) return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
    else if(syst_name.find("JES_Down") != std::string::npos) return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
    else return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.03 && hh.fatJet2PNetXbb() < 0.950; 
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
	return hh.fatJet1Tau3OverTau2() < 0.46 && hh.fatJet2Tau3OverTau2() < 0.46 &&  hh.fatJet1HasBJetCSVLoose() && hh.fatJet2HasBJetCSVLoose();},   [&]() {// && hh.fatJet1Pt() >450. && hh.fatJet2Pt() > 450.;},   [&]() {
   
    if(syst_name.find("JES_Up") != std::string::npos){
        return isTTJets ? (TopTagSF("0.46", year_, hh.fatJet1Pt_JESUp()) * TopTagSF("0.46", year_, hh.fatJet2Pt_JESUp())) : 1.0;
    }
    else if(syst_name.find("JES_Down") != std::string::npos){
        return isTTJets ? (TopTagSF("0.46", year_, hh.fatJet1Pt_JESDown()) * TopTagSF("0.46", year_, hh.fatJet2Pt_JESDown())) : 1.0;
    }
    else if(syst_name.find("JER_Up") != std::string::npos){
        return isTTJets ? (TopTagSF("0.46", year_, hh.fatJet1Pt_JERUp()) * TopTagSF("0.46", year_, hh.fatJet2Pt_JERUp())) : 1.0;
    }
    else if(syst_name.find("JER_Down") != std::string::npos){
        return isTTJets ? (TopTagSF("0.46", year_, hh.fatJet1Pt_JERDown()) * TopTagSF("0.46", year_, hh.fatJet2Pt_JERDown())) : 1.0;
    }
    else{
        return isTTJets ? (TopTagSF("0.46", year_, hh.fatJet1Pt()) * TopTagSF("0.46", year_, hh.fatJet2Pt())) : 1.0;
    } 
    } );

    cutflow.getCut("TTBarCR");
    cutflow.addCutToLastActiveCut("TTBarCRBDT1v8p2",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() < 0.00008 && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.000; },   UNITY);
    cutflow.getCut("TTBarCR");
    cutflow.addCutToLastActiveCut("TTBarCRBDT2v8p2",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() < 0.0002 && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.00008; },   UNITY);
    cutflow.getCut("TTBarCR");
    cutflow.addCutToLastActiveCut("TTBarCRBDT3v8p2",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() < 0.0004 && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.0002; },   UNITY);
    cutflow.getCut("TTBarCR");
    cutflow.addCutToLastActiveCut("TTBarCRBDT4v8p2",       [&](){ return  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() < 1.0 && hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2() > 0.0004; },   UNITY);

    //cutflow.getCut("CutfatJetsXbb");
    cutflow.getCut("CutfatJetsXbb");
    cutflow.addCutToLastActiveCut("TTBarCRTight",       [&](){ return hh.fatJet1Tau3OverTau2() < 0.39 && hh.fatJet2Tau3OverTau2() < 0.39 &&  hh.fatJet1HasBJetCSVLoose() && hh.fatJet2HasBJetCSVLoose(); },   [&]() {return isTTJets ? (TopTagSF("0.40", year_, hh.fatJet1Pt()) * TopTagSF("0.40", year_, hh.fatJet2Pt())) : 1.0;} );  
    }
}

/****Systematics******/
if(doSystematics && (outputFileName.find("qcd") == std::string::npos ) && (outputFileName.find("data") == std::string::npos ) && (syst_name.find("nominal") != std::string::npos) )
{    
    //pdf uncertainty for the HH signal acceptance
    if(isHH && (outputFileName.find("VBF") == std::string::npos)){
        
    cutflow.addWgtSyst("LHEPDFEigenv0",  [&](){return isHH ?  hh.LHEPdfWeight()[0] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv1",  [&](){return isHH ?  hh.LHEPdfWeight()[1] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv2",  [&](){return isHH ?  hh.LHEPdfWeight()[2] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv3",  [&](){return isHH ?  hh.LHEPdfWeight()[3] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv4",  [&](){return isHH ?  hh.LHEPdfWeight()[4] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv5",  [&](){return isHH ?  hh.LHEPdfWeight()[5] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv6",  [&](){return isHH ?  hh.LHEPdfWeight()[6] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv7",  [&](){return isHH ?  hh.LHEPdfWeight()[7] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv8",  [&](){return isHH ?  hh.LHEPdfWeight()[8] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv9",  [&](){return isHH ?  hh.LHEPdfWeight()[9] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv10",  [&](){return isHH ?  hh.LHEPdfWeight()[10] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv11",  [&](){return isHH ?  hh.LHEPdfWeight()[11] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv12",  [&](){return isHH ?  hh.LHEPdfWeight()[12] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv13",  [&](){return isHH ?  hh.LHEPdfWeight()[13] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv14",  [&](){return isHH ?  hh.LHEPdfWeight()[14] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv15",  [&](){return isHH ?  hh.LHEPdfWeight()[15] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv16",  [&](){return isHH ?  hh.LHEPdfWeight()[16] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv17",  [&](){return isHH ?  hh.LHEPdfWeight()[17] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv18",  [&](){return isHH ?  hh.LHEPdfWeight()[18] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv19",  [&](){return isHH ?  hh.LHEPdfWeight()[19] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv20",  [&](){return isHH ?  hh.LHEPdfWeight()[20] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv21",  [&](){return isHH ?  hh.LHEPdfWeight()[21] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv22",  [&](){return isHH ?  hh.LHEPdfWeight()[22] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv23",  [&](){return isHH ?  hh.LHEPdfWeight()[23] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv24",  [&](){return isHH ?  hh.LHEPdfWeight()[24] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv25",  [&](){return isHH ?  hh.LHEPdfWeight()[25] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv26",  [&](){return isHH ?  hh.LHEPdfWeight()[26] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv27",  [&](){return isHH ?  hh.LHEPdfWeight()[27] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv28",  [&](){return isHH ?  hh.LHEPdfWeight()[28] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv29",  [&](){return isHH ?  hh.LHEPdfWeight()[29] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv30",  [&](){return isHH ?  hh.LHEPdfWeight()[30] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv31",  [&](){return isHH ?  hh.LHEPdfWeight()[31] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv32",  [&](){return isHH ?  hh.LHEPdfWeight()[32] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv33",  [&](){return isHH ?  hh.LHEPdfWeight()[33] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv34",  [&](){return isHH ?  hh.LHEPdfWeight()[34] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv35",  [&](){return isHH ?  hh.LHEPdfWeight()[35] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv36",  [&](){return isHH ?  hh.LHEPdfWeight()[36] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv37",  [&](){return isHH ?  hh.LHEPdfWeight()[37] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv38",  [&](){return isHH ?  hh.LHEPdfWeight()[38] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv39",  [&](){return isHH ?  hh.LHEPdfWeight()[39] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv40",  [&](){return isHH ?  hh.LHEPdfWeight()[40] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv41",  [&](){return isHH ?  hh.LHEPdfWeight()[41] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv42",  [&](){return isHH ?  hh.LHEPdfWeight()[42] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv43",  [&](){return isHH ?  hh.LHEPdfWeight()[43] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv44",  [&](){return isHH ?  hh.LHEPdfWeight()[44] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv45",  [&](){return isHH ?  hh.LHEPdfWeight()[45] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv46",  [&](){return isHH ?  hh.LHEPdfWeight()[46] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv47",  [&](){return isHH ?  hh.LHEPdfWeight()[47] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv48",  [&](){return isHH ?  hh.LHEPdfWeight()[48] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv49",  [&](){return isHH ?  hh.LHEPdfWeight()[49] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv50",  [&](){return isHH ?  hh.LHEPdfWeight()[50] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv51",  [&](){return isHH ?  hh.LHEPdfWeight()[51] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv52",  [&](){return isHH ?  hh.LHEPdfWeight()[52] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv53",  [&](){return isHH ?  hh.LHEPdfWeight()[53] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv54",  [&](){return isHH ?  hh.LHEPdfWeight()[54] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv55",  [&](){return isHH ?  hh.LHEPdfWeight()[55] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv56",  [&](){return isHH ?  hh.LHEPdfWeight()[56] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv57",  [&](){return isHH ?  hh.LHEPdfWeight()[57] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv58",  [&](){return isHH ?  hh.LHEPdfWeight()[58] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv59",  [&](){return isHH ?  hh.LHEPdfWeight()[59] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv60",  [&](){return isHH ?  hh.LHEPdfWeight()[60] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv61",  [&](){return isHH ?  hh.LHEPdfWeight()[61] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv62",  [&](){return isHH ?  hh.LHEPdfWeight()[62] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv63",  [&](){return isHH ?  hh.LHEPdfWeight()[63] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv64",  [&](){return isHH ?  hh.LHEPdfWeight()[64] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv65",  [&](){return isHH ?  hh.LHEPdfWeight()[65] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv66",  [&](){return isHH ?  hh.LHEPdfWeight()[66] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv67",  [&](){return isHH ?  hh.LHEPdfWeight()[67] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv68",  [&](){return isHH ?  hh.LHEPdfWeight()[68] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv69",  [&](){return isHH ?  hh.LHEPdfWeight()[69] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv70",  [&](){return isHH ?  hh.LHEPdfWeight()[70] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv71",  [&](){return isHH ?  hh.LHEPdfWeight()[71] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv72",  [&](){return isHH ?  hh.LHEPdfWeight()[72] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv73",  [&](){return isHH ?  hh.LHEPdfWeight()[73] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv74",  [&](){return isHH ?  hh.LHEPdfWeight()[74] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv75",  [&](){return isHH ?  hh.LHEPdfWeight()[75] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv76",  [&](){return isHH ?  hh.LHEPdfWeight()[76] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv77",  [&](){return isHH ?  hh.LHEPdfWeight()[77] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv78",  [&](){return isHH ?  hh.LHEPdfWeight()[78] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv79",  [&](){return isHH ?  hh.LHEPdfWeight()[79] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv80",  [&](){return isHH ?  hh.LHEPdfWeight()[80] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv81",  [&](){return isHH ?  hh.LHEPdfWeight()[81] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv82",  [&](){return isHH ?  hh.LHEPdfWeight()[82] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv83",  [&](){return isHH ?  hh.LHEPdfWeight()[83] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv84",  [&](){return isHH ?  hh.LHEPdfWeight()[84] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv85",  [&](){return isHH ?  hh.LHEPdfWeight()[85] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv86",  [&](){return isHH ?  hh.LHEPdfWeight()[86] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv87",  [&](){return isHH ?  hh.LHEPdfWeight()[87] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv88",  [&](){return isHH ?  hh.LHEPdfWeight()[88] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv89",  [&](){return isHH ?  hh.LHEPdfWeight()[89] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv90",  [&](){return isHH ?  hh.LHEPdfWeight()[90] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv91",  [&](){return isHH ?  hh.LHEPdfWeight()[91] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv92",  [&](){return isHH ?  hh.LHEPdfWeight()[92] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv93",  [&](){return isHH ?  hh.LHEPdfWeight()[93] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv94",  [&](){return isHH ?  hh.LHEPdfWeight()[94] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv95",  [&](){return isHH ?  hh.LHEPdfWeight()[95] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv96",  [&](){return isHH ?  hh.LHEPdfWeight()[96] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv97",  [&](){return isHH ?  hh.LHEPdfWeight()[97] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv98",  [&](){return isHH ?  hh.LHEPdfWeight()[98] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv99",  [&](){return isHH ?  hh.LHEPdfWeight()[99] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv100",  [&](){return isHH ?  hh.LHEPdfWeight()[100] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv101",  [&](){return isHH ?  hh.LHEPdfWeight()[101] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv102",  [&](){return isHH ?  hh.LHEPdfWeight()[102] : 1.0;});

    /*
    cutflow.addWgtSyst("LHEPDFEigenv0",  [&](){return isHH ?  hh.LHEPdfWeight()[0]*hh.LHEPdfWeightNorm()[0] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv1",  [&](){return isHH ?  hh.LHEPdfWeight()[1]*hh.LHEPdfWeightNorm()[1] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv2",  [&](){return isHH ?  hh.LHEPdfWeight()[2]*hh.LHEPdfWeightNorm()[2] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv3",  [&](){return isHH ?  hh.LHEPdfWeight()[3]*hh.LHEPdfWeightNorm()[3] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv4",  [&](){return isHH ?  hh.LHEPdfWeight()[4]*hh.LHEPdfWeightNorm()[4] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv5",  [&](){return isHH ?  hh.LHEPdfWeight()[5]*hh.LHEPdfWeightNorm()[5] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv6",  [&](){return isHH ?  hh.LHEPdfWeight()[6]*hh.LHEPdfWeightNorm()[6] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv7",  [&](){return isHH ?  hh.LHEPdfWeight()[7]*hh.LHEPdfWeightNorm()[7] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv8",  [&](){return isHH ?  hh.LHEPdfWeight()[8]*hh.LHEPdfWeightNorm()[8] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv9",  [&](){return isHH ?  hh.LHEPdfWeight()[9]*hh.LHEPdfWeightNorm()[9] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv10",  [&](){return isHH ?  hh.LHEPdfWeight()[10]*hh.LHEPdfWeightNorm()[10] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv11",  [&](){return isHH ?  hh.LHEPdfWeight()[11]*hh.LHEPdfWeightNorm()[11] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv12",  [&](){return isHH ?  hh.LHEPdfWeight()[12]*hh.LHEPdfWeightNorm()[12] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv13",  [&](){return isHH ?  hh.LHEPdfWeight()[13]*hh.LHEPdfWeightNorm()[13] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv14",  [&](){return isHH ?  hh.LHEPdfWeight()[14]*hh.LHEPdfWeightNorm()[14] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv15",  [&](){return isHH ?  hh.LHEPdfWeight()[15]*hh.LHEPdfWeightNorm()[15] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv16",  [&](){return isHH ?  hh.LHEPdfWeight()[16]*hh.LHEPdfWeightNorm()[16] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv17",  [&](){return isHH ?  hh.LHEPdfWeight()[17]*hh.LHEPdfWeightNorm()[17] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv18",  [&](){return isHH ?  hh.LHEPdfWeight()[18]*hh.LHEPdfWeightNorm()[18] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv19",  [&](){return isHH ?  hh.LHEPdfWeight()[19]*hh.LHEPdfWeightNorm()[19] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv20",  [&](){return isHH ?  hh.LHEPdfWeight()[20]*hh.LHEPdfWeightNorm()[20] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv21",  [&](){return isHH ?  hh.LHEPdfWeight()[21]*hh.LHEPdfWeightNorm()[21] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv22",  [&](){return isHH ?  hh.LHEPdfWeight()[22]*hh.LHEPdfWeightNorm()[22] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv23",  [&](){return isHH ?  hh.LHEPdfWeight()[23]*hh.LHEPdfWeightNorm()[23] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv24",  [&](){return isHH ?  hh.LHEPdfWeight()[24]*hh.LHEPdfWeightNorm()[24] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv25",  [&](){return isHH ?  hh.LHEPdfWeight()[25]*hh.LHEPdfWeightNorm()[25] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv26",  [&](){return isHH ?  hh.LHEPdfWeight()[26]*hh.LHEPdfWeightNorm()[26] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv27",  [&](){return isHH ?  hh.LHEPdfWeight()[27]*hh.LHEPdfWeightNorm()[27] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv28",  [&](){return isHH ?  hh.LHEPdfWeight()[28]*hh.LHEPdfWeightNorm()[28] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv29",  [&](){return isHH ?  hh.LHEPdfWeight()[29]*hh.LHEPdfWeightNorm()[29] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv30",  [&](){return isHH ?  hh.LHEPdfWeight()[30]*hh.LHEPdfWeightNorm()[30] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv31",  [&](){return isHH ?  hh.LHEPdfWeight()[31]*hh.LHEPdfWeightNorm()[31] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv32",  [&](){return isHH ?  hh.LHEPdfWeight()[32]*hh.LHEPdfWeightNorm()[32] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv33",  [&](){return isHH ?  hh.LHEPdfWeight()[33]*hh.LHEPdfWeightNorm()[33] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv34",  [&](){return isHH ?  hh.LHEPdfWeight()[34]*hh.LHEPdfWeightNorm()[34] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv35",  [&](){return isHH ?  hh.LHEPdfWeight()[35]*hh.LHEPdfWeightNorm()[35] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv36",  [&](){return isHH ?  hh.LHEPdfWeight()[36]*hh.LHEPdfWeightNorm()[36] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv37",  [&](){return isHH ?  hh.LHEPdfWeight()[37]*hh.LHEPdfWeightNorm()[37] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv38",  [&](){return isHH ?  hh.LHEPdfWeight()[38]*hh.LHEPdfWeightNorm()[38] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv39",  [&](){return isHH ?  hh.LHEPdfWeight()[39]*hh.LHEPdfWeightNorm()[39] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv40",  [&](){return isHH ?  hh.LHEPdfWeight()[40]*hh.LHEPdfWeightNorm()[40] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv41",  [&](){return isHH ?  hh.LHEPdfWeight()[41]*hh.LHEPdfWeightNorm()[41] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv42",  [&](){return isHH ?  hh.LHEPdfWeight()[42]*hh.LHEPdfWeightNorm()[42] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv43",  [&](){return isHH ?  hh.LHEPdfWeight()[43]*hh.LHEPdfWeightNorm()[43] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv44",  [&](){return isHH ?  hh.LHEPdfWeight()[44]*hh.LHEPdfWeightNorm()[44] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv45",  [&](){return isHH ?  hh.LHEPdfWeight()[45]*hh.LHEPdfWeightNorm()[45] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv46",  [&](){return isHH ?  hh.LHEPdfWeight()[46]*hh.LHEPdfWeightNorm()[46] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv47",  [&](){return isHH ?  hh.LHEPdfWeight()[47]*hh.LHEPdfWeightNorm()[47] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv48",  [&](){return isHH ?  hh.LHEPdfWeight()[48]*hh.LHEPdfWeightNorm()[48] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv49",  [&](){return isHH ?  hh.LHEPdfWeight()[49]*hh.LHEPdfWeightNorm()[49] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv50",  [&](){return isHH ?  hh.LHEPdfWeight()[50]*hh.LHEPdfWeightNorm()[50] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv51",  [&](){return isHH ?  hh.LHEPdfWeight()[51]*hh.LHEPdfWeightNorm()[51] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv52",  [&](){return isHH ?  hh.LHEPdfWeight()[52]*hh.LHEPdfWeightNorm()[52] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv53",  [&](){return isHH ?  hh.LHEPdfWeight()[53]*hh.LHEPdfWeightNorm()[53] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv54",  [&](){return isHH ?  hh.LHEPdfWeight()[54]*hh.LHEPdfWeightNorm()[54] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv55",  [&](){return isHH ?  hh.LHEPdfWeight()[55]*hh.LHEPdfWeightNorm()[55] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv56",  [&](){return isHH ?  hh.LHEPdfWeight()[56]*hh.LHEPdfWeightNorm()[56] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv57",  [&](){return isHH ?  hh.LHEPdfWeight()[57]*hh.LHEPdfWeightNorm()[57] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv58",  [&](){return isHH ?  hh.LHEPdfWeight()[58]*hh.LHEPdfWeightNorm()[58] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv59",  [&](){return isHH ?  hh.LHEPdfWeight()[59]*hh.LHEPdfWeightNorm()[59] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv60",  [&](){return isHH ?  hh.LHEPdfWeight()[60]*hh.LHEPdfWeightNorm()[60] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv61",  [&](){return isHH ?  hh.LHEPdfWeight()[61]*hh.LHEPdfWeightNorm()[61] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv62",  [&](){return isHH ?  hh.LHEPdfWeight()[62]*hh.LHEPdfWeightNorm()[62] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv63",  [&](){return isHH ?  hh.LHEPdfWeight()[63]*hh.LHEPdfWeightNorm()[63] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv64",  [&](){return isHH ?  hh.LHEPdfWeight()[64]*hh.LHEPdfWeightNorm()[64] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv65",  [&](){return isHH ?  hh.LHEPdfWeight()[65]*hh.LHEPdfWeightNorm()[65] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv66",  [&](){return isHH ?  hh.LHEPdfWeight()[66]*hh.LHEPdfWeightNorm()[66] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv67",  [&](){return isHH ?  hh.LHEPdfWeight()[67]*hh.LHEPdfWeightNorm()[67] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv68",  [&](){return isHH ?  hh.LHEPdfWeight()[68]*hh.LHEPdfWeightNorm()[68] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv69",  [&](){return isHH ?  hh.LHEPdfWeight()[69]*hh.LHEPdfWeightNorm()[69] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv70",  [&](){return isHH ?  hh.LHEPdfWeight()[70]*hh.LHEPdfWeightNorm()[70] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv71",  [&](){return isHH ?  hh.LHEPdfWeight()[71]*hh.LHEPdfWeightNorm()[71] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv72",  [&](){return isHH ?  hh.LHEPdfWeight()[72]*hh.LHEPdfWeightNorm()[72] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv73",  [&](){return isHH ?  hh.LHEPdfWeight()[73]*hh.LHEPdfWeightNorm()[73] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv74",  [&](){return isHH ?  hh.LHEPdfWeight()[74]*hh.LHEPdfWeightNorm()[74] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv75",  [&](){return isHH ?  hh.LHEPdfWeight()[75]*hh.LHEPdfWeightNorm()[75] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv76",  [&](){return isHH ?  hh.LHEPdfWeight()[76]*hh.LHEPdfWeightNorm()[76] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv77",  [&](){return isHH ?  hh.LHEPdfWeight()[77]*hh.LHEPdfWeightNorm()[77] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv78",  [&](){return isHH ?  hh.LHEPdfWeight()[78]*hh.LHEPdfWeightNorm()[78] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv79",  [&](){return isHH ?  hh.LHEPdfWeight()[79]*hh.LHEPdfWeightNorm()[79] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv80",  [&](){return isHH ?  hh.LHEPdfWeight()[80]*hh.LHEPdfWeightNorm()[80] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv81",  [&](){return isHH ?  hh.LHEPdfWeight()[81]*hh.LHEPdfWeightNorm()[81] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv82",  [&](){return isHH ?  hh.LHEPdfWeight()[82]*hh.LHEPdfWeightNorm()[82] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv83",  [&](){return isHH ?  hh.LHEPdfWeight()[83]*hh.LHEPdfWeightNorm()[83] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv84",  [&](){return isHH ?  hh.LHEPdfWeight()[84]*hh.LHEPdfWeightNorm()[84] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv85",  [&](){return isHH ?  hh.LHEPdfWeight()[85]*hh.LHEPdfWeightNorm()[85] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv86",  [&](){return isHH ?  hh.LHEPdfWeight()[86]*hh.LHEPdfWeightNorm()[86] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv87",  [&](){return isHH ?  hh.LHEPdfWeight()[87]*hh.LHEPdfWeightNorm()[87] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv88",  [&](){return isHH ?  hh.LHEPdfWeight()[88]*hh.LHEPdfWeightNorm()[88] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv89",  [&](){return isHH ?  hh.LHEPdfWeight()[89]*hh.LHEPdfWeightNorm()[89] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv90",  [&](){return isHH ?  hh.LHEPdfWeight()[90]*hh.LHEPdfWeightNorm()[90] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv91",  [&](){return isHH ?  hh.LHEPdfWeight()[91]*hh.LHEPdfWeightNorm()[91] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv92",  [&](){return isHH ?  hh.LHEPdfWeight()[92]*hh.LHEPdfWeightNorm()[92] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv93",  [&](){return isHH ?  hh.LHEPdfWeight()[93]*hh.LHEPdfWeightNorm()[93] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv94",  [&](){return isHH ?  hh.LHEPdfWeight()[94]*hh.LHEPdfWeightNorm()[94] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv95",  [&](){return isHH ?  hh.LHEPdfWeight()[95]*hh.LHEPdfWeightNorm()[95] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv96",  [&](){return isHH ?  hh.LHEPdfWeight()[96]*hh.LHEPdfWeightNorm()[96] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv97",  [&](){return isHH ?  hh.LHEPdfWeight()[97]*hh.LHEPdfWeightNorm()[97] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv98",  [&](){return isHH ?  hh.LHEPdfWeight()[98]*hh.LHEPdfWeightNorm()[98] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv99",  [&](){return isHH ?  hh.LHEPdfWeight()[99]*hh.LHEPdfWeightNorm()[99] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv100",  [&](){return isHH ?  hh.LHEPdfWeight()[100]*hh.LHEPdfWeightNorm()[100] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv101",  [&](){return isHH ?  hh.LHEPdfWeight()[101]*hh.LHEPdfWeightNorm()[101] : 1.0;});
    cutflow.addWgtSyst("LHEPDFEigenv102",  [&](){return isHH ?  hh.LHEPdfWeight()[102]*hh.LHEPdfWeightNorm()[102] : 1.0;});
    */
        
    //QCD scale uncertainty for the HH signal acceptance: take the envelope for weights [0,1,3,4,5,8,9]
    /*
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Factorization_and_renormalizatio
        ['LHE scale variation weights (w_var / w_nominal)',
        ' [0] is renscfact=0.5d0 facscfact=0.5d0 ',
        ' [1] is renscfact=0.5d0 facscfact=1d0 ',
        ' [2] is renscfact=0.5d0 facscfact=2d0 ',
        ' [3] is renscfact=1d0 facscfact=0.5d0 ',
        ' [4] is renscfact=1d0 facscfact=1d0 ',
        ' [5] is renscfact=1d0 facscfact=2d0 ',
        ' [6] is renscfact=2d0 facscfact=0.5d0 ',
        ' [7] is renscfact=2d0 facscfact=1d0 ',
        ' [8] is renscfact=2d0 facscfact=2d0 ']
    */
        
    cutflow.addWgtSyst("QCDscale0",  [&](){return isHH ?  hh.LHEScaleWeight()[0] : 1.0;});
    cutflow.addWgtSyst("QCDscale1",  [&](){return isHH ?  hh.LHEScaleWeight()[1] : 1.0;});
    cutflow.addWgtSyst("QCDscale2",  [&](){return isHH ?  hh.LHEScaleWeight()[2] : 1.0;});
    cutflow.addWgtSyst("QCDscale3",  [&](){return isHH ?  hh.LHEScaleWeight()[3] : 1.0;});
    cutflow.addWgtSyst("QCDscale4",  [&](){return isHH ?  hh.LHEScaleWeight()[4] : 1.0;});
    cutflow.addWgtSyst("QCDscale5",  [&](){return isHH ?  hh.LHEScaleWeight()[5] : 1.0;});
    cutflow.addWgtSyst("QCDscale6",  [&](){return isHH ?  hh.LHEScaleWeight()[6] : 1.0;});
    cutflow.addWgtSyst("QCDscale7",  [&](){return isHH ?  hh.LHEScaleWeight()[7] : 1.0;});
    cutflow.addWgtSyst("QCDscale8",  [&](){return isHH ?  hh.LHEScaleWeight()[8] : 1.0;});
        
    /*   
    cutflow.addWgtSyst("QCDscale0",  [&](){return isHH ?  hh.LHEScaleWeight()[0]*hh.LHEScaleWeightNorm()[0] : 1.0;});
    cutflow.addWgtSyst("QCDscale1",  [&](){return isHH ?  hh.LHEScaleWeight()[1]*hh.LHEScaleWeightNorm()[1] : 1.0;});
    cutflow.addWgtSyst("QCDscale2",  [&](){return isHH ?  hh.LHEScaleWeight()[2]*hh.LHEScaleWeightNorm()[2] : 1.0;});
    cutflow.addWgtSyst("QCDscale3",  [&](){return isHH ?  hh.LHEScaleWeight()[3]*hh.LHEScaleWeightNorm()[3] : 1.0;});
    cutflow.addWgtSyst("QCDscale4",  [&](){return isHH ?  hh.LHEScaleWeight()[4]*hh.LHEScaleWeightNorm()[4] : 1.0;});
    cutflow.addWgtSyst("QCDscale5",  [&](){return isHH ?  hh.LHEScaleWeight()[5]*hh.LHEScaleWeightNorm()[5] : 1.0;});
    cutflow.addWgtSyst("QCDscale6",  [&](){return isHH ?  hh.LHEScaleWeight()[6]*hh.LHEScaleWeightNorm()[6] : 1.0;});
    cutflow.addWgtSyst("QCDscale7",  [&](){return isHH ?  hh.LHEScaleWeight()[7]*hh.LHEScaleWeightNorm()[7] : 1.0;});
    cutflow.addWgtSyst("QCDscale8",  [&](){return isHH ?  hh.LHEScaleWeight()[8]*hh.LHEScaleWeightNorm()[8] : 1.0;});
    */  
    }
 
    if(isHH){    
        //trigger eff correction uncertainty
        cutflow.addWgtSyst("trigCorrHH"+year_+"Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1MassSD(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2MassSD(), hh.fatJet2PNetXbb(), 0, 0,  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"up");});
        cutflow.addWgtSyst("trigCorrHH"+year_+"Down", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1MassSD(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2MassSD(), hh.fatJet2PNetXbb(), 0, 0, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"down");});   
        //parton shower uncertainty for the gluon fusion and VBF HH signal yield (ISR and FSR)
        cutflow.addWgtSyst("ISRPartonShowerUp",  [&](){ return (isHH && fabs(hh.PSWeight()[0])< 20) ? hh.PSWeight()[0] : 1.0;}); //ISR=2 FSR=1
        cutflow.addWgtSyst("ISRPartonShowerDown",  [&](){ return (isHH && fabs(hh.PSWeight()[2])< 20) ? hh.PSWeight()[2] : 1.0;}); //ISR=0.5 FSR=1
        cutflow.addWgtSyst("FSRPartonShowerUp",  [&](){ return (isHH && fabs(hh.PSWeight()[1])< 20) ? hh.PSWeight()[1] : 1.0;}); //ISR=1 FSR=2
        cutflow.addWgtSyst("FSRPartonShowerDown",  [&](){ return (isHH && fabs(hh.PSWeight()[3])< 20) ? hh.PSWeight()[3] : 1.0;}); //ISR=1 FSR=0.5
    }
    
    //BDT modeling uncertainty for ttbar
    cutflow.addWgtSyst("BDTv8p2ShapeUp",  [&](){return isTTJets ? ( hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2()  < 0.00008 ? 0.96 : ( hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2()  <  0.0002 ?  0.91 : (hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2()  < 0.0004 ? 0.90 :  1.12))) : 1.0;});
    cutflow.addWgtSyst("BDTv8p2ShapeDown",  [&](){return isTTJets ? ( hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2()  < 0.00008 ? 1.04 : ( hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2()  <  0.0002 ?  1.09 : (hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2()  < 0.0004 ? 1.10 :  0.88))) : 1.0;});
    
    //ttbar recoil correction uncertainty for ttbar
    cutflow.addWgtSyst("ttJetsCorrUp",  [&](){return isTTJets ?  ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 1)/ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 0) : 1.0;});
    cutflow.addWgtSyst("ttJetsCorrDown",  [&](){return isTTJets ? ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), -1)/ttjets_sf.getScaleFactorsFit(year_, hh.hh_pt(), 0) : 1.0;});
    
    //PNet shape uncertainty for ttbar 
    cutflow.addWgtSyst("PNetShapeUp",  [&](){return isTTJets ?  ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 1)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 1)/(ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)) : 1.0;});
    cutflow.addWgtSyst("PNetShapeDown",  [&](){return isTTJets ?  ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), -1)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), -1)/(ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet1PNetXbb(), 0)*ttjets_sf.getPNetXbbShapeScaleFactors(year_, hh.fatJet2PNetXbb(), 0)) : 1.0;});
    
    //pileup reweighting uncertainty
    cutflow.addWgtSyst("pileupWeightUp",  [&](){return hh.puWeight()!=0 ? hh.puWeightUp()/hh.puWeight() : 1.0;});
    cutflow.addWgtSyst("pileupWeightDown",  [&](){return hh.puWeight()!=0 ? hh.puWeightDown()/hh.puWeight() : 1.0;});
    
    //PNetHbbScaleFactors uncertainty   
    if(doPNetSFsys){
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin11Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 1, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 1, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 1, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 1, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin11Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 1, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 1, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 1, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 1, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin12Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 2, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 2, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 2, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 2, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin12Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 2, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 2, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 2, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 2, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin13Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 3, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 3, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 3, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 3, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin13Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 3, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 3, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 3, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 3, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin14Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 4, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 4, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 4, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 4, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin14Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 4, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 4, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 4, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 4, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin15Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 5, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 5, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 5, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 5, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin15Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 5, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 5, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 1, 5, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 1, 5, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin21Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 1, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 1, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 1, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 1, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin21Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 1, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 1, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 1, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 1, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin22Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 2, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 2, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 2, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 2, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin22Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 2, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 2, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 2, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 2, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin23Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 3, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 3, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 3, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 3, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin23Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 3, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 3, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 3, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 3, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin24Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 4, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 4, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 4, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 4, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin24Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 4, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 4, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 4, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 4, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin25Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 5, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 5, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 5, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 5, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin25Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 5, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 5, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 2, 5, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 2, 5, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin31Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 1, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 1, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 1, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 1, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin31Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 1, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 1, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 1, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 1, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin32Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 2, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 2, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 2, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 2, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin32Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 2, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 2, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 2, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 2, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin33Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 3, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 3, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 3, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 3, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin33Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 3, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 3, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 3, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 3, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin34Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 4, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 4, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 4, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 4, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin34Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 4, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 4, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 4, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 4, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin35Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 5, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 5, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 5, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 5, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin35Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 5, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 5, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 3, 5, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 3, 5, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin41Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 1, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 1, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 1, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 1, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin41Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 1, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 1, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 1, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 1, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin42Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 2, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 2, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 2, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 2, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin42Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 2, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 2, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 2, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 2, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin43Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 3, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 3, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 3, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 3, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin43Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 3, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 3, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 3, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 3, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin44Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 4, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 4, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 4, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 4, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin44Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 4, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 4, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 4, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 4, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin45Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 5, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 5, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 5, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 5, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin45Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 5, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 5, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 4, 5, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 4, 5, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin51Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 1, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 1, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 1, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 1, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin51Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 1, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 1, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 1, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 1, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin52Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 2, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 2, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 2, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 2, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin52Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 2, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 2, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 2, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 2, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin53Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 3, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 3, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 3, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 3, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin53Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 3, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 3, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 3, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 3, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin54Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 4, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 4, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 4, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 4, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin54Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 4, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 4, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 4, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 4, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin55Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 5, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 5, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 5, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 5, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin55Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 5, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 5, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 5, 5, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 5, 5, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin61Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 1, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 1, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 1, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 1, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin61Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 1, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 1, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 1, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 1, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin62Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 2, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 2, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 2, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 2, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin62Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 2, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 2, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 2, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 2, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin63Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 3, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 3, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 3, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 3, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin63Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 3, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 3, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 3, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 3, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin64Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 4, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 4, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 4, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 4, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin64Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 4, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 4, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 4, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 4, -1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin65Up",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 5, 1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 5, 1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 5, 1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 5, 1)/nominal_weight;
        }
        return varied_weight;
    });
    cutflow.addWgtSyst("PNetHbbScaleFactors"+year_+"bin65Down",  [&](){
        float nominal_weight =  1.0;
        float varied_weight = 1.0;
        if(isH){
            float fatJet1_genH_dR = DeltaR(hh.fatJet1Eta(), hh.fatJet1Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            float fatJet2_genH_dR = DeltaR(hh.fatJet2Eta(), hh.fatJet2Phi(), hh.genHiggs1Eta(), hh.genHiggs1Phi());
            if( (fatJet1_genH_dR < fatJet2_genH_dR) && (fatJet1_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0);
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 5, -1)/nominal_weight;
            }  
            else if( (fatJet1_genH_dR > fatJet2_genH_dR) && (fatJet2_genH_dR < 0.4) ){
                nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0); 
                if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 5, -1)/nominal_weight;
            } 
        }
        else if(isHH){
            nominal_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 0, 0, 0)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 0, 0, 0);
            if(nominal_weight!=0) varied_weight = PNet_sf.getPNetHbbScaleFactors(hh.fatJet1Pt(), hh.fatJet1PNetXbb(), 6, 5, -1)*PNet_sf.getPNetHbbScaleFactors(hh.fatJet2Pt(), hh.fatJet2PNetXbb(), 6, 5, -1)/nominal_weight;
        }
        return varied_weight;
    });
        
    }
    
    if(dotrigsys){

    //trigger eff SF uncertainty
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin3Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 3,  hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin3Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 3, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin4Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 4, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin4Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 4, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin5Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 5, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin5Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 5, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin6Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 6, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin6Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 6, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin7Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 7, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin7Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 7, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin8Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 8, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin8Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 8, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin9Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 9, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin9Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 9, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin10Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 10, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin10Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 10, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin11Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 11, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin11Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 11, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin12Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 12, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin12Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 12, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin15Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 15, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin15Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 15, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin16Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 16, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin16Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 16, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin17Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 17, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin17Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 17, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin18Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 18, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin18Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 18, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin19Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 19, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin19Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 19, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin20Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 20, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin20Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 20, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin21Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 21, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin21Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 21, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin22Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 22, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin22Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 22, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin23Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 23, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin23Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 23, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin24Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 24, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin24Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 24, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin27Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 27, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin27Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 27, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin28Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 28, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin28Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 28, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin29Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 29, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin29Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 29, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin30Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 30, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin30Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 30, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin31Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 31, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin31Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 31, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin32Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 32, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin32Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 32, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin33Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 33, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin33Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 33, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin34Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 34, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin34Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 34, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin35Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 35, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin35Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 35, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin36Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 36, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin36Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 36, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin39Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 39, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin39Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 39, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin40Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 40, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin40Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 40, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin41Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 41, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin41Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 41, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin42Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 42, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin42Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 42, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin43Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 43, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin43Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 43, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin44Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 44, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin44Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 44, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin45Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 45, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin45Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 45, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin46Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 46, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin46Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 46, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin47Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 47, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin47Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 47, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin48Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 48, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin48Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 48, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin51Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 51, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin51Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 51, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin52Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 52, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin52Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 52, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin53Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 53, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin53Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 53, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin54Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 54, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin54Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 54, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin55Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 55, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin55Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 55, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin56Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 56, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin56Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 56, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin57Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 57, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin57Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 57, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin58Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 58, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin58Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 58, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin59Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 59, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin59Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 59, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin60Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 60, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin60Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 60, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin63Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 63, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin63Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 63, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin64Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 64, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin64Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 64, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin65Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 65, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin65Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 65, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin66Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 66, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin66Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 66, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin67Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 67, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin67Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 67, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin68Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 68, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin68Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 68, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin69Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 69, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin69Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 69, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin70Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 70, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin70Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 70, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin71Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 71, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin71Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 71, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin72Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 72, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin72Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 72, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin75Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 75, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin75Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 75, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin76Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 76, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin76Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 76, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin77Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 77, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin77Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 77, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin78Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 78, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin78Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 78, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin79Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 79, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin79Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 79, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin80Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 80, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin80Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 80, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin81Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 81, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin81Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 81, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin82Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 82, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin82Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 82, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin83Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 83, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin83Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 83, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin84Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 84, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin84Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 84, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin87Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 87, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin87Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 87, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin88Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 88, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin88Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 88, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin89Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 89, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin89Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 89, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin90Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 90, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin90Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 90, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin91Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 91, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin91Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 91, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin92Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 92, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin92Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 92, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin93Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 93, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin93Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 93, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin94Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 94, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin94Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 94, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin95Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 95, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin95Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 95, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin96Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 96, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin96Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 96, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin99Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 99, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin99Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 99, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin100Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 100, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin100Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 100, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin101Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 101, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin101Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 101, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin102Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 102, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin102Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 102, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin103Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 103, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin103Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 103, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin104Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 104, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin104Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 104, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin105Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 105, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin105Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 105, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin106Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 106, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin106Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 106, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin107Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 107, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin107Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 107, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin108Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 108, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin108Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 108, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin111Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 111, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin111Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 111, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin112Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 112, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin112Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 112, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin113Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 113, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin113Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 113, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin114Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 114, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin114Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 114, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin115Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 115, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin115Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 115, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin116Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 116, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin116Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 116, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin117Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 117, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin117Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 117, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin118Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 118, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin118Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 118, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin119Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 119, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin119Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 119, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin120Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 120, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin120Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 120, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin123Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 123, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin123Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 123, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin124Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 124, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin124Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 124, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin125Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 125, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin125Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 125, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin126Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 126, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin126Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 126, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin127Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 127, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin127Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 127, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin128Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 128, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin128Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 128, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin129Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 129, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin129Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 129, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin130Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 130, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin130Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 130, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin131Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 131, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin131Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 131, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin132Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 132, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin132Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 132, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin135Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 135, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin135Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 135, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin136Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 136, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin136Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 136, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin137Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 137, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin137Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 137, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin138Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 138, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin138Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 138, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin139Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 139, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin139Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 139, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin140Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 140, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin140Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 140, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin141Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 141, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin141Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 141, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin142Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 142, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin142Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 142, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin143Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 143, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin143Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 143, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin144Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 144, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin144Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 144, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin147Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 147, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin147Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 147, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin148Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 148, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin148Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 148, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin149Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 149, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin149Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 149, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin150Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 150, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin150Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 150, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin151Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 151, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin151Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 151, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin152Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 152, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin152Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 152, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin153Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 153, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin153Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 153, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin154Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 154, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin154Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 154, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin155Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 155, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin155Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 155, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin156Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 156, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin156Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 156, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin159Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 159, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin159Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 159, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin160Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 160, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin160Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 160, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin161Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 161, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin161Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 161, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin162Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 162, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin162Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 162, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin163Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 163, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin163Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 163, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin164Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 164, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin164Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 164, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin165Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 165, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin165Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 165, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin166Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 166, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin166Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 166, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin167Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 167, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin167Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 167, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin168Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 168, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin168Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 168, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin171Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 171, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin171Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 171, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin172Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 172, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin172Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 172, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin173Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 173, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin173Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 173, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin174Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 174, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin174Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 174, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin175Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 175, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin175Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 175, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin176Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 176, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin176Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 176, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin177Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 177, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin177Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 177, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin178Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 178, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin178Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 178, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin179Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 179, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin179Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 179, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin180Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 180, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin180Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 180, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin183Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 183, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin183Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 183, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin184Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 184, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin184Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 184, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin185Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 185, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin185Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 185, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin186Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 186, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin186Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 186, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin187Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 187, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin187Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 187, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin188Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 188, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin188Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 188, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin189Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 189, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin189Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 189, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin190Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 190, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin190Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 190, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin191Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 191, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin191Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 191, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin192Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 192, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin192Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 192, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin195Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 195, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin195Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 195, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin196Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 196, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin196Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 196, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin197Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 197, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin197Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 197, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin198Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 198, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin198Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 198, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin199Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 199, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin199Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 199, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin200Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 200, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin200Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 200, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin201Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 201, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin201Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 201, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin202Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 202, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin202Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 202, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin203Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 203, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin203Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 203, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin204Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 204, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin204Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 204, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin207Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 207, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin207Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 207, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin208Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 208, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin208Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 208, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin209Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 209, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin209Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 209, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin210Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 210, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin210Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 210, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin211Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 211, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin211Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 211, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin212Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 212, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin212Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 212, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin213Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 213, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin213Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 213, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin214Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 214, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin214Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 214, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin215Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 215, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin215Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 215, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin216Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 216, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin216Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 216, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin219Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 219, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin219Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 219, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin220Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 220, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin220Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 220, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin221Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 221, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin221Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 221, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin222Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 222, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin222Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 222, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin223Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 223, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin223Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 223, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin224Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 224, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin224Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 224, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin225Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 225, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin225Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 225, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin226Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 226, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin226Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 226, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin227Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 227, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin227Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 227, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin228Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 228, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin228Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 228, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin231Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 231, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin231Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 231, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin232Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 232, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin232Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 232, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin233Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 233, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin233Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 233, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin234Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 234, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin234Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 234, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin235Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 235, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin235Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 235, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin236Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 236, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin236Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 236, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin237Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 237, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin237Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 237, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin238Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 238, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin238Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 238, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin239Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 239, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin239Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 239, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin240Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 240, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin240Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 240, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin243Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 243, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin243Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 243, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin244Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 244, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin244Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 244, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin245Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 245, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin245Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 245, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin246Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 246, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin246Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 246, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin247Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 247, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin247Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 247, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin248Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 248, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin248Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 248, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin249Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 249, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin249Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 249, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin250Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 250, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin250Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 250, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin251Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 251, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin251Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 251, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin252Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 252, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin252Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 252, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin255Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 255, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin255Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 255, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin256Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 256, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin256Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 256, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin257Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 257, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin257Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 257, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin258Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 258, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin258Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 258, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin259Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 259, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin259Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 259, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin260Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 260, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin260Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 260, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin261Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 261, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin261Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 261, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin262Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 262, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin262Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 262, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin263Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 263, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin263Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 263, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin264Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 264, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin264Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 264, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin267Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 267, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin267Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 267, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin268Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 268, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin268Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 268, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin269Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 269, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin269Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 269, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin270Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 270, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin270Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 270, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin271Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 271, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin271Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 271, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin272Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 272, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin272Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 272, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin273Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 273, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin273Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 273, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin274Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 274, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin274Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 274, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin275Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 275, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin275Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 275, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin276Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 276, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin276Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 276, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin279Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 279, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin279Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 279, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin280Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 280, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin280Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 280, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin281Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 281, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin281Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 281, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin282Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 282, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin282Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 282, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin283Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 283, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin283Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 283, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin284Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 284, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin284Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 284, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin285Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 285, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin285Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 285, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin286Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 286, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin286Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 286, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin287Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 287, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin287Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 287, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin288Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 288, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin288Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 288, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin291Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 291, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin291Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 291, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin292Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 292, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin292Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 292, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin293Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 293, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin293Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 293, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin294Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 294, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin294Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 294, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin295Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 295, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin295Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 295, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin296Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 296, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin296Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 296, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin297Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 297, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin297Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 297, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin298Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 298, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin298Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 298, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin299Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 299, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin299Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 299, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin300Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 300, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin300Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 300, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin303Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 303, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin303Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 303, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin304Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 304, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin304Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 304, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin305Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 305, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin305Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 305, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin306Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 306, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin306Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 306, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin307Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 307, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin307Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 307, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin308Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 308, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin308Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 308, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin309Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 309, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin309Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 309, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin310Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 310, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin310Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 310, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin311Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 311, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin311Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 311, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin312Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 312, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin312Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 312, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin315Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 315, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin315Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 315, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin316Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 316, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin316Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 316, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin317Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 317, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin317Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 317, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin318Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 318, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin318Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 318, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin319Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 319, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin319Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 319, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin320Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 320, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin320Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 320, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin321Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 321, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin321Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 321, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin322Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 322, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin322Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 322, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin323Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 323, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin323Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 323, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin324Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 324, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin324Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 324, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin327Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 327, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin327Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 327, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin328Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 328, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin328Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 328, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin329Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 329, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin329Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 329, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin330Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 330, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin330Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 330, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin331Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 331, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin331Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 331, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin332Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 332, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin332Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 332, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin333Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 333, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin333Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 333, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin334Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 334, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin334Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 334, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin335Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 335, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin335Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 335, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin336Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 336, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin336Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 336, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin339Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 339, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin339Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 339, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin340Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 340, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin340Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 340, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin341Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 341, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin341Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 341, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin342Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 342, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin342Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 342, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin343Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 343, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin343Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 343, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin344Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 344, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin344Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 344, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin345Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 345, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin345Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 345, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin346Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 346, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin346Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 346, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin347Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 347, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin347Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 347, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin348Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 348, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin348Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 348, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin351Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 351, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin351Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 351, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin352Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 352, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin352Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 352, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin353Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 353, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin353Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 353, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin354Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 354, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin354Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 354, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin355Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 355, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin355Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 355, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin356Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 356, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin356Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 356, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin357Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 357, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin357Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 357, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin358Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 358, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin358Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 358, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin359Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 359, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin359Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 359, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin360Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 360, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin360Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 360, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin363Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 363, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin363Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 363, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin364Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 364, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin364Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 364, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin365Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 365, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin365Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 365, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin366Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 366, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin366Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 366, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin367Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 367, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin367Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 367, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin368Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 368, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin368Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 368, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin369Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 369, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin369Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 369, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin370Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 370, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin370Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 370, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin371Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 371, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin371Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 371, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin372Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 372, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin372Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 372, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin375Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 375, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin375Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 375, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin376Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 376, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin376Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 376, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin377Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 377, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin377Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 377, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin378Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 378, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin378Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 378, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin379Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 379, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin379Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 379, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin380Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 380, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin380Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 380, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin381Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 381, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin381Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 381, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin382Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 382, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin382Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 382, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin383Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 383, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin383Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 383, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin384Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 384, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin384Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 384, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin387Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 387, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin387Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 387, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin388Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 388, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin388Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 388, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin389Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 389, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin389Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 389, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin390Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 390, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin390Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 390, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin391Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 391, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin391Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 391, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin392Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 392, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin392Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 392, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin393Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 393, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin393Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 393, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin394Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 394, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin394Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 394, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin395Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 395, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin395Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 395, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin396Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 396, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin396Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 396, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin399Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 399, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin399Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 399, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin400Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 400, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin400Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 400, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin401Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 401, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin401Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 401, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin402Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 402, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin402Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 402, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin403Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 403, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin403Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 403, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin404Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 404, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin404Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 404, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin405Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 405, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin405Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 405, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin406Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 406, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin406Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 406, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin407Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 407, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin407Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 407, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin408Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 408, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin408Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 408, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin411Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 411, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin411Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 411, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin412Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 412, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin412Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 412, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin413Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 413, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin413Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 413, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin414Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 414, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin414Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 414, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin415Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 415, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin415Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 415, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin416Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 416, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin416Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 416, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin417Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 417, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin417Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 417, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin418Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 418, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin418Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 418, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin419Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 419, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin419Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 419, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin420Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 420, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin420Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 420, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin423Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 423, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin423Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 423, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin424Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 424, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin424Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 424, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin425Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 425, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin425Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 425, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin426Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 426, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin426Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 426, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin427Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 427, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin427Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 427, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin428Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 428, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin428Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 428, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin429Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 429, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin429Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 429, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin430Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 430, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin430Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 430, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin431Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 431, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin431Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 431, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin432Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 432, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin432Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 432, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin435Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 435, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin435Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 435, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin436Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 436, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin436Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 436, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin437Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 437, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin437Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 437, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin438Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 438, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin438Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 438, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin439Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 439, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin439Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 439, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin440Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 440, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin440Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 440, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin441Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 441, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin441Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 441, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin442Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 442, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin442Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 442, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin443Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 443, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin443Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 443, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin444Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 444, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin444Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 444, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin447Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 447, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin447Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 447, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin448Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 448, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin448Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 448, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin449Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 449, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin449Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 449, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin450Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 450, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin450Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 450, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin451Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 451, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin451Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 451, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin452Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 452, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin452Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 452, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin453Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 453, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin453Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 453, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin454Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 454, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin454Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 454, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin455Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 455, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin455Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 455, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin456Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 456, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin456Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 456, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin459Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 459, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin459Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 459, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin460Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 460, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin460Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 460, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin461Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 461, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin461Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 461, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin462Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 462, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin462Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 462, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin463Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 463, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin463Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 463, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin464Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 464, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin464Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 464, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin465Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 465, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin465Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 465, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin466Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 466, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin466Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 466, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin467Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 467, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin467Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 467, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin468Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 468, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin468Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 468, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin471Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 471, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin471Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 471, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin472Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 472, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin472Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 472, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin473Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 473, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin473Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 473, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin474Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 474, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin474Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 474, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin475Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 475, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin475Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 475, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin476Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 476, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin476Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 476, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin477Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 477, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin477Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 477, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin478Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 478, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin478Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 478, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin479Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 479, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin479Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 479, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin480Up", [&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), 1, 480, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    cutflow.addWgtSyst("triggerEffSF"+year_+"bin480Down",[&](){return trig_sf.get_unc_ratio(hh.fatJet1Pt(), hh.fatJet1Mass(), hh.fatJet1PNetXbb(), hh.fatJet2Pt(), hh.fatJet2Mass(), hh.fatJet2PNetXbb(), -1, 480, hh.disc_qcd_and_ttbar_Run2_enhanced_v8p2(), year_, isHH,"nom");});
    }
}

//book histograms for cuts
if(not doSystematics) cutflow.bookHistogramsForCutAndBelow(histograms, "CutWeight");
else
{
    if(input.find("Tau32TopSkim") != std::string::npos) // this is ttbar enriched input
    {
        cutflow.bookHistogramsForCut(histograms, "TTBarCR");
        //cutflow.bookHistogramsForCut(histograms, "TTBarCRTight");
        //cutflow.bookHistogramsForCut(histograms, "TTBarCRBDT1v8p2");
        //cutflow.bookHistogramsForCut(histograms, "TTBarCRBDT2v8p2");
        //cutflow.bookHistogramsForCut(histograms, "TTBarCRBDT3v8p2");
        //cutflow.bookHistogramsForCut(histograms, "TTBarCRBDT4v8p2");
    }
    else // normal  signal enriched input
    {
        //cutflow.bookHistogramsForCut(histograms, "SRv24Bin1");
        //cutflow.bookHistogramsForCut(histograms, "SRv24Bin2");
        //cutflow.bookHistogramsForCut(histograms, "SRv24Bin3");
        //cutflow.bookHistogramsForCut(histograms, "SRv24Bin4");
        //cutflow.bookHistogramsForCut(histograms, "FailSRv24");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandv24Bin1");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandv24Bin2");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandv24Bin3");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandv24Bin4");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandFailSRv24");
        //cutflow.bookHistogramsForCut(histograms, "FitCRv24");
        //cutflow.bookHistogramsForCut(histograms, "FailFitCRv24");

        cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1");

        //cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1PNetp95");
        cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1PNetp9");
        //cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1PNetp8");
        //cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1PNetp5");
        //cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1PNetp2");

        //cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin1PNetp0");
        cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin2");
        cutflow.bookHistogramsForCut(histograms, "SRv8p2Bin3");
        cutflow.bookHistogramsForCut(histograms, "FailSRv8p2");

        //check preselection 
        //cutflow.bookHistogramsForCut(histograms, "CutfatJetsMassSD");

        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandv8p2Bin1");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandv8p2Bin2");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandv8p2Bin3");
        //cutflow.bookHistogramsForCut(histograms, "J2MassSideBandFailSRv8p2");
        //cutflow.bookHistogramsForCut(histograms, "FitCRv8p2");
        //cutflow.bookHistogramsForCut(histograms, "FailFitCRv8p2");

    }
}
cutflow.bookCutflows();

int iEntry = 0;

if(saveSkim) outfile_skim->cd();
TTree *tree_out;

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
tree_out->Branch("weight", &weight, "weight/F");
tree_out->Branch("fatJet1MassSD", &fatJet1MassSD, "fatJet1MassSD/F");
 tree_out->Branch("fatJet1MassRegressed", &fatJet1MassRegressed, "fatJet1MassRegressed/F");
tree_out->Branch("fatJet1PNetXbb", &fatJet1PNetXbb, "fatJet1PNetXbb/F");
tree_out->Branch("fatJet1Pt", &fatJet1Pt, "fatJet1Pt/F");
tree_out->Branch("fatJet1Eta", &fatJet1Eta, "fatJet1Eta/F");
tree_out->Branch("fatJet1Phi", &fatJet1Phi, "fatJet1Phi/F");
tree_out->Branch("fatJet1PtOverMHH", &fatJet1PtOverMHH, "fatJet1PtOverMHH/F");
tree_out->Branch("fatJet1PtOverMSD", &fatJet1PtOverMSD, "fatJet1PtOverMSD/F");
 tree_out->Branch("fatJet1PtOverMRegressed", &fatJet1PtOverMRegressed, "fatJet1PtOverMRegressed/F");
tree_out->Branch("fatJet2MassSD", &fatJet2MassSD, "fatJet2MassSD/F");
tree_out->Branch("fatJet2MassRegressed", &fatJet2MassRegressed, "fatJet2MassRegressed/F");
tree_out->Branch("fatJet2PNetXbb", &fatJet2PNetXbb, "fatJet2PNetXbb/F");
tree_out->Branch("fatJet2Pt", &fatJet2Pt, "fatJet2Pt/F");
tree_out->Branch("fatJet2Eta", &fatJet2Eta, "fatJet2Eta/F");
tree_out->Branch("fatJet2Phi", &fatJet2Phi, "fatJet2Phi/F");
tree_out->Branch("fatJet2PtOverMHH", &fatJet2PtOverMHH, "fatJet2PtOverMHH/F");
tree_out->Branch("fatJet2PtOverMSD", &fatJet2PtOverMSD, "fatJet2PtOverMSD/F");
 tree_out->Branch("fatJet2PtOverMRegressed", &fatJet2PtOverMRegressed, "fatJet2PtOverMRegressed/F");
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
  tree_this = (TTree*)file_this->Get("Events");
  hh.Init(tree_this);
  int nEntries_this = tree_this->GetEntries();
  for(int iEntry_this=0; iEntry_this<nEntries_this; iEntry_this++)
  {
	hh.GetEntry(iEntry_this);
        if(saveSkim) outfile->cd();
	cutflow.fill();
	//if(saveSkim && cutflow.getCut("TwofatJets").pass)
	if(saveSkim && cutflow.getCut("BDTTrainPreSelection").pass)
	{
	  outfile_skim->cd();	
	  weight = isData ?  1.0 : lumi*hh.weight()*hh.triggerEff3DWeight()*hh.puWeight()* hh.l1PreFiringWeight() * hh.xsecWeight()*hh.genWeight();
	  fatJet1MassSD = hh.fatJet1MassSD();
	  fatJet1MassRegressed = hh.fatJet1MassRegressed();
	  fatJet1PNetXbb = hh.fatJet1PNetXbb();
	  fatJet1Pt = hh.fatJet1Pt();
	  fatJet1Eta = hh.fatJet1Eta();
	  fatJet1Phi = hh.fatJet1Phi();
	  fatJet1PtOverMHH = hh.fatJet1PtOverMHH();
	  fatJet1PtOverMSD = hh.fatJet1PtOverMSD();
	  fatJet1PtOverMRegressed = hh.fatJet1PtOverMRegressed();
	  fatJet2MassSD = hh.fatJet2MassSD();
	  fatJet2MassRegressed = hh.fatJet2MassRegressed();
	  fatJet2PNetXbb = hh.fatJet2PNetXbb();
	  fatJet2Pt = hh.fatJet2Pt();
	  fatJet2Eta = hh.fatJet2Eta();
	  fatJet2Phi = hh.fatJet2Phi();
	  fatJet2PtOverMHH = hh.fatJet2PtOverMHH();
	  fatJet2PtOverMSD = hh.fatJet2PtOverMSD();
	  fatJet2PtOverMRegressed = hh.fatJet2PtOverMRegressed();

	  abs_dEta_j1j2 = fabs(hh.fatJet1Eta() - hh.fatJet2Eta());
	  abs_dPhi_j1j2 = fabs(hh.fatJet1Phi() - hh.fatJet2Phi());
	  abs_dR_j1j2 = sqrt((hh.fatJet1Eta() - hh.fatJet2Eta())*(hh.fatJet1Eta() - hh.fatJet2Eta())  + (hh.fatJet1Phi() - hh.fatJet2Phi())*(hh.fatJet1Phi() - hh.fatJet2Phi()));
	  ptj2_over_ptj1 = hh.fatJet2Pt() / hh.fatJet1Pt();
	  mj2_over_mj1 = fatJet2MassSD/fatJet1MassSD;
	  mregj2_over_mregj1 = fatJet2MassRegressed/fatJet1MassRegressed;
	  hh_pt = hh.hh_pt();
	  hh_eta = hh.hh_eta();
	  hh_phi = hh.hh_phi();
	  hh_mass = hh.hh_mass();
	  TLorentzVector gh1,gh2;
	  gh1.SetPtEtaPhiM(hh.genHiggs1Pt(), hh.genHiggs1Eta(),hh.genHiggs1Phi(),125.0);
	  gh2.SetPtEtaPhiM(hh.genHiggs2Pt(), hh.genHiggs2Eta(),hh.genHiggs2Phi(),125.0);
	  gen_hh_mass = (gh1+gh2).M();
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
