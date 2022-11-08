#ifndef reweight_HH_h
#define reweight_HH_h
#include <TBranch.h> 
#include <TTree.h> 
#include <TH1F.h> 
#include <TFile.h> 
#include <TBits.h> 
#include <vector> 
#include <unistd.h> 
#include <sstream>
using namespace std; 

vector<double> get_CMS_EFT_benchmarks( string name, string year, bool cms_fake = false ){
vector<double> answer;
 if( name == "SM" or name == "sm" ) answer = {1, 1, 0, 0, 0};
 if( name ==  "1" ) answer = {7.5, 1, -1, 0, 0};
 if( name ==  "2" ) answer = {1.0, 1.0, 0.5, -0.8, 0.6};
 if( name ==  "3" ) answer = {1.0, 1.0, -1.5, 0.0, -0.8};
 if( name ==  "4" ) answer = {-3.5, 1.5, -3.0, 0.0, 0.0};
 if( name ==  "5" ) answer = {1.0, 1.0, 0.0, 0.8, -1.0};
 if( name ==  "6" ) answer = {2.4, 1.0, 0.0, 0.2, -0.2};
 if( name ==  "7" ) answer = {5.0, 1.0, 0.0, 0.2, -0.2};
 if( name ==  "8" ) answer = {15.0, 1.0, 0.0, -1.0, 1.0};
 if( name ==  "9" ) answer = {1.0, 1.0, 1.0, -0.6, 0.6};
 if( name ==  "10" ) answer = {10.0, 1.5, -1.0, 0.0, 0.0};
 if( name ==  "11" ) answer = {2.4, 1.0, 0.0, 1.0, -1.0};
 if( name ==  "12" ) answer = {15.0, 1.0, 1.0, 0.0, 0.0};

 if( (year == "2017" or year == "2018") and answer.size() and cms_fake ){
   answer[4] = 1.0;
   double tr = answer[0];
   answer[0] = answer[1];
   answer[1] = tr;
 }
 if( (year == "2016") and answer.size() and cms_fake ){
   double tr = answer[0];
   answer[0] = answer[1];
   answer[1] = tr;
 }

 if( name == "BOX" or name == "box" ) answer = {0, 1, 0, 0, 0};
 if( name ==  "cHHH0" ) answer = {0.0,  1.0, 0.0, 0.0, 0.0};
 if( name ==  "cHHH1" ) answer = {1.0,  1.0, 0.0, 0.0, 0.0};
 if( name ==  "cHHH2" ) answer = {2.45, 1.0, 0.0, 0.0, 0.0};
 if( name ==  "cHHH5" ) answer = {5.0,  1.0, 0.0, 0.0, 0.0};
  
 // 8a from https://link.springer.com/article/10.1007/JHEP09(2018)057
 if( name ==  "8a" ) answer = {1.0, 1.0, 0.5, 0.8/3, 0.0};

 // New benchmarks from https://arxiv.org/pdf/1908.08923.pdf
 if( name == "JHEP03(2020)091_1" or name == "1b") answer = {  3.94, 0.94, -1./3.,  0.5 * 1.5,  1./3. * (-3.) };
 if( name == "JHEP03(2020)091_2" or name == "2b") answer = {  6.84, 0.61,  1./3.,  0.0 * 1.5, -1./3. * (-3.) };
 if( name == "JHEP03(2020)091_3" or name == "3b") answer = {  2.21, 1.05, -1./3.,  0.5 * 1.5,   0.5  * (-3.) };
 if( name == "JHEP03(2020)091_4" or name == "4b") answer = {  2.79, 0.61,  1./3., -0.5 * 1.5,  1./6. * (-3.) };
 if( name == "JHEP03(2020)091_5" or name == "5b") answer = {  3.95, 1.17, -1./3., 1./6.* 1.5,  -0.5  * (-3.) };
 if( name == "JHEP03(2020)091_6" or name == "6b") answer = {  5.68, 0.83,  1./3., -0.5 * 1.5,  1./3. * (-3.) };
 if( name == "JHEP03(2020)091_7" or name == "7b") answer = { -0.10, 0.94,     1., 1./6.* 1.5, -1./6. * (-3.) };

 return answer;
}

double get_eft_xsec_13TeV(vector<double> kappas, string order, string uncertantie=""){
  double answer_xsection = 0;
  double kl  = kappas[0];
  double kt  = kappas[1];
  double c2  = kappas[2];
  double cg  = kappas[3];
  double c2g = kappas[4];
  vector<double> couplings = { 
    // LO
    pow(kt,4), pow(c2,2), pow(kt,2)*pow(kl,2), pow(cg,2)*pow(kl,2),
    pow(c2g, 2), c2*pow(kt, 2), kl*pow(kt, 3), kt*kl*c2, cg*kl*c2, c2*c2g,
    cg*kl*pow(kt, 2), c2g*pow(kt, 2), 
    pow(kl,2)*cg*kt, c2g*kt*kl, cg*c2g*kl,
    // NLO
    pow(kt,3)*cg, kt*c2*cg, kt*pow(cg,2)*kl, cg*kt*c2g, 
    pow(kt*cg,2), c2*pow(cg,2), pow(cg,3)*kl, pow(cg,2)*c2g };

  if( order == "lo" ){
    vector<double> A_values_lo    = {35.0111,169.908,4.72866,2.38523,22.3288,-142.521,-22.996,47.2901,28.0101,-82.3576,-13.1345,31.2217,6.37158,-13.9821,-10.8268,};
    if( uncertantie=="muR_UP")A_values_lo    = {29.339989,142.828652,3.944472,1.987152,18.841645,-119.547333,-19.223360,39.572711,23.399797,-69.040228,-10.957083,26.085625,5.308951,-11.681204,-9.055525,};
    if( uncertantie=="muR_DN")A_values_lo    = {42.541636,205.687054,5.778218,2.919222,26.909724,-172.982077,-28.026622,57.564385,34.165875,-100.032567,-16.048800,38.075269,7.796551,-17.054354,-13.187661,};
    if( uncertantie=="muF_UP")A_values_lo    = {33.214886,159.640761,4.551589,2.305399,20.727429,-134.819116,-21.986524,45.071311,26.839155,-78.057553,-12.641570,29.899841,6.155038,-13.395502,-10.335953,};
    if( uncertantie=="muF_DN")A_values_lo    = {36.805367,180.705794,4.885372,2.452094,24.101164,-150.349033,-23.950697,49.443752,29.097654,-86.678003,-13.570472,32.448121,6.554115,-14.528044,-11.299634,};
    if( uncertantie=="muRF_UP")A_values_lo    = {27.832322,134.173116,3.796522,1.920497,17.483549,-113.072653,-18.378290,37.711670,22.418837,-65.428160,-10.544753,24.981035,5.128215,-11.189653,-8.643234,};
    if( uncertantie=="muRF_DN")A_values_lo    = {44.714854,218.701521,5.968957,3.000666,29.029593,-182.448410,-29.186181,60.175085,35.486974,-105.262740,-16.579853,39.569843,8.019070,-17.717366,-13.759752,};
    if( uncertantie=="PDF4LHC15_nlo_30_pdfas_UP")A_values_lo    = {36.208893,176.168489,4.879759,2.459933,23.284642,-147.500053,-23.752158,48.877622,28.922822,-85.183450,-13.550707,32.220345,6.571496,-14.439991,-11.195454,};
    if( uncertantie=="PDF4LHC15_nlo_30_pdfas_DN")A_values_lo    = {33.812066,163.645533,4.577366,2.310462,21.372602,-137.539060,-22.238752,45.701219,27.097318,-79.530065,-12.718469,30.222196,6.171466,-13.524284,-10.458108,};
    if( uncertantie=="PDF4LHC15_nlo_30_pdfas_as_UP")A_values_lo    = {34.557787,167.736680,4.665718,2.353319,22.038584,-140.685346,-22.693751,46.672758,27.642437,-81.296870,-12.961223,30.811765,6.286366,-13.797880,-10.684346,};
    if( uncertantie=="PDF4LHC15_nlo_30_pdfas_as_DN")A_values_lo    = {35.481751,172.242144,4.791610,2.416866,22.658916,-144.446484,-23.303070,47.923962,28.382176,-83.461994,-13.307699,31.634385,6.456062,-14.168891,-10.973232,};

    for(int i = 0; i < A_values_lo.size(); i++)
      answer_xsection += A_values_lo.at(i) * couplings.at(i);
    return answer_xsection;
  }
  if( order == "lo_MadGraph" ){ // pm_mg_LO-Ais-13TeV.txt
    vector<double> A_values_lo_mg = {27.7684,134.174,3.79328,1.35023,17.4785,-112.933,-18.3424,37.7066,19.5523,-65.5234,-9.08617,24.9774,4.31592,-11.2505,-7.59558,}; // pm_mg_LO-Ais-13TeV.txt
    for(int i = 0; i < A_values_lo_mg.size(); i++)
      answer_xsection += A_values_lo_mg.at(i) * couplings.at(i);
    return answer_xsection;
  }

  // TODO add NLO
  vector<double> A_values_nlo  = {62.5088,345.604,9.63451,4.34841,39.0143,-268.644,-44.2924,96.5595,53.515,-155.793,-23.678,54.5601,12.2273,-26.8654,-19.3723,-0.0904439,0.321092,0.452381,-0.0190758,-0.607163,1.27408,0.364487,-0.499263,};
  for(int i = 0; i < A_values_nlo.size(); i++)
    answer_xsection += A_values_nlo.at(i) * couplings.at(i);
  return answer_xsection;
}

double get_eft_xsec_13TeV(string mark, string order, string year = "2016", bool cms_fake=false){
  vector<double> answer = get_CMS_EFT_benchmarks(mark, year, cms_fake);
  return  get_eft_xsec_13TeV(answer, order);
}


class ReweightMandrik {
 public:
  vector< vector<double>> A_values_lo, A_values_nlo;

  ReweightMandrik(string input_lo="pm_pw_LO-Ais-13TeV_V2.txt", string input_nlo="pm_pw_NLO_Ais_13TeV_V2.txt"){
    LoadData( input_lo, A_values_lo );
    LoadData( input_nlo, A_values_nlo );
  }

  vector<double> GetEFTBenchmark(string name, string year="2016", bool cms_fake = false){
    vector<double> couplings = get_CMS_EFT_benchmarks( name, year, cms_fake );
    return couplings;
  }

  void LoadData( string input_name, vector< vector<double>> & answer ){
    std::string line;
    std::ifstream infile( input_name );
    if(!infile)cout<<"Reweight coefficient input file not found!\n";
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      answer.push_back( vector<double>() );
      string value;
      while (std::getline(iss, value, ',')) {
        answer[ answer.size()-1 ].push_back( atof(value.c_str()) );
      }
    }
  }
  
  vector< pair<double, double> > GetBinCenters(string order="nlo"){
    vector< vector<double> > * A_map = & A_values_nlo;
    if(order=="lo") A_map            = & A_values_lo;

    vector< pair<double, double> > answer;
    for(int i = 0; i < A_map->size(); i++){
      vector<double> & values = A_map->at(i);
      double mass_bin_center = (values.at(0)+values.at(1))*0.5;
      double cos_bin_center  = (values.at(2)+values.at(3))*0.5;
      answer.push_back( make_pair(mass_bin_center, cos_bin_center) );
    }
    return answer;
  }

  double GetDiffXsection(double m_hh, double cos_theta, const vector<double> & couplings, const vector< vector<double> > * A_map){
    for(int i = 0; i < A_map->size(); i++){
      const vector<double> & values = A_map->at(i);
      double mass_bin_end = values.at(1);
      double cos_bin_end  = values.at(3);
      if(m_hh > mass_bin_end or cos_theta > cos_bin_end) continue;
      double dXsec = 0;
      for( int j = 0; j < couplings.size(); j++ ){
	dXsec += values.at(j+4) * couplings.at(j);
      }
      return dXsec;
    }
    return 0;
  }

  double GetDiffXsection(double m_hh, double cos_theta, const vector<double> & eft_parameters, string order="nlo"){
    if(m_hh < 250) return 0;
    double kl, kt, c2, cg, c2g;
    kl = eft_parameters[0];
    kt = eft_parameters[1];
    c2 = eft_parameters[2];
    cg = eft_parameters[3];
    c2g = eft_parameters[4];

    vector<double> couplings;
    if(order == "lo") couplings = { // LO
	pow(kt,4), pow(c2,2), pow(kt,2)*pow(kl,2), pow(cg,2)*pow(kl,2),
	pow(c2g, 2), c2*pow(kt, 2), kl*pow(kt, 3), kt*kl*c2, cg*kl*c2, c2*c2g,
	cg*kl*pow(kt, 2), c2g*pow(kt, 2), 
	pow(kl,2)*cg*kt, c2g*kt*kl, cg*c2g*kl };
    if(order == "nlo") couplings = { // LO
	pow(kt,4), pow(c2,2), pow(kt,2)*pow(kl,2), pow(cg,2)*pow(kl,2),
	pow(c2g, 2), c2*pow(kt, 2), kl*pow(kt, 3), kt*kl*c2, cg*kl*c2, c2*c2g,
	cg*kl*pow(kt, 2), c2g*pow(kt, 2), 
	pow(kl,2)*cg*kt, c2g*kt*kl, cg*c2g*kl,
	// NLO
	pow(kt,3)*cg, kt*c2*cg, kt*pow(cg,2)*kl, cg*kt*c2g, 
	pow(kt*cg,2), c2*pow(cg,2), pow(cg,3)*kl, pow(cg,2)*c2g };
    
    vector< vector<double> > * A_map = & A_values_lo;
    if(order=="nlo")           A_map = & A_values_nlo;
    return GetDiffXsection(m_hh, cos_theta, couplings, A_map);
  }

  double GetA(int index, double m_hh, double cos_theta, string order="nlo"){
    vector< vector<double> > * A_map = & A_values_lo;
    if(order=="nlo")           A_map = & A_values_nlo;
    
    for(int i = 0; i < A_map->size(); i++){
      const vector<double> & values = A_map->at(i);
      double mass_bin_end = values.at(1);
      double cos_bin_end  = values.at(3);
      if(m_hh > mass_bin_end or cos_theta > cos_bin_end) continue;
      return values.at(index+4);
    }
    return 0;
  }

};

double make_reweight_prediction(string BMpoint, double event_mHH=0., double event_costhetaHH=0.){
  // Define the reweighter object                                                 
                                                     
  ReweightMandrik rm_pw_NLO = ReweightMandrik("data/scale_factor/HHEFT_reweight_histo/pm_pw_LO-Ais-13TeV_V2.txt", "data/scale_factor/HHEFT_reweight_histo/pm_pw_NLO_Ais_13TeV_V2.txt");

  // Load the 2D distribution of the input events                                                                   
                                                       
  TFile* inputfile =  new TFile("data/scale_factor/HHEFT_reweight_histo/2016/HHEFT_2016_skim.root");
  TH2F* histo_Nev = (TH2F*) inputfile->Get("Nev");
  double Nevtot = histo_Nev->Integral();

  // The target benchmark is for example the number 10                                                              
                                                       
  vector<double> couplings = rm_pw_NLO.GetEFTBenchmark(BMpoint);
  double XStot = get_eft_xsec_13TeV(BMpoint,"nlo");
                                                        
  // Compute the reweight
  double Nev = histo_Nev->GetBinContent( histo_Nev->FindBin(event_mHH, event_costhetaHH) );
  double XS = rm_pw_NLO.GetDiffXsection(event_mHH, event_costhetaHH, couplings, "nlo") / 1000.;//diff XS in [fb]
  
  // before using the differential XS, scale it by the bin area to properly compare it with histo_Nev conten
                                                     
  int ibinmhh = histo_Nev->GetXaxis()->FindBin(event_mHH);
  int ibincosthetaHH = histo_Nev->GetYaxis()->FindBin(event_costhetaHH);
  double Noutputev = XS * histo_Nev->GetXaxis()->GetBinWidth(ibinmhh) * histo_Nev->GetYaxis()->GetBinWidth(ibincosthetaHH);
  double reweight = Noutputev/Nev * Nevtot/XStot;
  //cout<<Noutputev<<" "<<Nev<<" "<<Nevtot<<" "<<XS<<" "<<XStot<<" "<<reweight<<endl;
  return reweight;
}
#endif
