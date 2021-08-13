#ifndef hhtree_HH
#include <TBranch.h> 
#include <TTree.h> 
#include <TH1F.h> 
#include <TFile.h> 
#include <TBits.h> 
#include <vector> 
#include <unistd.h> 
using namespace std; 

class hhtree{
 private:
 protected:
   unsigned int index;
   unsigned int run_;
   TBranch *run_branch;
   bool run_isLoaded;
   unsigned int luminosityBlock_;
   TBranch *luminosityBlock_branch;
   bool luminosityBlock_isLoaded;
   unsigned long int event_;
   TBranch *event_branch;
   bool event_isLoaded;
   float genWeight_;
   TBranch *genWeight_branch;
   bool genWeight_isLoaded;
   int nLHEPdfWeight_;
   TBranch *nLHEPdfWeight_branch;
   bool nLHEPdfWeight_isLoaded;
   float LHEPdfWeight_[103];
   TBranch *LHEPdfWeight_branch;
   bool LHEPdfWeight_isLoaded;
   int nLHEScaleWeight_;
   TBranch *nLHEScaleWeight_branch;
   bool nLHEScaleWeight_isLoaded;
   float LHEScaleWeight_[9];
   TBranch *LHEScaleWeight_branch;
   bool LHEScaleWeight_isLoaded;
   int nPSWeight_;
   TBranch *nPSWeight_branch;
   bool nPSWeight_isLoaded;
   float PSWeight_[4];
   TBranch *PSWeight_branch;
   bool PSWeight_isLoaded;

  bool HLT_Ele27_WPTight_Gsf_;
  TBranch *HLT_Ele27_WPTight_Gsf_branch;
  bool HLT_Ele27_WPTight_Gsf_isLoaded;
  bool HLT_Ele28_WPTight_Gsf_;
  TBranch *HLT_Ele28_WPTight_Gsf_branch;
  bool HLT_Ele28_WPTight_Gsf_isLoaded;
  bool HLT_Ele30_WPTight_Gsf_;
  TBranch *HLT_Ele30_WPTight_Gsf_branch;
  bool HLT_Ele30_WPTight_Gsf_isLoaded;
  bool HLT_Ele32_WPTight_Gsf_;
  TBranch *HLT_Ele32_WPTight_Gsf_branch;
  bool HLT_Ele32_WPTight_Gsf_isLoaded;
  bool HLT_Ele35_WPTight_Gsf_;
  TBranch *HLT_Ele35_WPTight_Gsf_branch;
  bool HLT_Ele35_WPTight_Gsf_isLoaded;
  bool HLT_Ele38_WPTight_Gsf_;
  TBranch *HLT_Ele38_WPTight_Gsf_branch;
  bool HLT_Ele38_WPTight_Gsf_isLoaded;
  bool HLT_Ele40_WPTight_Gsf_;
  TBranch *HLT_Ele40_WPTight_Gsf_branch;
  bool HLT_Ele40_WPTight_Gsf_isLoaded;
  bool HLT_IsoMu20_;
  TBranch *HLT_IsoMu20_branch;
  bool HLT_IsoMu20_isLoaded;
  bool HLT_IsoMu24_;
  TBranch *HLT_IsoMu24_branch;
  bool HLT_IsoMu24_isLoaded;
  bool HLT_IsoMu24_eta2p1_;
  TBranch *HLT_IsoMu24_eta2p1_branch;
  bool HLT_IsoMu24_eta2p1_isLoaded;
  bool HLT_IsoMu27_;
  TBranch *HLT_IsoMu27_branch;
  bool HLT_IsoMu27_isLoaded;
  bool HLT_IsoMu30_;
  TBranch *HLT_IsoMu30_branch;
  bool HLT_IsoMu30_isLoaded;
  bool HLT_Mu50_;
  TBranch *HLT_Mu50_branch;
  bool HLT_Mu50_isLoaded;
  bool HLT_Mu55_;
  TBranch *HLT_Mu55_branch;
  bool HLT_Mu55_isLoaded;
  bool HLT_Photon175_;
  TBranch *HLT_Photon175_branch;
  bool HLT_Photon175_isLoaded;
  bool HLT_PFHT780_;
  TBranch *HLT_PFHT780_branch;
  bool HLT_PFHT780_isLoaded;
  bool HLT_PFHT890_;
  TBranch *HLT_PFHT890_branch;
  bool HLT_PFHT890_isLoaded;
  bool HLT_PFHT1050_;
  TBranch *HLT_PFHT1050_branch;
  bool HLT_PFHT1050_isLoaded;
  bool HLT_AK8PFJet360_TrimMass30_=false;
  TBranch *HLT_AK8PFJet360_TrimMass30_branch;
  bool HLT_AK8PFJet360_TrimMass30_isLoaded;
  bool HLT_AK8PFJet380_TrimMass30_=false;
  TBranch *HLT_AK8PFJet380_TrimMass30_branch;
  bool HLT_AK8PFJet380_TrimMass30_isLoaded;
  bool HLT_AK8PFJet400_TrimMass30_=false;
  TBranch *HLT_AK8PFJet400_TrimMass30_branch;
  bool HLT_AK8PFJet400_TrimMass30_isLoaded;
  bool HLT_AK8PFJet420_TrimMass30_;
  TBranch *HLT_AK8PFJet420_TrimMass30_branch;
  bool HLT_AK8PFJet420_TrimMass30_isLoaded;
  bool HLT_AK8PFHT750_TrimMass50_=false;
  TBranch *HLT_AK8PFHT750_TrimMass50_branch;
  bool HLT_AK8PFHT750_TrimMass50_isLoaded;
  bool HLT_AK8PFHT800_TrimMass50_=false;
  TBranch *HLT_AK8PFHT800_TrimMass50_branch;
  bool HLT_AK8PFHT800_TrimMass50_isLoaded;
  bool HLT_AK8PFHT850_TrimMass50_;
  TBranch *HLT_AK8PFHT850_TrimMass50_branch;
  bool HLT_AK8PFHT850_TrimMass50_isLoaded;
  bool HLT_AK8PFHT900_TrimMass50_;
  TBranch *HLT_AK8PFHT900_TrimMass50_branch;
  bool HLT_AK8PFHT900_TrimMass50_isLoaded;
  bool HLT_PFJet450_;
  TBranch *HLT_PFJet450_branch;
  bool HLT_PFJet450_isLoaded;
  bool HLT_PFJet500_;
  TBranch *HLT_PFJet500_branch;
  bool HLT_PFJet500_isLoaded;
  bool HLT_PFJet550_;
  TBranch *HLT_PFJet550_branch;
  bool HLT_PFJet550_isLoaded;
  bool HLT_AK8PFJet400_;
  TBranch *HLT_AK8PFJet400_branch;
  bool HLT_AK8PFJet400_isLoaded;
  bool HLT_AK8PFJet450_;
  TBranch *HLT_AK8PFJet450_branch;
  bool HLT_AK8PFJet450_isLoaded;
  bool HLT_AK8PFJet500_;
  TBranch *HLT_AK8PFJet500_branch;
  bool HLT_AK8PFJet500_isLoaded;
  bool HLT_AK8PFJet550_;
  TBranch *HLT_AK8PFJet550_branch;
  bool HLT_AK8PFJet550_isLoaded;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_;
  TBranch *HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_branch;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_isLoaded;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_;
  TBranch *HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_branch;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_isLoaded;
  bool HLT_AK8PFJet330_PFAK8BTagCSV_p17_=false;
  TBranch *HLT_AK8PFJet330_PFAK8BTagCSV_p17_branch;
  bool HLT_AK8PFJet330_PFAK8BTagCSV_p17_isLoaded;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_;
  TBranch *HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_branch;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_isLoaded;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_;
  TBranch *HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_branch;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_isLoaded;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_=false;
  TBranch *HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_branch;
  bool HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_isLoaded;
  bool HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_;
  TBranch *HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_branch;
  bool HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_isLoaded;
  bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_;
  TBranch *HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_branch;
  bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_isLoaded;
  bool HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_;
  TBranch *HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_branch;
  bool HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_isLoaded;
  bool HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_;
  TBranch *HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_branch;
  bool HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_isLoaded;
  bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_;
  TBranch *HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_branch;
  bool HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_isLoaded;
  bool HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_;
  TBranch *HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_branch;
  bool HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_isLoaded;
   float weight_;
   TBranch *weight_branch;
   bool weight_isLoaded;
   float met_;
   TBranch *met_branch;
   bool met_isLoaded;
   float metphi_;
   TBranch *metphi_branch;
   bool metphi_isLoaded;
   float ht_;
   TBranch *ht_branch;
   bool ht_isLoaded;
   bool passmetfilters_;
   TBranch *passmetfilters_branch;
   bool passmetfilters_isLoaded;
   float l1PreFiringWeight_;
   TBranch *l1PreFiringWeight_branch;
   bool l1PreFiringWeight_isLoaded;
   float l1PreFiringWeightUp_;
   TBranch *l1PreFiringWeightUp_branch;
   bool l1PreFiringWeightUp_isLoaded;
   float l1PreFiringWeightDown_;
   TBranch *l1PreFiringWeightDown_branch;
   bool l1PreFiringWeightDown_isLoaded;
   float triggerEffWeight_;
   TBranch *triggerEffWeight_branch;
   bool triggerEffWeight_isLoaded;
   float triggerEff3DWeight_;
   TBranch *triggerEff3DWeight_branch;
   bool triggerEff3DWeight_isLoaded;
   float triggerEffMCWeight_;
   TBranch *triggerEffMCWeight_branch;
   bool triggerEffMCWeight_isLoaded;
   float triggerEffMC3DWeight_;
   TBranch *triggerEffMC3DWeight_branch;
   bool triggerEffMC3DWeight_isLoaded;
   float fatJet1Pt_;
   TBranch *fatJet1Pt_branch;
   bool fatJet1Pt_isLoaded;
   float fatJet1Eta_;
   TBranch *fatJet1Eta_branch;
   bool fatJet1Eta_isLoaded;
   float fatJet1Phi_;
   TBranch *fatJet1Phi_branch;
   bool fatJet1Phi_isLoaded;
   float fatJet1Mass_;
   TBranch *fatJet1Mass_branch;
   bool fatJet1Mass_isLoaded;
   float fatJet1MassSD_;
   TBranch *fatJet1MassSD_branch;
   bool fatJet1MassSD_isLoaded;
   float fatJet1MassSD_UnCorrected_;
   TBranch *fatJet1MassSD_UnCorrected_branch;
   bool fatJet1MassSD_UnCorrected_isLoaded;
   float fatJet1MassRegressed_;
   TBranch *fatJet1MassRegressed_branch;
   bool fatJet1MassRegressed_isLoaded;
   float fatJet1MassRegressed_UnCorrected_;
   TBranch *fatJet1MassRegressed_UnCorrected_branch;
   bool fatJet1MassRegressed_UnCorrected_isLoaded;
   float fatJet1PNetXbb_;
   TBranch *fatJet1PNetXbb_branch;
   bool fatJet1PNetXbb_isLoaded;
   float fatJet1PNetQCDb_;
   TBranch *fatJet1PNetQCDb_branch;
   bool fatJet1PNetQCDb_isLoaded;
   float fatJet1PNetQCDbb_;
   TBranch *fatJet1PNetQCDbb_branch;
   bool fatJet1PNetQCDbb_isLoaded;
   float fatJet1PNetQCDc_;
   TBranch *fatJet1PNetQCDc_branch;
   bool fatJet1PNetQCDc_isLoaded;
   float fatJet1PNetQCDcc_;
   TBranch *fatJet1PNetQCDcc_branch;
   bool fatJet1PNetQCDcc_isLoaded;
   float fatJet1PNetQCDothers_;
   TBranch *fatJet1PNetQCDothers_branch;
   bool fatJet1PNetQCDothers_isLoaded;
   float fatJet1Tau3OverTau2_;
   TBranch *fatJet1Tau3OverTau2_branch;
   bool fatJet1Tau3OverTau2_isLoaded;
   int fatJet1GenMatchIndex_;
   TBranch *fatJet1GenMatchIndex_branch;
   bool fatJet1GenMatchIndex_isLoaded;
   bool fatJet1HasMuon_;
   TBranch *fatJet1HasMuon_branch;
   bool fatJet1HasMuon_isLoaded;
   bool fatJet1HasElectron_;
   TBranch *fatJet1HasElectron_branch;
   bool fatJet1HasElectron_isLoaded;
   bool fatJet1HasBJetCSVLoose_;
   TBranch *fatJet1HasBJetCSVLoose_branch;
   bool fatJet1HasBJetCSVLoose_isLoaded;
   bool fatJet1HasBJetCSVMedium_;
   TBranch *fatJet1HasBJetCSVMedium_branch;
   bool fatJet1HasBJetCSVMedium_isLoaded;
   bool fatJet1HasBJetCSVTight_;
   TBranch *fatJet1HasBJetCSVTight_branch;
   bool fatJet1HasBJetCSVTight_isLoaded;
   bool fatJet1OppositeHemisphereHasBJet_;
   TBranch *fatJet1OppositeHemisphereHasBJet_branch;
   bool fatJet1OppositeHemisphereHasBJet_isLoaded;
   float fatJet1PtOverMHH_;
   TBranch *fatJet1PtOverMHH_branch;
   bool fatJet1PtOverMHH_isLoaded;
   float fatJet1PtOverMSD_;
   TBranch *fatJet1PtOverMSD_branch;
   bool fatJet1PtOverMSD_isLoaded;
   float fatJet1PtOverMRegressed_;
   TBranch *fatJet1PtOverMRegressed_branch;
   bool fatJet1PtOverMRegressed_isLoaded;
   float fatJet1MassSD_JMS_Down_;
   TBranch *fatJet1MassSD_JMS_Down_branch;
   bool fatJet1MassSD_JMS_Down_isLoaded;
   float fatJet1MassSD_JMS_Up_;
   TBranch *fatJet1MassSD_JMS_Up_branch;
   bool fatJet1MassSD_JMS_Up_isLoaded;
   float fatJet1MassSD_JMR_Down_;
   TBranch *fatJet1MassSD_JMR_Down_branch;
   bool fatJet1MassSD_JMR_Down_isLoaded;
   float fatJet1MassSD_JMR_Up_;
   TBranch *fatJet1MassSD_JMR_Up_branch;
   bool fatJet1MassSD_JMR_Up_isLoaded;
   float fatJet1MassRegressed_JMS_Down_;
   TBranch *fatJet1MassRegressed_JMS_Down_branch;
   bool fatJet1MassRegressed_JMS_Down_isLoaded;
   float fatJet1MassRegressed_JMS_Up_;
   TBranch *fatJet1MassRegressed_JMS_Up_branch;
   bool fatJet1MassRegressed_JMS_Up_isLoaded;
   float fatJet1MassRegressed_JMR_Down_;
   TBranch *fatJet1MassRegressed_JMR_Down_branch;
   bool fatJet1MassRegressed_JMR_Down_isLoaded;
   float fatJet1MassRegressed_JMR_Up_;
   TBranch *fatJet1MassRegressed_JMR_Up_branch;
   bool fatJet1MassRegressed_JMR_Up_isLoaded;
   float fatJet1PtOverMHH_JMS_Down_;
   TBranch *fatJet1PtOverMHH_JMS_Down_branch;
   bool fatJet1PtOverMHH_JMS_Down_isLoaded;
   float fatJet1PtOverMHH_JMS_Up_;
   TBranch *fatJet1PtOverMHH_JMS_Up_branch;
   bool fatJet1PtOverMHH_JMS_Up_isLoaded;
   float fatJet1PtOverMHH_JMR_Down_;
   TBranch *fatJet1PtOverMHH_JMR_Down_branch;
   bool fatJet1PtOverMHH_JMR_Down_isLoaded;
   float fatJet1PtOverMHH_JMR_Up_;
   TBranch *fatJet1PtOverMHH_JMR_Up_branch;
   bool fatJet1PtOverMHH_JMR_Up_isLoaded;
   float fatJet1Pt_JERUp_;
   TBranch *fatJet1Pt_JERUp_branch;
   bool fatJet1Pt_JERUp_isLoaded;
   float fatJet1PtOverMHH_JERUp_;
   TBranch *fatJet1PtOverMHH_JERUp_branch;
   bool fatJet1PtOverMHH_JERUp_isLoaded;
   float fatJet1Pt_JERDown_;
   TBranch *fatJet1Pt_JERDown_branch;
   bool fatJet1Pt_JERDown_isLoaded;
   float fatJet1PtOverMHH_JERDown_;
   TBranch *fatJet1PtOverMHH_JERDown_branch;
   bool fatJet1PtOverMHH_JERDown_isLoaded;
   float fatJet1Pt_JESUp_;
   TBranch *fatJet1Pt_JESUp_branch;
   bool fatJet1Pt_JESUp_isLoaded;
   float fatJet1PtOverMHH_JESUp_;
   TBranch *fatJet1PtOverMHH_JESUp_branch;
   bool fatJet1PtOverMHH_JESUp_isLoaded;
   float fatJet1Pt_JESDown_;
   TBranch *fatJet1Pt_JESDown_branch;
   bool fatJet1Pt_JESDown_isLoaded;
   float fatJet1PtOverMHH_JESDown_;
   TBranch *fatJet1PtOverMHH_JESDown_branch;
   bool fatJet1PtOverMHH_JESDown_isLoaded;
   float fatJet2Pt_;
   TBranch *fatJet2Pt_branch;
   bool fatJet2Pt_isLoaded;
   float fatJet2Eta_;
   TBranch *fatJet2Eta_branch;
   bool fatJet2Eta_isLoaded;
   float fatJet2Phi_;
   TBranch *fatJet2Phi_branch;
   bool fatJet2Phi_isLoaded;
   float fatJet2Mass_;
   TBranch *fatJet2Mass_branch;
   bool fatJet2Mass_isLoaded;
   float fatJet2MassSD_;
   TBranch *fatJet2MassSD_branch;
   bool fatJet2MassSD_isLoaded;
   float fatJet2MassSD_UnCorrected_;
   TBranch *fatJet2MassSD_UnCorrected_branch;
   bool fatJet2MassSD_UnCorrected_isLoaded;
   float fatJet2MassRegressed_;
   TBranch *fatJet2MassRegressed_branch;
   bool fatJet2MassRegressed_isLoaded;
   float fatJet2MassRegressed_UnCorrected_;
   TBranch *fatJet2MassRegressed_UnCorrected_branch;
   bool fatJet2MassRegressed_UnCorrected_isLoaded;
   float fatJet2PNetXbb_;
   TBranch *fatJet2PNetXbb_branch;
   bool fatJet2PNetXbb_isLoaded;
   float fatJet2PNetQCDb_;
   TBranch *fatJet2PNetQCDb_branch;
   bool fatJet2PNetQCDb_isLoaded;
   float fatJet2PNetQCDbb_;
   TBranch *fatJet2PNetQCDbb_branch;
   bool fatJet2PNetQCDbb_isLoaded;
   float fatJet2PNetQCDc_;
   TBranch *fatJet2PNetQCDc_branch;
   bool fatJet2PNetQCDc_isLoaded;
   float fatJet2PNetQCDcc_;
   TBranch *fatJet2PNetQCDcc_branch;
   bool fatJet2PNetQCDcc_isLoaded;
   float fatJet2PNetQCDothers_;
   TBranch *fatJet2PNetQCDothers_branch;
   bool fatJet2PNetQCDothers_isLoaded;
   float fatJet2Tau3OverTau2_;
   TBranch *fatJet2Tau3OverTau2_branch;
   bool fatJet2Tau3OverTau2_isLoaded;
   int fatJet2GenMatchIndex_;
   TBranch *fatJet2GenMatchIndex_branch;
   bool fatJet2GenMatchIndex_isLoaded;
   bool fatJet2HasMuon_;
   TBranch *fatJet2HasMuon_branch;
   bool fatJet2HasMuon_isLoaded;
   bool fatJet2HasElectron_;
   TBranch *fatJet2HasElectron_branch;
   bool fatJet2HasElectron_isLoaded;
   bool fatJet2HasBJetCSVLoose_;
   TBranch *fatJet2HasBJetCSVLoose_branch;
   bool fatJet2HasBJetCSVLoose_isLoaded;
   bool fatJet2HasBJetCSVMedium_;
   TBranch *fatJet2HasBJetCSVMedium_branch;
   bool fatJet2HasBJetCSVMedium_isLoaded;
   bool fatJet2HasBJetCSVTight_;
   TBranch *fatJet2HasBJetCSVTight_branch;
   bool fatJet2HasBJetCSVTight_isLoaded;
   bool fatJet2OppositeHemisphereHasBJet_;
   TBranch *fatJet2OppositeHemisphereHasBJet_branch;
   bool fatJet2OppositeHemisphereHasBJet_isLoaded;
   float fatJet2PtOverMHH_;
   TBranch *fatJet2PtOverMHH_branch;
   bool fatJet2PtOverMHH_isLoaded;
   float fatJet2PtOverMSD_;
   TBranch *fatJet2PtOverMSD_branch;
   bool fatJet2PtOverMSD_isLoaded;
   float fatJet2PtOverMRegressed_;
   TBranch *fatJet2PtOverMRegressed_branch;
   bool fatJet2PtOverMRegressed_isLoaded;
   float fatJet2MassSD_JMS_Down_;
   TBranch *fatJet2MassSD_JMS_Down_branch;
   bool fatJet2MassSD_JMS_Down_isLoaded;
   float fatJet2MassSD_JMS_Up_;
   TBranch *fatJet2MassSD_JMS_Up_branch;
   bool fatJet2MassSD_JMS_Up_isLoaded;
   float fatJet2MassSD_JMR_Down_;
   TBranch *fatJet2MassSD_JMR_Down_branch;
   bool fatJet2MassSD_JMR_Down_isLoaded;
   float fatJet2MassSD_JMR_Up_;
   TBranch *fatJet2MassSD_JMR_Up_branch;
   bool fatJet2MassSD_JMR_Up_isLoaded;
   float fatJet2MassRegressed_JMS_Down_;
   TBranch *fatJet2MassRegressed_JMS_Down_branch;
   bool fatJet2MassRegressed_JMS_Down_isLoaded;
   float fatJet2MassRegressed_JMS_Up_;
   TBranch *fatJet2MassRegressed_JMS_Up_branch;
   bool fatJet2MassRegressed_JMS_Up_isLoaded;
   float fatJet2MassRegressed_JMR_Down_;
   TBranch *fatJet2MassRegressed_JMR_Down_branch;
   bool fatJet2MassRegressed_JMR_Down_isLoaded;
   float fatJet2MassRegressed_JMR_Up_;
   TBranch *fatJet2MassRegressed_JMR_Up_branch;
   bool fatJet2MassRegressed_JMR_Up_isLoaded;
   float fatJet2PtOverMHH_JMS_Down_;
   TBranch *fatJet2PtOverMHH_JMS_Down_branch;
   bool fatJet2PtOverMHH_JMS_Down_isLoaded;
   float fatJet2PtOverMHH_JMS_Up_;
   TBranch *fatJet2PtOverMHH_JMS_Up_branch;
   bool fatJet2PtOverMHH_JMS_Up_isLoaded;
   float fatJet2PtOverMHH_JMR_Down_;
   TBranch *fatJet2PtOverMHH_JMR_Down_branch;
   bool fatJet2PtOverMHH_JMR_Down_isLoaded;
   float fatJet2PtOverMHH_JMR_Up_;
   TBranch *fatJet2PtOverMHH_JMR_Up_branch;
   bool fatJet2PtOverMHH_JMR_Up_isLoaded;
   float fatJet2Pt_JERUp_;
   TBranch *fatJet2Pt_JERUp_branch;
   bool fatJet2Pt_JERUp_isLoaded;
   float fatJet2PtOverMHH_JERUp_;
   TBranch *fatJet2PtOverMHH_JERUp_branch;
   bool fatJet2PtOverMHH_JERUp_isLoaded;
   float fatJet2Pt_JERDown_;
   TBranch *fatJet2Pt_JERDown_branch;
   bool fatJet2Pt_JERDown_isLoaded;
   float fatJet2PtOverMHH_JERDown_;
   TBranch *fatJet2PtOverMHH_JERDown_branch;
   bool fatJet2PtOverMHH_JERDown_isLoaded;
   float fatJet2Pt_JESUp_;
   TBranch *fatJet2Pt_JESUp_branch;
   bool fatJet2Pt_JESUp_isLoaded;
   float fatJet2PtOverMHH_JESUp_;
   TBranch *fatJet2PtOverMHH_JESUp_branch;
   bool fatJet2PtOverMHH_JESUp_isLoaded;
   float fatJet2Pt_JESDown_;
   TBranch *fatJet2Pt_JESDown_branch;
   bool fatJet2Pt_JESDown_isLoaded;
   float fatJet2PtOverMHH_JESDown_;
   TBranch *fatJet2PtOverMHH_JESDown_branch;
   bool fatJet2PtOverMHH_JESDown_isLoaded;
   float hh_pt_;
   TBranch *hh_pt_branch;
   bool hh_pt_isLoaded;
   float hh_eta_;
   TBranch *hh_eta_branch;
   bool hh_eta_isLoaded;
   float hh_phi_;
   TBranch *hh_phi_branch;
   bool hh_phi_isLoaded;
   float hh_mass_;
   TBranch *hh_mass_branch;
   bool hh_mass_isLoaded;
   float hh_pt_JMS_Down_;
   TBranch *hh_pt_JMS_Down_branch;
   bool hh_pt_JMS_Down_isLoaded;
   float hh_pt_JMS_Up_;
   TBranch *hh_pt_JMS_Up_branch;
   bool hh_pt_JMS_Up_isLoaded;
   float hh_eta_JMS_Down_;
   TBranch *hh_eta_JMS_Down_branch;
   bool hh_eta_JMS_Down_isLoaded;
   float hh_eta_JMS_Up_;
   TBranch *hh_eta_JMS_Up_branch;
   bool hh_eta_JMS_Up_isLoaded;
   float hh_mass_JMS_Down_;
   TBranch *hh_mass_JMS_Down_branch;
   bool hh_mass_JMS_Down_isLoaded;
   float hh_mass_JMS_Up_;
   TBranch *hh_mass_JMS_Up_branch;
   bool hh_mass_JMS_Up_isLoaded;
   float hh_pt_JMR_Down_;
   TBranch *hh_pt_JMR_Down_branch;
   bool hh_pt_JMR_Down_isLoaded;
   float hh_pt_JMR_Up_;
   TBranch *hh_pt_JMR_Up_branch;
   bool hh_pt_JMR_Up_isLoaded;
   float hh_eta_JMR_Down_;
   TBranch *hh_eta_JMR_Down_branch;
   bool hh_eta_JMR_Down_isLoaded;
   float hh_eta_JMR_Up_;
   TBranch *hh_eta_JMR_Up_branch;
   bool hh_eta_JMR_Up_isLoaded;
   float hh_mass_JMR_Down_;
   TBranch *hh_mass_JMR_Down_branch;
   bool hh_mass_JMR_Down_isLoaded;
   float hh_mass_JMR_Up_;
   TBranch *hh_mass_JMR_Up_branch;
   bool hh_mass_JMR_Up_isLoaded;
   float hh_pt_JERUp_;
   TBranch *hh_pt_JERUp_branch;
   bool hh_pt_JERUp_isLoaded;
   float hh_eta_JERUp_;
   TBranch *hh_eta_JERUp_branch;
   bool hh_eta_JERUp_isLoaded;
   float hh_mass_JERUp_;
   TBranch *hh_mass_JERUp_branch;
   bool hh_mass_JERUp_isLoaded;
   float hh_pt_JERDown_;
   TBranch *hh_pt_JERDown_branch;
   bool hh_pt_JERDown_isLoaded;
   float hh_eta_JERDown_;
   TBranch *hh_eta_JERDown_branch;
   bool hh_eta_JERDown_isLoaded;
   float hh_mass_JERDown_;
   TBranch *hh_mass_JERDown_branch;
   bool hh_mass_JERDown_isLoaded;
   float hh_pt_JESUp_;
   TBranch *hh_pt_JESUp_branch;
   bool hh_pt_JESUp_isLoaded;
   float hh_eta_JESUp_;
   TBranch *hh_eta_JESUp_branch;
   bool hh_eta_JESUp_isLoaded;
   float hh_mass_JESUp_;
   TBranch *hh_mass_JESUp_branch;
   bool hh_mass_JESUp_isLoaded;
   float hh_pt_JESDown_;
   TBranch *hh_pt_JESDown_branch;
   bool hh_pt_JESDown_isLoaded;
   float hh_eta_JESDown_;
   TBranch *hh_eta_JESDown_branch;
   bool hh_eta_JESDown_isLoaded;
   float hh_mass_JESDown_;
   TBranch *hh_mass_JESDown_branch;
   bool hh_mass_JESDown_isLoaded;
   float deltaEta_j1j2_;
   TBranch *deltaEta_j1j2_branch;
   bool deltaEta_j1j2_isLoaded;
   float deltaPhi_j1j2_;
   TBranch *deltaPhi_j1j2_branch;
   bool deltaPhi_j1j2_isLoaded;
   float deltaR_j1j2_;
   TBranch *deltaR_j1j2_branch;
   bool deltaR_j1j2_isLoaded;
   float ptj2_over_ptj1_;
   TBranch *ptj2_over_ptj1_branch;
   bool ptj2_over_ptj1_isLoaded;
   float mj2_over_mj1_;
   TBranch *mj2_over_mj1_branch;
   bool mj2_over_mj1_isLoaded;
   int isVBFtag_;
   TBranch *isVBFtag_branch;
   bool isVBFtag_isLoaded;
   float dijetmass_;
   TBranch *dijetmass_branch;
   bool dijetmass_isLoaded;
   float vbfjet1Pt_;
   TBranch *vbfjet1Pt_branch;
   bool vbfjet1Pt_isLoaded;
   float vbfjet1Eta_;
   TBranch *vbfjet1Eta_branch;
   bool vbfjet1Eta_isLoaded;
   float vbfjet1Phi_;
   TBranch *vbfjet1Phi_branch;
   bool vbfjet1Phi_isLoaded;
   float vbfjet1Mass_;
   TBranch *vbfjet1Mass_branch;
   bool vbfjet1Mass_isLoaded;
   float vbffatJet1Pt_;
   TBranch *vbffatJet1Pt_branch;
   bool vbffatJet1Pt_isLoaded;
   float vbffatJet1Eta_;
   TBranch *vbffatJet1Eta_branch;
   bool vbffatJet1Eta_isLoaded;
   float vbffatJet1Phi_;
   TBranch *vbffatJet1Phi_branch;
   bool vbffatJet1Phi_isLoaded;
   float vbffatJet1PNetXbb_;
   TBranch *vbffatJet1PNetXbb_branch;
   bool vbffatJet1PNetXbb_isLoaded;
   float vbfjet2Pt_;
   TBranch *vbfjet2Pt_branch;
   bool vbfjet2Pt_isLoaded;
   float vbfjet2Eta_;
   TBranch *vbfjet2Eta_branch;
   bool vbfjet2Eta_isLoaded;
   float vbfjet2Phi_;
   TBranch *vbfjet2Phi_branch;
   bool vbfjet2Phi_isLoaded;
   float vbfjet2Mass_;
   TBranch *vbfjet2Mass_branch;
   bool vbfjet2Mass_isLoaded;
   float vbffatJet2Pt_;
   TBranch *vbffatJet2Pt_branch;
   bool vbffatJet2Pt_isLoaded;
   float vbffatJet2Eta_;
   TBranch *vbffatJet2Eta_branch;
   bool vbffatJet2Eta_isLoaded;
   float vbffatJet2Phi_;
   TBranch *vbffatJet2Phi_branch;
   bool vbffatJet2Phi_isLoaded;
   float vbffatJet2PNetXbb_;
   TBranch *vbffatJet2PNetXbb_branch;
   bool vbffatJet2PNetXbb_isLoaded;
   float jet1Pt_;
   TBranch *jet1Pt_branch;
   bool jet1Pt_isLoaded;
   float jet1Eta_;
   TBranch *jet1Eta_branch;
   bool jet1Eta_isLoaded;
   float jet1Phi_;
   TBranch *jet1Phi_branch;
   bool jet1Phi_isLoaded;
   float jet2Pt_;
   TBranch *jet2Pt_branch;
   bool jet2Pt_isLoaded;
   float jet2Eta_;
   TBranch *jet2Eta_branch;
   bool jet2Eta_isLoaded;
   float jet2Phi_;
   TBranch *jet2Phi_branch;
   bool jet2Phi_isLoaded;
   float jet3Pt_;
   TBranch *jet3Pt_branch;
   bool jet3Pt_isLoaded;
   float jet3Eta_;
   TBranch *jet3Eta_branch;
   bool jet3Eta_isLoaded;
   float jet3Phi_;
   TBranch *jet3Phi_branch;
   bool jet3Phi_isLoaded;
   float jet4Pt_;
   TBranch *jet4Pt_branch;
   bool jet4Pt_isLoaded;
   float jet4Eta_;
   TBranch *jet4Eta_branch;
   bool jet4Eta_isLoaded;
   float jet4Phi_;
   TBranch *jet4Phi_branch;
   bool jet4Phi_isLoaded;
   float lep1Pt_;
   TBranch *lep1Pt_branch;
   bool lep1Pt_isLoaded;
   float lep1Eta_;
   TBranch *lep1Eta_branch;
   bool lep1Eta_isLoaded;
   float lep1Phi_;
   TBranch *lep1Phi_branch;
   bool lep1Phi_isLoaded;
   int lep1Id_;
   TBranch *lep1Id_branch;
   bool lep1Id_isLoaded;
   float lep2Pt_;
   TBranch *lep2Pt_branch;
   bool lep2Pt_isLoaded;
   float lep2Eta_;
   TBranch *lep2Eta_branch;
   bool lep2Eta_isLoaded;
   float lep2Phi_;
   TBranch *lep2Phi_branch;
   bool lep2Phi_isLoaded;
   int lep2Id_;
   TBranch *lep2Id_branch;
   bool lep2Id_isLoaded;
   float genHiggs1Pt_;
   TBranch *genHiggs1Pt_branch;
   bool genHiggs1Pt_isLoaded;
   float genHiggs1Eta_;
   TBranch *genHiggs1Eta_branch;
   bool genHiggs1Eta_isLoaded;
   float genHiggs1Phi_;
   TBranch *genHiggs1Phi_branch;
   bool genHiggs1Phi_isLoaded;
   float genHiggs2Pt_;
   TBranch *genHiggs2Pt_branch;
   bool genHiggs2Pt_isLoaded;
   float genHiggs2Eta_;
   TBranch *genHiggs2Eta_branch;
   bool genHiggs2Eta_isLoaded;
   float genHiggs2Phi_;
   TBranch *genHiggs2Phi_branch;
   bool genHiggs2Phi_isLoaded;
   float puWeight_;
   TBranch *puWeight_branch;
   bool puWeight_isLoaded;
   float puWeightUp_;
   TBranch *puWeightUp_branch;
   bool puWeightUp_isLoaded;
   float puWeightDown_;
   TBranch *puWeightDown_branch;
   bool puWeightDown_isLoaded;
   float xsecWeight_;
   TBranch *xsecWeight_branch;
   bool xsecWeight_isLoaded;
   float LHEScaleWeightNorm_[9];
   TBranch *LHEScaleWeightNorm_branch;
   bool LHEScaleWeightNorm_isLoaded;
   float LHEPdfWeightNorm_[103];
   TBranch *LHEPdfWeightNorm_branch;
   bool LHEPdfWeightNorm_isLoaded;
   float disc_qcd_and_ttbar_Run2_enhanced_v8p2_;
   TBranch *disc_qcd_and_ttbar_Run2_enhanced_v8p2_branch;
   bool disc_qcd_and_ttbar_Run2_enhanced_v8p2_isLoaded;
   float disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_;
   TBranch *disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_branch;
   bool disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_isLoaded;
   float disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_;
   TBranch *disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_branch;
   bool disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_isLoaded;
   float disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_;
   TBranch *disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_branch;
   bool disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_isLoaded;
   float disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_;
   TBranch *disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_branch;
   bool disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_isLoaded;
   float disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_;
   TBranch *disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_branch;
   bool disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_isLoaded;
   float disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_;
   TBranch *disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_branch;
   bool disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_isLoaded;
 
 
 public:
   void Init(TTree * tree);
   void GetEntry(unsigned int idx);
   void LoadAllBranches();
   const unsigned int &run();
   const unsigned int &luminosityBlock();
   const unsigned long int &event();
   const float &genWeight();
   const int &nLHEPdfWeight();
   const float *LHEPdfWeight();
   const int &nLHEScaleWeight();
   const float *LHEScaleWeight();
   const int &nPSWeight();
   const float *PSWeight();
  const bool &HLT_Ele27_WPTight_Gsf();
  const bool &HLT_Ele28_WPTight_Gsf();
  const bool &HLT_Ele30_WPTight_Gsf();
  const bool &HLT_Ele32_WPTight_Gsf();
  const bool &HLT_Ele35_WPTight_Gsf();
  const bool &HLT_Ele38_WPTight_Gsf();
  const bool &HLT_Ele40_WPTight_Gsf();
  const bool &HLT_IsoMu20();
  const bool &HLT_IsoMu24();
  const bool &HLT_IsoMu24_eta2p1();
  const bool &HLT_IsoMu27();
  const bool &HLT_IsoMu30();
  const bool &HLT_Mu50();
  const bool &HLT_Mu55();
  const bool &HLT_Photon175();
  const bool &HLT_PFHT780();
  const bool &HLT_PFHT890();
  const bool &HLT_PFHT1050();
  const bool &HLT_AK8PFJet360_TrimMass30();
  const bool &HLT_AK8PFJet380_TrimMass30();
  const bool &HLT_AK8PFJet400_TrimMass30();
  const bool &HLT_AK8PFJet420_TrimMass30();
  const bool &HLT_AK8PFHT750_TrimMass50();
  const bool &HLT_AK8PFHT800_TrimMass50();
  const bool &HLT_AK8PFHT850_TrimMass50();
  const bool &HLT_AK8PFHT900_TrimMass50();
  const bool &HLT_PFJet450();
  const bool &HLT_PFJet500();
  const bool &HLT_PFJet550();
  const bool &HLT_AK8PFJet400();
  const bool &HLT_AK8PFJet450();
  const bool &HLT_AK8PFJet500();
  const bool &HLT_AK8PFJet550();
  const bool &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17();
  const bool &HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1();
  const bool &HLT_AK8PFJet330_PFAK8BTagCSV_p17();
  const bool &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02();
  const bool &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2();
  const bool &HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4();
  const bool &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20();
  const bool &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087();
  const bool &HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087();
  const bool &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20();
  const bool &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20();
  const bool &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20();
   const float &weight();
   const float &met();
   const float &metphi();
   const float &ht();
   const bool &passmetfilters();
   const float &l1PreFiringWeight();
   const float &l1PreFiringWeightUp();
   const float &l1PreFiringWeightDown();
   const float &triggerEffWeight();
   const float &triggerEff3DWeight();
   const float &triggerEffMCWeight();
   const float &triggerEffMC3DWeight();
   const float &fatJet1Pt();
   const float &fatJet1Eta();
   const float &fatJet1Phi();
   const float &fatJet1Mass();
   const float &fatJet1MassSD();
   const float &fatJet1MassSD_UnCorrected();
   const float &fatJet1MassRegressed();
   const float &fatJet1MassRegressed_UnCorrected();
   const float &fatJet1PNetXbb();
   const float &fatJet1PNetQCDb();
   const float &fatJet1PNetQCDbb();
   const float &fatJet1PNetQCDc();
   const float &fatJet1PNetQCDcc();
   const float &fatJet1PNetQCDothers();
   const float &fatJet1Tau3OverTau2();
   const int &fatJet1GenMatchIndex();
   const bool &fatJet1HasMuon();
   const bool &fatJet1HasElectron();
   const bool &fatJet1HasBJetCSVLoose();
   const bool &fatJet1HasBJetCSVMedium();
   const bool &fatJet1HasBJetCSVTight();
   const bool &fatJet1OppositeHemisphereHasBJet();
   const float &fatJet1PtOverMHH();
   const float &fatJet1PtOverMSD();
   const float &fatJet1PtOverMRegressed();
   const float &fatJet1MassSD_JMS_Down();
   const float &fatJet1MassSD_JMS_Up();
   const float &fatJet1MassSD_JMR_Down();
   const float &fatJet1MassSD_JMR_Up();
   const float &fatJet1MassRegressed_JMS_Down();
   const float &fatJet1MassRegressed_JMS_Up();
   const float &fatJet1MassRegressed_JMR_Down();
   const float &fatJet1MassRegressed_JMR_Up();
   const float &fatJet1PtOverMHH_JMS_Down();
   const float &fatJet1PtOverMHH_JMS_Up();
   const float &fatJet1PtOverMHH_JMR_Down();
   const float &fatJet1PtOverMHH_JMR_Up();
   const float &fatJet1Pt_JERUp();
   const float &fatJet1PtOverMHH_JERUp();
   const float &fatJet1Pt_JERDown();
   const float &fatJet1PtOverMHH_JERDown();
   const float &fatJet1Pt_JESUp();
   const float &fatJet1PtOverMHH_JESUp();
   const float &fatJet1Pt_JESDown();
   const float &fatJet1PtOverMHH_JESDown();
   const float &fatJet2Pt();
   const float &fatJet2Eta();
   const float &fatJet2Phi();
   const float &fatJet2Mass();
   const float &fatJet2MassSD();
   const float &fatJet2MassSD_UnCorrected();
   const float &fatJet2MassRegressed();
   const float &fatJet2MassRegressed_UnCorrected();
   const float &fatJet2PNetXbb();
   const float &fatJet2PNetQCDb();
   const float &fatJet2PNetQCDbb();
   const float &fatJet2PNetQCDc();
   const float &fatJet2PNetQCDcc();
   const float &fatJet2PNetQCDothers();
   const float &fatJet2Tau3OverTau2();
   const int &fatJet2GenMatchIndex();
   const bool &fatJet2HasMuon();
   const bool &fatJet2HasElectron();
   const bool &fatJet2HasBJetCSVLoose();
   const bool &fatJet2HasBJetCSVMedium();
   const bool &fatJet2HasBJetCSVTight();
   const bool &fatJet2OppositeHemisphereHasBJet();
   const float &fatJet2PtOverMHH();
   const float &fatJet2PtOverMSD();
   const float &fatJet2PtOverMRegressed();
   const float &fatJet2MassSD_JMS_Down();
   const float &fatJet2MassSD_JMS_Up();
   const float &fatJet2MassSD_JMR_Down();
   const float &fatJet2MassSD_JMR_Up();
   const float &fatJet2MassRegressed_JMS_Down();
   const float &fatJet2MassRegressed_JMS_Up();
   const float &fatJet2MassRegressed_JMR_Down();
   const float &fatJet2MassRegressed_JMR_Up();
   const float &fatJet2PtOverMHH_JMS_Down();
   const float &fatJet2PtOverMHH_JMS_Up();
   const float &fatJet2PtOverMHH_JMR_Down();
   const float &fatJet2PtOverMHH_JMR_Up();
   const float &fatJet2Pt_JERUp();
   const float &fatJet2PtOverMHH_JERUp();
   const float &fatJet2Pt_JERDown();
   const float &fatJet2PtOverMHH_JERDown();
   const float &fatJet2Pt_JESUp();
   const float &fatJet2PtOverMHH_JESUp();
   const float &fatJet2Pt_JESDown();
   const float &fatJet2PtOverMHH_JESDown();
   const float &hh_pt();
   const float &hh_eta();
   const float &hh_phi();
   const float &hh_mass();
   const float &hh_pt_JMS_Down();
   const float &hh_pt_JMS_Up();
   const float &hh_eta_JMS_Down();
   const float &hh_eta_JMS_Up();
   const float &hh_mass_JMS_Down();
   const float &hh_mass_JMS_Up();
   const float &hh_pt_JMR_Down();
   const float &hh_pt_JMR_Up();
   const float &hh_eta_JMR_Down();
   const float &hh_eta_JMR_Up();
   const float &hh_mass_JMR_Down();
   const float &hh_mass_JMR_Up();
   const float &hh_pt_JERUp();
   const float &hh_eta_JERUp();
   const float &hh_mass_JERUp();
   const float &hh_pt_JERDown();
   const float &hh_eta_JERDown();
   const float &hh_mass_JERDown();
   const float &hh_pt_JESUp();
   const float &hh_eta_JESUp();
   const float &hh_mass_JESUp();
   const float &hh_pt_JESDown();
   const float &hh_eta_JESDown();
   const float &hh_mass_JESDown();
   const float &deltaEta_j1j2();
   const float &deltaPhi_j1j2();
   const float &deltaR_j1j2();
   const float &ptj2_over_ptj1();
   const float &mj2_over_mj1();
   const int &isVBFtag();
   const float &dijetmass();
   const float &vbfjet1Pt();
   const float &vbfjet1Eta();
   const float &vbfjet1Phi();
   const float &vbfjet1Mass();
   const float &vbffatJet1Pt();
   const float &vbffatJet1Eta();
   const float &vbffatJet1Phi();
   const float &vbffatJet1PNetXbb();
   const float &vbfjet2Pt();
   const float &vbfjet2Eta();
   const float &vbfjet2Phi();
   const float &vbfjet2Mass();
   const float &vbffatJet2Pt();
   const float &vbffatJet2Eta();
   const float &vbffatJet2Phi();
   const float &vbffatJet2PNetXbb();
   const float &jet1Pt();
   const float &jet1Eta();
   const float &jet1Phi();
   const float &jet2Pt();
   const float &jet2Eta();
   const float &jet2Phi();
   const float &jet3Pt();
   const float &jet3Eta();
   const float &jet3Phi();
   const float &jet4Pt();
   const float &jet4Eta();
   const float &jet4Phi();
   const float &lep1Pt();
   const float &lep1Eta();
   const float &lep1Phi();
   const int &lep1Id();
   const float &lep2Pt();
   const float &lep2Eta();
   const float &lep2Phi();
   const int &lep2Id();
   const float &genHiggs1Pt();
   const float &genHiggs1Eta();
   const float &genHiggs1Phi();
   const float &genHiggs2Pt();
   const float &genHiggs2Eta();
   const float &genHiggs2Phi();
   const float &puWeight();
   const float &puWeightUp();
   const float &puWeightDown();
   const float &xsecWeight();
   const float *LHEScaleWeightNorm();
   const float *LHEPdfWeightNorm();
   const float &disc_qcd_and_ttbar_Run2_enhanced_v8p2();
   const float &disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp();
   const float &disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown();
   const float &disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp();
   const float &disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown();
   const float &disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp();
   const float &disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown();
};   

#ifndef __CINT__   
extern hhtree hh;   
#endif   

#endif   
