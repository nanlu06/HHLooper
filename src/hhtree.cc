#include "hhtree.hh"
hhtree hh; 

void hhtree::Init(TTree *tree) 
{

 run_branch = tree->GetBranch("run");
 if(run_branch) run_branch->SetAddress(&run_);
 luminosityBlock_branch = tree->GetBranch("luminosityBlock");
 if(luminosityBlock_branch) luminosityBlock_branch->SetAddress(&luminosityBlock_);
 event_branch = tree->GetBranch("event");
 if(event_branch) event_branch->SetAddress(&event_);
 genWeight_branch = tree->GetBranch("genWeight");
 if(genWeight_branch) genWeight_branch->SetAddress(&genWeight_);
 nLHEPdfWeight_branch = tree->GetBranch("nLHEPdfWeight");
 if(nLHEPdfWeight_branch) nLHEPdfWeight_branch->SetAddress(&nLHEPdfWeight_);
 LHEPdfWeight_branch = tree->GetBranch("LHEPdfWeight");
 if(LHEPdfWeight_branch) LHEPdfWeight_branch->SetAddress(&LHEPdfWeight_);
 nLHEScaleWeight_branch = tree->GetBranch("nLHEScaleWeight");
 if(nLHEScaleWeight_branch) nLHEScaleWeight_branch->SetAddress(&nLHEScaleWeight_);
 LHEScaleWeight_branch = tree->GetBranch("LHEScaleWeight");
 if(LHEScaleWeight_branch) LHEScaleWeight_branch->SetAddress(LHEScaleWeight_);
 nPSWeight_branch = tree->GetBranch("nPSWeight");
 if(nPSWeight_branch) nPSWeight_branch->SetAddress(&nPSWeight_);
 PSWeight_branch = tree->GetBranch("PSWeight");
 if(PSWeight_branch) PSWeight_branch->SetAddress(PSWeight_);
 HLT_Ele27_WPTight_Gsf_branch = tree->GetBranch("HLT_Ele27_WPTight_Gsf");
 if(HLT_Ele27_WPTight_Gsf_branch) HLT_Ele27_WPTight_Gsf_branch->SetAddress(&HLT_Ele27_WPTight_Gsf_);
 HLT_Ele28_WPTight_Gsf_branch = tree->GetBranch("HLT_Ele28_WPTight_Gsf");
 if(HLT_Ele28_WPTight_Gsf_branch) HLT_Ele28_WPTight_Gsf_branch->SetAddress(&HLT_Ele28_WPTight_Gsf_);
 HLT_Ele30_WPTight_Gsf_branch = tree->GetBranch("HLT_Ele30_WPTight_Gsf");
 if(HLT_Ele30_WPTight_Gsf_branch) HLT_Ele30_WPTight_Gsf_branch->SetAddress(&HLT_Ele30_WPTight_Gsf_);
 HLT_Ele32_WPTight_Gsf_branch = tree->GetBranch("HLT_Ele32_WPTight_Gsf");
 if(HLT_Ele32_WPTight_Gsf_branch) HLT_Ele32_WPTight_Gsf_branch->SetAddress(&HLT_Ele32_WPTight_Gsf_);
 HLT_Ele32_WPTight_Gsf_L1DoubleEG_branch = tree->GetBranch("HLT_Ele32_WPTight_Gsf_L1DoubleEG");
 if(HLT_Ele32_WPTight_Gsf_L1DoubleEG_branch) HLT_Ele32_WPTight_Gsf_L1DoubleEG_branch->SetAddress(&HLT_Ele32_WPTight_Gsf_L1DoubleEG_);

 HLT_Ele35_WPTight_Gsf_branch = tree->GetBranch("HLT_Ele35_WPTight_Gsf");
 if(HLT_Ele35_WPTight_Gsf_branch) HLT_Ele35_WPTight_Gsf_branch->SetAddress(&HLT_Ele35_WPTight_Gsf_);
 HLT_Ele38_WPTight_Gsf_branch = tree->GetBranch("HLT_Ele38_WPTight_Gsf");
 if(HLT_Ele38_WPTight_Gsf_branch) HLT_Ele38_WPTight_Gsf_branch->SetAddress(&HLT_Ele38_WPTight_Gsf_);
 HLT_Ele40_WPTight_Gsf_branch = tree->GetBranch("HLT_Ele40_WPTight_Gsf");
 if(HLT_Ele40_WPTight_Gsf_branch) HLT_Ele40_WPTight_Gsf_branch->SetAddress(&HLT_Ele40_WPTight_Gsf_);
 HLT_IsoMu20_branch = tree->GetBranch("HLT_IsoMu20");
 if(HLT_IsoMu20_branch) HLT_IsoMu20_branch->SetAddress(&HLT_IsoMu20_);
 HLT_IsoMu24_branch = tree->GetBranch("HLT_IsoMu24");
 if(HLT_IsoMu24_branch) HLT_IsoMu24_branch->SetAddress(&HLT_IsoMu24_);
 HLT_IsoMu24_eta2p1_branch = tree->GetBranch("HLT_IsoMu24_eta2p1");
 if(HLT_IsoMu24_eta2p1_branch) HLT_IsoMu24_eta2p1_branch->SetAddress(&HLT_IsoMu24_eta2p1_);
 HLT_IsoMu27_branch = tree->GetBranch("HLT_IsoMu27");
 if(HLT_IsoMu27_branch) HLT_IsoMu27_branch->SetAddress(&HLT_IsoMu27_);
 HLT_IsoMu30_branch = tree->GetBranch("HLT_IsoMu30");
 if(HLT_IsoMu30_branch) HLT_IsoMu30_branch->SetAddress(&HLT_IsoMu30_);
 HLT_Mu50_branch = tree->GetBranch("HLT_Mu50");
 if(HLT_Mu50_branch) HLT_Mu50_branch->SetAddress(&HLT_Mu50_);
 HLT_Mu55_branch = tree->GetBranch("HLT_Mu55");
 if(HLT_Mu55_branch) HLT_Mu55_branch->SetAddress(&HLT_Mu55_);
 HLT_Photon175_branch = tree->GetBranch("HLT_Photon175");
 if(HLT_Photon175_branch) HLT_Photon175_branch->SetAddress(&HLT_Photon175_);
 HLT_PFHT780_branch = tree->GetBranch("HLT_PFHT780");
 if(HLT_PFHT780_branch) HLT_PFHT780_branch->SetAddress(&HLT_PFHT780_);
 HLT_PFHT890_branch = tree->GetBranch("HLT_PFHT890");
 if(HLT_PFHT890_branch) HLT_PFHT890_branch->SetAddress(&HLT_PFHT890_);
 HLT_PFHT1050_branch = tree->GetBranch("HLT_PFHT1050");
 if(HLT_PFHT1050_branch) HLT_PFHT1050_branch->SetAddress(&HLT_PFHT1050_);
 HLT_AK8PFJet360_TrimMass30_branch = tree->GetBranch("HLT_AK8PFJet360_TrimMass30");
 if(HLT_AK8PFJet360_TrimMass30_branch) HLT_AK8PFJet360_TrimMass30_branch->SetAddress(&HLT_AK8PFJet360_TrimMass30_);
 HLT_AK8PFJet380_TrimMass30_branch = tree->GetBranch("HLT_AK8PFJet380_TrimMass30");
 if(HLT_AK8PFJet380_TrimMass30_branch) HLT_AK8PFJet380_TrimMass30_branch->SetAddress(&HLT_AK8PFJet380_TrimMass30_);
 HLT_AK8PFJet400_TrimMass30_branch = tree->GetBranch("HLT_AK8PFJet400_TrimMass30");
 if(HLT_AK8PFJet400_TrimMass30_branch) HLT_AK8PFJet400_TrimMass30_branch->SetAddress(&HLT_AK8PFJet400_TrimMass30_);
 HLT_AK8PFJet420_TrimMass30_branch = tree->GetBranch("HLT_AK8PFJet420_TrimMass30");
 if(HLT_AK8PFJet420_TrimMass30_branch) HLT_AK8PFJet420_TrimMass30_branch->SetAddress(&HLT_AK8PFJet420_TrimMass30_);
 HLT_AK8PFHT750_TrimMass50_branch = tree->GetBranch("HLT_AK8PFHT750_TrimMass50");
 if(HLT_AK8PFHT750_TrimMass50_branch) HLT_AK8PFHT750_TrimMass50_branch->SetAddress(&HLT_AK8PFHT750_TrimMass50_);
 HLT_AK8PFHT800_TrimMass50_branch = tree->GetBranch("HLT_AK8PFHT800_TrimMass50");
 if(HLT_AK8PFHT800_TrimMass50_branch) HLT_AK8PFHT800_TrimMass50_branch->SetAddress(&HLT_AK8PFHT800_TrimMass50_);
 HLT_AK8PFHT850_TrimMass50_branch = tree->GetBranch("HLT_AK8PFHT850_TrimMass50");
 if(HLT_AK8PFHT850_TrimMass50_branch) HLT_AK8PFHT850_TrimMass50_branch->SetAddress(&HLT_AK8PFHT850_TrimMass50_);
 HLT_AK8PFHT900_TrimMass50_branch = tree->GetBranch("HLT_AK8PFHT900_TrimMass50");
 if(HLT_AK8PFHT900_TrimMass50_branch) HLT_AK8PFHT900_TrimMass50_branch->SetAddress(&HLT_AK8PFHT900_TrimMass50_);
 HLT_PFJet450_branch = tree->GetBranch("HLT_PFJet450");
 if(HLT_PFJet450_branch) HLT_PFJet450_branch->SetAddress(&HLT_PFJet450_);
 HLT_PFJet500_branch = tree->GetBranch("HLT_PFJet500");
 if(HLT_PFJet500_branch) HLT_PFJet500_branch->SetAddress(&HLT_PFJet500_);
 HLT_PFJet550_branch = tree->GetBranch("HLT_PFJet550");
 if(HLT_PFJet550_branch) HLT_PFJet550_branch->SetAddress(&HLT_PFJet550_);
 HLT_AK8PFJet400_branch = tree->GetBranch("HLT_AK8PFJet400");
 if(HLT_AK8PFJet400_branch) HLT_AK8PFJet400_branch->SetAddress(&HLT_AK8PFJet400_);
 HLT_AK8PFJet450_branch = tree->GetBranch("HLT_AK8PFJet450");
 if(HLT_AK8PFJet450_branch) HLT_AK8PFJet450_branch->SetAddress(&HLT_AK8PFJet450_);
 HLT_AK8PFJet500_branch = tree->GetBranch("HLT_AK8PFJet500");
 if(HLT_AK8PFJet500_branch) HLT_AK8PFJet500_branch->SetAddress(&HLT_AK8PFJet500_);
 HLT_AK8PFJet550_branch = tree->GetBranch("HLT_AK8PFJet550");
 if(HLT_AK8PFJet550_branch) HLT_AK8PFJet550_branch->SetAddress(&HLT_AK8PFJet550_);
 HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_branch = tree->GetBranch("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17");
 if(HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_branch) HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_branch->SetAddress(&HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_);
 HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_branch = tree->GetBranch("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1");
 if(HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_branch) HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_branch->SetAddress(&HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_);
 HLT_AK8PFJet330_PFAK8BTagCSV_p17_branch = tree->GetBranch("HLT_AK8PFJet330_PFAK8BTagCSV_p17");
 if(HLT_AK8PFJet330_PFAK8BTagCSV_p17_branch) HLT_AK8PFJet330_PFAK8BTagCSV_p17_branch->SetAddress(&HLT_AK8PFJet330_PFAK8BTagCSV_p17_);
 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_branch = tree->GetBranch("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02");
 if(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_branch) HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_branch->SetAddress(&HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_);
 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_branch = tree->GetBranch("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2");
 if(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_branch) HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_branch->SetAddress(&HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_);
 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_branch = tree->GetBranch("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4");
 if(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_branch) HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_branch->SetAddress(&HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_);
 HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_branch = tree->GetBranch("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20");
 if(HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_branch) HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_branch->SetAddress(&HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_);
 HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_branch = tree->GetBranch("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087");
 if(HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_branch) HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_branch->SetAddress(&HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_);
 HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_branch = tree->GetBranch("HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087");
 if(HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_branch) HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_branch->SetAddress(&HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_);
 HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_branch = tree->GetBranch("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20");
 if(HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_branch) HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_branch->SetAddress(&HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_);
 HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_branch = tree->GetBranch("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20");
 if(HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_branch) HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_branch->SetAddress(&HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_);
 HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_branch = tree->GetBranch("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20");
 if(HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_branch) HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_branch->SetAddress(&HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_);
 weight_branch = tree->GetBranch("weight");
 if(weight_branch) weight_branch->SetAddress(&weight_);
 met_branch = tree->GetBranch("met");
 if(met_branch) met_branch->SetAddress(&met_);
 metphi_branch = tree->GetBranch("metphi");
 if(metphi_branch) metphi_branch->SetAddress(&metphi_);
 ht_branch = tree->GetBranch("ht");
 if(ht_branch) ht_branch->SetAddress(&ht_);
 passmetfilters_branch = tree->GetBranch("passmetfilters");
 if(passmetfilters_branch) passmetfilters_branch->SetAddress(&passmetfilters_);
 l1PreFiringWeight_branch = tree->GetBranch("l1PreFiringWeight");
 if(l1PreFiringWeight_branch) l1PreFiringWeight_branch->SetAddress(&l1PreFiringWeight_);
 l1PreFiringWeightUp_branch = tree->GetBranch("l1PreFiringWeightUp");
 if(l1PreFiringWeightUp_branch) l1PreFiringWeightUp_branch->SetAddress(&l1PreFiringWeightUp_);
 l1PreFiringWeightDown_branch = tree->GetBranch("l1PreFiringWeightDown");
 if(l1PreFiringWeightDown_branch) l1PreFiringWeightDown_branch->SetAddress(&l1PreFiringWeightDown_);
 triggerEffWeight_branch = tree->GetBranch("triggerEffWeight");
 if(triggerEffWeight_branch) triggerEffWeight_branch->SetAddress(&triggerEffWeight_);
 triggerEff3DWeight_branch = tree->GetBranch("triggerEff3DWeight");
 if(triggerEff3DWeight_branch) triggerEff3DWeight_branch->SetAddress(&triggerEff3DWeight_);
 triggerEffMCWeight_branch = tree->GetBranch("triggerEffMCWeight");
 if(triggerEffMCWeight_branch) triggerEffMCWeight_branch->SetAddress(&triggerEffMCWeight_);
 triggerEffMC3DWeight_branch = tree->GetBranch("triggerEffMC3DWeight");
 if(triggerEffMC3DWeight_branch) triggerEffMC3DWeight_branch->SetAddress(&triggerEffMC3DWeight_);
 fatJet1Pt_branch = tree->GetBranch("fatJet1Pt");
 if(fatJet1Pt_branch) fatJet1Pt_branch->SetAddress(&fatJet1Pt_);
 fatJet1Eta_branch = tree->GetBranch("fatJet1Eta");
 if(fatJet1Eta_branch) fatJet1Eta_branch->SetAddress(&fatJet1Eta_);
 fatJet1Phi_branch = tree->GetBranch("fatJet1Phi");
 if(fatJet1Phi_branch) fatJet1Phi_branch->SetAddress(&fatJet1Phi_);
 fatJet1Mass_branch = tree->GetBranch("fatJet1Mass");
 if(fatJet1Mass_branch) fatJet1Mass_branch->SetAddress(&fatJet1Mass_);
 fatJet1MassSD_branch = tree->GetBranch("fatJet1MassSD");
 if(fatJet1MassSD_branch) fatJet1MassSD_branch->SetAddress(&fatJet1MassSD_);
 fatJet1MassSD_UnCorrected_branch = tree->GetBranch("fatJet1MassSD_UnCorrected");
 if(fatJet1MassSD_UnCorrected_branch) fatJet1MassSD_UnCorrected_branch->SetAddress(&fatJet1MassSD_UnCorrected_);
 fatJet1MassRegressed_branch = tree->GetBranch("fatJet1MassRegressed");
 if(fatJet1MassRegressed_branch) fatJet1MassRegressed_branch->SetAddress(&fatJet1MassRegressed_);
 fatJet1MassRegressed_UnCorrected_branch = tree->GetBranch("fatJet1MassRegressed_UnCorrected");
 if(fatJet1MassRegressed_UnCorrected_branch) fatJet1MassRegressed_UnCorrected_branch->SetAddress(&fatJet1MassRegressed_UnCorrected_);
 fatJet1PNetXbb_branch = tree->GetBranch("fatJet1PNetXbb");
 if(fatJet1PNetXbb_branch) fatJet1PNetXbb_branch->SetAddress(&fatJet1PNetXbb_);
 fatJet1PNetQCDb_branch = tree->GetBranch("fatJet1PNetQCDb");
 if(fatJet1PNetQCDb_branch) fatJet1PNetQCDb_branch->SetAddress(&fatJet1PNetQCDb_);
 fatJet1PNetQCDbb_branch = tree->GetBranch("fatJet1PNetQCDbb");
 if(fatJet1PNetQCDbb_branch) fatJet1PNetQCDbb_branch->SetAddress(&fatJet1PNetQCDbb_);
 fatJet1PNetQCDc_branch = tree->GetBranch("fatJet1PNetQCDc");
 if(fatJet1PNetQCDc_branch) fatJet1PNetQCDc_branch->SetAddress(&fatJet1PNetQCDc_);
 fatJet1PNetQCDcc_branch = tree->GetBranch("fatJet1PNetQCDcc");
 if(fatJet1PNetQCDcc_branch) fatJet1PNetQCDcc_branch->SetAddress(&fatJet1PNetQCDcc_);
 fatJet1PNetQCDothers_branch = tree->GetBranch("fatJet1PNetQCDothers");
 if(fatJet1PNetQCDothers_branch) fatJet1PNetQCDothers_branch->SetAddress(&fatJet1PNetQCDothers_);
 fatJet1Tau3OverTau2_branch = tree->GetBranch("fatJet1Tau3OverTau2");
 if(fatJet1Tau3OverTau2_branch) fatJet1Tau3OverTau2_branch->SetAddress(&fatJet1Tau3OverTau2_);
 fatJet1GenMatchIndex_branch = tree->GetBranch("fatJet1GenMatchIndex");
 if(fatJet1GenMatchIndex_branch) fatJet1GenMatchIndex_branch->SetAddress(&fatJet1GenMatchIndex_);
 fatJet1HasMuon_branch = tree->GetBranch("fatJet1HasMuon");
 if(fatJet1HasMuon_branch) fatJet1HasMuon_branch->SetAddress(&fatJet1HasMuon_);
 fatJet1HasElectron_branch = tree->GetBranch("fatJet1HasElectron");
 if(fatJet1HasElectron_branch) fatJet1HasElectron_branch->SetAddress(&fatJet1HasElectron_);
 fatJet1HasBJetCSVLoose_branch = tree->GetBranch("fatJet1HasBJetCSVLoose");
 if(fatJet1HasBJetCSVLoose_branch) fatJet1HasBJetCSVLoose_branch->SetAddress(&fatJet1HasBJetCSVLoose_);
 fatJet1HasBJetCSVMedium_branch = tree->GetBranch("fatJet1HasBJetCSVMedium");
 if(fatJet1HasBJetCSVMedium_branch) fatJet1HasBJetCSVMedium_branch->SetAddress(&fatJet1HasBJetCSVMedium_);
 fatJet1HasBJetCSVTight_branch = tree->GetBranch("fatJet1HasBJetCSVTight");
 if(fatJet1HasBJetCSVTight_branch) fatJet1HasBJetCSVTight_branch->SetAddress(&fatJet1HasBJetCSVTight_);
 fatJet1OppositeHemisphereHasBJet_branch = tree->GetBranch("fatJet1OppositeHemisphereHasBJet");
 if(fatJet1OppositeHemisphereHasBJet_branch) fatJet1OppositeHemisphereHasBJet_branch->SetAddress(&fatJet1OppositeHemisphereHasBJet_);
 fatJet1PtOverMHH_branch = tree->GetBranch("fatJet1PtOverMHH");
 if(fatJet1PtOverMHH_branch) fatJet1PtOverMHH_branch->SetAddress(&fatJet1PtOverMHH_);
 fatJet1PtOverMSD_branch = tree->GetBranch("fatJet1PtOverMSD");
 if(fatJet1PtOverMSD_branch) fatJet1PtOverMSD_branch->SetAddress(&fatJet1PtOverMSD_);
 fatJet1PtOverMRegressed_branch = tree->GetBranch("fatJet1PtOverMRegressed");
 if(fatJet1PtOverMRegressed_branch) fatJet1PtOverMRegressed_branch->SetAddress(&fatJet1PtOverMRegressed_);
 fatJet1MassSD_JMS_Down_branch = tree->GetBranch("fatJet1MassSD_JMS_Down");
 if(fatJet1MassSD_JMS_Down_branch) fatJet1MassSD_JMS_Down_branch->SetAddress(&fatJet1MassSD_JMS_Down_);
 fatJet1MassSD_JMS_Up_branch = tree->GetBranch("fatJet1MassSD_JMS_Up");
 if(fatJet1MassSD_JMS_Up_branch) fatJet1MassSD_JMS_Up_branch->SetAddress(&fatJet1MassSD_JMS_Up_);
 fatJet1MassSD_JMR_Down_branch = tree->GetBranch("fatJet1MassSD_JMR_Down");
 if(fatJet1MassSD_JMR_Down_branch) fatJet1MassSD_JMR_Down_branch->SetAddress(&fatJet1MassSD_JMR_Down_);
 fatJet1MassSD_JMR_Up_branch = tree->GetBranch("fatJet1MassSD_JMR_Up");
 if(fatJet1MassSD_JMR_Up_branch) fatJet1MassSD_JMR_Up_branch->SetAddress(&fatJet1MassSD_JMR_Up_);
 fatJet1MassRegressed_JMS_Down_branch = tree->GetBranch("fatJet1MassRegressed_JMS_Down");
 if(fatJet1MassRegressed_JMS_Down_branch) fatJet1MassRegressed_JMS_Down_branch->SetAddress(&fatJet1MassRegressed_JMS_Down_);
 fatJet1MassRegressed_JMS_Up_branch = tree->GetBranch("fatJet1MassRegressed_JMS_Up");
 if(fatJet1MassRegressed_JMS_Up_branch) fatJet1MassRegressed_JMS_Up_branch->SetAddress(&fatJet1MassRegressed_JMS_Up_);
 fatJet1MassRegressed_JMR_Down_branch = tree->GetBranch("fatJet1MassRegressed_JMR_Down");
 if(fatJet1MassRegressed_JMR_Down_branch) fatJet1MassRegressed_JMR_Down_branch->SetAddress(&fatJet1MassRegressed_JMR_Down_);
 fatJet1MassRegressed_JMR_Up_branch = tree->GetBranch("fatJet1MassRegressed_JMR_Up");
 if(fatJet1MassRegressed_JMR_Up_branch) fatJet1MassRegressed_JMR_Up_branch->SetAddress(&fatJet1MassRegressed_JMR_Up_);
 fatJet1PtOverMHH_JMS_Down_branch = tree->GetBranch("fatJet1PtOverMHH_JMS_Down");
 if(fatJet1PtOverMHH_JMS_Down_branch) fatJet1PtOverMHH_JMS_Down_branch->SetAddress(&fatJet1PtOverMHH_JMS_Down_);
 fatJet1PtOverMHH_JMS_Up_branch = tree->GetBranch("fatJet1PtOverMHH_JMS_Up");
 if(fatJet1PtOverMHH_JMS_Up_branch) fatJet1PtOverMHH_JMS_Up_branch->SetAddress(&fatJet1PtOverMHH_JMS_Up_);
 fatJet1PtOverMHH_JMR_Down_branch = tree->GetBranch("fatJet1PtOverMHH_JMR_Down");
 if(fatJet1PtOverMHH_JMR_Down_branch) fatJet1PtOverMHH_JMR_Down_branch->SetAddress(&fatJet1PtOverMHH_JMR_Down_);
 fatJet1PtOverMHH_JMR_Up_branch = tree->GetBranch("fatJet1PtOverMHH_JMR_Up");
 if(fatJet1PtOverMHH_JMR_Up_branch) fatJet1PtOverMHH_JMR_Up_branch->SetAddress(&fatJet1PtOverMHH_JMR_Up_);
 fatJet1Pt_JERUp_branch = tree->GetBranch("fatJet1Pt_JERUp");
 if(fatJet1Pt_JERUp_branch) fatJet1Pt_JERUp_branch->SetAddress(&fatJet1Pt_JERUp_);
 fatJet1PtOverMHH_JERUp_branch = tree->GetBranch("fatJet1PtOverMHH_JERUp");
 if(fatJet1PtOverMHH_JERUp_branch) fatJet1PtOverMHH_JERUp_branch->SetAddress(&fatJet1PtOverMHH_JERUp_);
 fatJet1Pt_JERDown_branch = tree->GetBranch("fatJet1Pt_JERDown");
 if(fatJet1Pt_JERDown_branch) fatJet1Pt_JERDown_branch->SetAddress(&fatJet1Pt_JERDown_);
 fatJet1PtOverMHH_JERDown_branch = tree->GetBranch("fatJet1PtOverMHH_JERDown");
 if(fatJet1PtOverMHH_JERDown_branch) fatJet1PtOverMHH_JERDown_branch->SetAddress(&fatJet1PtOverMHH_JERDown_);
 fatJet1Pt_JESUp_branch = tree->GetBranch("fatJet1Pt_JESUp");
 if(fatJet1Pt_JESUp_branch) fatJet1Pt_JESUp_branch->SetAddress(&fatJet1Pt_JESUp_);
 fatJet1PtOverMHH_JESUp_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp");
 if(fatJet1PtOverMHH_JESUp_branch) fatJet1PtOverMHH_JESUp_branch->SetAddress(&fatJet1PtOverMHH_JESUp_);
 fatJet1Pt_JESDown_branch = tree->GetBranch("fatJet1Pt_JESDown");
 if(fatJet1Pt_JESDown_branch) fatJet1Pt_JESDown_branch->SetAddress(&fatJet1Pt_JESDown_);
 fatJet1PtOverMHH_JESDown_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown");
 if(fatJet1PtOverMHH_JESDown_branch) fatJet1PtOverMHH_JESDown_branch->SetAddress(&fatJet1PtOverMHH_JESDown_);
  fatJet1Pt_JESUp_Abs_branch = tree->GetBranch("fatJet1Pt_JESUp_Abs");
  if(fatJet1Pt_JESUp_Abs_branch) fatJet1Pt_JESUp_Abs_branch->SetAddress(&fatJet1Pt_JESUp_Abs_);
  fatJet1PtOverMHH_JESUp_Abs_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_Abs");
  if(fatJet1PtOverMHH_JESUp_Abs_branch) fatJet1PtOverMHH_JESUp_Abs_branch->SetAddress(&fatJet1PtOverMHH_JESUp_Abs_);
  fatJet1Pt_JESDown_Abs_branch = tree->GetBranch("fatJet1Pt_JESDown_Abs");
  if(fatJet1Pt_JESDown_Abs_branch) fatJet1Pt_JESDown_Abs_branch->SetAddress(&fatJet1Pt_JESDown_Abs_);
  fatJet1PtOverMHH_JESDown_Abs_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_Abs");
  if(fatJet1PtOverMHH_JESDown_Abs_branch) fatJet1PtOverMHH_JESDown_Abs_branch->SetAddress(&fatJet1PtOverMHH_JESDown_Abs_);
  fatJet1Pt_JESUp_Abs_2016_branch = tree->GetBranch("fatJet1Pt_JESUp_Abs_2016");
  if(fatJet1Pt_JESUp_Abs_2016_branch) fatJet1Pt_JESUp_Abs_2016_branch->SetAddress(&fatJet1Pt_JESUp_Abs_2016_);
  fatJet1PtOverMHH_JESUp_Abs_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_Abs_2016");
  if(fatJet1PtOverMHH_JESUp_Abs_2016_branch) fatJet1PtOverMHH_JESUp_Abs_2016_branch->SetAddress(&fatJet1PtOverMHH_JESUp_Abs_2016_);
  fatJet1Pt_JESDowb_Abs_2016_branch = tree->GetBranch("fatJet1Pt_JESDowb_Abs_2016");
  if(fatJet1Pt_JESDowb_Abs_2016_branch) fatJet1Pt_JESDowb_Abs_2016_branch->SetAddress(&fatJet1Pt_JESDowb_Abs_2016_);
  fatJet1PtOverMHH_JESDowb_Abs_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESDowb_Abs_2016");
  if(fatJet1PtOverMHH_JESDowb_Abs_2016_branch) fatJet1PtOverMHH_JESDowb_Abs_2016_branch->SetAddress(&fatJet1PtOverMHH_JESDowb_Abs_2016_);
  fatJet1Pt_JESUp_Abs_2017_branch = tree->GetBranch("fatJet1Pt_JESUp_Abs_2017");
  if(fatJet1Pt_JESUp_Abs_2017_branch) fatJet1Pt_JESUp_Abs_2017_branch->SetAddress(&fatJet1Pt_JESUp_Abs_2017_);
  fatJet1PtOverMHH_JESUp_Abs_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_Abs_2017");
  if(fatJet1PtOverMHH_JESUp_Abs_2017_branch) fatJet1PtOverMHH_JESUp_Abs_2017_branch->SetAddress(&fatJet1PtOverMHH_JESUp_Abs_2017_);
  fatJet1Pt_JESDowb_Abs_2017_branch = tree->GetBranch("fatJet1Pt_JESDowb_Abs_2017");
  if(fatJet1Pt_JESDowb_Abs_2017_branch) fatJet1Pt_JESDowb_Abs_2017_branch->SetAddress(&fatJet1Pt_JESDowb_Abs_2017_);
  fatJet1PtOverMHH_JESDowb_Abs_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESDowb_Abs_2017");
  if(fatJet1PtOverMHH_JESDowb_Abs_2017_branch) fatJet1PtOverMHH_JESDowb_Abs_2017_branch->SetAddress(&fatJet1PtOverMHH_JESDowb_Abs_2017_);
  fatJet1Pt_JESUp_Abs_2018_branch = tree->GetBranch("fatJet1Pt_JESUp_Abs_2018");
  if(fatJet1Pt_JESUp_Abs_2018_branch) fatJet1Pt_JESUp_Abs_2018_branch->SetAddress(&fatJet1Pt_JESUp_Abs_2018_);
  fatJet1PtOverMHH_JESUp_Abs_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_Abs_2018");
  if(fatJet1PtOverMHH_JESUp_Abs_2018_branch) fatJet1PtOverMHH_JESUp_Abs_2018_branch->SetAddress(&fatJet1PtOverMHH_JESUp_Abs_2018_);
  fatJet1Pt_JESDowb_Abs_2018_branch = tree->GetBranch("fatJet1Pt_JESDowb_Abs_2018");
  if(fatJet1Pt_JESDowb_Abs_2018_branch) fatJet1Pt_JESDowb_Abs_2018_branch->SetAddress(&fatJet1Pt_JESDowb_Abs_2018_);
  fatJet1PtOverMHH_JESDowb_Abs_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESDowb_Abs_2018");
  if(fatJet1PtOverMHH_JESDowb_Abs_2018_branch) fatJet1PtOverMHH_JESDowb_Abs_2018_branch->SetAddress(&fatJet1PtOverMHH_JESDowb_Abs_2018_);
  fatJet1Pt_JESUp_BBEC1_branch = tree->GetBranch("fatJet1Pt_JESUp_BBEC1");
  if(fatJet1Pt_JESUp_BBEC1_branch) fatJet1Pt_JESUp_BBEC1_branch->SetAddress(&fatJet1Pt_JESUp_BBEC1_);
  fatJet1PtOverMHH_JESUp_BBEC1_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_BBEC1");
  if(fatJet1PtOverMHH_JESUp_BBEC1_branch) fatJet1PtOverMHH_JESUp_BBEC1_branch->SetAddress(&fatJet1PtOverMHH_JESUp_BBEC1_);
  fatJet1Pt_JESDown_BBEC1_branch = tree->GetBranch("fatJet1Pt_JESDown_BBEC1");
  if(fatJet1Pt_JESDown_BBEC1_branch) fatJet1Pt_JESDown_BBEC1_branch->SetAddress(&fatJet1Pt_JESDown_BBEC1_);
  fatJet1PtOverMHH_JESDown_BBEC1_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_BBEC1");
  if(fatJet1PtOverMHH_JESDown_BBEC1_branch) fatJet1PtOverMHH_JESDown_BBEC1_branch->SetAddress(&fatJet1PtOverMHH_JESDown_BBEC1_);
  fatJet1Pt_JESUp_BBEC1_2016_branch = tree->GetBranch("fatJet1Pt_JESUp_BBEC1_2016");
  if(fatJet1Pt_JESUp_BBEC1_2016_branch) fatJet1Pt_JESUp_BBEC1_2016_branch->SetAddress(&fatJet1Pt_JESUp_BBEC1_2016_);
  fatJet1PtOverMHH_JESUp_BBEC1_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_BBEC1_2016");
  if(fatJet1PtOverMHH_JESUp_BBEC1_2016_branch) fatJet1PtOverMHH_JESUp_BBEC1_2016_branch->SetAddress(&fatJet1PtOverMHH_JESUp_BBEC1_2016_);
  fatJet1Pt_JESDown_BBEC1_2016_branch = tree->GetBranch("fatJet1Pt_JESDown_BBEC1_2016");
  if(fatJet1Pt_JESDown_BBEC1_2016_branch) fatJet1Pt_JESDown_BBEC1_2016_branch->SetAddress(&fatJet1Pt_JESDown_BBEC1_2016_);
  fatJet1PtOverMHH_JESDown_BBEC1_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_BBEC1_2016");
  if(fatJet1PtOverMHH_JESDown_BBEC1_2016_branch) fatJet1PtOverMHH_JESDown_BBEC1_2016_branch->SetAddress(&fatJet1PtOverMHH_JESDown_BBEC1_2016_);
  fatJet1Pt_JESUp_BBEC1_2017_branch = tree->GetBranch("fatJet1Pt_JESUp_BBEC1_2017");
  if(fatJet1Pt_JESUp_BBEC1_2017_branch) fatJet1Pt_JESUp_BBEC1_2017_branch->SetAddress(&fatJet1Pt_JESUp_BBEC1_2017_);
  fatJet1PtOverMHH_JESUp_BBEC1_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_BBEC1_2017");
  if(fatJet1PtOverMHH_JESUp_BBEC1_2017_branch) fatJet1PtOverMHH_JESUp_BBEC1_2017_branch->SetAddress(&fatJet1PtOverMHH_JESUp_BBEC1_2017_);
  fatJet1Pt_JESDown_BBEC1_2017_branch = tree->GetBranch("fatJet1Pt_JESDown_BBEC1_2017");
  if(fatJet1Pt_JESDown_BBEC1_2017_branch) fatJet1Pt_JESDown_BBEC1_2017_branch->SetAddress(&fatJet1Pt_JESDown_BBEC1_2017_);
  fatJet1PtOverMHH_JESDown_BBEC1_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_BBEC1_2017");
  if(fatJet1PtOverMHH_JESDown_BBEC1_2017_branch) fatJet1PtOverMHH_JESDown_BBEC1_2017_branch->SetAddress(&fatJet1PtOverMHH_JESDown_BBEC1_2017_);
  fatJet1Pt_JESUp_BBEC1_2018_branch = tree->GetBranch("fatJet1Pt_JESUp_BBEC1_2018");
  if(fatJet1Pt_JESUp_BBEC1_2018_branch) fatJet1Pt_JESUp_BBEC1_2018_branch->SetAddress(&fatJet1Pt_JESUp_BBEC1_2018_);
  fatJet1PtOverMHH_JESUp_BBEC1_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_BBEC1_2018");
  if(fatJet1PtOverMHH_JESUp_BBEC1_2018_branch) fatJet1PtOverMHH_JESUp_BBEC1_2018_branch->SetAddress(&fatJet1PtOverMHH_JESUp_BBEC1_2018_);
  fatJet1Pt_JESDown_BBEC1_2018_branch = tree->GetBranch("fatJet1Pt_JESDown_BBEC1_2018");
  if(fatJet1Pt_JESDown_BBEC1_2018_branch) fatJet1Pt_JESDown_BBEC1_2018_branch->SetAddress(&fatJet1Pt_JESDown_BBEC1_2018_);
  fatJet1PtOverMHH_JESDown_BBEC1_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_BBEC1_2018");
  if(fatJet1PtOverMHH_JESDown_BBEC1_2018_branch) fatJet1PtOverMHH_JESDown_BBEC1_2018_branch->SetAddress(&fatJet1PtOverMHH_JESDown_BBEC1_2018_);
  fatJet1Pt_JESUp_EC2_branch = tree->GetBranch("fatJet1Pt_JESUp_EC2");
  if(fatJet1Pt_JESUp_EC2_branch) fatJet1Pt_JESUp_EC2_branch->SetAddress(&fatJet1Pt_JESUp_EC2_);
  fatJet1PtOverMHH_JESUp_EC2_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_EC2");
  if(fatJet1PtOverMHH_JESUp_EC2_branch) fatJet1PtOverMHH_JESUp_EC2_branch->SetAddress(&fatJet1PtOverMHH_JESUp_EC2_);
  fatJet1Pt_JESDown_EC2_branch = tree->GetBranch("fatJet1Pt_JESDown_EC2");
  if(fatJet1Pt_JESDown_EC2_branch) fatJet1Pt_JESDown_EC2_branch->SetAddress(&fatJet1Pt_JESDown_EC2_);
  fatJet1PtOverMHH_JESDown_EC2_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_EC2");
  if(fatJet1PtOverMHH_JESDown_EC2_branch) fatJet1PtOverMHH_JESDown_EC2_branch->SetAddress(&fatJet1PtOverMHH_JESDown_EC2_);
  fatJet1Pt_JESUp_EC2_2016_branch = tree->GetBranch("fatJet1Pt_JESUp_EC2_2016");
  if(fatJet1Pt_JESUp_EC2_2016_branch) fatJet1Pt_JESUp_EC2_2016_branch->SetAddress(&fatJet1Pt_JESUp_EC2_2016_);
  fatJet1PtOverMHH_JESUp_EC2_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_EC2_2016");
  if(fatJet1PtOverMHH_JESUp_EC2_2016_branch) fatJet1PtOverMHH_JESUp_EC2_2016_branch->SetAddress(&fatJet1PtOverMHH_JESUp_EC2_2016_);
  fatJet1Pt_JESDown_EC2_2016_branch = tree->GetBranch("fatJet1Pt_JESDown_EC2_2016");
  if(fatJet1Pt_JESDown_EC2_2016_branch) fatJet1Pt_JESDown_EC2_2016_branch->SetAddress(&fatJet1Pt_JESDown_EC2_2016_);
  fatJet1PtOverMHH_JESDown_EC2_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_EC2_2016");
  if(fatJet1PtOverMHH_JESDown_EC2_2016_branch) fatJet1PtOverMHH_JESDown_EC2_2016_branch->SetAddress(&fatJet1PtOverMHH_JESDown_EC2_2016_);
  fatJet1Pt_JESUp_EC2_2017_branch = tree->GetBranch("fatJet1Pt_JESUp_EC2_2017");
  if(fatJet1Pt_JESUp_EC2_2017_branch) fatJet1Pt_JESUp_EC2_2017_branch->SetAddress(&fatJet1Pt_JESUp_EC2_2017_);
  fatJet1PtOverMHH_JESUp_EC2_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_EC2_2017");
  if(fatJet1PtOverMHH_JESUp_EC2_2017_branch) fatJet1PtOverMHH_JESUp_EC2_2017_branch->SetAddress(&fatJet1PtOverMHH_JESUp_EC2_2017_);
  fatJet1Pt_JESDown_EC2_2017_branch = tree->GetBranch("fatJet1Pt_JESDown_EC2_2017");
  if(fatJet1Pt_JESDown_EC2_2017_branch) fatJet1Pt_JESDown_EC2_2017_branch->SetAddress(&fatJet1Pt_JESDown_EC2_2017_);
  fatJet1PtOverMHH_JESDown_EC2_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_EC2_2017");
  if(fatJet1PtOverMHH_JESDown_EC2_2017_branch) fatJet1PtOverMHH_JESDown_EC2_2017_branch->SetAddress(&fatJet1PtOverMHH_JESDown_EC2_2017_);
  fatJet1Pt_JESUp_EC2_2018_branch = tree->GetBranch("fatJet1Pt_JESUp_EC2_2018");
  if(fatJet1Pt_JESUp_EC2_2018_branch) fatJet1Pt_JESUp_EC2_2018_branch->SetAddress(&fatJet1Pt_JESUp_EC2_2018_);
  fatJet1PtOverMHH_JESUp_EC2_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_EC2_2018");
  if(fatJet1PtOverMHH_JESUp_EC2_2018_branch) fatJet1PtOverMHH_JESUp_EC2_2018_branch->SetAddress(&fatJet1PtOverMHH_JESUp_EC2_2018_);
  fatJet1Pt_JESDown_EC2_2018_branch = tree->GetBranch("fatJet1Pt_JESDown_EC2_2018");
  if(fatJet1Pt_JESDown_EC2_2018_branch) fatJet1Pt_JESDown_EC2_2018_branch->SetAddress(&fatJet1Pt_JESDown_EC2_2018_);
  fatJet1PtOverMHH_JESDown_EC2_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_EC2_2018");
  if(fatJet1PtOverMHH_JESDown_EC2_2018_branch) fatJet1PtOverMHH_JESDown_EC2_2018_branch->SetAddress(&fatJet1PtOverMHH_JESDown_EC2_2018_);
  fatJet1Pt_JESUp_FlavQCD_branch = tree->GetBranch("fatJet1Pt_JESUp_FlavQCD");
  if(fatJet1Pt_JESUp_FlavQCD_branch) fatJet1Pt_JESUp_FlavQCD_branch->SetAddress(&fatJet1Pt_JESUp_FlavQCD_);
  fatJet1PtOverMHH_JESUp_FlavQCD_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_FlavQCD");
  if(fatJet1PtOverMHH_JESUp_FlavQCD_branch) fatJet1PtOverMHH_JESUp_FlavQCD_branch->SetAddress(&fatJet1PtOverMHH_JESUp_FlavQCD_);
  fatJet1Pt_JESDown_FlavQCD_branch = tree->GetBranch("fatJet1Pt_JESDown_FlavQCD");
  if(fatJet1Pt_JESDown_FlavQCD_branch) fatJet1Pt_JESDown_FlavQCD_branch->SetAddress(&fatJet1Pt_JESDown_FlavQCD_);
  fatJet1PtOverMHH_JESDown_FlavQCD_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_FlavQCD");
  if(fatJet1PtOverMHH_JESDown_FlavQCD_branch) fatJet1PtOverMHH_JESDown_FlavQCD_branch->SetAddress(&fatJet1PtOverMHH_JESDown_FlavQCD_);
  fatJet1Pt_JESUp_HF_branch = tree->GetBranch("fatJet1Pt_JESUp_HF");
  if(fatJet1Pt_JESUp_HF_branch) fatJet1Pt_JESUp_HF_branch->SetAddress(&fatJet1Pt_JESUp_HF_);
  fatJet1PtOverMHH_JESUp_HF_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_HF");
  if(fatJet1PtOverMHH_JESUp_HF_branch) fatJet1PtOverMHH_JESUp_HF_branch->SetAddress(&fatJet1PtOverMHH_JESUp_HF_);
  fatJet1Pt_JESDown_HF_branch = tree->GetBranch("fatJet1Pt_JESDown_HF");
  if(fatJet1Pt_JESDown_HF_branch) fatJet1Pt_JESDown_HF_branch->SetAddress(&fatJet1Pt_JESDown_HF_);
  fatJet1PtOverMHH_JESDown_HF_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_HF");
  if(fatJet1PtOverMHH_JESDown_HF_branch) fatJet1PtOverMHH_JESDown_HF_branch->SetAddress(&fatJet1PtOverMHH_JESDown_HF_);
  fatJet1Pt_JESUp_HF_2016_branch = tree->GetBranch("fatJet1Pt_JESUp_HF_2016");
  if(fatJet1Pt_JESUp_HF_2016_branch) fatJet1Pt_JESUp_HF_2016_branch->SetAddress(&fatJet1Pt_JESUp_HF_2016_);//fatJet1Pt_JESUp_HF_2016);
  fatJet1PtOverMHH_JESUp_HF_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_HF_2016");
  if(fatJet1PtOverMHH_JESUp_HF_2016_branch) fatJet1PtOverMHH_JESUp_HF_2016_branch->SetAddress(&fatJet1PtOverMHH_JESUp_HF_2016_);//fatJet1PtOverMHH_JESUp_HF_2016);
  fatJet1Pt_JESDown_HF_2016_branch = tree->GetBranch("fatJet1Pt_JESDown_HF_2016");
  if(fatJet1Pt_JESDown_HF_2016_branch) fatJet1Pt_JESDown_HF_2016_branch->SetAddress(&fatJet1Pt_JESDown_HF_2016_);//fatJet1Pt_JESDown_HF_2016);
  fatJet1PtOverMHH_JESDown_HF_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_HF_2016");
  if(fatJet1PtOverMHH_JESDown_HF_2016_branch) fatJet1PtOverMHH_JESDown_HF_2016_branch->SetAddress(&fatJet1PtOverMHH_JESDown_HF_2016_);//fatJet1PtOverMHH_JESDown_HF_2016);
  fatJet1Pt_JESUp_HF_2017_branch = tree->GetBranch("fatJet1Pt_JESUp_HF_2017");
  if(fatJet1Pt_JESUp_HF_2017_branch) fatJet1Pt_JESUp_HF_2017_branch->SetAddress(&fatJet1Pt_JESUp_HF_2017_);//fatJet1Pt_JESUp_HF_2017);
  fatJet1PtOverMHH_JESUp_HF_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_HF_2017");
  if(fatJet1PtOverMHH_JESUp_HF_2017_branch) fatJet1PtOverMHH_JESUp_HF_2017_branch->SetAddress(&fatJet1PtOverMHH_JESUp_HF_2017_);//fatJet1PtOverMHH_JESUp_HF_2017);
  fatJet1Pt_JESDown_HF_2017_branch = tree->GetBranch("fatJet1Pt_JESDown_HF_2017");
  if(fatJet1Pt_JESDown_HF_2017_branch) fatJet1Pt_JESDown_HF_2017_branch->SetAddress(&fatJet1Pt_JESDown_HF_2017_);//fatJet1Pt_JESDown_HF_2017);
  fatJet1PtOverMHH_JESDown_HF_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_HF_2017");
  if(fatJet1PtOverMHH_JESDown_HF_2017_branch) fatJet1PtOverMHH_JESDown_HF_2017_branch->SetAddress(&fatJet1PtOverMHH_JESDown_HF_2017_);//fatJet1PtOverMHH_JESDown_HF_2017);
  fatJet1Pt_JESUp_HF_2018_branch = tree->GetBranch("fatJet1Pt_JESUp_HF_2018");
  if(fatJet1Pt_JESUp_HF_2018_branch) fatJet1Pt_JESUp_HF_2018_branch->SetAddress(&fatJet1Pt_JESUp_HF_2018_);//fatJet1Pt_JESUp_HF_2018);
  fatJet1PtOverMHH_JESUp_HF_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_HF_2018");
  if(fatJet1PtOverMHH_JESUp_HF_2018_branch) fatJet1PtOverMHH_JESUp_HF_2018_branch->SetAddress(&fatJet1PtOverMHH_JESUp_HF_2018_);//fatJet1PtOverMHH_JESUp_HF_2018);
  fatJet1Pt_JESDown_HF_2018_branch = tree->GetBranch("fatJet1Pt_JESDown_HF_2018");
  if(fatJet1Pt_JESDown_HF_2018_branch) fatJet1Pt_JESDown_HF_2018_branch->SetAddress(&fatJet1Pt_JESDown_HF_2018_);//fatJet1Pt_JESDown_HF_2018);
  fatJet1PtOverMHH_JESDown_HF_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_HF_2018");
  if(fatJet1PtOverMHH_JESDown_HF_2018_branch) fatJet1PtOverMHH_JESDown_HF_2018_branch->SetAddress(&fatJet1PtOverMHH_JESDown_HF_2018_);//fatJet1PtOverMHH_JESDown_HF_2018);
  fatJet1Pt_JESUp_RelBal_branch = tree->GetBranch("fatJet1Pt_JESUp_RelBal");
  if(fatJet1Pt_JESUp_RelBal_branch) fatJet1Pt_JESUp_RelBal_branch->SetAddress(&fatJet1Pt_JESUp_RelBal_);//fatJet1Pt_JESUp_RelBal);
  fatJet1PtOverMHH_JESUp_RelBal_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_RelBal");
  if(fatJet1PtOverMHH_JESUp_RelBal_branch) fatJet1PtOverMHH_JESUp_RelBal_branch->SetAddress(&fatJet1PtOverMHH_JESUp_RelBal_);//fatJet1PtOverMHH_JESUp_RelBal);
  fatJet1Pt_JESDown_RelBal_branch = tree->GetBranch("fatJet1Pt_JESDown_RelBal");
  if(fatJet1Pt_JESDown_RelBal_branch) fatJet1Pt_JESDown_RelBal_branch->SetAddress(&fatJet1Pt_JESDown_RelBal_);//fatJet1Pt_JESDown_RelBal);
 fatJet1PtOverMHH_JESDown_RelBal_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_RelBal");
 if(fatJet1PtOverMHH_JESDown_RelBal_branch) fatJet1PtOverMHH_JESDown_RelBal_branch->SetAddress(&fatJet1PtOverMHH_JESDown_RelBal_);//fatJet1PtOverMHH_JESDown_RelBal);
 fatJet1Pt_JESUp_RelSample_2016_branch = tree->GetBranch("fatJet1Pt_JESUp_RelSample_2016");
 if(fatJet1Pt_JESUp_RelSample_2016_branch) fatJet1Pt_JESUp_RelSample_2016_branch->SetAddress(&fatJet1Pt_JESUp_RelSample_2016_);//fatJet1Pt_JESUp_RelSample_2016);
 fatJet1PtOverMHH_JESUp_RelSample_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_RelSample_2016");
 if(fatJet1PtOverMHH_JESUp_RelSample_2016_branch) fatJet1PtOverMHH_JESUp_RelSample_2016_branch->SetAddress(&fatJet1PtOverMHH_JESUp_RelSample_2016_);//fatJet1PtOverMHH_JESUp_RelSample_2016);
 fatJet1Pt_JESDown_RelSample_2016_branch = tree->GetBranch("fatJet1Pt_JESDown_RelSample_2016");
 if(fatJet1Pt_JESDown_RelSample_2016_branch) fatJet1Pt_JESDown_RelSample_2016_branch->SetAddress(&fatJet1Pt_JESDown_RelSample_2016_);//fatJet1Pt_JESDown_RelSample_2016);
 fatJet1PtOverMHH_JESDown_RelSample_2016_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_RelSample_2016");
 if(fatJet1PtOverMHH_JESDown_RelSample_2016_branch) fatJet1PtOverMHH_JESDown_RelSample_2016_branch->SetAddress(&fatJet1PtOverMHH_JESDown_RelSample_2016_);//fatJet1PtOverMHH_JESDown_RelSample_2016);
 fatJet1Pt_JESUp_RelSample_2017_branch = tree->GetBranch("fatJet1Pt_JESUp_RelSample_2017");
 if(fatJet1Pt_JESUp_RelSample_2017_branch) fatJet1Pt_JESUp_RelSample_2017_branch->SetAddress(&fatJet1Pt_JESUp_RelSample_2017_);//fatJet1Pt_JESUp_RelSample_2017);
 fatJet1PtOverMHH_JESUp_RelSample_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_RelSample_2017");
 if(fatJet1PtOverMHH_JESUp_RelSample_2017_branch) fatJet1PtOverMHH_JESUp_RelSample_2017_branch->SetAddress(&fatJet1PtOverMHH_JESUp_RelSample_2017_);//fatJet1PtOverMHH_JESUp_RelSample_2017);
 fatJet1Pt_JESDown_RelSample_2017_branch = tree->GetBranch("fatJet1Pt_JESDown_RelSample_2017");
 if(fatJet1Pt_JESDown_RelSample_2017_branch) fatJet1Pt_JESDown_RelSample_2017_branch->SetAddress(&fatJet1Pt_JESDown_RelSample_2017_);//fatJet1Pt_JESDown_RelSample_2017);
 fatJet1PtOverMHH_JESDown_RelSample_2017_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_RelSample_2017");
 if(fatJet1PtOverMHH_JESDown_RelSample_2017_branch) fatJet1PtOverMHH_JESDown_RelSample_2017_branch->SetAddress(&fatJet1PtOverMHH_JESDown_RelSample_2017_);//fatJet1PtOverMHH_JESDown_RelSample_2017);
  fatJet1Pt_JESUp_RelSample_2018_branch = tree->GetBranch("fatJet1Pt_JESUp_RelSample_2018");
 if(fatJet1Pt_JESUp_RelSample_2018_branch) fatJet1Pt_JESUp_RelSample_2018_branch->SetAddress(&fatJet1Pt_JESUp_RelSample_2018_);//fatJet1Pt_JESUp_RelSample_2018);
 fatJet1PtOverMHH_JESUp_RelSample_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESUp_RelSample_2018");
 if(fatJet1PtOverMHH_JESUp_RelSample_2018_branch) fatJet1PtOverMHH_JESUp_RelSample_2018_branch->SetAddress(&fatJet1PtOverMHH_JESUp_RelSample_2018_);//fatJet1PtOverMHH_JESUp_RelSample_2018);
 fatJet1Pt_JESDown_RelSample_2018_branch = tree->GetBranch("fatJet1Pt_JESDown_RelSample_2018");
 if(fatJet1Pt_JESDown_RelSample_2018_branch) fatJet1Pt_JESDown_RelSample_2018_branch->SetAddress(&fatJet1Pt_JESDown_RelSample_2018_);//fatJet1Pt_JESDown_RelSample_2018);
 fatJet1PtOverMHH_JESDown_RelSample_2018_branch = tree->GetBranch("fatJet1PtOverMHH_JESDown_RelSample_2018");
 if(fatJet1PtOverMHH_JESDown_RelSample_2018_branch) fatJet1PtOverMHH_JESDown_RelSample_2018_branch->SetAddress(&fatJet1PtOverMHH_JESDown_RelSample_2018_);//fatJet1PtOverMHH_JESDown_RelSample_2018);


 fatJet2Pt_branch = tree->GetBranch("fatJet2Pt");
 if(fatJet2Pt_branch) fatJet2Pt_branch->SetAddress(&fatJet2Pt_);//
 fatJet2Eta_branch = tree->GetBranch("fatJet2Eta");
 if(fatJet2Eta_branch) fatJet2Eta_branch->SetAddress(&fatJet2Eta_);//
 fatJet2Phi_branch = tree->GetBranch("fatJet2Phi");
 if(fatJet2Phi_branch) fatJet2Phi_branch->SetAddress(&fatJet2Phi_);//
 fatJet2Mass_branch = tree->GetBranch("fatJet2Mass");
 if(fatJet2Mass_branch) fatJet2Mass_branch->SetAddress(&fatJet2Mass_);//
 fatJet2MassSD_branch = tree->GetBranch("fatJet2MassSD");
 if(fatJet2MassSD_branch) fatJet2MassSD_branch->SetAddress(&fatJet2MassSD_);//
 fatJet2MassSD_UnCorrected_branch = tree->GetBranch("fatJet2MassSD_UnCorrected");
 if(fatJet2MassSD_UnCorrected_branch) fatJet2MassSD_UnCorrected_branch->SetAddress(&fatJet2MassSD_UnCorrected_);//
 fatJet2MassRegressed_branch = tree->GetBranch("fatJet2MassRegressed");
 if(fatJet2MassRegressed_branch) fatJet2MassRegressed_branch->SetAddress(&fatJet2MassRegressed_);//
 fatJet2MassRegressed_UnCorrected_branch = tree->GetBranch("fatJet2MassRegressed_UnCorrected");
 if(fatJet2MassRegressed_UnCorrected_branch) fatJet2MassRegressed_UnCorrected_branch->SetAddress(&fatJet2MassRegressed_UnCorrected_);//
 fatJet2PNetXbb_branch = tree->GetBranch("fatJet2PNetXbb");
 if(fatJet2PNetXbb_branch) fatJet2PNetXbb_branch->SetAddress(&fatJet2PNetXbb_);//
 fatJet2PNetQCDb_branch = tree->GetBranch("fatJet2PNetQCDb");
 if(fatJet2PNetQCDb_branch) fatJet2PNetQCDb_branch->SetAddress(&fatJet2PNetQCDb_);//
 fatJet2PNetQCDbb_branch = tree->GetBranch("fatJet2PNetQCDbb");
 if(fatJet2PNetQCDbb_branch) fatJet2PNetQCDbb_branch->SetAddress(&fatJet2PNetQCDbb_);//
 fatJet2PNetQCDc_branch = tree->GetBranch("fatJet2PNetQCDc");
 if(fatJet2PNetQCDc_branch) fatJet2PNetQCDc_branch->SetAddress(&fatJet2PNetQCDc_);//
 fatJet2PNetQCDcc_branch = tree->GetBranch("fatJet2PNetQCDcc");
 if(fatJet2PNetQCDcc_branch) fatJet2PNetQCDcc_branch->SetAddress(&fatJet2PNetQCDcc_);//
 fatJet2PNetQCDothers_branch = tree->GetBranch("fatJet2PNetQCDothers");
 if(fatJet2PNetQCDothers_branch) fatJet2PNetQCDothers_branch->SetAddress(&fatJet2PNetQCDothers_);//
 fatJet2Tau3OverTau2_branch = tree->GetBranch("fatJet2Tau3OverTau2");
 if(fatJet2Tau3OverTau2_branch) fatJet2Tau3OverTau2_branch->SetAddress(&fatJet2Tau3OverTau2_);//
 fatJet2GenMatchIndex_branch = tree->GetBranch("fatJet2GenMatchIndex");
 if(fatJet2GenMatchIndex_branch) fatJet2GenMatchIndex_branch->SetAddress(&fatJet2GenMatchIndex_);//
 fatJet2HasMuon_branch = tree->GetBranch("fatJet2HasMuon");
 if(fatJet2HasMuon_branch) fatJet2HasMuon_branch->SetAddress(&fatJet2HasMuon_);//
 fatJet2HasElectron_branch = tree->GetBranch("fatJet2HasElectron");
 if(fatJet2HasElectron_branch) fatJet2HasElectron_branch->SetAddress(&fatJet2HasElectron_);//
 fatJet2HasBJetCSVLoose_branch = tree->GetBranch("fatJet2HasBJetCSVLoose");
 if(fatJet2HasBJetCSVLoose_branch) fatJet2HasBJetCSVLoose_branch->SetAddress(&fatJet2HasBJetCSVLoose_);//
 fatJet2HasBJetCSVMedium_branch = tree->GetBranch("fatJet2HasBJetCSVMedium");
 if(fatJet2HasBJetCSVMedium_branch) fatJet2HasBJetCSVMedium_branch->SetAddress(&fatJet2HasBJetCSVMedium_);//
 fatJet2HasBJetCSVTight_branch = tree->GetBranch("fatJet2HasBJetCSVTight");
 if(fatJet2HasBJetCSVTight_branch) fatJet2HasBJetCSVTight_branch->SetAddress(&fatJet2HasBJetCSVTight_);//
 fatJet2OppositeHemisphereHasBJet_branch = tree->GetBranch("fatJet2OppositeHemisphereHasBJet");
 if(fatJet2OppositeHemisphereHasBJet_branch) fatJet2OppositeHemisphereHasBJet_branch->SetAddress(&fatJet2OppositeHemisphereHasBJet_);//
 fatJet2PtOverMHH_branch = tree->GetBranch("fatJet2PtOverMHH");
 if(fatJet2PtOverMHH_branch) fatJet2PtOverMHH_branch->SetAddress(&fatJet2PtOverMHH_);//
 fatJet2PtOverMSD_branch = tree->GetBranch("fatJet2PtOverMSD");
 if(fatJet2PtOverMSD_branch) fatJet2PtOverMSD_branch->SetAddress(&fatJet2PtOverMSD_);//
 fatJet2PtOverMRegressed_branch = tree->GetBranch("fatJet2PtOverMRegressed");
 if(fatJet2PtOverMRegressed_branch) fatJet2PtOverMRegressed_branch->SetAddress(&fatJet2PtOverMRegressed_);//
 fatJet2MassSD_JMS_Down_branch = tree->GetBranch("fatJet2MassSD_JMS_Down");
 if(fatJet2MassSD_JMS_Down_branch) fatJet2MassSD_JMS_Down_branch->SetAddress(&fatJet2MassSD_JMS_Down_);//
 fatJet2MassSD_JMS_Up_branch = tree->GetBranch("fatJet2MassSD_JMS_Up");
 if(fatJet2MassSD_JMS_Up_branch) fatJet2MassSD_JMS_Up_branch->SetAddress(&fatJet2MassSD_JMS_Up_);//
 fatJet2MassSD_JMR_Down_branch = tree->GetBranch("fatJet2MassSD_JMR_Down");
 if(fatJet2MassSD_JMR_Down_branch) fatJet2MassSD_JMR_Down_branch->SetAddress(&fatJet2MassSD_JMR_Down_);//
 fatJet2MassSD_JMR_Up_branch = tree->GetBranch("fatJet2MassSD_JMR_Up");
 if(fatJet2MassSD_JMR_Up_branch) fatJet2MassSD_JMR_Up_branch->SetAddress(&fatJet2MassSD_JMR_Up_);//
 fatJet2MassRegressed_JMS_Down_branch = tree->GetBranch("fatJet2MassRegressed_JMS_Down");
 if(fatJet2MassRegressed_JMS_Down_branch) fatJet2MassRegressed_JMS_Down_branch->SetAddress(&fatJet2MassRegressed_JMS_Down_);//
 fatJet2MassRegressed_JMS_Up_branch = tree->GetBranch("fatJet2MassRegressed_JMS_Up");
 if(fatJet2MassRegressed_JMS_Up_branch) fatJet2MassRegressed_JMS_Up_branch->SetAddress(&fatJet2MassRegressed_JMS_Up_);//
 fatJet2MassRegressed_JMR_Down_branch = tree->GetBranch("fatJet2MassRegressed_JMR_Down");
 if(fatJet2MassRegressed_JMR_Down_branch) fatJet2MassRegressed_JMR_Down_branch->SetAddress(&fatJet2MassRegressed_JMR_Down_);//
 fatJet2MassRegressed_JMR_Up_branch = tree->GetBranch("fatJet2MassRegressed_JMR_Up");
 if(fatJet2MassRegressed_JMR_Up_branch) fatJet2MassRegressed_JMR_Up_branch->SetAddress(&fatJet2MassRegressed_JMR_Up_);//
 fatJet2PtOverMHH_JMS_Down_branch = tree->GetBranch("fatJet2PtOverMHH_JMS_Down");
 if(fatJet2PtOverMHH_JMS_Down_branch) fatJet2PtOverMHH_JMS_Down_branch->SetAddress(&fatJet2PtOverMHH_JMS_Down_);//
 fatJet2PtOverMHH_JMS_Up_branch = tree->GetBranch("fatJet2PtOverMHH_JMS_Up");
 if(fatJet2PtOverMHH_JMS_Up_branch) fatJet2PtOverMHH_JMS_Up_branch->SetAddress(&fatJet2PtOverMHH_JMS_Up_);//
 fatJet2PtOverMHH_JMR_Down_branch = tree->GetBranch("fatJet2PtOverMHH_JMR_Down");
 if(fatJet2PtOverMHH_JMR_Down_branch) fatJet2PtOverMHH_JMR_Down_branch->SetAddress(&fatJet2PtOverMHH_JMR_Down_);//
 fatJet2PtOverMHH_JMR_Up_branch = tree->GetBranch("fatJet2PtOverMHH_JMR_Up");
 if(fatJet2PtOverMHH_JMR_Up_branch) fatJet2PtOverMHH_JMR_Up_branch->SetAddress(&fatJet2PtOverMHH_JMR_Up_);//
 fatJet2Pt_JERUp_branch = tree->GetBranch("fatJet2Pt_JERUp");
 if(fatJet2Pt_JERUp_branch) fatJet2Pt_JERUp_branch->SetAddress(&fatJet2Pt_JERUp_);//
 fatJet2PtOverMHH_JERUp_branch = tree->GetBranch("fatJet2PtOverMHH_JERUp");
 if(fatJet2PtOverMHH_JERUp_branch) fatJet2PtOverMHH_JERUp_branch->SetAddress(&fatJet2PtOverMHH_JERUp_);//
 fatJet2Pt_JERDown_branch = tree->GetBranch("fatJet2Pt_JERDown");
 if(fatJet2Pt_JERDown_branch) fatJet2Pt_JERDown_branch->SetAddress(&fatJet2Pt_JERDown_);//
 fatJet2PtOverMHH_JERDown_branch = tree->GetBranch("fatJet2PtOverMHH_JERDown");
 if(fatJet2PtOverMHH_JERDown_branch) fatJet2PtOverMHH_JERDown_branch->SetAddress(&fatJet2PtOverMHH_JERDown_);//
 fatJet2Pt_JESUp_branch = tree->GetBranch("fatJet2Pt_JESUp");
 if(fatJet2Pt_JESUp_branch) fatJet2Pt_JESUp_branch->SetAddress(&fatJet2Pt_JESUp_);//
 fatJet2PtOverMHH_JESUp_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp");
 if(fatJet2PtOverMHH_JESUp_branch) fatJet2PtOverMHH_JESUp_branch->SetAddress(&fatJet2PtOverMHH_JESUp_);//
 fatJet2Pt_JESDown_branch = tree->GetBranch("fatJet2Pt_JESDown");
 if(fatJet2Pt_JESDown_branch) fatJet2Pt_JESDown_branch->SetAddress(&fatJet2Pt_JESDown_);//
 fatJet2PtOverMHH_JESDown_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown");
 if(fatJet2PtOverMHH_JESDown_branch) fatJet2PtOverMHH_JESDown_branch->SetAddress(&fatJet2PtOverMHH_JESDown_);//
 fatJet2Pt_JESUp_Abs_branch = tree->GetBranch("fatJet2Pt_JESUp_Abs");
 if(fatJet2Pt_JESUp_Abs_branch) fatJet2Pt_JESUp_Abs_branch->SetAddress(&fatJet2Pt_JESUp_Abs_);//
 fatJet2PtOverMHH_JESUp_Abs_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_Abs");
 if(fatJet2PtOverMHH_JESUp_Abs_branch) fatJet2PtOverMHH_JESUp_Abs_branch->SetAddress(&fatJet2PtOverMHH_JESUp_Abs_);//fatJet2PtOverMHH_JESUp_Abs);
 fatJet2Pt_JESDown_Abs_branch = tree->GetBranch("fatJet2Pt_JESDown_Abs");
 if(fatJet2Pt_JESDown_Abs_branch) fatJet2Pt_JESDown_Abs_branch->SetAddress(&fatJet2Pt_JESDown_Abs_);//fatJet2Pt_JESDown_Abs);
 fatJet2PtOverMHH_JESDown_Abs_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_Abs");
 if(fatJet2PtOverMHH_JESDown_Abs_branch) fatJet2PtOverMHH_JESDown_Abs_branch->SetAddress(&fatJet2PtOverMHH_JESDown_Abs_);//fatJet2PtOverMHH_JESDown_Abs);
 fatJet2Pt_JESUp_Abs_2016_branch = tree->GetBranch("fatJet2Pt_JESUp_Abs_2016");
 if(fatJet2Pt_JESUp_Abs_2016_branch) fatJet2Pt_JESUp_Abs_2016_branch->SetAddress(&fatJet2Pt_JESUp_Abs_2016_);//fatJet2Pt_JESUp_Abs_2016);
 fatJet2PtOverMHH_JESUp_Abs_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_Abs_2016");
 if(fatJet2PtOverMHH_JESUp_Abs_2016_branch) fatJet2PtOverMHH_JESUp_Abs_2016_branch->SetAddress(&fatJet2PtOverMHH_JESUp_Abs_2016_);//fatJet2PtOverMHH_JESUp_Abs_2016);
 fatJet2Pt_JESDowb_Abs_2016_branch = tree->GetBranch("fatJet2Pt_JESDowb_Abs_2016");
 if(fatJet2Pt_JESDowb_Abs_2016_branch) fatJet2Pt_JESDowb_Abs_2016_branch->SetAddress(&fatJet2Pt_JESDowb_Abs_2016_);//fatJet2Pt_JESDowb_Abs_2016);
 fatJet2PtOverMHH_JESDowb_Abs_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESDowb_Abs_2016");
 if(fatJet2PtOverMHH_JESDowb_Abs_2016_branch) fatJet2PtOverMHH_JESDowb_Abs_2016_branch->SetAddress(&fatJet2PtOverMHH_JESDowb_Abs_2016_);//fatJet2PtOverMHH_JESDowb_Abs_2016);
 fatJet2Pt_JESUp_Abs_2017_branch = tree->GetBranch("fatJet2Pt_JESUp_Abs_2017");
  if(fatJet2Pt_JESUp_Abs_2017_branch) fatJet2Pt_JESUp_Abs_2017_branch->SetAddress(&fatJet2Pt_JESUp_Abs_2017_);
  fatJet2PtOverMHH_JESUp_Abs_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_Abs_2017");
  if(fatJet2PtOverMHH_JESUp_Abs_2017_branch) fatJet2PtOverMHH_JESUp_Abs_2017_branch->SetAddress(&fatJet2PtOverMHH_JESUp_Abs_2017_);
  fatJet2Pt_JESDowb_Abs_2017_branch = tree->GetBranch("fatJet2Pt_JESDowb_Abs_2017");
  if(fatJet2Pt_JESDowb_Abs_2017_branch) fatJet2Pt_JESDowb_Abs_2017_branch->SetAddress(&fatJet2Pt_JESDowb_Abs_2017_);
  fatJet2PtOverMHH_JESDowb_Abs_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESDowb_Abs_2017");
  if(fatJet2PtOverMHH_JESDowb_Abs_2017_branch) fatJet2PtOverMHH_JESDowb_Abs_2017_branch->SetAddress(&fatJet2PtOverMHH_JESDowb_Abs_2017_);
  fatJet2Pt_JESUp_Abs_2018_branch = tree->GetBranch("fatJet2Pt_JESUp_Abs_2018");
  if(fatJet2Pt_JESUp_Abs_2018_branch) fatJet2Pt_JESUp_Abs_2018_branch->SetAddress(&fatJet2Pt_JESUp_Abs_2018_);
  fatJet2PtOverMHH_JESUp_Abs_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_Abs_2018");
  if(fatJet2PtOverMHH_JESUp_Abs_2018_branch) fatJet2PtOverMHH_JESUp_Abs_2018_branch->SetAddress(&fatJet2PtOverMHH_JESUp_Abs_2018_);
  fatJet2Pt_JESDowb_Abs_2018_branch = tree->GetBranch("fatJet2Pt_JESDowb_Abs_2018");
  if(fatJet2Pt_JESDowb_Abs_2018_branch) fatJet2Pt_JESDowb_Abs_2018_branch->SetAddress(&fatJet2Pt_JESDowb_Abs_2018_);
  fatJet2PtOverMHH_JESDowb_Abs_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESDowb_Abs_2018");
  if(fatJet2PtOverMHH_JESDowb_Abs_2018_branch) fatJet2PtOverMHH_JESDowb_Abs_2018_branch->SetAddress(&fatJet2PtOverMHH_JESDowb_Abs_2018_);
 fatJet2Pt_JESUp_BBEC1_branch = tree->GetBranch("fatJet2Pt_JESUp_BBEC1");
 if(fatJet2Pt_JESUp_BBEC1_branch) fatJet2Pt_JESUp_BBEC1_branch->SetAddress(&fatJet2Pt_JESUp_BBEC1_);//fatJet2Pt_JESUp_BBEC1);
 fatJet2PtOverMHH_JESUp_BBEC1_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_BBEC1");
 if(fatJet2PtOverMHH_JESUp_BBEC1_branch) fatJet2PtOverMHH_JESUp_BBEC1_branch->SetAddress(&fatJet2PtOverMHH_JESUp_BBEC1_);//fatJet2PtOverMHH_JESUp_BBEC1);
 fatJet2Pt_JESDown_BBEC1_branch = tree->GetBranch("fatJet2Pt_JESDown_BBEC1");
 if(fatJet2Pt_JESDown_BBEC1_branch) fatJet2Pt_JESDown_BBEC1_branch->SetAddress(&fatJet2Pt_JESDown_BBEC1_);//fatJet2Pt_JESDown_BBEC1);
 fatJet2PtOverMHH_JESDown_BBEC1_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_BBEC1");
 if(fatJet2PtOverMHH_JESDown_BBEC1_branch) fatJet2PtOverMHH_JESDown_BBEC1_branch->SetAddress(&fatJet2PtOverMHH_JESDown_BBEC1_);//fatJet2PtOverMHH_JESDown_BBEC1);
 fatJet2Pt_JESUp_BBEC1_2016_branch = tree->GetBranch("fatJet2Pt_JESUp_BBEC1_2016");
 if(fatJet2Pt_JESUp_BBEC1_2016_branch) fatJet2Pt_JESUp_BBEC1_2016_branch->SetAddress(&fatJet2Pt_JESUp_BBEC1_2016_);//fatJet2Pt_JESUp_BBEC1_2016);
 fatJet2PtOverMHH_JESUp_BBEC1_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_BBEC1_2016");
 if(fatJet2PtOverMHH_JESUp_BBEC1_2016_branch) fatJet2PtOverMHH_JESUp_BBEC1_2016_branch->SetAddress(&fatJet2PtOverMHH_JESUp_BBEC1_2016_);//fatJet2PtOverMHH_JESUp_BBEC1_2016);
 fatJet2Pt_JESDown_BBEC1_2016_branch = tree->GetBranch("fatJet2Pt_JESDown_BBEC1_2016");
 if(fatJet2Pt_JESDown_BBEC1_2016_branch) fatJet2Pt_JESDown_BBEC1_2016_branch->SetAddress(&fatJet2Pt_JESDown_BBEC1_2016_);//fatJet2Pt_JESDown_BBEC1_2016);
 fatJet2PtOverMHH_JESDown_BBEC1_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_BBEC1_2016");
 if(fatJet2PtOverMHH_JESDown_BBEC1_2016_branch) fatJet2PtOverMHH_JESDown_BBEC1_2016_branch->SetAddress(&fatJet2PtOverMHH_JESDown_BBEC1_2016_);//fatJet2PtOverMHH_JESDown_BBEC1_2016);
 fatJet2Pt_JESUp_BBEC1_2017_branch = tree->GetBranch("fatJet2Pt_JESUp_BBEC1_2017");
  if(fatJet2Pt_JESUp_BBEC1_2017_branch) fatJet2Pt_JESUp_BBEC1_2017_branch->SetAddress(&fatJet2Pt_JESUp_BBEC1_2017_);
  fatJet2PtOverMHH_JESUp_BBEC1_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_BBEC1_2017");
  if(fatJet2PtOverMHH_JESUp_BBEC1_2017_branch) fatJet2PtOverMHH_JESUp_BBEC1_2017_branch->SetAddress(&fatJet2PtOverMHH_JESUp_BBEC1_2017_);
  fatJet2Pt_JESDown_BBEC1_2017_branch = tree->GetBranch("fatJet2Pt_JESDown_BBEC1_2017");
  if(fatJet2Pt_JESDown_BBEC1_2017_branch) fatJet2Pt_JESDown_BBEC1_2017_branch->SetAddress(&fatJet2Pt_JESDown_BBEC1_2017_);
  fatJet2PtOverMHH_JESDown_BBEC1_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_BBEC1_2017");
  if(fatJet2PtOverMHH_JESDown_BBEC1_2017_branch) fatJet2PtOverMHH_JESDown_BBEC1_2017_branch->SetAddress(&fatJet2PtOverMHH_JESDown_BBEC1_2017_);
  fatJet2Pt_JESUp_BBEC1_2018_branch = tree->GetBranch("fatJet2Pt_JESUp_BBEC1_2018");
  if(fatJet2Pt_JESUp_BBEC1_2018_branch) fatJet2Pt_JESUp_BBEC1_2018_branch->SetAddress(&fatJet2Pt_JESUp_BBEC1_2018_);
  fatJet2PtOverMHH_JESUp_BBEC1_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_BBEC1_2018");
  if(fatJet2PtOverMHH_JESUp_BBEC1_2018_branch) fatJet2PtOverMHH_JESUp_BBEC1_2018_branch->SetAddress(&fatJet2PtOverMHH_JESUp_BBEC1_2018_);
  fatJet2Pt_JESDown_BBEC1_2018_branch = tree->GetBranch("fatJet2Pt_JESDown_BBEC1_2018");
  if(fatJet2Pt_JESDown_BBEC1_2018_branch) fatJet2Pt_JESDown_BBEC1_2018_branch->SetAddress(&fatJet2Pt_JESDown_BBEC1_2018_);
  fatJet2PtOverMHH_JESDown_BBEC1_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_BBEC1_2018");
  if(fatJet2PtOverMHH_JESDown_BBEC1_2018_branch) fatJet2PtOverMHH_JESDown_BBEC1_2018_branch->SetAddress(&fatJet2PtOverMHH_JESDown_BBEC1_2018_);
 fatJet2Pt_JESUp_EC2_branch = tree->GetBranch("fatJet2Pt_JESUp_EC2");
 if(fatJet2Pt_JESUp_EC2_branch) fatJet2Pt_JESUp_EC2_branch->SetAddress(&fatJet2Pt_JESUp_EC2_);//fatJet2Pt_JESUp_EC2);
 fatJet2PtOverMHH_JESUp_EC2_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_EC2");
 if(fatJet2PtOverMHH_JESUp_EC2_branch) fatJet2PtOverMHH_JESUp_EC2_branch->SetAddress(&fatJet2PtOverMHH_JESUp_EC2_);//fatJet2PtOverMHH_JESUp_EC2);
 fatJet2Pt_JESDown_EC2_branch = tree->GetBranch("fatJet2Pt_JESDown_EC2");
 if(fatJet2Pt_JESDown_EC2_branch) fatJet2Pt_JESDown_EC2_branch->SetAddress(&fatJet2Pt_JESDown_EC2_);//fatJet2Pt_JESDown_EC2);
 fatJet2PtOverMHH_JESDown_EC2_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_EC2");
 if(fatJet2PtOverMHH_JESDown_EC2_branch) fatJet2PtOverMHH_JESDown_EC2_branch->SetAddress(&fatJet2PtOverMHH_JESDown_EC2_);//fatJet2PtOverMHH_JESDown_EC2);
 fatJet2Pt_JESUp_EC2_2016_branch = tree->GetBranch("fatJet2Pt_JESUp_EC2_2016");
 if(fatJet2Pt_JESUp_EC2_2016_branch) fatJet2Pt_JESUp_EC2_2016_branch->SetAddress(&fatJet2Pt_JESUp_EC2_2016_);//fatJet2Pt_JESUp_EC2_2016);
 fatJet2PtOverMHH_JESUp_EC2_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_EC2_2016");
 if(fatJet2PtOverMHH_JESUp_EC2_2016_branch) fatJet2PtOverMHH_JESUp_EC2_2016_branch->SetAddress(&fatJet2PtOverMHH_JESUp_EC2_2016_);//fatJet2PtOverMHH_JESUp_EC2_2016);
 fatJet2Pt_JESDown_EC2_2016_branch = tree->GetBranch("fatJet2Pt_JESDown_EC2_2016");
 if(fatJet2Pt_JESDown_EC2_2016_branch) fatJet2Pt_JESDown_EC2_2016_branch->SetAddress(&fatJet2Pt_JESDown_EC2_2016_);//fatJet2Pt_JESDown_EC2_2016);
 fatJet2PtOverMHH_JESDown_EC2_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_EC2_2016");
 if(fatJet2PtOverMHH_JESDown_EC2_2016_branch) fatJet2PtOverMHH_JESDown_EC2_2016_branch->SetAddress(&fatJet2PtOverMHH_JESDown_EC2_2016_);//fatJet2PtOverMHH_JESDown_EC2_2016);
 fatJet2Pt_JESUp_EC2_2017_branch = tree->GetBranch("fatJet2Pt_JESUp_EC2_2017");
  if(fatJet2Pt_JESUp_EC2_2017_branch) fatJet2Pt_JESUp_EC2_2017_branch->SetAddress(&fatJet2Pt_JESUp_EC2_2017_);
  fatJet2PtOverMHH_JESUp_EC2_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_EC2_2017");
  if(fatJet2PtOverMHH_JESUp_EC2_2017_branch) fatJet2PtOverMHH_JESUp_EC2_2017_branch->SetAddress(&fatJet2PtOverMHH_JESUp_EC2_2017_);
  fatJet2Pt_JESDown_EC2_2017_branch = tree->GetBranch("fatJet2Pt_JESDown_EC2_2017");
  if(fatJet2Pt_JESDown_EC2_2017_branch) fatJet2Pt_JESDown_EC2_2017_branch->SetAddress(&fatJet2Pt_JESDown_EC2_2017_);
  fatJet2PtOverMHH_JESDown_EC2_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_EC2_2017");
  if(fatJet2PtOverMHH_JESDown_EC2_2017_branch) fatJet2PtOverMHH_JESDown_EC2_2017_branch->SetAddress(&fatJet2PtOverMHH_JESDown_EC2_2017_);
  fatJet2Pt_JESUp_EC2_2018_branch = tree->GetBranch("fatJet2Pt_JESUp_EC2_2018");
  if(fatJet2Pt_JESUp_EC2_2018_branch) fatJet2Pt_JESUp_EC2_2018_branch->SetAddress(&fatJet2Pt_JESUp_EC2_2018_);
  fatJet2PtOverMHH_JESUp_EC2_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_EC2_2018");
  if(fatJet2PtOverMHH_JESUp_EC2_2018_branch) fatJet2PtOverMHH_JESUp_EC2_2018_branch->SetAddress(&fatJet2PtOverMHH_JESUp_EC2_2018_);
  fatJet2Pt_JESDown_EC2_2018_branch = tree->GetBranch("fatJet2Pt_JESDown_EC2_2018");
  if(fatJet2Pt_JESDown_EC2_2018_branch) fatJet2Pt_JESDown_EC2_2018_branch->SetAddress(&fatJet2Pt_JESDown_EC2_2018_);
  fatJet2PtOverMHH_JESDown_EC2_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_EC2_2018");
  if(fatJet2PtOverMHH_JESDown_EC2_2018_branch) fatJet2PtOverMHH_JESDown_EC2_2018_branch->SetAddress(&fatJet2PtOverMHH_JESDown_EC2_2018_);
 fatJet2Pt_JESUp_FlavQCD_branch = tree->GetBranch("fatJet2Pt_JESUp_FlavQCD");
 if(fatJet2Pt_JESUp_FlavQCD_branch) fatJet2Pt_JESUp_FlavQCD_branch->SetAddress(&fatJet2Pt_JESUp_FlavQCD_);//fatJet2Pt_JESUp_FlavQCD);
 fatJet2PtOverMHH_JESUp_FlavQCD_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_FlavQCD");
 if(fatJet2PtOverMHH_JESUp_FlavQCD_branch) fatJet2PtOverMHH_JESUp_FlavQCD_branch->SetAddress(&fatJet2PtOverMHH_JESUp_FlavQCD_);//fatJet2PtOverMHH_JESUp_FlavQCD);
 fatJet2Pt_JESDown_FlavQCD_branch = tree->GetBranch("fatJet2Pt_JESDown_FlavQCD");
 if(fatJet2Pt_JESDown_FlavQCD_branch) fatJet2Pt_JESDown_FlavQCD_branch->SetAddress(&fatJet2Pt_JESDown_FlavQCD_);//fatJet2Pt_JESDown_FlavQCD);
 fatJet2PtOverMHH_JESDown_FlavQCD_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_FlavQCD");
 if(fatJet2PtOverMHH_JESDown_FlavQCD_branch) fatJet2PtOverMHH_JESDown_FlavQCD_branch->SetAddress(&fatJet2PtOverMHH_JESDown_FlavQCD_);//fatJet2PtOverMHH_JESDown_FlavQCD);
 fatJet2Pt_JESUp_HF_branch = tree->GetBranch("fatJet2Pt_JESUp_HF");
 if(fatJet2Pt_JESUp_HF_branch) fatJet2Pt_JESUp_HF_branch->SetAddress(&fatJet2Pt_JESUp_HF_);//fatJet2Pt_JESUp_HF);
 fatJet2PtOverMHH_JESUp_HF_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_HF");
 if(fatJet2PtOverMHH_JESUp_HF_branch) fatJet2PtOverMHH_JESUp_HF_branch->SetAddress(&fatJet2PtOverMHH_JESUp_HF_);//fatJet2PtOverMHH_JESUp_HF);
 fatJet2Pt_JESDown_HF_branch = tree->GetBranch("fatJet2Pt_JESDown_HF");
 if(fatJet2Pt_JESDown_HF_branch) fatJet2Pt_JESDown_HF_branch->SetAddress(&fatJet2Pt_JESDown_HF_);//fatJet2Pt_JESDown_HF);
 fatJet2PtOverMHH_JESDown_HF_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_HF");
 if(fatJet2PtOverMHH_JESDown_HF_branch) fatJet2PtOverMHH_JESDown_HF_branch->SetAddress(&fatJet2PtOverMHH_JESDown_HF_);//fatJet2PtOverMHH_JESDown_HF);
 fatJet2Pt_JESUp_HF_2016_branch = tree->GetBranch("fatJet2Pt_JESUp_HF_2016");
 if(fatJet2Pt_JESUp_HF_2016_branch) fatJet2Pt_JESUp_HF_2016_branch->SetAddress(&fatJet2Pt_JESUp_HF_2016_);//fatJet2Pt_JESUp_HF_2016);
 fatJet2PtOverMHH_JESUp_HF_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_HF_2016");
 if(fatJet2PtOverMHH_JESUp_HF_2016_branch) fatJet2PtOverMHH_JESUp_HF_2016_branch->SetAddress(&fatJet2PtOverMHH_JESUp_HF_2016_);//fatJet2PtOverMHH_JESUp_HF_2016);
 fatJet2Pt_JESDown_HF_2016_branch = tree->GetBranch("fatJet2Pt_JESDown_HF_2016");
 if(fatJet2Pt_JESDown_HF_2016_branch) fatJet2Pt_JESDown_HF_2016_branch->SetAddress(&fatJet2Pt_JESDown_HF_2016_);//fatJet2Pt_JESDown_HF_2016);
 fatJet2PtOverMHH_JESDown_HF_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_HF_2016");
 if(fatJet2PtOverMHH_JESDown_HF_2016_branch) fatJet2PtOverMHH_JESDown_HF_2016_branch->SetAddress(&fatJet2PtOverMHH_JESDown_HF_2016_);//fatJet2PtOverMHH_JESDown_HF_2016);
 fatJet2Pt_JESUp_HF_2017_branch = tree->GetBranch("fatJet2Pt_JESUp_HF_2017");
  if(fatJet2Pt_JESUp_HF_2017_branch) fatJet2Pt_JESUp_HF_2017_branch->SetAddress(&fatJet2Pt_JESUp_HF_2017_);//fatJet2Pt_JESUp_HF_2017);
  fatJet2PtOverMHH_JESUp_HF_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_HF_2017");
  if(fatJet2PtOverMHH_JESUp_HF_2017_branch) fatJet2PtOverMHH_JESUp_HF_2017_branch->SetAddress(&fatJet2PtOverMHH_JESUp_HF_2017_);//fatJet2PtOverMHH_JESUp_HF_2017);
  fatJet2Pt_JESDown_HF_2017_branch = tree->GetBranch("fatJet2Pt_JESDown_HF_2017");
  if(fatJet2Pt_JESDown_HF_2017_branch) fatJet2Pt_JESDown_HF_2017_branch->SetAddress(&fatJet2Pt_JESDown_HF_2017_);//fatJet2Pt_JESDown_HF_2017);
  fatJet2PtOverMHH_JESDown_HF_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_HF_2017");
  if(fatJet2PtOverMHH_JESDown_HF_2017_branch) fatJet2PtOverMHH_JESDown_HF_2017_branch->SetAddress(&fatJet2PtOverMHH_JESDown_HF_2017_);//fatJet2PtOverMHH_JESDown_HF_2017);
  fatJet2Pt_JESUp_HF_2018_branch = tree->GetBranch("fatJet2Pt_JESUp_HF_2018");
  if(fatJet2Pt_JESUp_HF_2018_branch) fatJet2Pt_JESUp_HF_2018_branch->SetAddress(&fatJet2Pt_JESUp_HF_2018_);//fatJet2Pt_JESUp_HF_2018);
  fatJet2PtOverMHH_JESUp_HF_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_HF_2018");
  if(fatJet2PtOverMHH_JESUp_HF_2018_branch) fatJet2PtOverMHH_JESUp_HF_2018_branch->SetAddress(&fatJet2PtOverMHH_JESUp_HF_2018_);//fatJet2PtOverMHH_JESUp_HF_2018);
  fatJet2Pt_JESDown_HF_2018_branch = tree->GetBranch("fatJet2Pt_JESDown_HF_2018");
  if(fatJet2Pt_JESDown_HF_2018_branch) fatJet2Pt_JESDown_HF_2018_branch->SetAddress(&fatJet2Pt_JESDown_HF_2018_);//fatJet2Pt_JESDown_HF_2018);
  fatJet2PtOverMHH_JESDown_HF_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_HF_2018");
  if(fatJet2PtOverMHH_JESDown_HF_2018_branch) fatJet2PtOverMHH_JESDown_HF_2018_branch->SetAddress(&fatJet2PtOverMHH_JESDown_HF_2018_);//fatJet2PtOverMHH_JESDown_HF_2018);
 fatJet2Pt_JESUp_RelBal_branch = tree->GetBranch("fatJet2Pt_JESUp_RelBal");
 if(fatJet2Pt_JESUp_RelBal_branch) fatJet2Pt_JESUp_RelBal_branch->SetAddress(&fatJet2Pt_JESUp_RelBal_);//fatJet2Pt_JESUp_RelBal);
 fatJet2PtOverMHH_JESUp_RelBal_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_RelBal");
 if(fatJet2PtOverMHH_JESUp_RelBal_branch) fatJet2PtOverMHH_JESUp_RelBal_branch->SetAddress(&fatJet2PtOverMHH_JESUp_RelBal_);//fatJet2PtOverMHH_JESUp_RelBal);
 fatJet2Pt_JESDown_RelBal_branch = tree->GetBranch("fatJet2Pt_JESDown_RelBal");
 if(fatJet2Pt_JESDown_RelBal_branch) fatJet2Pt_JESDown_RelBal_branch->SetAddress(&fatJet2Pt_JESDown_RelBal_);//fatJet2Pt_JESDown_RelBal);
 fatJet2PtOverMHH_JESDown_RelBal_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_RelBal");
 if(fatJet2PtOverMHH_JESDown_RelBal_branch) fatJet2PtOverMHH_JESDown_RelBal_branch->SetAddress(&fatJet2PtOverMHH_JESDown_RelBal_);//fatJet2PtOverMHH_JESDown_RelBal);
 fatJet2Pt_JESUp_RelSample_2016_branch = tree->GetBranch("fatJet2Pt_JESUp_RelSample_2016");
 if(fatJet2Pt_JESUp_RelSample_2016_branch) fatJet2Pt_JESUp_RelSample_2016_branch->SetAddress(&fatJet2Pt_JESUp_RelSample_2016_);//fatJet2Pt_JESUp_RelSample_2016);
 fatJet2PtOverMHH_JESUp_RelSample_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_RelSample_2016");
 if(fatJet2PtOverMHH_JESUp_RelSample_2016_branch) fatJet2PtOverMHH_JESUp_RelSample_2016_branch->SetAddress(&fatJet2PtOverMHH_JESUp_RelSample_2016_);//fatJet2PtOverMHH_JESUp_RelSample_2016);
 fatJet2Pt_JESDown_RelSample_2016_branch = tree->GetBranch("fatJet2Pt_JESDown_RelSample_2016");
 if(fatJet2Pt_JESDown_RelSample_2016_branch) fatJet2Pt_JESDown_RelSample_2016_branch->SetAddress(&fatJet2Pt_JESDown_RelSample_2016_);//fatJet2Pt_JESDown_RelSample_2016);
 fatJet2PtOverMHH_JESDown_RelSample_2016_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_RelSample_2016");
 if(fatJet2PtOverMHH_JESDown_RelSample_2016_branch) fatJet2PtOverMHH_JESDown_RelSample_2016_branch->SetAddress(&fatJet2PtOverMHH_JESDown_RelSample_2016_);//fatJet2PtOverMHH_JESDown_RelSample_2016);
fatJet2Pt_JESUp_RelSample_2017_branch = tree->GetBranch("fatJet2Pt_JESUp_RelSample_2017");
 if(fatJet2Pt_JESUp_RelSample_2017_branch) fatJet2Pt_JESUp_RelSample_2017_branch->SetAddress(&fatJet2Pt_JESUp_RelSample_2017_);//fatJet2Pt_JESUp_RelSample_2017);
 fatJet2PtOverMHH_JESUp_RelSample_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_RelSample_2017");
 if(fatJet2PtOverMHH_JESUp_RelSample_2017_branch) fatJet2PtOverMHH_JESUp_RelSample_2017_branch->SetAddress(&fatJet2PtOverMHH_JESUp_RelSample_2017_);//fatJet2PtOverMHH_JESUp_RelSample_2017);
 fatJet2Pt_JESDown_RelSample_2017_branch = tree->GetBranch("fatJet2Pt_JESDown_RelSample_2017");
 if(fatJet2Pt_JESDown_RelSample_2017_branch) fatJet2Pt_JESDown_RelSample_2017_branch->SetAddress(&fatJet2Pt_JESDown_RelSample_2017_);//fatJet2Pt_JESDown_RelSample_2017);
 fatJet2PtOverMHH_JESDown_RelSample_2017_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_RelSample_2017");
 if(fatJet2PtOverMHH_JESDown_RelSample_2017_branch) fatJet2PtOverMHH_JESDown_RelSample_2017_branch->SetAddress(&fatJet2PtOverMHH_JESDown_RelSample_2017_);//fatJet2PtOverMHH_JESDown_RelSample_2017);
  fatJet2Pt_JESUp_RelSample_2018_branch = tree->GetBranch("fatJet2Pt_JESUp_RelSample_2018");
 if(fatJet2Pt_JESUp_RelSample_2018_branch) fatJet2Pt_JESUp_RelSample_2018_branch->SetAddress(&fatJet2Pt_JESUp_RelSample_2018_);//fatJet2Pt_JESUp_RelSample_2018);
 fatJet2PtOverMHH_JESUp_RelSample_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESUp_RelSample_2018");
 if(fatJet2PtOverMHH_JESUp_RelSample_2018_branch) fatJet2PtOverMHH_JESUp_RelSample_2018_branch->SetAddress(&fatJet2PtOverMHH_JESUp_RelSample_2018_);//fatJet2PtOverMHH_JESUp_RelSample_2018);
 fatJet2Pt_JESDown_RelSample_2018_branch = tree->GetBranch("fatJet2Pt_JESDown_RelSample_2018");
 if(fatJet2Pt_JESDown_RelSample_2018_branch) fatJet2Pt_JESDown_RelSample_2018_branch->SetAddress(&fatJet2Pt_JESDown_RelSample_2018_);//fatJet2Pt_JESDown_RelSample_2018);
 fatJet2PtOverMHH_JESDown_RelSample_2018_branch = tree->GetBranch("fatJet2PtOverMHH_JESDown_RelSample_2018");
 if(fatJet2PtOverMHH_JESDown_RelSample_2018_branch) fatJet2PtOverMHH_JESDown_RelSample_2018_branch->SetAddress(&fatJet2PtOverMHH_JESDown_RelSample_2018_);//fatJet2PtOverMHH_JESDown_RelSample_2018);



 hh_pt_branch = tree->GetBranch("hh_pt");
 if(hh_pt_branch) hh_pt_branch->SetAddress(&hh_pt_);//
 hh_eta_branch = tree->GetBranch("hh_eta");
 if(hh_eta_branch) hh_eta_branch->SetAddress(&hh_eta_);//
 hh_phi_branch = tree->GetBranch("hh_phi");
 if(hh_phi_branch) hh_phi_branch->SetAddress(&hh_phi_);//
 hh_mass_branch = tree->GetBranch("hh_mass");
 if(hh_mass_branch) hh_mass_branch->SetAddress(&hh_mass_);//
 hh_pt_JMS_Down_branch = tree->GetBranch("hh_pt_JMS_Down");
 if(hh_pt_JMS_Down_branch) hh_pt_JMS_Down_branch->SetAddress(&hh_pt_JMS_Down_);//
 hh_pt_JMS_Up_branch = tree->GetBranch("hh_pt_JMS_Up");
 if(hh_pt_JMS_Up_branch) hh_pt_JMS_Up_branch->SetAddress(&hh_pt_JMS_Up_);//
 hh_eta_JMS_Down_branch = tree->GetBranch("hh_eta_JMS_Down");
 if(hh_eta_JMS_Down_branch) hh_eta_JMS_Down_branch->SetAddress(&hh_eta_JMS_Down_);//
 hh_eta_JMS_Up_branch = tree->GetBranch("hh_eta_JMS_Up");
 if(hh_eta_JMS_Up_branch) hh_eta_JMS_Up_branch->SetAddress(&hh_eta_JMS_Up_);//
 hh_mass_JMS_Down_branch = tree->GetBranch("hh_mass_JMS_Down");
 if(hh_mass_JMS_Down_branch) hh_mass_JMS_Down_branch->SetAddress(&hh_mass_JMS_Down_);//
 hh_mass_JMS_Up_branch = tree->GetBranch("hh_mass_JMS_Up");
 if(hh_mass_JMS_Up_branch) hh_mass_JMS_Up_branch->SetAddress(&hh_mass_JMS_Up_);//
 hh_pt_JMR_Down_branch = tree->GetBranch("hh_pt_JMR_Down");
 if(hh_pt_JMR_Down_branch) hh_pt_JMR_Down_branch->SetAddress(&hh_pt_JMR_Down_);//
 hh_pt_JMR_Up_branch = tree->GetBranch("hh_pt_JMR_Up");
 if(hh_pt_JMR_Up_branch) hh_pt_JMR_Up_branch->SetAddress(&hh_pt_JMR_Up_);//
 hh_eta_JMR_Down_branch = tree->GetBranch("hh_eta_JMR_Down");
 if(hh_eta_JMR_Down_branch) hh_eta_JMR_Down_branch->SetAddress(&hh_eta_JMR_Down_);//
 hh_eta_JMR_Up_branch = tree->GetBranch("hh_eta_JMR_Up");
 if(hh_eta_JMR_Up_branch) hh_eta_JMR_Up_branch->SetAddress(&hh_eta_JMR_Up_);//
 hh_mass_JMR_Down_branch = tree->GetBranch("hh_mass_JMR_Down");
 if(hh_mass_JMR_Down_branch) hh_mass_JMR_Down_branch->SetAddress(&hh_mass_JMR_Down_);//
 hh_mass_JMR_Up_branch = tree->GetBranch("hh_mass_JMR_Up");
 if(hh_mass_JMR_Up_branch) hh_mass_JMR_Up_branch->SetAddress(&hh_mass_JMR_Up_);//
 hh_pt_JERUp_branch = tree->GetBranch("hh_pt_JERUp");
 if(hh_pt_JERUp_branch) hh_pt_JERUp_branch->SetAddress(&hh_pt_JERUp_);//
 hh_eta_JERUp_branch = tree->GetBranch("hh_eta_JERUp");
 if(hh_eta_JERUp_branch) hh_eta_JERUp_branch->SetAddress(&hh_eta_JERUp_);//
 hh_mass_JERUp_branch = tree->GetBranch("hh_mass_JERUp");
 if(hh_mass_JERUp_branch) hh_mass_JERUp_branch->SetAddress(&hh_mass_JERUp_);//
 hh_pt_JERDown_branch = tree->GetBranch("hh_pt_JERDown");
 if(hh_pt_JERDown_branch) hh_pt_JERDown_branch->SetAddress(&hh_pt_JERDown_);//
 hh_eta_JERDown_branch = tree->GetBranch("hh_eta_JERDown");
 if(hh_eta_JERDown_branch) hh_eta_JERDown_branch->SetAddress(&hh_eta_JERDown_);//
 hh_mass_JERDown_branch = tree->GetBranch("hh_mass_JERDown");
 if(hh_mass_JERDown_branch) hh_mass_JERDown_branch->SetAddress(&hh_mass_JERDown_);//
 hh_pt_JESUp_branch = tree->GetBranch("hh_pt_JESUp");
 if(hh_pt_JESUp_branch) hh_pt_JESUp_branch->SetAddress(&hh_pt_JESUp_);//
 hh_eta_JESUp_branch = tree->GetBranch("hh_eta_JESUp");
 if(hh_eta_JESUp_branch) hh_eta_JESUp_branch->SetAddress(&hh_eta_JESUp_);//
 hh_mass_JESUp_branch = tree->GetBranch("hh_mass_JESUp");
 if(hh_mass_JESUp_branch) hh_mass_JESUp_branch->SetAddress(&hh_mass_JESUp_);//
 hh_pt_JESDown_branch = tree->GetBranch("hh_pt_JESDown");
 if(hh_pt_JESDown_branch) hh_pt_JESDown_branch->SetAddress(&hh_pt_JESDown_);//
 hh_eta_JESDown_branch = tree->GetBranch("hh_eta_JESDown");
 if(hh_eta_JESDown_branch) hh_eta_JESDown_branch->SetAddress(&hh_eta_JESDown_);//
 hh_mass_JESDown_branch = tree->GetBranch("hh_mass_JESDown");
 if(hh_mass_JESDown_branch) hh_mass_JESDown_branch->SetAddress(&hh_mass_JESDown_);//
 hh_pt_JESUp_Abs_branch = tree->GetBranch("hh_pt_JESUp_Abs");
 if(hh_pt_JESUp_Abs_branch) hh_pt_JESUp_Abs_branch->SetAddress(&hh_pt_JESUp_Abs_);//hh_pt_JESUp_Abs);
 hh_eta_JESUp_Abs_branch = tree->GetBranch("hh_eta_JESUp_Abs");
 if(hh_eta_JESUp_Abs_branch) hh_eta_JESUp_Abs_branch->SetAddress(&hh_eta_JESUp_Abs_);//hh_eta_JESUp_Abs);
 hh_mass_JESUp_Abs_branch = tree->GetBranch("hh_mass_JESUp_Abs");
 if(hh_mass_JESUp_Abs_branch) hh_mass_JESUp_Abs_branch->SetAddress(&hh_mass_JESUp_Abs_);//hh_mass_JESUp_Abs);
 hh_pt_JESDown_Abs_branch = tree->GetBranch("hh_pt_JESDown_Abs");
 if(hh_pt_JESDown_Abs_branch) hh_pt_JESDown_Abs_branch->SetAddress(&hh_pt_JESDown_Abs_);//hh_pt_JESDown_Abs);
 hh_eta_JESDown_Abs_branch = tree->GetBranch("hh_eta_JESDown_Abs");
 if(hh_eta_JESDown_Abs_branch) hh_eta_JESDown_Abs_branch->SetAddress(&hh_eta_JESDown_Abs_);//hh_eta_JESDown_Abs);
 hh_mass_JESDown_Abs_branch = tree->GetBranch("hh_mass_JESDown_Abs");
 if(hh_mass_JESDown_Abs_branch) hh_mass_JESDown_Abs_branch->SetAddress(&hh_mass_JESDown_Abs_);//hh_mass_JESDown_Abs);
 hh_pt_JESUp_Abs_2016_branch = tree->GetBranch("hh_pt_JESUp_Abs_2016");
 if(hh_pt_JESUp_Abs_2016_branch) hh_pt_JESUp_Abs_2016_branch->SetAddress(&hh_pt_JESUp_Abs_2016_);//hh_pt_JESUp_Abs_2016);
 hh_eta_JESUp_Abs_2016_branch = tree->GetBranch("hh_eta_JESUp_Abs_2016");
 if(hh_eta_JESUp_Abs_2016_branch) hh_eta_JESUp_Abs_2016_branch->SetAddress(&hh_eta_JESUp_Abs_2016_);//hh_eta_JESUp_Abs_2016);
 hh_mass_JESUp_Abs_2016_branch = tree->GetBranch("hh_mass_JESUp_Abs_2016");
 if(hh_mass_JESUp_Abs_2016_branch) hh_mass_JESUp_Abs_2016_branch->SetAddress(&hh_mass_JESUp_Abs_2016_);//hh_mass_JESUp_Abs_2016);
 hh_pt_JESDowb_Abs_2016_branch = tree->GetBranch("hh_pt_JESDowb_Abs_2016");
 if(hh_pt_JESDowb_Abs_2016_branch) hh_pt_JESDowb_Abs_2016_branch->SetAddress(&hh_pt_JESDowb_Abs_2016_);//hh_pt_JESDowb_Abs_2016);
 hh_eta_JESDowb_Abs_2016_branch = tree->GetBranch("hh_eta_JESDowb_Abs_2016");
 if(hh_eta_JESDowb_Abs_2016_branch) hh_eta_JESDowb_Abs_2016_branch->SetAddress(&hh_eta_JESDowb_Abs_2016_);//hh_eta_JESDowb_Abs_2016);
 hh_mass_JESDowb_Abs_2016_branch = tree->GetBranch("hh_mass_JESDowb_Abs_2016");
 if(hh_mass_JESDowb_Abs_2016_branch) hh_mass_JESDowb_Abs_2016_branch->SetAddress(&hh_mass_JESDowb_Abs_2016_);//hh_mass_JESDowb_Abs_2016);
 hh_pt_JESUp_Abs_2017_branch = tree->GetBranch("hh_pt_JESUp_Abs_2017");
 if(hh_pt_JESUp_Abs_2017_branch) hh_pt_JESUp_Abs_2017_branch->SetAddress(&hh_pt_JESUp_Abs_2017_);//hh_pt_JESUp_Abs_2017);
 hh_eta_JESUp_Abs_2017_branch = tree->GetBranch("hh_eta_JESUp_Abs_2017");
 if(hh_eta_JESUp_Abs_2017_branch) hh_eta_JESUp_Abs_2017_branch->SetAddress(&hh_eta_JESUp_Abs_2017_);//hh_eta_JESUp_Abs_2017);
 hh_mass_JESUp_Abs_2017_branch = tree->GetBranch("hh_mass_JESUp_Abs_2017");
 if(hh_mass_JESUp_Abs_2017_branch) hh_mass_JESUp_Abs_2017_branch->SetAddress(&hh_mass_JESUp_Abs_2017_);//hh_mass_JESUp_Abs_2017);
 hh_pt_JESDowb_Abs_2017_branch = tree->GetBranch("hh_pt_JESDowb_Abs_2017");
 if(hh_pt_JESDowb_Abs_2017_branch) hh_pt_JESDowb_Abs_2017_branch->SetAddress(&hh_pt_JESDowb_Abs_2017_);//hh_pt_JESDowb_Abs_2017);
 hh_eta_JESDowb_Abs_2017_branch = tree->GetBranch("hh_eta_JESDowb_Abs_2017");
 if(hh_eta_JESDowb_Abs_2017_branch) hh_eta_JESDowb_Abs_2017_branch->SetAddress(&hh_eta_JESDowb_Abs_2017_);//hh_eta_JESDowb_Abs_2017);
 hh_mass_JESDowb_Abs_2017_branch = tree->GetBranch("hh_mass_JESDowb_Abs_2017");
 if(hh_mass_JESDowb_Abs_2017_branch) hh_mass_JESDowb_Abs_2017_branch->SetAddress(&hh_mass_JESDowb_Abs_2017_);//hh_mass_JESDowb_Abs_2017);
 hh_pt_JESUp_Abs_2018_branch = tree->GetBranch("hh_pt_JESUp_Abs_2018");
 if(hh_pt_JESUp_Abs_2018_branch) hh_pt_JESUp_Abs_2018_branch->SetAddress(&hh_pt_JESUp_Abs_2018_);//hh_pt_JESUp_Abs_2018);
 hh_eta_JESUp_Abs_2018_branch = tree->GetBranch("hh_eta_JESUp_Abs_2018");
 if(hh_eta_JESUp_Abs_2018_branch) hh_eta_JESUp_Abs_2018_branch->SetAddress(&hh_eta_JESUp_Abs_2018_);//hh_eta_JESUp_Abs_2018);
 hh_mass_JESUp_Abs_2018_branch = tree->GetBranch("hh_mass_JESUp_Abs_2018");
 if(hh_mass_JESUp_Abs_2018_branch) hh_mass_JESUp_Abs_2018_branch->SetAddress(&hh_mass_JESUp_Abs_2018_);//hh_mass_JESUp_Abs_2018);
 hh_pt_JESDowb_Abs_2018_branch = tree->GetBranch("hh_pt_JESDowb_Abs_2018");
 if(hh_pt_JESDowb_Abs_2018_branch) hh_pt_JESDowb_Abs_2018_branch->SetAddress(&hh_pt_JESDowb_Abs_2018_);//hh_pt_JESDowb_Abs_2018);
 hh_eta_JESDowb_Abs_2018_branch = tree->GetBranch("hh_eta_JESDowb_Abs_2018");
 if(hh_eta_JESDowb_Abs_2018_branch) hh_eta_JESDowb_Abs_2018_branch->SetAddress(&hh_eta_JESDowb_Abs_2018_);//hh_eta_JESDowb_Abs_2018);
 hh_mass_JESDowb_Abs_2018_branch = tree->GetBranch("hh_mass_JESDowb_Abs_2018");
 if(hh_mass_JESDowb_Abs_2018_branch) hh_mass_JESDowb_Abs_2018_branch->SetAddress(&hh_mass_JESDowb_Abs_2018_);//hh_mass_JESDowb_Abs_2018);
 hh_pt_JESUp_BBEC1_branch = tree->GetBranch("hh_pt_JESUp_BBEC1");
 if(hh_pt_JESUp_BBEC1_branch) hh_pt_JESUp_BBEC1_branch->SetAddress(&hh_pt_JESUp_BBEC1_);//hh_pt_JESUp_BBEC1);
 hh_eta_JESUp_BBEC1_branch = tree->GetBranch("hh_eta_JESUp_BBEC1");
 if(hh_eta_JESUp_BBEC1_branch) hh_eta_JESUp_BBEC1_branch->SetAddress(&hh_eta_JESUp_BBEC1_);//hh_eta_JESUp_BBEC1);
 hh_eta_JESUp_BBEC1_branch = tree->GetBranch("hh_eta_JESUp_BBEC1");
 if(hh_eta_JESUp_BBEC1_branch) hh_eta_JESUp_BBEC1_branch->SetAddress(&hh_mass_JESUp_BBEC1_);//hh_mass_JESUp_BBEC1);
 hh_pt_JESDown_BBEC1_branch = tree->GetBranch("hh_pt_JESDown_BBEC1");
 if(hh_pt_JESDown_BBEC1_branch) hh_pt_JESDown_BBEC1_branch->SetAddress(&hh_pt_JESDown_BBEC1_);//hh_pt_JESDown_BBEC1);
 hh_eta_JESDown_BBEC1_branch = tree->GetBranch("hh_eta_JESDown_BBEC1");
 if(hh_eta_JESDown_BBEC1_branch) hh_eta_JESDown_BBEC1_branch->SetAddress(&hh_eta_JESDown_BBEC1_);//hh_eta_JESDown_BBEC1);
 hh_mass_JESDown_BBEC1_branch = tree->GetBranch("hh_mass_JESDown_BBEC1");
 if(hh_mass_JESDown_BBEC1_branch) hh_mass_JESDown_BBEC1_branch->SetAddress(&hh_mass_JESDown_BBEC1_);//hh_mass_JESDown_BBEC1);
 hh_pt_JESUp_BBEC1_2016_branch = tree->GetBranch("hh_pt_JESUp_BBEC1_2016");
 if(hh_pt_JESUp_BBEC1_2016_branch) hh_pt_JESUp_BBEC1_2016_branch->SetAddress(&hh_pt_JESUp_BBEC1_2016_);//hh_pt_JESUp_BBEC1_2016);
 hh_eta_JESUp_BBEC1_2016_branch = tree->GetBranch("hh_eta_JESUp_BBEC1_2016");
 if(hh_eta_JESUp_BBEC1_2016_branch) hh_eta_JESUp_BBEC1_2016_branch->SetAddress(&hh_eta_JESUp_BBEC1_2016_);//hh_eta_JESUp_BBEC1_2016);
 hh_mass_JESUp_BBEC1_2016_branch = tree->GetBranch("hh_mass_JESUp_BBEC1_2016");
 if(hh_mass_JESUp_BBEC1_2016_branch) hh_mass_JESUp_BBEC1_2016_branch->SetAddress(&hh_mass_JESUp_BBEC1_2016_);//hh_mass_JESUp_BBEC1_2016);
 hh_pt_JESDown_BBEC1_2016_branch = tree->GetBranch("hh_pt_JESDown_BBEC1_2016");
 if(hh_pt_JESDown_BBEC1_2016_branch) hh_pt_JESDown_BBEC1_2016_branch->SetAddress(&hh_pt_JESDown_BBEC1_2016_);//hh_pt_JESDown_BBEC1_2016);
 hh_eta_JESDown_BBEC1_2016_branch = tree->GetBranch("hh_eta_JESDown_BBEC1_2016");
 if(hh_eta_JESDown_BBEC1_2016_branch) hh_eta_JESDown_BBEC1_2016_branch->SetAddress(&hh_eta_JESDown_BBEC1_2016_);//hh_eta_JESDown_BBEC1_2016);
 hh_mass_JESDown_BBEC1_2016_branch = tree->GetBranch("hh_mass_JESDown_BBEC1_2016");
 if(hh_mass_JESDown_BBEC1_2016_branch) hh_mass_JESDown_BBEC1_2016_branch->SetAddress(&hh_mass_JESDown_BBEC1_2016_);//hh_mass_JESDown_BBEC1_2016);
 hh_pt_JESUp_BBEC1_2017_branch = tree->GetBranch("hh_pt_JESUp_BBEC1_2017");
 if(hh_pt_JESUp_BBEC1_2017_branch) hh_pt_JESUp_BBEC1_2017_branch->SetAddress(&hh_pt_JESUp_BBEC1_2017_);//hh_pt_JESUp_BBEC1_2017);
 hh_eta_JESUp_BBEC1_2017_branch = tree->GetBranch("hh_eta_JESUp_BBEC1_2017");
 if(hh_eta_JESUp_BBEC1_2017_branch) hh_eta_JESUp_BBEC1_2017_branch->SetAddress(&hh_eta_JESUp_BBEC1_2017_);//hh_eta_JESUp_BBEC1_2017);
 hh_mass_JESUp_BBEC1_2017_branch = tree->GetBranch("hh_mass_JESUp_BBEC1_2017");
 if(hh_mass_JESUp_BBEC1_2017_branch) hh_mass_JESUp_BBEC1_2017_branch->SetAddress(&hh_mass_JESUp_BBEC1_2017_);//hh_mass_JESUp_BBEC1_2017);
 hh_pt_JESDown_BBEC1_2017_branch = tree->GetBranch("hh_pt_JESDown_BBEC1_2017");
 if(hh_pt_JESDown_BBEC1_2017_branch) hh_pt_JESDown_BBEC1_2017_branch->SetAddress(&hh_pt_JESDown_BBEC1_2017_);//hh_pt_JESDown_BBEC1_2017);
 hh_eta_JESDown_BBEC1_2017_branch = tree->GetBranch("hh_eta_JESDown_BBEC1_2017");
 if(hh_eta_JESDown_BBEC1_2017_branch) hh_eta_JESDown_BBEC1_2017_branch->SetAddress(&hh_eta_JESDown_BBEC1_2017_);//hh_eta_JESDown_BBEC1_2017);
 hh_mass_JESDown_BBEC1_2017_branch = tree->GetBranch("hh_mass_JESDown_BBEC1_2017");
 if(hh_mass_JESDown_BBEC1_2017_branch) hh_mass_JESDown_BBEC1_2017_branch->SetAddress(&hh_mass_JESDown_BBEC1_2017_);//hh_mass_JESDown_BBEC1_2017);
 hh_pt_JESUp_BBEC1_2018_branch = tree->GetBranch("hh_pt_JESUp_BBEC1_2018");
 if(hh_pt_JESUp_BBEC1_2018_branch) hh_pt_JESUp_BBEC1_2018_branch->SetAddress(&hh_pt_JESUp_BBEC1_2018_);//hh_pt_JESUp_BBEC1_2018);
 hh_eta_JESUp_BBEC1_2018_branch = tree->GetBranch("hh_eta_JESUp_BBEC1_2018");
 if(hh_eta_JESUp_BBEC1_2018_branch) hh_eta_JESUp_BBEC1_2018_branch->SetAddress(&hh_eta_JESUp_BBEC1_2018_);//hh_eta_JESUp_BBEC1_2018);
 hh_mass_JESUp_BBEC1_2018_branch = tree->GetBranch("hh_mass_JESUp_BBEC1_2018");
 if(hh_mass_JESUp_BBEC1_2018_branch) hh_mass_JESUp_BBEC1_2018_branch->SetAddress(&hh_mass_JESUp_BBEC1_2018_);//hh_mass_JESUp_BBEC1_2018);
 hh_pt_JESDown_BBEC1_2018_branch = tree->GetBranch("hh_pt_JESDown_BBEC1_2018");
 if(hh_pt_JESDown_BBEC1_2018_branch) hh_pt_JESDown_BBEC1_2018_branch->SetAddress(&hh_pt_JESDown_BBEC1_2018_);//hh_pt_JESDown_BBEC1_2018);
 hh_eta_JESDown_BBEC1_2018_branch = tree->GetBranch("hh_eta_JESDown_BBEC1_2018");
 if(hh_eta_JESDown_BBEC1_2018_branch) hh_eta_JESDown_BBEC1_2018_branch->SetAddress(&hh_eta_JESDown_BBEC1_2018_);//hh_eta_JESDown_BBEC1_2018);
 hh_mass_JESDown_BBEC1_2018_branch = tree->GetBranch("hh_mass_JESDown_BBEC1_2018");
 if(hh_mass_JESDown_BBEC1_2018_branch) hh_mass_JESDown_BBEC1_2018_branch->SetAddress(&hh_mass_JESDown_BBEC1_2018_);//hh_mass_JESDown_BBEC1_2018);
 hh_pt_JESUp_EC2_branch = tree->GetBranch("hh_pt_JESUp_EC2");
 if(hh_pt_JESUp_EC2_branch) hh_pt_JESUp_EC2_branch->SetAddress(&hh_pt_JESUp_EC2_);//hh_pt_JESUp_EC2);
 hh_eta_JESUp_EC2_branch = tree->GetBranch("hh_eta_JESUp_EC2");
 if(hh_eta_JESUp_EC2_branch) hh_eta_JESUp_EC2_branch->SetAddress(&hh_eta_JESUp_EC2_);//hh_eta_JESUp_EC2);
 hh_mass_JESUp_EC2_branch = tree->GetBranch("hh_mass_JESUp_EC2");
 if(hh_mass_JESUp_EC2_branch) hh_mass_JESUp_EC2_branch->SetAddress(&hh_mass_JESUp_EC2_);//hh_mass_JESUp_EC2);
 hh_pt_JESDown_EC2_branch = tree->GetBranch("hh_pt_JESDown_EC2");
 if(hh_pt_JESDown_EC2_branch) hh_pt_JESDown_EC2_branch->SetAddress(&hh_pt_JESDown_EC2_);//hh_pt_JESDown_EC2);
 hh_eta_JESDown_EC2_branch = tree->GetBranch("hh_eta_JESDown_EC2");
 if(hh_eta_JESDown_EC2_branch) hh_eta_JESDown_EC2_branch->SetAddress(&hh_eta_JESDown_EC2_);//hh_eta_JESDown_EC2);
 hh_mass_JESDown_EC2_branch = tree->GetBranch("hh_mass_JESDown_EC2");
 if(hh_mass_JESDown_EC2_branch) hh_mass_JESDown_EC2_branch->SetAddress(&hh_mass_JESDown_EC2_);//hh_mass_JESDown_EC2);
 hh_pt_JESUp_EC2_2016_branch = tree->GetBranch("hh_pt_JESUp_EC2_2016");
 if(hh_pt_JESUp_EC2_2016_branch) hh_pt_JESUp_EC2_2016_branch->SetAddress(&hh_pt_JESUp_EC2_2016_);//hh_pt_JESUp_EC2_2016);
 hh_eta_JESUp_EC2_2016_branch = tree->GetBranch("hh_eta_JESUp_EC2_2016");
 if(hh_eta_JESUp_EC2_2016_branch) hh_eta_JESUp_EC2_2016_branch->SetAddress(&hh_eta_JESUp_EC2_2016_);//hh_eta_JESUp_EC2_2016);
 hh_mass_JESUp_EC2_2016_branch = tree->GetBranch("hh_mass_JESUp_EC2_2016");
 if(hh_mass_JESUp_EC2_2016_branch) hh_mass_JESUp_EC2_2016_branch->SetAddress(&hh_mass_JESUp_EC2_2016_);//hh_mass_JESUp_EC2_2016);
 hh_pt_JESDown_EC2_2016_branch = tree->GetBranch("hh_pt_JESDown_EC2_2016");
 if(hh_pt_JESDown_EC2_2016_branch) hh_pt_JESDown_EC2_2016_branch->SetAddress(&hh_pt_JESDown_EC2_2016_);//hh_pt_JESDown_EC2_2016);
 hh_eta_JESDown_EC2_2016_branch = tree->GetBranch("hh_eta_JESDown_EC2_2016");
 if(hh_eta_JESDown_EC2_2016_branch) hh_eta_JESDown_EC2_2016_branch->SetAddress(&hh_eta_JESDown_EC2_2016_);//hh_eta_JESDown_EC2_2016);
 hh_mass_JESDown_EC2_2016_branch = tree->GetBranch("hh_mass_JESDown_EC2_2016");
 if(hh_mass_JESDown_EC2_2016_branch) hh_mass_JESDown_EC2_2016_branch->SetAddress(&hh_mass_JESDown_EC2_2016_);//hh_mass_JESDown_EC2_2016);
 hh_pt_JESUp_EC2_2017_branch = tree->GetBranch("hh_pt_JESUp_EC2_2017");
 if(hh_pt_JESUp_EC2_2017_branch) hh_pt_JESUp_EC2_2017_branch->SetAddress(&hh_pt_JESUp_EC2_2017_);//hh_pt_JESUp_EC2_2017);
 hh_eta_JESUp_EC2_2017_branch = tree->GetBranch("hh_eta_JESUp_EC2_2017");
 if(hh_eta_JESUp_EC2_2017_branch) hh_eta_JESUp_EC2_2017_branch->SetAddress(&hh_eta_JESUp_EC2_2017_);//hh_eta_JESUp_EC2_2017);
 hh_mass_JESUp_EC2_2017_branch = tree->GetBranch("hh_mass_JESUp_EC2_2017");
 if(hh_mass_JESUp_EC2_2017_branch) hh_mass_JESUp_EC2_2017_branch->SetAddress(&hh_mass_JESUp_EC2_2017_);//hh_mass_JESUp_EC2_2017);
 hh_pt_JESDown_EC2_2017_branch = tree->GetBranch("hh_pt_JESDown_EC2_2017");
 if(hh_pt_JESDown_EC2_2017_branch) hh_pt_JESDown_EC2_2017_branch->SetAddress(&hh_pt_JESDown_EC2_2017_);//hh_pt_JESDown_EC2_2017);
 hh_eta_JESDown_EC2_2017_branch = tree->GetBranch("hh_eta_JESDown_EC2_2017");
 if(hh_eta_JESDown_EC2_2017_branch) hh_eta_JESDown_EC2_2017_branch->SetAddress(&hh_eta_JESDown_EC2_2017_);//hh_eta_JESDown_EC2_2017);
 hh_mass_JESDown_EC2_2017_branch = tree->GetBranch("hh_mass_JESDown_EC2_2017");
 if(hh_mass_JESDown_EC2_2017_branch) hh_mass_JESDown_EC2_2017_branch->SetAddress(&hh_mass_JESDown_EC2_2017_);//hh_mass_JESDown_EC2_2017);
 hh_pt_JESUp_EC2_2018_branch = tree->GetBranch("hh_pt_JESUp_EC2_2018");
 if(hh_pt_JESUp_EC2_2018_branch) hh_pt_JESUp_EC2_2018_branch->SetAddress(&hh_pt_JESUp_EC2_2018_);//hh_pt_JESUp_EC2_2018);
 hh_eta_JESUp_EC2_2018_branch = tree->GetBranch("hh_eta_JESUp_EC2_2018");
 if(hh_eta_JESUp_EC2_2018_branch) hh_eta_JESUp_EC2_2018_branch->SetAddress(&hh_eta_JESUp_EC2_2018_);//hh_eta_JESUp_EC2_2018);
 hh_mass_JESUp_EC2_2018_branch = tree->GetBranch("hh_mass_JESUp_EC2_2018");
 if(hh_mass_JESUp_EC2_2018_branch) hh_mass_JESUp_EC2_2018_branch->SetAddress(&hh_mass_JESUp_EC2_2018_);//hh_mass_JESUp_EC2_2018);
 hh_pt_JESDown_EC2_2018_branch = tree->GetBranch("hh_pt_JESDown_EC2_2018");
 if(hh_pt_JESDown_EC2_2018_branch) hh_pt_JESDown_EC2_2018_branch->SetAddress(&hh_pt_JESDown_EC2_2018_);//hh_pt_JESDown_EC2_2018);
 hh_eta_JESDown_EC2_2018_branch = tree->GetBranch("hh_eta_JESDown_EC2_2018");
 if(hh_eta_JESDown_EC2_2018_branch) hh_eta_JESDown_EC2_2018_branch->SetAddress(&hh_eta_JESDown_EC2_2018_);//hh_eta_JESDown_EC2_2018);
 hh_mass_JESDown_EC2_2018_branch = tree->GetBranch("hh_mass_JESDown_EC2_2018");
 if(hh_mass_JESDown_EC2_2018_branch) hh_mass_JESDown_EC2_2018_branch->SetAddress(&hh_mass_JESDown_EC2_2018_);//hh_mass_JESDown_EC2_2018);
 hh_pt_JESUp_FlavQCD_branch = tree->GetBranch("hh_pt_JESUp_FlavQCD");
 if(hh_pt_JESUp_FlavQCD_branch) hh_pt_JESUp_FlavQCD_branch->SetAddress(&hh_pt_JESUp_FlavQCD_);//hh_pt_JESUp_FlavQCD);
 hh_eta_JESUp_FlavQCD_branch = tree->GetBranch("hh_eta_JESUp_FlavQCD");
 if(hh_eta_JESUp_FlavQCD_branch) hh_eta_JESUp_FlavQCD_branch->SetAddress(&hh_eta_JESUp_FlavQCD_);//hh_eta_JESUp_FlavQCD);
 hh_mass_JESUp_FlavQCD_branch = tree->GetBranch("hh_mass_JESUp_FlavQCD");
 if(hh_mass_JESUp_FlavQCD_branch) hh_mass_JESUp_FlavQCD_branch->SetAddress(&hh_mass_JESUp_FlavQCD_);//hh_mass_JESUp_FlavQCD);
 hh_pt_JESDown_FlavQCD_branch = tree->GetBranch("hh_pt_JESDown_FlavQCD");
 if(hh_pt_JESDown_FlavQCD_branch) hh_pt_JESDown_FlavQCD_branch->SetAddress(&hh_pt_JESDown_FlavQCD_);//hh_pt_JESDown_FlavQCD);
 hh_eta_JESDown_FlavQCD_branch = tree->GetBranch("hh_eta_JESDown_FlavQCD");
 if(hh_eta_JESDown_FlavQCD_branch) hh_eta_JESDown_FlavQCD_branch->SetAddress(&hh_eta_JESDown_FlavQCD_);//hh_eta_JESDown_FlavQCD);
 hh_mass_JESDown_FlavQCD_branch = tree->GetBranch("hh_mass_JESDown_FlavQCD");
 if(hh_mass_JESDown_FlavQCD_branch) hh_mass_JESDown_FlavQCD_branch->SetAddress(&hh_mass_JESDown_FlavQCD_);//hh_mass_JESDown_FlavQCD);
 hh_pt_JESUp_HF_branch = tree->GetBranch("hh_pt_JESUp_HF");
 if(hh_pt_JESUp_HF_branch) hh_pt_JESUp_HF_branch->SetAddress(&hh_pt_JESUp_HF_);//hh_pt_JESUp_HF);
 hh_eta_JESUp_HF_branch = tree->GetBranch("hh_eta_JESUp_HF");
 if(hh_eta_JESUp_HF_branch) hh_eta_JESUp_HF_branch->SetAddress(&hh_eta_JESUp_HF_);//hh_eta_JESUp_HF);
 hh_mass_JESUp_HF_branch = tree->GetBranch("hh_mass_JESUp_HF");
 if(hh_mass_JESUp_HF_branch) hh_mass_JESUp_HF_branch->SetAddress(&hh_mass_JESUp_HF_);//hh_mass_JESUp_HF);
 hh_pt_JESDown_HF_branch = tree->GetBranch("hh_pt_JESDown_HF");
 if(hh_pt_JESDown_HF_branch) hh_pt_JESDown_HF_branch->SetAddress(&hh_pt_JESDown_HF_);//hh_pt_JESDown_HF);
 hh_eta_JESDown_HF_branch = tree->GetBranch("hh_eta_JESDown_HF");
 if(hh_eta_JESDown_HF_branch) hh_eta_JESDown_HF_branch->SetAddress(&hh_eta_JESDown_HF_);//hh_eta_JESDown_HF);
 hh_mass_JESDown_HF_branch = tree->GetBranch("hh_mass_JESDown_HF");
 if(hh_mass_JESDown_HF_branch) hh_mass_JESDown_HF_branch->SetAddress(&hh_mass_JESDown_HF_);//hh_mass_JESDown_HF);
 hh_pt_JESUp_HF_2016_branch = tree->GetBranch("hh_pt_JESUp_HF_2016");
 if(hh_pt_JESUp_HF_2016_branch) hh_pt_JESUp_HF_2016_branch->SetAddress(&hh_pt_JESUp_HF_2016_);//hh_pt_JESUp_HF_2016);
 hh_eta_JESUp_HF_2016_branch = tree->GetBranch("hh_eta_JESUp_HF_2016");
 if(hh_eta_JESUp_HF_2016_branch) hh_eta_JESUp_HF_2016_branch->SetAddress(&hh_eta_JESUp_HF_2016_);//hh_eta_JESUp_HF_2016);
 hh_mass_JESUp_HF_2016_branch = tree->GetBranch("hh_mass_JESUp_HF_2016");
if(hh_mass_JESUp_HF_2016_branch) hh_mass_JESUp_HF_2016_branch->SetAddress(&hh_mass_JESUp_HF_2016_);//hh_mass_JESUp_HF_2016);
 hh_pt_JESDown_HF_2016_branch = tree->GetBranch("hh_pt_JESDown_HF_2016");
 if(hh_pt_JESDown_HF_2016_branch) hh_pt_JESDown_HF_2016_branch->SetAddress(&hh_pt_JESDown_HF_2016_);//hh_pt_JESDown_HF_2016);
 hh_eta_JESDown_HF_2016_branch = tree->GetBranch("hh_eta_JESDown_HF_2016");
 if(hh_eta_JESDown_HF_2016_branch) hh_eta_JESDown_HF_2016_branch->SetAddress(&hh_eta_JESDown_HF_2016_);//hh_eta_JESDown_HF_2016);
 hh_mass_JESDown_HF_2016_branch = tree->GetBranch("hh_mass_JESDown_HF_2016");
 if(hh_mass_JESDown_HF_2016_branch) hh_mass_JESDown_HF_2016_branch->SetAddress(&hh_mass_JESDown_HF_2016_);//hh_mass_JESDown_HF_2016);
 hh_pt_JESUp_HF_2017_branch = tree->GetBranch("hh_pt_JESUp_HF_2017");
 if(hh_pt_JESUp_HF_2017_branch) hh_pt_JESUp_HF_2017_branch->SetAddress(&hh_pt_JESUp_HF_2017_);//hh_pt_JESUp_HF_2017);
 hh_eta_JESUp_HF_2017_branch = tree->GetBranch("hh_eta_JESUp_HF_2017");
 if(hh_eta_JESUp_HF_2017_branch) hh_eta_JESUp_HF_2017_branch->SetAddress(&hh_eta_JESUp_HF_2017_);//hh_eta_JESUp_HF_2017);
 hh_mass_JESUp_HF_2017_branch = tree->GetBranch("hh_mass_JESUp_HF_2017");
if(hh_mass_JESUp_HF_2017_branch) hh_mass_JESUp_HF_2017_branch->SetAddress(&hh_mass_JESUp_HF_2017_);//hh_mass_JESUp_HF_2017);
 hh_pt_JESDown_HF_2017_branch = tree->GetBranch("hh_pt_JESDown_HF_2017");
 if(hh_pt_JESDown_HF_2017_branch) hh_pt_JESDown_HF_2017_branch->SetAddress(&hh_pt_JESDown_HF_2017_);//hh_pt_JESDown_HF_2017);
 hh_eta_JESDown_HF_2017_branch = tree->GetBranch("hh_eta_JESDown_HF_2017");
 if(hh_eta_JESDown_HF_2017_branch) hh_eta_JESDown_HF_2017_branch->SetAddress(&hh_eta_JESDown_HF_2017_);//hh_eta_JESDown_HF_2017);
 hh_mass_JESDown_HF_2017_branch = tree->GetBranch("hh_mass_JESDown_HF_2017");
 if(hh_mass_JESDown_HF_2017_branch) hh_mass_JESDown_HF_2017_branch->SetAddress(&hh_mass_JESDown_HF_2017_);//hh_mass_JESDown_HF_2017);
 hh_pt_JESUp_HF_2018_branch = tree->GetBranch("hh_pt_JESUp_HF_2018");
 if(hh_pt_JESUp_HF_2018_branch) hh_pt_JESUp_HF_2018_branch->SetAddress(&hh_pt_JESUp_HF_2018_);//hh_pt_JESUp_HF_2018);
 hh_eta_JESUp_HF_2018_branch = tree->GetBranch("hh_eta_JESUp_HF_2018");
 if(hh_eta_JESUp_HF_2018_branch) hh_eta_JESUp_HF_2018_branch->SetAddress(&hh_eta_JESUp_HF_2018_);//hh_eta_JESUp_HF_2018);
 hh_mass_JESUp_HF_2018_branch = tree->GetBranch("hh_mass_JESUp_HF_2018");
if(hh_mass_JESUp_HF_2018_branch) hh_mass_JESUp_HF_2018_branch->SetAddress(&hh_mass_JESUp_HF_2018_);//hh_mass_JESUp_HF_2018);
 hh_pt_JESDown_HF_2018_branch = tree->GetBranch("hh_pt_JESDown_HF_2018");
 if(hh_pt_JESDown_HF_2018_branch) hh_pt_JESDown_HF_2018_branch->SetAddress(&hh_pt_JESDown_HF_2018_);//hh_pt_JESDown_HF_2018);
 hh_eta_JESDown_HF_2018_branch = tree->GetBranch("hh_eta_JESDown_HF_2018");
 if(hh_eta_JESDown_HF_2018_branch) hh_eta_JESDown_HF_2018_branch->SetAddress(&hh_eta_JESDown_HF_2018_);//hh_eta_JESDown_HF_2018);
 hh_mass_JESDown_HF_2018_branch = tree->GetBranch("hh_mass_JESDown_HF_2018");
 if(hh_mass_JESDown_HF_2018_branch) hh_mass_JESDown_HF_2018_branch->SetAddress(&hh_mass_JESDown_HF_2018_);//hh_mass_JESDown_HF_2018);
 hh_pt_JESUp_RelBal_branch = tree->GetBranch("hh_pt_JESUp_RelBal");
 if(hh_pt_JESUp_RelBal_branch) hh_pt_JESUp_RelBal_branch->SetAddress(&hh_pt_JESUp_RelBal_);//hh_pt_JESUp_RelBal);
 hh_eta_JESUp_RelBal_branch = tree->GetBranch("hh_eta_JESUp_RelBal");
 if(hh_eta_JESUp_RelBal_branch) hh_eta_JESUp_RelBal_branch->SetAddress(&hh_eta_JESUp_RelBal_);//hh_eta_JESUp_RelBal);
 hh_mass_JESUp_RelBal_branch = tree->GetBranch("hh_mass_JESUp_RelBal");
 if(hh_mass_JESUp_RelBal_branch) hh_mass_JESUp_RelBal_branch->SetAddress(&hh_mass_JESUp_RelBal_);//hh_mass_JESUp_RelBal);
 hh_pt_JESDown_RelBal_branch = tree->GetBranch("hh_pt_JESDown_RelBal");
 if(hh_pt_JESDown_RelBal_branch) hh_pt_JESDown_RelBal_branch->SetAddress(&hh_pt_JESDown_RelBal_);//hh_pt_JESDown_RelBal);
 hh_eta_JESDown_RelBal_branch = tree->GetBranch("hh_eta_JESDown_RelBal");
 if(hh_eta_JESDown_RelBal_branch) hh_eta_JESDown_RelBal_branch->SetAddress(&hh_eta_JESDown_RelBal_);//hh_eta_JESDown_RelBal);
 hh_mass_JESDown_RelBal_branch = tree->GetBranch("hh_mass_JESDown_RelBal");
 if(hh_mass_JESDown_RelBal_branch) hh_mass_JESDown_RelBal_branch->SetAddress(&hh_mass_JESDown_RelBal_);//hh_mass_JESDown_RelBal);
 hh_pt_JESUp_RelSample_2016_branch = tree->GetBranch("hh_pt_JESUp_RelSample_2016");
 if(hh_pt_JESUp_RelSample_2016_branch) hh_pt_JESUp_RelSample_2016_branch->SetAddress(&hh_pt_JESUp_RelSample_2016_);//hh_pt_JESUp_RelSample_2016);
 hh_eta_JESUp_RelSample_2016_branch = tree->GetBranch("hh_eta_JESUp_RelSample_2016");
 if(hh_eta_JESUp_RelSample_2016_branch) hh_eta_JESUp_RelSample_2016_branch->SetAddress(&hh_eta_JESUp_RelSample_2016_);//hh_eta_JESUp_RelSample_2016);
 hh_mass_JESUp_RelSample_2016_branch = tree->GetBranch("hh_mass_JESUp_RelSample_2016");
 if(hh_mass_JESUp_RelSample_2016_branch) hh_mass_JESUp_RelSample_2016_branch->SetAddress(&hh_mass_JESUp_RelSample_2016_);//hh_mass_JESUp_RelSample_2016);
 hh_pt_JESDown_RelSample_2016_branch = tree->GetBranch("hh_pt_JESDown_RelSample_2016");
 if(hh_pt_JESDown_RelSample_2016_branch) hh_pt_JESDown_RelSample_2016_branch->SetAddress(&hh_pt_JESDown_RelSample_2016_);//hh_pt_JESDown_RelSample_2016);
 hh_eta_JESDown_RelSample_2016_branch = tree->GetBranch("hh_eta_JESDown_RelSample_2016");
 if(hh_eta_JESDown_RelSample_2016_branch) hh_eta_JESDown_RelSample_2016_branch->SetAddress(&hh_eta_JESDown_RelSample_2016_);//hh_eta_JESDown_RelSample_2016);
 hh_mass_JESDown_RelSample_2016_branch = tree->GetBranch("hh_mass_JESDown_RelSample_2016");
 if(hh_mass_JESDown_RelSample_2016_branch) hh_mass_JESDown_RelSample_2016_branch->SetAddress(&hh_mass_JESDown_RelSample_2016_);//hh_mass_JESDown_RelSample_2016);
 hh_pt_JESUp_RelSample_2017_branch = tree->GetBranch("hh_pt_JESUp_RelSample_2017");
 if(hh_pt_JESUp_RelSample_2017_branch) hh_pt_JESUp_RelSample_2017_branch->SetAddress(&hh_pt_JESUp_RelSample_2017_);//hh_pt_JESUp_RelSample_2017);
 hh_eta_JESUp_RelSample_2017_branch = tree->GetBranch("hh_eta_JESUp_RelSample_2017");
 if(hh_eta_JESUp_RelSample_2017_branch) hh_eta_JESUp_RelSample_2017_branch->SetAddress(&hh_eta_JESUp_RelSample_2017_);//hh_eta_JESUp_RelSample_2017);
 hh_mass_JESUp_RelSample_2017_branch = tree->GetBranch("hh_mass_JESUp_RelSample_2017");
 if(hh_mass_JESUp_RelSample_2017_branch) hh_mass_JESUp_RelSample_2017_branch->SetAddress(&hh_mass_JESUp_RelSample_2017_);//hh_mass_JESUp_RelSample_2017);
 hh_pt_JESDown_RelSample_2017_branch = tree->GetBranch("hh_pt_JESDown_RelSample_2017");
 if(hh_pt_JESDown_RelSample_2017_branch) hh_pt_JESDown_RelSample_2017_branch->SetAddress(&hh_pt_JESDown_RelSample_2017_);//hh_pt_JESDown_RelSample_2017);
 hh_eta_JESDown_RelSample_2017_branch = tree->GetBranch("hh_eta_JESDown_RelSample_2017");
 if(hh_eta_JESDown_RelSample_2017_branch) hh_eta_JESDown_RelSample_2017_branch->SetAddress(&hh_eta_JESDown_RelSample_2017_);//hh_eta_JESDown_RelSample_2017);
 hh_mass_JESDown_RelSample_2017_branch = tree->GetBranch("hh_mass_JESDown_RelSample_2017");
 if(hh_mass_JESDown_RelSample_2017_branch) hh_mass_JESDown_RelSample_2017_branch->SetAddress(&hh_mass_JESDown_RelSample_2017_);//hh_mass_JESDown_RelSample_2017);
 hh_pt_JESUp_RelSample_2018_branch = tree->GetBranch("hh_pt_JESUp_RelSample_2018");
 if(hh_pt_JESUp_RelSample_2018_branch) hh_pt_JESUp_RelSample_2018_branch->SetAddress(&hh_pt_JESUp_RelSample_2018_);//hh_pt_JESUp_RelSample_2018);
 hh_eta_JESUp_RelSample_2018_branch = tree->GetBranch("hh_eta_JESUp_RelSample_2018");
 if(hh_eta_JESUp_RelSample_2018_branch) hh_eta_JESUp_RelSample_2018_branch->SetAddress(&hh_eta_JESUp_RelSample_2018_);//hh_eta_JESUp_RelSample_2018);
 hh_mass_JESUp_RelSample_2018_branch = tree->GetBranch("hh_mass_JESUp_RelSample_2018");
 if(hh_mass_JESUp_RelSample_2018_branch) hh_mass_JESUp_RelSample_2018_branch->SetAddress(&hh_mass_JESUp_RelSample_2018_);//hh_mass_JESUp_RelSample_2018);
 hh_pt_JESDown_RelSample_2018_branch = tree->GetBranch("hh_pt_JESDown_RelSample_2018");
 if(hh_pt_JESDown_RelSample_2018_branch) hh_pt_JESDown_RelSample_2018_branch->SetAddress(&hh_pt_JESDown_RelSample_2018_);//hh_pt_JESDown_RelSample_2018);
 hh_eta_JESDown_RelSample_2018_branch = tree->GetBranch("hh_eta_JESDown_RelSample_2018");
 if(hh_eta_JESDown_RelSample_2018_branch) hh_eta_JESDown_RelSample_2018_branch->SetAddress(&hh_eta_JESDown_RelSample_2018_);//hh_eta_JESDown_RelSample_2018);
 hh_mass_JESDown_RelSample_2018_branch = tree->GetBranch("hh_mass_JESDown_RelSample_2018");
 if(hh_mass_JESDown_RelSample_2018_branch) hh_mass_JESDown_RelSample_2018_branch->SetAddress(&hh_mass_JESDown_RelSample_2018_);//hh_mass_JESDown_RelSample_2018);


 deltaEta_j1j2_branch = tree->GetBranch("deltaEta_j1j2");
 if(deltaEta_j1j2_branch) deltaEta_j1j2_branch->SetAddress(&deltaEta_j1j2_);//
 deltaPhi_j1j2_branch = tree->GetBranch("deltaPhi_j1j2");
 if(deltaPhi_j1j2_branch) deltaPhi_j1j2_branch->SetAddress(&deltaPhi_j1j2_);//
 deltaR_j1j2_branch = tree->GetBranch("deltaR_j1j2");
 if(deltaR_j1j2_branch) deltaR_j1j2_branch->SetAddress(&deltaR_j1j2_);//
 ptj2_over_ptj1_branch = tree->GetBranch("ptj2_over_ptj1");
 if(ptj2_over_ptj1_branch) ptj2_over_ptj1_branch->SetAddress(&ptj2_over_ptj1_);//
 mj2_over_mj1_branch = tree->GetBranch("mj2_over_mj1");
 if(mj2_over_mj1_branch) mj2_over_mj1_branch->SetAddress(&mj2_over_mj1_);//
 isVBFtag_branch = tree->GetBranch("isVBFtag");
 if(isVBFtag_branch) isVBFtag_branch->SetAddress(&isVBFtag_);//
 isVBFtag_JERUp_branch = tree->GetBranch("isVBFtag_JERUp");
 if(isVBFtag_JERUp_branch) isVBFtag_JERUp_branch->SetAddress(&isVBFtag_JERUp_);//isVBFtag_JERUp);
 isVBFtag_JERDown_branch = tree->GetBranch("isVBFtag_JERDown");
 if(isVBFtag_JERDown_branch) isVBFtag_JERDown_branch->SetAddress(&isVBFtag_JERDown_);//isVBFtag_JERDown);
 isVBFtag_JESUp_branch = tree->GetBranch("isVBFtag_JESUp");
 if(isVBFtag_JESUp_branch) isVBFtag_JESUp_branch->SetAddress(&isVBFtag_JESUp_);//isVBFtag_JESUp);
 isVBFtag_JESDown_branch = tree->GetBranch("isVBFtag_JESDown");
 if(isVBFtag_JESDown_branch) isVBFtag_JESDown_branch->SetAddress(&isVBFtag_JESDown_);//isVBFtag_JESDown);
 isVBFtag_JESUp_Abs_branch = tree->GetBranch("isVBFtag_JESUp_Abs");
 if(isVBFtag_JESUp_Abs_branch) isVBFtag_JESUp_Abs_branch->SetAddress(&isVBFtag_JESUp_Abs_);//isVBFtag_JESUp_Abs);
 isVBFtag_JESDown_Abs_branch = tree->GetBranch("isVBFtag_JESDown_Abs");
 if(isVBFtag_JESDown_Abs_branch) isVBFtag_JESDown_Abs_branch->SetAddress(&isVBFtag_JESDown_Abs_);//isVBFtag_JESDown_Abs);
 isVBFtag_JESUp_Abs_2016_branch = tree->GetBranch("isVBFtag_JESUp_Abs_2016");
 if(isVBFtag_JESUp_Abs_2016_branch) isVBFtag_JESUp_Abs_2016_branch->SetAddress(&isVBFtag_JESUp_Abs_2016_);//isVBFtag_JESUp_Abs_2016);
 isVBFtag_JESDowb_Abs_2016_branch = tree->GetBranch("isVBFtag_JESDowb_Abs_2016");
 if(isVBFtag_JESDowb_Abs_2016_branch) isVBFtag_JESDowb_Abs_2016_branch->SetAddress(&isVBFtag_JESDowb_Abs_2016_);//isVBFtag_JESDowb_Abs_2016);
 isVBFtag_JESUp_Abs_2017_branch = tree->GetBranch("isVBFtag_JESUp_Abs_2017");
 if(isVBFtag_JESUp_Abs_2017_branch) isVBFtag_JESUp_Abs_2017_branch->SetAddress(&isVBFtag_JESUp_Abs_2017_);//isVBFtag_JESUp_Abs_2017);
 isVBFtag_JESDowb_Abs_2017_branch = tree->GetBranch("isVBFtag_JESDowb_Abs_2017");
 if(isVBFtag_JESDowb_Abs_2017_branch) isVBFtag_JESDowb_Abs_2017_branch->SetAddress(&isVBFtag_JESDowb_Abs_2017_);//isVBFtag_JESDowb_Abs_2017);
 isVBFtag_JESUp_Abs_2018_branch = tree->GetBranch("isVBFtag_JESUp_Abs_2018");
 if(isVBFtag_JESUp_Abs_2018_branch) isVBFtag_JESUp_Abs_2018_branch->SetAddress(&isVBFtag_JESUp_Abs_2018_);//isVBFtag_JESUp_Abs_2018);
 isVBFtag_JESDowb_Abs_2018_branch = tree->GetBranch("isVBFtag_JESDowb_Abs_2018");
 if(isVBFtag_JESDowb_Abs_2018_branch) isVBFtag_JESDowb_Abs_2018_branch->SetAddress(&isVBFtag_JESDowb_Abs_2018_);//isVBFtag_JESDowb_Abs_2018);
 isVBFtag_JESUp_BBEC1_branch = tree->GetBranch("isVBFtag_JESUp_BBEC1");
 if(isVBFtag_JESUp_BBEC1_branch) isVBFtag_JESUp_BBEC1_branch->SetAddress(&isVBFtag_JESUp_BBEC1_);//isVBFtag_JESUp_BBEC1);
 isVBFtag_JESDown_BBEC1_branch = tree->GetBranch("isVBFtag_JESDown_BBEC1");
 if(isVBFtag_JESDown_BBEC1_branch) isVBFtag_JESDown_BBEC1_branch->SetAddress(&isVBFtag_JESDown_BBEC1_);//isVBFtag_JESDown_BBEC1);
 isVBFtag_JESUp_BBEC1_2016_branch = tree->GetBranch("isVBFtag_JESUp_BBEC1_2016");
 if(isVBFtag_JESUp_BBEC1_2016_branch) isVBFtag_JESUp_BBEC1_2016_branch->SetAddress(&isVBFtag_JESUp_BBEC1_2016_);//isVBFtag_JESUp_BBEC1_2016);
 isVBFtag_JESDown_BBEC1_2016_branch = tree->GetBranch("isVBFtag_JESDown_BBEC1_2016");
 if(isVBFtag_JESDown_BBEC1_2016_branch) isVBFtag_JESDown_BBEC1_2016_branch->SetAddress(&isVBFtag_JESDown_BBEC1_2016_);//isVBFtag_JESDown_BBEC1_2016);
 isVBFtag_JESUp_BBEC1_2017_branch = tree->GetBranch("isVBFtag_JESUp_BBEC1_2017");
 if(isVBFtag_JESUp_BBEC1_2017_branch) isVBFtag_JESUp_BBEC1_2017_branch->SetAddress(&isVBFtag_JESUp_BBEC1_2017_);//isVBFtag_JESUp_BBEC1_2017);
 isVBFtag_JESDown_BBEC1_2017_branch = tree->GetBranch("isVBFtag_JESDown_BBEC1_2017");
 if(isVBFtag_JESDown_BBEC1_2017_branch) isVBFtag_JESDown_BBEC1_2017_branch->SetAddress(&isVBFtag_JESDown_BBEC1_2017_);//isVBFtag_JESDown_BBEC1_2017);
 isVBFtag_JESUp_BBEC1_2018_branch = tree->GetBranch("isVBFtag_JESUp_BBEC1_2018");
 if(isVBFtag_JESUp_BBEC1_2018_branch) isVBFtag_JESUp_BBEC1_2018_branch->SetAddress(&isVBFtag_JESUp_BBEC1_2018_);//isVBFtag_JESUp_BBEC1_2018);
 isVBFtag_JESDown_BBEC1_2018_branch = tree->GetBranch("isVBFtag_JESDown_BBEC1_2018");
 if(isVBFtag_JESDown_BBEC1_2018_branch) isVBFtag_JESDown_BBEC1_2018_branch->SetAddress(&isVBFtag_JESDown_BBEC1_2018_);//isVBFtag_JESDown_BBEC1_2018);
 isVBFtag_JESUp_EC2_branch = tree->GetBranch("isVBFtag_JESUp_EC2");
 if(isVBFtag_JESUp_EC2_branch) isVBFtag_JESUp_EC2_branch->SetAddress(&isVBFtag_JESUp_EC2_);//isVBFtag_JESUp_EC2);
 isVBFtag_JESDown_EC2_branch = tree->GetBranch("isVBFtag_JESDown_EC2");
 if(isVBFtag_JESDown_EC2_branch) isVBFtag_JESDown_EC2_branch->SetAddress(&isVBFtag_JESDown_EC2_);//isVBFtag_JESDown_EC2);
 isVBFtag_JESUp_EC2_2016_branch = tree->GetBranch("isVBFtag_JESUp_EC2_2016");
 if(isVBFtag_JESUp_EC2_2016_branch) isVBFtag_JESUp_EC2_2016_branch->SetAddress(&isVBFtag_JESUp_EC2_2016_);//isVBFtag_JESUp_EC2_2016);
 isVBFtag_JESDown_EC2_2016_branch = tree->GetBranch("isVBFtag_JESDown_EC2_2016");
 if(isVBFtag_JESDown_EC2_2016_branch) isVBFtag_JESDown_EC2_2016_branch->SetAddress(&isVBFtag_JESDown_EC2_2016_);//isVBFtag_JESDown_EC2_2016);
 isVBFtag_JESUp_EC2_2017_branch = tree->GetBranch("isVBFtag_JESUp_EC2_2017");
 if(isVBFtag_JESUp_EC2_2017_branch) isVBFtag_JESUp_EC2_2017_branch->SetAddress(&isVBFtag_JESUp_EC2_2017_);//isVBFtag_JESUp_EC2_2017);
 isVBFtag_JESDown_EC2_2017_branch = tree->GetBranch("isVBFtag_JESDown_EC2_2017");
 if(isVBFtag_JESDown_EC2_2017_branch) isVBFtag_JESDown_EC2_2017_branch->SetAddress(&isVBFtag_JESDown_EC2_2017_);//isVBFtag_JESDown_EC2_2017);
 isVBFtag_JESUp_EC2_2018_branch = tree->GetBranch("isVBFtag_JESUp_EC2_2018");
 if(isVBFtag_JESUp_EC2_2018_branch) isVBFtag_JESUp_EC2_2018_branch->SetAddress(&isVBFtag_JESUp_EC2_2018_);//isVBFtag_JESUp_EC2_2018);
 isVBFtag_JESDown_EC2_2018_branch = tree->GetBranch("isVBFtag_JESDown_EC2_2018");
 if(isVBFtag_JESDown_EC2_2018_branch) isVBFtag_JESDown_EC2_2018_branch->SetAddress(&isVBFtag_JESDown_EC2_2018_);//isVBFtag_JESDown_EC2_2018);
 isVBFtag_JESUp_FlavQCD_branch = tree->GetBranch("isVBFtag_JESUp_FlavQCD");
 if(isVBFtag_JESUp_FlavQCD_branch) isVBFtag_JESUp_FlavQCD_branch->SetAddress(&isVBFtag_JESUp_FlavQCD_);//isVBFtag_JESUp_FlavQCD);
 isVBFtag_JESDown_FlavQCD_branch = tree->GetBranch("isVBFtag_JESDown_FlavQCD");
 if(isVBFtag_JESDown_FlavQCD_branch) isVBFtag_JESDown_FlavQCD_branch->SetAddress(&isVBFtag_JESDown_FlavQCD_);//isVBFtag_JESDown_FlavQCD);
 isVBFtag_JESUp_HF_branch = tree->GetBranch("isVBFtag_JESUp_HF");
 if(isVBFtag_JESUp_HF_branch) isVBFtag_JESUp_HF_branch->SetAddress(&isVBFtag_JESUp_HF_);//isVBFtag_JESUp_HF);
 isVBFtag_JESDown_HF_branch = tree->GetBranch("isVBFtag_JESDown_HF");
 if(isVBFtag_JESDown_HF_branch) isVBFtag_JESDown_HF_branch->SetAddress(&isVBFtag_JESDown_HF_);//isVBFtag_JESDown_HF);
 isVBFtag_JESUp_HF_2016_branch = tree->GetBranch("isVBFtag_JESUp_HF_2016");
 if(isVBFtag_JESUp_HF_2016_branch) isVBFtag_JESUp_HF_2016_branch->SetAddress(&isVBFtag_JESUp_HF_2016_);//isVBFtag_JESUp_HF_2016);
 isVBFtag_JESDown_HF_2016_branch = tree->GetBranch("isVBFtag_JESDown_HF_2016");
 if(isVBFtag_JESDown_HF_2016_branch) isVBFtag_JESDown_HF_2016_branch->SetAddress(&isVBFtag_JESDown_HF_2016_);//isVBFtag_JESDown_HF_2016);
 isVBFtag_JESUp_HF_2017_branch = tree->GetBranch("isVBFtag_JESUp_HF_2017");
 if(isVBFtag_JESUp_HF_2017_branch) isVBFtag_JESUp_HF_2017_branch->SetAddress(&isVBFtag_JESUp_HF_2017_);//isVBFtag_JESUp_HF_2017);
 isVBFtag_JESDown_HF_2017_branch = tree->GetBranch("isVBFtag_JESDown_HF_2017");
 if(isVBFtag_JESDown_HF_2017_branch) isVBFtag_JESDown_HF_2017_branch->SetAddress(&isVBFtag_JESDown_HF_2017_);//isVBFtag_JESDown_HF_2017);
 isVBFtag_JESUp_HF_2018_branch = tree->GetBranch("isVBFtag_JESUp_HF_2018");
 if(isVBFtag_JESUp_HF_2018_branch) isVBFtag_JESUp_HF_2018_branch->SetAddress(&isVBFtag_JESUp_HF_2018_);//isVBFtag_JESUp_HF_2018);
 isVBFtag_JESDown_HF_2018_branch = tree->GetBranch("isVBFtag_JESDown_HF_2018");
 if(isVBFtag_JESDown_HF_2018_branch) isVBFtag_JESDown_HF_2018_branch->SetAddress(&isVBFtag_JESDown_HF_2018_);//isVBFtag_JESDown_HF_2018);
 isVBFtag_JESUp_RelBal_branch = tree->GetBranch("isVBFtag_JESUp_RelBal");
 if(isVBFtag_JESUp_RelBal_branch) isVBFtag_JESUp_RelBal_branch->SetAddress(&isVBFtag_JESUp_RelBal_);//isVBFtag_JESUp_RelBal);
 isVBFtag_JESDown_RelBal_branch = tree->GetBranch("isVBFtag_JESDown_RelBal");
 if(isVBFtag_JESDown_RelBal_branch) isVBFtag_JESDown_RelBal_branch->SetAddress(&isVBFtag_JESDown_RelBal_);//isVBFtag_JESDown_RelBal);
 isVBFtag_JESUp_RelSample_2016_branch = tree->GetBranch("isVBFtag_JESUp_RelSample_2016");
 if(isVBFtag_JESUp_RelSample_2016_branch) isVBFtag_JESUp_RelSample_2016_branch->SetAddress(&isVBFtag_JESUp_RelSample_2016_);//isVBFtag_JESUp_RelSample_2016);
 isVBFtag_JESDown_RelSample_2016_branch = tree->GetBranch("isVBFtag_JESDown_RelSample_2016");
 if(isVBFtag_JESDown_RelSample_2016_branch) isVBFtag_JESDown_RelSample_2016_branch->SetAddress(&isVBFtag_JESDown_RelSample_2016_);//isVBFtag_JESDown_RelSample_2016);
 isVBFtag_JESUp_RelSample_2017_branch = tree->GetBranch("isVBFtag_JESUp_RelSample_2017");
 if(isVBFtag_JESUp_RelSample_2017_branch) isVBFtag_JESUp_RelSample_2017_branch->SetAddress(&isVBFtag_JESUp_RelSample_2017_);//isVBFtag_JESUp_RelSample_2017);
 isVBFtag_JESDown_RelSample_2017_branch = tree->GetBranch("isVBFtag_JESDown_RelSample_2017");
 if(isVBFtag_JESDown_RelSample_2017_branch) isVBFtag_JESDown_RelSample_2017_branch->SetAddress(&isVBFtag_JESDown_RelSample_2017_);//isVBFtag_JESDown_RelSample_2017);
 isVBFtag_JESUp_RelSample_2018_branch = tree->GetBranch("isVBFtag_JESUp_RelSample_2018");
 if(isVBFtag_JESUp_RelSample_2018_branch) isVBFtag_JESUp_RelSample_2018_branch->SetAddress(&isVBFtag_JESUp_RelSample_2018_);//isVBFtag_JESUp_RelSample_2018);
 isVBFtag_JESDown_RelSample_2018_branch = tree->GetBranch("isVBFtag_JESDown_RelSample_2018");
 if(isVBFtag_JESDown_RelSample_2018_branch) isVBFtag_JESDown_RelSample_2018_branch->SetAddress(&isVBFtag_JESDown_RelSample_2018_);//isVBFtag_JESDown_RelSample_2018);




 dijetmass_branch = tree->GetBranch("dijetmass");
 if(dijetmass_branch) dijetmass_branch->SetAddress(&dijetmass_);//
 vbfjet1Pt_branch = tree->GetBranch("vbfjet1Pt");
 if(vbfjet1Pt_branch) vbfjet1Pt_branch->SetAddress(&vbfjet1Pt_);//
 vbfjet1Eta_branch = tree->GetBranch("vbfjet1Eta");
 if(vbfjet1Eta_branch) vbfjet1Eta_branch->SetAddress(&vbfjet1Eta_);//
 vbfjet1Phi_branch = tree->GetBranch("vbfjet1Phi");
 if(vbfjet1Phi_branch) vbfjet1Phi_branch->SetAddress(&vbfjet1Phi_);//
 vbfjet1Mass_branch = tree->GetBranch("vbfjet1Mass");
 if(vbfjet1Mass_branch) vbfjet1Mass_branch->SetAddress(&vbfjet1Mass_);//
 vbffatJet1Pt_branch = tree->GetBranch("vbffatJet1Pt");
 if(vbffatJet1Pt_branch) vbffatJet1Pt_branch->SetAddress(&vbffatJet1Pt_);//
 vbffatJet1Eta_branch = tree->GetBranch("vbffatJet1Eta");
 if(vbffatJet1Eta_branch) vbffatJet1Eta_branch->SetAddress(&vbffatJet1Eta_);//
 vbffatJet1Phi_branch = tree->GetBranch("vbffatJet1Phi");
 if(vbffatJet1Phi_branch) vbffatJet1Phi_branch->SetAddress(&vbffatJet1Phi_);//
 vbffatJet1PNetXbb_branch = tree->GetBranch("vbffatJet1PNetXbb");
 if(vbffatJet1PNetXbb_branch) vbffatJet1PNetXbb_branch->SetAddress(&vbffatJet1PNetXbb_);//
 vbfjet2Pt_branch = tree->GetBranch("vbfjet2Pt");
 if(vbfjet2Pt_branch) vbfjet2Pt_branch->SetAddress(&vbfjet2Pt_);//
 vbfjet2Eta_branch = tree->GetBranch("vbfjet2Eta");
 if(vbfjet2Eta_branch) vbfjet2Eta_branch->SetAddress(&vbfjet2Eta_);//
 vbfjet2Phi_branch = tree->GetBranch("vbfjet2Phi");
 if(vbfjet2Phi_branch) vbfjet2Phi_branch->SetAddress(&vbfjet2Phi_);//
 vbfjet2Mass_branch = tree->GetBranch("vbfjet2Mass");
 if(vbfjet2Mass_branch) vbfjet2Mass_branch->SetAddress(&vbfjet2Mass_);//
 vbffatJet2Pt_branch = tree->GetBranch("vbffatJet2Pt");
 if(vbffatJet2Pt_branch) vbffatJet2Pt_branch->SetAddress(&vbffatJet2Pt_);//
 vbffatJet2Eta_branch = tree->GetBranch("vbffatJet2Eta");
 if(vbffatJet2Eta_branch) vbffatJet2Eta_branch->SetAddress(&vbffatJet2Eta_);//
 vbffatJet2Phi_branch = tree->GetBranch("vbffatJet2Phi");
 if(vbffatJet2Phi_branch) vbffatJet2Phi_branch->SetAddress(&vbffatJet2Phi_);//
 vbffatJet2PNetXbb_branch = tree->GetBranch("vbffatJet2PNetXbb");
 if(vbffatJet2PNetXbb_branch) vbffatJet2PNetXbb_branch->SetAddress(&vbffatJet2PNetXbb_);//
 jet1Pt_branch = tree->GetBranch("jet1Pt");
 if(jet1Pt_branch) jet1Pt_branch->SetAddress(&jet1Pt_);//
 jet1Eta_branch = tree->GetBranch("jet1Eta");
 if(jet1Eta_branch) jet1Eta_branch->SetAddress(&jet1Eta_);//
 jet1Phi_branch = tree->GetBranch("jet1Phi");
 if(jet1Phi_branch) jet1Phi_branch->SetAddress(&jet1Phi_);//
 jet2Pt_branch = tree->GetBranch("jet2Pt");
 if(jet2Pt_branch) jet2Pt_branch->SetAddress(&jet2Pt_);//
 jet2Eta_branch = tree->GetBranch("jet2Eta");
 if(jet2Eta_branch) jet2Eta_branch->SetAddress(&jet2Eta_);//
 jet2Phi_branch = tree->GetBranch("jet2Phi");
 if(jet2Phi_branch) jet2Phi_branch->SetAddress(&jet2Phi_);//
 jet3Pt_branch = tree->GetBranch("jet3Pt");
 if(jet3Pt_branch) jet3Pt_branch->SetAddress(&jet3Pt_);//
 jet3Eta_branch = tree->GetBranch("jet3Eta");
 if(jet3Eta_branch) jet3Eta_branch->SetAddress(&jet3Eta_);//
 jet3Phi_branch = tree->GetBranch("jet3Phi");
 if(jet3Phi_branch) jet3Phi_branch->SetAddress(&jet3Phi_);//
 jet4Pt_branch = tree->GetBranch("jet4Pt");
 if(jet4Pt_branch) jet4Pt_branch->SetAddress(&jet4Pt_);//
 jet4Eta_branch = tree->GetBranch("jet4Eta");
 if(jet4Eta_branch) jet4Eta_branch->SetAddress(&jet4Eta_);//
 jet4Phi_branch = tree->GetBranch("jet4Phi");
 if(jet4Phi_branch) jet4Phi_branch->SetAddress(&jet4Phi_);//
 lep1Pt_branch = tree->GetBranch("lep1Pt");
 if(lep1Pt_branch) lep1Pt_branch->SetAddress(&lep1Pt_);//
 lep1Eta_branch = tree->GetBranch("lep1Eta");
 if(lep1Eta_branch) lep1Eta_branch->SetAddress(&lep1Eta_);//
 lep1Phi_branch = tree->GetBranch("lep1Phi");
 if(lep1Phi_branch) lep1Phi_branch->SetAddress(&lep1Phi_);//
 lep1Id_branch = tree->GetBranch("lep1Id");
 if(lep1Id_branch) lep1Id_branch->SetAddress(&lep1Id_);//
 lep2Pt_branch = tree->GetBranch("lep2Pt");
 if(lep2Pt_branch) lep2Pt_branch->SetAddress(&lep2Pt_);//
 lep2Eta_branch = tree->GetBranch("lep2Eta");
 if(lep2Eta_branch) lep2Eta_branch->SetAddress(&lep2Eta_);//
 lep2Phi_branch = tree->GetBranch("lep2Phi");
 if(lep2Phi_branch) lep2Phi_branch->SetAddress(&lep2Phi_);//
 lep2Id_branch = tree->GetBranch("lep2Id");
 if(lep2Id_branch) lep2Id_branch->SetAddress(&lep2Id_);//
 genHiggs1Pt_branch = tree->GetBranch("genHiggs1Pt");
 if(genHiggs1Pt_branch) genHiggs1Pt_branch->SetAddress(&genHiggs1Pt_);//
 genHiggs1Eta_branch = tree->GetBranch("genHiggs1Eta");
 if(genHiggs1Eta_branch) genHiggs1Eta_branch->SetAddress(&genHiggs1Eta_);//
 genHiggs1Phi_branch = tree->GetBranch("genHiggs1Phi");
 if(genHiggs1Phi_branch) genHiggs1Phi_branch->SetAddress(&genHiggs1Phi_);//
 genHiggs2Pt_branch = tree->GetBranch("genHiggs2Pt");
 if(genHiggs2Pt_branch) genHiggs2Pt_branch->SetAddress(&genHiggs2Pt_);//
 genHiggs2Eta_branch = tree->GetBranch("genHiggs2Eta");
 if(genHiggs2Eta_branch) genHiggs2Eta_branch->SetAddress(&genHiggs2Eta_);//
 genHiggs2Phi_branch = tree->GetBranch("genHiggs2Phi");
 if(genHiggs2Phi_branch) genHiggs2Phi_branch->SetAddress(&genHiggs2Phi_);//
 puWeight_branch = tree->GetBranch("puWeight");
 if(puWeight_branch) puWeight_branch->SetAddress(&puWeight_);//
 puWeightUp_branch = tree->GetBranch("puWeightUp");
 if(puWeightUp_branch) puWeightUp_branch->SetAddress(&puWeightUp_);//
 puWeightDown_branch = tree->GetBranch("puWeightDown");
 if(puWeightDown_branch) puWeightDown_branch->SetAddress(&puWeightDown_);//
 xsecWeight_branch = tree->GetBranch("xsecWeight");
 if(xsecWeight_branch) xsecWeight_branch->SetAddress(&xsecWeight_);//
 LHEScaleWeightNormNew_branch = tree->GetBranch("LHEScaleWeightNormNew");
 if(LHEScaleWeightNormNew_branch) LHEScaleWeightNormNew_branch->SetAddress(LHEScaleWeightNormNew_);//
 LHEPdfWeightNorm_branch = tree->GetBranch("LHEPdfWeightNorm");
 if(LHEPdfWeightNorm_branch) LHEPdfWeightNorm_branch->SetAddress(LHEPdfWeightNorm_);//
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_);//
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_);//
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_);//
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_);//
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_);//
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_);//
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_);//


 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017);
  disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017);
  disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016);
disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017);
disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018);
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_branch = tree->GetBranch("disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018");
 if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_branch) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_branch->SetAddress(&disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_);//disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018);




}

void hhtree::GetEntry(unsigned int idx)
{

 index = idx;
 run_isLoaded = false;
 luminosityBlock_isLoaded = false;
 event_isLoaded = false;
 genWeight_isLoaded = false;
 nLHEPdfWeight_isLoaded = false;
 LHEPdfWeight_isLoaded = false;
 nLHEScaleWeight_isLoaded = false;
 LHEScaleWeight_isLoaded = false;
 nPSWeight_isLoaded = false;
 PSWeight_isLoaded = false;
 HLT_Ele27_WPTight_Gsf_isLoaded = false;
 HLT_Ele28_WPTight_Gsf_isLoaded = false;
 HLT_Ele30_WPTight_Gsf_isLoaded = false;
 HLT_Ele32_WPTight_Gsf_isLoaded = false;
 HLT_Ele32_WPTight_Gsf_L1DoubleEG_isLoaded = false;
 HLT_Ele35_WPTight_Gsf_isLoaded = false;
 HLT_Ele38_WPTight_Gsf_isLoaded = false;
 HLT_Ele40_WPTight_Gsf_isLoaded = false;
 HLT_IsoMu20_isLoaded = false;
 HLT_IsoMu24_isLoaded = false;
 HLT_IsoMu24_eta2p1_isLoaded = false;
 HLT_IsoMu27_isLoaded = false;
 HLT_IsoMu30_isLoaded = false;
 HLT_Mu50_isLoaded = false;
 HLT_Mu55_isLoaded = false;
 HLT_Photon175_isLoaded = false;
 HLT_PFHT780_isLoaded = false;
 HLT_PFHT890_isLoaded = false;
 HLT_PFHT1050_isLoaded = false;
 HLT_AK8PFJet360_TrimMass30_isLoaded = false;
 HLT_AK8PFJet380_TrimMass30_isLoaded = false;
 HLT_AK8PFJet400_TrimMass30_isLoaded = false;
 HLT_AK8PFJet420_TrimMass30_isLoaded = false;
 HLT_AK8PFHT750_TrimMass50_isLoaded = false;
 HLT_AK8PFHT800_TrimMass50_isLoaded = false;
 HLT_AK8PFHT850_TrimMass50_isLoaded = false;
 HLT_AK8PFHT900_TrimMass50_isLoaded = false;
 HLT_PFJet450_isLoaded = false;
 HLT_PFJet500_isLoaded = false;
 HLT_PFJet550_isLoaded = false;
 HLT_AK8PFJet400_isLoaded = false;
 HLT_AK8PFJet450_isLoaded = false;
 HLT_AK8PFJet500_isLoaded = false;
 HLT_AK8PFJet550_isLoaded = false;
 HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_isLoaded = false;
 HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_isLoaded = false;
 HLT_AK8PFJet330_PFAK8BTagCSV_p17_isLoaded = false;
 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_isLoaded = false;
 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_isLoaded = false;
 HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_isLoaded = false;
 HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_isLoaded = false;
 HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_isLoaded = false;
 HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_isLoaded = false;
 HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_isLoaded = false;
 HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_isLoaded = false;
 HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_isLoaded = false;
 weight_isLoaded = false;
 met_isLoaded = false;
 metphi_isLoaded = false;
 ht_isLoaded = false;
 passmetfilters_isLoaded = false;
 l1PreFiringWeight_isLoaded = false;
 l1PreFiringWeightUp_isLoaded = false;
 l1PreFiringWeightDown_isLoaded = false;
 triggerEffWeight_isLoaded = false;
 triggerEff3DWeight_isLoaded = false;
 triggerEffMCWeight_isLoaded = false;
 triggerEffMC3DWeight_isLoaded = false;
 fatJet1Pt_isLoaded = false;
 fatJet1Eta_isLoaded = false;
 fatJet1Phi_isLoaded = false;
 fatJet1Mass_isLoaded = false;
 fatJet1MassSD_isLoaded = false;
 fatJet1MassSD_UnCorrected_isLoaded = false;
 fatJet1MassRegressed_isLoaded = false;
 fatJet1MassRegressed_UnCorrected_isLoaded = false;
 fatJet1PNetXbb_isLoaded = false;
 fatJet1PNetQCDb_isLoaded = false;
 fatJet1PNetQCDbb_isLoaded = false;
 fatJet1PNetQCDc_isLoaded = false;
 fatJet1PNetQCDcc_isLoaded = false;
 fatJet1PNetQCDothers_isLoaded = false;
 fatJet1Tau3OverTau2_isLoaded = false;
 fatJet1GenMatchIndex_isLoaded = false;
 fatJet1HasMuon_isLoaded = false;
 fatJet1HasElectron_isLoaded = false;
 fatJet1HasBJetCSVLoose_isLoaded = false;
 fatJet1HasBJetCSVMedium_isLoaded = false;
 fatJet1HasBJetCSVTight_isLoaded = false;
 fatJet1OppositeHemisphereHasBJet_isLoaded = false;
 fatJet1PtOverMHH_isLoaded = false;
 fatJet1PtOverMSD_isLoaded = false;
 fatJet1PtOverMRegressed_isLoaded = false;
 fatJet1MassSD_JMS_Down_isLoaded = false;
 fatJet1MassSD_JMS_Up_isLoaded = false;
 fatJet1MassSD_JMR_Down_isLoaded = false;
 fatJet1MassSD_JMR_Up_isLoaded = false;
 fatJet1MassRegressed_JMS_Down_isLoaded = false;
 fatJet1MassRegressed_JMS_Up_isLoaded = false;
 fatJet1MassRegressed_JMR_Down_isLoaded = false;
 fatJet1MassRegressed_JMR_Up_isLoaded = false;
 fatJet1PtOverMHH_JMS_Down_isLoaded = false;
 fatJet1PtOverMHH_JMS_Up_isLoaded = false;
 fatJet1PtOverMHH_JMR_Down_isLoaded = false;
 fatJet1PtOverMHH_JMR_Up_isLoaded = false;
 fatJet1Pt_JERUp_isLoaded = false;
 fatJet1PtOverMHH_JERUp_isLoaded = false;
 fatJet1Pt_JERDown_isLoaded = false;
 fatJet1PtOverMHH_JERDown_isLoaded = false;
 fatJet1Pt_JESUp_isLoaded = false;
 fatJet1PtOverMHH_JESUp_isLoaded = false;
 fatJet1Pt_JESDown_isLoaded = false;
 fatJet1PtOverMHH_JESDown_isLoaded = false;
 fatJet1Pt_JESUp_Abs_isLoaded = false;
 fatJet1PtOverMHH_JESUp_Abs_isLoaded = false;
 fatJet1Pt_JESDown_Abs_isLoaded = false;
 fatJet1PtOverMHH_JESDown_Abs_isLoaded = false;
 fatJet1Pt_JESUp_Abs_2016_isLoaded = false;
 fatJet1PtOverMHH_JESUp_Abs_2016_isLoaded = false;
 fatJet1Pt_JESDowb_Abs_2016_isLoaded = false;
 fatJet1PtOverMHH_JESDowb_Abs_2016_isLoaded = false;
 fatJet1Pt_JESUp_Abs_2017_isLoaded = false;
 fatJet1PtOverMHH_JESUp_Abs_2017_isLoaded = false;
 fatJet1Pt_JESDowb_Abs_2017_isLoaded = false;
 fatJet1PtOverMHH_JESDowb_Abs_2017_isLoaded = false;
 fatJet1Pt_JESUp_Abs_2018_isLoaded = false;
 fatJet1PtOverMHH_JESUp_Abs_2018_isLoaded = false;
 fatJet1Pt_JESDowb_Abs_2018_isLoaded = false;
 fatJet1PtOverMHH_JESDowb_Abs_2018_isLoaded = false;
 fatJet1Pt_JESUp_BBEC1_isLoaded = false;
 fatJet1PtOverMHH_JESUp_BBEC1_isLoaded = false;
 fatJet1Pt_JESDown_BBEC1_isLoaded = false;
 fatJet1PtOverMHH_JESDown_BBEC1_isLoaded = false;
 fatJet1Pt_JESUp_BBEC1_2016_isLoaded = false;
 fatJet1PtOverMHH_JESUp_BBEC1_2017_isLoaded = false;
 fatJet1Pt_JESDown_BBEC1_2017_isLoaded = false;
 fatJet1PtOverMHH_JESDown_BBEC1_2017_isLoaded = false;
 fatJet1Pt_JESUp_BBEC1_2017_isLoaded = false;
 fatJet1PtOverMHH_JESUp_BBEC1_2017_isLoaded = false;
 fatJet1Pt_JESDown_BBEC1_2017_isLoaded = false;
 fatJet1PtOverMHH_JESDown_BBEC1_2017_isLoaded = false;
 fatJet1PtOverMHH_JESUp_BBEC1_2018_isLoaded = false;
 fatJet1Pt_JESDown_BBEC1_2018_isLoaded = false;
 fatJet1PtOverMHH_JESDown_BBEC1_2018_isLoaded = false;
 fatJet1Pt_JESUp_BBEC1_2018_isLoaded = false;
 fatJet1PtOverMHH_JESUp_BBEC1_2018_isLoaded = false;
 fatJet1Pt_JESDown_BBEC1_2018_isLoaded = false;
 fatJet1PtOverMHH_JESDown_BBEC1_2018_isLoaded = false;
 fatJet1Pt_JESUp_EC2_isLoaded = false;
 
 fatJet1PtOverMHH_JESUp_EC2_isLoaded = false;
 fatJet1Pt_JESDown_EC2_isLoaded = false;
 fatJet1PtOverMHH_JESDown_EC2_isLoaded = false;
 fatJet1Pt_JESUp_EC2_2016_isLoaded = false;
 fatJet1PtOverMHH_JESUp_EC2_2016_isLoaded = false;
 fatJet1Pt_JESDown_EC2_2016_isLoaded = false;
 fatJet1PtOverMHH_JESDown_EC2_2016_isLoaded = false;
 fatJet1Pt_JESUp_EC2_2017_isLoaded = false;
 fatJet1PtOverMHH_JESUp_EC2_2017_isLoaded = false;
 fatJet1Pt_JESDown_EC2_2017_isLoaded = false;
 fatJet1PtOverMHH_JESDown_EC2_2017_isLoaded = false;
 fatJet1Pt_JESUp_EC2_2018_isLoaded = false;
 fatJet1PtOverMHH_JESUp_EC2_2018_isLoaded = false;
 fatJet1Pt_JESDown_EC2_2018_isLoaded = false;
 fatJet1PtOverMHH_JESDown_EC2_2018_isLoaded = false;
 fatJet1Pt_JESUp_FlavQCD_isLoaded = false;
 fatJet1PtOverMHH_JESUp_FlavQCD_isLoaded = false;
 fatJet1Pt_JESDown_FlavQCD_isLoaded = false;
 fatJet1PtOverMHH_JESDown_FlavQCD_isLoaded = false;
 
 fatJet1Pt_JESUp_HF_isLoaded = false;
 fatJet1PtOverMHH_JESUp_HF_isLoaded = false;
 fatJet1Pt_JESDown_HF_isLoaded = false;
 fatJet1PtOverMHH_JESDown_HF_isLoaded = false;
 fatJet1Pt_JESUp_HF_2016_isLoaded = false;
 fatJet1PtOverMHH_JESUp_HF_2016_isLoaded = false;
 fatJet1Pt_JESDown_HF_2016_isLoaded = false;
 fatJet1PtOverMHH_JESDown_HF_2016_isLoaded = false;
 fatJet1Pt_JESUp_HF_2017_isLoaded = false;
 fatJet1PtOverMHH_JESUp_HF_2017_isLoaded = false;
 fatJet1Pt_JESDown_HF_2017_isLoaded = false;
 fatJet1PtOverMHH_JESDown_HF_2017_isLoaded = false;
 fatJet1Pt_JESUp_HF_2018_isLoaded = false;
 fatJet1PtOverMHH_JESUp_HF_2018_isLoaded = false;
 fatJet1Pt_JESDown_HF_2018_isLoaded = false;
 fatJet1PtOverMHH_JESDown_HF_2018_isLoaded = false;
 fatJet1Pt_JESUp_RelBal_isLoaded = false;
 fatJet1PtOverMHH_JESUp_RelBal_isLoaded = false;
 fatJet1Pt_JESDown_RelBal_isLoaded = false;
 fatJet1PtOverMHH_JESDown_RelBal_isLoaded = false;
 fatJet1Pt_JESUp_RelSample_2016_isLoaded = false;
 fatJet1PtOverMHH_JESUp_RelSample_2016_isLoaded = false;
 fatJet1Pt_JESDown_RelSample_2016_isLoaded = false;
 fatJet1PtOverMHH_JESDown_RelSample_2016_isLoaded = false;
 fatJet1Pt_JESUp_RelSample_2017_isLoaded = false;
 fatJet1PtOverMHH_JESUp_RelSample_2017_isLoaded = false;
 fatJet1Pt_JESDown_RelSample_2017_isLoaded = false;
 fatJet1PtOverMHH_JESDown_RelSample_2017_isLoaded = false;
 fatJet1Pt_JESUp_RelSample_2018_isLoaded = false;
 fatJet1PtOverMHH_JESUp_RelSample_2018_isLoaded = false;
 fatJet1Pt_JESDown_RelSample_2018_isLoaded = false;
 fatJet1PtOverMHH_JESDown_RelSample_2018_isLoaded = false;
 fatJet2Pt_isLoaded = false;
 fatJet2Eta_isLoaded = false;
 fatJet2Phi_isLoaded = false;
 fatJet2Mass_isLoaded = false;
 fatJet2MassSD_isLoaded = false;
 fatJet2MassSD_UnCorrected_isLoaded = false;
 fatJet2MassRegressed_isLoaded = false;
 fatJet2MassRegressed_UnCorrected_isLoaded = false;
 fatJet2PNetXbb_isLoaded = false;
 fatJet2PNetQCDb_isLoaded = false;
 fatJet2PNetQCDbb_isLoaded = false;
 fatJet2PNetQCDc_isLoaded = false;
 fatJet2PNetQCDcc_isLoaded = false;
 fatJet2PNetQCDothers_isLoaded = false;
 fatJet2Tau3OverTau2_isLoaded = false;
 fatJet2GenMatchIndex_isLoaded = false;
 fatJet2HasMuon_isLoaded = false;
 fatJet2HasElectron_isLoaded = false;
 fatJet2HasBJetCSVLoose_isLoaded = false;
 fatJet2HasBJetCSVMedium_isLoaded = false;
 fatJet2HasBJetCSVTight_isLoaded = false;
 fatJet2OppositeHemisphereHasBJet_isLoaded = false;
 fatJet2PtOverMHH_isLoaded = false;
 fatJet2PtOverMSD_isLoaded = false;
 fatJet2PtOverMRegressed_isLoaded = false;
 fatJet2MassSD_JMS_Down_isLoaded = false;
 fatJet2MassSD_JMS_Up_isLoaded = false;
 fatJet2MassSD_JMR_Down_isLoaded = false;
 fatJet2MassSD_JMR_Up_isLoaded = false;
 fatJet2MassRegressed_JMS_Down_isLoaded = false;
 fatJet2MassRegressed_JMS_Up_isLoaded = false;
 fatJet2MassRegressed_JMR_Down_isLoaded = false;
 fatJet2MassRegressed_JMR_Up_isLoaded = false;
 fatJet2PtOverMHH_JMS_Down_isLoaded = false;
 fatJet2PtOverMHH_JMS_Up_isLoaded = false;
 fatJet2PtOverMHH_JMR_Down_isLoaded = false;
 fatJet2PtOverMHH_JMR_Up_isLoaded = false;
 fatJet2Pt_JERUp_isLoaded = false;
 fatJet2PtOverMHH_JERUp_isLoaded = false;
 fatJet2Pt_JERDown_isLoaded = false;
 fatJet2PtOverMHH_JERDown_isLoaded = false;
 fatJet2Pt_JESUp_isLoaded = false;
 fatJet2PtOverMHH_JESUp_isLoaded = false;
 fatJet2Pt_JESDown_isLoaded = false;
 fatJet2PtOverMHH_JESDown_isLoaded = false;

 fatJet2Pt_JESUp_Abs_isLoaded = false;
 fatJet2PtOverMHH_JESUp_Abs_isLoaded = false;
 fatJet2Pt_JESDown_Abs_isLoaded = false;
 fatJet2PtOverMHH_JESDown_Abs_isLoaded = false;
 fatJet2Pt_JESUp_Abs_2016_isLoaded = false;
 fatJet2PtOverMHH_JESUp_Abs_2016_isLoaded = false;
 fatJet2Pt_JESDowb_Abs_2016_isLoaded = false;
 fatJet2PtOverMHH_JESDowb_Abs_2016_isLoaded = false;
 fatJet2Pt_JESUp_Abs_2017_isLoaded = false;
 fatJet2PtOverMHH_JESUp_Abs_2017_isLoaded = false;
 fatJet2Pt_JESDowb_Abs_2017_isLoaded = false;
 fatJet2PtOverMHH_JESDowb_Abs_2017_isLoaded = false;
 fatJet2Pt_JESUp_Abs_2018_isLoaded = false;
 fatJet2PtOverMHH_JESUp_Abs_2018_isLoaded = false;
 fatJet2Pt_JESDowb_Abs_2018_isLoaded = false;
 fatJet2PtOverMHH_JESDowb_Abs_2018_isLoaded = false;
 fatJet2Pt_JESUp_BBEC1_isLoaded = false;
 fatJet2PtOverMHH_JESUp_BBEC1_isLoaded = false;
 fatJet2Pt_JESDown_BBEC1_isLoaded = false;
 fatJet2PtOverMHH_JESDown_BBEC1_isLoaded = false;
 fatJet2Pt_JESUp_BBEC1_2016_isLoaded = false;
 fatJet2PtOverMHH_JESUp_BBEC1_2016_isLoaded = false;
 fatJet2Pt_JESDown_BBEC1_2016_isLoaded = false;
 fatJet2PtOverMHH_JESDown_BBEC1_2016_isLoaded = false;
 fatJet2Pt_JESUp_BBEC1_2017_isLoaded = false;
 fatJet2PtOverMHH_JESUp_BBEC1_2017_isLoaded = false;
 fatJet2Pt_JESDown_BBEC1_2017_isLoaded = false;
 fatJet2PtOverMHH_JESDown_BBEC1_2017_isLoaded = false;
 fatJet2Pt_JESUp_BBEC1_2018_isLoaded = false;
 fatJet2PtOverMHH_JESUp_BBEC1_2018_isLoaded = false;
 fatJet2Pt_JESDown_BBEC1_2018_isLoaded = false;
 fatJet2PtOverMHH_JESDown_BBEC1_2018_isLoaded = false;
 fatJet2Pt_JESUp_EC2_isLoaded = false;
 fatJet2PtOverMHH_JESUp_EC2_isLoaded = false;
 fatJet2Pt_JESDown_EC2_isLoaded = false;
 fatJet2PtOverMHH_JESDown_EC2_isLoaded = false;
 fatJet2Pt_JESUp_EC2_2016_isLoaded = false;
 fatJet2PtOverMHH_JESUp_EC2_2016_isLoaded = false;
 fatJet2Pt_JESDown_EC2_2016_isLoaded = false;
 fatJet2PtOverMHH_JESDown_EC2_2016_isLoaded = false;
 fatJet2Pt_JESUp_EC2_2017_isLoaded = false;
 fatJet2PtOverMHH_JESUp_EC2_2017_isLoaded = false;
 fatJet2Pt_JESDown_EC2_2017_isLoaded = false;
 fatJet2PtOverMHH_JESDown_EC2_2017_isLoaded = false;
 fatJet2Pt_JESUp_EC2_2018_isLoaded = false;
 fatJet2PtOverMHH_JESUp_EC2_2018_isLoaded = false;
 fatJet2Pt_JESDown_EC2_2018_isLoaded = false;
 fatJet2PtOverMHH_JESDown_EC2_2018_isLoaded = false;
 fatJet2Pt_JESUp_FlavQCD_isLoaded = false;
 fatJet2PtOverMHH_JESUp_FlavQCD_isLoaded = false;
 fatJet2Pt_JESDown_FlavQCD_isLoaded = false;
 fatJet2PtOverMHH_JESDown_FlavQCD_isLoaded = false;
 fatJet2Pt_JESUp_HF_isLoaded = false;
 fatJet2PtOverMHH_JESUp_HF_isLoaded = false;
 fatJet2Pt_JESDown_HF_isLoaded = false;
 fatJet2PtOverMHH_JESDown_HF_isLoaded = false;
 fatJet2Pt_JESUp_HF_2016_isLoaded = false;
 fatJet2PtOverMHH_JESUp_HF_2016_isLoaded = false;
 fatJet2Pt_JESDown_HF_2016_isLoaded = false;
 fatJet2PtOverMHH_JESDown_HF_2016_isLoaded = false;
 fatJet2Pt_JESUp_HF_2017_isLoaded = false;
 fatJet2PtOverMHH_JESUp_HF_2017_isLoaded = false;
 fatJet2Pt_JESDown_HF_2017_isLoaded = false;
 fatJet2PtOverMHH_JESDown_HF_2017_isLoaded = false;
 fatJet2Pt_JESUp_HF_2018_isLoaded = false;
 fatJet2PtOverMHH_JESUp_HF_2018_isLoaded = false;
 fatJet2Pt_JESDown_HF_2018_isLoaded = false;
 fatJet2PtOverMHH_JESDown_HF_2018_isLoaded = false;
 fatJet2Pt_JESUp_RelBal_isLoaded = false;
 fatJet2PtOverMHH_JESUp_RelBal_isLoaded = false;
 fatJet2Pt_JESDown_RelBal_isLoaded = false;
 fatJet2PtOverMHH_JESDown_RelBal_isLoaded = false;
 fatJet2Pt_JESUp_RelSample_2016_isLoaded = false;
 fatJet2PtOverMHH_JESUp_RelSample_2016_isLoaded = false;
 fatJet2Pt_JESDown_RelSample_2016_isLoaded = false;
 fatJet2PtOverMHH_JESDown_RelSample_2016_isLoaded = false;
 fatJet2Pt_JESUp_RelSample_2017_isLoaded = false;
 fatJet2PtOverMHH_JESUp_RelSample_2017_isLoaded = false;
 fatJet2Pt_JESDown_RelSample_2017_isLoaded = false;
 fatJet2PtOverMHH_JESDown_RelSample_2017_isLoaded = false;
 fatJet2Pt_JESUp_RelSample_2018_isLoaded = false;
 fatJet2PtOverMHH_JESUp_RelSample_2018_isLoaded = false;
 fatJet2Pt_JESDown_RelSample_2018_isLoaded = false;
 fatJet2PtOverMHH_JESDown_RelSample_2018_isLoaded = false;
 hh_pt_isLoaded = false;
 hh_eta_isLoaded = false;
 hh_phi_isLoaded = false;
 hh_mass_isLoaded = false;
 hh_pt_JMS_Down_isLoaded = false;
 hh_pt_JMS_Up_isLoaded = false;
 hh_eta_JMS_Down_isLoaded = false;
 hh_eta_JMS_Up_isLoaded = false;
 hh_mass_JMS_Down_isLoaded = false;
 hh_mass_JMS_Up_isLoaded = false;
 hh_pt_JMR_Down_isLoaded = false;
 hh_pt_JMR_Up_isLoaded = false;
 hh_eta_JMR_Down_isLoaded = false;
 hh_eta_JMR_Up_isLoaded = false;
 hh_mass_JMR_Down_isLoaded = false;
 hh_mass_JMR_Up_isLoaded = false;
 hh_pt_JERUp_isLoaded = false;
 hh_eta_JERUp_isLoaded = false;
 hh_mass_JERUp_isLoaded = false;
 hh_pt_JERDown_isLoaded = false;
 hh_eta_JERDown_isLoaded = false;
 hh_mass_JERDown_isLoaded = false;
 hh_pt_JESUp_isLoaded = false;
 hh_eta_JESUp_isLoaded = false;
 hh_mass_JESUp_isLoaded = false;
 hh_pt_JESDown_isLoaded = false;
 hh_eta_JESDown_isLoaded = false;
 hh_mass_JESDown_isLoaded = false;

 hh_pt_JESUp_Abs_isLoaded = false;
 hh_eta_JESUp_Abs_isLoaded = false;
 hh_mass_JESUp_Abs_isLoaded = false;
 hh_pt_JESDown_Abs_isLoaded = false;
 hh_eta_JESDown_Abs_isLoaded = false;
 hh_mass_JESDown_Abs_isLoaded = false;
 hh_pt_JESUp_Abs_2016_isLoaded = false;
 hh_eta_JESUp_Abs_2016_isLoaded = false;
 hh_mass_JESUp_Abs_2016_isLoaded = false;
 hh_pt_JESDowb_Abs_2016_isLoaded = false;
 hh_eta_JESDowb_Abs_2016_isLoaded = false;
 hh_mass_JESDowb_Abs_2016_isLoaded = false;
 hh_pt_JESUp_Abs_2017_isLoaded = false;
 hh_eta_JESUp_Abs_2017_isLoaded = false;
 hh_mass_JESUp_Abs_2017_isLoaded = false;
 hh_pt_JESDowb_Abs_2017_isLoaded = false;
 hh_eta_JESDowb_Abs_2017_isLoaded = false;
 hh_mass_JESDowb_Abs_2017_isLoaded = false;
 hh_pt_JESUp_Abs_2018_isLoaded = false;
 hh_eta_JESUp_Abs_2018_isLoaded = false;
 hh_mass_JESUp_Abs_2018_isLoaded = false;
 hh_pt_JESDowb_Abs_2018_isLoaded = false;
 hh_eta_JESDowb_Abs_2018_isLoaded = false;
 hh_mass_JESDowb_Abs_2018_isLoaded = false;
 hh_pt_JESUp_BBEC1_isLoaded = false;
 hh_eta_JESUp_BBEC1_isLoaded = false;
 hh_mass_JESUp_BBEC1_isLoaded = false;
 hh_pt_JESDown_BBEC1_isLoaded = false;
 hh_eta_JESDown_BBEC1_isLoaded = false;
 hh_mass_JESDown_BBEC1_isLoaded = false;
 hh_pt_JESUp_BBEC1_2016_isLoaded = false;
 hh_eta_JESUp_BBEC1_2016_isLoaded = false;
 hh_mass_JESUp_BBEC1_2016_isLoaded = false;
 hh_pt_JESDown_BBEC1_2016_isLoaded = false;
 hh_eta_JESDown_BBEC1_2016_isLoaded = false;
 hh_mass_JESDown_BBEC1_2016_isLoaded = false;
 hh_pt_JESUp_BBEC1_2017_isLoaded = false;
 hh_eta_JESUp_BBEC1_2017_isLoaded = false;
 hh_mass_JESUp_BBEC1_2017_isLoaded = false;
 hh_pt_JESDown_BBEC1_2017_isLoaded = false;
 hh_eta_JESDown_BBEC1_2017_isLoaded = false;
 hh_mass_JESDown_BBEC1_2017_isLoaded = false;
 hh_pt_JESUp_BBEC1_2018_isLoaded = false;
 hh_eta_JESUp_BBEC1_2018_isLoaded = false;
 hh_mass_JESUp_BBEC1_2018_isLoaded = false;
 hh_pt_JESDown_BBEC1_2018_isLoaded = false;
 hh_eta_JESDown_BBEC1_2018_isLoaded = false;
 hh_mass_JESDown_BBEC1_2018_isLoaded = false;
 hh_pt_JESUp_EC2_isLoaded = false;
 hh_eta_JESUp_EC2_isLoaded = false;
 hh_mass_JESUp_EC2_isLoaded = false;
 hh_pt_JESDown_EC2_isLoaded = false;
 hh_eta_JESDown_EC2_isLoaded = false;
 hh_mass_JESDown_EC2_isLoaded = false;
 hh_pt_JESUp_EC2_2016_isLoaded = false;
 hh_eta_JESUp_EC2_2016_isLoaded = false;
 hh_mass_JESUp_EC2_2016_isLoaded = false;
 hh_pt_JESDown_EC2_2016_isLoaded = false;
 hh_eta_JESDown_EC2_2016_isLoaded = false;
 hh_mass_JESDown_EC2_2016_isLoaded = false;
 hh_pt_JESUp_EC2_2017_isLoaded = false;
 hh_eta_JESUp_EC2_2017_isLoaded = false;
 hh_mass_JESUp_EC2_2017_isLoaded = false;
 hh_pt_JESDown_EC2_2017_isLoaded = false;
 hh_eta_JESDown_EC2_2017_isLoaded = false;
 hh_mass_JESDown_EC2_2017_isLoaded = false;
 hh_pt_JESUp_EC2_2018_isLoaded = false;
 hh_eta_JESUp_EC2_2018_isLoaded = false;
 hh_mass_JESUp_EC2_2018_isLoaded = false;
 hh_pt_JESDown_EC2_2018_isLoaded = false;
 hh_eta_JESDown_EC2_2018_isLoaded = false;
 hh_mass_JESDown_EC2_2018_isLoaded = false;
 hh_pt_JESUp_FlavQCD_isLoaded = false;
 hh_eta_JESUp_FlavQCD_isLoaded = false;
 hh_mass_JESUp_FlavQCD_isLoaded = false;
 hh_pt_JESDown_FlavQCD_isLoaded = false;
 hh_eta_JESDown_FlavQCD_isLoaded = false;
 hh_mass_JESDown_FlavQCD_isLoaded = false;
 hh_pt_JESUp_HF_isLoaded = false;
 hh_eta_JESUp_HF_isLoaded = false;
 hh_mass_JESUp_HF_isLoaded = false;
 hh_pt_JESDown_HF_isLoaded = false;
 hh_eta_JESDown_HF_isLoaded = false;
 hh_mass_JESDown_HF_isLoaded = false;
 hh_pt_JESUp_HF_2016_isLoaded = false;
 hh_eta_JESUp_HF_2016_isLoaded = false;
 hh_mass_JESUp_HF_2016_isLoaded = false;
 hh_pt_JESDown_HF_2016_isLoaded = false;
 hh_eta_JESDown_HF_2016_isLoaded = false;
 hh_mass_JESDown_HF_2016_isLoaded = false;
 hh_pt_JESUp_HF_2017_isLoaded = false;
 hh_eta_JESUp_HF_2017_isLoaded = false;
 hh_mass_JESUp_HF_2017_isLoaded = false;
 hh_pt_JESDown_HF_2017_isLoaded = false;
 hh_eta_JESDown_HF_2017_isLoaded = false;
 hh_mass_JESDown_HF_2017_isLoaded = false;
 hh_pt_JESUp_HF_2018_isLoaded = false;
 hh_eta_JESUp_HF_2018_isLoaded = false;
 hh_mass_JESUp_HF_2018_isLoaded = false;
 hh_pt_JESDown_HF_2018_isLoaded = false;
 hh_eta_JESDown_HF_2018_isLoaded = false;
 hh_mass_JESDown_HF_2018_isLoaded = false;
 hh_pt_JESUp_RelBal_isLoaded = false;
 hh_eta_JESUp_RelBal_isLoaded = false;
 hh_mass_JESUp_RelBal_isLoaded = false;
 hh_pt_JESDown_RelBal_isLoaded = false;
 hh_eta_JESDown_RelBal_isLoaded = false;
 hh_mass_JESDown_RelBal_isLoaded = false;
 hh_pt_JESUp_RelSample_2016_isLoaded = false;
 hh_eta_JESUp_RelSample_2016_isLoaded = false;
 hh_mass_JESUp_RelSample_2016_isLoaded = false;
 hh_pt_JESDown_RelSample_2016_isLoaded = false;
 hh_eta_JESDown_RelSample_2016_isLoaded = false;
 hh_mass_JESDown_RelSample_2016_isLoaded = false;
 hh_pt_JESUp_RelSample_2017_isLoaded = false;
 hh_eta_JESUp_RelSample_2017_isLoaded = false;
 hh_mass_JESUp_RelSample_2017_isLoaded = false;
 hh_pt_JESDown_RelSample_2017_isLoaded = false;
 hh_eta_JESDown_RelSample_2017_isLoaded = false;
 hh_mass_JESDown_RelSample_2017_isLoaded = false;
 hh_pt_JESUp_RelSample_2018_isLoaded = false;
 hh_eta_JESUp_RelSample_2018_isLoaded = false;
 hh_mass_JESUp_RelSample_2018_isLoaded = false;
 hh_pt_JESDown_RelSample_2018_isLoaded = false;
 hh_eta_JESDown_RelSample_2018_isLoaded = false;
 hh_mass_JESDown_RelSample_2018_isLoaded = false;
 deltaEta_j1j2_isLoaded = false;
 deltaPhi_j1j2_isLoaded = false;
 deltaR_j1j2_isLoaded = false;
 ptj2_over_ptj1_isLoaded = false;
 mj2_over_mj1_isLoaded = false;
 isVBFtag_isLoaded = false;
 isVBFtag_JERUp_isLoaded = false;
 isVBFtag_JERDown_isLoaded = false;
 isVBFtag_JESUp_isLoaded = false;
 isVBFtag_JESDown_isLoaded = false;
 isVBFtag_JESUp_Abs_isLoaded = false;
 isVBFtag_JESDown_Abs_isLoaded = false;
 isVBFtag_JESUp_Abs_2016_isLoaded = false;
 isVBFtag_JESDowb_Abs_2016_isLoaded = false;
 isVBFtag_JESUp_Abs_2017_isLoaded = false;
 isVBFtag_JESDowb_Abs_2017_isLoaded = false;
isVBFtag_JESUp_Abs_2018_isLoaded = false;
 isVBFtag_JESDowb_Abs_2018_isLoaded = false;
 isVBFtag_JESUp_BBEC1_isLoaded = false;
 isVBFtag_JESDown_BBEC1_isLoaded = false;
 isVBFtag_JESUp_BBEC1_2016_isLoaded = false;
 isVBFtag_JESDown_BBEC1_2016_isLoaded = false;
 isVBFtag_JESUp_BBEC1_2017_isLoaded = false;
 isVBFtag_JESDown_BBEC1_2017_isLoaded = false;
isVBFtag_JESUp_BBEC1_2018_isLoaded = false;
 isVBFtag_JESDown_BBEC1_2018_isLoaded = false;
 isVBFtag_JESUp_EC2_isLoaded = false;
 isVBFtag_JESDown_EC2_isLoaded = false;
 isVBFtag_JESUp_EC2_2016_isLoaded = false;
 isVBFtag_JESDown_EC2_2016_isLoaded = false;
 isVBFtag_JESUp_EC2_2017_isLoaded = false;
 isVBFtag_JESDown_EC2_2017_isLoaded = false;
 isVBFtag_JESUp_EC2_2018_isLoaded = false;
 isVBFtag_JESDown_EC2_2018_isLoaded = false;
 isVBFtag_JESUp_FlavQCD_isLoaded = false;
 isVBFtag_JESDown_FlavQCD_isLoaded = false;
 isVBFtag_JESUp_HF_isLoaded = false;
 isVBFtag_JESDown_HF_isLoaded = false;
 isVBFtag_JESUp_HF_2016_isLoaded = false;
 isVBFtag_JESDown_HF_2016_isLoaded = false;
 isVBFtag_JESUp_HF_2017_isLoaded = false;
 isVBFtag_JESDown_HF_2017_isLoaded = false;
 isVBFtag_JESUp_HF_2018_isLoaded = false;
 isVBFtag_JESDown_HF_2018_isLoaded = false;
 isVBFtag_JESUp_RelBal_isLoaded = false;
 isVBFtag_JESDown_RelBal_isLoaded = false;
 isVBFtag_JESUp_RelSample_2016_isLoaded = false;
 isVBFtag_JESDown_RelSample_2016_isLoaded = false;
 isVBFtag_JESUp_RelSample_2017_isLoaded = false;
 isVBFtag_JESDown_RelSample_2017_isLoaded = false;
 isVBFtag_JESUp_RelSample_2018_isLoaded = false;
 isVBFtag_JESDown_RelSample_2018_isLoaded = false;

 dijetmass_isLoaded = false;
 vbfjet1Pt_isLoaded = false;
 vbfjet1Eta_isLoaded = false;
 vbfjet1Phi_isLoaded = false;
 vbfjet1Mass_isLoaded = false;
 vbffatJet1Pt_isLoaded = false;
 vbffatJet1Eta_isLoaded = false;
 vbffatJet1Phi_isLoaded = false;
 vbffatJet1PNetXbb_isLoaded = false;
 vbfjet2Pt_isLoaded = false;
 vbfjet2Eta_isLoaded = false;
 vbfjet2Phi_isLoaded = false;
 vbfjet2Mass_isLoaded = false;
 vbffatJet2Pt_isLoaded = false;
 vbffatJet2Eta_isLoaded = false;
 vbffatJet2Phi_isLoaded = false;
 vbffatJet2PNetXbb_isLoaded = false;
 jet1Pt_isLoaded = false;
 jet1Eta_isLoaded = false;
 jet1Phi_isLoaded = false;
 jet2Pt_isLoaded = false;
 jet2Eta_isLoaded = false;
 jet2Phi_isLoaded = false;
 jet3Pt_isLoaded = false;
 jet3Eta_isLoaded = false;
 jet3Phi_isLoaded = false;
 jet4Pt_isLoaded = false;
 jet4Eta_isLoaded = false;
 jet4Phi_isLoaded = false;
 lep1Pt_isLoaded = false;
 lep1Eta_isLoaded = false;
 lep1Phi_isLoaded = false;
 lep1Id_isLoaded = false;
 lep2Pt_isLoaded = false;
 lep2Eta_isLoaded = false;
 lep2Phi_isLoaded = false;
 lep2Id_isLoaded = false;
 genHiggs1Pt_isLoaded = false;
 genHiggs1Eta_isLoaded = false;
 genHiggs1Phi_isLoaded = false;
 genHiggs2Pt_isLoaded = false;
 genHiggs2Eta_isLoaded = false;
 genHiggs2Phi_isLoaded = false;
 puWeight_isLoaded = false;
 puWeightUp_isLoaded = false;
 puWeightDown_isLoaded = false;
 xsecWeight_isLoaded = false;
 LHEScaleWeightNormNew_isLoaded = false;
 LHEPdfWeightNorm_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_isLoaded = false;

 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_isLoaded = false;
disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_isLoaded = false;
disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_isLoaded = false;
 disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_isLoaded = false;
}

const unsigned int &hhtree::run() 
{
  if(not run_isLoaded)
    {
      if(run_branch != 0) run_branch->GetEntry(index);
      else
	{
	  printf("branch run_branch does not exist!\n");
	  //exit(1);
	}
      run_isLoaded = true;
    }
  return run_;
}

const unsigned int &hhtree::luminosityBlock() 
{
  if(not luminosityBlock_isLoaded)
    {
      if(luminosityBlock_branch != 0) luminosityBlock_branch->GetEntry(index);
      else
	{
	  printf("branch luminosityBlock_branch does not exist!\n");
	  //exit(1);
	}
      luminosityBlock_isLoaded = true;
    }
  return luminosityBlock_;
}

const unsigned long int &hhtree::event() 
{
  if(not event_isLoaded)
    {
      if(event_branch != 0) event_branch->GetEntry(index);
      else
	{
	  printf("branch event_branch does not exist!\n");
	  //exit(1);
	}
      event_isLoaded = true;
    }
  return event_;
}

const float &hhtree::genWeight() 
{
  if(not genWeight_isLoaded)
    {
      if(genWeight_branch != 0) genWeight_branch->GetEntry(index);
      else
	{
	  printf("branch genWeight_branch does not exist!\n");
	  //exit(1);
	}
      genWeight_isLoaded = true;
    }
  return genWeight_;
}
const int &hhtree::nLHEPdfWeight()
{
  if(not nLHEPdfWeight_isLoaded)
    {
      if(nLHEPdfWeight_branch != 0){
	nLHEPdfWeight_branch->GetEntry(index);
	nLHEPdfWeight_isLoaded = true;
      } 
      //else
      //{
      //printf("branch nLHEPdfWeight_branch does not exist!\n");
      //exit(1);
      //}  
    }
  return nLHEPdfWeight_;
}
const float * hhtree::LHEPdfWeight()
{
  if(not LHEPdfWeight_isLoaded)
    {
      if(LHEPdfWeight_branch != 0){
	LHEPdfWeight_branch->GetEntry(index);
	LHEPdfWeight_isLoaded = true;
      } 
      //else
      //{
      //printf("branch LHEPdfWeight_branch does not exist!\n");
      //exit(1);
      //}
    }
  return LHEPdfWeight_;
}
const int &hhtree::nLHEScaleWeight()
{
  if(not nLHEScaleWeight_isLoaded)
    {
      if(nLHEScaleWeight_branch != 0){
	nLHEScaleWeight_branch->GetEntry(index);
	nLHEScaleWeight_isLoaded = true;
      } 
      //else
      //{
      //printf("branch nLHEScaleWeight_branch does not exist!\n");
      //exit(1);
      //}
    }
  return nLHEScaleWeight_;
}

const float * hhtree::LHEScaleWeight()
{
  if(not LHEScaleWeight_isLoaded)
    {
      if(LHEScaleWeight_branch != 0){
	LHEScaleWeight_branch->GetEntry(index);
	LHEScaleWeight_isLoaded = true;
      } 
      //else
      //{
      //printf("branch LHEScaleWeight_branch does not exist!\n");
      //exit(1);
      //}  
    }
  return LHEScaleWeight_;
}
const int &hhtree::nPSWeight()
{
  if(not nPSWeight_isLoaded)
    {
      if(nPSWeight_branch != 0){
	nPSWeight_branch->GetEntry(index);
	nPSWeight_isLoaded = true;
      }
      //else
      //{
      //printf("branch nPSWeight_branch does not exist!\n");
      //exit(1);
      //}  
    }
  return nPSWeight_;
}
const float * hhtree::PSWeight()
{
  if(not PSWeight_isLoaded)
    {
      if(PSWeight_branch != 0){
	PSWeight_branch->GetEntry(index);
	PSWeight_isLoaded = true;
      } 
      //else
      //{
      //printf("branch PSWeight_branch does not exist!\n");
      //exit(1);
      //}   
    }
  return PSWeight_;
}
const bool &hhtree::HLT_Ele27_WPTight_Gsf() 
{
  if(not HLT_Ele27_WPTight_Gsf_isLoaded)
    {
      if(HLT_Ele27_WPTight_Gsf_branch != 0) HLT_Ele27_WPTight_Gsf_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele27_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele27_WPTight_Gsf_isLoaded = true;
    }
  return HLT_Ele27_WPTight_Gsf_;
}

const bool &hhtree::HLT_Ele28_WPTight_Gsf() 
{
  if(not HLT_Ele28_WPTight_Gsf_isLoaded)
    {
      if(HLT_Ele28_WPTight_Gsf_branch != 0) HLT_Ele28_WPTight_Gsf_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele28_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele28_WPTight_Gsf_isLoaded = true;
    }
  return HLT_Ele28_WPTight_Gsf_;
}

const bool &hhtree::HLT_Ele30_WPTight_Gsf() 
{
  if(not HLT_Ele30_WPTight_Gsf_isLoaded)
    {
      if(HLT_Ele30_WPTight_Gsf_branch != 0) HLT_Ele30_WPTight_Gsf_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele30_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele30_WPTight_Gsf_isLoaded = true;
    }
  return HLT_Ele30_WPTight_Gsf_;
}

const bool &hhtree::HLT_Ele32_WPTight_Gsf() 
{
  if(not HLT_Ele32_WPTight_Gsf_isLoaded)
    {
      if(HLT_Ele32_WPTight_Gsf_branch != 0) HLT_Ele32_WPTight_Gsf_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele32_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele32_WPTight_Gsf_isLoaded = true;
    }
  return HLT_Ele32_WPTight_Gsf_;
}
const bool &hhtree::HLT_Ele32_WPTight_Gsf_L1DoubleEG() 
{
  if(not HLT_Ele32_WPTight_Gsf_L1DoubleEG_isLoaded)
    {
      if(HLT_Ele32_WPTight_Gsf_L1DoubleEG_branch != 0) HLT_Ele32_WPTight_Gsf_L1DoubleEG_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele32_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele32_WPTight_Gsf_L1DoubleEG_isLoaded = true;
    }
  return HLT_Ele32_WPTight_Gsf_L1DoubleEG_;
}
const bool &hhtree::HLT_Ele35_WPTight_Gsf() 
{
  if(not HLT_Ele35_WPTight_Gsf_isLoaded)
    {
      if(HLT_Ele35_WPTight_Gsf_branch != 0) HLT_Ele35_WPTight_Gsf_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele35_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele35_WPTight_Gsf_isLoaded = true;
    }
  return HLT_Ele35_WPTight_Gsf_;
}

const bool &hhtree::HLT_Ele38_WPTight_Gsf() 
{
  if(not HLT_Ele38_WPTight_Gsf_isLoaded)
    {
      if(HLT_Ele38_WPTight_Gsf_branch != 0) HLT_Ele38_WPTight_Gsf_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele38_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele38_WPTight_Gsf_isLoaded = true;
    }
  return HLT_Ele38_WPTight_Gsf_;
}

const bool &hhtree::HLT_Ele40_WPTight_Gsf() 
{
  if(not HLT_Ele40_WPTight_Gsf_isLoaded)
    {
      if(HLT_Ele40_WPTight_Gsf_branch != 0) HLT_Ele40_WPTight_Gsf_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Ele40_WPTight_Gsf_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Ele40_WPTight_Gsf_isLoaded = true;
    }
  return HLT_Ele40_WPTight_Gsf_;
}

const bool &hhtree::HLT_IsoMu20() 
{
  if(not HLT_IsoMu20_isLoaded)
    {
      if(HLT_IsoMu20_branch != 0) HLT_IsoMu20_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_IsoMu20_branch does not exist!\n");
	  //exit(1);
	}
      HLT_IsoMu20_isLoaded = true;
    }
  return HLT_IsoMu20_;
}

const bool &hhtree::HLT_IsoMu24() 
{
  if(not HLT_IsoMu24_isLoaded)
    {
      if(HLT_IsoMu24_branch != 0) HLT_IsoMu24_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_IsoMu24_branch does not exist!\n");
	  //exit(1);
	}
      HLT_IsoMu24_isLoaded = true;
    }
  return HLT_IsoMu24_;
}

const bool &hhtree::HLT_IsoMu24_eta2p1() 
{
  if(not HLT_IsoMu24_eta2p1_isLoaded)
    {
      if(HLT_IsoMu24_eta2p1_branch != 0) HLT_IsoMu24_eta2p1_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_IsoMu24_eta2p1_branch does not exist!\n");
	  //exit(1);
	}
      HLT_IsoMu24_eta2p1_isLoaded = true;
    }
  return HLT_IsoMu24_eta2p1_;
}

const bool &hhtree::HLT_IsoMu27() 
{
  if(not HLT_IsoMu27_isLoaded)
    {
      if(HLT_IsoMu27_branch != 0) HLT_IsoMu27_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_IsoMu27_branch does not exist!\n");
	  //exit(1);
	}
      HLT_IsoMu27_isLoaded = true;
    }
  return HLT_IsoMu27_;
}

const bool &hhtree::HLT_IsoMu30() 
{
  if(not HLT_IsoMu30_isLoaded)
    {
      if(HLT_IsoMu30_branch != 0) HLT_IsoMu30_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_IsoMu30_branch does not exist!\n");
	  //exit(1);
	}
      HLT_IsoMu30_isLoaded = true;
    }
  return HLT_IsoMu30_;
}

const bool &hhtree::HLT_Mu50() 
{
  if(not HLT_Mu50_isLoaded)
    {
      if(HLT_Mu50_branch != 0) HLT_Mu50_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Mu50_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Mu50_isLoaded = true;
    }
  return HLT_Mu50_;
}

const bool &hhtree::HLT_Mu55() 
{
  if(not HLT_Mu55_isLoaded)
    {
      if(HLT_Mu55_branch != 0) HLT_Mu55_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Mu55_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Mu55_isLoaded = true;
    }
  return HLT_Mu55_;
}

const bool &hhtree::HLT_Photon175() 
{
  if(not HLT_Photon175_isLoaded)
    {
      if(HLT_Photon175_branch != 0) HLT_Photon175_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_Photon175_branch does not exist!\n");
	  //exit(1);
	}
      HLT_Photon175_isLoaded = true;
    }
  return HLT_Photon175_;
}

const bool &hhtree::HLT_PFHT780() 
{
  if(not HLT_PFHT780_isLoaded)
    {
      if(HLT_PFHT780_branch != 0) HLT_PFHT780_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_PFHT780_branch does not exist!\n");
	  //exit(1);
	}
      HLT_PFHT780_isLoaded = true;
    }
  return HLT_PFHT780_;
}

const bool &hhtree::HLT_PFHT890() 
{
  if(not HLT_PFHT890_isLoaded)
    {
      if(HLT_PFHT890_branch != 0) HLT_PFHT890_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_PFHT890_branch does not exist!\n");
	  //exit(1);
	}
      HLT_PFHT890_isLoaded = true;
    }
  return HLT_PFHT890_;
}

const bool &hhtree::HLT_PFHT1050() 
{
  if(not HLT_PFHT1050_isLoaded)
    {
      if(HLT_PFHT1050_branch != 0) HLT_PFHT1050_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_PFHT1050_branch does not exist!\n");
	  //exit(1);
	}
      HLT_PFHT1050_isLoaded = true;
    }
  return HLT_PFHT1050_;
}

const bool &hhtree::HLT_AK8PFJet360_TrimMass30() 
{
  if(not HLT_AK8PFJet360_TrimMass30_isLoaded)
    {
      if(HLT_AK8PFJet360_TrimMass30_branch != 0){
	HLT_AK8PFJet360_TrimMass30_branch->GetEntry(index);
	//else
	//{
	//  printf("branch HLT_AK8PFJet360_TrimMass30_branch does not exist!\n");
	//exit(1);
	//}
	HLT_AK8PFJet360_TrimMass30_isLoaded = true;
      }
    }
  return HLT_AK8PFJet360_TrimMass30_;
}

const bool &hhtree::HLT_AK8PFJet380_TrimMass30() 
{
  if(not HLT_AK8PFJet380_TrimMass30_isLoaded)
    {
      if(HLT_AK8PFJet380_TrimMass30_branch != 0){
	HLT_AK8PFJet380_TrimMass30_branch->GetEntry(index);
	//else
	//{
	//  printf("branch HLT_AK8PFJet380_TrimMass30_branch does not exist!\n");
	//exit(1);
	//}
	HLT_AK8PFJet380_TrimMass30_isLoaded = true;
      }
    }
  return HLT_AK8PFJet380_TrimMass30_;
}

const bool &hhtree::HLT_AK8PFJet400_TrimMass30() 
{
  if(not HLT_AK8PFJet400_TrimMass30_isLoaded)
    {
      if(HLT_AK8PFJet400_TrimMass30_branch != 0){
	HLT_AK8PFJet400_TrimMass30_branch->GetEntry(index);
	//else
	//{
	//  printf("branch HLT_AK8PFJet400_TrimMass30_branch does not exist!\n");
	//exit(1);
	//}
	HLT_AK8PFJet400_TrimMass30_isLoaded = true;
      }
    }
  return HLT_AK8PFJet400_TrimMass30_;
}

const bool &hhtree::HLT_AK8PFJet420_TrimMass30() 
{
  if(not HLT_AK8PFJet420_TrimMass30_isLoaded)
    {
      if(HLT_AK8PFJet420_TrimMass30_branch != 0) HLT_AK8PFJet420_TrimMass30_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet420_TrimMass30_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet420_TrimMass30_isLoaded = true;
    }
  return HLT_AK8PFJet420_TrimMass30_;
}

const bool &hhtree::HLT_AK8PFHT750_TrimMass50() 
{
  if(not HLT_AK8PFHT750_TrimMass50_isLoaded)
    {
      if(HLT_AK8PFHT750_TrimMass50_branch != 0) {
	HLT_AK8PFHT750_TrimMass50_branch->GetEntry(index);
	//else
	//{
	//  printf("branch HLT_AK8PFHT750_TrimMass50_branch does not exist!\n");
	//exit(1);
	//}
	HLT_AK8PFHT750_TrimMass50_isLoaded = true;
      }
    }
  return HLT_AK8PFHT750_TrimMass50_;
}

const bool &hhtree::HLT_AK8PFHT800_TrimMass50() 
{
  if(not HLT_AK8PFHT800_TrimMass50_isLoaded)
    {
      if(HLT_AK8PFHT800_TrimMass50_branch != 0){
	HLT_AK8PFHT800_TrimMass50_branch->GetEntry(index);
	//else
	//{
	//  printf("branch HLT_AK8PFHT800_TrimMass50_branch does not exist!\n");
	//exit(1);
	//}
	HLT_AK8PFHT800_TrimMass50_isLoaded = true;
      }
    }
  return HLT_AK8PFHT800_TrimMass50_;
}

const bool &hhtree::HLT_AK8PFHT850_TrimMass50() 
{
  if(not HLT_AK8PFHT850_TrimMass50_isLoaded)
    {
      if(HLT_AK8PFHT850_TrimMass50_branch != 0) HLT_AK8PFHT850_TrimMass50_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFHT850_TrimMass50_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFHT850_TrimMass50_isLoaded = true;
    }
  return HLT_AK8PFHT850_TrimMass50_;
}

const bool &hhtree::HLT_AK8PFHT900_TrimMass50() 
{
  if(not HLT_AK8PFHT900_TrimMass50_isLoaded)
    {
      if(HLT_AK8PFHT900_TrimMass50_branch != 0) HLT_AK8PFHT900_TrimMass50_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFHT900_TrimMass50_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFHT900_TrimMass50_isLoaded = true;
    }
  return HLT_AK8PFHT900_TrimMass50_;
}

const bool &hhtree::HLT_PFJet450() 
{
  if(not HLT_PFJet450_isLoaded)
    {
      if(HLT_PFJet450_branch != 0) HLT_PFJet450_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_PFJet450_branch does not exist!\n");
	  //exit(1);
	}
      HLT_PFJet450_isLoaded = true;
    }
  return HLT_PFJet450_;
}

const bool &hhtree::HLT_PFJet500() 
{
  if(not HLT_PFJet500_isLoaded)
    {
      if(HLT_PFJet500_branch != 0) HLT_PFJet500_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_PFJet500_branch does not exist!\n");
	  //exit(1);
	}
      HLT_PFJet500_isLoaded = true;
    }
  return HLT_PFJet500_;
}

const bool &hhtree::HLT_PFJet550() 
{
  if(not HLT_PFJet550_isLoaded)
    {
      if(HLT_PFJet550_branch != 0) HLT_PFJet550_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_PFJet550_branch does not exist!\n");
	  //exit(1);
	}
      HLT_PFJet550_isLoaded = true;
    }
  return HLT_PFJet550_;
}

const bool &hhtree::HLT_AK8PFJet400() 
{
  if(not HLT_AK8PFJet400_isLoaded)
    {
      if(HLT_AK8PFJet400_branch != 0) HLT_AK8PFJet400_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet400_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet400_isLoaded = true;
    }
  return HLT_AK8PFJet400_;
}

const bool &hhtree::HLT_AK8PFJet450() 
{
  if(not HLT_AK8PFJet450_isLoaded)
    {
      if(HLT_AK8PFJet450_branch != 0) HLT_AK8PFJet450_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet450_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet450_isLoaded = true;
    }
  return HLT_AK8PFJet450_;
}

const bool &hhtree::HLT_AK8PFJet500() 
{
  if(not HLT_AK8PFJet500_isLoaded)
    {
      if(HLT_AK8PFJet500_branch != 0) HLT_AK8PFJet500_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet500_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet500_isLoaded = true;
    }
  return HLT_AK8PFJet500_;
}

const bool &hhtree::HLT_AK8PFJet550() 
{
  if(not HLT_AK8PFJet550_isLoaded)
    {
      if(HLT_AK8PFJet550_branch != 0) HLT_AK8PFJet550_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet550_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet550_isLoaded = true;
    }
  return HLT_AK8PFJet550_;
}

const bool &hhtree::HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17() 
{
  if(not HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_isLoaded)
    {
      if(HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_branch != 0) HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_isLoaded = true;
    }
  return HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_;
}

const bool &hhtree::HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1() 
{
  if(not HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_isLoaded)
    {
      if(HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_branch != 0) HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_isLoaded = true;
    }
  return HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_;
}

const bool &hhtree::HLT_AK8PFJet330_PFAK8BTagCSV_p17() 
{
  if(not HLT_AK8PFJet330_PFAK8BTagCSV_p17_isLoaded)
    {
      if(HLT_AK8PFJet330_PFAK8BTagCSV_p17_branch != 0){
	HLT_AK8PFJet330_PFAK8BTagCSV_p17_branch->GetEntry(index);
	HLT_AK8PFJet330_PFAK8BTagCSV_p17_isLoaded = true;
      }
      //else
      //{
      //printf("branch HLT_AK8PFJet330_PFAK8BTagCSV_p17_branch does not exist!\n");
      //exit(1);
      //}
    }
  return HLT_AK8PFJet330_PFAK8BTagCSV_p17_;
}

const bool &hhtree::HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02() 
{
  if(not HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_isLoaded)
    {
      if(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_branch != 0) HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_isLoaded = true;
    }
  return HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_;
}

const bool &hhtree::HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2() 
{
  if(not HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_isLoaded)
    {
      if(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_branch != 0) HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_isLoaded = true;
    }
  return HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_;
}

const bool &hhtree::HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4() 
{
  if(not HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_isLoaded)
    {
      if(HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_branch != 0){
	HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_branch->GetEntry(index);
	//else
	//{
	//printf("branch HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_branch does not exist!\n");
	//exit(1);
	//}
	HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_isLoaded = true;
      }
    }
  return HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_;
}

const bool &hhtree::HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20() 
{
  if(not HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_isLoaded)
    {
      if(HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_branch != 0) HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_isLoaded = true;
    }
  return HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p20_;
}

const bool &hhtree::HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087() 
{
  if(not HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_isLoaded)
    {
      if(HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_branch != 0) HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_isLoaded = true;
    }
  return HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p087_;
}

const bool &hhtree::HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087() 
{
  if(not HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_isLoaded)
    {
      if(HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_branch != 0) HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_isLoaded = true;
    }
  return HLT_AK8DiPFJet300_200_TrimMass30_BTagCSV_p087_;
}

const bool &hhtree::HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20() 
{
  if(not HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_isLoaded)
    {
      if(HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_branch != 0) HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_isLoaded = true;
    }
  return HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_;
}

const bool &hhtree::HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20() 
{
  if(not HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_isLoaded)
    {
      if(HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_branch != 0) HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_isLoaded = true;
    }
  return HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_;
}

const bool &hhtree::HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20() 
{
  if(not HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_isLoaded)
    {
      if(HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_branch != 0) HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_branch->GetEntry(index);
      else
	{
	  //printf("branch HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_branch does not exist!\n");
	  //exit(1);
	}
      HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_isLoaded = true;
    }
  return HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_;
}

const float &hhtree::weight() 
{
  if(not weight_isLoaded)
    {
      if(weight_branch != 0) weight_branch->GetEntry(index);
      else
	{
	  printf("branch weight_branch does not exist!\n");
	  //exit(1);
	}
      weight_isLoaded = true;
    }
  return weight_;
}

const float &hhtree::met() 
{
  if(not met_isLoaded)
    {
      if(met_branch != 0) met_branch->GetEntry(index);
      else
	{
	  printf("branch met_branch does not exist!\n");
	  //exit(1);
	}
      met_isLoaded = true;
    }
  return met_;
}

const float &hhtree::metphi() 
{
  if(not metphi_isLoaded)
    {
      if(metphi_branch != 0) metphi_branch->GetEntry(index);
      else
	{
	  printf("branch metphi_branch does not exist!\n");
	  //exit(1);
	}
      metphi_isLoaded = true;
    }
  return metphi_;
}

const float &hhtree::ht() 
{
  if(not ht_isLoaded)
    {
      if(ht_branch != 0) ht_branch->GetEntry(index);
      else
	{
	  printf("branch ht_branch does not exist!\n");
	  //exit(1);
	}
      ht_isLoaded = true;
    }
  return ht_;
}

const bool &hhtree::passmetfilters() 
{
  if(not passmetfilters_isLoaded)
    {
      if(passmetfilters_branch != 0) passmetfilters_branch->GetEntry(index);
      else
	{
	  printf("branch passmetfilters_branch does not exist!\n");
	  //exit(1);
	}
      passmetfilters_isLoaded = true;
    }
  return passmetfilters_;
}

const float &hhtree::l1PreFiringWeight() 
{
  if(not l1PreFiringWeight_isLoaded)
    {
      if(l1PreFiringWeight_branch != 0) l1PreFiringWeight_branch->GetEntry(index);
      else
	{
	  printf("branch l1PreFiringWeight_branch does not exist!\n");
	  //exit(1);
	}
      l1PreFiringWeight_isLoaded = true;
    }
  return l1PreFiringWeight_;
}

const float &hhtree::l1PreFiringWeightUp() 
{
  if(not l1PreFiringWeightUp_isLoaded)
    {
      if(l1PreFiringWeightUp_branch != 0) l1PreFiringWeightUp_branch->GetEntry(index);
      else
	{
	  printf("branch l1PreFiringWeightUp_branch does not exist!\n");
	  //exit(1);
	}
      l1PreFiringWeightUp_isLoaded = true;
    }
  return l1PreFiringWeightUp_;
}

const float &hhtree::l1PreFiringWeightDown() 
{
  if(not l1PreFiringWeightDown_isLoaded)
    {
      if(l1PreFiringWeightDown_branch != 0) l1PreFiringWeightDown_branch->GetEntry(index);
      else
	{
	  printf("branch l1PreFiringWeightDown_branch does not exist!\n");
	  //exit(1);
	}
      l1PreFiringWeightDown_isLoaded = true;
    }
  return l1PreFiringWeightDown_;
}

const float &hhtree::triggerEffWeight() 
{
  if(not triggerEffWeight_isLoaded)
    {
      if(triggerEffWeight_branch != 0) triggerEffWeight_branch->GetEntry(index);
      else
	{
	  printf("branch triggerEffWeight_branch does not exist!\n");
	  //exit(1);
	}
      triggerEffWeight_isLoaded = true;
    }
  return triggerEffWeight_;
}

const float &hhtree::triggerEff3DWeight() 
{
  if(not triggerEff3DWeight_isLoaded)
    {
      if(triggerEff3DWeight_branch != 0) triggerEff3DWeight_branch->GetEntry(index);
      else
	{
	  printf("branch triggerEff3DWeight_branch does not exist!\n");
	  //exit(1);
	}
      triggerEff3DWeight_isLoaded = true;
    }
  return triggerEff3DWeight_;
}

const float &hhtree::triggerEffMCWeight() 
{
  if(not triggerEffMCWeight_isLoaded)
    {
      if(triggerEffMCWeight_branch != 0) triggerEffMCWeight_branch->GetEntry(index);
      else
	{
	  printf("branch triggerEffMCWeight_branch does not exist!\n");
	  //exit(1);
	}
      triggerEffMCWeight_isLoaded = true;
    }
  return triggerEffMCWeight_;
}

const float &hhtree::triggerEffMC3DWeight() 
{
  if(not triggerEffMC3DWeight_isLoaded)
    {
      if(triggerEffMC3DWeight_branch != 0) triggerEffMC3DWeight_branch->GetEntry(index);
      else
	{
	  printf("branch triggerEffMC3DWeight_branch does not exist!\n");
	  //exit(1);
	}
      triggerEffMC3DWeight_isLoaded = true;
    }
  return triggerEffMC3DWeight_;
}

const float &hhtree::fatJet1Pt() 
{
  if(not fatJet1Pt_isLoaded)
    {
      if(fatJet1Pt_branch != 0) fatJet1Pt_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1Pt_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1Pt_isLoaded = true;
    }
  return fatJet1Pt_;
}

const float &hhtree::fatJet1Eta() 
{
  if(not fatJet1Eta_isLoaded)
    {
      if(fatJet1Eta_branch != 0) fatJet1Eta_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1Eta_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1Eta_isLoaded = true;
    }
  return fatJet1Eta_;
}

const float &hhtree::fatJet1Phi() 
{
  if(not fatJet1Phi_isLoaded)
    {
      if(fatJet1Phi_branch != 0) fatJet1Phi_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1Phi_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1Phi_isLoaded = true;
    }
  return fatJet1Phi_;
}

const float &hhtree::fatJet1Mass() 
{
  if(not fatJet1Mass_isLoaded)
    {
      if(fatJet1Mass_branch != 0) fatJet1Mass_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1Mass_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1Mass_isLoaded = true;
    }
  return fatJet1Mass_;
}

const float &hhtree::fatJet1MassSD() 
{
  if(not fatJet1MassSD_isLoaded)
    {
      if(fatJet1MassSD_branch != 0) fatJet1MassSD_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassSD_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassSD_isLoaded = true;
    }
  return fatJet1MassSD_;
}

const float &hhtree::fatJet1MassSD_UnCorrected() 
{
  if(not fatJet1MassSD_UnCorrected_isLoaded)
    {
      if(fatJet1MassSD_UnCorrected_branch != 0) fatJet1MassSD_UnCorrected_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassSD_UnCorrected_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassSD_UnCorrected_isLoaded = true;
    }
  return fatJet1MassSD_UnCorrected_;
}

const float &hhtree::fatJet1MassRegressed() 
{
  if(not fatJet1MassRegressed_isLoaded)
    {
      if(fatJet1MassRegressed_branch != 0) fatJet1MassRegressed_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassRegressed_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassRegressed_isLoaded = true;
    }
  return fatJet1MassRegressed_;
}

const float &hhtree::fatJet1MassRegressed_UnCorrected() 
{
  if(not fatJet1MassRegressed_UnCorrected_isLoaded)
    {
      if(fatJet1MassRegressed_UnCorrected_branch != 0) fatJet1MassRegressed_UnCorrected_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassRegressed_UnCorrected_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassRegressed_UnCorrected_isLoaded = true;
    }
  return fatJet1MassRegressed_UnCorrected_;
}

const float &hhtree::fatJet1PNetXbb() 
{
  if(not fatJet1PNetXbb_isLoaded)
    {
      if(fatJet1PNetXbb_branch != 0) fatJet1PNetXbb_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PNetXbb_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PNetXbb_isLoaded = true;
    }
  return fatJet1PNetXbb_;
}

const float &hhtree::fatJet1PNetQCDb() 
{
  if(not fatJet1PNetQCDb_isLoaded)
    {
      if(fatJet1PNetQCDb_branch != 0) fatJet1PNetQCDb_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PNetQCDb_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PNetQCDb_isLoaded = true;
    }
  return fatJet1PNetQCDb_;
}

const float &hhtree::fatJet1PNetQCDbb() 
{
  if(not fatJet1PNetQCDbb_isLoaded)
    {
      if(fatJet1PNetQCDbb_branch != 0) fatJet1PNetQCDbb_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PNetQCDbb_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PNetQCDbb_isLoaded = true;
    }
  return fatJet1PNetQCDbb_;
}

const float &hhtree::fatJet1PNetQCDc() 
{
  if(not fatJet1PNetQCDc_isLoaded)
    {
      if(fatJet1PNetQCDc_branch != 0) fatJet1PNetQCDc_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PNetQCDc_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PNetQCDc_isLoaded = true;
    }
  return fatJet1PNetQCDc_;
}

const float &hhtree::fatJet1PNetQCDcc() 
{
  if(not fatJet1PNetQCDcc_isLoaded)
    {
      if(fatJet1PNetQCDcc_branch != 0) fatJet1PNetQCDcc_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PNetQCDcc_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PNetQCDcc_isLoaded = true;
    }
  return fatJet1PNetQCDcc_;
}

const float &hhtree::fatJet1PNetQCDothers() 
{
  if(not fatJet1PNetQCDothers_isLoaded)
    {
      if(fatJet1PNetQCDothers_branch != 0) fatJet1PNetQCDothers_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PNetQCDothers_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PNetQCDothers_isLoaded = true;
    }
  return fatJet1PNetQCDothers_;
}

const float &hhtree::fatJet1Tau3OverTau2() 
{
  if(not fatJet1Tau3OverTau2_isLoaded)
    {
      if(fatJet1Tau3OverTau2_branch != 0) fatJet1Tau3OverTau2_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1Tau3OverTau2_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1Tau3OverTau2_isLoaded = true;
    }
  return fatJet1Tau3OverTau2_;
}

const int &hhtree::fatJet1GenMatchIndex() 
{
  if(not fatJet1GenMatchIndex_isLoaded)
    {
      if(fatJet1GenMatchIndex_branch != 0) fatJet1GenMatchIndex_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1GenMatchIndex_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1GenMatchIndex_isLoaded = true;
    }
  return fatJet1GenMatchIndex_;
}

const bool &hhtree::fatJet1HasMuon() 
{
  if(not fatJet1HasMuon_isLoaded)
    {
      if(fatJet1HasMuon_branch != 0) fatJet1HasMuon_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1HasMuon_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1HasMuon_isLoaded = true;
    }
  return fatJet1HasMuon_;
}

const bool &hhtree::fatJet1HasElectron() 
{
  if(not fatJet1HasElectron_isLoaded)
    {
      if(fatJet1HasElectron_branch != 0) fatJet1HasElectron_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1HasElectron_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1HasElectron_isLoaded = true;
    }
  return fatJet1HasElectron_;
}

const bool &hhtree::fatJet1HasBJetCSVLoose() 
{
  if(not fatJet1HasBJetCSVLoose_isLoaded)
    {
      if(fatJet1HasBJetCSVLoose_branch != 0) fatJet1HasBJetCSVLoose_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1HasBJetCSVLoose_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1HasBJetCSVLoose_isLoaded = true;
    }
  return fatJet1HasBJetCSVLoose_;
}

const bool &hhtree::fatJet1HasBJetCSVMedium() 
{
  if(not fatJet1HasBJetCSVMedium_isLoaded)
    {
      if(fatJet1HasBJetCSVMedium_branch != 0) fatJet1HasBJetCSVMedium_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1HasBJetCSVMedium_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1HasBJetCSVMedium_isLoaded = true;
    }
  return fatJet1HasBJetCSVMedium_;
}

const bool &hhtree::fatJet1HasBJetCSVTight() 
{
  if(not fatJet1HasBJetCSVTight_isLoaded)
    {
      if(fatJet1HasBJetCSVTight_branch != 0) fatJet1HasBJetCSVTight_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1HasBJetCSVTight_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1HasBJetCSVTight_isLoaded = true;
    }
  return fatJet1HasBJetCSVTight_;
}

const bool &hhtree::fatJet1OppositeHemisphereHasBJet() 
{
  if(not fatJet1OppositeHemisphereHasBJet_isLoaded)
    {
      if(fatJet1OppositeHemisphereHasBJet_branch != 0) fatJet1OppositeHemisphereHasBJet_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1OppositeHemisphereHasBJet_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1OppositeHemisphereHasBJet_isLoaded = true;
    }
  return fatJet1OppositeHemisphereHasBJet_;
}

const float &hhtree::fatJet1PtOverMHH() 
{
  if(not fatJet1PtOverMHH_isLoaded)
    {
      if(fatJet1PtOverMHH_branch != 0) fatJet1PtOverMHH_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PtOverMHH_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PtOverMHH_isLoaded = true;
    }
  return fatJet1PtOverMHH_;
}

const float &hhtree::fatJet1PtOverMSD() 
{
  if(not fatJet1PtOverMSD_isLoaded)
    {
      if(fatJet1PtOverMSD_branch != 0) fatJet1PtOverMSD_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PtOverMSD_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PtOverMSD_isLoaded = true;
    }
  return fatJet1PtOverMSD_;
}

const float &hhtree::fatJet1PtOverMRegressed() 
{
  if(not fatJet1PtOverMRegressed_isLoaded)
    {
      if(fatJet1PtOverMRegressed_branch != 0) fatJet1PtOverMRegressed_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1PtOverMRegressed_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1PtOverMRegressed_isLoaded = true;
    }
  return fatJet1PtOverMRegressed_;
}

const float &hhtree::fatJet1MassSD_JMS_Down() 
{
  if(not fatJet1MassSD_JMS_Down_isLoaded)
    {
      if(fatJet1MassSD_JMS_Down_branch != 0) fatJet1MassSD_JMS_Down_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassSD_JMS_Down_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassSD_JMS_Down_isLoaded = true;
    }
  return fatJet1MassSD_JMS_Down_;
}

const float &hhtree::fatJet1MassSD_JMS_Up() 
{
  if(not fatJet1MassSD_JMS_Up_isLoaded)
    {
      if(fatJet1MassSD_JMS_Up_branch != 0) fatJet1MassSD_JMS_Up_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassSD_JMS_Up_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassSD_JMS_Up_isLoaded = true;
    }
  return fatJet1MassSD_JMS_Up_;
}

const float &hhtree::fatJet1MassSD_JMR_Down() 
{
  if(not fatJet1MassSD_JMR_Down_isLoaded)
    {
      if(fatJet1MassSD_JMR_Down_branch != 0) fatJet1MassSD_JMR_Down_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassSD_JMR_Down_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassSD_JMR_Down_isLoaded = true;
    }
  return fatJet1MassSD_JMR_Down_;
}

const float &hhtree::fatJet1MassSD_JMR_Up() 
{
  if(not fatJet1MassSD_JMR_Up_isLoaded)
    {
      if(fatJet1MassSD_JMR_Up_branch != 0) fatJet1MassSD_JMR_Up_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassSD_JMR_Up_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassSD_JMR_Up_isLoaded = true;
    }
  return fatJet1MassSD_JMR_Up_;
}

const float &hhtree::fatJet1MassRegressed_JMS_Down() 
{
  if(not fatJet1MassRegressed_JMS_Down_isLoaded)
    {
      if(fatJet1MassRegressed_JMS_Down_branch != 0) fatJet1MassRegressed_JMS_Down_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassRegressed_JMS_Down_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassRegressed_JMS_Down_isLoaded = true;
    }
  return fatJet1MassRegressed_JMS_Down_;
}

const float &hhtree::fatJet1MassRegressed_JMS_Up() 
{
  if(not fatJet1MassRegressed_JMS_Up_isLoaded)
    {
      if(fatJet1MassRegressed_JMS_Up_branch != 0) fatJet1MassRegressed_JMS_Up_branch->GetEntry(index);
      else
	{
	  printf("branch fatJet1MassRegressed_JMS_Up_branch does not exist!\n");
	  //exit(1);
	}
      fatJet1MassRegressed_JMS_Up_isLoaded = true;
    }
  return fatJet1MassRegressed_JMS_Up_;
}

const float &hhtree::fatJet1MassRegressed_JMR_Down() 
{
 if(not fatJet1MassRegressed_JMR_Down_isLoaded)
 {
   if(fatJet1MassRegressed_JMR_Down_branch != 0) fatJet1MassRegressed_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1MassRegressed_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet1MassRegressed_JMR_Down_isLoaded = true;
 }
 return fatJet1MassRegressed_JMR_Down_;
}

const float &hhtree::fatJet1MassRegressed_JMR_Up() 
{
 if(not fatJet1MassRegressed_JMR_Up_isLoaded)
 {
   if(fatJet1MassRegressed_JMR_Up_branch != 0) fatJet1MassRegressed_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1MassRegressed_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet1MassRegressed_JMR_Up_isLoaded = true;
 }
 return fatJet1MassRegressed_JMR_Up_;
}

const float &hhtree::fatJet1PtOverMHH_JMS_Down() 
{
 if(not fatJet1PtOverMHH_JMS_Down_isLoaded)
 {
   if(fatJet1PtOverMHH_JMS_Down_branch != 0) fatJet1PtOverMHH_JMS_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JMS_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JMS_Down_isLoaded = true;
 }
 return fatJet1PtOverMHH_JMS_Down_;
}

const float &hhtree::fatJet1PtOverMHH_JMS_Up() 
{
 if(not fatJet1PtOverMHH_JMS_Up_isLoaded)
 {
   if(fatJet1PtOverMHH_JMS_Up_branch != 0) fatJet1PtOverMHH_JMS_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JMS_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JMS_Up_isLoaded = true;
 }
 return fatJet1PtOverMHH_JMS_Up_;
}

const float &hhtree::fatJet1PtOverMHH_JMR_Down() 
{
 if(not fatJet1PtOverMHH_JMR_Down_isLoaded)
 {
   if(fatJet1PtOverMHH_JMR_Down_branch != 0) fatJet1PtOverMHH_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JMR_Down_isLoaded = true;
 }
 return fatJet1PtOverMHH_JMR_Down_;
}

const float &hhtree::fatJet1PtOverMHH_JMR_Up() 
{
 if(not fatJet1PtOverMHH_JMR_Up_isLoaded)
 {
   if(fatJet1PtOverMHH_JMR_Up_branch != 0) fatJet1PtOverMHH_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JMR_Up_isLoaded = true;
 }
 return fatJet1PtOverMHH_JMR_Up_;
}

const float &hhtree::fatJet1Pt_JERUp() 
{
 if(not fatJet1Pt_JERUp_isLoaded)
 {
   if(fatJet1Pt_JERUp_branch != 0) fatJet1Pt_JERUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1Pt_JERUp_branch does not exist!\n");
//exit(1);
   }
   fatJet1Pt_JERUp_isLoaded = true;
 }
 return fatJet1Pt_JERUp_;
}

const float &hhtree::fatJet1PtOverMHH_JERUp() 
{
 if(not fatJet1PtOverMHH_JERUp_isLoaded)
 {
   if(fatJet1PtOverMHH_JERUp_branch != 0) fatJet1PtOverMHH_JERUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JERUp_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JERUp_isLoaded = true;
 }
 return fatJet1PtOverMHH_JERUp_;
}

const float &hhtree::fatJet1Pt_JERDown() 
{
 if(not fatJet1Pt_JERDown_isLoaded)
 {
   if(fatJet1Pt_JERDown_branch != 0) fatJet1Pt_JERDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1Pt_JERDown_branch does not exist!\n");
//exit(1);
   }
   fatJet1Pt_JERDown_isLoaded = true;
 }
 return fatJet1Pt_JERDown_;
}

const float &hhtree::fatJet1PtOverMHH_JERDown() 
{
 if(not fatJet1PtOverMHH_JERDown_isLoaded)
 {
   if(fatJet1PtOverMHH_JERDown_branch != 0) fatJet1PtOverMHH_JERDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JERDown_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JERDown_isLoaded = true;
 }
 return fatJet1PtOverMHH_JERDown_;
}

const float &hhtree::fatJet1Pt_JESUp() 
{
 if(not fatJet1Pt_JESUp_isLoaded)
 {
   if(fatJet1Pt_JESUp_branch != 0) fatJet1Pt_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1Pt_JESUp_branch does not exist!\n");
//exit(1);
   }
   fatJet1Pt_JESUp_isLoaded = true;
 }
 return fatJet1Pt_JESUp_;
}

const float &hhtree::fatJet1PtOverMHH_JESUp() 
{
 if(not fatJet1PtOverMHH_JESUp_isLoaded)
 {
   if(fatJet1PtOverMHH_JESUp_branch != 0) fatJet1PtOverMHH_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JESUp_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JESUp_isLoaded = true;
 }
 return fatJet1PtOverMHH_JESUp_;
}

const float &hhtree::fatJet1Pt_JESDown() 
{
 if(not fatJet1Pt_JESDown_isLoaded)
 {
   if(fatJet1Pt_JESDown_branch != 0) fatJet1Pt_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1Pt_JESDown_branch does not exist!\n");
//exit(1);
   }
   fatJet1Pt_JESDown_isLoaded = true;
 }
 return fatJet1Pt_JESDown_;
}

const float &hhtree::fatJet1PtOverMHH_JESDown() 
{
 if(not fatJet1PtOverMHH_JESDown_isLoaded)
 {
   if(fatJet1PtOverMHH_JESDown_branch != 0) fatJet1PtOverMHH_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet1PtOverMHH_JESDown_branch does not exist!\n");
//exit(1);
   }
   fatJet1PtOverMHH_JESDown_isLoaded = true;
 }
 return fatJet1PtOverMHH_JESDown_;
}

const float &hhtree::fatJet1Pt_JESUp_Abs(){
	if(not fatJet1Pt_JESUp_Abs_isLoaded){
		if(fatJet1Pt_JESUp_Abs_branch != 0) fatJet1Pt_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_Abs_branch does not exist!\n");}
		fatJet1Pt_JESUp_Abs_isLoaded = true;
}
return fatJet1Pt_JESUp_Abs_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_Abs(){
	if(not fatJet1PtOverMHH_JESUp_Abs_isLoaded){
		if(fatJet1PtOverMHH_JESUp_Abs_branch != 0) fatJet1PtOverMHH_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_Abs_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_Abs_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_Abs_;
}
const float &hhtree::fatJet1Pt_JESDown_Abs(){
	if(not fatJet1Pt_JESDown_Abs_isLoaded){
		if(fatJet1Pt_JESDown_Abs_branch != 0) fatJet1Pt_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_Abs_branch does not exist!\n");}
		fatJet1Pt_JESDown_Abs_isLoaded = true;
}
return fatJet1Pt_JESDown_Abs_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_Abs(){
	if(not fatJet1PtOverMHH_JESDown_Abs_isLoaded){
		if(fatJet1PtOverMHH_JESDown_Abs_branch != 0) fatJet1PtOverMHH_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_Abs_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_Abs_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_Abs_;
}
const float &hhtree::fatJet1Pt_JESUp_Abs_2016(){
	if(not fatJet1Pt_JESUp_Abs_2016_isLoaded){
		if(fatJet1Pt_JESUp_Abs_2016_branch != 0) fatJet1Pt_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_Abs_2016_branch does not exist!\n");}
		fatJet1Pt_JESUp_Abs_2016_isLoaded = true;
}
return fatJet1Pt_JESUp_Abs_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_Abs_2016(){
	if(not fatJet1PtOverMHH_JESUp_Abs_2016_isLoaded){
		if(fatJet1PtOverMHH_JESUp_Abs_2016_branch != 0) fatJet1PtOverMHH_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_Abs_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_Abs_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_Abs_2016_;
}
const float &hhtree::fatJet1Pt_JESDowb_Abs_2016(){
	if(not fatJet1Pt_JESDowb_Abs_2016_isLoaded){
		if(fatJet1Pt_JESDowb_Abs_2016_branch != 0) fatJet1Pt_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDowb_Abs_2016_branch does not exist!\n");}
		fatJet1Pt_JESDowb_Abs_2016_isLoaded = true;
}
return fatJet1Pt_JESDowb_Abs_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESDowb_Abs_2016(){
	if(not fatJet1PtOverMHH_JESDowb_Abs_2016_isLoaded){
		if(fatJet1PtOverMHH_JESDowb_Abs_2016_branch != 0) fatJet1PtOverMHH_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDowb_Abs_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDowb_Abs_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESDowb_Abs_2016_;
}
const float &hhtree::fatJet1Pt_JESUp_BBEC1(){
	if(not fatJet1Pt_JESUp_BBEC1_isLoaded){
		if(fatJet1Pt_JESUp_BBEC1_branch != 0) fatJet1Pt_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_BBEC1_branch does not exist!\n");}
		fatJet1Pt_JESUp_BBEC1_isLoaded = true;
}
return fatJet1Pt_JESUp_BBEC1_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_BBEC1(){
	if(not fatJet1PtOverMHH_JESUp_BBEC1_isLoaded){
		if(fatJet1PtOverMHH_JESUp_BBEC1_branch != 0) fatJet1PtOverMHH_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_BBEC1_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_BBEC1_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_BBEC1_;
}
const float &hhtree::fatJet1Pt_JESDown_BBEC1(){
	if(not fatJet1Pt_JESDown_BBEC1_isLoaded){
		if(fatJet1Pt_JESDown_BBEC1_branch != 0) fatJet1Pt_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_BBEC1_branch does not exist!\n");}
		fatJet1Pt_JESDown_BBEC1_isLoaded = true;
}
return fatJet1Pt_JESDown_BBEC1_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_BBEC1(){
	if(not fatJet1PtOverMHH_JESDown_BBEC1_isLoaded){
		if(fatJet1PtOverMHH_JESDown_BBEC1_branch != 0) fatJet1PtOverMHH_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_BBEC1_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_BBEC1_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_BBEC1_;
}
const float &hhtree::fatJet1Pt_JESUp_BBEC1_2016(){
	if(not fatJet1Pt_JESUp_BBEC1_2016_isLoaded){
		if(fatJet1Pt_JESUp_BBEC1_2016_branch != 0) fatJet1Pt_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_BBEC1_2016_branch does not exist!\n");}
		fatJet1Pt_JESUp_BBEC1_2016_isLoaded = true;
}
return fatJet1Pt_JESUp_BBEC1_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_BBEC1_2016(){
	if(not fatJet1PtOverMHH_JESUp_BBEC1_2016_isLoaded){
		if(fatJet1PtOverMHH_JESUp_BBEC1_2016_branch != 0) fatJet1PtOverMHH_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_BBEC1_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_BBEC1_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_BBEC1_2016_;
}
const float &hhtree::fatJet1Pt_JESDown_BBEC1_2016(){
	if(not fatJet1Pt_JESDown_BBEC1_2016_isLoaded){
		if(fatJet1Pt_JESDown_BBEC1_2016_branch != 0) fatJet1Pt_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_BBEC1_2016_branch does not exist!\n");}
		fatJet1Pt_JESDown_BBEC1_2016_isLoaded = true;
}
return fatJet1Pt_JESDown_BBEC1_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_BBEC1_2016(){
	if(not fatJet1PtOverMHH_JESDown_BBEC1_2016_isLoaded){
		if(fatJet1PtOverMHH_JESDown_BBEC1_2016_branch != 0) fatJet1PtOverMHH_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_BBEC1_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_BBEC1_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_BBEC1_2016_;
}
const float &hhtree::fatJet1Pt_JESUp_EC2(){
	if(not fatJet1Pt_JESUp_EC2_isLoaded){
		if(fatJet1Pt_JESUp_EC2_branch != 0) fatJet1Pt_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_EC2_branch does not exist!\n");}
		fatJet1Pt_JESUp_EC2_isLoaded = true;
}
return fatJet1Pt_JESUp_EC2_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_EC2(){
	if(not fatJet1PtOverMHH_JESUp_EC2_isLoaded){
		if(fatJet1PtOverMHH_JESUp_EC2_branch != 0) fatJet1PtOverMHH_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_EC2_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_EC2_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_EC2_;
}
const float &hhtree::fatJet1Pt_JESDown_EC2(){
	if(not fatJet1Pt_JESDown_EC2_isLoaded){
		if(fatJet1Pt_JESDown_EC2_branch != 0) fatJet1Pt_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_EC2_branch does not exist!\n");}
		fatJet1Pt_JESDown_EC2_isLoaded = true;
}
return fatJet1Pt_JESDown_EC2_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_EC2(){
	if(not fatJet1PtOverMHH_JESDown_EC2_isLoaded){
		if(fatJet1PtOverMHH_JESDown_EC2_branch != 0) fatJet1PtOverMHH_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_EC2_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_EC2_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_EC2_;
}
const float &hhtree::fatJet1Pt_JESUp_EC2_2016(){
	if(not fatJet1Pt_JESUp_EC2_2016_isLoaded){
		if(fatJet1Pt_JESUp_EC2_2016_branch != 0) fatJet1Pt_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_EC2_2016_branch does not exist!\n");}
		fatJet1Pt_JESUp_EC2_2016_isLoaded = true;
}
return fatJet1Pt_JESUp_EC2_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_EC2_2016(){
	if(not fatJet1PtOverMHH_JESUp_EC2_2016_isLoaded){
		if(fatJet1PtOverMHH_JESUp_EC2_2016_branch != 0) fatJet1PtOverMHH_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_EC2_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_EC2_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_EC2_2016_;
}
const float &hhtree::fatJet1Pt_JESDown_EC2_2016(){
	if(not fatJet1Pt_JESDown_EC2_2016_isLoaded){
		if(fatJet1Pt_JESDown_EC2_2016_branch != 0) fatJet1Pt_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_EC2_2016_branch does not exist!\n");}
		fatJet1Pt_JESDown_EC2_2016_isLoaded = true;
}
return fatJet1Pt_JESDown_EC2_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_EC2_2016(){
	if(not fatJet1PtOverMHH_JESDown_EC2_2016_isLoaded){
		if(fatJet1PtOverMHH_JESDown_EC2_2016_branch != 0) fatJet1PtOverMHH_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_EC2_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_EC2_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_EC2_2016_;
}
const float &hhtree::fatJet1Pt_JESUp_FlavQCD(){
	if(not fatJet1Pt_JESUp_FlavQCD_isLoaded){
		if(fatJet1Pt_JESUp_FlavQCD_branch != 0) fatJet1Pt_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_FlavQCD_branch does not exist!\n");}
		fatJet1Pt_JESUp_FlavQCD_isLoaded = true;
}
return fatJet1Pt_JESUp_FlavQCD_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_FlavQCD(){
	if(not fatJet1PtOverMHH_JESUp_FlavQCD_isLoaded){
		if(fatJet1PtOverMHH_JESUp_FlavQCD_branch != 0) fatJet1PtOverMHH_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_FlavQCD_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_FlavQCD_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_FlavQCD_;
}
const float &hhtree::fatJet1Pt_JESDown_FlavQCD(){
	if(not fatJet1Pt_JESDown_FlavQCD_isLoaded){
		if(fatJet1Pt_JESDown_FlavQCD_branch != 0) fatJet1Pt_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_FlavQCD_branch does not exist!\n");}
		fatJet1Pt_JESDown_FlavQCD_isLoaded = true;
}
return fatJet1Pt_JESDown_FlavQCD_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_FlavQCD(){
	if(not fatJet1PtOverMHH_JESDown_FlavQCD_isLoaded){
		if(fatJet1PtOverMHH_JESDown_FlavQCD_branch != 0) fatJet1PtOverMHH_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_FlavQCD_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_FlavQCD_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_FlavQCD_;
}
const float &hhtree::fatJet1Pt_JESUp_HF(){
	if(not fatJet1Pt_JESUp_HF_isLoaded){
		if(fatJet1Pt_JESUp_HF_branch != 0) fatJet1Pt_JESUp_HF_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_HF_branch does not exist!\n");}
		fatJet1Pt_JESUp_HF_isLoaded = true;
}
return fatJet1Pt_JESUp_HF_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_HF(){
	if(not fatJet1PtOverMHH_JESUp_HF_isLoaded){
		if(fatJet1PtOverMHH_JESUp_HF_branch != 0) fatJet1PtOverMHH_JESUp_HF_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_HF_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_HF_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_HF_;
}
const float &hhtree::fatJet1Pt_JESDown_HF(){
	if(not fatJet1Pt_JESDown_HF_isLoaded){
		if(fatJet1Pt_JESDown_HF_branch != 0) fatJet1Pt_JESDown_HF_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_HF_branch does not exist!\n");}
		fatJet1Pt_JESDown_HF_isLoaded = true;
}
return fatJet1Pt_JESDown_HF_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_HF(){
	if(not fatJet1PtOverMHH_JESDown_HF_isLoaded){
		if(fatJet1PtOverMHH_JESDown_HF_branch != 0) fatJet1PtOverMHH_JESDown_HF_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_HF_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_HF_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_HF_;
}
const float &hhtree::fatJet1Pt_JESUp_HF_2016(){
	if(not fatJet1Pt_JESUp_HF_2016_isLoaded){
		if(fatJet1Pt_JESUp_HF_2016_branch != 0) fatJet1Pt_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_HF_2016_branch does not exist!\n");}
		fatJet1Pt_JESUp_HF_2016_isLoaded = true;
}
return fatJet1Pt_JESUp_HF_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_HF_2016(){
	if(not fatJet1PtOverMHH_JESUp_HF_2016_isLoaded){
		if(fatJet1PtOverMHH_JESUp_HF_2016_branch != 0) fatJet1PtOverMHH_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_HF_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_HF_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_HF_2016_;
}
const float &hhtree::fatJet1Pt_JESDown_HF_2016(){
	if(not fatJet1Pt_JESDown_HF_2016_isLoaded){
		if(fatJet1Pt_JESDown_HF_2016_branch != 0) fatJet1Pt_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_HF_2016_branch does not exist!\n");}
		fatJet1Pt_JESDown_HF_2016_isLoaded = true;
}
return fatJet1Pt_JESDown_HF_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_HF_2016(){
	if(not fatJet1PtOverMHH_JESDown_HF_2016_isLoaded){
		if(fatJet1PtOverMHH_JESDown_HF_2016_branch != 0) fatJet1PtOverMHH_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_HF_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_HF_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_HF_2016_;
}
const float &hhtree::fatJet1Pt_JESUp_RelBal(){
	if(not fatJet1Pt_JESUp_RelBal_isLoaded){
		if(fatJet1Pt_JESUp_RelBal_branch != 0) fatJet1Pt_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_RelBal_branch does not exist!\n");}
		fatJet1Pt_JESUp_RelBal_isLoaded = true;
}
return fatJet1Pt_JESUp_RelBal_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_RelBal(){
	if(not fatJet1PtOverMHH_JESUp_RelBal_isLoaded){
		if(fatJet1PtOverMHH_JESUp_RelBal_branch != 0) fatJet1PtOverMHH_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_RelBal_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_RelBal_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_RelBal_;
}
const float &hhtree::fatJet1Pt_JESDown_RelBal(){
	if(not fatJet1Pt_JESDown_RelBal_isLoaded){
		if(fatJet1Pt_JESDown_RelBal_branch != 0) fatJet1Pt_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_RelBal_branch does not exist!\n");}
		fatJet1Pt_JESDown_RelBal_isLoaded = true;
}
return fatJet1Pt_JESDown_RelBal_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_RelBal(){
	if(not fatJet1PtOverMHH_JESDown_RelBal_isLoaded){
		if(fatJet1PtOverMHH_JESDown_RelBal_branch != 0) fatJet1PtOverMHH_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_RelBal_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_RelBal_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_RelBal_;
}
const float &hhtree::fatJet1Pt_JESUp_RelSample_2016(){
	if(not fatJet1Pt_JESUp_RelSample_2016_isLoaded){
		if(fatJet1Pt_JESUp_RelSample_2016_branch != 0) fatJet1Pt_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_RelSample_2016_branch does not exist!\n");}
		fatJet1Pt_JESUp_RelSample_2016_isLoaded = true;
}
return fatJet1Pt_JESUp_RelSample_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_RelSample_2016(){
	if(not fatJet1PtOverMHH_JESUp_RelSample_2016_isLoaded){
		if(fatJet1PtOverMHH_JESUp_RelSample_2016_branch != 0) fatJet1PtOverMHH_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_RelSample_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_RelSample_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_RelSample_2016_;
}
const float &hhtree::fatJet1Pt_JESDown_RelSample_2016(){
	if(not fatJet1Pt_JESDown_RelSample_2016_isLoaded){
		if(fatJet1Pt_JESDown_RelSample_2016_branch != 0) fatJet1Pt_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_RelSample_2016_branch does not exist!\n");}
		fatJet1Pt_JESDown_RelSample_2016_isLoaded = true;
}
return fatJet1Pt_JESDown_RelSample_2016_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_RelSample_2016(){
	if(not fatJet1PtOverMHH_JESDown_RelSample_2016_isLoaded){
		if(fatJet1PtOverMHH_JESDown_RelSample_2016_branch != 0) fatJet1PtOverMHH_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_RelSample_2016_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_RelSample_2016_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_RelSample_2016_;
}



const float &hhtree::fatJet2Pt() 
{
 if(not fatJet2Pt_isLoaded)
 {
   if(fatJet2Pt_branch != 0) fatJet2Pt_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Pt_branch does not exist!\n");
//exit(1);
   }
   fatJet2Pt_isLoaded = true;
 }
 return fatJet2Pt_;
}

const float &hhtree::fatJet2Eta() 
{
 if(not fatJet2Eta_isLoaded)
 {
   if(fatJet2Eta_branch != 0) fatJet2Eta_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Eta_branch does not exist!\n");
//exit(1);
   }
   fatJet2Eta_isLoaded = true;
 }
 return fatJet2Eta_;
}

const float &hhtree::fatJet2Phi() 
{
 if(not fatJet2Phi_isLoaded)
 {
   if(fatJet2Phi_branch != 0) fatJet2Phi_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Phi_branch does not exist!\n");
//exit(1);
   }
   fatJet2Phi_isLoaded = true;
 }
 return fatJet2Phi_;
}

const float &hhtree::fatJet2Mass() 
{
 if(not fatJet2Mass_isLoaded)
 {
   if(fatJet2Mass_branch != 0) fatJet2Mass_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Mass_branch does not exist!\n");
//exit(1);
   }
   fatJet2Mass_isLoaded = true;
 }
 return fatJet2Mass_;
}

const float &hhtree::fatJet2MassSD() 
{
 if(not fatJet2MassSD_isLoaded)
 {
   if(fatJet2MassSD_branch != 0) fatJet2MassSD_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassSD_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassSD_isLoaded = true;
 }
 return fatJet2MassSD_;
}

const float &hhtree::fatJet2MassSD_UnCorrected() 
{
 if(not fatJet2MassSD_UnCorrected_isLoaded)
 {
   if(fatJet2MassSD_UnCorrected_branch != 0) fatJet2MassSD_UnCorrected_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassSD_UnCorrected_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassSD_UnCorrected_isLoaded = true;
 }
 return fatJet2MassSD_UnCorrected_;
}

const float &hhtree::fatJet2MassRegressed() 
{
 if(not fatJet2MassRegressed_isLoaded)
 {
   if(fatJet2MassRegressed_branch != 0) fatJet2MassRegressed_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassRegressed_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassRegressed_isLoaded = true;
 }
 return fatJet2MassRegressed_;
}

const float &hhtree::fatJet2MassRegressed_UnCorrected() 
{
 if(not fatJet2MassRegressed_UnCorrected_isLoaded)
 {
   if(fatJet2MassRegressed_UnCorrected_branch != 0) fatJet2MassRegressed_UnCorrected_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassRegressed_UnCorrected_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassRegressed_UnCorrected_isLoaded = true;
 }
 return fatJet2MassRegressed_UnCorrected_;
}

const float &hhtree::fatJet2PNetXbb() 
{
 if(not fatJet2PNetXbb_isLoaded)
 {
   if(fatJet2PNetXbb_branch != 0) fatJet2PNetXbb_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PNetXbb_branch does not exist!\n");
//exit(1);
   }
   fatJet2PNetXbb_isLoaded = true;
 }
 return fatJet2PNetXbb_;
}

const float &hhtree::fatJet2PNetQCDb() 
{
 if(not fatJet2PNetQCDb_isLoaded)
 {
   if(fatJet2PNetQCDb_branch != 0) fatJet2PNetQCDb_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PNetQCDb_branch does not exist!\n");
//exit(1);
   }
   fatJet2PNetQCDb_isLoaded = true;
 }
 return fatJet2PNetQCDb_;
}

const float &hhtree::fatJet2PNetQCDbb() 
{
 if(not fatJet2PNetQCDbb_isLoaded)
 {
   if(fatJet2PNetQCDbb_branch != 0) fatJet2PNetQCDbb_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PNetQCDbb_branch does not exist!\n");
//exit(1);
   }
   fatJet2PNetQCDbb_isLoaded = true;
 }
 return fatJet2PNetQCDbb_;
}

const float &hhtree::fatJet2PNetQCDc() 
{
 if(not fatJet2PNetQCDc_isLoaded)
 {
   if(fatJet2PNetQCDc_branch != 0) fatJet2PNetQCDc_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PNetQCDc_branch does not exist!\n");
//exit(1);
   }
   fatJet2PNetQCDc_isLoaded = true;
 }
 return fatJet2PNetQCDc_;
}

const float &hhtree::fatJet2PNetQCDcc() 
{
 if(not fatJet2PNetQCDcc_isLoaded)
 {
   if(fatJet2PNetQCDcc_branch != 0) fatJet2PNetQCDcc_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PNetQCDcc_branch does not exist!\n");
//exit(1);
   }
   fatJet2PNetQCDcc_isLoaded = true;
 }
 return fatJet2PNetQCDcc_;
}

const float &hhtree::fatJet2PNetQCDothers() 
{
 if(not fatJet2PNetQCDothers_isLoaded)
 {
   if(fatJet2PNetQCDothers_branch != 0) fatJet2PNetQCDothers_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PNetQCDothers_branch does not exist!\n");
//exit(1);
   }
   fatJet2PNetQCDothers_isLoaded = true;
 }
 return fatJet2PNetQCDothers_;
}

const float &hhtree::fatJet2Tau3OverTau2() 
{
 if(not fatJet2Tau3OverTau2_isLoaded)
 {
   if(fatJet2Tau3OverTau2_branch != 0) fatJet2Tau3OverTau2_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Tau3OverTau2_branch does not exist!\n");
//exit(1);
   }
   fatJet2Tau3OverTau2_isLoaded = true;
 }
 return fatJet2Tau3OverTau2_;
}

const int &hhtree::fatJet2GenMatchIndex() 
{
 if(not fatJet2GenMatchIndex_isLoaded)
 {
   if(fatJet2GenMatchIndex_branch != 0) fatJet2GenMatchIndex_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2GenMatchIndex_branch does not exist!\n");
//exit(1);
   }
   fatJet2GenMatchIndex_isLoaded = true;
 }
 return fatJet2GenMatchIndex_;
}

const bool &hhtree::fatJet2HasMuon() 
{
 if(not fatJet2HasMuon_isLoaded)
 {
   if(fatJet2HasMuon_branch != 0) fatJet2HasMuon_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2HasMuon_branch does not exist!\n");
//exit(1);
   }
   fatJet2HasMuon_isLoaded = true;
 }
 return fatJet2HasMuon_;
}

const bool &hhtree::fatJet2HasElectron() 
{
 if(not fatJet2HasElectron_isLoaded)
 {
   if(fatJet2HasElectron_branch != 0) fatJet2HasElectron_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2HasElectron_branch does not exist!\n");
//exit(1);
   }
   fatJet2HasElectron_isLoaded = true;
 }
 return fatJet2HasElectron_;
}

const bool &hhtree::fatJet2HasBJetCSVLoose() 
{
 if(not fatJet2HasBJetCSVLoose_isLoaded)
 {
   if(fatJet2HasBJetCSVLoose_branch != 0) fatJet2HasBJetCSVLoose_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2HasBJetCSVLoose_branch does not exist!\n");
//exit(1);
   }
   fatJet2HasBJetCSVLoose_isLoaded = true;
 }
 return fatJet2HasBJetCSVLoose_;
}

const bool &hhtree::fatJet2HasBJetCSVMedium() 
{
 if(not fatJet2HasBJetCSVMedium_isLoaded)
 {
   if(fatJet2HasBJetCSVMedium_branch != 0) fatJet2HasBJetCSVMedium_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2HasBJetCSVMedium_branch does not exist!\n");
//exit(1);
   }
   fatJet2HasBJetCSVMedium_isLoaded = true;
 }
 return fatJet2HasBJetCSVMedium_;
}

const bool &hhtree::fatJet2HasBJetCSVTight() 
{
 if(not fatJet2HasBJetCSVTight_isLoaded)
 {
   if(fatJet2HasBJetCSVTight_branch != 0) fatJet2HasBJetCSVTight_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2HasBJetCSVTight_branch does not exist!\n");
//exit(1);
   }
   fatJet2HasBJetCSVTight_isLoaded = true;
 }
 return fatJet2HasBJetCSVTight_;
}

const bool &hhtree::fatJet2OppositeHemisphereHasBJet() 
{
 if(not fatJet2OppositeHemisphereHasBJet_isLoaded)
 {
   if(fatJet2OppositeHemisphereHasBJet_branch != 0) fatJet2OppositeHemisphereHasBJet_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2OppositeHemisphereHasBJet_branch does not exist!\n");
//exit(1);
   }
   fatJet2OppositeHemisphereHasBJet_isLoaded = true;
 }
 return fatJet2OppositeHemisphereHasBJet_;
}

const float &hhtree::fatJet2PtOverMHH() 
{
 if(not fatJet2PtOverMHH_isLoaded)
 {
   if(fatJet2PtOverMHH_branch != 0) fatJet2PtOverMHH_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_isLoaded = true;
 }
 return fatJet2PtOverMHH_;
}

const float &hhtree::fatJet2PtOverMSD() 
{
 if(not fatJet2PtOverMSD_isLoaded)
 {
   if(fatJet2PtOverMSD_branch != 0) fatJet2PtOverMSD_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMSD_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMSD_isLoaded = true;
 }
 return fatJet2PtOverMSD_;
}

const float &hhtree::fatJet2PtOverMRegressed() 
{
 if(not fatJet2PtOverMRegressed_isLoaded)
 {
   if(fatJet2PtOverMRegressed_branch != 0) fatJet2PtOverMRegressed_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMRegressed_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMRegressed_isLoaded = true;
 }
 return fatJet2PtOverMRegressed_;
}

const float &hhtree::fatJet2MassSD_JMS_Down() 
{
 if(not fatJet2MassSD_JMS_Down_isLoaded)
 {
   if(fatJet2MassSD_JMS_Down_branch != 0) fatJet2MassSD_JMS_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassSD_JMS_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassSD_JMS_Down_isLoaded = true;
 }
 return fatJet2MassSD_JMS_Down_;
}

const float &hhtree::fatJet2MassSD_JMS_Up() 
{
 if(not fatJet2MassSD_JMS_Up_isLoaded)
 {
   if(fatJet2MassSD_JMS_Up_branch != 0) fatJet2MassSD_JMS_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassSD_JMS_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassSD_JMS_Up_isLoaded = true;
 }
 return fatJet2MassSD_JMS_Up_;
}

const float &hhtree::fatJet2MassSD_JMR_Down() 
{
 if(not fatJet2MassSD_JMR_Down_isLoaded)
 {
   if(fatJet2MassSD_JMR_Down_branch != 0) fatJet2MassSD_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassSD_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassSD_JMR_Down_isLoaded = true;
 }
 return fatJet2MassSD_JMR_Down_;
}

const float &hhtree::fatJet2MassSD_JMR_Up() 
{
 if(not fatJet2MassSD_JMR_Up_isLoaded)
 {
   if(fatJet2MassSD_JMR_Up_branch != 0) fatJet2MassSD_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassSD_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassSD_JMR_Up_isLoaded = true;
 }
 return fatJet2MassSD_JMR_Up_;
}

const float &hhtree::fatJet2MassRegressed_JMS_Down() 
{
 if(not fatJet2MassRegressed_JMS_Down_isLoaded)
 {
   if(fatJet2MassRegressed_JMS_Down_branch != 0) fatJet2MassRegressed_JMS_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassRegressed_JMS_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassRegressed_JMS_Down_isLoaded = true;
 }
 return fatJet2MassRegressed_JMS_Down_;
}

const float &hhtree::fatJet2MassRegressed_JMS_Up() 
{
 if(not fatJet2MassRegressed_JMS_Up_isLoaded)
 {
   if(fatJet2MassRegressed_JMS_Up_branch != 0) fatJet2MassRegressed_JMS_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassRegressed_JMS_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassRegressed_JMS_Up_isLoaded = true;
 }
 return fatJet2MassRegressed_JMS_Up_;
}

const float &hhtree::fatJet2MassRegressed_JMR_Down() 
{
 if(not fatJet2MassRegressed_JMR_Down_isLoaded)
 {
   if(fatJet2MassRegressed_JMR_Down_branch != 0) fatJet2MassRegressed_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassRegressed_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassRegressed_JMR_Down_isLoaded = true;
 }
 return fatJet2MassRegressed_JMR_Down_;
}

const float &hhtree::fatJet2MassRegressed_JMR_Up() 
{
 if(not fatJet2MassRegressed_JMR_Up_isLoaded)
 {
   if(fatJet2MassRegressed_JMR_Up_branch != 0) fatJet2MassRegressed_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2MassRegressed_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet2MassRegressed_JMR_Up_isLoaded = true;
 }
 return fatJet2MassRegressed_JMR_Up_;
}

const float &hhtree::fatJet2PtOverMHH_JMS_Down() 
{
 if(not fatJet2PtOverMHH_JMS_Down_isLoaded)
 {
   if(fatJet2PtOverMHH_JMS_Down_branch != 0) fatJet2PtOverMHH_JMS_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JMS_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JMS_Down_isLoaded = true;
 }
 return fatJet2PtOverMHH_JMS_Down_;
}

const float &hhtree::fatJet2PtOverMHH_JMS_Up() 
{
 if(not fatJet2PtOverMHH_JMS_Up_isLoaded)
 {
   if(fatJet2PtOverMHH_JMS_Up_branch != 0) fatJet2PtOverMHH_JMS_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JMS_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JMS_Up_isLoaded = true;
 }
 return fatJet2PtOverMHH_JMS_Up_;
}

const float &hhtree::fatJet2PtOverMHH_JMR_Down() 
{
 if(not fatJet2PtOverMHH_JMR_Down_isLoaded)
 {
   if(fatJet2PtOverMHH_JMR_Down_branch != 0) fatJet2PtOverMHH_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JMR_Down_isLoaded = true;
 }
 return fatJet2PtOverMHH_JMR_Down_;
}

const float &hhtree::fatJet2PtOverMHH_JMR_Up() 
{
 if(not fatJet2PtOverMHH_JMR_Up_isLoaded)
 {
   if(fatJet2PtOverMHH_JMR_Up_branch != 0) fatJet2PtOverMHH_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JMR_Up_isLoaded = true;
 }
 return fatJet2PtOverMHH_JMR_Up_;
}

const float &hhtree::fatJet2Pt_JERUp() 
{
 if(not fatJet2Pt_JERUp_isLoaded)
 {
   if(fatJet2Pt_JERUp_branch != 0) fatJet2Pt_JERUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Pt_JERUp_branch does not exist!\n");
//exit(1);
   }
   fatJet2Pt_JERUp_isLoaded = true;
 }
 return fatJet2Pt_JERUp_;
}

const float &hhtree::fatJet2PtOverMHH_JERUp() 
{
 if(not fatJet2PtOverMHH_JERUp_isLoaded)
 {
   if(fatJet2PtOverMHH_JERUp_branch != 0) fatJet2PtOverMHH_JERUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JERUp_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JERUp_isLoaded = true;
 }
 return fatJet2PtOverMHH_JERUp_;
}

const float &hhtree::fatJet2Pt_JERDown() 
{
 if(not fatJet2Pt_JERDown_isLoaded)
 {
   if(fatJet2Pt_JERDown_branch != 0) fatJet2Pt_JERDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Pt_JERDown_branch does not exist!\n");
//exit(1);
   }
   fatJet2Pt_JERDown_isLoaded = true;
 }
 return fatJet2Pt_JERDown_;
}

const float &hhtree::fatJet2PtOverMHH_JERDown() 
{
 if(not fatJet2PtOverMHH_JERDown_isLoaded)
 {
   if(fatJet2PtOverMHH_JERDown_branch != 0) fatJet2PtOverMHH_JERDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JERDown_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JERDown_isLoaded = true;
 }
 return fatJet2PtOverMHH_JERDown_;
}

const float &hhtree::fatJet2Pt_JESUp() 
{
 if(not fatJet2Pt_JESUp_isLoaded)
 {
   if(fatJet2Pt_JESUp_branch != 0) fatJet2Pt_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Pt_JESUp_branch does not exist!\n");
//exit(1);
   }
   fatJet2Pt_JESUp_isLoaded = true;
 }
 return fatJet2Pt_JESUp_;
}

const float &hhtree::fatJet2PtOverMHH_JESUp() 
{
 if(not fatJet2PtOverMHH_JESUp_isLoaded)
 {
   if(fatJet2PtOverMHH_JESUp_branch != 0) fatJet2PtOverMHH_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JESUp_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JESUp_isLoaded = true;
 }
 return fatJet2PtOverMHH_JESUp_;
}

const float &hhtree::fatJet2Pt_JESDown() 
{
 if(not fatJet2Pt_JESDown_isLoaded)
 {
   if(fatJet2Pt_JESDown_branch != 0) fatJet2Pt_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2Pt_JESDown_branch does not exist!\n");
//exit(1);
   }
   fatJet2Pt_JESDown_isLoaded = true;
 }
 return fatJet2Pt_JESDown_;
}

const float &hhtree::fatJet2PtOverMHH_JESDown() 
{
 if(not fatJet2PtOverMHH_JESDown_isLoaded)
 {
   if(fatJet2PtOverMHH_JESDown_branch != 0) fatJet2PtOverMHH_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch fatJet2PtOverMHH_JESDown_branch does not exist!\n");
//exit(1);
   }
   fatJet2PtOverMHH_JESDown_isLoaded = true;
 }
 return fatJet2PtOverMHH_JESDown_;
}

const float &hhtree::fatJet2Pt_JESUp_Abs(){
	if(not fatJet2Pt_JESUp_Abs_isLoaded){
		if(fatJet2Pt_JESUp_Abs_branch != 0) fatJet2Pt_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_Abs_branch does not exist!\n");}
		fatJet2Pt_JESUp_Abs_isLoaded = true;
}
return fatJet2Pt_JESUp_Abs_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_Abs(){
	if(not fatJet2PtOverMHH_JESUp_Abs_isLoaded){
		if(fatJet2PtOverMHH_JESUp_Abs_branch != 0) fatJet2PtOverMHH_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_Abs_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_Abs_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_Abs_;
}
const float &hhtree::fatJet2Pt_JESDown_Abs(){
	if(not fatJet2Pt_JESDown_Abs_isLoaded){
		if(fatJet2Pt_JESDown_Abs_branch != 0) fatJet2Pt_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_Abs_branch does not exist!\n");}
		fatJet2Pt_JESDown_Abs_isLoaded = true;
}
return fatJet2Pt_JESDown_Abs_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_Abs(){
	if(not fatJet2PtOverMHH_JESDown_Abs_isLoaded){
		if(fatJet2PtOverMHH_JESDown_Abs_branch != 0) fatJet2PtOverMHH_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_Abs_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_Abs_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_Abs_;
}
const float &hhtree::fatJet2Pt_JESUp_Abs_2016(){
	if(not fatJet2Pt_JESUp_Abs_2016_isLoaded){
		if(fatJet2Pt_JESUp_Abs_2016_branch != 0) fatJet2Pt_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_Abs_2016_branch does not exist!\n");}
		fatJet2Pt_JESUp_Abs_2016_isLoaded = true;
}
return fatJet2Pt_JESUp_Abs_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_Abs_2016(){
	if(not fatJet2PtOverMHH_JESUp_Abs_2016_isLoaded){
		if(fatJet2PtOverMHH_JESUp_Abs_2016_branch != 0) fatJet2PtOverMHH_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_Abs_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_Abs_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_Abs_2016_;
}
const float &hhtree::fatJet2Pt_JESDowb_Abs_2016(){
	if(not fatJet2Pt_JESDowb_Abs_2016_isLoaded){
		if(fatJet2Pt_JESDowb_Abs_2016_branch != 0) fatJet2Pt_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDowb_Abs_2016_branch does not exist!\n");}
		fatJet2Pt_JESDowb_Abs_2016_isLoaded = true;
}
return fatJet2Pt_JESDowb_Abs_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESDowb_Abs_2016(){
	if(not fatJet2PtOverMHH_JESDowb_Abs_2016_isLoaded){
		if(fatJet2PtOverMHH_JESDowb_Abs_2016_branch != 0) fatJet2PtOverMHH_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDowb_Abs_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDowb_Abs_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESDowb_Abs_2016_;
}
const float &hhtree::fatJet2Pt_JESUp_BBEC1(){
	if(not fatJet2Pt_JESUp_BBEC1_isLoaded){
		if(fatJet2Pt_JESUp_BBEC1_branch != 0) fatJet2Pt_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_BBEC1_branch does not exist!\n");}
		fatJet2Pt_JESUp_BBEC1_isLoaded = true;
}
return fatJet2Pt_JESUp_BBEC1_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_BBEC1(){
	if(not fatJet2PtOverMHH_JESUp_BBEC1_isLoaded){
		if(fatJet2PtOverMHH_JESUp_BBEC1_branch != 0) fatJet2PtOverMHH_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_BBEC1_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_BBEC1_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_BBEC1_;
}
const float &hhtree::fatJet2Pt_JESDown_BBEC1(){
	if(not fatJet2Pt_JESDown_BBEC1_isLoaded){
		if(fatJet2Pt_JESDown_BBEC1_branch != 0) fatJet2Pt_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_BBEC1_branch does not exist!\n");}
		fatJet2Pt_JESDown_BBEC1_isLoaded = true;
}
return fatJet2Pt_JESDown_BBEC1_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_BBEC1(){
	if(not fatJet2PtOverMHH_JESDown_BBEC1_isLoaded){
		if(fatJet2PtOverMHH_JESDown_BBEC1_branch != 0) fatJet2PtOverMHH_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_BBEC1_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_BBEC1_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_BBEC1_;
}
const float &hhtree::fatJet2Pt_JESUp_BBEC1_2016(){
	if(not fatJet2Pt_JESUp_BBEC1_2016_isLoaded){
		if(fatJet2Pt_JESUp_BBEC1_2016_branch != 0) fatJet2Pt_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_BBEC1_2016_branch does not exist!\n");}
		fatJet2Pt_JESUp_BBEC1_2016_isLoaded = true;
}
return fatJet2Pt_JESUp_BBEC1_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_BBEC1_2016(){
	if(not fatJet2PtOverMHH_JESUp_BBEC1_2016_isLoaded){
		if(fatJet2PtOverMHH_JESUp_BBEC1_2016_branch != 0) fatJet2PtOverMHH_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_BBEC1_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_BBEC1_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_BBEC1_2016_;
}
const float &hhtree::fatJet2Pt_JESDown_BBEC1_2016(){
	if(not fatJet2Pt_JESDown_BBEC1_2016_isLoaded){
		if(fatJet2Pt_JESDown_BBEC1_2016_branch != 0) fatJet2Pt_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_BBEC1_2016_branch does not exist!\n");}
		fatJet2Pt_JESDown_BBEC1_2016_isLoaded = true;
}
return fatJet2Pt_JESDown_BBEC1_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_BBEC1_2016(){
	if(not fatJet2PtOverMHH_JESDown_BBEC1_2016_isLoaded){
		if(fatJet2PtOverMHH_JESDown_BBEC1_2016_branch != 0) fatJet2PtOverMHH_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_BBEC1_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_BBEC1_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_BBEC1_2016_;
}
const float &hhtree::fatJet2Pt_JESUp_EC2(){
	if(not fatJet2Pt_JESUp_EC2_isLoaded){
		if(fatJet2Pt_JESUp_EC2_branch != 0) fatJet2Pt_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_EC2_branch does not exist!\n");}
		fatJet2Pt_JESUp_EC2_isLoaded = true;
}
return fatJet2Pt_JESUp_EC2_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_EC2(){
	if(not fatJet2PtOverMHH_JESUp_EC2_isLoaded){
		if(fatJet2PtOverMHH_JESUp_EC2_branch != 0) fatJet2PtOverMHH_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_EC2_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_EC2_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_EC2_;
}
const float &hhtree::fatJet2Pt_JESDown_EC2(){
	if(not fatJet2Pt_JESDown_EC2_isLoaded){
		if(fatJet2Pt_JESDown_EC2_branch != 0) fatJet2Pt_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_EC2_branch does not exist!\n");}
		fatJet2Pt_JESDown_EC2_isLoaded = true;
}
return fatJet2Pt_JESDown_EC2_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_EC2(){
	if(not fatJet2PtOverMHH_JESDown_EC2_isLoaded){
		if(fatJet2PtOverMHH_JESDown_EC2_branch != 0) fatJet2PtOverMHH_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_EC2_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_EC2_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_EC2_;
}
const float &hhtree::fatJet2Pt_JESUp_EC2_2016(){
	if(not fatJet2Pt_JESUp_EC2_2016_isLoaded){
		if(fatJet2Pt_JESUp_EC2_2016_branch != 0) fatJet2Pt_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_EC2_2016_branch does not exist!\n");}
		fatJet2Pt_JESUp_EC2_2016_isLoaded = true;
}
return fatJet2Pt_JESUp_EC2_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_EC2_2016(){
	if(not fatJet2PtOverMHH_JESUp_EC2_2016_isLoaded){
		if(fatJet2PtOverMHH_JESUp_EC2_2016_branch != 0) fatJet2PtOverMHH_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_EC2_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_EC2_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_EC2_2016_;
}
const float &hhtree::fatJet2Pt_JESDown_EC2_2016(){
	if(not fatJet2Pt_JESDown_EC2_2016_isLoaded){
		if(fatJet2Pt_JESDown_EC2_2016_branch != 0) fatJet2Pt_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_EC2_2016_branch does not exist!\n");}
		fatJet2Pt_JESDown_EC2_2016_isLoaded = true;
}
return fatJet2Pt_JESDown_EC2_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_EC2_2016(){
	if(not fatJet2PtOverMHH_JESDown_EC2_2016_isLoaded){
		if(fatJet2PtOverMHH_JESDown_EC2_2016_branch != 0) fatJet2PtOverMHH_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_EC2_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_EC2_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_EC2_2016_;
}
const float &hhtree::fatJet2Pt_JESUp_FlavQCD(){
	if(not fatJet2Pt_JESUp_FlavQCD_isLoaded){
		if(fatJet2Pt_JESUp_FlavQCD_branch != 0) fatJet2Pt_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_FlavQCD_branch does not exist!\n");}
		fatJet2Pt_JESUp_FlavQCD_isLoaded = true;
}
return fatJet2Pt_JESUp_FlavQCD_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_FlavQCD(){
	if(not fatJet2PtOverMHH_JESUp_FlavQCD_isLoaded){
		if(fatJet2PtOverMHH_JESUp_FlavQCD_branch != 0) fatJet2PtOverMHH_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_FlavQCD_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_FlavQCD_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_FlavQCD_;
}
const float &hhtree::fatJet2Pt_JESDown_FlavQCD(){
	if(not fatJet2Pt_JESDown_FlavQCD_isLoaded){
		if(fatJet2Pt_JESDown_FlavQCD_branch != 0) fatJet2Pt_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_FlavQCD_branch does not exist!\n");}
		fatJet2Pt_JESDown_FlavQCD_isLoaded = true;
}
return fatJet2Pt_JESDown_FlavQCD_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_FlavQCD(){
	if(not fatJet2PtOverMHH_JESDown_FlavQCD_isLoaded){
		if(fatJet2PtOverMHH_JESDown_FlavQCD_branch != 0) fatJet2PtOverMHH_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_FlavQCD_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_FlavQCD_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_FlavQCD_;
}
const float &hhtree::fatJet2Pt_JESUp_HF(){
	if(not fatJet2Pt_JESUp_HF_isLoaded){
		if(fatJet2Pt_JESUp_HF_branch != 0) fatJet2Pt_JESUp_HF_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_HF_branch does not exist!\n");}
		fatJet2Pt_JESUp_HF_isLoaded = true;
}
return fatJet2Pt_JESUp_HF_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_HF(){
	if(not fatJet2PtOverMHH_JESUp_HF_isLoaded){
		if(fatJet2PtOverMHH_JESUp_HF_branch != 0) fatJet2PtOverMHH_JESUp_HF_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_HF_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_HF_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_HF_;
}
const float &hhtree::fatJet2Pt_JESDown_HF(){
	if(not fatJet2Pt_JESDown_HF_isLoaded){
		if(fatJet2Pt_JESDown_HF_branch != 0) fatJet2Pt_JESDown_HF_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_HF_branch does not exist!\n");}
		fatJet2Pt_JESDown_HF_isLoaded = true;
}
return fatJet2Pt_JESDown_HF_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_HF(){
	if(not fatJet2PtOverMHH_JESDown_HF_isLoaded){
		if(fatJet2PtOverMHH_JESDown_HF_branch != 0) fatJet2PtOverMHH_JESDown_HF_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_HF_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_HF_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_HF_;
}
const float &hhtree::fatJet2Pt_JESUp_HF_2016(){
	if(not fatJet2Pt_JESUp_HF_2016_isLoaded){
		if(fatJet2Pt_JESUp_HF_2016_branch != 0) fatJet2Pt_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_HF_2016_branch does not exist!\n");}
		fatJet2Pt_JESUp_HF_2016_isLoaded = true;
}
return fatJet2Pt_JESUp_HF_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_HF_2016(){
	if(not fatJet2PtOverMHH_JESUp_HF_2016_isLoaded){
		if(fatJet2PtOverMHH_JESUp_HF_2016_branch != 0) fatJet2PtOverMHH_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_HF_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_HF_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_HF_2016_;
}
const float &hhtree::fatJet2Pt_JESDown_HF_2016(){
	if(not fatJet2Pt_JESDown_HF_2016_isLoaded){
		if(fatJet2Pt_JESDown_HF_2016_branch != 0) fatJet2Pt_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_HF_2016_branch does not exist!\n");}
		fatJet2Pt_JESDown_HF_2016_isLoaded = true;
}
return fatJet2Pt_JESDown_HF_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_HF_2016(){
	if(not fatJet2PtOverMHH_JESDown_HF_2016_isLoaded){
		if(fatJet2PtOverMHH_JESDown_HF_2016_branch != 0) fatJet2PtOverMHH_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_HF_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_HF_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_HF_2016_;
}
const float &hhtree::fatJet2Pt_JESUp_RelBal(){
	if(not fatJet2Pt_JESUp_RelBal_isLoaded){
		if(fatJet2Pt_JESUp_RelBal_branch != 0) fatJet2Pt_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_RelBal_branch does not exist!\n");}
		fatJet2Pt_JESUp_RelBal_isLoaded = true;
}
return fatJet2Pt_JESUp_RelBal_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_RelBal(){
	if(not fatJet2PtOverMHH_JESUp_RelBal_isLoaded){
		if(fatJet2PtOverMHH_JESUp_RelBal_branch != 0) fatJet2PtOverMHH_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_RelBal_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_RelBal_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_RelBal_;
}
const float &hhtree::fatJet2Pt_JESDown_RelBal(){
	if(not fatJet2Pt_JESDown_RelBal_isLoaded){
		if(fatJet2Pt_JESDown_RelBal_branch != 0) fatJet2Pt_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_RelBal_branch does not exist!\n");}
		fatJet2Pt_JESDown_RelBal_isLoaded = true;
}
return fatJet2Pt_JESDown_RelBal_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_RelBal(){
	if(not fatJet2PtOverMHH_JESDown_RelBal_isLoaded){
		if(fatJet2PtOverMHH_JESDown_RelBal_branch != 0) fatJet2PtOverMHH_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_RelBal_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_RelBal_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_RelBal_;
}
const float &hhtree::fatJet2Pt_JESUp_RelSample_2016(){
	if(not fatJet2Pt_JESUp_RelSample_2016_isLoaded){
		if(fatJet2Pt_JESUp_RelSample_2016_branch != 0) fatJet2Pt_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_RelSample_2016_branch does not exist!\n");}
		fatJet2Pt_JESUp_RelSample_2016_isLoaded = true;
}
return fatJet2Pt_JESUp_RelSample_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_RelSample_2016(){
	if(not fatJet2PtOverMHH_JESUp_RelSample_2016_isLoaded){
		if(fatJet2PtOverMHH_JESUp_RelSample_2016_branch != 0) fatJet2PtOverMHH_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_RelSample_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_RelSample_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_RelSample_2016_;
}
const float &hhtree::fatJet2Pt_JESDown_RelSample_2016(){
	if(not fatJet2Pt_JESDown_RelSample_2016_isLoaded){
		if(fatJet2Pt_JESDown_RelSample_2016_branch != 0) fatJet2Pt_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_RelSample_2016_branch does not exist!\n");}
		fatJet2Pt_JESDown_RelSample_2016_isLoaded = true;
}
return fatJet2Pt_JESDown_RelSample_2016_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_RelSample_2016(){
	if(not fatJet2PtOverMHH_JESDown_RelSample_2016_isLoaded){
		if(fatJet2PtOverMHH_JESDown_RelSample_2016_branch != 0) fatJet2PtOverMHH_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_RelSample_2016_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_RelSample_2016_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_RelSample_2016_;
}


const float &hhtree::hh_pt() 
{
 if(not hh_pt_isLoaded)
 {
   if(hh_pt_branch != 0) hh_pt_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_branch does not exist!\n");
//exit(1);
   }
   hh_pt_isLoaded = true;
 }
 return hh_pt_;
}

const float &hhtree::hh_eta() 
{
 if(not hh_eta_isLoaded)
 {
   if(hh_eta_branch != 0) hh_eta_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_branch does not exist!\n");
//exit(1);
   }
   hh_eta_isLoaded = true;
 }
 return hh_eta_;
}

const float &hhtree::hh_phi() 
{
 if(not hh_phi_isLoaded)
 {
   if(hh_phi_branch != 0) hh_phi_branch->GetEntry(index);
   else
   {
     printf("branch hh_phi_branch does not exist!\n");
//exit(1);
   }
   hh_phi_isLoaded = true;
 }
 return hh_phi_;
}

const float &hhtree::hh_mass() 
{
 if(not hh_mass_isLoaded)
 {
   if(hh_mass_branch != 0) hh_mass_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_branch does not exist!\n");
//exit(1);
   }
   hh_mass_isLoaded = true;
 }
 return hh_mass_;
}

const float &hhtree::hh_pt_JMS_Down() 
{
 if(not hh_pt_JMS_Down_isLoaded)
 {
   if(hh_pt_JMS_Down_branch != 0) hh_pt_JMS_Down_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JMS_Down_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JMS_Down_isLoaded = true;
 }
 return hh_pt_JMS_Down_;
}

const float &hhtree::hh_pt_JMS_Up() 
{
 if(not hh_pt_JMS_Up_isLoaded)
 {
   if(hh_pt_JMS_Up_branch != 0) hh_pt_JMS_Up_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JMS_Up_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JMS_Up_isLoaded = true;
 }
 return hh_pt_JMS_Up_;
}

const float &hhtree::hh_eta_JMS_Down() 
{
 if(not hh_eta_JMS_Down_isLoaded)
 {
   if(hh_eta_JMS_Down_branch != 0) hh_eta_JMS_Down_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JMS_Down_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JMS_Down_isLoaded = true;
 }
 return hh_eta_JMS_Down_;
}

const float &hhtree::hh_eta_JMS_Up() 
{
 if(not hh_eta_JMS_Up_isLoaded)
 {
   if(hh_eta_JMS_Up_branch != 0) hh_eta_JMS_Up_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JMS_Up_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JMS_Up_isLoaded = true;
 }
 return hh_eta_JMS_Up_;
}

const float &hhtree::hh_mass_JMS_Down() 
{
 if(not hh_mass_JMS_Down_isLoaded)
 {
   if(hh_mass_JMS_Down_branch != 0) hh_mass_JMS_Down_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JMS_Down_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JMS_Down_isLoaded = true;
 }
 return hh_mass_JMS_Down_;
}

const float &hhtree::hh_mass_JMS_Up() 
{
 if(not hh_mass_JMS_Up_isLoaded)
 {
   if(hh_mass_JMS_Up_branch != 0) hh_mass_JMS_Up_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JMS_Up_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JMS_Up_isLoaded = true;
 }
 return hh_mass_JMS_Up_;
}

const float &hhtree::hh_pt_JMR_Down() 
{
 if(not hh_pt_JMR_Down_isLoaded)
 {
   if(hh_pt_JMR_Down_branch != 0) hh_pt_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JMR_Down_isLoaded = true;
 }
 return hh_pt_JMR_Down_;
}

const float &hhtree::hh_pt_JMR_Up() 
{
 if(not hh_pt_JMR_Up_isLoaded)
 {
   if(hh_pt_JMR_Up_branch != 0) hh_pt_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JMR_Up_isLoaded = true;
 }
 return hh_pt_JMR_Up_;
}

const float &hhtree::hh_eta_JMR_Down() 
{
 if(not hh_eta_JMR_Down_isLoaded)
 {
   if(hh_eta_JMR_Down_branch != 0) hh_eta_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JMR_Down_isLoaded = true;
 }
 return hh_eta_JMR_Down_;
}

const float &hhtree::hh_eta_JMR_Up() 
{
 if(not hh_eta_JMR_Up_isLoaded)
 {
   if(hh_eta_JMR_Up_branch != 0) hh_eta_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JMR_Up_isLoaded = true;
 }
 return hh_eta_JMR_Up_;
}

const float &hhtree::hh_mass_JMR_Down() 
{
 if(not hh_mass_JMR_Down_isLoaded)
 {
   if(hh_mass_JMR_Down_branch != 0) hh_mass_JMR_Down_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JMR_Down_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JMR_Down_isLoaded = true;
 }
 return hh_mass_JMR_Down_;
}

const float &hhtree::hh_mass_JMR_Up() 
{
 if(not hh_mass_JMR_Up_isLoaded)
 {
   if(hh_mass_JMR_Up_branch != 0) hh_mass_JMR_Up_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JMR_Up_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JMR_Up_isLoaded = true;
 }
 return hh_mass_JMR_Up_;
}

const float &hhtree::hh_pt_JERUp() 
{
 if(not hh_pt_JERUp_isLoaded)
 {
   if(hh_pt_JERUp_branch != 0) hh_pt_JERUp_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JERUp_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JERUp_isLoaded = true;
 }
 return hh_pt_JERUp_;
}

const float &hhtree::hh_eta_JERUp() 
{
 if(not hh_eta_JERUp_isLoaded)
 {
   if(hh_eta_JERUp_branch != 0) hh_eta_JERUp_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JERUp_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JERUp_isLoaded = true;
 }
 return hh_eta_JERUp_;
}

const float &hhtree::hh_mass_JERUp() 
{
 if(not hh_mass_JERUp_isLoaded)
 {
   if(hh_mass_JERUp_branch != 0) hh_mass_JERUp_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JERUp_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JERUp_isLoaded = true;
 }
 return hh_mass_JERUp_;
}

const float &hhtree::hh_pt_JERDown() 
{
 if(not hh_pt_JERDown_isLoaded)
 {
   if(hh_pt_JERDown_branch != 0) hh_pt_JERDown_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JERDown_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JERDown_isLoaded = true;
 }
 return hh_pt_JERDown_;
}

const float &hhtree::hh_eta_JERDown() 
{
 if(not hh_eta_JERDown_isLoaded)
 {
   if(hh_eta_JERDown_branch != 0) hh_eta_JERDown_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JERDown_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JERDown_isLoaded = true;
 }
 return hh_eta_JERDown_;
}

const float &hhtree::hh_mass_JERDown() 
{
 if(not hh_mass_JERDown_isLoaded)
 {
   if(hh_mass_JERDown_branch != 0) hh_mass_JERDown_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JERDown_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JERDown_isLoaded = true;
 }
 return hh_mass_JERDown_;
}

const float &hhtree::hh_pt_JESUp() 
{
 if(not hh_pt_JESUp_isLoaded)
 {
   if(hh_pt_JESUp_branch != 0) hh_pt_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JESUp_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JESUp_isLoaded = true;
 }
 return hh_pt_JESUp_;
}

const float &hhtree::hh_eta_JESUp() 
{
 if(not hh_eta_JESUp_isLoaded)
 {
   if(hh_eta_JESUp_branch != 0) hh_eta_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JESUp_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JESUp_isLoaded = true;
 }
 return hh_eta_JESUp_;
}

const float &hhtree::hh_mass_JESUp() 
{
 if(not hh_mass_JESUp_isLoaded)
 {
   if(hh_mass_JESUp_branch != 0) hh_mass_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JESUp_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JESUp_isLoaded = true;
 }
 return hh_mass_JESUp_;
}

const float &hhtree::hh_pt_JESDown() 
{
 if(not hh_pt_JESDown_isLoaded)
 {
   if(hh_pt_JESDown_branch != 0) hh_pt_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch hh_pt_JESDown_branch does not exist!\n");
//exit(1);
   }
   hh_pt_JESDown_isLoaded = true;
 }
 return hh_pt_JESDown_;
}

const float &hhtree::hh_eta_JESDown() 
{
 if(not hh_eta_JESDown_isLoaded)
 {
   if(hh_eta_JESDown_branch != 0) hh_eta_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch hh_eta_JESDown_branch does not exist!\n");
//exit(1);
   }
   hh_eta_JESDown_isLoaded = true;
 }
 return hh_eta_JESDown_;
}

const float &hhtree::hh_mass_JESDown() 
{
 if(not hh_mass_JESDown_isLoaded)
 {
   if(hh_mass_JESDown_branch != 0) hh_mass_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch hh_mass_JESDown_branch does not exist!\n");
//exit(1);
   }
   hh_mass_JESDown_isLoaded = true;
 }
 return hh_mass_JESDown_;
}
const float &hhtree::hh_pt_JESUp_Abs(){
	if(not hh_pt_JESUp_Abs_isLoaded){
		if(hh_pt_JESUp_Abs_branch != 0) hh_pt_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_Abs_branch does not exist!\n");}
		hh_pt_JESUp_Abs_isLoaded = true;
}
return hh_pt_JESUp_Abs_;
}
const float &hhtree::hh_eta_JESUp_Abs(){
	if(not hh_eta_JESUp_Abs_isLoaded){
		if(hh_eta_JESUp_Abs_branch != 0) hh_eta_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_Abs_branch does not exist!\n");}
		hh_eta_JESUp_Abs_isLoaded = true;
}
return hh_eta_JESUp_Abs_;
}
const float &hhtree::hh_mass_JESUp_Abs(){
	if(not hh_mass_JESUp_Abs_isLoaded){
		if(hh_mass_JESUp_Abs_branch != 0) hh_mass_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_Abs_branch does not exist!\n");}
		hh_mass_JESUp_Abs_isLoaded = true;
}
return hh_mass_JESUp_Abs_;
}
const float &hhtree::hh_pt_JESDown_Abs(){
	if(not hh_pt_JESDown_Abs_isLoaded){
		if(hh_pt_JESDown_Abs_branch != 0) hh_pt_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_Abs_branch does not exist!\n");}
		hh_pt_JESDown_Abs_isLoaded = true;
}
return hh_pt_JESDown_Abs_;
}
const float &hhtree::hh_eta_JESDown_Abs(){
	if(not hh_eta_JESDown_Abs_isLoaded){
		if(hh_eta_JESDown_Abs_branch != 0) hh_eta_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_Abs_branch does not exist!\n");}
		hh_eta_JESDown_Abs_isLoaded = true;
}
return hh_eta_JESDown_Abs_;
}
const float &hhtree::hh_mass_JESDown_Abs(){
	if(not hh_mass_JESDown_Abs_isLoaded){
		if(hh_mass_JESDown_Abs_branch != 0) hh_mass_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_Abs_branch does not exist!\n");}
		hh_mass_JESDown_Abs_isLoaded = true;
}
return hh_mass_JESDown_Abs_;
}
const float &hhtree::hh_pt_JESUp_Abs_2016(){
	if(not hh_pt_JESUp_Abs_2016_isLoaded){
		if(hh_pt_JESUp_Abs_2016_branch != 0) hh_pt_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_Abs_2016_branch does not exist!\n");}
		hh_pt_JESUp_Abs_2016_isLoaded = true;
}
return hh_pt_JESUp_Abs_2016_;
}
const float &hhtree::hh_eta_JESUp_Abs_2016(){
	if(not hh_eta_JESUp_Abs_2016_isLoaded){
		if(hh_eta_JESUp_Abs_2016_branch != 0) hh_eta_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_Abs_2016_branch does not exist!\n");}
		hh_eta_JESUp_Abs_2016_isLoaded = true;
}
return hh_eta_JESUp_Abs_2016_;
}
const float &hhtree::hh_mass_JESUp_Abs_2016(){
	if(not hh_mass_JESUp_Abs_2016_isLoaded){
		if(hh_mass_JESUp_Abs_2016_branch != 0) hh_mass_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_Abs_2016_branch does not exist!\n");}
		hh_mass_JESUp_Abs_2016_isLoaded = true;
}
return hh_mass_JESUp_Abs_2016_;
}
const float &hhtree::hh_pt_JESDowb_Abs_2016(){
	if(not hh_pt_JESDowb_Abs_2016_isLoaded){
		if(hh_pt_JESDowb_Abs_2016_branch != 0) hh_pt_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDowb_Abs_2016_branch does not exist!\n");}
		hh_pt_JESDowb_Abs_2016_isLoaded = true;
}
return hh_pt_JESDowb_Abs_2016_;
}
const float &hhtree::hh_eta_JESDowb_Abs_2016(){
	if(not hh_eta_JESDowb_Abs_2016_isLoaded){
		if(hh_eta_JESDowb_Abs_2016_branch != 0) hh_eta_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDowb_Abs_2016_branch does not exist!\n");}
		hh_eta_JESDowb_Abs_2016_isLoaded = true;
}
return hh_eta_JESDowb_Abs_2016_;
}
const float &hhtree::hh_mass_JESDowb_Abs_2016(){
	if(not hh_mass_JESDowb_Abs_2016_isLoaded){
		if(hh_mass_JESDowb_Abs_2016_branch != 0) hh_mass_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDowb_Abs_2016_branch does not exist!\n");}
		hh_mass_JESDowb_Abs_2016_isLoaded = true;
}
return hh_mass_JESDowb_Abs_2016_;
}
const float &hhtree::hh_pt_JESUp_BBEC1(){
	if(not hh_pt_JESUp_BBEC1_isLoaded){
		if(hh_pt_JESUp_BBEC1_branch != 0) hh_pt_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_BBEC1_branch does not exist!\n");}
		hh_pt_JESUp_BBEC1_isLoaded = true;
}
return hh_pt_JESUp_BBEC1_;
}
const float &hhtree::hh_eta_JESUp_BBEC1(){
	if(not hh_eta_JESUp_BBEC1_isLoaded){
		if(hh_eta_JESUp_BBEC1_branch != 0) hh_eta_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_BBEC1_branch does not exist!\n");}
		hh_eta_JESUp_BBEC1_isLoaded = true;
}
return hh_eta_JESUp_BBEC1_;
}
const float &hhtree::hh_mass_JESUp_BBEC1(){
	if(not hh_mass_JESUp_BBEC1_isLoaded){
		if(hh_mass_JESUp_BBEC1_branch != 0) hh_mass_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_BBEC1_branch does not exist!\n");}
		hh_mass_JESUp_BBEC1_isLoaded = true;
}
return hh_mass_JESUp_BBEC1_;
}
const float &hhtree::hh_pt_JESDown_BBEC1(){
	if(not hh_pt_JESDown_BBEC1_isLoaded){
		if(hh_pt_JESDown_BBEC1_branch != 0) hh_pt_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_BBEC1_branch does not exist!\n");}
		hh_pt_JESDown_BBEC1_isLoaded = true;
}
return hh_pt_JESDown_BBEC1_;
}
const float &hhtree::hh_eta_JESDown_BBEC1(){
	if(not hh_eta_JESDown_BBEC1_isLoaded){
		if(hh_eta_JESDown_BBEC1_branch != 0) hh_eta_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_BBEC1_branch does not exist!\n");}
		hh_eta_JESDown_BBEC1_isLoaded = true;
}
return hh_eta_JESDown_BBEC1_;
}
const float &hhtree::hh_mass_JESDown_BBEC1(){
	if(not hh_mass_JESDown_BBEC1_isLoaded){
		if(hh_mass_JESDown_BBEC1_branch != 0) hh_mass_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_BBEC1_branch does not exist!\n");}
		hh_mass_JESDown_BBEC1_isLoaded = true;
}
return hh_mass_JESDown_BBEC1_;
}
const float &hhtree::hh_pt_JESUp_BBEC1_2016(){
	if(not hh_pt_JESUp_BBEC1_2016_isLoaded){
		if(hh_pt_JESUp_BBEC1_2016_branch != 0) hh_pt_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_BBEC1_2016_branch does not exist!\n");}
		hh_pt_JESUp_BBEC1_2016_isLoaded = true;
}
return hh_pt_JESUp_BBEC1_2016_;
}
const float &hhtree::hh_eta_JESUp_BBEC1_2016(){
	if(not hh_eta_JESUp_BBEC1_2016_isLoaded){
		if(hh_eta_JESUp_BBEC1_2016_branch != 0) hh_eta_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_BBEC1_2016_branch does not exist!\n");}
		hh_eta_JESUp_BBEC1_2016_isLoaded = true;
}
return hh_eta_JESUp_BBEC1_2016_;
}
const float &hhtree::hh_mass_JESUp_BBEC1_2016(){
	if(not hh_mass_JESUp_BBEC1_2016_isLoaded){
		if(hh_mass_JESUp_BBEC1_2016_branch != 0) hh_mass_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_BBEC1_2016_branch does not exist!\n");}
		hh_mass_JESUp_BBEC1_2016_isLoaded = true;
}
return hh_mass_JESUp_BBEC1_2016_;
}
const float &hhtree::hh_pt_JESDown_BBEC1_2016(){
	if(not hh_pt_JESDown_BBEC1_2016_isLoaded){
		if(hh_pt_JESDown_BBEC1_2016_branch != 0) hh_pt_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_BBEC1_2016_branch does not exist!\n");}
		hh_pt_JESDown_BBEC1_2016_isLoaded = true;
}
return hh_pt_JESDown_BBEC1_2016_;
}
const float &hhtree::hh_eta_JESDown_BBEC1_2016(){
	if(not hh_eta_JESDown_BBEC1_2016_isLoaded){
		if(hh_eta_JESDown_BBEC1_2016_branch != 0) hh_eta_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_BBEC1_2016_branch does not exist!\n");}
		hh_eta_JESDown_BBEC1_2016_isLoaded = true;
}
return hh_eta_JESDown_BBEC1_2016_;
}
const float &hhtree::hh_mass_JESDown_BBEC1_2016(){
	if(not hh_mass_JESDown_BBEC1_2016_isLoaded){
		if(hh_mass_JESDown_BBEC1_2016_branch != 0) hh_mass_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_BBEC1_2016_branch does not exist!\n");}
		hh_mass_JESDown_BBEC1_2016_isLoaded = true;
}
return hh_mass_JESDown_BBEC1_2016_;
}
const float &hhtree::hh_pt_JESUp_EC2(){
	if(not hh_pt_JESUp_EC2_isLoaded){
		if(hh_pt_JESUp_EC2_branch != 0) hh_pt_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_EC2_branch does not exist!\n");}
		hh_pt_JESUp_EC2_isLoaded = true;
}
return hh_pt_JESUp_EC2_;
}
const float &hhtree::hh_eta_JESUp_EC2(){
	if(not hh_eta_JESUp_EC2_isLoaded){
		if(hh_eta_JESUp_EC2_branch != 0) hh_eta_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_EC2_branch does not exist!\n");}
		hh_eta_JESUp_EC2_isLoaded = true;
}
return hh_eta_JESUp_EC2_;
}
const float &hhtree::hh_mass_JESUp_EC2(){
	if(not hh_mass_JESUp_EC2_isLoaded){
		if(hh_mass_JESUp_EC2_branch != 0) hh_mass_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_EC2_branch does not exist!\n");}
		hh_mass_JESUp_EC2_isLoaded = true;
}
return hh_mass_JESUp_EC2_;
}
const float &hhtree::hh_pt_JESDown_EC2(){
	if(not hh_pt_JESDown_EC2_isLoaded){
		if(hh_pt_JESDown_EC2_branch != 0) hh_pt_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_EC2_branch does not exist!\n");}
		hh_pt_JESDown_EC2_isLoaded = true;
}
return hh_pt_JESDown_EC2_;
}
const float &hhtree::hh_eta_JESDown_EC2(){
	if(not hh_eta_JESDown_EC2_isLoaded){
		if(hh_eta_JESDown_EC2_branch != 0) hh_eta_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_EC2_branch does not exist!\n");}
		hh_eta_JESDown_EC2_isLoaded = true;
}
return hh_eta_JESDown_EC2_;
}
const float &hhtree::hh_mass_JESDown_EC2(){
	if(not hh_mass_JESDown_EC2_isLoaded){
		if(hh_mass_JESDown_EC2_branch != 0) hh_mass_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_EC2_branch does not exist!\n");}
		hh_mass_JESDown_EC2_isLoaded = true;
}
return hh_mass_JESDown_EC2_;
}
const float &hhtree::hh_pt_JESUp_EC2_2016(){
	if(not hh_pt_JESUp_EC2_2016_isLoaded){
		if(hh_pt_JESUp_EC2_2016_branch != 0) hh_pt_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_EC2_2016_branch does not exist!\n");}
		hh_pt_JESUp_EC2_2016_isLoaded = true;
}
return hh_pt_JESUp_EC2_2016_;
}
const float &hhtree::hh_eta_JESUp_EC2_2016(){
	if(not hh_eta_JESUp_EC2_2016_isLoaded){
		if(hh_eta_JESUp_EC2_2016_branch != 0) hh_eta_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_EC2_2016_branch does not exist!\n");}
		hh_eta_JESUp_EC2_2016_isLoaded = true;
}
return hh_eta_JESUp_EC2_2016_;
}
const float &hhtree::hh_mass_JESUp_EC2_2016(){
	if(not hh_mass_JESUp_EC2_2016_isLoaded){
		if(hh_mass_JESUp_EC2_2016_branch != 0) hh_mass_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_EC2_2016_branch does not exist!\n");}
		hh_mass_JESUp_EC2_2016_isLoaded = true;
}
return hh_mass_JESUp_EC2_2016_;
}
const float &hhtree::hh_pt_JESDown_EC2_2016(){
	if(not hh_pt_JESDown_EC2_2016_isLoaded){
		if(hh_pt_JESDown_EC2_2016_branch != 0) hh_pt_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_EC2_2016_branch does not exist!\n");}
		hh_pt_JESDown_EC2_2016_isLoaded = true;
}
return hh_pt_JESDown_EC2_2016_;
}
const float &hhtree::hh_eta_JESDown_EC2_2016(){
	if(not hh_eta_JESDown_EC2_2016_isLoaded){
		if(hh_eta_JESDown_EC2_2016_branch != 0) hh_eta_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_EC2_2016_branch does not exist!\n");}
		hh_eta_JESDown_EC2_2016_isLoaded = true;
}
return hh_eta_JESDown_EC2_2016_;
}
const float &hhtree::hh_mass_JESDown_EC2_2016(){
	if(not hh_mass_JESDown_EC2_2016_isLoaded){
		if(hh_mass_JESDown_EC2_2016_branch != 0) hh_mass_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_EC2_2016_branch does not exist!\n");}
		hh_mass_JESDown_EC2_2016_isLoaded = true;
}
return hh_mass_JESDown_EC2_2016_;
}
const float &hhtree::hh_pt_JESUp_FlavQCD(){
	if(not hh_pt_JESUp_FlavQCD_isLoaded){
		if(hh_pt_JESUp_FlavQCD_branch != 0) hh_pt_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_FlavQCD_branch does not exist!\n");}
		hh_pt_JESUp_FlavQCD_isLoaded = true;
}
return hh_pt_JESUp_FlavQCD_;
}
const float &hhtree::hh_eta_JESUp_FlavQCD(){
	if(not hh_eta_JESUp_FlavQCD_isLoaded){
		if(hh_eta_JESUp_FlavQCD_branch != 0) hh_eta_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_FlavQCD_branch does not exist!\n");}
		hh_eta_JESUp_FlavQCD_isLoaded = true;
}
return hh_eta_JESUp_FlavQCD_;
}
const float &hhtree::hh_mass_JESUp_FlavQCD(){
	if(not hh_mass_JESUp_FlavQCD_isLoaded){
		if(hh_mass_JESUp_FlavQCD_branch != 0) hh_mass_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_FlavQCD_branch does not exist!\n");}
		hh_mass_JESUp_FlavQCD_isLoaded = true;
}
return hh_mass_JESUp_FlavQCD_;
}
const float &hhtree::hh_pt_JESDown_FlavQCD(){
	if(not hh_pt_JESDown_FlavQCD_isLoaded){
		if(hh_pt_JESDown_FlavQCD_branch != 0) hh_pt_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_FlavQCD_branch does not exist!\n");}
		hh_pt_JESDown_FlavQCD_isLoaded = true;
}
return hh_pt_JESDown_FlavQCD_;
}
const float &hhtree::hh_eta_JESDown_FlavQCD(){
	if(not hh_eta_JESDown_FlavQCD_isLoaded){
		if(hh_eta_JESDown_FlavQCD_branch != 0) hh_eta_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_FlavQCD_branch does not exist!\n");}
		hh_eta_JESDown_FlavQCD_isLoaded = true;
}
return hh_eta_JESDown_FlavQCD_;
}
const float &hhtree::hh_mass_JESDown_FlavQCD(){
	if(not hh_mass_JESDown_FlavQCD_isLoaded){
		if(hh_mass_JESDown_FlavQCD_branch != 0) hh_mass_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_FlavQCD_branch does not exist!\n");}
		hh_mass_JESDown_FlavQCD_isLoaded = true;
}
return hh_mass_JESDown_FlavQCD_;
}
const float &hhtree::hh_pt_JESUp_HF(){
	if(not hh_pt_JESUp_HF_isLoaded){
		if(hh_pt_JESUp_HF_branch != 0) hh_pt_JESUp_HF_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_HF_branch does not exist!\n");}
		hh_pt_JESUp_HF_isLoaded = true;
}
return hh_pt_JESUp_HF_;
}
const float &hhtree::hh_eta_JESUp_HF(){
	if(not hh_eta_JESUp_HF_isLoaded){
		if(hh_eta_JESUp_HF_branch != 0) hh_eta_JESUp_HF_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_HF_branch does not exist!\n");}
		hh_eta_JESUp_HF_isLoaded = true;
}
return hh_eta_JESUp_HF_;
}
const float &hhtree::hh_mass_JESUp_HF(){
	if(not hh_mass_JESUp_HF_isLoaded){
		if(hh_mass_JESUp_HF_branch != 0) hh_mass_JESUp_HF_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_HF_branch does not exist!\n");}
		hh_mass_JESUp_HF_isLoaded = true;
}
return hh_mass_JESUp_HF_;
}
const float &hhtree::hh_pt_JESDown_HF(){
	if(not hh_pt_JESDown_HF_isLoaded){
		if(hh_pt_JESDown_HF_branch != 0) hh_pt_JESDown_HF_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_HF_branch does not exist!\n");}
		hh_pt_JESDown_HF_isLoaded = true;
}
return hh_pt_JESDown_HF_;
}
const float &hhtree::hh_eta_JESDown_HF(){
	if(not hh_eta_JESDown_HF_isLoaded){
		if(hh_eta_JESDown_HF_branch != 0) hh_eta_JESDown_HF_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_HF_branch does not exist!\n");}
		hh_eta_JESDown_HF_isLoaded = true;
}
return hh_eta_JESDown_HF_;
}
const float &hhtree::hh_mass_JESDown_HF(){
	if(not hh_mass_JESDown_HF_isLoaded){
		if(hh_mass_JESDown_HF_branch != 0) hh_mass_JESDown_HF_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_HF_branch does not exist!\n");}
		hh_mass_JESDown_HF_isLoaded = true;
}
return hh_mass_JESDown_HF_;
}
const float &hhtree::hh_pt_JESUp_HF_2016(){
	if(not hh_pt_JESUp_HF_2016_isLoaded){
		if(hh_pt_JESUp_HF_2016_branch != 0) hh_pt_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_HF_2016_branch does not exist!\n");}
		hh_pt_JESUp_HF_2016_isLoaded = true;
}
return hh_pt_JESUp_HF_2016_;
}
const float &hhtree::hh_eta_JESUp_HF_2016(){
	if(not hh_eta_JESUp_HF_2016_isLoaded){
		if(hh_eta_JESUp_HF_2016_branch != 0) hh_eta_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_HF_2016_branch does not exist!\n");}
		hh_eta_JESUp_HF_2016_isLoaded = true;
}
return hh_eta_JESUp_HF_2016_;
}
const float &hhtree::hh_mass_JESUp_HF_2016(){
	if(not hh_mass_JESUp_HF_2016_isLoaded){
		if(hh_mass_JESUp_HF_2016_branch != 0) hh_mass_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_HF_2016_branch does not exist!\n");}
		hh_mass_JESUp_HF_2016_isLoaded = true;
}
return hh_mass_JESUp_HF_2016_;
}
const float &hhtree::hh_pt_JESDown_HF_2016(){
	if(not hh_pt_JESDown_HF_2016_isLoaded){
		if(hh_pt_JESDown_HF_2016_branch != 0) hh_pt_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_HF_2016_branch does not exist!\n");}
		hh_pt_JESDown_HF_2016_isLoaded = true;
}
return hh_pt_JESDown_HF_2016_;
}
const float &hhtree::hh_eta_JESDown_HF_2016(){
	if(not hh_eta_JESDown_HF_2016_isLoaded){
		if(hh_eta_JESDown_HF_2016_branch != 0) hh_eta_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_HF_2016_branch does not exist!\n");}
		hh_eta_JESDown_HF_2016_isLoaded = true;
}
return hh_eta_JESDown_HF_2016_;
}
const float &hhtree::hh_mass_JESDown_HF_2016(){
	if(not hh_mass_JESDown_HF_2016_isLoaded){
		if(hh_mass_JESDown_HF_2016_branch != 0) hh_mass_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_HF_2016_branch does not exist!\n");}
		hh_mass_JESDown_HF_2016_isLoaded = true;
}
return hh_mass_JESDown_HF_2016_;
}
const float &hhtree::hh_pt_JESUp_RelBal(){
	if(not hh_pt_JESUp_RelBal_isLoaded){
		if(hh_pt_JESUp_RelBal_branch != 0) hh_pt_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_RelBal_branch does not exist!\n");}
		hh_pt_JESUp_RelBal_isLoaded = true;
}
return hh_pt_JESUp_RelBal_;
}
const float &hhtree::hh_eta_JESUp_RelBal(){
	if(not hh_eta_JESUp_RelBal_isLoaded){
		if(hh_eta_JESUp_RelBal_branch != 0) hh_eta_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_RelBal_branch does not exist!\n");}
		hh_eta_JESUp_RelBal_isLoaded = true;
}
return hh_eta_JESUp_RelBal_;
}
const float &hhtree::hh_mass_JESUp_RelBal(){
	if(not hh_mass_JESUp_RelBal_isLoaded){
		if(hh_mass_JESUp_RelBal_branch != 0) hh_mass_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_RelBal_branch does not exist!\n");}
		hh_mass_JESUp_RelBal_isLoaded = true;
}
return hh_mass_JESUp_RelBal_;
}
const float &hhtree::hh_pt_JESDown_RelBal(){
	if(not hh_pt_JESDown_RelBal_isLoaded){
		if(hh_pt_JESDown_RelBal_branch != 0) hh_pt_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_RelBal_branch does not exist!\n");}
		hh_pt_JESDown_RelBal_isLoaded = true;
}
return hh_pt_JESDown_RelBal_;
}
const float &hhtree::hh_eta_JESDown_RelBal(){
	if(not hh_eta_JESDown_RelBal_isLoaded){
		if(hh_eta_JESDown_RelBal_branch != 0) hh_eta_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_RelBal_branch does not exist!\n");}
		hh_eta_JESDown_RelBal_isLoaded = true;
}
return hh_eta_JESDown_RelBal_;
}
const float &hhtree::hh_mass_JESDown_RelBal(){
	if(not hh_mass_JESDown_RelBal_isLoaded){
		if(hh_mass_JESDown_RelBal_branch != 0) hh_mass_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_RelBal_branch does not exist!\n");}
		hh_mass_JESDown_RelBal_isLoaded = true;
}
return hh_mass_JESDown_RelBal_;
}
const float &hhtree::hh_pt_JESUp_RelSample_2016(){
	if(not hh_pt_JESUp_RelSample_2016_isLoaded){
		if(hh_pt_JESUp_RelSample_2016_branch != 0) hh_pt_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_RelSample_2016_branch does not exist!\n");}
		hh_pt_JESUp_RelSample_2016_isLoaded = true;
}
return hh_pt_JESUp_RelSample_2016_;
}
const float &hhtree::hh_eta_JESUp_RelSample_2016(){
	if(not hh_eta_JESUp_RelSample_2016_isLoaded){
		if(hh_eta_JESUp_RelSample_2016_branch != 0) hh_eta_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_RelSample_2016_branch does not exist!\n");}
		hh_eta_JESUp_RelSample_2016_isLoaded = true;
}
return hh_eta_JESUp_RelSample_2016_;
}
const float &hhtree::hh_mass_JESUp_RelSample_2016(){
	if(not hh_mass_JESUp_RelSample_2016_isLoaded){
		if(hh_mass_JESUp_RelSample_2016_branch != 0) hh_mass_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_RelSample_2016_branch does not exist!\n");}
		hh_mass_JESUp_RelSample_2016_isLoaded = true;
}
return hh_mass_JESUp_RelSample_2016_;
}
const float &hhtree::hh_pt_JESDown_RelSample_2016(){
	if(not hh_pt_JESDown_RelSample_2016_isLoaded){
		if(hh_pt_JESDown_RelSample_2016_branch != 0) hh_pt_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_RelSample_2016_branch does not exist!\n");}
		hh_pt_JESDown_RelSample_2016_isLoaded = true;
}
return hh_pt_JESDown_RelSample_2016_;
}
const float &hhtree::hh_eta_JESDown_RelSample_2016(){
	if(not hh_eta_JESDown_RelSample_2016_isLoaded){
		if(hh_eta_JESDown_RelSample_2016_branch != 0) hh_eta_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_RelSample_2016_branch does not exist!\n");}
		hh_eta_JESDown_RelSample_2016_isLoaded = true;
}
return hh_eta_JESDown_RelSample_2016_;
}
const float &hhtree::hh_mass_JESDown_RelSample_2016(){
	if(not hh_mass_JESDown_RelSample_2016_isLoaded){
		if(hh_mass_JESDown_RelSample_2016_branch != 0) hh_mass_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_RelSample_2016_branch does not exist!\n");}
		hh_mass_JESDown_RelSample_2016_isLoaded = true;
}
return hh_mass_JESDown_RelSample_2016_;
}

const float &hhtree::deltaEta_j1j2() 
{
 if(not deltaEta_j1j2_isLoaded)
 {
   if(deltaEta_j1j2_branch != 0) deltaEta_j1j2_branch->GetEntry(index);
   else
   {
     printf("branch deltaEta_j1j2_branch does not exist!\n");
//exit(1);
   }
   deltaEta_j1j2_isLoaded = true;
 }
 return deltaEta_j1j2_;
}

const float &hhtree::deltaPhi_j1j2() 
{
 if(not deltaPhi_j1j2_isLoaded)
 {
   if(deltaPhi_j1j2_branch != 0) deltaPhi_j1j2_branch->GetEntry(index);
   else
   {
     printf("branch deltaPhi_j1j2_branch does not exist!\n");
//exit(1);
   }
   deltaPhi_j1j2_isLoaded = true;
 }
 return deltaPhi_j1j2_;
}

const float &hhtree::deltaR_j1j2() 
{
 if(not deltaR_j1j2_isLoaded)
 {
   if(deltaR_j1j2_branch != 0) deltaR_j1j2_branch->GetEntry(index);
   else
   {
     printf("branch deltaR_j1j2_branch does not exist!\n");
//exit(1);
   }
   deltaR_j1j2_isLoaded = true;
 }
 return deltaR_j1j2_;
}

const float &hhtree::ptj2_over_ptj1() 
{
 if(not ptj2_over_ptj1_isLoaded)
 {
   if(ptj2_over_ptj1_branch != 0) ptj2_over_ptj1_branch->GetEntry(index);
   else
   {
     printf("branch ptj2_over_ptj1_branch does not exist!\n");
//exit(1);
   }
   ptj2_over_ptj1_isLoaded = true;
 }
 return ptj2_over_ptj1_;
}

const float &hhtree::mj2_over_mj1() 
{
 if(not mj2_over_mj1_isLoaded)
 {
   if(mj2_over_mj1_branch != 0) mj2_over_mj1_branch->GetEntry(index);
   else
   {
     printf("branch mj2_over_mj1_branch does not exist!\n");
//exit(1);
   }
   mj2_over_mj1_isLoaded = true;
 }
 return mj2_over_mj1_;
}

const int &hhtree::isVBFtag() 
{
 if(not isVBFtag_isLoaded)
 {
   if(isVBFtag_branch != 0) isVBFtag_branch->GetEntry(index);
   else
   {
     printf("branch isVBFtag_branch does not exist!\n");
//exit(1);
   }
   isVBFtag_isLoaded = true;
 }
 return isVBFtag_;
}

const float &hhtree::isVBFtag_JERUp(){
	if(not isVBFtag_JERUp_isLoaded){
		if(isVBFtag_JERUp_branch != 0) isVBFtag_JERUp_branch->GetEntry(index);
		else {printf("branch isVBFtag_JERUp_branch does not exist!\n");}
		isVBFtag_JERUp_isLoaded = true;
}
return isVBFtag_JERUp_;
}
const float &hhtree::isVBFtag_JERDown(){
	if(not isVBFtag_JERDown_isLoaded){
		if(isVBFtag_JERDown_branch != 0) isVBFtag_JERDown_branch->GetEntry(index);
		else {printf("branch isVBFtag_JERDown_branch does not exist!\n");}
		isVBFtag_JERDown_isLoaded = true;
}
return isVBFtag_JERDown_;
}
const float &hhtree::isVBFtag_JESUp(){
	if(not isVBFtag_JESUp_isLoaded){
		if(isVBFtag_JESUp_branch != 0) isVBFtag_JESUp_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_branch does not exist!\n");}
		isVBFtag_JESUp_isLoaded = true;
}
return isVBFtag_JESUp_;
}
const float &hhtree::isVBFtag_JESDown(){
	if(not isVBFtag_JESDown_isLoaded){
		if(isVBFtag_JESDown_branch != 0) isVBFtag_JESDown_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_branch does not exist!\n");}
		isVBFtag_JESDown_isLoaded = true;
}
return isVBFtag_JESDown_;
}
const float &hhtree::isVBFtag_JESUp_Abs(){
	if(not isVBFtag_JESUp_Abs_isLoaded){
		if(isVBFtag_JESUp_Abs_branch != 0) isVBFtag_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_Abs_branch does not exist!\n");}
		isVBFtag_JESUp_Abs_isLoaded = true;
}
return isVBFtag_JESUp_Abs_;
}
const float &hhtree::isVBFtag_JESDown_Abs(){
	if(not isVBFtag_JESDown_Abs_isLoaded){
		if(isVBFtag_JESDown_Abs_branch != 0) isVBFtag_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_Abs_branch does not exist!\n");}
		isVBFtag_JESDown_Abs_isLoaded = true;
}
return isVBFtag_JESDown_Abs_;
}
const float &hhtree::isVBFtag_JESUp_Abs_2016(){
	if(not isVBFtag_JESUp_Abs_2016_isLoaded){
		if(isVBFtag_JESUp_Abs_2016_branch != 0) isVBFtag_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_Abs_2016_branch does not exist!\n");}
		isVBFtag_JESUp_Abs_2016_isLoaded = true;
}
return isVBFtag_JESUp_Abs_2016_;
}
const float &hhtree::isVBFtag_JESDowb_Abs_2016(){
	if(not isVBFtag_JESDowb_Abs_2016_isLoaded){
		if(isVBFtag_JESDowb_Abs_2016_branch != 0) isVBFtag_JESDowb_Abs_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDowb_Abs_2016_branch does not exist!\n");}
		isVBFtag_JESDowb_Abs_2016_isLoaded = true;
}
return isVBFtag_JESDowb_Abs_2016_;
}
const float &hhtree::isVBFtag_JESUp_BBEC1(){
	if(not isVBFtag_JESUp_BBEC1_isLoaded){
		if(isVBFtag_JESUp_BBEC1_branch != 0) isVBFtag_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_BBEC1_branch does not exist!\n");}
		isVBFtag_JESUp_BBEC1_isLoaded = true;
}
return isVBFtag_JESUp_BBEC1_;
}
const float &hhtree::isVBFtag_JESDown_BBEC1(){
	if(not isVBFtag_JESDown_BBEC1_isLoaded){
		if(isVBFtag_JESDown_BBEC1_branch != 0) isVBFtag_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_BBEC1_branch does not exist!\n");}
		isVBFtag_JESDown_BBEC1_isLoaded = true;
}
return isVBFtag_JESDown_BBEC1_;
}
const float &hhtree::isVBFtag_JESUp_BBEC1_2016(){
	if(not isVBFtag_JESUp_BBEC1_2016_isLoaded){
		if(isVBFtag_JESUp_BBEC1_2016_branch != 0) isVBFtag_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_BBEC1_2016_branch does not exist!\n");}
		isVBFtag_JESUp_BBEC1_2016_isLoaded = true;
}
return isVBFtag_JESUp_BBEC1_2016_;
}
const float &hhtree::isVBFtag_JESDown_BBEC1_2016(){
	if(not isVBFtag_JESDown_BBEC1_2016_isLoaded){
		if(isVBFtag_JESDown_BBEC1_2016_branch != 0) isVBFtag_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_BBEC1_2016_branch does not exist!\n");}
		isVBFtag_JESDown_BBEC1_2016_isLoaded = true;
}
return isVBFtag_JESDown_BBEC1_2016_;
}
const float &hhtree::isVBFtag_JESUp_EC2(){
	if(not isVBFtag_JESUp_EC2_isLoaded){
		if(isVBFtag_JESUp_EC2_branch != 0) isVBFtag_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_EC2_branch does not exist!\n");}
		isVBFtag_JESUp_EC2_isLoaded = true;
}
return isVBFtag_JESUp_EC2_;
}
const float &hhtree::isVBFtag_JESDown_EC2(){
	if(not isVBFtag_JESDown_EC2_isLoaded){
		if(isVBFtag_JESDown_EC2_branch != 0) isVBFtag_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_EC2_branch does not exist!\n");}
		isVBFtag_JESDown_EC2_isLoaded = true;
}
return isVBFtag_JESDown_EC2_;
}
const float &hhtree::isVBFtag_JESUp_EC2_2016(){
	if(not isVBFtag_JESUp_EC2_2016_isLoaded){
		if(isVBFtag_JESUp_EC2_2016_branch != 0) isVBFtag_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_EC2_2016_branch does not exist!\n");}
		isVBFtag_JESUp_EC2_2016_isLoaded = true;
}
return isVBFtag_JESUp_EC2_2016_;
}
const float &hhtree::isVBFtag_JESDown_EC2_2016(){
	if(not isVBFtag_JESDown_EC2_2016_isLoaded){
		if(isVBFtag_JESDown_EC2_2016_branch != 0) isVBFtag_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_EC2_2016_branch does not exist!\n");}
		isVBFtag_JESDown_EC2_2016_isLoaded = true;
}
return isVBFtag_JESDown_EC2_2016_;
}
const float &hhtree::isVBFtag_JESUp_FlavQCD(){
	if(not isVBFtag_JESUp_FlavQCD_isLoaded){
		if(isVBFtag_JESUp_FlavQCD_branch != 0) isVBFtag_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_FlavQCD_branch does not exist!\n");}
		isVBFtag_JESUp_FlavQCD_isLoaded = true;
}
return isVBFtag_JESUp_FlavQCD_;
}
const float &hhtree::isVBFtag_JESDown_FlavQCD(){
	if(not isVBFtag_JESDown_FlavQCD_isLoaded){
		if(isVBFtag_JESDown_FlavQCD_branch != 0) isVBFtag_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_FlavQCD_branch does not exist!\n");}
		isVBFtag_JESDown_FlavQCD_isLoaded = true;
}
return isVBFtag_JESDown_FlavQCD_;
}
const float &hhtree::isVBFtag_JESUp_HF(){
	if(not isVBFtag_JESUp_HF_isLoaded){
		if(isVBFtag_JESUp_HF_branch != 0) isVBFtag_JESUp_HF_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_HF_branch does not exist!\n");}
		isVBFtag_JESUp_HF_isLoaded = true;
}
return isVBFtag_JESUp_HF_;
}
const float &hhtree::isVBFtag_JESDown_HF(){
	if(not isVBFtag_JESDown_HF_isLoaded){
		if(isVBFtag_JESDown_HF_branch != 0) isVBFtag_JESDown_HF_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_HF_branch does not exist!\n");}
		isVBFtag_JESDown_HF_isLoaded = true;
}
return isVBFtag_JESDown_HF_;
}
const float &hhtree::isVBFtag_JESUp_HF_2016(){
	if(not isVBFtag_JESUp_HF_2016_isLoaded){
		if(isVBFtag_JESUp_HF_2016_branch != 0) isVBFtag_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_HF_2016_branch does not exist!\n");}
		isVBFtag_JESUp_HF_2016_isLoaded = true;
}
return isVBFtag_JESUp_HF_2016_;
}
const float &hhtree::isVBFtag_JESDown_HF_2016(){
	if(not isVBFtag_JESDown_HF_2016_isLoaded){
		if(isVBFtag_JESDown_HF_2016_branch != 0) isVBFtag_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_HF_2016_branch does not exist!\n");}
		isVBFtag_JESDown_HF_2016_isLoaded = true;
}
return isVBFtag_JESDown_HF_2016_;
}
const float &hhtree::isVBFtag_JESUp_RelBal(){
	if(not isVBFtag_JESUp_RelBal_isLoaded){
		if(isVBFtag_JESUp_RelBal_branch != 0) isVBFtag_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_RelBal_branch does not exist!\n");}
		isVBFtag_JESUp_RelBal_isLoaded = true;
}
return isVBFtag_JESUp_RelBal_;
}
const float &hhtree::isVBFtag_JESDown_RelBal(){
	if(not isVBFtag_JESDown_RelBal_isLoaded){
		if(isVBFtag_JESDown_RelBal_branch != 0) isVBFtag_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_RelBal_branch does not exist!\n");}
		isVBFtag_JESDown_RelBal_isLoaded = true;
}
return isVBFtag_JESDown_RelBal_;
}
const float &hhtree::isVBFtag_JESUp_RelSample_2016(){
	if(not isVBFtag_JESUp_RelSample_2016_isLoaded){
		if(isVBFtag_JESUp_RelSample_2016_branch != 0) isVBFtag_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_RelSample_2016_branch does not exist!\n");}
		isVBFtag_JESUp_RelSample_2016_isLoaded = true;
}
return isVBFtag_JESUp_RelSample_2016_;
}
const float &hhtree::isVBFtag_JESDown_RelSample_2016(){
	if(not isVBFtag_JESDown_RelSample_2016_isLoaded){
		if(isVBFtag_JESDown_RelSample_2016_branch != 0) isVBFtag_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_RelSample_2016_branch does not exist!\n");}
		isVBFtag_JESDown_RelSample_2016_isLoaded = true;
}
return isVBFtag_JESDown_RelSample_2016_;
}

const float &hhtree::dijetmass() 
{
 if(not dijetmass_isLoaded)
 {
   if(dijetmass_branch != 0) dijetmass_branch->GetEntry(index);
   else
   {
     printf("branch dijetmass_branch does not exist!\n");
//exit(1);
   }
   dijetmass_isLoaded = true;
 }
 return dijetmass_;
}

const float &hhtree::vbfjet1Pt() 
{
 if(not vbfjet1Pt_isLoaded)
 {
   if(vbfjet1Pt_branch != 0) vbfjet1Pt_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet1Pt_branch does not exist!\n");
//exit(1);
   }
   vbfjet1Pt_isLoaded = true;
 }
 return vbfjet1Pt_;
}

const float &hhtree::vbfjet1Eta() 
{
 if(not vbfjet1Eta_isLoaded)
 {
   if(vbfjet1Eta_branch != 0) vbfjet1Eta_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet1Eta_branch does not exist!\n");
//exit(1);
   }
   vbfjet1Eta_isLoaded = true;
 }
 return vbfjet1Eta_;
}

const float &hhtree::vbfjet1Phi() 
{
 if(not vbfjet1Phi_isLoaded)
 {
   if(vbfjet1Phi_branch != 0) vbfjet1Phi_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet1Phi_branch does not exist!\n");
//exit(1);
   }
   vbfjet1Phi_isLoaded = true;
 }
 return vbfjet1Phi_;
}

const float &hhtree::vbfjet1Mass() 
{
 if(not vbfjet1Mass_isLoaded)
 {
   if(vbfjet1Mass_branch != 0) vbfjet1Mass_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet1Mass_branch does not exist!\n");
//exit(1);
   }
   vbfjet1Mass_isLoaded = true;
 }
 return vbfjet1Mass_;
}

const float &hhtree::vbffatJet1Pt() 
{
 if(not vbffatJet1Pt_isLoaded)
 {
   if(vbffatJet1Pt_branch != 0) vbffatJet1Pt_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet1Pt_branch does not exist!\n");
//exit(1);
   }
   vbffatJet1Pt_isLoaded = true;
 }
 return vbffatJet1Pt_;
}

const float &hhtree::vbffatJet1Eta() 
{
 if(not vbffatJet1Eta_isLoaded)
 {
   if(vbffatJet1Eta_branch != 0) vbffatJet1Eta_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet1Eta_branch does not exist!\n");
//exit(1);
   }
   vbffatJet1Eta_isLoaded = true;
 }
 return vbffatJet1Eta_;
}

const float &hhtree::vbffatJet1Phi() 
{
 if(not vbffatJet1Phi_isLoaded)
 {
   if(vbffatJet1Phi_branch != 0) vbffatJet1Phi_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet1Phi_branch does not exist!\n");
//exit(1);
   }
   vbffatJet1Phi_isLoaded = true;
 }
 return vbffatJet1Phi_;
}

const float &hhtree::vbffatJet1PNetXbb() 
{
 if(not vbffatJet1PNetXbb_isLoaded)
 {
   if(vbffatJet1PNetXbb_branch != 0) vbffatJet1PNetXbb_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet1PNetXbb_branch does not exist!\n");
//exit(1);
   }
   vbffatJet1PNetXbb_isLoaded = true;
 }
 return vbffatJet1PNetXbb_;
}

const float &hhtree::vbfjet2Pt() 
{
 if(not vbfjet2Pt_isLoaded)
 {
   if(vbfjet2Pt_branch != 0) vbfjet2Pt_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet2Pt_branch does not exist!\n");
//exit(1);
   }
   vbfjet2Pt_isLoaded = true;
 }
 return vbfjet2Pt_;
}

const float &hhtree::vbfjet2Eta() 
{
 if(not vbfjet2Eta_isLoaded)
 {
   if(vbfjet2Eta_branch != 0) vbfjet2Eta_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet2Eta_branch does not exist!\n");
//exit(1);
   }
   vbfjet2Eta_isLoaded = true;
 }
 return vbfjet2Eta_;
}

const float &hhtree::vbfjet2Phi() 
{
 if(not vbfjet2Phi_isLoaded)
 {
   if(vbfjet2Phi_branch != 0) vbfjet2Phi_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet2Phi_branch does not exist!\n");
//exit(1);
   }
   vbfjet2Phi_isLoaded = true;
 }
 return vbfjet2Phi_;
}

const float &hhtree::vbfjet2Mass() 
{
 if(not vbfjet2Mass_isLoaded)
 {
   if(vbfjet2Mass_branch != 0) vbfjet2Mass_branch->GetEntry(index);
   else
   {
     printf("branch vbfjet2Mass_branch does not exist!\n");
//exit(1);
   }
   vbfjet2Mass_isLoaded = true;
 }
 return vbfjet2Mass_;
}

const float &hhtree::vbffatJet2Pt() 
{
 if(not vbffatJet2Pt_isLoaded)
 {
   if(vbffatJet2Pt_branch != 0) vbffatJet2Pt_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet2Pt_branch does not exist!\n");
//exit(1);
   }
   vbffatJet2Pt_isLoaded = true;
 }
 return vbffatJet2Pt_;
}

const float &hhtree::vbffatJet2Eta() 
{
 if(not vbffatJet2Eta_isLoaded)
 {
   if(vbffatJet2Eta_branch != 0) vbffatJet2Eta_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet2Eta_branch does not exist!\n");
//exit(1);
   }
   vbffatJet2Eta_isLoaded = true;
 }
 return vbffatJet2Eta_;
}

const float &hhtree::vbffatJet2Phi() 
{
 if(not vbffatJet2Phi_isLoaded)
 {
   if(vbffatJet2Phi_branch != 0) vbffatJet2Phi_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet2Phi_branch does not exist!\n");
//exit(1);
   }
   vbffatJet2Phi_isLoaded = true;
 }
 return vbffatJet2Phi_;
}

const float &hhtree::vbffatJet2PNetXbb() 
{
 if(not vbffatJet2PNetXbb_isLoaded)
 {
   if(vbffatJet2PNetXbb_branch != 0) vbffatJet2PNetXbb_branch->GetEntry(index);
   else
   {
     printf("branch vbffatJet2PNetXbb_branch does not exist!\n");
//exit(1);
   }
   vbffatJet2PNetXbb_isLoaded = true;
 }
 return vbffatJet2PNetXbb_;
}

const float &hhtree::jet1Pt() 
{
 if(not jet1Pt_isLoaded)
 {
   if(jet1Pt_branch != 0) jet1Pt_branch->GetEntry(index);
   else
   {
     printf("branch jet1Pt_branch does not exist!\n");
//exit(1);
   }
   jet1Pt_isLoaded = true;
 }
 return jet1Pt_;
}

const float &hhtree::jet1Eta() 
{
 if(not jet1Eta_isLoaded)
 {
   if(jet1Eta_branch != 0) jet1Eta_branch->GetEntry(index);
   else
   {
     printf("branch jet1Eta_branch does not exist!\n");
//exit(1);
   }
   jet1Eta_isLoaded = true;
 }
 return jet1Eta_;
}

const float &hhtree::jet1Phi() 
{
 if(not jet1Phi_isLoaded)
 {
   if(jet1Phi_branch != 0) jet1Phi_branch->GetEntry(index);
   else
   {
     printf("branch jet1Phi_branch does not exist!\n");
//exit(1);
   }
   jet1Phi_isLoaded = true;
 }
 return jet1Phi_;
}

const float &hhtree::jet2Pt() 
{
 if(not jet2Pt_isLoaded)
 {
   if(jet2Pt_branch != 0) jet2Pt_branch->GetEntry(index);
   else
   {
     printf("branch jet2Pt_branch does not exist!\n");
//exit(1);
   }
   jet2Pt_isLoaded = true;
 }
 return jet2Pt_;
}

const float &hhtree::jet2Eta() 
{
 if(not jet2Eta_isLoaded)
 {
   if(jet2Eta_branch != 0) jet2Eta_branch->GetEntry(index);
   else
   {
     printf("branch jet2Eta_branch does not exist!\n");
//exit(1);
   }
   jet2Eta_isLoaded = true;
 }
 return jet2Eta_;
}

const float &hhtree::jet2Phi() 
{
 if(not jet2Phi_isLoaded)
 {
   if(jet2Phi_branch != 0) jet2Phi_branch->GetEntry(index);
   else
   {
     printf("branch jet2Phi_branch does not exist!\n");
//exit(1);
   }
   jet2Phi_isLoaded = true;
 }
 return jet2Phi_;
}

const float &hhtree::jet3Pt() 
{
 if(not jet3Pt_isLoaded)
 {
   if(jet3Pt_branch != 0) jet3Pt_branch->GetEntry(index);
   else
   {
     printf("branch jet3Pt_branch does not exist!\n");
//exit(1);
   }
   jet3Pt_isLoaded = true;
 }
 return jet3Pt_;
}

const float &hhtree::jet3Eta() 
{
 if(not jet3Eta_isLoaded)
 {
   if(jet3Eta_branch != 0) jet3Eta_branch->GetEntry(index);
   else
   {
     printf("branch jet3Eta_branch does not exist!\n");
//exit(1);
   }
   jet3Eta_isLoaded = true;
 }
 return jet3Eta_;
}

const float &hhtree::jet3Phi() 
{
 if(not jet3Phi_isLoaded)
 {
   if(jet3Phi_branch != 0) jet3Phi_branch->GetEntry(index);
   else
   {
     printf("branch jet3Phi_branch does not exist!\n");
//exit(1);
   }
   jet3Phi_isLoaded = true;
 }
 return jet3Phi_;
}

const float &hhtree::jet4Pt() 
{
 if(not jet4Pt_isLoaded)
 {
   if(jet4Pt_branch != 0) jet4Pt_branch->GetEntry(index);
   else
   {
     printf("branch jet4Pt_branch does not exist!\n");
//exit(1);
   }
   jet4Pt_isLoaded = true;
 }
 return jet4Pt_;
}

const float &hhtree::jet4Eta() 
{
 if(not jet4Eta_isLoaded)
 {
   if(jet4Eta_branch != 0) jet4Eta_branch->GetEntry(index);
   else
   {
     printf("branch jet4Eta_branch does not exist!\n");
//exit(1);
   }
   jet4Eta_isLoaded = true;
 }
 return jet4Eta_;
}

const float &hhtree::jet4Phi() 
{
 if(not jet4Phi_isLoaded)
 {
   if(jet4Phi_branch != 0) jet4Phi_branch->GetEntry(index);
   else
   {
     printf("branch jet4Phi_branch does not exist!\n");
//exit(1);
   }
   jet4Phi_isLoaded = true;
 }
 return jet4Phi_;
}

const float &hhtree::lep1Pt() 
{
 if(not lep1Pt_isLoaded)
 {
   if(lep1Pt_branch != 0) lep1Pt_branch->GetEntry(index);
   else
   {
     printf("branch lep1Pt_branch does not exist!\n");
//exit(1);
   }
   lep1Pt_isLoaded = true;
 }
 return lep1Pt_;
}

const float &hhtree::lep1Eta() 
{
 if(not lep1Eta_isLoaded)
 {
   if(lep1Eta_branch != 0) lep1Eta_branch->GetEntry(index);
   else
   {
     printf("branch lep1Eta_branch does not exist!\n");
//exit(1);
   }
   lep1Eta_isLoaded = true;
 }
 return lep1Eta_;
}

const float &hhtree::lep1Phi() 
{
 if(not lep1Phi_isLoaded)
 {
   if(lep1Phi_branch != 0) lep1Phi_branch->GetEntry(index);
   else
   {
     printf("branch lep1Phi_branch does not exist!\n");
//exit(1);
   }
   lep1Phi_isLoaded = true;
 }
 return lep1Phi_;
}

const int &hhtree::lep1Id() 
{
 if(not lep1Id_isLoaded)
 {
   if(lep1Id_branch != 0) lep1Id_branch->GetEntry(index);
   else
   {
     printf("branch lep1Id_branch does not exist!\n");
//exit(1);
   }
   lep1Id_isLoaded = true;
 }
 return lep1Id_;
}

const float &hhtree::lep2Pt() 
{
 if(not lep2Pt_isLoaded)
 {
   if(lep2Pt_branch != 0) lep2Pt_branch->GetEntry(index);
   else
   {
     printf("branch lep2Pt_branch does not exist!\n");
//exit(1);
   }
   lep2Pt_isLoaded = true;
 }
 return lep2Pt_;
}

const float &hhtree::lep2Eta() 
{
 if(not lep2Eta_isLoaded)
 {
   if(lep2Eta_branch != 0) lep2Eta_branch->GetEntry(index);
   else
   {
     printf("branch lep2Eta_branch does not exist!\n");
//exit(1);
   }
   lep2Eta_isLoaded = true;
 }
 return lep2Eta_;
}

const float &hhtree::lep2Phi() 
{
 if(not lep2Phi_isLoaded)
 {
   if(lep2Phi_branch != 0) lep2Phi_branch->GetEntry(index);
   else
   {
     printf("branch lep2Phi_branch does not exist!\n");
//exit(1);
   }
   lep2Phi_isLoaded = true;
 }
 return lep2Phi_;
}

const int &hhtree::lep2Id() 
{
 if(not lep2Id_isLoaded)
 {
   if(lep2Id_branch != 0) lep2Id_branch->GetEntry(index);
   else
   {
     printf("branch lep2Id_branch does not exist!\n");
//exit(1);
   }
   lep2Id_isLoaded = true;
 }
 return lep2Id_;
}

const float &hhtree::genHiggs1Pt() 
{
 if(not genHiggs1Pt_isLoaded)
 {
   if(genHiggs1Pt_branch != 0) genHiggs1Pt_branch->GetEntry(index);
   else
   {
     printf("branch genHiggs1Pt_branch does not exist!\n");
//exit(1);
   }
   genHiggs1Pt_isLoaded = true;
 }
 return genHiggs1Pt_;
}

const float &hhtree::genHiggs1Eta() 
{
 if(not genHiggs1Eta_isLoaded)
 {
   if(genHiggs1Eta_branch != 0) genHiggs1Eta_branch->GetEntry(index);
   else
   {
     printf("branch genHiggs1Eta_branch does not exist!\n");
//exit(1);
   }
   genHiggs1Eta_isLoaded = true;
 }
 return genHiggs1Eta_;
}

const float &hhtree::genHiggs1Phi() 
{
 if(not genHiggs1Phi_isLoaded)
 {
   if(genHiggs1Phi_branch != 0) genHiggs1Phi_branch->GetEntry(index);
   else
   {
     printf("branch genHiggs1Phi_branch does not exist!\n");
//exit(1);
   }
   genHiggs1Phi_isLoaded = true;
 }
 return genHiggs1Phi_;
}

const float &hhtree::genHiggs2Pt() 
{
 if(not genHiggs2Pt_isLoaded)
 {
   if(genHiggs2Pt_branch != 0) genHiggs2Pt_branch->GetEntry(index);
   else
   {
     printf("branch genHiggs2Pt_branch does not exist!\n");
//exit(1);
   }
   genHiggs2Pt_isLoaded = true;
 }
 return genHiggs2Pt_;
}

const float &hhtree::genHiggs2Eta() 
{
 if(not genHiggs2Eta_isLoaded)
 {
   if(genHiggs2Eta_branch != 0) genHiggs2Eta_branch->GetEntry(index);
   else
   {
     printf("branch genHiggs2Eta_branch does not exist!\n");
//exit(1);
   }
   genHiggs2Eta_isLoaded = true;
 }
 return genHiggs2Eta_;
}

const float &hhtree::genHiggs2Phi() 
{
 if(not genHiggs2Phi_isLoaded)
 {
   if(genHiggs2Phi_branch != 0) genHiggs2Phi_branch->GetEntry(index);
   else
   {
     printf("branch genHiggs2Phi_branch does not exist!\n");
//exit(1);
   }
   genHiggs2Phi_isLoaded = true;
 }
 return genHiggs2Phi_;
}

const float &hhtree::puWeight() 
{
 if(not puWeight_isLoaded)
 {
   if(puWeight_branch != 0) puWeight_branch->GetEntry(index);
   else
   {
     printf("branch puWeight_branch does not exist!\n");
//exit(1);
   }
   puWeight_isLoaded = true;
 }
 return puWeight_;
}

const float &hhtree::puWeightUp() 
{
 if(not puWeightUp_isLoaded)
 {
   if(puWeightUp_branch != 0) puWeightUp_branch->GetEntry(index);
   else
   {
     printf("branch puWeightUp_branch does not exist!\n");
//exit(1);
   }
   puWeightUp_isLoaded = true;
 }
 return puWeightUp_;
}

const float &hhtree::puWeightDown() 
{
 if(not puWeightDown_isLoaded)
 {
   if(puWeightDown_branch != 0) puWeightDown_branch->GetEntry(index);
   else
   {
     printf("branch puWeightDown_branch does not exist!\n");
//exit(1);
   }
   puWeightDown_isLoaded = true;
 }
 return puWeightDown_;
}

const float &hhtree::xsecWeight() 
{
 if(not xsecWeight_isLoaded)
 {
   if(xsecWeight_branch != 0) xsecWeight_branch->GetEntry(index);
   else
   {
     printf("branch xsecWeight_branch does not exist!\n");
     exit(1);
   }
   xsecWeight_isLoaded = true;
 }
 return xsecWeight_;
}
const float * hhtree::LHEScaleWeightNormNew()
{
 if(not LHEScaleWeightNormNew_isLoaded)
 {
   if(LHEScaleWeightNormNew_branch != 0) LHEScaleWeightNormNew_branch->GetEntry(index);
   else
   {
     //printf("branch LHEScaleWeightNormNew_branch does not exist!\n");
//exit(1);
   }
   LHEScaleWeightNormNew_isLoaded = true;
 }
 return LHEScaleWeightNormNew_;
}

const float * hhtree::LHEPdfWeightNorm()
{
 if(not LHEPdfWeightNorm_isLoaded)
 {
   if(LHEPdfWeightNorm_branch != 0) LHEPdfWeightNorm_branch->GetEntry(index);
   else
   {
     //printf("branch LHEPdfWeightNorm_branch does not exist!\n");
//exit(1);
   }
   LHEPdfWeightNorm_isLoaded = true;
 }
 return LHEPdfWeightNorm_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2() 
{
 if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_isLoaded)
 {
   if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_branch->GetEntry(index);
   else
   {
     printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_branch does not exist!\n");
//exit(1);
   }
   disc_qcd_and_ttbar_Run2_enhanced_v8p2_isLoaded = true;
 }
 return disc_qcd_and_ttbar_Run2_enhanced_v8p2_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp() 
{
 if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_isLoaded)
 {
   if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_branch->GetEntry(index);
   else
   {
     printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_branch does not exist!\n");
//exit(1);
   }
   disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_isLoaded = true;
 }
 return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown() 
{
 if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_isLoaded)
 {
   if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_branch->GetEntry(index);
   else
   {
     printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_branch does not exist!\n");
//exit(1);
   }
   disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_isLoaded = true;
 }
 return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp() 
{
 if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_isLoaded)
 {
   if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_branch->GetEntry(index);
   else
   {
     printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_branch does not exist!\n");
//exit(1);
   }
   disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_isLoaded = true;
 }
 return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSUp_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown() 
{
 if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_isLoaded)
 {
   if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_branch->GetEntry(index);
   else
   {
     printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_branch does not exist!\n");
//exit(1);
   }
   disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_isLoaded = true;
 }
 return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMSDown_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp() 
{
 if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_isLoaded)
 {
   if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_branch->GetEntry(index);
   else
   {
     printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_branch does not exist!\n");
//exit(1);
   }
   disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_isLoaded = true;
 }
 return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRUp_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown() 
{
 if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_isLoaded)
 {
   if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_branch->GetEntry(index);
   else
   {
     printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_branch does not exist!\n");
//exit(1);
   }
   disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_isLoaded = true;
 }
 return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JMRDown_;
}

const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_FlavQCD_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_FlavQCD_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelBal_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelBal_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2016_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2016_;
}


const float &hhtree::fatJet1Pt_JESUp_Abs_2017(){
	if(not fatJet1Pt_JESUp_Abs_2017_isLoaded){
		if(fatJet1Pt_JESUp_Abs_2017_branch != 0) fatJet1Pt_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_Abs_2017_branch does not exist!\n");}
		fatJet1Pt_JESUp_Abs_2017_isLoaded = true;
}
return fatJet1Pt_JESUp_Abs_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_Abs_2017(){
	if(not fatJet1PtOverMHH_JESUp_Abs_2017_isLoaded){
		if(fatJet1PtOverMHH_JESUp_Abs_2017_branch != 0) fatJet1PtOverMHH_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_Abs_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_Abs_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_Abs_2017_;
}
const float &hhtree::fatJet1Pt_JESDowb_Abs_2017(){
	if(not fatJet1Pt_JESDowb_Abs_2017_isLoaded){
		if(fatJet1Pt_JESDowb_Abs_2017_branch != 0) fatJet1Pt_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDowb_Abs_2017_branch does not exist!\n");}
		fatJet1Pt_JESDowb_Abs_2017_isLoaded = true;
}
return fatJet1Pt_JESDowb_Abs_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESDowb_Abs_2017(){
	if(not fatJet1PtOverMHH_JESDowb_Abs_2017_isLoaded){
		if(fatJet1PtOverMHH_JESDowb_Abs_2017_branch != 0) fatJet1PtOverMHH_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDowb_Abs_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDowb_Abs_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESDowb_Abs_2017_;
}
const float &hhtree::fatJet1Pt_JESUp_BBEC1_2017(){
	if(not fatJet1Pt_JESUp_BBEC1_2017_isLoaded){
		if(fatJet1Pt_JESUp_BBEC1_2017_branch != 0) fatJet1Pt_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_BBEC1_2017_branch does not exist!\n");}
		fatJet1Pt_JESUp_BBEC1_2017_isLoaded = true;
}
return fatJet1Pt_JESUp_BBEC1_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_BBEC1_2017(){
	if(not fatJet1PtOverMHH_JESUp_BBEC1_2017_isLoaded){
		if(fatJet1PtOverMHH_JESUp_BBEC1_2017_branch != 0) fatJet1PtOverMHH_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_BBEC1_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_BBEC1_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_BBEC1_2017_;
}
const float &hhtree::fatJet1Pt_JESDown_BBEC1_2017(){
	if(not fatJet1Pt_JESDown_BBEC1_2017_isLoaded){
		if(fatJet1Pt_JESDown_BBEC1_2017_branch != 0) fatJet1Pt_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_BBEC1_2017_branch does not exist!\n");}
		fatJet1Pt_JESDown_BBEC1_2017_isLoaded = true;
}
return fatJet1Pt_JESDown_BBEC1_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_BBEC1_2017(){
	if(not fatJet1PtOverMHH_JESDown_BBEC1_2017_isLoaded){
		if(fatJet1PtOverMHH_JESDown_BBEC1_2017_branch != 0) fatJet1PtOverMHH_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_BBEC1_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_BBEC1_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_BBEC1_2017_;
}

const float &hhtree::fatJet1Pt_JESUp_EC2_2017(){
	if(not fatJet1Pt_JESUp_EC2_2017_isLoaded){
		if(fatJet1Pt_JESUp_EC2_2017_branch != 0) fatJet1Pt_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_EC2_2017_branch does not exist!\n");}
		fatJet1Pt_JESUp_EC2_2017_isLoaded = true;
}
return fatJet1Pt_JESUp_EC2_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_EC2_2017(){
	if(not fatJet1PtOverMHH_JESUp_EC2_2017_isLoaded){
		if(fatJet1PtOverMHH_JESUp_EC2_2017_branch != 0) fatJet1PtOverMHH_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_EC2_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_EC2_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_EC2_2017_;
}
const float &hhtree::fatJet1Pt_JESDown_EC2_2017(){
	if(not fatJet1Pt_JESDown_EC2_2017_isLoaded){
		if(fatJet1Pt_JESDown_EC2_2017_branch != 0) fatJet1Pt_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_EC2_2017_branch does not exist!\n");}
		fatJet1Pt_JESDown_EC2_2017_isLoaded = true;
}
return fatJet1Pt_JESDown_EC2_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_EC2_2017(){
	if(not fatJet1PtOverMHH_JESDown_EC2_2017_isLoaded){
		if(fatJet1PtOverMHH_JESDown_EC2_2017_branch != 0) fatJet1PtOverMHH_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_EC2_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_EC2_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_EC2_2017_;
}
const float &hhtree::fatJet1Pt_JESUp_HF_2017(){
	if(not fatJet1Pt_JESUp_HF_2017_isLoaded){
		if(fatJet1Pt_JESUp_HF_2017_branch != 0) fatJet1Pt_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_HF_2017_branch does not exist!\n");}
		fatJet1Pt_JESUp_HF_2017_isLoaded = true;
}
return fatJet1Pt_JESUp_HF_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_HF_2017(){
	if(not fatJet1PtOverMHH_JESUp_HF_2017_isLoaded){
		if(fatJet1PtOverMHH_JESUp_HF_2017_branch != 0) fatJet1PtOverMHH_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_HF_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_HF_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_HF_2017_;
}
const float &hhtree::fatJet1Pt_JESDown_HF_2017(){
	if(not fatJet1Pt_JESDown_HF_2017_isLoaded){
		if(fatJet1Pt_JESDown_HF_2017_branch != 0) fatJet1Pt_JESDown_HF_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_HF_2017_branch does not exist!\n");}
		fatJet1Pt_JESDown_HF_2017_isLoaded = true;
}
return fatJet1Pt_JESDown_HF_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_HF_2017(){
  if(not fatJet1PtOverMHH_JESDown_HF_2017_isLoaded){
    if(fatJet1PtOverMHH_JESDown_HF_2017_branch != 0) fatJet1PtOverMHH_JESDown_HF_2017_branch->GetEntry(index);
    else {printf("branch fatJet1PtOverMHH_JESDown_HF_2017_branch does not exist!\n");}
    fatJet1PtOverMHH_JESDown_HF_2017_isLoaded = true;
  }
  return fatJet1PtOverMHH_JESDown_HF_2017_;
}
const float &hhtree::fatJet1Pt_JESUp_RelSample_2017(){
	if(not fatJet1Pt_JESUp_RelSample_2017_isLoaded){
		if(fatJet1Pt_JESUp_RelSample_2017_branch != 0) fatJet1Pt_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_RelSample_2017_branch does not exist!\n");}
		fatJet1Pt_JESUp_RelSample_2017_isLoaded = true;
}
return fatJet1Pt_JESUp_RelSample_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_RelSample_2017(){
	if(not fatJet1PtOverMHH_JESUp_RelSample_2017_isLoaded){
		if(fatJet1PtOverMHH_JESUp_RelSample_2017_branch != 0) fatJet1PtOverMHH_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_RelSample_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_RelSample_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_RelSample_2017_;
}
const float &hhtree::fatJet1Pt_JESDown_RelSample_2017(){
	if(not fatJet1Pt_JESDown_RelSample_2017_isLoaded){
		if(fatJet1Pt_JESDown_RelSample_2017_branch != 0) fatJet1Pt_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_RelSample_2017_branch does not exist!\n");}
		fatJet1Pt_JESDown_RelSample_2017_isLoaded = true;
}
return fatJet1Pt_JESDown_RelSample_2017_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_RelSample_2017(){
	if(not fatJet1PtOverMHH_JESDown_RelSample_2017_isLoaded){
		if(fatJet1PtOverMHH_JESDown_RelSample_2017_branch != 0) fatJet1PtOverMHH_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_RelSample_2017_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_RelSample_2017_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_RelSample_2017_;
}

const float &hhtree::fatJet2Pt_JESUp_Abs_2017(){
	if(not fatJet2Pt_JESUp_Abs_2017_isLoaded){
		if(fatJet2Pt_JESUp_Abs_2017_branch != 0) fatJet2Pt_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_Abs_2017_branch does not exist!\n");}
		fatJet2Pt_JESUp_Abs_2017_isLoaded = true;
}
return fatJet2Pt_JESUp_Abs_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_Abs_2017(){
	if(not fatJet2PtOverMHH_JESUp_Abs_2017_isLoaded){
		if(fatJet2PtOverMHH_JESUp_Abs_2017_branch != 0) fatJet2PtOverMHH_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_Abs_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_Abs_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_Abs_2017_;
}
const float &hhtree::fatJet2Pt_JESDowb_Abs_2017(){
	if(not fatJet2Pt_JESDowb_Abs_2017_isLoaded){
		if(fatJet2Pt_JESDowb_Abs_2017_branch != 0) fatJet2Pt_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDowb_Abs_2017_branch does not exist!\n");}
		fatJet2Pt_JESDowb_Abs_2017_isLoaded = true;
}
return fatJet2Pt_JESDowb_Abs_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESDowb_Abs_2017(){
	if(not fatJet2PtOverMHH_JESDowb_Abs_2017_isLoaded){
		if(fatJet2PtOverMHH_JESDowb_Abs_2017_branch != 0) fatJet2PtOverMHH_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDowb_Abs_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDowb_Abs_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESDowb_Abs_2017_;
}
const float &hhtree::fatJet2Pt_JESUp_BBEC1_2017(){
	if(not fatJet2Pt_JESUp_BBEC1_2017_isLoaded){
		if(fatJet2Pt_JESUp_BBEC1_2017_branch != 0) fatJet2Pt_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_BBEC1_2017_branch does not exist!\n");}
		fatJet2Pt_JESUp_BBEC1_2017_isLoaded = true;
}
return fatJet2Pt_JESUp_BBEC1_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_BBEC1_2017(){
	if(not fatJet2PtOverMHH_JESUp_BBEC1_2017_isLoaded){
		if(fatJet2PtOverMHH_JESUp_BBEC1_2017_branch != 0) fatJet2PtOverMHH_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_BBEC1_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_BBEC1_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_BBEC1_2017_;
}
const float &hhtree::fatJet2Pt_JESDown_BBEC1_2017(){
	if(not fatJet2Pt_JESDown_BBEC1_2017_isLoaded){
		if(fatJet2Pt_JESDown_BBEC1_2017_branch != 0) fatJet2Pt_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_BBEC1_2017_branch does not exist!\n");}
		fatJet2Pt_JESDown_BBEC1_2017_isLoaded = true;
}
return fatJet2Pt_JESDown_BBEC1_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_BBEC1_2017(){
	if(not fatJet2PtOverMHH_JESDown_BBEC1_2017_isLoaded){
		if(fatJet2PtOverMHH_JESDown_BBEC1_2017_branch != 0) fatJet2PtOverMHH_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_BBEC1_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_BBEC1_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_BBEC1_2017_;
}

const float &hhtree::fatJet2Pt_JESUp_EC2_2017(){
	if(not fatJet2Pt_JESUp_EC2_2017_isLoaded){
		if(fatJet2Pt_JESUp_EC2_2017_branch != 0) fatJet2Pt_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_EC2_2017_branch does not exist!\n");}
		fatJet2Pt_JESUp_EC2_2017_isLoaded = true;
}
return fatJet2Pt_JESUp_EC2_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_EC2_2017(){
	if(not fatJet2PtOverMHH_JESUp_EC2_2017_isLoaded){
		if(fatJet2PtOverMHH_JESUp_EC2_2017_branch != 0) fatJet2PtOverMHH_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_EC2_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_EC2_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_EC2_2017_;
}
const float &hhtree::fatJet2Pt_JESDown_EC2_2017(){
	if(not fatJet2Pt_JESDown_EC2_2017_isLoaded){
		if(fatJet2Pt_JESDown_EC2_2017_branch != 0) fatJet2Pt_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_EC2_2017_branch does not exist!\n");}
		fatJet2Pt_JESDown_EC2_2017_isLoaded = true;
}
return fatJet2Pt_JESDown_EC2_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_EC2_2017(){
	if(not fatJet2PtOverMHH_JESDown_EC2_2017_isLoaded){
		if(fatJet2PtOverMHH_JESDown_EC2_2017_branch != 0) fatJet2PtOverMHH_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_EC2_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_EC2_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_EC2_2017_;
}
const float &hhtree::fatJet2Pt_JESUp_HF_2017(){
	if(not fatJet2Pt_JESUp_HF_2017_isLoaded){
		if(fatJet2Pt_JESUp_HF_2017_branch != 0) fatJet2Pt_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_HF_2017_branch does not exist!\n");}
		fatJet2Pt_JESUp_HF_2017_isLoaded = true;
}
return fatJet2Pt_JESUp_HF_2017_;
}

const float &hhtree::fatJet2PtOverMHH_JESUp_HF_2017(){
	if(not fatJet2PtOverMHH_JESUp_HF_2017_isLoaded){
		if(fatJet2PtOverMHH_JESUp_HF_2017_branch != 0) fatJet2PtOverMHH_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_HF_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_HF_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_HF_2017_;
}
const float &hhtree::fatJet2Pt_JESDown_HF_2017(){
	if(not fatJet2Pt_JESDown_HF_2017_isLoaded){
		if(fatJet2Pt_JESDown_HF_2017_branch != 0) fatJet2Pt_JESDown_HF_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_HF_2017_branch does not exist!\n");}
		fatJet2Pt_JESDown_HF_2017_isLoaded = true;
}
return fatJet2Pt_JESDown_HF_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_HF_2017(){
	if(not fatJet2PtOverMHH_JESDown_HF_2017_isLoaded){
		if(fatJet2PtOverMHH_JESDown_HF_2017_branch != 0) fatJet2PtOverMHH_JESDown_HF_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_HF_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_HF_2017_isLoaded = true;
}
	return fatJet2PtOverMHH_JESDown_HF_2017_;
}
const float &hhtree::fatJet2Pt_JESUp_RelSample_2017(){
	if(not fatJet2Pt_JESUp_RelSample_2017_isLoaded){
		if(fatJet2Pt_JESUp_RelSample_2017_branch != 0) fatJet2Pt_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_RelSample_2017_branch does not exist!\n");}
		fatJet2Pt_JESUp_RelSample_2017_isLoaded = true;
}
return fatJet2Pt_JESUp_RelSample_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_RelSample_2017(){
	if(not fatJet2PtOverMHH_JESUp_RelSample_2017_isLoaded){
		if(fatJet2PtOverMHH_JESUp_RelSample_2017_branch != 0) fatJet2PtOverMHH_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_RelSample_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_RelSample_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_RelSample_2017_;
}
const float &hhtree::fatJet2Pt_JESDown_RelSample_2017(){
	if(not fatJet2Pt_JESDown_RelSample_2017_isLoaded){
		if(fatJet2Pt_JESDown_RelSample_2017_branch != 0) fatJet2Pt_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_RelSample_2017_branch does not exist!\n");}
		fatJet2Pt_JESDown_RelSample_2017_isLoaded = true;
}
return fatJet2Pt_JESDown_RelSample_2017_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_RelSample_2017(){
	if(not fatJet2PtOverMHH_JESDown_RelSample_2017_isLoaded){
		if(fatJet2PtOverMHH_JESDown_RelSample_2017_branch != 0) fatJet2PtOverMHH_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_RelSample_2017_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_RelSample_2017_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_RelSample_2017_;
}
const float &hhtree::hh_pt_JESUp_Abs_2017(){
	if(not hh_pt_JESUp_Abs_2017_isLoaded){
		if(hh_pt_JESUp_Abs_2017_branch != 0) hh_pt_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_Abs_2017_branch does not exist!\n");}
		hh_pt_JESUp_Abs_2017_isLoaded = true;
}
return hh_pt_JESUp_Abs_2017_;
}
const float &hhtree::hh_eta_JESUp_Abs_2017(){
	if(not hh_eta_JESUp_Abs_2017_isLoaded){
		if(hh_eta_JESUp_Abs_2017_branch != 0) hh_eta_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_Abs_2017_branch does not exist!\n");}
		hh_eta_JESUp_Abs_2017_isLoaded = true;
}
return hh_eta_JESUp_Abs_2017_;
}
const float &hhtree::hh_mass_JESUp_Abs_2017(){
	if(not hh_mass_JESUp_Abs_2017_isLoaded){
		if(hh_mass_JESUp_Abs_2017_branch != 0) hh_mass_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_Abs_2017_branch does not exist!\n");}
		hh_mass_JESUp_Abs_2017_isLoaded = true;
}
return hh_mass_JESUp_Abs_2017_;
}
const float &hhtree::hh_pt_JESDowb_Abs_2017(){
	if(not hh_pt_JESDowb_Abs_2017_isLoaded){
		if(hh_pt_JESDowb_Abs_2017_branch != 0) hh_pt_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDowb_Abs_2017_branch does not exist!\n");}
		hh_pt_JESDowb_Abs_2017_isLoaded = true;
}
return hh_pt_JESDowb_Abs_2017_;
}
const float &hhtree::hh_eta_JESDowb_Abs_2017(){
	if(not hh_eta_JESDowb_Abs_2017_isLoaded){
		if(hh_eta_JESDowb_Abs_2017_branch != 0) hh_eta_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDowb_Abs_2017_branch does not exist!\n");}
		hh_eta_JESDowb_Abs_2017_isLoaded = true;
}
return hh_eta_JESDowb_Abs_2017_;
}
const float &hhtree::hh_mass_JESDowb_Abs_2017(){
	if(not hh_mass_JESDowb_Abs_2017_isLoaded){
		if(hh_mass_JESDowb_Abs_2017_branch != 0) hh_mass_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDowb_Abs_2017_branch does not exist!\n");}
		hh_mass_JESDowb_Abs_2017_isLoaded = true;
}
return hh_mass_JESDowb_Abs_2017_;
}
const float &hhtree::hh_pt_JESUp_BBEC1_2017(){
	if(not hh_pt_JESUp_BBEC1_2017_isLoaded){
		if(hh_pt_JESUp_BBEC1_2017_branch != 0) hh_pt_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_BBEC1_2017_branch does not exist!\n");}
		hh_pt_JESUp_BBEC1_2017_isLoaded = true;
}
return hh_pt_JESUp_BBEC1_2017_;
}
const float &hhtree::hh_eta_JESUp_BBEC1_2017(){
	if(not hh_eta_JESUp_BBEC1_2017_isLoaded){
		if(hh_eta_JESUp_BBEC1_2017_branch != 0) hh_eta_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_BBEC1_2017_branch does not exist!\n");}
		hh_eta_JESUp_BBEC1_2017_isLoaded = true;
}
return hh_eta_JESUp_BBEC1_2017_;
}
const float &hhtree::hh_mass_JESUp_BBEC1_2017(){
	if(not hh_mass_JESUp_BBEC1_2017_isLoaded){
		if(hh_mass_JESUp_BBEC1_2017_branch != 0) hh_mass_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_BBEC1_2017_branch does not exist!\n");}
		hh_mass_JESUp_BBEC1_2017_isLoaded = true;
}
return hh_mass_JESUp_BBEC1_2017_;
}
const float &hhtree::hh_pt_JESDown_BBEC1_2017(){
	if(not hh_pt_JESDown_BBEC1_2017_isLoaded){
		if(hh_pt_JESDown_BBEC1_2017_branch != 0) hh_pt_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_BBEC1_2017_branch does not exist!\n");}
		hh_pt_JESDown_BBEC1_2017_isLoaded = true;
}
return hh_pt_JESDown_BBEC1_2017_;
}
const float &hhtree::hh_eta_JESDown_BBEC1_2017(){
	if(not hh_eta_JESDown_BBEC1_2017_isLoaded){
		if(hh_eta_JESDown_BBEC1_2017_branch != 0) hh_eta_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_BBEC1_2017_branch does not exist!\n");}
		hh_eta_JESDown_BBEC1_2017_isLoaded = true;
}
return hh_eta_JESDown_BBEC1_2017_;
}
const float &hhtree::hh_mass_JESDown_BBEC1_2017(){
	if(not hh_mass_JESDown_BBEC1_2017_isLoaded){
		if(hh_mass_JESDown_BBEC1_2017_branch != 0) hh_mass_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_BBEC1_2017_branch does not exist!\n");}
		hh_mass_JESDown_BBEC1_2017_isLoaded = true;
}
return hh_mass_JESDown_BBEC1_2017_;
}
const float &hhtree::hh_pt_JESUp_EC2_2017(){
	if(not hh_pt_JESUp_EC2_2017_isLoaded){
		if(hh_pt_JESUp_EC2_2017_branch != 0) hh_pt_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_EC2_2017_branch does not exist!\n");}
		hh_pt_JESUp_EC2_2017_isLoaded = true;
}
return hh_pt_JESUp_EC2_2017_;
}
const float &hhtree::hh_eta_JESUp_EC2_2017(){
	if(not hh_eta_JESUp_EC2_2017_isLoaded){
		if(hh_eta_JESUp_EC2_2017_branch != 0) hh_eta_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_EC2_2017_branch does not exist!\n");}
		hh_eta_JESUp_EC2_2017_isLoaded = true;
}
return hh_eta_JESUp_EC2_2017_;
}
const float &hhtree::hh_mass_JESUp_EC2_2017(){
	if(not hh_mass_JESUp_EC2_2017_isLoaded){
		if(hh_mass_JESUp_EC2_2017_branch != 0) hh_mass_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_EC2_2017_branch does not exist!\n");}
		hh_mass_JESUp_EC2_2017_isLoaded = true;
}
return hh_mass_JESUp_EC2_2017_;
}
const float &hhtree::hh_pt_JESDown_EC2_2017(){
	if(not hh_pt_JESDown_EC2_2017_isLoaded){
		if(hh_pt_JESDown_EC2_2017_branch != 0) hh_pt_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_EC2_2017_branch does not exist!\n");}
		hh_pt_JESDown_EC2_2017_isLoaded = true;
}
return hh_pt_JESDown_EC2_2017_;
}
const float &hhtree::hh_eta_JESDown_EC2_2017(){
	if(not hh_eta_JESDown_EC2_2017_isLoaded){
		if(hh_eta_JESDown_EC2_2017_branch != 0) hh_eta_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_EC2_2017_branch does not exist!\n");}
		hh_eta_JESDown_EC2_2017_isLoaded = true;
}
return hh_eta_JESDown_EC2_2017_;
}
const float &hhtree::hh_mass_JESDown_EC2_2017(){
	if(not hh_mass_JESDown_EC2_2017_isLoaded){
		if(hh_mass_JESDown_EC2_2017_branch != 0) hh_mass_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_EC2_2017_branch does not exist!\n");}
		hh_mass_JESDown_EC2_2017_isLoaded = true;
}
return hh_mass_JESDown_EC2_2017_;
}
const float &hhtree::hh_pt_JESUp_HF_2017(){
	if(not hh_pt_JESUp_HF_2017_isLoaded){
		if(hh_pt_JESUp_HF_2017_branch != 0) hh_pt_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_HF_2017_branch does not exist!\n");}
		hh_pt_JESUp_HF_2017_isLoaded = true;
}
return hh_pt_JESUp_HF_2017_;
}
const float &hhtree::hh_eta_JESUp_HF_2017(){
	if(not hh_eta_JESUp_HF_2017_isLoaded){
		if(hh_eta_JESUp_HF_2017_branch != 0) hh_eta_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_HF_2017_branch does not exist!\n");}
		hh_eta_JESUp_HF_2017_isLoaded = true;
}
return hh_eta_JESUp_HF_2017_;
}
const float &hhtree::hh_mass_JESUp_HF_2017(){
	if(not hh_mass_JESUp_HF_2017_isLoaded){
		if(hh_mass_JESUp_HF_2017_branch != 0) hh_mass_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_HF_2017_branch does not exist!\n");}
		hh_mass_JESUp_HF_2017_isLoaded = true;
}
return hh_mass_JESUp_HF_2017_;
}
const float &hhtree::hh_pt_JESDown_HF_2017(){
	if(not hh_pt_JESDown_HF_2017_isLoaded){
		if(hh_pt_JESDown_HF_2017_branch != 0) hh_pt_JESDown_HF_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_HF_2017_branch does not exist!\n");}
		hh_pt_JESDown_HF_2017_isLoaded = true;
}
return hh_pt_JESDown_HF_2017_;
}
const float &hhtree::hh_eta_JESDown_HF_2017(){
	if(not hh_eta_JESDown_HF_2017_isLoaded){
		if(hh_eta_JESDown_HF_2017_branch != 0) hh_eta_JESDown_HF_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_HF_2017_branch does not exist!\n");}
		hh_eta_JESDown_HF_2017_isLoaded = true;
}
return hh_eta_JESDown_HF_2017_;
}
const float &hhtree::hh_mass_JESDown_HF_2017(){
	if(not hh_mass_JESDown_HF_2017_isLoaded){
		if(hh_mass_JESDown_HF_2017_branch != 0) hh_mass_JESDown_HF_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_HF_2017_branch does not exist!\n");}
		hh_mass_JESDown_HF_2017_isLoaded = true;
}
return hh_mass_JESDown_HF_2017_;
}
const float &hhtree::hh_pt_JESUp_RelSample_2017(){
	if(not hh_pt_JESUp_RelSample_2017_isLoaded){
		if(hh_pt_JESUp_RelSample_2017_branch != 0) hh_pt_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_RelSample_2017_branch does not exist!\n");}
		hh_pt_JESUp_RelSample_2017_isLoaded = true;
}
return hh_pt_JESUp_RelSample_2017_;
}
const float &hhtree::hh_eta_JESUp_RelSample_2017(){
	if(not hh_eta_JESUp_RelSample_2017_isLoaded){
		if(hh_eta_JESUp_RelSample_2017_branch != 0) hh_eta_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_RelSample_2017_branch does not exist!\n");}
		hh_eta_JESUp_RelSample_2017_isLoaded = true;
}
return hh_eta_JESUp_RelSample_2017_;
}
const float &hhtree::hh_mass_JESUp_RelSample_2017(){
	if(not hh_mass_JESUp_RelSample_2017_isLoaded){
		if(hh_mass_JESUp_RelSample_2017_branch != 0) hh_mass_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_RelSample_2017_branch does not exist!\n");}
		hh_mass_JESUp_RelSample_2017_isLoaded = true;
}
return hh_mass_JESUp_RelSample_2017_;
}
const float &hhtree::hh_pt_JESDown_RelSample_2017(){
	if(not hh_pt_JESDown_RelSample_2017_isLoaded){
		if(hh_pt_JESDown_RelSample_2017_branch != 0) hh_pt_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_RelSample_2017_branch does not exist!\n");}
		hh_pt_JESDown_RelSample_2017_isLoaded = true;
}
return hh_pt_JESDown_RelSample_2017_;
}
const float &hhtree::hh_eta_JESDown_RelSample_2017(){
	if(not hh_eta_JESDown_RelSample_2017_isLoaded){
		if(hh_eta_JESDown_RelSample_2017_branch != 0) hh_eta_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_RelSample_2017_branch does not exist!\n");}
		hh_eta_JESDown_RelSample_2017_isLoaded = true;
}
return hh_eta_JESDown_RelSample_2017_;
}
const float &hhtree::hh_mass_JESDown_RelSample_2017(){
	if(not hh_mass_JESDown_RelSample_2017_isLoaded){
		if(hh_mass_JESDown_RelSample_2017_branch != 0) hh_mass_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_RelSample_2017_branch does not exist!\n");}
		hh_mass_JESDown_RelSample_2017_isLoaded = true;
}
return hh_mass_JESDown_RelSample_2017_;
}
const float &hhtree::isVBFtag_JESUp_Abs_2017(){
	if(not isVBFtag_JESUp_Abs_2017_isLoaded){
		if(isVBFtag_JESUp_Abs_2017_branch != 0) isVBFtag_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_Abs_2017_branch does not exist!\n");}
		isVBFtag_JESUp_Abs_2017_isLoaded = true;
}
return isVBFtag_JESUp_Abs_2017_;
}
const float &hhtree::isVBFtag_JESDowb_Abs_2017(){
	if(not isVBFtag_JESDowb_Abs_2017_isLoaded){
		if(isVBFtag_JESDowb_Abs_2017_branch != 0) isVBFtag_JESDowb_Abs_2017_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDowb_Abs_2017_branch does not exist!\n");}
		isVBFtag_JESDowb_Abs_2017_isLoaded = true;
}
return isVBFtag_JESDowb_Abs_2017_;
}
const float &hhtree::isVBFtag_JESUp_EC2_2017(){
	if(not isVBFtag_JESUp_EC2_2017_isLoaded){
		if(isVBFtag_JESUp_EC2_2017_branch != 0) isVBFtag_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_EC2_2017_branch does not exist!\n");}
		isVBFtag_JESUp_EC2_2017_isLoaded = true;
}
return isVBFtag_JESUp_EC2_2017_;
}
const float &hhtree::isVBFtag_JESDown_EC2_2017(){
	if(not isVBFtag_JESDown_EC2_2017_isLoaded){
		if(isVBFtag_JESDown_EC2_2017_branch != 0) isVBFtag_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_EC2_2017_branch does not exist!\n");}
		isVBFtag_JESDown_EC2_2017_isLoaded = true;
}
return isVBFtag_JESDown_EC2_2017_;
}
const float &hhtree::isVBFtag_JESUp_RelSample_2017(){
	if(not isVBFtag_JESUp_RelSample_2017_isLoaded){
		if(isVBFtag_JESUp_RelSample_2017_branch != 0) isVBFtag_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_RelSample_2017_branch does not exist!\n");}
		isVBFtag_JESUp_RelSample_2017_isLoaded = true;
}
return isVBFtag_JESUp_RelSample_2017_;
}
const float &hhtree::isVBFtag_JESDown_RelSample_2017(){
	if(not isVBFtag_JESDown_RelSample_2017_isLoaded){
		if(isVBFtag_JESDown_RelSample_2017_branch != 0) isVBFtag_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_RelSample_2017_branch does not exist!\n");}
		isVBFtag_JESDown_RelSample_2017_isLoaded = true;
}
return isVBFtag_JESDown_RelSample_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2017_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2017_;
}



const float &hhtree::fatJet1Pt_JESUp_Abs_2018(){
	if(not fatJet1Pt_JESUp_Abs_2018_isLoaded){
		if(fatJet1Pt_JESUp_Abs_2018_branch != 0) fatJet1Pt_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_Abs_2018_branch does not exist!\n");}
		fatJet1Pt_JESUp_Abs_2018_isLoaded = true;
}
return fatJet1Pt_JESUp_Abs_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_Abs_2018(){
	if(not fatJet1PtOverMHH_JESUp_Abs_2018_isLoaded){
		if(fatJet1PtOverMHH_JESUp_Abs_2018_branch != 0) fatJet1PtOverMHH_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_Abs_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_Abs_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_Abs_2018_;
}
const float &hhtree::fatJet1Pt_JESDowb_Abs_2018(){
	if(not fatJet1Pt_JESDowb_Abs_2018_isLoaded){
		if(fatJet1Pt_JESDowb_Abs_2018_branch != 0) fatJet1Pt_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDowb_Abs_2018_branch does not exist!\n");}
		fatJet1Pt_JESDowb_Abs_2018_isLoaded = true;
}
return fatJet1Pt_JESDowb_Abs_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESDowb_Abs_2018(){
	if(not fatJet1PtOverMHH_JESDowb_Abs_2018_isLoaded){
		if(fatJet1PtOverMHH_JESDowb_Abs_2018_branch != 0) fatJet1PtOverMHH_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDowb_Abs_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDowb_Abs_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESDowb_Abs_2018_;
}
const float &hhtree::fatJet1Pt_JESUp_BBEC1_2018(){
	if(not fatJet1Pt_JESUp_BBEC1_2018_isLoaded){
		if(fatJet1Pt_JESUp_BBEC1_2018_branch != 0) fatJet1Pt_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_BBEC1_2018_branch does not exist!\n");}
		fatJet1Pt_JESUp_BBEC1_2018_isLoaded = true;
}
return fatJet1Pt_JESUp_BBEC1_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_BBEC1_2018(){
	if(not fatJet1PtOverMHH_JESUp_BBEC1_2018_isLoaded){
		if(fatJet1PtOverMHH_JESUp_BBEC1_2018_branch != 0) fatJet1PtOverMHH_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_BBEC1_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_BBEC1_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_BBEC1_2018_;
}
const float &hhtree::fatJet1Pt_JESDown_BBEC1_2018(){
	if(not fatJet1Pt_JESDown_BBEC1_2018_isLoaded){
		if(fatJet1Pt_JESDown_BBEC1_2018_branch != 0) fatJet1Pt_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_BBEC1_2018_branch does not exist!\n");}
		fatJet1Pt_JESDown_BBEC1_2018_isLoaded = true;
}
return fatJet1Pt_JESDown_BBEC1_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_BBEC1_2018(){
	if(not fatJet1PtOverMHH_JESDown_BBEC1_2018_isLoaded){
		if(fatJet1PtOverMHH_JESDown_BBEC1_2018_branch != 0) fatJet1PtOverMHH_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_BBEC1_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_BBEC1_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_BBEC1_2018_;
}

const float &hhtree::fatJet1Pt_JESUp_EC2_2018(){
	if(not fatJet1Pt_JESUp_EC2_2018_isLoaded){
		if(fatJet1Pt_JESUp_EC2_2018_branch != 0) fatJet1Pt_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_EC2_2018_branch does not exist!\n");}
		fatJet1Pt_JESUp_EC2_2018_isLoaded = true;
}
return fatJet1Pt_JESUp_EC2_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_EC2_2018(){
	if(not fatJet1PtOverMHH_JESUp_EC2_2018_isLoaded){
		if(fatJet1PtOverMHH_JESUp_EC2_2018_branch != 0) fatJet1PtOverMHH_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_EC2_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_EC2_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_EC2_2018_;
}
const float &hhtree::fatJet1Pt_JESDown_EC2_2018(){
	if(not fatJet1Pt_JESDown_EC2_2018_isLoaded){
		if(fatJet1Pt_JESDown_EC2_2018_branch != 0) fatJet1Pt_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_EC2_2018_branch does not exist!\n");}
		fatJet1Pt_JESDown_EC2_2018_isLoaded = true;
}
return fatJet1Pt_JESDown_EC2_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_EC2_2018(){
	if(not fatJet1PtOverMHH_JESDown_EC2_2018_isLoaded){
		if(fatJet1PtOverMHH_JESDown_EC2_2018_branch != 0) fatJet1PtOverMHH_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_EC2_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_EC2_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_EC2_2018_;
}
const float &hhtree::fatJet1Pt_JESUp_HF_2018(){
	if(not fatJet1Pt_JESUp_HF_2018_isLoaded){
		if(fatJet1Pt_JESUp_HF_2018_branch != 0) fatJet1Pt_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_HF_2018_branch does not exist!\n");}
		fatJet1Pt_JESUp_HF_2018_isLoaded = true;
}
return fatJet1Pt_JESUp_HF_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_HF_2018(){
	if(not fatJet1PtOverMHH_JESUp_HF_2018_isLoaded){
		if(fatJet1PtOverMHH_JESUp_HF_2018_branch != 0) fatJet1PtOverMHH_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_HF_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_HF_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_HF_2018_;
}
const float &hhtree::fatJet1Pt_JESDown_HF_2018(){
	if(not fatJet1Pt_JESDown_HF_2018_isLoaded){
		if(fatJet1Pt_JESDown_HF_2018_branch != 0) fatJet1Pt_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_HF_2018_branch does not exist!\n");}
		fatJet1Pt_JESDown_HF_2018_isLoaded = true;
}
return fatJet1Pt_JESDown_HF_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_HF_2018(){
	if(not fatJet1PtOverMHH_JESDown_HF_2018_isLoaded){
		if(fatJet1PtOverMHH_JESDown_HF_2018_branch != 0) fatJet1PtOverMHH_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_HF_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_HF_2018_isLoaded = true;
}
	return fatJet1PtOverMHH_JESDown_HF_2018_;
}
const float &hhtree::fatJet1Pt_JESUp_RelSample_2018(){
	if(not fatJet1Pt_JESUp_RelSample_2018_isLoaded){
		if(fatJet1Pt_JESUp_RelSample_2018_branch != 0) fatJet1Pt_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESUp_RelSample_2018_branch does not exist!\n");}
		fatJet1Pt_JESUp_RelSample_2018_isLoaded = true;
}
return fatJet1Pt_JESUp_RelSample_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESUp_RelSample_2018(){
	if(not fatJet1PtOverMHH_JESUp_RelSample_2018_isLoaded){
		if(fatJet1PtOverMHH_JESUp_RelSample_2018_branch != 0) fatJet1PtOverMHH_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESUp_RelSample_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESUp_RelSample_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESUp_RelSample_2018_;
}
const float &hhtree::fatJet1Pt_JESDown_RelSample_2018(){
	if(not fatJet1Pt_JESDown_RelSample_2018_isLoaded){
		if(fatJet1Pt_JESDown_RelSample_2018_branch != 0) fatJet1Pt_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet1Pt_JESDown_RelSample_2018_branch does not exist!\n");}
		fatJet1Pt_JESDown_RelSample_2018_isLoaded = true;
}
return fatJet1Pt_JESDown_RelSample_2018_;
}
const float &hhtree::fatJet1PtOverMHH_JESDown_RelSample_2018(){
	if(not fatJet1PtOverMHH_JESDown_RelSample_2018_isLoaded){
		if(fatJet1PtOverMHH_JESDown_RelSample_2018_branch != 0) fatJet1PtOverMHH_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet1PtOverMHH_JESDown_RelSample_2018_branch does not exist!\n");}
		fatJet1PtOverMHH_JESDown_RelSample_2018_isLoaded = true;
}
return fatJet1PtOverMHH_JESDown_RelSample_2018_;
}

const float &hhtree::fatJet2Pt_JESUp_Abs_2018(){
	if(not fatJet2Pt_JESUp_Abs_2018_isLoaded){
		if(fatJet2Pt_JESUp_Abs_2018_branch != 0) fatJet2Pt_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_Abs_2018_branch does not exist!\n");}
		fatJet2Pt_JESUp_Abs_2018_isLoaded = true;
}
return fatJet2Pt_JESUp_Abs_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_Abs_2018(){
	if(not fatJet2PtOverMHH_JESUp_Abs_2018_isLoaded){
		if(fatJet2PtOverMHH_JESUp_Abs_2018_branch != 0) fatJet2PtOverMHH_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_Abs_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_Abs_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_Abs_2018_;
}
const float &hhtree::fatJet2Pt_JESDowb_Abs_2018(){
	if(not fatJet2Pt_JESDowb_Abs_2018_isLoaded){
		if(fatJet2Pt_JESDowb_Abs_2018_branch != 0) fatJet2Pt_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDowb_Abs_2018_branch does not exist!\n");}
		fatJet2Pt_JESDowb_Abs_2018_isLoaded = true;
}
return fatJet2Pt_JESDowb_Abs_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESDowb_Abs_2018(){
	if(not fatJet2PtOverMHH_JESDowb_Abs_2018_isLoaded){
		if(fatJet2PtOverMHH_JESDowb_Abs_2018_branch != 0) fatJet2PtOverMHH_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDowb_Abs_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDowb_Abs_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESDowb_Abs_2018_;
}
const float &hhtree::fatJet2Pt_JESUp_BBEC1_2018(){
	if(not fatJet2Pt_JESUp_BBEC1_2018_isLoaded){
		if(fatJet2Pt_JESUp_BBEC1_2018_branch != 0) fatJet2Pt_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_BBEC1_2018_branch does not exist!\n");}
		fatJet2Pt_JESUp_BBEC1_2018_isLoaded = true;
}
return fatJet2Pt_JESUp_BBEC1_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_BBEC1_2018(){
	if(not fatJet2PtOverMHH_JESUp_BBEC1_2018_isLoaded){
		if(fatJet2PtOverMHH_JESUp_BBEC1_2018_branch != 0) fatJet2PtOverMHH_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_BBEC1_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_BBEC1_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_BBEC1_2018_;
}
const float &hhtree::fatJet2Pt_JESDown_BBEC1_2018(){
	if(not fatJet2Pt_JESDown_BBEC1_2018_isLoaded){
		if(fatJet2Pt_JESDown_BBEC1_2018_branch != 0) fatJet2Pt_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_BBEC1_2018_branch does not exist!\n");}
		fatJet2Pt_JESDown_BBEC1_2018_isLoaded = true;
}
return fatJet2Pt_JESDown_BBEC1_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_BBEC1_2018(){
	if(not fatJet2PtOverMHH_JESDown_BBEC1_2018_isLoaded){
		if(fatJet2PtOverMHH_JESDown_BBEC1_2018_branch != 0) fatJet2PtOverMHH_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_BBEC1_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_BBEC1_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_BBEC1_2018_;
}

const float &hhtree::fatJet2Pt_JESUp_EC2_2018(){
	if(not fatJet2Pt_JESUp_EC2_2018_isLoaded){
		if(fatJet2Pt_JESUp_EC2_2018_branch != 0) fatJet2Pt_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_EC2_2018_branch does not exist!\n");}
		fatJet2Pt_JESUp_EC2_2018_isLoaded = true;
}
return fatJet2Pt_JESUp_EC2_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_EC2_2018(){
	if(not fatJet2PtOverMHH_JESUp_EC2_2018_isLoaded){
		if(fatJet2PtOverMHH_JESUp_EC2_2018_branch != 0) fatJet2PtOverMHH_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_EC2_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_EC2_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_EC2_2018_;
}
const float &hhtree::fatJet2Pt_JESDown_EC2_2018(){
	if(not fatJet2Pt_JESDown_EC2_2018_isLoaded){
		if(fatJet2Pt_JESDown_EC2_2018_branch != 0) fatJet2Pt_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_EC2_2018_branch does not exist!\n");}
		fatJet2Pt_JESDown_EC2_2018_isLoaded = true;
}
return fatJet2Pt_JESDown_EC2_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_EC2_2018(){
	if(not fatJet2PtOverMHH_JESDown_EC2_2018_isLoaded){
		if(fatJet2PtOverMHH_JESDown_EC2_2018_branch != 0) fatJet2PtOverMHH_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_EC2_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_EC2_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_EC2_2018_;
}
const float &hhtree::fatJet2Pt_JESUp_HF_2018(){
	if(not fatJet2Pt_JESUp_HF_2018_isLoaded){
		if(fatJet2Pt_JESUp_HF_2018_branch != 0) fatJet2Pt_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_HF_2018_branch does not exist!\n");}
		fatJet2Pt_JESUp_HF_2018_isLoaded = true;
}
return fatJet2Pt_JESUp_HF_2018_;
}

const float &hhtree::fatJet2PtOverMHH_JESUp_HF_2018(){
	if(not fatJet2PtOverMHH_JESUp_HF_2018_isLoaded){
		if(fatJet2PtOverMHH_JESUp_HF_2018_branch != 0) fatJet2PtOverMHH_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_HF_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_HF_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_HF_2018_;
}
const float &hhtree::fatJet2Pt_JESDown_HF_2018(){
	if(not fatJet2Pt_JESDown_HF_2018_isLoaded){
		if(fatJet2Pt_JESDown_HF_2018_branch != 0) fatJet2Pt_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_HF_2018_branch does not exist!\n");}
		fatJet2Pt_JESDown_HF_2018_isLoaded = true;
}
return fatJet2Pt_JESDown_HF_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_HF_2018(){
	if(not fatJet2PtOverMHH_JESDown_HF_2018_isLoaded){
		if(fatJet2PtOverMHH_JESDown_HF_2018_branch != 0) fatJet2PtOverMHH_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_HF_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_HF_2018_isLoaded = true;
}
	return fatJet2PtOverMHH_JESDown_HF_2018_;
}
const float &hhtree::fatJet2Pt_JESUp_RelSample_2018(){
	if(not fatJet2Pt_JESUp_RelSample_2018_isLoaded){
		if(fatJet2Pt_JESUp_RelSample_2018_branch != 0) fatJet2Pt_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESUp_RelSample_2018_branch does not exist!\n");}
		fatJet2Pt_JESUp_RelSample_2018_isLoaded = true;
}
return fatJet2Pt_JESUp_RelSample_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESUp_RelSample_2018(){
	if(not fatJet2PtOverMHH_JESUp_RelSample_2018_isLoaded){
		if(fatJet2PtOverMHH_JESUp_RelSample_2018_branch != 0) fatJet2PtOverMHH_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESUp_RelSample_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESUp_RelSample_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESUp_RelSample_2018_;
}
const float &hhtree::fatJet2Pt_JESDown_RelSample_2018(){
	if(not fatJet2Pt_JESDown_RelSample_2018_isLoaded){
		if(fatJet2Pt_JESDown_RelSample_2018_branch != 0) fatJet2Pt_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet2Pt_JESDown_RelSample_2018_branch does not exist!\n");}
		fatJet2Pt_JESDown_RelSample_2018_isLoaded = true;
}
return fatJet2Pt_JESDown_RelSample_2018_;
}
const float &hhtree::fatJet2PtOverMHH_JESDown_RelSample_2018(){
	if(not fatJet2PtOverMHH_JESDown_RelSample_2018_isLoaded){
		if(fatJet2PtOverMHH_JESDown_RelSample_2018_branch != 0) fatJet2PtOverMHH_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch fatJet2PtOverMHH_JESDown_RelSample_2018_branch does not exist!\n");}
		fatJet2PtOverMHH_JESDown_RelSample_2018_isLoaded = true;
}
return fatJet2PtOverMHH_JESDown_RelSample_2018_;
}
const float &hhtree::hh_pt_JESUp_Abs_2018(){
	if(not hh_pt_JESUp_Abs_2018_isLoaded){
		if(hh_pt_JESUp_Abs_2018_branch != 0) hh_pt_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_Abs_2018_branch does not exist!\n");}
		hh_pt_JESUp_Abs_2018_isLoaded = true;
}
return hh_pt_JESUp_Abs_2018_;
}
const float &hhtree::hh_eta_JESUp_Abs_2018(){
	if(not hh_eta_JESUp_Abs_2018_isLoaded){
		if(hh_eta_JESUp_Abs_2018_branch != 0) hh_eta_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_Abs_2018_branch does not exist!\n");}
		hh_eta_JESUp_Abs_2018_isLoaded = true;
}
return hh_eta_JESUp_Abs_2018_;
}
const float &hhtree::hh_mass_JESUp_Abs_2018(){
	if(not hh_mass_JESUp_Abs_2018_isLoaded){
		if(hh_mass_JESUp_Abs_2018_branch != 0) hh_mass_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_Abs_2018_branch does not exist!\n");}
		hh_mass_JESUp_Abs_2018_isLoaded = true;
}
return hh_mass_JESUp_Abs_2018_;
}
const float &hhtree::hh_pt_JESDowb_Abs_2018(){
	if(not hh_pt_JESDowb_Abs_2018_isLoaded){
		if(hh_pt_JESDowb_Abs_2018_branch != 0) hh_pt_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDowb_Abs_2018_branch does not exist!\n");}
		hh_pt_JESDowb_Abs_2018_isLoaded = true;
}
return hh_pt_JESDowb_Abs_2018_;
}
const float &hhtree::hh_eta_JESDowb_Abs_2018(){
	if(not hh_eta_JESDowb_Abs_2018_isLoaded){
		if(hh_eta_JESDowb_Abs_2018_branch != 0) hh_eta_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDowb_Abs_2018_branch does not exist!\n");}
		hh_eta_JESDowb_Abs_2018_isLoaded = true;
}
return hh_eta_JESDowb_Abs_2018_;
}
const float &hhtree::hh_mass_JESDowb_Abs_2018(){
	if(not hh_mass_JESDowb_Abs_2018_isLoaded){
		if(hh_mass_JESDowb_Abs_2018_branch != 0) hh_mass_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDowb_Abs_2018_branch does not exist!\n");}
		hh_mass_JESDowb_Abs_2018_isLoaded = true;
}
return hh_mass_JESDowb_Abs_2018_;
}
const float &hhtree::hh_pt_JESUp_BBEC1_2018(){
	if(not hh_pt_JESUp_BBEC1_2018_isLoaded){
		if(hh_pt_JESUp_BBEC1_2018_branch != 0) hh_pt_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_BBEC1_2018_branch does not exist!\n");}
		hh_pt_JESUp_BBEC1_2018_isLoaded = true;
}
return hh_pt_JESUp_BBEC1_2018_;
}
const float &hhtree::hh_eta_JESUp_BBEC1_2018(){
	if(not hh_eta_JESUp_BBEC1_2018_isLoaded){
		if(hh_eta_JESUp_BBEC1_2018_branch != 0) hh_eta_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_BBEC1_2018_branch does not exist!\n");}
		hh_eta_JESUp_BBEC1_2018_isLoaded = true;
}
return hh_eta_JESUp_BBEC1_2018_;
}
const float &hhtree::hh_mass_JESUp_BBEC1_2018(){
	if(not hh_mass_JESUp_BBEC1_2018_isLoaded){
		if(hh_mass_JESUp_BBEC1_2018_branch != 0) hh_mass_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_BBEC1_2018_branch does not exist!\n");}
		hh_mass_JESUp_BBEC1_2018_isLoaded = true;
}
return hh_mass_JESUp_BBEC1_2018_;
}
const float &hhtree::hh_pt_JESDown_BBEC1_2018(){
	if(not hh_pt_JESDown_BBEC1_2018_isLoaded){
		if(hh_pt_JESDown_BBEC1_2018_branch != 0) hh_pt_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_BBEC1_2018_branch does not exist!\n");}
		hh_pt_JESDown_BBEC1_2018_isLoaded = true;
}
return hh_pt_JESDown_BBEC1_2018_;
}
const float &hhtree::hh_eta_JESDown_BBEC1_2018(){
	if(not hh_eta_JESDown_BBEC1_2018_isLoaded){
		if(hh_eta_JESDown_BBEC1_2018_branch != 0) hh_eta_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_BBEC1_2018_branch does not exist!\n");}
		hh_eta_JESDown_BBEC1_2018_isLoaded = true;
}
return hh_eta_JESDown_BBEC1_2018_;
}
const float &hhtree::hh_mass_JESDown_BBEC1_2018(){
	if(not hh_mass_JESDown_BBEC1_2018_isLoaded){
		if(hh_mass_JESDown_BBEC1_2018_branch != 0) hh_mass_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_BBEC1_2018_branch does not exist!\n");}
		hh_mass_JESDown_BBEC1_2018_isLoaded = true;
}
return hh_mass_JESDown_BBEC1_2018_;
}
const float &hhtree::hh_pt_JESUp_EC2_2018(){
	if(not hh_pt_JESUp_EC2_2018_isLoaded){
		if(hh_pt_JESUp_EC2_2018_branch != 0) hh_pt_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_EC2_2018_branch does not exist!\n");}
		hh_pt_JESUp_EC2_2018_isLoaded = true;
}
return hh_pt_JESUp_EC2_2018_;
}
const float &hhtree::hh_eta_JESUp_EC2_2018(){
	if(not hh_eta_JESUp_EC2_2018_isLoaded){
		if(hh_eta_JESUp_EC2_2018_branch != 0) hh_eta_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_EC2_2018_branch does not exist!\n");}
		hh_eta_JESUp_EC2_2018_isLoaded = true;
}
return hh_eta_JESUp_EC2_2018_;
}
const float &hhtree::hh_mass_JESUp_EC2_2018(){
	if(not hh_mass_JESUp_EC2_2018_isLoaded){
		if(hh_mass_JESUp_EC2_2018_branch != 0) hh_mass_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_EC2_2018_branch does not exist!\n");}
		hh_mass_JESUp_EC2_2018_isLoaded = true;
}
return hh_mass_JESUp_EC2_2018_;
}
const float &hhtree::hh_pt_JESDown_EC2_2018(){
	if(not hh_pt_JESDown_EC2_2018_isLoaded){
		if(hh_pt_JESDown_EC2_2018_branch != 0) hh_pt_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_EC2_2018_branch does not exist!\n");}
		hh_pt_JESDown_EC2_2018_isLoaded = true;
}
return hh_pt_JESDown_EC2_2018_;
}
const float &hhtree::hh_eta_JESDown_EC2_2018(){
	if(not hh_eta_JESDown_EC2_2018_isLoaded){
		if(hh_eta_JESDown_EC2_2018_branch != 0) hh_eta_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_EC2_2018_branch does not exist!\n");}
		hh_eta_JESDown_EC2_2018_isLoaded = true;
}
return hh_eta_JESDown_EC2_2018_;
}
const float &hhtree::hh_mass_JESDown_EC2_2018(){
	if(not hh_mass_JESDown_EC2_2018_isLoaded){
		if(hh_mass_JESDown_EC2_2018_branch != 0) hh_mass_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_EC2_2018_branch does not exist!\n");}
		hh_mass_JESDown_EC2_2018_isLoaded = true;
}
return hh_mass_JESDown_EC2_2018_;
}
const float &hhtree::hh_pt_JESUp_HF_2018(){
	if(not hh_pt_JESUp_HF_2018_isLoaded){
		if(hh_pt_JESUp_HF_2018_branch != 0) hh_pt_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_HF_2018_branch does not exist!\n");}
		hh_pt_JESUp_HF_2018_isLoaded = true;
}
return hh_pt_JESUp_HF_2018_;
}
const float &hhtree::hh_eta_JESUp_HF_2018(){
	if(not hh_eta_JESUp_HF_2018_isLoaded){
		if(hh_eta_JESUp_HF_2018_branch != 0) hh_eta_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_HF_2018_branch does not exist!\n");}
		hh_eta_JESUp_HF_2018_isLoaded = true;
}
return hh_eta_JESUp_HF_2018_;
}
const float &hhtree::hh_mass_JESUp_HF_2018(){
	if(not hh_mass_JESUp_HF_2018_isLoaded){
		if(hh_mass_JESUp_HF_2018_branch != 0) hh_mass_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_HF_2018_branch does not exist!\n");}
		hh_mass_JESUp_HF_2018_isLoaded = true;
}
return hh_mass_JESUp_HF_2018_;
}
const float &hhtree::hh_pt_JESDown_HF_2018(){
	if(not hh_pt_JESDown_HF_2018_isLoaded){
		if(hh_pt_JESDown_HF_2018_branch != 0) hh_pt_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_HF_2018_branch does not exist!\n");}
		hh_pt_JESDown_HF_2018_isLoaded = true;
}
return hh_pt_JESDown_HF_2018_;
}
const float &hhtree::hh_eta_JESDown_HF_2018(){
	if(not hh_eta_JESDown_HF_2018_isLoaded){
		if(hh_eta_JESDown_HF_2018_branch != 0) hh_eta_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_HF_2018_branch does not exist!\n");}
		hh_eta_JESDown_HF_2018_isLoaded = true;
}
return hh_eta_JESDown_HF_2018_;
}
const float &hhtree::hh_mass_JESDown_HF_2018(){
	if(not hh_mass_JESDown_HF_2018_isLoaded){
		if(hh_mass_JESDown_HF_2018_branch != 0) hh_mass_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_HF_2018_branch does not exist!\n");}
		hh_mass_JESDown_HF_2018_isLoaded = true;
}
return hh_mass_JESDown_HF_2018_;
}
const float &hhtree::hh_pt_JESUp_RelSample_2018(){
	if(not hh_pt_JESUp_RelSample_2018_isLoaded){
		if(hh_pt_JESUp_RelSample_2018_branch != 0) hh_pt_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESUp_RelSample_2018_branch does not exist!\n");}
		hh_pt_JESUp_RelSample_2018_isLoaded = true;
}
return hh_pt_JESUp_RelSample_2018_;
}
const float &hhtree::hh_eta_JESUp_RelSample_2018(){
	if(not hh_eta_JESUp_RelSample_2018_isLoaded){
		if(hh_eta_JESUp_RelSample_2018_branch != 0) hh_eta_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESUp_RelSample_2018_branch does not exist!\n");}
		hh_eta_JESUp_RelSample_2018_isLoaded = true;
}
return hh_eta_JESUp_RelSample_2018_;
}
const float &hhtree::hh_mass_JESUp_RelSample_2018(){
	if(not hh_mass_JESUp_RelSample_2018_isLoaded){
		if(hh_mass_JESUp_RelSample_2018_branch != 0) hh_mass_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESUp_RelSample_2018_branch does not exist!\n");}
		hh_mass_JESUp_RelSample_2018_isLoaded = true;
}
return hh_mass_JESUp_RelSample_2018_;
}
const float &hhtree::hh_pt_JESDown_RelSample_2018(){
	if(not hh_pt_JESDown_RelSample_2018_isLoaded){
		if(hh_pt_JESDown_RelSample_2018_branch != 0) hh_pt_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch hh_pt_JESDown_RelSample_2018_branch does not exist!\n");}
		hh_pt_JESDown_RelSample_2018_isLoaded = true;
}
return hh_pt_JESDown_RelSample_2018_;
}
const float &hhtree::hh_eta_JESDown_RelSample_2018(){
	if(not hh_eta_JESDown_RelSample_2018_isLoaded){
		if(hh_eta_JESDown_RelSample_2018_branch != 0) hh_eta_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch hh_eta_JESDown_RelSample_2018_branch does not exist!\n");}
		hh_eta_JESDown_RelSample_2018_isLoaded = true;
}
return hh_eta_JESDown_RelSample_2018_;
}
const float &hhtree::hh_mass_JESDown_RelSample_2018(){
	if(not hh_mass_JESDown_RelSample_2018_isLoaded){
		if(hh_mass_JESDown_RelSample_2018_branch != 0) hh_mass_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch hh_mass_JESDown_RelSample_2018_branch does not exist!\n");}
		hh_mass_JESDown_RelSample_2018_isLoaded = true;
}
return hh_mass_JESDown_RelSample_2018_;
}
const float &hhtree::isVBFtag_JESUp_Abs_2018(){
	if(not isVBFtag_JESUp_Abs_2018_isLoaded){
		if(isVBFtag_JESUp_Abs_2018_branch != 0) isVBFtag_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_Abs_2018_branch does not exist!\n");}
		isVBFtag_JESUp_Abs_2018_isLoaded = true;
}
return isVBFtag_JESUp_Abs_2018_;
}
const float &hhtree::isVBFtag_JESDowb_Abs_2018(){
	if(not isVBFtag_JESDowb_Abs_2018_isLoaded){
		if(isVBFtag_JESDowb_Abs_2018_branch != 0) isVBFtag_JESDowb_Abs_2018_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDowb_Abs_2018_branch does not exist!\n");}
		isVBFtag_JESDowb_Abs_2018_isLoaded = true;
}
return isVBFtag_JESDowb_Abs_2018_;
}
const float &hhtree::isVBFtag_JESUp_EC2_2018(){
	if(not isVBFtag_JESUp_EC2_2018_isLoaded){
		if(isVBFtag_JESUp_EC2_2018_branch != 0) isVBFtag_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_EC2_2018_branch does not exist!\n");}
		isVBFtag_JESUp_EC2_2018_isLoaded = true;
}
return isVBFtag_JESUp_EC2_2018_;
}
const float &hhtree::isVBFtag_JESDown_EC2_2018(){
	if(not isVBFtag_JESDown_EC2_2018_isLoaded){
		if(isVBFtag_JESDown_EC2_2018_branch != 0) isVBFtag_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_EC2_2018_branch does not exist!\n");}
		isVBFtag_JESDown_EC2_2018_isLoaded = true;
}
return isVBFtag_JESDown_EC2_2018_;
}
const float &hhtree::isVBFtag_JESUp_RelSample_2018(){
	if(not isVBFtag_JESUp_RelSample_2018_isLoaded){
		if(isVBFtag_JESUp_RelSample_2018_branch != 0) isVBFtag_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESUp_RelSample_2018_branch does not exist!\n");}
		isVBFtag_JESUp_RelSample_2018_isLoaded = true;
}
return isVBFtag_JESUp_RelSample_2018_;
}
const float &hhtree::isVBFtag_JESDown_RelSample_2018(){
	if(not isVBFtag_JESDown_RelSample_2018_isLoaded){
		if(isVBFtag_JESDown_RelSample_2018_branch != 0) isVBFtag_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch isVBFtag_JESDown_RelSample_2018_branch does not exist!\n");}
		isVBFtag_JESDown_RelSample_2018_isLoaded = true;
}
return isVBFtag_JESDown_RelSample_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_Abs_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_Abs_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_BBEC1_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_BBEC1_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_EC2_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_EC2_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_HF_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_HF_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESUp_RelSample_2018_;
}
const float &hhtree::disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018(){
	if(not disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_isLoaded){
		if(disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_branch != 0) disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_branch->GetEntry(index);
		else {printf("branch disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_branch does not exist!\n");}
		disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_isLoaded = true;
}
return disc_qcd_and_ttbar_Run2_enhanced_v8p2_JESDown_RelSample_2018_;
}


