#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TObject.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>
#include <TPaveText.h>

TStyle* DrawStyle()
{

TStyle *drawStyle = new TStyle("Plot","Plot style");
Int_t icol=0;
drawStyle->SetFrameBorderMode(icol);
drawStyle->SetFrameFillColor(icol);
drawStyle->SetCanvasBorderMode(icol);
drawStyle->SetCanvasColor(icol);
drawStyle->SetPadBorderMode(icol);
drawStyle->SetPadColor(icol);
drawStyle->SetStatColor(icol);

drawStyle->SetPadTopMargin(0.05);
drawStyle->SetPadRightMargin(0.05);
drawStyle->SetPadBottomMargin(0.16);
drawStyle->SetPadLeftMargin(0.16);

drawStyle->SetTitleXOffset(1.4);
drawStyle->SetTitleYOffset(1.4);

Int_t font=42;
Double_t tsize=0.05;
drawStyle->SetTextFont(font);

drawStyle->SetTextSize(tsize);
drawStyle->SetLabelFont(font,"x");
drawStyle->SetTitleFont(font,"x");
drawStyle->SetLabelFont(font,"y");
drawStyle->SetTitleFont(font,"y");
drawStyle->SetLabelFont(font,"z");
drawStyle->SetTitleFont(font,"z");

drawStyle->SetLabelSize(tsize,"x");
drawStyle->SetTitleSize(tsize,"x");
drawStyle->SetLabelSize(tsize,"y");
drawStyle->SetTitleSize(tsize,"y");
drawStyle->SetLabelSize(tsize,"z");
drawStyle->SetTitleSize(tsize,"z");

drawStyle->SetMarkerStyle(20);
drawStyle->SetMarkerSize(1.);
drawStyle->SetHistLineWidth(2.);
drawStyle->SetLineStyleString(2,"[12 12]");

drawStyle->SetEndErrorSize(0.);

drawStyle->SetOptTitle(0);
drawStyle->SetOptStat(0);
drawStyle->SetOptFit(1);

drawStyle->SetPadTickX(1);
drawStyle->SetPadTickY(1);

return drawStyle;

}


void SetPlotStyle()
{
    std::cout << "\nApplying plot style settings...\n" << std::endl ;
    TStyle* drawStyle = DrawStyle();
    gROOT->SetStyle("Plot");
    gROOT->ForceStyle();
}

void PNetcorrection() {
    
    gStyle->SetOptFit(0111);
    SetPlotStyle();
    
    TString year="2018";

    TFile *f = new TFile("../../hists/yield_AN_1Lttbar_0329pileupweight/"+year+"/TTBarLepJetCR__fatJet1PNetXbb_Bin2.root");
    TH1F *h_data = (TH1F*)f->Get("TTBarLepJetCR__fatJet1PNetXbb_Bin2_data");
    TH1F *h_bkg1 = (TH1F*)f->Get("TTBarLepJetCR__fatJet1PNetXbb_Bin2_bkg_1_stack_2_stack_2");
    TH1F *h_bkg2 = (TH1F*)f->Get("TTBarLepJetCR__fatJet1PNetXbb_Bin2_bkg_2_stack_3_stack_3");
    TH1F *h_sig2 = (TH1F*)f->Get("TTBarLepJetCR__fatJet1PNetXbb_Bin2_bkg_0_stack_1_stack_1");
    TH1F *h_sig1 = (TH1F*)f->Get("TTBarLepJetCR__fatJet1PNetXbb_Bin2_sig_0_stack_4_stack_4");
    
    TH1F *h_data_mbkg = (TH1F*)h_data->Clone("data_mbkg");
    h_data_mbkg->Add(h_bkg1,-1);
    h_data_mbkg->Add(h_bkg2,-1);
    h_sig1->Add(h_sig2);
    h_data_mbkg->Divide(h_sig1);
    TH1F *hpx = (TH1F*)f->Get("ratio_data_over_mc");
    TCanvas *cv = new TCanvas("c","c");
    hpx->Draw();
    h_data_mbkg->SetLineColor(kRed);
    h_data_mbkg->Draw("same");

    cout<< "result to be used Data-bkg/sig:"<<endl;
    h_data_mbkg->Print("all");

    cout<< "check Data/MC:"<<endl;
    hpx->Print("all");
    
    cv->SaveAs("ratio_comparison_PNet_ttbar.pdf");
}
