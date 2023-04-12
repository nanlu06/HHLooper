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

Double_t fitf(Double_t *x,Double_t *par) {
    Double_t fitval = par[0];
    if(x[0]<300){
        fitval = par[0]-300.0*par[1] + par[1]*x[0];
    }
    return fitval;
}

void RampModel() {
    TString year="2017";
    
    gStyle->SetOptFit(0111);
    SetPlotStyle();
    
    TF1 *func = new TF1("fit",fitf,-10,10,2);
    //TF1 *func = new TF1("fit",fitf,0,1000);
    func->SetParameters(0,0);
    func->SetParNames("Constant","Slope");
    

    TFile *f = new TFile("../../python/plots/yield_AN_ttbar_reg//"+year+"/TTBarCR__hh_pt.root");
    
    TH1F *h_data = (TH1F*)f->Get("TTBarCR__hh_pt_data");
    TH1F *h_bkg1 = (TH1F*)f->Get("TTBarCR__hh_pt_bkg_0_stack_1_stack_1");
    TH1F *h_bkg2 = (TH1F*)f->Get("TTBarCR__hh_pt_bkg_1_stack_2_stack_2");
    TH1F *h_sig = (TH1F*)f->Get("TTBarCR__hh_pt_bkg_2_stack_3_stack_3");
    
    TH1F *h_data_mbkg = (TH1F*)h_data->Clone("data_mbkg");
    h_data_mbkg->Add(h_bkg1,-1);
    h_data_mbkg->Add(h_bkg2,-1);
    h_data_mbkg->Divide(h_sig);
    h_data_mbkg->Fit("fit","E");
    
    //ratio of data/MC
    TH1F *hpx = (TH1F*)f->Get("ratio_data_over_mc");
    TCanvas *cv = new TCanvas("c","c");
    hpx->Fit("fit","E");
    hpx->Draw();
    
    cout<< "Data-bkg/sig:"<<endl;
    h_data_mbkg->Print("all");
    cout<< "Data/MC:"<<endl;
    hpx->Print("all");
    
    //perform the fit
    TF1 *fit = hpx->GetFunction("fit");
    cout <<"chi2 "<<fit->GetChisquare()<<" Constant: "<<fit->GetParameter(0)<<"+/-"<<fit->GetParError(0)<< " Slope: "<<fit->GetParameter(1)<<"+/-"<<fit->GetParError(1)<<endl;
    
    TF1 *fit_2 = h_data_mbkg->GetFunction("fit");
    cout <<"chi2 "<<fit_2->GetChisquare()<<" Constant: "<<fit_2->GetParameter(0)<<"+/-"<<fit_2->GetParError(0)<< " Slope: "<<fit_2->GetParameter(1)<<"+/-"<<fit_2->GetParError(1)<<endl;
    
    float par0 = fit->GetParameter(0);
    float err0 = fit->GetParError(0);
    float par1 = fit->GetParameter(1);
    float err1 = fit->GetParError(1);
    
    fit->SetLineColor(kRed);
    fit->Draw("same");
    
    cv->SaveAs(year+"_fit_piecewise_data_MC.png");

    TCanvas *cv2 = new TCanvas("c2","c2");
    h_data_mbkg->GetYaxis()->SetRangeUser(0,2);
    h_data_mbkg->Draw();
    fit_2->SetLineColor(kBlue);
    fit_2->Draw("same");
    cv2->SaveAs(year+"_fit_piecewise_data_mMC.png");
    
}
