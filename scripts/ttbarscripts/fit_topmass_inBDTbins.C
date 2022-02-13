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
drawStyle->SetPadRightMargin(0.10);
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

void fit_topmass_inBDTbins() {
    
    const Int_t nbdt = 4;
    
    Double_t x[nbdt]  = {0.00004, 0.00014, 0.0003, 0.0007};
    Double_t yData[nbdt]  = {0, 0, 0, 0};
    Double_t yMC[nbdt]  = {0, 0, 0, 0};

    Double_t ex[nbdt] = {0.00002, 0.00006, 0.0001, 0.0003};
    Double_t eyMC[nbdt] = {0, 0, 0, 0};
    Double_t eyData[nbdt] = {0, 0, 0, 0};

    
    for(int ibdt=1; ibdt< nbdt+1; ibdt++){

    TString ibdtst = Form("%d",ibdt);
    TString name="TTBarCRBDT"+ibdtst+"v8p2";
    
    gStyle->SetOptFit(0111);
    SetPlotStyle();
    
    TF1 *func = new TF1("fit","fitf",-10,10,2);
    func->SetParameters(0,0);
    func->SetParNames("Constant","Slope");
    
    TFile *f = new TFile("../../python/plots/yield_AN_ttbar_reg_cor/combine/"+name+"__fatJet2MassRegressed.root");
    //TFile *f = new TFile("../../python/plots/yield_AN_ttbar_reg_cor/combine/"+name+"__fatJet2MassSD.root");
    
    TH1F *h_data = (TH1F*)f->Get(name+"__fatJet2MassRegressed_data");
    TH1F *h_bkg1 = (TH1F*)f->Get(name+"__fatJet2MassRegressed_bkg_0_stack_1_stack_1");
    TH1F *h_bkg2 = (TH1F*)f->Get(name+"__fatJet2MassRegressed_bkg_1_stack_2_stack_2");
    TH1F *h_bkg3 = (TH1F*)f->Get(name+"__fatJet2MassRegressed_bkg_2_stack_3_stack_3");
    /*TH1F *h_data = (TH1F*)f->Get(name+"__fatJet2MassSD_data");
    TH1F *h_bkg1 = (TH1F*)f->Get(name+"__fatJet2MassSD_bkg_0_stack_1_stack_1");
    TH1F *h_bkg2 = (TH1F*)f->Get(name+"__fatJet2MassSD_bkg_1_stack_2_stack_2");
    TH1F *h_bkg3 = (TH1F*)f->Get(name+"__fatJet2MassSD_bkg_2_stack_3_stack_3");*/
    TF1 *f1 =  new TF1("f1","gaus",130,200);
    TH1F *h_mc = (TH1F*)h_bkg1->Clone("bkg");
    h_mc->Add(h_bkg2);
    h_mc->Add(h_bkg3);
    
    h_mc->Fit("f1","R");
  
    TCanvas *cv = new TCanvas("c","c");
    h_mc->Draw();
    cv->SaveAs(name+"_mcfit.pdf");
    //cv->SaveAs(name+"_mcfit_msd.pdf");
    cout <<"fit data"<<endl;
    h_data->Fit("f1","R");
    TCanvas *cv1 = new TCanvas("c1","c1");
    h_data->Draw();
    cv1->SaveAs(name+"_datafit.pdf");
    //cv1->SaveAs(name+"_datafit_msd.pdf");
    //perform the fit
    cout <<"mc fit result"<<endl;
    TF1 *fit = h_mc->GetFunction("f1");
    cout <<"chi2 "<<fit->GetChisquare()<<" Constant: "<<fit->GetParameter(0)<<"+/-"<<fit->GetParError(0)<< " Slope: "<<fit->GetParameter(1)<<"+/-"<<fit->GetParError(1)<<endl;
    yMC[ibdt-1] = fit->GetParameter(1);
    eyMC[ibdt-1] = fit->GetParError(1);

    cout <<"data fit result"<<endl;
    TF1 *fit_2 = h_data->GetFunction("f1");
    cout <<"chi2 "<<fit_2->GetChisquare()<<" Constant: "<<fit_2->GetParameter(0)<<"+/-"<<fit_2->GetParError(0)<< " Slope: "<<fit_2->GetParameter(1)<<"+/-"<<fit_2->GetParError(1)<<endl;
        
    yData[ibdt-1] = fit_2->GetParameter(1);
    eyData[ibdt-1] = fit_2->GetParError(1);
        
    }
    
    TCanvas *cresult = new TCanvas("cresult","test",200,10,700,500);
    TGraph *grMC = new TGraphErrors(nbdt,x,yMC,ex,eyMC);
    grMC->SetTitle("");
    grMC->GetXaxis()->SetTitle("Event BDT score");
    grMC->GetYaxis()->SetTitle("Fitted top mass [GeV]");
    grMC->GetYaxis()->SetRangeUser(150,190);
    grMC->GetXaxis()->SetRangeUser(0,1);
    grMC->SetMarkerColor(kBlue);
    grMC->SetLineColor(kBlue);
    grMC->SetMarkerStyle(21);
    grMC->SetLineWidth(2);
    grMC->Draw("AP");
    
    TGraph *grData = new TGraphErrors(nbdt,x,yData,ex,eyData);
    grData->SetTitle("");
    grData->SetMarkerColor(kRed);
    grData->SetLineColor(kRed);
    grData->SetMarkerStyle(21);
    grData->SetLineWidth(2);
    grData->Draw("P same");
    
    TLegend *legend = new TLegend(0.6,0.75,0.8,0.9);
    legend->AddEntry(grMC,"MC","l");
    legend->AddEntry(grData,"Data","l");
    legend->Draw();
    
    cresult->SaveAs("mass_vs_BDT_data_vs_mc_v8p2_mreg.pdf");
    //cresult->SaveAs("mass_vs_BDT_data_vs_mc_v8p2_msd.pdf");
    
}
