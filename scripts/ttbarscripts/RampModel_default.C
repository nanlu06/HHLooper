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
#include <TMatrixD.h>

const int npars = 2;
TMatrixD cov_data_mbkg_mc(npars, npars); 

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

Double_t fitfrf(Double_t *x,Double_t *par) {
    Double_t fitval = par[0];
    if(x[0]<300){
        fitval = par[0] + par[1]*x[0] ;
    }
    else{
        fitval = par[0] + 300.*par[1] - 300.*par[2] + par[2]*x[0] ;
    }
    return fitval;
}

Double_t fitf_up(Double_t *x,Double_t *par) {
    
    TMatrixD a(1,2);
    TArrayD adata(2);

    adata[0] = 1.;
    adata[1] = x[0]-300.;
    a.SetMatrixArray(adata.GetArray());
    
    TMatrixD aT(2,1);
    aT.Transpose(a);

    TMatrixD c(1,2);
    TArrayD cdata(2);

    cdata[0] = 1.;
    cdata[1] = 0.;
    c.SetMatrixArray(cdata.GetArray());

    TMatrixD cT(2,1);
    cT.Transpose(c);

    //cout <<"inside fit: cov_data_mbkg_mc"<<endl;
    //cov_data_mbkg_mc.Print();
    //cout <<" outcome: "<<endl;

    TMatrixD tmp1(1,2),err1(1,1);
    tmp1.Mult(c,cov_data_mbkg_mc);
    err1.Mult(tmp1,cT);
    //err1.Print();
    //cout << err1.Determinant() <<endl;
    
    TMatrixD tmp2(1,2),err2(1,1);
    tmp2.Mult(a,cov_data_mbkg_mc);
    err2.Mult(tmp2,aT);
    //err2.Print();

    //cov_data_mbkg_mc
    Double_t fitval = par[0] + sqrt(err1.Determinant());
    if(x[0]<300){
        fitval = par[0]-300.0*par[1] + par[1]*x[0] + sqrt(err2.Determinant());
    }
    return fitval;
}

Double_t fitf_dn(Double_t *x,Double_t *par) {
    
    TMatrixD a(1,2);
    TArrayD adata(2);

    adata[0] = 1.;
    adata[1] = x[0]-300.;
    a.SetMatrixArray(adata.GetArray());
    
    TMatrixD aT(2,1);
    aT.Transpose(a);

    TMatrixD c(1,2);
    TArrayD cdata(2);

    cdata[0] = 1.;
    cdata[1] = 0.;
    c.SetMatrixArray(cdata.GetArray());

    TMatrixD cT(2,1);
    cT.Transpose(c);

    //cout <<"inside fit: cov_data_mbkg_mc"<<endl;
    //cov_data_mbkg_mc.Print();
    //cout <<" outcome: "<<endl;

    TMatrixD tmp1(1,2),err1(1,1);
    tmp1.Mult(c,cov_data_mbkg_mc);
    err1.Mult(tmp1,cT);
    //err1.Print();
    //cout << err1.Determinant() <<endl;
    
    TMatrixD tmp2(1,2),err2(1,1);
    tmp2.Mult(a,cov_data_mbkg_mc);
    err2.Mult(tmp2,aT);
    //err2.Print();

    //cov_data_mbkg_mc
    Double_t fitval = par[0] - sqrt(err1.Determinant());
    if(x[0]<300){
        fitval = par[0]-300.0*par[1] + par[1]*x[0] - sqrt(err2.Determinant());
    }
    return fitval;
}

Double_t fitf(Double_t *x,Double_t *par) {
    Double_t fitval = par[0];
    if(x[0]<300){
        fitval = par[0]-300.0*par[1] + par[1]*x[0];
    }
    return fitval;
}

void RampModel_default() {
    TString year="2018";
    
    gStyle->SetOptFit(0111);
    SetPlotStyle();
    
    float xmin = 0;
    float xmax = 1000.;
    //TF1 *func = new TF1("fit",fitfrf,-10,10,npars);
    TF1 *func = new TF1("fit",fitf,xmin,xmax,npars);

    //TF1 *func = new TF1("fit",fitf,0,1000);
    //func->SetParameters(0,0);
    //func->SetParNames("Constant","Slope");
    
    if(npars==3){
        func->SetParameters(0,0,0);
        func->SetParNames("Constant1","Slope1","Slope2");
    }
    else{
        func->SetParameters(0,0);
        func->SetParNames("Constant","Slope");
    }
    
    TFile *f = new TFile("yield_AN_ttbar_0329pileupweight/"+year+"/TTBarCR__hh_pt.root");
    
    TH1F *h_data = (TH1F*)f->Get("TTBarCR__hh_pt_data");
    TH1F *h_bkg1 = (TH1F*)f->Get("TTBarCR__hh_pt_bkg_0_stack_1_stack_1");
    TH1F *h_bkg2 = (TH1F*)f->Get("TTBarCR__hh_pt_bkg_1_stack_2_stack_2");
    TH1F *h_sig = (TH1F*)f->Get("TTBarCR__hh_pt_bkg_2_stack_3_stack_3");
    
    //ratio of data-mbkg/MC
    TH1F *h_data_mbkg = (TH1F*)h_data->Clone("data_mbkg");
    h_data_mbkg->Add(h_bkg1,-1);
    h_data_mbkg->Add(h_bkg2,-1);
    h_data_mbkg->Divide(h_sig);
    TFitResultPtr r_data_mbkg_mc = h_data_mbkg->Fit("fit","ES");
    cout <<"Fit result of data-mbkg/MC -----"<<endl;
    r_data_mbkg_mc->Print();
    cov_data_mbkg_mc = r_data_mbkg_mc->GetCovarianceMatrix();
    cov_data_mbkg_mc.Print();
    
    //ratio of data/MC
    TH1F *hpx = (TH1F*)f->Get("ratio_data_over_mc");
    TFitResultPtr r_data_mc = hpx->Fit("fit","ES");
    cout <<"Fit result of data/MC -----"<<endl;
    r_data_mc->Print();
    //TMatrixD cov_data_mc = r_data_mc->GetCovarianceMatrix();
    //cov_data_mc.Print();
    
    //get parameters from the fit
    TF1 *fit = hpx->GetFunction("fit");
    cout <<"fit parameters for data/ttbar MC chi2 "<<fit->GetChisquare()<<endl;
    for(int k=0; k<npars; k++){
        cout <<"par"<<k<<" : "<<fit->GetParameter(k)<<"+/-"<<fit->GetParError(k)<<endl;
    }
    
    TCanvas *cv = new TCanvas("c","c");
    hpx->Draw();
    fit->SetLineColor(kRed);
    fit->Draw("same");
    
    cv->SaveAs(year+"_fit_piecewise_data_MC.png");

    //get parameters from the fit
    TF1 *fit_2 = h_data_mbkg->GetFunction("fit");
    cout <<"fit parameters for data-mbkg/ttbar MC chi2 "<<fit_2->GetChisquare()<<endl;
    for(int k=0; k<npars; k++){
        cout <<"par"<<k<<" : "<<fit_2->GetParameter(k)<<"+/-"<<fit_2->GetParError(k)<<endl;
    }
    
    TF1 *fit_2_up = new TF1("fit_2_up",fitf_up,xmin,xmax,npars);
    fit_2_up->SetLineColor(kBlue); 
    fit_2_up->SetLineStyle(2); 
    fit_2_up->SetParameters(fit_2->GetParameter(0), fit_2->GetParameter(1));
    for(int k=0; k<npars; k++){
        cout <<"par for fit_2_up "<<k<<" : "<<fit_2_up->GetParameter(k)<<"+/-"<<fit_2_up->GetParError(k)<<endl;
    }

    TF1 *fit_2_dn = new TF1("fit_2_down",fitf_dn,xmin,xmax,npars);
    fit_2_dn->SetLineColor(kBlue);
    fit_2_dn->SetLineStyle(2); 
    fit_2_dn->SetParameters(fit_2->GetParameter(0), fit_2->GetParameter(1));
    for(int k=0; k<npars; k++){
        cout <<"par for fit_2_dn "<<k<<" : "<<fit_2_dn->GetParameter(k)<<"+/-"<<fit_2_dn->GetParError(k)<<endl;
    }

    TCanvas *cv2 = new TCanvas("c2","c2");
    h_data_mbkg->GetYaxis()->SetRangeUser(0,2);
    h_data_mbkg->Draw();
    fit_2->SetLineColor(kBlue);
    fit_2->Draw("same");
    fit_2_up->Draw("same");
    fit_2_dn->Draw("same");
    cv2->SaveAs(year+"_fit_piecewise_data_mMC.png");
    
}
