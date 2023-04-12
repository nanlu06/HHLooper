from __future__ import print_function
import numpy as np
import math
import ROOT as r
import shlex
import sys
import os
from ROOT import TLorentzVector,TGraphAsymmErrors,TMath
from ROOT import Double
from colors  import *

leftMargin   = 0.17
rightMargin  = 0.04
topMargin    = 0.07
bottomMargin = 0.12
r.gStyle.SetOptStat(0)
r.gStyle.SetOptFit(111)
r.gStyle.SetStatX(0.99)
r.gStyle.SetStatY(0.9)
r.gStyle.SetStatW(0.2)
r.gStyle.SetStatH(0.15)
r.gROOT.SetBatch(1)

def get_tgraph(data, subtract=None, denom=None):
    alpha = 1-0.6827
    TGraph_data = TGraphAsymmErrors(data)
    for i in range(TGraph_data.GetN()):
        idenom = 1
        isubtract = 0
        if denom: 
            idenom = denom.GetBinContent(i+1)
        if subtract:
            isubtract = subtract.GetBinContent(i+1)
        N = TGraph_data.GetY()[i]
        if N>0:
            L = Math.gamma_quantile(alpha/2,N,1.)
        U = Math.gamma_quantile_c(alpha/2,N+1,1)
        TGraph_data.SetPointEYlow(i, (N-L)/idenom)
        TGraph_data.SetPointEYhigh(i, (U-N)/idenom)
        TGraph_data.SetPoint(i, TGraph_data.GetX()[i], (N-isubtract)/idenom)
        TGraph_data.SetPointEXlow(i, 0)
        TGraph_data.SetPointEXhigh(i, 0)
    return TGraph_data

def get_sig_bkg_data(vbdt, dirName, bdtbin, pnames_sig, pnames_bkg, region, debug=False):
    dirNameSig = "shapes_fit_s"
    fin = r.TFile("combine-hh/combined_cards_"+vbdt+"/fitDiagnostics"+region+".root", "READ")
    fin_sig = r.TFile("combine-hh/combined_cards_"+vbdt+"/fitDiagnostics"+region+".root", "READ")
    h1_sig = []
    h1_bkg = []
    h1_data = None

    if debug:
        print("debug: dirName", dirName)
    dirpre = fin.GetDirectory(dirName)
    dirthis = dirpre.GetDirectory(bdtbin)

    for idx_bkg in range(len(pnames_bkg)):
        if debug:
            print("debug pnames_bkg", pnames_bkg[idx_bkg], bdtbin)
        if dirthis.GetListOfKeys().Contains(pnames_bkg[idx_bkg]):
            h = fin.Get(dirName+"/"+bdtbin+"/"+pnames_bkg[idx_bkg])
            h.SetDirectory(0)
            h.Scale(10.0)
            h1_bkg.append(h)
            h1_data = h.Clone("h1_data")
            h1_data.SetDirectory(0)
        else:
            h = h1_bkg[0].Clone("h1_"+pnames_bkg[idx_bkg])
            h.SetDirectory(0)
            h.Scale(0.0)
            h1_bkg.append(h)
    for idx_sig in range(len(pnames_sig)):
        if debug:
            print("debug v1 pnames_sig[idx_sig]",pnames_sig[idx_sig])
        if dirthis.GetListOfKeys().Contains(pnames_sig[idx_sig]):
            if debug:
                print("debug: get histogram ",dirNameSig),
                print(bdtbin)
                print(pnames_sig[idx_sig])
            h = fin_sig.Get(dirNameSig+"/"+bdtbin+"/"+pnames_sig[idx_sig])
            h.SetDirectory(0)
            if debug:
                print("debug signal from fin_sig ",h.Integral()) 
            h.Scale(10.0)
            h1_sig.append(h)
        else:
            print("debug signal does not exist")
            h = h1_bkg[0].Clone("h1_"+pnames_sig[idx_sig])
            h.SetDirectory(0)
            h.Scale(0.0)
            h1_sig.append(h)
        
    nBins =  h1_data.GetNbinsX()
                    
    bdtbin_data = bdtbin

    dirthis_data = dirpre.GetDirectory(bdtbin_data)
    if dirthis_data.GetListOfKeys().Contains("data"):
    
        g = fin.Get(dirName+"/"+bdtbin_data+"/data")
        for idx in range(nBins):
            y = g.GetY()[idx]
            h1_data.SetBinContent(idx+1, y*10)
    h1_data.GetXaxis().SetTitle("j_{2} regressed mass (GeV)")

    return h1_sig, h1_bkg, h1_data

def makeplot_single(
    h1_sig=None,
    h1_bkg=None,
    h1_data=None,
    sig_legends_=None,
    bkg_legends_=None,
    sig_colors_=None,
    bkg_colors_=None,
    hist_name_=None,
    sig_scale_=1.0,
    dir_name_="plots",
    output_name_=None,
    extraoptions=None
    ):

    if h1_sig ==  None or h1_bkg == None:
        print("nothing to plot...")
        return
    os.system("mkdir -p "+dir_name_)
    os.system("cp index.php "+dir_name_)
    s_color = [632, 617, 839, 800, 1]
    b_color = [920, 2007, 2005, 2003, 2001, 2011]
    if sig_colors_:
        s_color = sig_colors_
    if bkg_colors_:
        b_color = bkg_colors_
    for idx in range(len(h1_sig)):
        h1_sig[idx].SetLineWidth(3)
        h1_sig[idx].SetLineColor(s_color[idx])
    for idx in range(len(h1_bkg)):
        h1_bkg[idx].SetLineWidth(2)
        h1_bkg[idx].SetLineColor(b_color[idx])
        h1_bkg[idx].SetFillColorAlpha(b_color[idx], 1)
    if h1_data:
        h1_data.SetBinErrorOption(1)
        h1_data.SetLineColor(1)
        h1_data.SetLineWidth(2)
        h1_data.SetMarkerColor(1)
        h1_data.SetMarkerStyle(20)

    myC = r.TCanvas("myC","myC", 600, 600)
    myC.SetTicky(1)
    pad1 = r.TPad("pad1","pad1", 0.05, 0.33,0.95, 0.97)
    pad1.SetBottomMargin(0.027)
    pad1.SetRightMargin( rightMargin )
    pad1.SetLeftMargin( leftMargin )
    pad2 = r.TPad("pad2","pad2", 0.05, 0.04, 0.95, 0.31)
    pad2.SetBottomMargin(0.4)
    pad2.SetTopMargin(0.05)
    pad2.SetRightMargin( rightMargin )
    pad2.SetLeftMargin( leftMargin )

    pad2.Draw()
    pad1.Draw()

    pad1.cd()

    #for idx in range(len(h1_sig)):
    #    print("before signal scaling",h1_sig[idx].Integral())
    #    h1_sig[idx].Scale(1.0)
    #    print("after signal scaling",h1_sig[idx].Integral())
        
    stack = r.THStack("stack", "stack")
    nS = np.zeros(h1_bkg[0].GetNbinsX())
    eS = np.zeros(h1_bkg[0].GetNbinsX())
    #hist_all is used to make the data/mc ratio. remove signal for the moment due to signal is scaled right now
    hist_all = h1_sig[0].Clone("hist_all")
    hist_all.Scale(0.0)
    hist_s = h1_sig[0].Clone("hist_s")
    hist_b = h1_bkg[0].Clone("hist_b")
    for idx in range(len(h1_bkg)):
        stack.Add(h1_bkg[idx])
        for ib in range(h1_bkg[0].GetNbinsX()):
            nS[ib] += h1_bkg[idx].GetBinContent(ib+1)
            eS[ib] = math.sqrt(eS[ib]*eS[ib] + h1_bkg[idx].GetBinError(ib+1)*h1_bkg[idx].GetBinError(ib+1))
        hist_all.Add(h1_bkg[idx]) 
        if idx > 0:
            hist_b.Add(h1_bkg[idx]) 
            
    for idx in range(len(h1_sig)):
        print("ggH signal yield: ", hist_s.Integral())
        if idx > 0:
            hist_temp = h1_sig[idx].Clone(h1_sig[idx].GetName()+"_temp")
            #hist_all.Add(hist_temp)
            hist_s.Add(h1_sig[idx])
        print("all signal yield: ", hist_s.Integral())

    stack.SetTitle("")
    
    maxY = 0.0
    if "stack_signal" in extraoptions and extraoptions["stack_signal"]:
        for idx in range(len(h1_sig)):
            h1_sig[idx].SetFillColorAlpha(s_color[idx], 1)
            stack.Add(h1_sig[idx])
            for ib in range(h1_bkg[0].GetNbinsX()):
                nS[ib] += h1_sig[idx].GetBinContent(ib+1)
                eS[ib] = math.sqrt(eS[ib]*eS[ib] + h1_sig[idx].GetBinError(ib+1)*h1_sig[idx].GetBinError(ib+1))
        if stack.GetMaximum() > maxY:
            maxY = stack.GetMaximum()
        #if "SR" in h.GetTitle(): 
        stack.Draw("hist")
    else:
        stack.Draw("hist")
        if stack.GetMaximum() > maxY:
            maxY = stack.GetMaximum()
        for idx in range(len(h1_sig)):
            if h1_sig[idx].GetMaximum() > maxY:
                maxY = h1_sig[idx].GetMaximum()
            if "SR" in h1_bkg[0].GetTitle():
                #h1_sig[idx].Draw("samehist")
                hist_s.Draw("samehist")

    ##draw  stack total unc on top of total histogram
    box = r.TBox(0,0,1,1,)
    box.SetFillStyle(3002)
    box.SetLineWidth(0)
    box.SetFillColor(r.kBlack)
    for idx in range(h1_bkg[0].GetNbinsX()):
        box.DrawBox(h1_bkg[0].GetBinCenter(idx+1)-0.5*h1_bkg[0].GetBinWidth(idx+1), nS[idx]-eS[idx], h1_bkg[0].GetBinCenter(idx+1)+0.5*h1_bkg[0].GetBinWidth(idx+1), nS[idx]+eS[idx])

    if h1_data:
        if h1_data.GetMaximum() > maxY:
            maxY = h1_data.GetMaximum()+np.sqrt(h1_data.GetMaximum())
        #if not "SR" in h1_data.GetTitle() or "fail" in h1_data.GetTitle():  
        if True:
            #print("debug h1_data.GetName()",h1_data.GetName(), h1_data.GetTitle())           
            TGraph_data = TGraphAsymmErrors(h1_data)
            for i in range(TGraph_data.GetN()):
                #data point
                var_x, var_y = Double(0.), Double(0.)
                TGraph_data.GetPoint(i,var_x,var_y)    
                if np.fabs(var_y) < 1e-5:
                    TGraph_data.SetPoint(i,var_x,-1.0)
                    TGraph_data.SetPointEYlow(i,-1)
                    TGraph_data.SetPointEYhigh(i,-1)
                    #print("zero bins in the data TGraph: bin",i+1)
                else:
                    TGraph_data.SetPoint(i,var_x,var_y)
                    err_low = var_y - (0.5*TMath.ChisquareQuantile(0.1586555,2.*var_y))
                    TGraph_data.SetPointEYlow(i, var_y - (0.5*TMath.ChisquareQuantile(0.1586555,2.*var_y)))
                    TGraph_data.SetPointEYhigh(i, (0.5*TMath.ChisquareQuantile(1.-0.1586555,2.*(var_y+1))) - var_y)
        
            TGraph_data.SetMarkerColor(1)
            TGraph_data.SetMarkerSize(1)
            TGraph_data.SetMarkerStyle(20)
            TGraph_data.Draw("same P")

    stack.GetYaxis().SetTitle("Events")
    stack.GetYaxis().SetTitleOffset(1.05)
    stack.GetYaxis().SetTitleSize(0.08)
    stack.GetYaxis().SetLabelSize(0.06)
    #stack.GetYaxis().CenterTitle()
    stack.GetXaxis().SetLabelSize(0.)
    #stack.GetXaxis().SetLabelOffset(0.013)
    #if "xaxis_range" in extraoptions:
    #    stack.GetXaxis().SetRangeUser(float(extraoptions["xaxis_range"][0]),float(extraoptions["xaxis_range"][1]))

    leg = r.TLegend(0.2, 0.60, 0.9, 0.88)
    leg.SetNColumns(3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.05)
    for idx in range(len(h1_bkg)):
        leg.AddEntry(h1_bkg[idx], bkg_legends_[idx], "F")
    if "SR" in hist_s.GetTitle():
        leg.AddEntry(hist_s, 'HH #times {:1.2}'.format(sig_scale_), "L")

    leg.AddEntry(box, "Total  unc", "F")
    if h1_data:
        leg.AddEntry(h1_data, "Data", "ep")
    leg.Draw()

    pad2.cd()
    pad2.SetGridy(1)
    
    ratio = None
    ratio_Low  = 0.0
    ratio_High  = 4
    
    if h1_data:      
        ratio = TGraphAsymmErrors(h1_data)
        for i in range(ratio.GetN()):
            
            #bkg prediction
            imc = Double(hist_all.GetBinContent(i+1))
            #data point
            var_x, var_y = Double(0.), Double(0.)
            #if not ("SR" in h1_data.GetTitle() and (i>5 and i<9)):
            #if not ("SR" in h1_data.GetTitle()):
            ratio.GetPoint(i,var_x,var_y)    
            if var_y == 0.:
                ratio.SetPoint(i,var_x,-1.0)
                ratio.SetPointEYlow(i,-1)
                ratio.SetPointEYhigh(i,-1)
                continue
            ratio.SetPoint(i,var_x,var_y/imc)
            err_low = (var_y - (0.5*TMath.ChisquareQuantile(0.1586555,2.*var_y)))/imc
            err_high = ((0.5*TMath.ChisquareQuantile(1.-0.1586555,2.*(var_y+1))) - var_y)/imc
            ratio.SetPointEYlow(i, err_low)
            ratio.SetPointEYhigh(i, err_high)
        
        ratio.SetMarkerColor(1)
        ratio.SetMarkerSize(1)
        ratio.SetMarkerStyle(20)
        ratio.GetXaxis().SetTitle("j_{2} regressed mass [GeV]")
        #myC.Update()
        
        if "ratio_range" in extraoptions:
            ratio_Low = extraoptions["ratio_range"][0]
            ratio_High = extraoptions["ratio_range"][1]
        ratio.GetYaxis().SetTitle("data/mc")
        ratio.GetYaxis().SetRangeUser(ratio_Low, ratio_High)
        ratio.GetXaxis().SetRangeUser(50, 220)
        ratio.SetTitle("")
        ratio.Draw("same AP")
        pad2.Update()
        
        print(ratio.GetTitle(),ratio.GetName(),"debug")
    else:
        ratio = h1_sig[0].Clone("ratio")
        ratio_High = 0.0
        for ibin in range(1,ratio.GetNbinsX()+1):
            s = hist_s.GetBinContent(ibin) 
            b = hist_b.GetBinContent(ibin)
            L = 0.0
            if b > 0.0:
                L = s/math.sqrt(b)
                if L > ratio_High:
                    ratio_High = L
            ratio.SetBinContent(ibin, L)
        if ratio_High > 1.0:
            ratio_High = 1.0
        ratio.GetYaxis().SetRangeUser(ratio_Low, ratio_High*1.2)
        ratio.GetYaxis().SetTitle("S/#sqrt{B}")
        ratio.Draw("samehist")
    ratio.SetLineColor(1)
    ratio.SetLineWidth(2)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerColor(1)
    ratio.SetFillColorAlpha(1, 0)
    ratio.GetXaxis().SetTitleOffset(0.94)
    ratio.GetXaxis().SetTitleSize(0.18)
    ratio.GetXaxis().SetLabelSize(0.12)
    ratio.GetXaxis().SetLabelOffset(0.013)
    ratio.GetYaxis().SetTitleOffset(0.40)
    ratio.GetYaxis().SetTitleSize(0.17)
    ratio.GetYaxis().SetLabelSize(0.13)
    ratio.GetYaxis().SetTickLength(0.01)
    ratio.GetYaxis().SetNdivisions(505)
    #if "xaxis_range" in extraoptions:
    #    ratio.GetXaxis().SetRangeUser(float(extraoptions["xaxis_range"][0]),float(extraoptions["xaxis_range"][1]))

    #draw  stack total unc on the ratio plot to present the background uncertainty
    box_ratio = r.TBox(0,0,1,1,)
    box_ratio.SetFillStyle(3002)
    box_ratio.SetLineWidth(0)
    box_ratio.SetFillColor(r.kBlack)
    for idx in range(h1_bkg[0].GetNbinsX()):
        if np.fabs(nS[idx])> 1e-06: 
            box_ratio.DrawBox(h1_bkg[0].GetBinCenter(idx+1)-0.5*h1_bkg[0].GetBinWidth(idx+1), (nS[idx]-eS[idx])/nS[idx], h1_bkg[0].GetBinCenter(idx+1)+0.5*h1_bkg[0].GetBinWidth(idx+1), (nS[idx]+eS[idx])/nS[idx])
        else:
            print("blinded Higgs peak region") 
    
    if "xaxis_label" in extraoptions and extraoptions["xaxis_label"] != None:
        x_title = extraoptions["xaxis_label"]
        ratio.GetXaxis().SetTitle(x_title)
    ratio.GetYaxis().CenterTitle()

    ##########draw CMS preliminary
    pad1.cd()
    tex1 = r.TLatex(leftMargin, 0.91, "CMS")
    tex1.SetNDC()
    tex1.SetTextFont(61)
    tex1.SetTextSize(0.070)
    tex1.SetLineWidth(2)
    tex1.Draw()
    tex2 = r.TLatex(leftMargin+0.12,0.912,"Internal")
    tex2.SetNDC()
    tex2.SetTextFont(52)
    tex2.SetTextSize(0.055)
    tex2.SetLineWidth(2)
    tex2.Draw()

    lumi_value = 138
    if "lumi_value" in extraoptions:
        lumi_value = extraoptions["lumi_value"]
    tex3 = r.TLatex(0.72,0.912,"%d"%lumi_value+" fb^{-1} (13 TeV)")
    tex3.SetNDC()
    tex3.SetTextFont(42)
    tex3.SetTextSize(0.055)
    tex3.SetLineWidth(2)
    tex3.Draw()
    outFile = dir_name_
    if output_name_:
        outFile = outFile + "/" +output_name_
    else:
        outFile = outFile + "/" + hist_name_

    #print("maxY = "+str(maxY))
    stack.SetMaximum(maxY*1.7)

    #print everything into txt file
    text_file = open(outFile+"_linY.txt", "w")
    text_file.write("bin    |   x    ")
    for idx in range(len(h1_bkg)):
        text_file.write(" | %21s"%bkg_legends_[idx])
    text_file.write(" | %21s"%("total B"))
    for idx in range(len(sig_legends_)):
        text_file.write(" | %25s"%sig_legends_[idx])
    if h1_data:
        text_file.write(" | data | data/mc")
    text_file.write("\n-------------")
    for idx in range(24*(len(h1_bkg) + 1)+ 29*len(sig_legends_)):
        text_file.write("-")
    if h1_data:
        text_file.write("-------")
    text_file.write("\n")
    for ibin in range(0,h1_sig[0].GetNbinsX()+1):
        text_file.write("%3d"%ibin+"   ")
        text_file.write(" | %6.3f"%h1_data.GetBinCenter(ibin)+" ")
        for idx in range(len(h1_bkg)):
            text_file.write(" | %7.3f "%h1_bkg[idx].GetBinContent(ibin)+"$\\pm$"+ " %7.3f"%h1_bkg[idx].GetBinError(ibin))
        text_file.write(" | %7.3f "%hist_b.GetBinContent(ibin)+"$\\pm$"+ " %7.3f"%hist_b.GetBinError(ibin))
        for idx in range(len(sig_legends_)):
            text_file.write(" | %9.3f "%h1_sig[idx].GetBinContent(ibin)+"$\\pm$"+ " %9.3f"%h1_sig[idx].GetBinError(ibin))
        if h1_data:
            text_file.write(" | %d"%h1_data.GetBinContent(ibin) +  " | %7.3f "%h1_data.GetBinContent(ibin) +"$\\pm$"+ " %7.3f"%h1_data.GetBinError(ibin))
        text_file.write("\n\n")
        
    #print yield table for AN
    text_file.write("print yield table for AN\n")
    bkg_all = 0
    bkg_all_errsq = 0
    for idx in range(len(h1_bkg)):
        bkg_tmp = h1_bkg[idx].GetBinContent(7)+h1_bkg[idx].GetBinContent(8)+h1_bkg[idx].GetBinContent(9)
        bkg_errsq_tmp = h1_bkg[idx].GetBinError(7)*h1_bkg[idx].GetBinError(7)+h1_bkg[idx].GetBinError(8)*h1_bkg[idx].GetBinError(8)+h1_bkg[idx].GetBinError(9)*h1_bkg[idx].GetBinError(9)
        bkg_all += bkg_tmp
        bkg_all_errsq += bkg_errsq_tmp
        text_file.write("%s"%(bkg_legends_[idx])+"& %7.2f"%(bkg_tmp)+"$\\pm$"+ "%7.2f"%np.sqrt(bkg_errsq_tmp)+"\n")
    text_file.write("total background & %7.2f"%(bkg_all)+"$\\pm$"+ "%7.2f"%np.sqrt(bkg_all_errsq)+"\n")
    
    text_file.write("\ggHH SM ($\kapl=1$) & %7.2f"%((h1_sig[0].GetBinContent(7)+h1_sig[0].GetBinContent(8)+h1_sig[0].GetBinContent(9))/sig_scale_)+"$\\pm$"+ "%7.1f"%(sig_scale_*np.sqrt(h1_sig[0].GetBinError(7)*h1_sig[0].GetBinError(7)+h1_sig[0].GetBinError(8)*h1_sig[0].GetBinError(8)+h1_sig[0].GetBinError(9)*h1_sig[0].GetBinError(9)))+"\n")
    text_file.write("\VBFHH SM ($\kapl=1$) & %7.2f"%((h1_sig[1].GetBinContent(7)+h1_sig[1].GetBinContent(8)+h1_sig[1].GetBinContent(9))/sig_scale_)+"$\\pm$"+ "%7.1f"%(sig_scale_*np.sqrt(h1_sig[1].GetBinError(7)*h1_sig[1].GetBinError(7)+h1_sig[1].GetBinError(8)*h1_sig[1].GetBinError(8)+h1_sig[1].GetBinError(9)*h1_sig[1].GetBinError(9)))+"\n")
    
    text_file.write("HH bin 8 value %s"%h1_sig[0].GetBinContent(8)+"\n")
    text_file.write("HH bin 9 value %s"%h1_sig[0].GetBinContent(9)+"\n")
    text_file.write("HH bin 7 value %s"%h1_sig[0].GetBinContent(7)+"\n")

    text_file.write("HH bin 8 error %s"%h1_sig[0].GetBinError(8)+"\n")
    text_file.write("HH bin 9 error %s"%h1_sig[0].GetBinError(9)+"\n")
    text_file.write("HH bin 7 error %s"%h1_sig[0].GetBinError(7)+"\n")
    
    text_file.write("total & %7.2f"%(bkg_all+(h1_sig[0].GetBinContent(7)+h1_sig[0].GetBinContent(8)+h1_sig[0].GetBinContent(9)+h1_sig[1].GetBinContent(7)+h1_sig[1].GetBinContent(8)+h1_sig[1].GetBinContent(9))/sig_scale_)+"$\\pm$"+ "%7.2f"%(np.sqrt((h1_sig[0].GetBinError(7)*h1_sig[0].GetBinError(7)+h1_sig[0].GetBinError(8)*h1_sig[0].GetBinError(8)+h1_sig[0].GetBinError(9)*h1_sig[0].GetBinError(9))/(sig_scale_*sig_scale_)+(h1_sig[1].GetBinError(7)*h1_sig[1].GetBinError(7)+h1_sig[1].GetBinError(8)*h1_sig[1].GetBinError(8)+h1_sig[1].GetBinError(9)*h1_sig[1].GetBinError(9))/(sig_scale_*sig_scale_)+bkg_all_errsq))+"\n")
    
    text_file.close()
    os.system("cp "+outFile+"_linY.txt "+outFile+"_logY.txt")

    pad1.RedrawAxis()
    myC.SaveAs(outFile+"_linY.png")
    myC.SaveAs(outFile+"_linY.pdf")
    myC.SaveAs(outFile+"_linY.C")
    pad1.cd()
    stack.SetMaximum(maxY*100.0)
    stack.SetMinimum(0.5)
    pad1.SetLogy()
    pad1.RedrawAxis()
    myC.SaveAs(outFile+"_logY.png")
    myC.SaveAs(outFile+"_logY.pdf")
    myC.SaveAs(outFile+"_logY.C")
    #save histogram and ratio to root file
    outFile_root = r.TFile(outFile+".root", "recreate")
    outFile_root.cd()
    for idx in range(len(h1_bkg)):
        h1_bkg[idx].Write()
    for idx in range(len(sig_legends_)):
        h1_sig[idx].Write()
    if  h1_data:
        h1_data.Write()
        ratio.Write()
    #outFile_root.Write()
    outFile_root.Close()

def main(vbdt, HH_limit):
    
    dirNameSig = "shapes_fit_s"
    
    debug = True;

    for dirName in ["shapes_fit_s",  "shapes_fit_b"]:
    
        pnames_sig = ["ggHH_kl_1_kt_1_hbbhbb", "qqHH_CV_1_C2V_1_kl_1_hbbhbb"]
        pnames_bkg = ["ttH_hbb", "VH_hbb", "bbbb_boosted_ggf_others", "ttbar", "bbbb_boosted_ggf_qcd_datadriven"]
        bkg_legends = ["t#bar{t}H", "VH", "V+jets,VV", "t#bar{t}+jets", "QCD+ggH+VBFH"]
        sig_legends = ["ggHH", "VBFHH"]
        pname_data = "data"
        sig_colors = [2, 839, 800, 1, 632]
        bkg_colors = [2003, 2011, 2001, 2005, 2007, 800, 839]
        
        #color scheme for paper draft
        #sig_colors = [2, 4, 800, 1, 632]
        #bkg_colors = [616, 2011, 2001, 601, 797, 2007, 839]

        region = "SBplusfail"
        #the first three are the sideband region, the next three is the full AK8 jet 2 mass region, the last one is the common fail region
        bdtbins = ["SRBin1", "SRBin2", "SRBin3","fitfail"]
        for bdtbin in bdtbins:

            h1_sig, h1_bkg, h1_data = get_sig_bkg_data(vbdt, dirName, bdtbin, pnames_sig, pnames_bkg, region, debug=debug)

            makeplot_single(
                h1_sig=h1_sig,
                h1_bkg=h1_bkg,
                h1_data=h1_data,
                sig_legends_=sig_legends,
                bkg_legends_=bkg_legends,
                sig_colors_=sig_colors,
                bkg_colors_=bkg_colors,
                sig_scale_=HH_limit,
                output_name_=dirName+"_"+bdtbin+"_BDTv8p2",
                dir_name_= "plots/postfitplots/output_all_"+vbdt+"/",
                extraoptions={"stack_signal": False}
                )

if __name__ == "__main__":
    vbdt = "v8p2yield_AN_sr_sys_230222v1"
    #"v8p2yield_AN_sr_sys_0830_fix2017trigSF0908_SDv1"
    #"v8p2yield_AN_sr_sys_0830_fix2017trigSF0908"
    HH_limit = 3.8;
    main(vbdt,HH_limit)
