from __future__ import print_function
import numpy as np
import math
import ROOT as r
import shlex
import sys
import os

leftMargin   = 0.13
rightMargin  = 0.05
topMargin    = 0.07
bottomMargin = 0.12
r.gStyle.SetOptStat(0)
r.gStyle.SetOptFit(111)
r.gStyle.SetStatX(0.99)
r.gStyle.SetStatY(0.9)
r.gStyle.SetStatW(0.2)
r.gStyle.SetStatH(0.15)
r.gROOT.SetBatch(1)

def rebin(hists, nbin):
    for hist in hists:
        if not hist: continue
        currnbin = hist.GetNbinsX()
        fac = currnbin / nbin
        if float(fac).is_integer() and fac > 0:
            hist.Rebin(int(fac))
def remove_underflow(hists):
    def func(hist):
        hist.SetBinContent(0, 0)
        hist.SetBinError(0, 0)
    if isinstance(hists, list):
        for hist in hists:
            func(hist)
    else:
        func(hists)
    return hists
def remove_overflow(hists):
    def func(hist):
        hist.SetBinContent(hist.GetNbinsX()+1, 0)
        hist.SetBinError(hist.GetNbinsX()+1, 0)
    if isinstance(hists, list):
        for hist in hists:
            func(hist)
    else:
        func(hists)
    return hists

def add_underflow(hists):
    def func(hist):
        hist.SetBinContent(1, hist.GetBinContent(0)+hist.GetBinContent(1))
        hist.SetBinError(1, math.sqrt(hist.GetBinError(0)*hist.GetBinError(0)+hist.GetBinError(1)*hist.GetBinError(1)))
    if isinstance(hists, list):
        for hist in hists:
            func(hist)
    else:
        func(hists)
    return hists
def add_overflow(hists):
    def func(hist):
        nx=hist.GetNbinsX()
        hist.SetBinContent(nx, hist.GetBinContent(nx)+hist.GetBinContent(nx+1))
        hist.SetBinError(nx, math.sqrt(hist.GetBinError(nx)*hist.GetBinError(nx)+hist.GetBinError(nx+1)*hist.GetBinError(nx+1)))
    if isinstance(hists, list):
        for hist in hists:
            func(hist)
    else:
        func(hists)
    return hists

def blind_data(hists):
    def func(hist):
        nx=hist.FindBin(125.0)
        hist.SetBinContent(nx, 0.0)
        hist.SetBinError(nx, 0.0)
        hist.SetBinContent(nx-1, 0.0)
        hist.SetBinError(nx-1, 0.0)
        hist.SetBinContent(nx+1, 0.0)
        hist.SetBinError(nx+1, 0.0)
    if isinstance(hists, list):
        for hist in hists:
            func(hist)
    else:
        func(hists)
    return hists

def add_underflow(hists):
    def func(hist):
        hist.SetBinContent(1, hist.GetBinContent(1)+hist.GetBinContent(0))
        hist.SetBinError(1, math.sqrt(hist.GetBinError(1)*hist.GetBinError(1)+hist.GetBinError(0)*hist.GetBinError(0)))
    if isinstance(hists, list):
        for hist in hists:
            func(hist)
    else:
        func(hists)
    return hists



def makeplot_single_2d(
    sig_fnames_=None,
    bkg_fnames_=None,
    data_fname_=None,
    sig_legends_=None,
    bkg_legends_=None,
    sig_colors_=None,
    bkg_colors_=None,
    hist_name_=None,
    dir_name_="plots",
    output_name_=None,
    extraoptions=None
    ):

    if sig_fnames_ == None or bkg_fnames_ == None or hist_name_ == None or sig_legends_ == None or bkg_legends_ == None:
        print("nothing to plot....")
        return
    print("making plot for "+hist_name_)
    for idx in range(len(sig_legends_)):
        sig_legends_[idx] = sig_legends_[idx].split(" x ")[0]
    for idx in range(len(bkg_legends_)):
        bkg_legends_[idx] = bkg_legends_[idx].split(" x ")[0]
   
    s_color = [632, 617, 839, 800, 1]
    b_color = [920, 2007, 2005, 2003, 2001, 2011]
    if sig_colors_:
        s_color = sig_colors_
    if bkg_colors_:
        b_color = bkg_colors_
    
    tfs_sig = {}
    h2_sig = []
    maxY = 0.0
    for idx in range(len(sig_fnames_)): 
        fn = sig_fnames_[idx]
        n = os.path.basename(fn.replace(".root", ""))
        tfs_sig[n] = r.TFile(fn)
        h2 = tfs_sig[n].Get(hist_name_)
        h2.SetName(hist_name_+"_sig_"+str(idx))
        h2.SetMarkerColor(s_color[idx])
        h2.SetMarkerStyle(7)
        #h2.SetMarkerSize(0.1)
        #print(h2.Integral())
        h2_sig.append(h2)
     
    tfs_bkg = {}
    h2_bkg = []
    for idx in range(len(bkg_fnames_)): 
        fn = bkg_fnames_[idx]
        n = os.path.basename(fn.replace(".root", ""))
        tfs_bkg[n] = r.TFile(fn)
        h2 = tfs_bkg[n].Get(hist_name_)
        h2.SetName(hist_name_+"_bkg_"+str(idx))
        h2.SetMarkerColor(b_color[idx])
        h2.SetMarkerStyle(7)
        #h2.SetMarkerSize(0.1)
        #print(h2.Integral())
        h2_bkg.append(h2)
    h2_data = None 

    tfs_data = {}
    if data_fname_:
        tfs_data["fdata"] = r.TFile(data_fname_)
        h2_data = tfs_data["fdata"].Get(hist_name_)
        h2_data.SetName(hist_name_+"_data")
        h2_data.SetMarkerColor(1)
        h2_data.SetMarkerStyle(7)
        #print(h2_data.Integral())
    myC = r.TCanvas("myC","myC", 600, 600)
    myC.SetTicky(1)
    myC.SetRightMargin( 0.05 )
    myC.SetLeftMargin( leftMargin ) 

    
    h2_sig[0].Draw("scat=10.0")
    for idx in range(1, len(h2_sig)):
        h2_sig[idx].Draw("same scat=10.0")
    for idx in range(len(h2_bkg)):
        h2_bkg[idx].Draw("same scat=10.0")

    h2_sig[0].SetTitle("")
   
    if h2_data:
        h2_data.Draw("same scat=10.0")
    #h2_sig[0].GetZaxis().SetTitle("Events")
    h2_sig[0].GetYaxis().SetTitleOffset(1.1)
    h2_sig[0].GetYaxis().SetTitleSize(0.055)
    h2_sig[0].GetYaxis().SetLabelSize(0.045)
    h2_sig[0].GetYaxis().CenterTitle()

    leg = r.TLegend(leftMargin, 0.7, leftMargin+0.20, 0.9)
    #leg.SetNColumns(2)
    #leg.SetFillStyle(1)
    leg.SetFillColor(0)
    leg.SetBorderSize(1)
    leg.SetTextFont(32)
    for idx in range(len(h2_sig)):
        leg.AddEntry(h2_sig[idx], sig_legends_[idx], "p")
    for idx in range(len(h2_bkg)):
        leg.AddEntry(h2_bkg[idx], bkg_legends_[idx], "p")
    if h2_data:
        leg.AddEntry(h2_data, "Data", "p")
    leg.Draw()
    h2_sig[0].GetXaxis().SetTitleOffset(0.8)
    h2_sig[0].GetXaxis().SetTitleSize(0.055)
    h2_sig[0].GetXaxis().SetLabelSize(0.045)
    h2_sig[0].GetXaxis().CenterTitle()

    if "xaxis_label" in extraoptions and extraoptions["xaxis_label"] != None:
        x_title = extraoptions["xaxis_label"]
        h2_sig[0].GetXaxis().SetTitle(x_title)
    if "yaxis_label" in extraoptions and extraoptions["yaxis_label"] != None:
        y_title = extraoptions["yaxis_label"]
        h2_sig[0].GetYaxis().SetTitle(y_title)

    ##########draw CMS preliminary
    tex1 = r.TLatex(leftMargin, 0.91, "CMS")
    tex1.SetNDC()
    tex1.SetTextFont(61)
    tex1.SetTextSize(0.070)
    tex1.SetLineWidth(2)
    tex1.Draw()
    tex2 = r.TLatex(leftMargin+0.15,0.912,"Preliminary")
    tex2.SetNDC()
    tex2.SetTextFont(52)
    tex2.SetTextSize(0.055)
    tex2.SetLineWidth(2)
    tex2.Draw()
    
    lumi_value = 137
    if "lumi_value" in extraoptions:
        lumi_value = extraoptions["lumi_value"]
    tex3 = r.TLatex(0.6,0.912,"%d"%lumi_value+" fb^{-1} (13 TeV)")
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

    
    myC.SaveAs(outFile+".png")
    myC.SaveAs(outFile+".pdf")
    myC.SaveAs(outFile+".C")
 
def makeplot_fitmass(
    data_fname_=None,
    hist_name_=None,
    dir_name_="plots",
    output_name_=None,
    extraoptions=None
    ):

    if data_fname_ == None or hist_name_ == None: 
        print("nothing to plot....")
        return
    if "msoft" not in hist_name_:
        return
    print("making plot for "+hist_name_)
    print(data_fname_)
   
    tfs_data = r.TFile(data_fname_)
    h1_data = tfs_data.Get(hist_name_)
    h1_data.SetName(hist_name_+"_data")
    h1_data.SetBinErrorOption(1)
    h1_data.SetLineColor(1)
    h1_data.SetLineWidth(2)
    h1_data.SetMarkerColor(1)
    h1_data.SetMarkerStyle(20)
    #print(h1_data.Integral())
    if "nbins" in extraoptions:
        rebin([h1_data], extraoptions["nbins"])   
    if "remove_underflow" in extraoptions:
        if extraoptions["remove_underflow"]:
            remove_underflow([h1_data])
    if "remove_overflow" in extraoptions:
        if extraoptions["remove_overflow"]:
            remove_overflow([h1_data])
    myC = r.TCanvas("myC","myC", 600, 600)
    myC.SetTicky(1)
    myC.SetRightMargin( rightMargin )
    myC.SetLeftMargin( leftMargin ) 
    myC.SetBottomMargin(0.14)

   
    h1_data.SetTitle("")
    maxY = h1_data.GetMaximum()
    h1_data.Draw("PEX0")
    h1_data.GetYaxis().SetTitle("Events")
    h1_data.GetYaxis().SetTitleOffset(0.85)
    h1_data.GetYaxis().SetTitleSize(0.07)
    h1_data.GetYaxis().SetLabelSize(0.045)
    h1_data.GetYaxis().CenterTitle()

    h1_data.GetXaxis().SetTitleOffset(0.94)
    h1_data.GetXaxis().SetTitleSize(0.06)
    h1_data.GetXaxis().SetLabelSize(0.045)
    h1_data.GetXaxis().SetLabelOffset(0.013)

    if "xaxis_label" in extraoptions and extraoptions["xaxis_label"] != None:
        x_title = extraoptions["xaxis_label"]
        h1_data.GetXaxis().SetTitle(x_title)

    f1 = r.TF1("gaus", "gaus(0)", 50, 170)
    f1.SetLineColor(r.kRed)
    h1_data.Fit(f1, "", "", 50, 170)
   
    N_bkg = f1.Integral(95.0, 135.0)/h1_data.GetBinWidth(2)
    tex_nb = r.TLatex(0.20,0.80,"N_{b} = #int_{95}^{135} = %6.2f"%N_bkg)
    tex_nb.SetNDC()
    tex_nb.SetTextFont(42)
    tex_nb.SetTextSize(0.035)
    tex_nb.SetLineWidth(2)
    tex_nb.Draw()

    ##########draw CMS preliminary
    tex1 = r.TLatex(leftMargin, 0.91, "CMS")
    tex1.SetNDC()
    tex1.SetTextFont(61)
    tex1.SetTextSize(0.070)
    tex1.SetLineWidth(2)
    tex1.Draw()
    tex2 = r.TLatex(leftMargin+0.15,0.912,"Preliminary")
    tex2.SetNDC()
    tex2.SetTextFont(52)
    tex2.SetTextSize(0.055)
    tex2.SetLineWidth(2)
    tex2.Draw()
    
    lumi_value = 137
    if "lumi_value" in extraoptions:
        lumi_value = extraoptions["lumi_value"]
    tex3 = r.TLatex(0.63,0.912,"%d"%lumi_value+" fb^{-1} (13 TeV)")
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
    h1_data.SetMaximum(maxY*1.5)
    
    myC.SaveAs(outFile+"_fit.png")
    myC.SaveAs(outFile+"_fit.pdf")
    myC.SaveAs(outFile+"_fit.C")
  
def makeplot_single(
    sig_fnames_=None,
    bkg_fnames_=None,
    data_fname_=None,
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

    if sig_fnames_ == None or bkg_fnames_ == None or hist_name_ == None or sig_legends_ == None or bkg_legends_ == None:
        print("nothing to plot....")
        return

    if "_v_" in hist_name_:
        makeplot_single_2d(sig_fnames_=sig_fnames_, bkg_fnames_ = bkg_fnames_, data_fname_=data_fname_, sig_legends_=sig_legends_, bkg_legends_=bkg_legends_, sig_colors_=sig_colors_, bkg_colors_=bkg_colors_, hist_name_=hist_name_, dir_name_=dir_name_, extraoptions=extraoptions)
        return

    print("making plot for "+hist_name_)
    print(sig_fnames_)
    print(bkg_fnames_)
   
    s_color = [632, 617, 839, 800, 1]
    b_color = [920, 2007, 2005, 2003, 2001, 2011]
    if sig_colors_:
        s_color = sig_colors_
    if bkg_colors_:
        b_color = bkg_colors_
    
    tfs_sig = {}
    h1_sig = []
    maxY = 0.0
    if data_fname_ != None:
        extraoptions["stack_signal"]=True
        sig_scale_ = 1.0
    if "stack_signal" in extraoptions and extraoptions["stack_signal"]:
        sig_scale_ = 1.0
        for idx in range(len(sig_legends_)):
            sig_legends_[idx] = sig_legends_[idx].split(" x ")[0]
    for idx in range(len(sig_fnames_)): 
        fn = sig_fnames_[idx]
        n = os.path.basename(fn.replace(".root", ""))
        tfs_sig[n] = r.TFile(fn)
        h1 = tfs_sig[n].Get(hist_name_)
        h1.SetName(hist_name_+"_sig_"+str(idx))
        h1.Scale(sig_scale_)
        h1.SetLineWidth(2)
        h1.SetLineColor(s_color[idx])
        print("signal index and sum:",idx, h1.Integral())
        h1_sig.append(h1)
     
    tfs_bkg = {}
    h1_bkg = []
    for idx in range(len(bkg_fnames_)): 
        fn = bkg_fnames_[idx]
        n = os.path.basename(fn.replace(".root", ""))
        tfs_bkg[n] = r.TFile(fn)
        print("fn hist_name_",fn,hist_name_)
        h1 = tfs_bkg[n].Get(hist_name_)
        print("sum 1:",hist_name_,h1.Integral())
        h1.SetName(hist_name_+"_bkg_"+str(idx))
        h1.SetLineWidth(2)
        h1.SetLineColor(b_color[idx])
        h1.SetFillColorAlpha(b_color[idx], 1)
        print("sum 2:",hist_name_,h1.Integral())
        h1_bkg.append(h1)
    h1_data = None 

    tfs_data = {}
    if data_fname_:
        tfs_data["fdata"] = r.TFile(data_fname_)
        h1_data = tfs_data["fdata"].Get(hist_name_)
        h1_data.SetName(hist_name_+"_data")
        h1_data.SetBinErrorOption(1)
        h1_data.SetLineColor(1)
        h1_data.SetLineWidth(2)
        h1_data.SetMarkerColor(1)
        h1_data.SetMarkerStyle(20)
        print("sum data: ",h1_data.Integral())
    if "nbins" in extraoptions:
        rebin(h1_sig, extraoptions["nbins"])   
        rebin(h1_bkg, extraoptions["nbins"])   
        if h1_data:
            rebin([h1_data], extraoptions["nbins"])   
    if "remove_underflow" in extraoptions:
        if extraoptions["remove_underflow"]:
            remove_underflow(h1_sig)
            remove_underflow(h1_bkg)
            if h1_data:
                remove_underflow([h1_data])
    if "remove_overflow" in extraoptions:
        if extraoptions["remove_overflow"]:
            remove_overflow(h1_sig)
            remove_overflow(h1_bkg)
            if h1_data:
                remove_overflow([h1_data])
    if "add_overflow" in extraoptions:
        if extraoptions["add_overflow"]:
            add_overflow(h1_sig)
            add_overflow(h1_bkg)
            if h1_data:
                add_overflow([h1_data])
    if "add_underflow" in extraoptions:
        if extraoptions["add_underflow"]:
            add_underflow(h1_sig)
            add_underflow(h1_bkg)
            if h1_data:
                add_underflow([h1_data])
    if "blind_data" in extraoptions:
        if extraoptions["blind_data"]:
            #blind_data(h1_sig)
            #blind_data(h1_bkg)
            if h1_data:
                blind_data([h1_data])

    myC = r.TCanvas("myC","myC", 600, 600)
    myC.SetTicky(1)
    pad1 = r.TPad("pad1","pad1", 0.05, 0.33,0.95, 0.97) 
    pad1.SetBottomMargin(0)
    pad1.SetRightMargin( rightMargin )
    pad1.SetLeftMargin( leftMargin ) 
    pad2 = r.TPad("pad2","pad2", 0.05, 0.045, 0.95, 0.33)
    pad2.SetBottomMargin(0.4)
    pad2.SetTopMargin(0.008)
    pad2.SetRightMargin( rightMargin )
    pad2.SetLeftMargin( leftMargin )

    pad2.Draw()
    pad1.Draw()

    pad1.cd()

    stack = r.THStack("stack", "stack")
    #signal process exist
    #hist_all = h1_sig[0].Clone("hist_all")
    #hist_all.Scale(1.0/sig_scale_)
    hist_all = h1_bkg[0].Clone("hist_all")
    hist_all.Scale(0.0)

    #no signal process
    hist_s = hist_all.Clone("hist_s")
    #signal process exist
    #hist_s = h1_sig[0].Clone("hist_s")

    hist_b = h1_bkg[0].Clone("hist_b")
    for idx in range(len(h1_bkg)):
        print("stack add:", h1_bkg[idx].Integral())
        stack.Add(h1_bkg[idx])
        hist_all.Add(h1_bkg[idx])
        if idx > 0:
            hist_b.Add(h1_bkg[idx])
    for idx in range(len(h1_sig)):
        if idx > -1:
            hist_temp = h1_sig[idx].Clone(h1_sig[idx].GetName()+"_temp")
            hist_temp.Scale(1.0/sig_scale_)
            hist_all.Add(hist_temp)
            hist_s.Add(h1_sig[idx])

    stack.SetTitle("")
    if "stack_signal" in extraoptions and extraoptions["stack_signal"]:
        for idx in range(len(h1_sig)):
            h1_sig[idx].SetFillColorAlpha(s_color[idx], 1)
            stack.Add(h1_sig[idx])
        if stack.GetMaximum() > maxY:
            maxY = stack.GetMaximum()
        stack.Draw("hist")
    else:
        stack.Draw("hist")
        if stack.GetMaximum() > maxY:
            maxY = stack.GetMaximum()
        for idx in range(len(h1_sig)):
            if h1_sig[idx].GetMaximum() > maxY:
                maxY = h1_sig[idx].GetMaximum()
            h1_sig[idx].Draw("samehist")
    
    if h1_data:
        if h1_data.GetMaximum() > maxY:
            maxY = h1_data.GetMaximum()
        h1_data.Draw("samePEX0")
    stack.GetYaxis().SetTitle("Events")
    stack.GetYaxis().SetTitleOffset(0.88)
    stack.GetYaxis().SetTitleSize(0.08)
    stack.GetYaxis().SetLabelSize(0.045)
    stack.GetYaxis().CenterTitle()
    #if "xaxis_range" in extraoptions:
    #    stack.GetXaxis().SetRangeUser(float(extraoptions["xaxis_range"][0]),float(extraoptions["xaxis_range"][1]))

    leg = r.TLegend(0.16, 0.70, 0.97, 0.88)
    leg.SetNColumns(4)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.07)
    
    for idx in range(len(h1_bkg)):
        leg.AddEntry(h1_bkg[idx], bkg_legends_[idx], "F")
    for idx in range(len(h1_sig)):
        if "stack_signal" in extraoptions and extraoptions["stack_signal"]:
            leg.AddEntry(h1_sig[idx], sig_legends_[idx], "F")
        else:
            leg.AddEntry(h1_sig[idx], sig_legends_[idx], "L")
    if h1_data:
        leg.AddEntry(h1_data, "Data", "ep")

    leg.Draw()
    pad2.cd()
    pad2.SetGridy(1)
    ratio = None
    ratio_Low  = 0.0
    ratio_High  = 2.0
    if h1_data:
        ratio = h1_data.Clone("ratio_data_over_mc")
        ratio.Divide(hist_all)
        if "ratio_range" in extraoptions:
            ratio_Low = extraoptions["ratio_range"][0]
            ratio_High = extraoptions["ratio_range"][1]
        ratio.GetYaxis().SetTitle("data/mc")
        ratio.GetYaxis().SetRangeUser(ratio_Low, ratio_High)
        ratio.Draw("same PE")
    else:
        ratio = h1_sig[0].Clone("ratio")
        ratio_High = 0.0
        for ibin in range(1,ratio.GetNbinsX()+1):
            s = hist_s.GetBinContent(ibin) / sig_scale_
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
    #if "xaxis_range" in extraoptions:
    #    ratio.GetXaxis().SetRangeUser(float(extraoptions["xaxis_range"][0]),float(extraoptions["xaxis_range"][1]))
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

    if "cutflow" in hist_name_:
        ratio.GetXaxis().SetTitle("")
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
    tex2 = r.TLatex(leftMargin+0.12,0.912,"Preliminary")
    tex2.SetNDC()
    tex2.SetTextFont(52)
    tex2.SetTextSize(0.055)
    tex2.SetLineWidth(2)
    tex2.Draw()
    
    lumi_value = 137
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
    #Y maximum range
    stack.SetMaximum(maxY*1.7)

    #print everything into txt file
    text_file = open(outFile+"_linY.txt", "w")
    text_file.write("bin    |   x    ")
    for idx in range(len(bkg_legends_)):
        text_file.write(" | %21s"%bkg_legends_[idx])
    text_file.write(" | %21s"%("total B"))
    for idx in range(len(sig_legends_)):
        text_file.write(" | %25s"%sig_legends_[idx])
    if h1_data:
        text_file.write(" | data | data/mc")
    text_file.write("\n-------------")
    for idx in range(24*(len(bkg_legends_) + 1)+ 29*len(sig_legends_)):
        text_file.write("-")
    if h1_data:
        text_file.write("-------")
    text_file.write("\n")
    for ibin in range(0,h1_bkg[0].GetNbinsX()+1):
        text_file.write("%3d"%ibin+"   ")
        text_file.write(" | %6.3f"%ratio.GetBinCenter(ibin)+" ")
        for idx in range(len(bkg_legends_)):
            text_file.write(" | %7.3f "%h1_bkg[idx].GetBinContent(ibin)+"$\\pm$"+ " %7.3f"%h1_bkg[idx].GetBinError(ibin))
        text_file.write(" | %7.3f "%hist_b.GetBinContent(ibin)+"$\\pm$"+ " %7.3f"%hist_b.GetBinError(ibin))
        for idx in range(len(sig_legends_)):
            text_file.write(" | %9.3f "%h1_sig[idx].GetBinContent(ibin)+"$\\pm$"+ " %9.3f"%h1_sig[idx].GetBinError(ibin))
        if h1_data:
            text_file.write(" | %d"%h1_data.GetBinContent(ibin) +  " | %7.3f "%ratio.GetBinContent(ibin) +"$\\pm$"+ " %7.3f"%ratio.GetBinError(ibin))
        text_file.write("\n")
    text_file.close()
    os.system("cp "+outFile+"_linY.txt "+outFile+"_logY.txt")

    
    myC.SaveAs(outFile+"_linY.png")
    myC.SaveAs(outFile+"_linY.pdf")
    myC.SaveAs(outFile+"_linY.C")
    pad1.cd()
    stack.SetMaximum(maxY*100.0)
    stack.SetMinimum(0.5)
    pad1.SetLogy()
    myC.SaveAs(outFile+"_logY.png")
    myC.SaveAs(outFile+"_logY.pdf")
    myC.SaveAs(outFile+"_logY.C")
    #save histogram and ratio to root file
    outFile_root = r.TFile(outFile+".root", "recreate")
    outFile_root.cd()
    for idx in range(len(bkg_legends_)):
        h1_bkg[idx].Write()
    for idx in range(len(sig_legends_)):
        h1_sig[idx].Write()
    if  h1_data:
        h1_data.Write()
        ratio.Write()
    #outFile_root.Write()
    outFile_root.Close()

def makeplot_all(
    sig_fnames_=None,
    bkg_fnames_=None,
    data_fname_=None,
    sig_legends_=None,
    bkg_legends_=None,
    sig_colors_=None,
    bkg_colors_=None,
    regions_name_=None,
    sig_scale_=1.0,
    dir_name_="plots",
    output_name_=None,
    extraoptions=None
    ):

    if regions_name_ == None or sig_fnames_ == None:
        print("nothing to plot....")
        return
    print("making all plots in the following regions: ")
    print(regions_name_)

    hist_names = []
    for region in regions_name_:
        tf = r.TFile(bkg_fnames_[0])
        for key in tf.GetListOfKeys():
            if region in key.GetName():
                hist_names.append(key.GetName())
    print("find the following histgrams to plot:")
    print(hist_names)
    for hist_name in hist_names:
        if "_v_" in hist_name:
            makeplot_single_2d(sig_fnames_=sig_fnames_, bkg_fnames_ = bkg_fnames_, data_fname_=data_fname_, sig_legends_=sig_legends_, bkg_legends_=bkg_legends_, sig_colors_=sig_colors_, bkg_colors_=bkg_colors_, hist_name_=hist_name, dir_name_=dir_name_, extraoptions=extraoptions)
        else:
            makeplot_single(sig_fnames_=sig_fnames_, bkg_fnames_ = bkg_fnames_, data_fname_=data_fname_, sig_legends_=sig_legends_, bkg_legends_=bkg_legends_, sig_colors_=sig_colors_, bkg_colors_=bkg_colors_, hist_name_=hist_name, sig_scale_=sig_scale_, dir_name_=dir_name_, extraoptions=extraoptions)

def makeplot_cutOptimize(
    sig_fnames_=None,
    bkg_fnames_=None,
    sig_legends_=None,
    bkg_legends_=None,
    sig_colors_=None,
    bkg_colors_=None,
    hist_name_=None,
    sig_scale_=1.0,
    bkg_scale_=1.0,
    dir_name_="plots",
    output_name_=None,
    right_hand_=False,
    extraoptions=None
    ):

    if sig_fnames_ == None or bkg_fnames_ == None or hist_name_ == None or sig_legends_ == None or bkg_legends_ == None:
        print("nothing to plot....")
        return
   
    s_color = [632, 617, 839, 800, 1]
    b_color = [920, 2007, 2005, 2003, 2001, 2011]
    if sig_colors_:
        s_color = sig_colors_
    if bkg_colors_:
        b_color = bkg_colors_
    
    tfs_sig = {}
    h1_sig = []
    for idx in range(len(sig_fnames_)): 
        fn = sig_fnames_[idx]
        n = os.path.basename(fn.replace(".root", ""))
        tfs_sig[n] = r.TFile(fn)
        h1 = tfs_sig[n].Get(hist_name_)
        h1.SetName(hist_name_+"_sig_"+str(idx))
        h1.Scale(sig_scale_)
        h1.SetLineWidth(2)
        h1.SetLineColor(s_color[idx])
        h1_sig.append(h1)

    tfs_bkg = {}
    h1_bkg = []
    for idx in range(len(bkg_fnames_)): 
        fn = bkg_fnames_[idx]
        n = os.path.basename(fn.replace(".root", ""))
        tfs_bkg[n] = r.TFile(fn)
        h1 = tfs_bkg[n].Get(hist_name_)
        h1.SetName(hist_name_+"_bkg_"+str(idx))
        h1.Scale(bkg_scale_)
        h1.SetLineWidth(2)
        h1.SetLineColor(b_color[idx])
        h1.SetFillColorAlpha(b_color[idx], 1)
        h1_bkg.append(h1)
   
    if "remove_underflow" in extraoptions:
        if extraoptions["remove_underflow"]:
            remove_underflow(h1_sig)
            remove_underflow(h1_bkg)
    if "remove_overflow" in extraoptions:
        if extraoptions["remove_overflow"]:
            remove_overflow(h1_sig)
            remove_overflow(h1_bkg)

    myC = r.TCanvas("myC","myC", 600, 600)
    myC.SetTicky(1)
    pad1 = r.TPad("pad1","pad1", 0.05, 0.33,0.95, 0.97)
    pad1.SetBottomMargin(0)
    pad1.SetRightMargin( rightMargin )
    pad1.SetLeftMargin( leftMargin )
    pad2 = r.TPad("pad2","pad2", 0.05, 0.045, 0.95, 0.33)
    pad2.SetBottomMargin(0.4)
    pad2.SetTopMargin(0.008)
    pad2.SetRightMargin( rightMargin )
    pad2.SetLeftMargin( leftMargin )

    pad2.Draw()
    pad1.Draw()

    pad1.cd()


    hist_noCut_s = h1_sig[0].Clone("hist_noCut_s")
    hist_noCut_b = h1_bkg[0].Clone("hist_noCut_b")
    for idx in range(len(h1_bkg)):
        if idx > 0:
            hist_noCut_b.Add(h1_bkg[idx])
    for idx in range(len(h1_sig)):
        if idx > 0:
            hist_noCut_s.Add(h1_sig[idx])
    nbins = hist_noCut_s.GetNbinsX()
    print("nbins: "+str(nbins))
    print("nS before cut: "+ str(hist_noCut_s.Integral()/sig_scale_))
    print("nB before cut: "+ str(hist_noCut_b.Integral()))
    ratio = hist_noCut_s.Clone("ratio")
    hist_b = hist_noCut_b.Clone("hist_b")
    hist_s = hist_noCut_b.Clone("hist_s")

    maxY1 = 0.0
    maxY2 = 0.0
    bestCut = 0 
    for ibin in range(0, nbins+1):
        s = 0.0
        b = 0.0
        L = 0.0
        if right_hand_: #<=
            s = hist_noCut_s.Integral(0, ibin) 
            b = hist_noCut_b.Integral(0, ibin)
            for idx in range(len(h1_sig)):
                s_this = h1_sig[idx].Integral(0, nbins-ibin+1)
                h1_sig[idx].SetBinContent(nbins-ibin+1, s_this)
            for idx in range(len(h1_bkg)):
                b_this = h1_bkg[idx].Integral(0, nbins-ibin+1)
                h1_bkg[idx].SetBinContent(nbins-ibin+1, b_this)
        else: #>=
            s = hist_noCut_s.Integral(ibin, nbins+1)
            b = hist_noCut_b.Integral(ibin, nbins+1)
            for idx in range(len(h1_sig)):
                s_this = h1_sig[idx].Integral(ibin, nbins+1)
                h1_sig[idx].SetBinContent(ibin, s_this)
            for idx in range(len(h1_bkg)):
                b_this = h1_bkg[idx].Integral(ibin, nbins+1)
                h1_bkg[idx].SetBinContent(ibin, b_this)
        if b > 0.0:
            L = s/(sig_scale_*math.sqrt(b))
        if s > maxY1:
            maxY1 = s
        if b > maxY1:
            maxY1 = b
        if L > maxY2:
            maxY2 = L
            bestCut = ibin
        ratio.SetBinContent(ibin, L)
        hist_b.SetBinContent(ibin, b)
        hist_s.SetBinContent(ibin, s)
    
    x_title = "cut (>=) on "+hist_noCut_s.GetXaxis().GetTitle()
    if right_hand_:
        x_title = "cut (<=) on "+hist_noCut_s.GetXaxis().GetTitle()
    if "xaxis_label" in extraoptions and extraoptions["xaxis_label"] != None:
        x_title = extraoptions["xaxis_label"]

    stack = r.THStack("stack", "stack")
    for idx in range(len(h1_bkg)):
        stack.Add(h1_bkg[idx])
    stack.SetTitle("")
    stack.Draw("hist")
    for idx in range(len(h1_sig)):
        h1_sig[idx].Draw("samehist")
    stack.GetYaxis().SetTitle("Events")
    stack.GetYaxis().SetTitleOffset(0.88)
    stack.GetYaxis().SetTitleSize(0.08)
    stack.GetYaxis().SetLabelSize(0.045)
    stack.GetYaxis().CenterTitle()

    leg = r.TLegend(0.15, 0.70, 0.95, 0.88)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(28)
    for idx in range(len(h1_bkg)):
        leg.AddEntry(h1_bkg[idx], bkg_legends_[idx], "F")
    for idx in range(len(h1_sig)):
        leg.AddEntry(h1_sig[idx], sig_legends_[idx], "L")
    leg.Draw()

    pad2.cd()
    pad2.SetGridy(1)
    
    if maxY2 > 1.0:
        maxY2 = 1.0
    ratio.GetYaxis().SetRangeUser(0.0, 1.2*maxY2)
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
    ratio.GetXaxis().SetTitle(x_title)

    ##########draw CMS preliminary
    pad1.cd()
    tex1 = r.TLatex(leftMargin, 0.91, "CMS")
    tex1.SetNDC()
    tex1.SetTextFont(61)
    tex1.SetTextSize(0.070)
    tex1.SetLineWidth(2)
    tex1.Draw()
    tex2 = r.TLatex(leftMargin+0.12,0.912,"Preliminary")
    tex2.SetNDC()
    tex2.SetTextFont(52)
    tex2.SetTextSize(0.055)
    tex2.SetLineWidth(2)
    tex2.Draw()

    lumi_value = 137
    if "lumi_value" in extraoptions:
        lumi_value = extraoptions["lumi_value"]
    tex3 = r.TLatex(0.72,0.912,"%d"%lumi_value+" fb^{-1} (13 TeV)")
    tex3.SetNDC()
    tex3.SetTextFont(42)
    tex3.SetTextSize(0.055)
    tex3.SetLineWidth(2)
    tex3.Draw()

    ##print out best cut and s, b, L
    c_best = ratio.GetBinCenter(bestCut)
    s_best = hist_s.GetBinContent(bestCut) / sig_scale_
    b_best = hist_b.GetBinContent(bestCut)
    L_best = ratio.GetBinContent(bestCut)
    print("best cut: %6.3f"%c_best)
    print("nS at best cut: %7.3f"%s_best)
    print("nB at best cut: %7.3f"%b_best)
    print("s/sqrt(b) at best cut: %7.3f"%L_best)

    tex4 = r.TLatex(leftMargin+0.05, 0.85, "best cut = %6.2f"%c_best)
    tex4.SetNDC()
    tex4.SetTextFont(42)
    tex4.SetTextSize(0.055)
    tex4.SetLineWidth(2)
    tex4.Draw()

    tex5 = r.TLatex(leftMargin+0.05, 0.80, "s = %6.2f"%s_best)
    tex5.SetNDC()
    tex5.SetTextFont(42)
    tex5.SetTextSize(0.055)
    tex5.SetLineWidth(2)
    tex5.Draw()

    tex6 = r.TLatex(leftMargin+0.05, 0.75, "b = %6.2f"%b_best)
    tex6.SetNDC()
    tex6.SetTextFont(42)
    tex6.SetTextSize(0.055)
    tex6.SetLineWidth(2)
    tex6.Draw()

    tex7 = r.TLatex(leftMargin+0.05, 0.68, "s/#sqrt{b} = %6.2f"%L_best)
    tex7.SetNDC()
    tex7.SetTextFont(42)
    tex7.SetTextSize(0.055)
    tex7.SetLineWidth(2)
    tex7.Draw()

    #draw vertical line indcating best cut

    pad1.cd()
    line1 = r.TLine(c_best,0.001,c_best,0.9*maxY1)
    line1.SetLineStyle(7)
    line1.Draw("same")
    pad2.cd()
    line2 = r.TLine(c_best,0.0,c_best,1.2*maxY2)
    line2.SetLineStyle(7)
    line2.Draw("same")

    outFile = dir_name_
    if output_name_:
        outFile = outFile + "/" +output_name_
    else:
        outFile = outFile + "/" + hist_name_

    stack.SetMaximum(maxY1*1.3)

    #print everything into txt file
    text_file = open(outFile+"_linY.txt", "w")
    text_file.write("bin    |   x    ")
    for idx in range(len(bkg_legends_)):
        text_file.write(" | %21s"%bkg_legends_[idx])
    text_file.write(" | %21s"%("total B"))
    for idx in range(len(sig_legends_)):
        text_file.write(" | %25s"%sig_legends_[idx])
    text_file.write(" | %21s"%("S/sqrt(B)"))
    text_file.write("\n-------------")
    for idx in range(24*(len(bkg_legends_) + 2)+ 29*len(sig_legends_)):
        text_file.write("-")
    text_file.write("\n")
    for ibin in range(0,h1_sig[0].GetNbinsX()+1):
        text_file.write("%3d"%ibin+"   ")
        text_file.write(" | %6.3f"%ratio.GetBinCenter(ibin)+" ")
        for idx in range(len(bkg_legends_)):
            text_file.write(" | %7.3f "%h1_bkg[idx].GetBinContent(ibin)+"$\\pm$"+ " %7.3f"%h1_bkg[idx].GetBinError(ibin))
        text_file.write(" | %7.3f "%hist_b.GetBinContent(ibin)+"$\\pm$"+ " %7.3f"%hist_b.GetBinError(ibin))
        for idx in range(len(sig_legends_)):
            text_file.write(" | %9.3f "%h1_sig[idx].GetBinContent(ibin)+"$\\pm$"+ " %9.3f"%h1_sig[idx].GetBinError(ibin))
        text_file.write(" | %9.3f "%ratio.GetBinContent(ibin))
        text_file.write("\n")
    text_file.close()
    os.system("cp "+outFile+"_linY.txt "+outFile+"_logY.txt")

    myC.SaveAs(outFile+"_linY.png")
    myC.SaveAs(outFile+"_linY.pdf")
    myC.SaveAs(outFile+"_linY.C")
    pad1.cd()
    stack.SetMaximum(maxY1*1000.0)
    stack.SetMinimum(0.01)
    pad1.SetLogy()
    myC.SaveAs(outFile+"_logY.png")
    myC.SaveAs(outFile+"_logY.pdf")
    myC.SaveAs(outFile+"_logY.C")





    

