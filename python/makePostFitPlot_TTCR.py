from __future__ import print_function
import numpy as np
import math
import ROOT as r
import shlex
import sys
import os
from colors  import *

leftMargin   = 0.13
rightMargin  = 0.01
topMargin    = 0.07
bottomMargin = 0.12
r.gStyle.SetOptStat(0)
r.gStyle.SetOptFit(111)
r.gStyle.SetStatX(0.99)
r.gStyle.SetStatY(0.9)
r.gStyle.SetStatW(0.2)
r.gStyle.SetStatH(0.15)
r.gROOT.SetBatch(1)


def makeplot_single(
    h1_sig=None,
    h1_bkg=None,
    h1_data=None,
    sig_legends_=None,
    bkg_legends_=None,
    sig_colors_=None,
    bkg_colors_=None,
    lumi_=137,
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
        h1_sig[idx].SetLineWidth(2)
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
    nS = np.zeros(h1_bkg[0].GetNbinsX())
    eS = np.zeros(h1_bkg[0].GetNbinsX())
    hist_all = h1_sig[0].Clone("hist_all")
    hist_all.Scale(1.0/sig_scale_)
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

    #print("nS after h_bkg ", nS)
    #print("nS after h_bkg ", eS)

    for idx in range(len(h1_sig)):
        print("idx signal ",idx)
        for ib in range(h1_sig[0].GetNbinsX()):
            #print("before ",nS[ib], h1_sig[idx].GetBinContent(ib+1), h1_sig[idx].GetBinError(ib+1))
            nS[ib] += h1_sig[idx].GetBinContent(ib+1)
            #print("after ",nS[ib])
            eS[ib] = math.sqrt(eS[ib]*eS[ib] + h1_sig[idx].GetBinError(ib+1)*h1_sig[idx].GetBinError(ib+1))

    #print("nS after h_sig ", nS)
    #print("nS after h_sig ", eS)

    for idx in range(len(h1_sig)):
        if idx > 0:
            hist_temp = h1_sig[idx].Clone(h1_sig[idx].GetName()+"_temp")
            hist_temp.Scale(1.0/sig_scale_)
            hist_all.Add(hist_temp)
            hist_s.Add(h1_sig[idx])

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
        stack.Draw("hist")
    else:
        stack.Draw("hist")
        if stack.GetMaximum() > maxY:
            maxY = stack.GetMaximum()
        for idx in range(len(h1_sig)):
            if h1_sig[idx].GetMaximum() > maxY:
                maxY = h1_sig[idx].GetMaximum()
            h1_sig[idx].Draw("samehist")
    ##draw  stack total unc
    box = r.TBox(0,0,1,1,)
    box.SetFillStyle(3001)
    box.SetLineWidth(0)
    box.SetFillColor(r.kBlack)
    for idx in range(h1_bkg[0].GetNbinsX()):
        box.DrawBox(h1_bkg[0].GetBinCenter(idx+1)-0.5*h1_bkg[0].GetBinWidth(idx+1), nS[idx]-eS[idx], h1_bkg[0].GetBinCenter(idx+1)+0.5*h1_bkg[0].GetBinWidth(idx+1), nS[idx]+eS[idx])

    if h1_data:
        if h1_data.GetMaximum() > maxY:
            maxY = h1_data.GetMaximum()
        h1_data.Draw("samePEX0")

    #for idx in range(h1_sig[0].GetNbinsX()):
    #    print("%d"%idx+", %6.3f"%eS[idx])
    stack.GetYaxis().SetTitle("Events")
    stack.GetYaxis().SetTitleOffset(0.88)
    stack.GetYaxis().SetTitleSize(0.08)
    stack.GetYaxis().SetLabelSize(0.045)
    stack.GetYaxis().CenterTitle()
    #if "xaxis_range" in extraoptions:
    #    stack.GetXaxis().SetRangeUser(float(extraoptions["xaxis_range"][0]),float(extraoptions["xaxis_range"][1]))

    leg = r.TLegend(0.16, 0.60, 0.97, 0.88)
    leg.SetNColumns(3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.07)
    for idx in range(len(h1_bkg)):
        leg.AddEntry(h1_bkg[idx], bkg_legends_[idx], "F")
    #for idx in range(len(h1_sig)):
    #    if "stack_signal" in extraoptions and extraoptions["stack_signal"]:
    #        leg.AddEntry(h1_sig[idx], sig_legends_[idx], "F")
    #    else:
    #        leg.AddEntry(h1_sig[idx], sig_legends_[idx], "L")
    leg.AddEntry(box, "Total  unc", "F")
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
        ratio.SetTitle("")
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

    if "xaxis_label" in extraoptions and extraoptions["xaxis_label"] != None:
        x_title = extraoptions["xaxis_label"]
        ratio.GetXaxis().SetTitle(x_title)
    ratio.GetYaxis().CenterTitle()

    box_ratio = r.TBox(0,0,1,1,)
    box_ratio.SetFillStyle(3002)
    box_ratio.SetLineWidth(0)
    box_ratio.SetFillColor(r.kBlack)
    for idx in range(h1_bkg[0].GetNbinsX()):
        if np.fabs(nS[idx])> 1e-06:
            box_ratio.DrawBox(h1_bkg[0].GetBinCenter(idx+1)-0.5*h1_bkg[0].GetBinWidth(idx+1), (nS[idx]-eS[idx])/nS[idx], h1_bkg[0].GetBinCenter(idx+1)+0.5*h1_bkg[0].GetBinWidth(idx+1), (nS[idx]+eS[idx])/nS[idx])

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

    lumi_value = lumi_
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
    stack.SetMaximum(maxY*1.6)

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
        text_file.write(" | %6.3f"%ratio.GetBinCenter(ibin)+" ")
        for idx in range(len(h1_bkg)):
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
    for idx in range(len(h1_bkg)):
        h1_bkg[idx].Write()
    for idx in range(len(sig_legends_)):
        h1_sig[idx].Write()
    if  h1_data:
        h1_data.Write()
        ratio.Write()
    #outFile_root.Write()
    outFile_root.Close()

##from this the main code start
if __name__ == "__main__":
    inputfile = sys.argv[1]

lumi = 137
if "2016" in inputfile:
    lumi = 35
elif "2017" in inputfile:
    lumi = 41
elif "2018" in inputfile:
    lumi = 59
    
vbdt = "v8p2"
for dirName in ["shapes_prefit",  "shapes_fit_s"]:

    pnames_sig = ["HH"]
    pnames_bkg = ["others", "QCD", "TTJets"]
    bkg_legends = ["others", "QCD", "t#bar{t}+jets"]
    sig_legends = ["HH"]
    pname_data = "data"
    sig_colors = [2, 839, 800, 1, 632]
    bkg_colors = [2001, 2007, 2005, 800, 839]

    bdtbins = ["TTBarCR"]
        
    for bdtbin in bdtbins:
        fin = r.TFile(inputfile, "READ")
        h1_sig = []
        h1_bkg = []
        h1_data = None

        dirpre = fin.GetDirectory(dirName)
        dirthis = dirpre.GetDirectory("TTBarCR")

        for idx_bkg in range(len(pnames_bkg)):
            if dirthis.GetListOfKeys().Contains(pnames_bkg[idx_bkg]):
                h = fin.Get(dirName+"/TTBarCR/"+pnames_bkg[idx_bkg])
                h.Scale(10.0)
                h1_bkg.append(h)
                print("bkg: ",pnames_bkg[idx_bkg],h.Integral())
                h1_data = h.Clone("h1_data")
            else:
                h = h1_bkg[0].Clone("h1_"+pnames_bkg[idx_bkg])
                h.Scale(0.0)
                h1_bkg.append(h)
        for idx_sig in range(len(pnames_sig)):
            if dirthis.GetListOfKeys().Contains(pnames_sig[idx_sig]):
                h = fin.Get(dirName+"/TTBarCR/"+pnames_sig[idx_sig])
                h1_sig.append(h)
            else:
                h = h1_bkg[0].Clone("h1_"+pnames_sig[idx_sig])
                h.Scale(0.0)
                h1_sig.append(h)

        nBins =  h1_data.GetNbinsX()
        if dirthis.GetListOfKeys().Contains("data"):
            g = fin.Get(dirName+"/TTBarCR/data")
            x,y = r.Double(0), r.Double(0)
            for idx in range(nBins):
                g.GetPoint(idx, x, y)
                h1_data.SetBinContent(idx+1, y*10)
        #h1_data.GetXaxis().SetTitle("j_{2} soft drop mass (GeV)")
        h1_data.GetXaxis().SetTitle("p_{T}^{j_{1}} (GeV)")

        makeplot_single(
                h1_sig=h1_sig,
                h1_bkg=h1_bkg,
                h1_data=h1_data,
                sig_legends_=sig_legends,
                bkg_legends_=bkg_legends,
                sig_colors_=sig_colors,
                bkg_colors_=bkg_colors,
                lumi_=lumi,
                output_name_=dirName+str(lumi)+"_"+bdtbin+"_BDT"+vbdt,
                #dir_name_= "cards_shapes_v8p2_TTBarCR/cards_shapes_TTBarCR/HHModel/",
                extraoptions={"stack_signal": False}
                )

