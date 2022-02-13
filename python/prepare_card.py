import ROOT  as r
import sys



if __name__ == "__main__":
    tag = sys.argv[1]
    vbdt = sys.argv[2]

    proc  =      ["Data", "QCD", "TTJets", "H", "VH",  "ttH", "others",  "HH", "HHc0",  "HHc5", "HHc2p45"]
    proc_file  = ["data", "qcd", "ttbar", "Higgs", "VH",  "ttH", "others",  "HHc1", "HHc0",  "HHc5", "HHc2p45"]

    systs = ["BDT"+vbdt+"MassShape", "ttJetsCorr", "BDT"+vbdt+"Shape", "PNetShape"]
    outName = "HHTo4BPlots_Run2_BDT"+vbdt+".root"
    if "ttbar" in tag:
        outName = "HHTo4BPlots_Run2_ttbarSkim_BDT"+vbdt+".root"
    outFile =  r.TFile(outName, "recreate")

    for idx in range(len(proc)):
        inFile_this = r.TFile("../hists/"+tag+"/combine/"+proc_file[idx]+".root",  "READ")
        print("read file "+"../hists/"+tag+"/combine/"+proc_file[idx]+".root")
        region_list = []
        if "ttbar" in tag:
            region_list =  ["TTBarCR", "TTBarCRTight"]
        else:
            if vbdt == "v24":
                region_list = ["SRv24Bin1", "SRv24Bin2",  "SRv24Bin3", "SRv24Bin4", "FailSRv24", "FitCRv24", "FailFitCRv24"]
            else:
                region_list = ["SRv8p2Bin1", "SRv8p2Bin2",  "SRv8p2Bin3", "FailSRv8p2", "FitCRv8p2", "FailFitCRv8p2"]
        for region in region_list:
            inFile_this.cd()
            hist_nominal = inFile_this.Get(region+"__fatJet2MassSD")
            outBinName=region.replace("SR",  "").replace("Fail", "fail").replace(vbdt, "")
            hist_nominal.SetName("histJet2Mass_"+outBinName+"_"+proc[idx])

            hists_sys = []
            for sys in systs:
                hist_Up = inFile_this.Get(region+sys+"Up__fatJet2MassSD")
                hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                hists_sys.append(hist_Up)
                hist_Down = inFile_this.Get(region+sys+"Down__fatJet2MassSD")
                hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                hists_sys.append(hist_Down)

            hist_JMSUp =  inFile_this.Get(region+"__fatJet2MassSD_JMSUp")
            hist_JMSUp.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMSUp")
            hists_sys.append(hist_JMSUp)
            hist_JMSDown =  inFile_this.Get(region+"__fatJet2MassSD_JMSDown")
            hist_JMSDown.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMSDown")
            hists_sys.append(hist_JMSDown)

            hist_JMRUp =  inFile_this.Get(region+"__fatJet2MassSD_JMRUp")
            hist_JMRUp.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMRUp")
            hists_sys.append(hist_JMRUp)
            hist_JMRDown =  inFile_this.Get(region+"__fatJet2MassSD_JMRDown")
            hist_JMRDown.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMRDown")
            hists_sys.append(hist_JMRDown)

            outFile.cd()
            hist_nominal.Write()
            hist_nominal_Blind =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassBlind"))
            nx = hist_nominal_Blind.FindBin(125.0)
            hist_nominal_Blind.SetBinContent(nx,0.0)
            hist_nominal_Blind.SetBinError(nx,0.0)
            hist_nominal_Blind.SetBinContent(nx+1,0.0)
            hist_nominal_Blind.SetBinError(nx+1,0.0)
            hist_nominal_Blind.SetBinContent(nx-1,0.0)
            hist_nominal_Blind.SetBinError(nx-1,0.0)
            hist_nominal_Blind.Write()
            
            hist_nominal_SR =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR"))
            
            #make a separate histogram for the signal region
            if region != "FailSRv8p2" :
                #data other than in Fail region remove all
                if idx==0 :
                    for ibin in range(hist_nominal.GetNbinsX()):
                        hist_nominal_SR.SetBinContent(ibin+1,0)
                        hist_nominal_SR.SetBinError(ibin+1,0)
                #only SR for signal and background
                else:
                    for ibin in range(hist_nominal.GetNbinsX()):
                        if (ibin<nx-2) or (ibin>nx):
                            hist_nominal_SR.SetBinContent(ibin+1,0)
                            hist_nominal_SR.SetBinError(ibin+1,0)
                        
                #hist_nominal_SR = r.TH1F(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR"), hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR"), 3, hist_nominal.GetBinLowEdge(nx-1), hist_nominal.GetBinLowEdge(nx+1) )
                #hist_nominal_SR.SetBinContent(2,hist_nominal.GetBinContent(nx))
                #hist_nominal_SR.SetBinError(2,hist_nominal.GetBinError(nx))
                #hist_nominal_SR.SetBinContent(3,hist_nominal.GetBinContent(nx+1))
                #hist_nominal_SR.SetBinError(3,hist_nominal.GetBinError(nx+1))
                #hist_nominal_SR.SetBinContent(1,hist_nominal.GetBinContent(nx-1))
                #hist_nominal_SR.SetBinError(1,hist_nominal.GetBinError(nx-1))
                
            hist_nominal_SR.Write() 
            
            for  hist in hists_sys:
                hist.Write()

                hist_Blind =  hist.Clone(hist.GetName().replace("histJet2Mass", "histJet2MassBlind"))
                nx = hist_Blind.FindBin(125.0)
                hist_Blind.SetBinContent(nx,0.0)
                hist_Blind.SetBinError(nx,0.0)
                hist_Blind.SetBinContent(nx+1,0.0)
                hist_Blind.SetBinError(nx+1,0.0)
                hist_Blind.SetBinContent(nx-1,0.0)
                hist_Blind.SetBinError(nx-1,0.0)
                hist_Blind.Write()
                
                hist_SR =  hist.Clone(hist.GetName().replace("histJet2Mass", "histJet2MassSR"))
                
                if region != "FailSRv8p2" :
                    #data other than in Fail region remove all
                    if idx==0 :
                        for ibin in range(hist.GetNbinsX()):
                            hist_SR.SetBinContent(ibin+1,0)
                            hist_SR.SetBinError(ibin+1,0)
                    #only SR for signal and background
                    else:
                        for ibin in range(hist.GetNbinsX()):
                            if (ibin<nx-2) or (ibin>nx):
                                hist_SR.SetBinContent(ibin+1,0)
                                hist_SR.SetBinError(ibin+1,0)
                            
                hist_SR.Write() 

        inFile_this.Close()
    outFile.Close()

