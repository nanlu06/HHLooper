import ROOT  as r
import sys



if __name__ == "__main__":
    tag = sys.argv[1]
    vbdt = sys.argv[2]

    proc  =      ["Data", "QCD", "TTJets", "others"]#, "VH",  "ttH", "others",  "HH", "HHc0",  "HHc5", "HHc2p45"]
    proc_file  = ["data", "qcd", "ttbar", "others"]#, "VH",  "ttH", "others",  "HHc1", "HHc0",  "HHc5", "HHc2p45"]

    #systs = ["BDT"+vbdt+"MassShape", "ttJetsCorr", "BDT"+vbdt+"Shape", "PNetShape"]
    systs= ["pileupWeight", "ttJetsCorr", "BDTv8p2Shape", "PNetShape"]
    outName = "HHTo4BPlots_Run2_BDT"+vbdt+"_ttbar.root"
    if "ttbar" in tag:
        outName = "HHTo4BPlots_Run2_ttbarSkim_BDT"+vbdt+".root"
    outFile =  r.TFile(outName, "recreate")

    for idx in range(len(proc)):
        inFile_this = r.TFile("../hists/"+tag+"_wsys/combine/"+proc_file[idx]+".root",  "READ")
        print("read file "+"../hists/"+tag+"_wsys/combine/"+proc_file[idx]+".root")
        
        inFile_JERUp_this = r.TFile("../hists/"+tag+"_JER_Up/combine/"+proc_file[idx]+".root",  "READ")
        inFile_JERDown_this = r.TFile("../hists/"+tag+"_JER_Down/combine/"+proc_file[idx]+".root",  "READ")

        inFile_JMSUp_this = r.TFile("../hists/"+tag+"_JMS_Up/combine/"+proc_file[idx]+".root",  "READ")
        inFile_JMSDown_this = r.TFile("../hists/"+tag+"_JMS_Down/combine/"+proc_file[idx]+".root",  "READ")
        
        inFile_JMRUp_this = r.TFile("../hists/"+tag+"_JMR_Up/combine/"+proc_file[idx]+".root",  "READ")
        inFile_JMRDown_this = r.TFile("../hists/"+tag+"_JMR_Down/combine/"+proc_file[idx]+".root",  "READ")

        inFile_JESUp_this = r.TFile("../hists/"+tag+"_JES_Up/combine/"+proc_file[idx]+".root",  "READ")
        inFile_JESDown_this = r.TFile("../hists/"+tag+"_JES_Down/combine/"+proc_file[idx]+".root",  "READ")
        
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
            hist_nominal = inFile_this.Get(region+"__fatJet2MassRegressed")
            outBinName=region.replace("SR",  "").replace("Fail", "fail").replace(vbdt, "")
            hist_nominal.SetName("histJet2Mass_"+outBinName+"_"+proc[idx])

            hists_sys = []
            for sys in systs:
                if idx==0: 
                    continue
                print("hist_Up ",  region+sys+"Up__fatJet2MassRegressed")
                if idx==1:
                    hist_Up = inFile_this.Get(region+"__fatJet2MassRegressed")
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                    hist_Down = inFile_this.Get(region+"__fatJet2MassRegressed")
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                    
                else:
                    hist_Up = inFile_this.Get(region+sys+"Up__fatJet2MassRegressed")
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                    hist_Down = inFile_this.Get(region+sys+"Down__fatJet2MassRegressed")
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                
                hists_sys.append(hist_Up)               
                hists_sys.append(hist_Down)
      
            if idx!=0:

                hist_JERUp =  inFile_JERUp_this.Get(region+"__fatJet2MassRegressed")
                hist_JERUp.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JERUp")
                hists_sys.append(hist_JERUp)
                hist_JERDown =  inFile_JERDown_this.Get(region+"__fatJet2MassRegressed")
                hist_JERDown.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JERDown")
                hists_sys.append(hist_JERDown)

                hist_JMSUp =  inFile_JMSUp_this.Get(region+"__fatJet2MassRegressed")
                hist_JMSUp.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMSUp")
                hists_sys.append(hist_JMSUp)
                hist_JMSDown =  inFile_JMSDown_this.Get(region+"__fatJet2MassRegressed")
                hist_JMSDown.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMSDown")
                hists_sys.append(hist_JMSDown)

                hist_JMRUp =  inFile_JMRUp_this.Get(region+"__fatJet2MassRegressed")
                hist_JMRUp.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMRUp")
                hists_sys.append(hist_JMRUp)
                hist_JMRDown =  inFile_JMRDown_this.Get(region+"__fatJet2MassRegressed")
                hist_JMRDown.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JMRDown")
                hists_sys.append(hist_JMRDown)

                hist_JESUp =  inFile_JESUp_this.Get(region+"__fatJet2MassRegressed")
                hist_JESUp.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JESUp")
                hists_sys.append(hist_JESUp)
                hist_JESDown =  inFile_JESDown_this.Get(region+"__fatJet2MassRegressed")
                hist_JESDown.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_JESDown")
                hists_sys.append(hist_JESDown)

            outFile.cd()
            hist_nominal.Write()
            
            for  hist in hists_sys:
                hist.Write()

        inFile_this.Close()
    outFile.Close()

