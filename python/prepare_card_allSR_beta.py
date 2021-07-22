import ROOT  as r
import sys
import math
import numpy as np

if __name__ == "__main__":
    tag = sys.argv[1]
    tag_trig = sys.argv[2]
    tag_ttbar_sys = sys.argv[3]
    vbdt = sys.argv[4]

    #proc  =  ["Data", "QCD", "TTJets", "VH",  "ttH", "others",  "HH", "HHVBFSM"]
    #proc_file  = ["data", "QCDggHVBF", "ttbar", "VH",  "ttH", "others",  "HHc1", "HHVBFSM"]
    
    proc_file  =  ["data", "QCDggHVBF", "ttbar", "VH",  "ttH", "others",  "HHc1", "HHc0",  "HHc5", "HHc2p45", "HHVBFSM", "HHVBF0p511", "HHVBF1p511", "HHVBF110", "HHVBF112", "HHVBF121", "HHVBF101"]
    proc  = ["Data", "QCD", "TTJets", "VH",  "ttH", "others",  "ggHH_kl_1_kt_1_boost4b", "ggHH_kl_0_kt_1_boost4b",  
                  "ggHH_kl_5_kt_1_boost4b", "ggHH_kl_2p45_kt_1_boost4b", "qqHH_CV_1_C2V_1_kl_1_boost4b",
                  "qqHH_CV_0p5_C2V_1_kl_1_boost4b", "qqHH_CV_1p5_C2V_1_kl_1_boost4b", "qqHH_CV_1_C2V_1_kl_0_boost4b", 
                  "qqHH_CV_1_C2V_1_kl_2_boost4b", "qqHH_CV_1_C2V_2_kl_1_boost4b", "qqHH_CV_1_C2V_0_kl_1_boost4b"]

    #source of weight systematics 
    systs_weight = ["pileupWeight","PNetShape","ttJetsCorr","BDT"+vbdt+"Shape","triggerEffSF","PNetHbbScaleFactors"]
    
    #source of shape systematics 
    systs_shape = ["JES","JMS","JMR"]
    
    outName = "HHTo4BPlots_Run2_BDT"+vbdt+".root"
    
    ttbar_weight_syst_index_start = 0 
    if "ttbar" in tag:
        outName = "HHTo4BPlots_Run2_ttbarSkim_BDT"+vbdt+".root"
    outFile =  r.TFile(outName, "recreate")

    #spectial treatment of nominal shape and JES/JMS/JMR shape systematics for the ttbar bkg
    inFile_this_ttbar_nominal = r.TFile("../hists/"+tag+"_nominal/combine/ttbar.root",  "READ")
    inFile_this_ttbar_loose = r.TFile("../hists/"+tag_ttbar_sys+"_nominal_nosys/combine/ttbar.root",  "READ")
    inFile_this_ttbar_loose_JES_Up = r.TFile("../hists/"+tag_ttbar_sys+"_JES_Up/combine/ttbar.root",  "READ")
    inFile_this_ttbar_loose_JES_Down = r.TFile("../hists/"+tag_ttbar_sys+"_JES_Down/combine/ttbar.root",  "READ")
    inFile_this_ttbar_loose_JMR_Up = r.TFile("../hists/"+tag_ttbar_sys+"_JMR_Up/combine/ttbar.root",  "READ")
    inFile_this_ttbar_loose_JMR_Down = r.TFile("../hists/"+tag_ttbar_sys+"_JMR_Down/combine/ttbar.root",  "READ")
    inFile_this_ttbar_loose_JMS_Up = r.TFile("../hists/"+tag_ttbar_sys+"_JMS_Up/combine/ttbar.root",  "READ")
    inFile_this_ttbar_loose_JMS_Down = r.TFile("../hists/"+tag_ttbar_sys+"_JMS_Down/combine/ttbar.root",  "READ")        
    
    for idx in range(len(proc)):
        print("study process: ",proc_file[idx])
        
        #read histogram with nominal and weight syst
        inFile_this = r.TFile("../hists/"+tag+"_nominal/combine/"+proc_file[idx]+".root",  "READ")
        print("read file "+"../hists/"+tag+"_nominal/combine/"+proc_file[idx]+".root")
        
        #read histogram for shape syst (jet related syst)
        inFile_systs_shape = []
        for isysts_shape in systs_shape:
            inFile_systs_shape.append(r.TFile("../hists/"+tag+"_"+isysts_shape+"_Up/combine/"+proc_file[idx]+".root",  "READ"))
            inFile_systs_shape.append(r.TFile("../hists/"+tag+"_"+isysts_shape+"_Down/combine/"+proc_file[idx]+".root",  "READ"))
         
        #names of the analysis categories (3 BDT bins + 1 fail region for QCD fit)
        region_list = []
        if "ttbar" in tag:
            region_list =  ["TTBarCR", "TTBarCRTight"]
        else:
            region_list = ["SRv8p2Bin1", "SRv8p2Bin2",  "SRv8p2Bin3", "FailSRv8p2"] #, "FitCRv8p2", "FailFitCRv8p2"]
        
        #loop over analysis categories
        for region in region_list:
            inFile_this.cd()
            hist_nominal = inFile_this.Get(region+"__fatJet2MassSD")
            outBinName=region.replace("SR",  "").replace("Fail", "fail").replace(vbdt, "")
            hist_nominal.SetName("histJet2Mass_"+outBinName+"_"+proc[idx])

            print("yields nominal test 0 ", hist_nominal.Integral(), "histJet2Mass_"+outBinName+"_"+proc[idx])
            hists_sys = []
            
            #loop over weight systematics
            for sys in systs_weight:
                print("studying sys for proc:",sys,proc[idx])
                
                if proc[idx]=="Data" or proc[idx]=="QCD":
                    print("proc: ",proc[idx], "continue")
                    continue
                    
                if(sys == "PNetHbbScaleFactors"):
                    inFile2016_PNet = r.TFile("../hists/"+tag+"_nominal/2016/"+proc_file[idx]+".root",  "READ")
                    inFile2017_PNet = r.TFile("../hists/"+tag+"_nominal/2017/"+proc_file[idx]+".root",  "READ")
                    inFile2018_PNet = r.TFile("../hists/"+tag+"_nominal/2018/"+proc_file[idx]+".root",  "READ")
                  
                    hist_2016_PNet = inFile2016_PNet.Get(region+"__fatJet2MassSD")
                    hist_2017_PNet = inFile2017_PNet.Get(region+"__fatJet2MassSD")
                    hist_2018_PNet = inFile2018_PNet.Get(region+"__fatJet2MassSD")
                    hist_Up = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                    hist_Down = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                    
                    for ibin in range(hist_Up.GetNbinsX()):                      
                        tmp_bin_up_sq = 0.0
                        tmp_bin_down_sq = 0.0
                        #6 bins in pT and 5 bins in PNet score
                        for ii in range(6):
                            for jj in range(5):                                               
                                hist_Up_2016_PNet = inFile2016_PNet.Get(region+sys+"2016bin"+str(ii+1)+str(jj+1)+"Up__fatJet2MassSD")
                                hist_Up_2017_PNet = inFile2017_PNet.Get(region+sys+"2017bin"+str(ii+1)+str(jj+1)+"Up__fatJet2MassSD")
                                hist_Up_2018_PNet = inFile2018_PNet.Get(region+sys+"2018bin"+str(ii+1)+str(jj+1)+"Up__fatJet2MassSD")
                                hist_Down_2016_PNet = inFile2016_PNet.Get(region+sys+"2016bin"+str(ii+1)+str(jj+1)+"Down__fatJet2MassSD")
                                hist_Down_2017_PNet = inFile2017_PNet.Get(region+sys+"2017bin"+str(ii+1)+str(jj+1)+"Down__fatJet2MassSD")
                                hist_Down_2018_PNet = inFile2018_PNet.Get(region+sys+"2018bin"+str(ii+1)+str(jj+1)+"Down__fatJet2MassSD")
                                #print("debug PNet SF unc: ", ii, jj, hist_Up_2016.Integral(), hist_Up_2017.Integral(), hist_Up_2018.Integral())                             
                                up_all = 1.0*(hist_Up_2016_PNet.GetBinContent(ibin+1) - hist_2016_PNet.GetBinContent(ibin+1)) ** 2 
                                + 1.0*(hist_Up_2017_PNet.GetBinContent(ibin+1) - hist_2017_PNet.GetBinContent(ibin+1)) ** 2
                                + 1.0*(hist_Up_2018_PNet.GetBinContent(ibin+1) - hist_2018_PNet.GetBinContent(ibin+1)) ** 2                                
                                tmp_bin_up_sq += up_all
                                                               
                        hist_Up.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)+np.sqrt(tmp_bin_up_sq)) 
                        hist_Down.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)-np.sqrt(tmp_bin_up_sq))
                        
                    #print("debug PNet SF unc before append to hists_sys: ",hist_Up.Integral(), hist_Down.Integral(), hist_nominal.Integral(), hist_2016.Integral(), hist_2017.Integral(), hist_2018.Integral())
                    hists_sys.append(hist_Up)
                    hists_sys.append(hist_Down) 
                    
                #trigger syst
                elif(sys == "triggerEffSF"):
                    inFile2016_this = r.TFile("../hists/"+tag_trig+"_nominal/2016/"+proc_file[idx]+".root",  "READ")
                    inFile2017_this = r.TFile("../hists/"+tag_trig+"_nominal/2017/"+proc_file[idx]+".root",  "READ")
                    inFile2018_this = r.TFile("../hists/"+tag_trig+"_nominal/2018/"+proc_file[idx]+".root",  "READ")
                  
                    hist_2016 = inFile2016_this.Get(region+"__fatJet2MassSD")
                    hist_2017 = inFile2017_this.Get(region+"__fatJet2MassSD")
                    hist_2018 = inFile2018_this.Get(region+"__fatJet2MassSD")
                    hist_Up = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                    hist_Down = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                
                    for ibin in range(hist_Up.GetNbinsX()):                      
                        tmp_bin_up_sq = 0.0
                        #print("type",type(tmp_bin_up_sq))
                        tmp_bin_down_sq = 0.0
                        #10x12x4 = 480 each year has 480 bins
                        for itrig in range(480):
                            if((itrig+1)%12!=1 and (itrig+1)%12!=2):                             
                                hist_Up_2016 = inFile2016_this.Get(region+sys+"2016bin"+str(itrig+1)+"Up__fatJet2MassSD")
                                hist_Up_2017 = inFile2017_this.Get(region+sys+"2017bin"+str(itrig+1)+"Up__fatJet2MassSD")
                                hist_Up_2018 = inFile2018_this.Get(region+sys+"2018bin"+str(itrig+1)+"Up__fatJet2MassSD")
                                hist_Down_2016 = inFile2016_this.Get(region+sys+"2016bin"+str(itrig+1)+"Down__fatJet2MassSD")
                                hist_Down_2017 = inFile2017_this.Get(region+sys+"2017bin"+str(itrig+1)+"Down__fatJet2MassSD")
                                hist_Down_2018 = inFile2018_this.Get(region+sys+"2018bin"+str(itrig+1)+"Down__fatJet2MassSD")
                                                                
                                up_all = (hist_Up_2016.GetBinContent(ibin+1) - hist_2016.GetBinContent(ibin+1)) ** 2 
                                + (hist_Up_2017.GetBinContent(ibin+1) - hist_2017.GetBinContent(ibin+1)) ** 2
                                + (hist_Up_2018.GetBinContent(ibin+1) - hist_2018.GetBinContent(ibin+1)) ** 2                                
                                tmp_bin_up_sq += up_all
                                
                        hist_Up.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)+np.sqrt(tmp_bin_up_sq)) 
                        hist_Down.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)-np.sqrt(tmp_bin_up_sq))
                    hists_sys.append(hist_Up)
                    hists_sys.append(hist_Down)
                    
                #other weight sysetmatics like pileup weights unc
                else:
                    #up variation
                    hist_Up = inFile_this.Get(region+sys+"Up__fatJet2MassSD")
                    #print(region+sys+"Up__fatJet2MassSD")
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                    #check nan numbers
                    for ibin in range(hist_Up.GetNbinsX()):
                        if math.isnan(hist_Up.GetBinContent(ibin+1)):
                            print("bug isNaN")
                            #hist_Up.SetBinContent(ibin+1,0)
                            #hist_Up.SetBinError(ibin+1,0)     
                    hists_sys.append(hist_Up)
                    #down variation
                    hist_Down = inFile_this.Get(region+sys+"Down__fatJet2MassSD")
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                    for ibin in range(hist_Down.GetNbinsX()):
                        if math.isnan(hist_Down.GetBinContent(ibin+1)):
                            print("bug isNaN")
                            #hist_Down.SetBinContent(ibin+1,0)
                            #hist_Down.SetBinError(ibin+1,0) 
                    hists_sys.append(hist_Down)
                print("yields nominal test 0.1", hist_nominal.Integral())

            #loop over shape syst unc
            for isysts_shape in range(len(systs_shape)):
                if proc_file[idx]=="data" or proc_file[idx]=="QCDggHVBF":
                    continue
                print(systs_shape[isysts_shape], proc_file[idx])
                print("yields nominal test 1", hist_nominal.Integral())

                hist_Up =  inFile_systs_shape[isysts_shape*2].Get(region+"__fatJet2MassSD")
                hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Up")
                
                hist_Down =  inFile_systs_shape[isysts_shape*2+1].Get(region+"__fatJet2MassSD")                 
                hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Down")
                
                print("yields nominal", hist_nominal.Integral())
                
                if (proc_file[idx]=="ttbar") and (region=="SRv8p2Bin1"):
                    #hist_loose_nominal = inFile_this_ttbar_loose.Get(region+"__fatJet2MassSD")
                    print("ttbar yields nominal", hist_nominal.Integral())
                    ratio_yield_ttbar = hist_nominal.Integral()/inFile_this_ttbar_loose.Get(region+"PNetp0__fatJet2MassSD").Integral()
                    print("ratio for ttbar yield", ratio_yield_ttbar)
                    if(systs_shape[isysts_shape]=="JMS"):
                        hist_Up = inFile_this_ttbar_loose_JMS_Up.Get(region+"PNetp0__fatJet2MassSD")
                        hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Up")
                        hist_Down = inFile_this_ttbar_loose_JMS_Down.Get(region+"PNetp0__fatJet2MassSD")
                        hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Down")
                        print("debug for JMS hist_loose_Up: ", hist_Up.Integral())
                        print("debug for JMS hist_loose_Down: ", hist_Down.Integral())
                    elif(systs_shape[isysts_shape]=="JMR"):
                        hist_Up = inFile_this_ttbar_loose_JMR_Up.Get(region+"PNetp0__fatJet2MassSD")
                        hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Up")
                        hist_Down = inFile_this_ttbar_loose_JMR_Down.Get(region+"PNetp0__fatJet2MassSD")
                        hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Down")
                        print("debug for JMR hist_loose_Up: ", hist_Up.Integral())
                        print("debug for JMR hist_loose_Down: ", hist_Down.Integral())
                    elif(systs_shape[isysts_shape]=="JES"):
                        hist_Up = inFile_this_ttbar_loose_JES_Up.Get(region+"PNetp0__fatJet2MassSD")
                        hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Up")
                        hist_Down = inFile_this_ttbar_loose_JES_Down.Get(region+"PNetp0__fatJet2MassSD")
                        hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Down")
                        print("debug for JES hist_loose_Up: ", hist_Up.Integral())
                        print("debug for JES hist_loose_Down: ", hist_Down.Integral())
                    
                    if(isysts_shape==0):
                        hist_nominal_save = hist_nominal.Clone("nominal_ttbar_Bin1_clone")
                    
                    hist_nominal = inFile_this_ttbar_loose.Get(region+"PNetp0__fatJet2MassSD")
                    outBinName=region.replace("SR",  "").replace("Fail", "fail").replace(vbdt, "")
                    hist_nominal.SetName("histJet2Mass_"+outBinName+"_"+proc[idx])
                    hist_nominal.Scale(ratio_yield_ttbar)
                    hist_Up.Scale(ratio_yield_ttbar);
                    print("ttbar hist up integral ",hist_Up.Integral())
                    hist_Down.Scale(ratio_yield_ttbar); 
                    print("ttbar hist down integral ",hist_Down.Integral())
                    
                    #now apply the weight systematic uncertainty to the new nominal ttbar shape
                    if(isysts_shape==0):
                        for isys_weight in range(len(systs_weight)):
                            print("debug: apply the weight shape sys for ttbar new shape: ", systs_weight[isys_weight])
                            print(ttbar_weight_syst_index_start+2*isys_weight, " len(hists_sys): ",len(hists_sys))
                            hists_sys[ttbar_weight_syst_index_start+2*isys_weight].Divide(hist_nominal_save)
                            hists_sys[ttbar_weight_syst_index_start+2*isys_weight].Multiply(hist_nominal) 
                            for ibin in range(hist_nominal.GetNbinsX()):
			    	hists_sys[ttbar_weight_syst_index_start+2*isys_weight].SetBinError(ibin+1, hist_nominal.GetBinError(ibin+1))
                            hists_sys[ttbar_weight_syst_index_start+1+2*isys_weight].Divide(hist_nominal_save)
                            hists_sys[ttbar_weight_syst_index_start+1+2*isys_weight].Multiply(hist_nominal) 
                            for ibin in range(hist_nominal.GetNbinsX()):
                                hists_sys[ttbar_weight_syst_index_start+1+2*isys_weight].SetBinError(ibin+1, hist_nominal.GetBinError(ibin+1))
                
                hists_sys.append(hist_Up)
                hists_sys.append(hist_Down)
                                  

            outFile.cd()
            if (idx==0) and (region != "FailSRv8p2") and (region != "FailSRv24") :
                print("debug blind data: ",region, " hist name ", hist_nominal.GetName())
            	nx = hist_nominal.FindBin(125.0)
                hist_nominal.SetBinContent(nx,0.0)
                hist_nominal.SetBinError(nx,0.0)
                hist_nominal.SetBinContent(nx+1,0.0)
                hist_nominal.SetBinError(nx+1,0.0)
                hist_nominal.SetBinContent(nx-1,0.0)
                hist_nominal.SetBinError(nx-1,0.0)
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
            
            if((region == "FailSRv8p2") or (region == "FailSRv24")):
                hist_nominal_fit =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2Massfit"))
                hist_nominal_fit.Write() 
                for  hist in hists_sys:
                    hist_fit =  hist.Clone(hist.GetName().replace("histJet2Mass", "histJet2Massfit"))
                    hist_fit.Write()            

            #hist_nominal_SR = r.TH1F(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR"), hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR"), 3, hist_nominal.GetBinLowEdge(nx-1), hist_nominal.GetBinLowEdge(nx+1) )
            #hist_nominal_SR.SetBinContent(2,hist_nominal.GetBinContent(nx))
            #hist_nominal_SR.SetBinError(2,hist_nominal.GetBinError(nx))
            #hist_nominal_SR.SetBinContent(3,hist_nominal.GetBinContent(nx+1))
            #hist_nominal_SR.SetBinError(3,hist_nominal.GetBinError(nx+1))
            #hist_nominal_SR.SetBinContent(1,hist_nominal.GetBinContent(nx-1))
            #hist_nominal_SR.SetBinError(1,hist_nominal.GetBinError(nx-1))
            
            #make a separate histogram for the signal region
            if (idx==0) and (region != "FailSRv8p2") and (region != "FailSRv24") :
                print("debug 1:",idx,region)
                #data in SR is set to zero 
                for ibin in range(hist_nominal.GetNbinsX()):
                    hist_nominal_SR.SetBinContent(ibin+1,0)
                    hist_nominal_SR.SetBinError(ibin+1,0)
                 
            #elif (region != "FailSRv8p2") and (region != "FailSRv24"):
            else:
                print("debug 2 Nan:",idx,region)
                #SR MC only include Higgs peak region
                for ibin in range(hist_nominal.GetNbinsX()):
                    if ((ibin<nx-2) or (ibin>nx)):
                        hist_nominal_SR.SetBinContent(ibin+1,0)
                        hist_nominal_SR.SetBinError(ibin+1,0)
                
            hist_nominal_SR.Write() 
            
            ij=0
            for  hist in hists_sys:
                print("for  hist in hists_sys ij ",ij)
                hist.Write()
                ij=ij+1

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
                
                if (idx==0) and (region != "FailSRv8p2") and (region != "FailSRv24") :
                    print("debug syst 1:",idx,region)
                    for ibin in range(hist.GetNbinsX()):
                        hist_SR.SetBinContent(ibin+1,0)
                        hist_SR.SetBinError(ibin+1,0)
                #elif (region != "FailSRv8p2") and (region != "FailSRv24"):
                else:
                    print("debug syst 2:",idx,region)
                    for ibin in range(hist.GetNbinsX()):
                        if ((ibin<nx-2) or (ibin>nx)):
                            hist_SR.SetBinContent(ibin+1,0)
                            hist_SR.SetBinError(ibin+1,0)
                            
                hist_SR.Write() 

        inFile_this.Close()
    outFile.Close()

