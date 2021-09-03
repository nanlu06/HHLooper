import ROOT  as r
import sys
import math
import numpy as np

def checkbin(hist):
    #nx = hist.FindBin(125.0)

    for ibin in range(hist.GetNbinsX()):
        ibinv = hist.GetBinContent(ibin+1)
        if math.isnan(ibinv):
            print("bug isNaN")
        elif(ibinv==0 and ("Data" not in hist.GetName())): #  and (ibin<nx-1 or ibin>nx+1)):
            #print("empty bin, set to 0.000001")
            hist.SetBinContent(ibin+1, 0.000001)
    return;            
                
if __name__ == "__main__":
    tag = sys.argv[1]
    tag_trig = sys.argv[2]
    vbdt = sys.argv[3]
    histdir = sys.argv[4]
    ttbar_bin1 = sys.argv[5]
    ttbar_bin1_up = "PNetp95"
    ttbar_bin1_down = "PNetp2"

    obs = "__fatJet2MassSD"
    #obs = "__fatJet2MassRegressed" #Regressed
    #print out information for debugging
    debug = False
    
    proc_file  =  ["data", "QCDggHVBF", "ttbar", "VH",  "ttH", "others",  "HHc1", "HHc0",  "HHc5", "HHc2p45", "HHVBFSM", "HHVBF0p511", "HHVBF1p511", "HHVBF110", "HHVBF112", "HHVBF121", "HHVBF101"]
    proc  = ["Data", "QCD", "TTJets", "VH",  "ttH", "others",  "ggHH_kl_1_kt_1_boost4b", "ggHH_kl_0_kt_1_boost4b",  
                  "ggHH_kl_5_kt_1_boost4b", "ggHH_kl_2p45_kt_1_boost4b", "qqHH_CV_1_C2V_1_kl_1_boost4b",
                  "qqHH_CV_0p5_C2V_1_kl_1_boost4b", "qqHH_CV_1p5_C2V_1_kl_1_boost4b", "qqHH_CV_1_C2V_1_kl_0_boost4b", 
                  "qqHH_CV_1_C2V_1_kl_2_boost4b", "qqHH_CV_1_C2V_2_kl_1_boost4b", "qqHH_CV_1_C2V_0_kl_1_boost4b"]

    #source of weight systematics name here should match that in the histogram name
    systs_weight = ["pileupWeight","PNetShape","ttJetsCorr","BDT"+vbdt+"Shape", #"triggerEffSF",
    "PNetHbbScaleFactors","FSRPartonShower","ISRPartonShower","ggHHPDFacc","ggHHQCDacc"]
    
    #source of shape systematics 
    systs_shape = ["JER","JES","JMS","JMR","ttbarBin1Jet2PNetCut"]
    
    outName = "HHTo4BPlots_Run2_BDT"+vbdt+tag+".root"
    
    outFile =  r.TFile(outName, "recreate")

    #ttbar yield: scale PNet9 yield to the actual yield with PNet>0.98
    inFile_this_ttbar_loose = r.TFile(histdir+tag+"_nominal/combine/ttbar.root",  "READ")
    ttbar_bin1_yield = inFile_this_ttbar_loose.Get("SRv8p2Bin1"+obs).Integral()
    print("ttbar yields nominal", ttbar_bin1_yield)
    print("ttbar yields loose nominal", inFile_this_ttbar_loose.Get("SRv8p2Bin1"+ttbar_bin1+obs).Integral())
    ratio_yield_ttbar = ttbar_bin1_yield/inFile_this_ttbar_loose.Get("SRv8p2Bin1"+ttbar_bin1+obs).Integral()
    print("ratio for ttbar yield", ratio_yield_ttbar)
                
    for idx in range(len(proc)):
        print("study process: ",proc_file[idx])
        
        #read histogram with nominal and weight syst
        inFile_this = r.TFile(histdir+tag+"_nominal/combine/"+proc_file[idx]+".root",  "READ")
        print("read file "+histdir+tag+"_nominal/combine/"+proc_file[idx]+".root")
        
        #read histogram for shape syst (jet related syst)
        inFile_systs_shape = []
        for isysts_shape in systs_shape:
            if isysts_shape != "ttbarBin1Jet2PNetCut":
                inFile_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Up/combine/"+proc_file[idx]+".root",  "READ"))
                inFile_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Down/combine/"+proc_file[idx]+".root",  "READ"))
            else:
                inFile_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/combine/"+proc_file[idx]+".root",  "READ"))
                inFile_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/combine/"+proc_file[idx]+".root",  "READ"))
         
        #names of the analysis categories (3 BDT bins + 1 fail region for QCD fit)
        region_list = ["SRv8p2Bin1", "SRv8p2Bin2",  "SRv8p2Bin3", "FailSRv8p2"] #, "FitCRv8p2", "FailFitCRv8p2"]
        
        #loop over analysis categories
        for iregion in region_list:
                
            inFile_this.cd()
            
            if (proc[idx]=='TTJets') and ('Bin1' in iregion):
                region = iregion+ttbar_bin1
            else:
                region = iregion
                
            hist_nominal = inFile_this.Get(region+obs)
            outBinName=region.replace("SR",  "").replace("Fail", "fail").replace(vbdt, "").replace(ttbar_bin1,"")
            hist_nominal.SetName("histJet2Mass_"+outBinName+"_"+proc[idx])

            if(debug): print("yields nominal test 0 ", hist_nominal.Integral(), "histJet2Mass_"+outBinName+"_"+proc[idx])
            hists_sys = []
            
            #loop over weight systematics
            for sys in systs_weight:
                print("studying systs_weight for proc in region:",sys,proc[idx],region)
                
                if proc[idx]=="Data" or proc[idx]=="QCD":
                    if(debug): print("proc: ",proc[idx], "continue")
                    continue
                    
                if(sys == "PNetHbbScaleFactors"):
                    inFile2016_PNet = r.TFile(histdir+tag+"_nominal/2016/"+proc_file[idx]+".root",  "READ")
                    inFile2017_PNet = r.TFile(histdir+tag+"_nominal/2017/"+proc_file[idx]+".root",  "READ")
                    inFile2018_PNet = r.TFile(histdir+tag+"_nominal/2018/"+proc_file[idx]+".root",  "READ")
                  
                    hist_2016_PNet = inFile2016_PNet.Get(region+obs)
                    hist_2017_PNet = inFile2017_PNet.Get(region+obs)
                    hist_2018_PNet = inFile2018_PNet.Get(region+obs)
                    hist_Up = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up") 

                    hist_Down = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down") 

                    for ibin in range(hist_Up.GetNbinsX()):                      
                        tmp_bin_up_sq = 0.0
                        tmp_bin_down_sq = 0.0
                        #6 bins in pT and 5 bins in PNet score
                        for ii in range(6):
                            for jj in range(5):                                               
                                hist_Up_2016_PNet = inFile2016_PNet.Get(region+sys+"2016bin"+str(ii+1)+str(jj+1)+"Up"+obs)
                                hist_Up_2017_PNet = inFile2017_PNet.Get(region+sys+"2017bin"+str(ii+1)+str(jj+1)+"Up"+obs)
                                hist_Up_2018_PNet = inFile2018_PNet.Get(region+sys+"2018bin"+str(ii+1)+str(jj+1)+"Up"+obs)
                                hist_Down_2016_PNet = inFile2016_PNet.Get(region+sys+"2016bin"+str(ii+1)+str(jj+1)+"Down"+obs)
                                hist_Down_2017_PNet = inFile2017_PNet.Get(region+sys+"2017bin"+str(ii+1)+str(jj+1)+"Down"+obs)
                                hist_Down_2018_PNet = inFile2018_PNet.Get(region+sys+"2018bin"+str(ii+1)+str(jj+1)+"Down"+obs)
                                #print("debug PNet SF unc: ", ii, jj, hist_Up_2016.Integral(), hist_Up_2017.Integral(), hist_Up_2018.Integral())                             
                                up_all = 1.0*(hist_Up_2016_PNet.GetBinContent(ibin+1) - hist_2016_PNet.GetBinContent(ibin+1)) ** 2 
                                + 1.0*(hist_Up_2017_PNet.GetBinContent(ibin+1) - hist_2017_PNet.GetBinContent(ibin+1)) ** 2
                                + 1.0*(hist_Up_2018_PNet.GetBinContent(ibin+1) - hist_2018_PNet.GetBinContent(ibin+1)) ** 2                                
                                tmp_bin_up_sq += up_all
                                                            
                        hist_Up.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)+np.sqrt(tmp_bin_up_sq)) 
                        hist_Down.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)-np.sqrt(tmp_bin_up_sq))
                            
                    #hist_Up = hist_nominal.Clone(region+sys+"Up"+obs)
                    #hist_Down = hist_nominal.Clone(region+sys+"Down"+obs)                        
                    
                #trigger syst
                elif(sys == "triggerEffSF"):
                    inFile2016_this = r.TFile(histdir+tag_trig+"_nominal/2016/"+proc_file[idx]+".root",  "READ")
                    inFile2017_this = r.TFile(histdir+tag_trig+"_nominal/2017/"+proc_file[idx]+".root",  "READ")
                    inFile2018_this = r.TFile(histdir+tag_trig+"_nominal/2018/"+proc_file[idx]+".root",  "READ")
                  
                    hist_2016 = inFile2016_this.Get(region+obs)
                    hist_2017 = inFile2017_this.Get(region+obs)
                    hist_2018 = inFile2018_this.Get(region+obs)
                    hist_Up = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up") 

                    hist_Down = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down") 

                    for ibin in range(hist_Up.GetNbinsX()):                      
                        tmp_bin_up_sq = 0.0
                        #print("type",type(tmp_bin_up_sq))
                        tmp_bin_down_sq = 0.0
                        #10x12x4 = 480 each year has 480 bins
                        for itrig in range(480):
                            if((itrig+1)%12!=1 and (itrig+1)%12!=2):                             
                                hist_Up_2016 = inFile2016_this.Get(region+sys+"2016bin"+str(itrig+1)+"Up"+obs)
                                hist_Up_2017 = inFile2017_this.Get(region+sys+"2017bin"+str(itrig+1)+"Up"+obs)
                                hist_Up_2018 = inFile2018_this.Get(region+sys+"2018bin"+str(itrig+1)+"Up"+obs)
                                hist_Down_2016 = inFile2016_this.Get(region+sys+"2016bin"+str(itrig+1)+"Down"+obs)
                                hist_Down_2017 = inFile2017_this.Get(region+sys+"2017bin"+str(itrig+1)+"Down"+obs)
                                hist_Down_2018 = inFile2018_this.Get(region+sys+"2018bin"+str(itrig+1)+"Down"+obs)
                                                                
                                up_all = (hist_Up_2016.GetBinContent(ibin+1) - hist_2016.GetBinContent(ibin+1)) ** 2 
                                + (hist_Up_2017.GetBinContent(ibin+1) - hist_2017.GetBinContent(ibin+1)) ** 2
                                + (hist_Up_2018.GetBinContent(ibin+1) - hist_2018.GetBinContent(ibin+1)) ** 2                                
                                tmp_bin_up_sq += up_all
                                
                        hist_Up.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)+np.sqrt(tmp_bin_up_sq)) 
                        hist_Down.SetBinContent(ibin+1,hist_nominal.GetBinContent(ibin+1)-np.sqrt(tmp_bin_up_sq))
                    
                elif "PartonShower" in sys:
                    if "HH" in proc[idx]:
                        hist_Up = inFile_this.Get(region+sys+"Up"+obs)
                        hist_Down = inFile_this.Get(region+sys+"Down"+obs)   
                    else:
                        #print("debug PartonShower:", sys,proc[idx])
                        hist_Up = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                        hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up") 

                        hist_Down = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")               
                        hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down") 
  
                elif "ggHHPDFacc" in sys:            
                    print("starting to cal ggHHPDFacc unc for ", proc[idx])
                    hist_Up = hist_nominal.Clone(region+sys+"Up"+obs)
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up") 

                    hist_Down = hist_nominal.Clone(region+sys+"Down"+obs)
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down") 

                    if ("HH" in proc[idx]) and ("qqHH" not in proc[idx]):         
                        for ibin in range(hist_nominal.GetNbinsX()):
                            nom_bin_content = hist_nominal.GetBinContent(ibin)
                            sq_tmp = 0;
                            for ipdf in range(103):
                                sq_tmp = sq_tmp + pow(inFile_this.Get(region+"LHEPDFEigenv"+str(ipdf)+obs).GetBinContent(ibin)-hist_nominal.GetBinContent(ibin),2);
                                
                            hist_Up.SetBinContent(ibin, nom_bin_content+np.sqrt(sq_tmp));
                            hist_Down.SetBinContent(ibin, nom_bin_content-np.sqrt(sq_tmp));                  
                            
                elif "ggHHQCDacc" in sys:
                    print("starting to cal ggHHQCDacc unc for ", proc[idx])
                    hist_Up = hist_nominal.Clone(region+sys+"Up"+obs)
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up") 

                    hist_Down = hist_nominal.Clone(region+sys+"Down"+obs)
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down") 

                    if ("HH" in proc[idx]) and ("qqHH" not in proc[idx]): 

                        for ibin in range(hist_nominal.GetNbinsX()):
                            nom_bin_content = hist_nominal.GetBinContent(ibin)
                            up_bin_content = nom_bin_content
                            down_bin_content = nom_bin_content
                            #take the envelope for weights [0,1,3,4,5,8,9]
                            for iqcd in range(9):
                                if iqcd == 2 or iqcd == 6:
                                    continue
                                else:
                                    tmp = inFile_this.Get(region+"QCDscale"+str(iqcd)+obs).GetBinContent(ibin)
                                    if tmp > up_bin_content:
                                        up_bin_content = tmp
                                    elif tmp < down_bin_content:
                                        down_bin_content = tmp
                            hist_Up.SetBinContent(ibin, up_bin_content)
                            hist_Down.SetBinContent(ibin, down_bin_content)                       
                    
                #weight sys which are not PNet, trig eff SF, HH theory acceptance uncertainties
                else:   
                    hist_Up = inFile_this.Get(region+sys+"Up"+obs)
                    hist_Down = inFile_this.Get(region+sys+"Down"+obs)

                if proc[idx]=="TTJets" and ("SRv8p2Bin1" in region):
                    hist_Up.Scale(ratio_yield_ttbar)
                    hist_Down.Scale(ratio_yield_ttbar) 
                        
                hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                hists_sys.append(hist_Up)
                hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down") 
                hists_sys.append(hist_Down)

            #loop over shape syst unc
            for isysts_shape in range(len(systs_shape)):
                if proc_file[idx]=="data" or proc_file[idx]=="QCDggHVBF":
                    continue
                if(debug): print("studying systs_shape for proc:",systs_shape[isysts_shape], proc_file[idx])
                    
                #ttbar bin 1 PNet cut uncertainty for the ttbar Jet 2 mass shape
                if (proc[idx]=='TTJets') and ('Bin1' in region) and systs_shape[isysts_shape] == "ttbarBin1Jet2PNetCut": 
                    print("debug ttbarBin1Jet2PNetCut")
                    hist_Up =  inFile_systs_shape[isysts_shape*2].Get(iregion+ttbar_bin1_up+obs)
                    hist_Down =  inFile_systs_shape[isysts_shape*2+1].Get(iregion+ttbar_bin1_down+obs)     
                else:
                    hist_Up =  inFile_systs_shape[isysts_shape*2].Get(region+obs)
                    hist_Down =  inFile_systs_shape[isysts_shape*2+1].Get(region+obs) 
                    
                hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Up")
                hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Down")
                
                if(debug): print("yields nominal", hist_nominal.Integral())
                
                if proc[idx]=="TTJets" and ("SRv8p2Bin1" in region):
                    
                    hist_up_total = hist_Up.Integral();
                    hist_Up.Scale(ttbar_bin1_yield/hist_up_total);
                    print("ttbar hist up integral ",hist_Up.Integral())
                    
                    hist_down_total = hist_Down.Integral();
                    hist_Down.Scale(ttbar_bin1_yield/hist_down_total); 
                    print("ttbar hist down integral ",hist_Down.Integral())
                
                    hist_nominal_total = hist_nominal.Integral();
                    hist_nominal.Scale(ttbar_bin1_yield/hist_nominal_total)
                    
                hists_sys.append(hist_Up)
                hists_sys.append(hist_Down)
                                  

            outFile.cd()
            checkbin(hist_nominal)
            if (idx==0) and (region != "FailSRv8p2") and (region != "FailSRv24") :
                if(debug): print("debug blind data: ",region, " hist name ", hist_nominal.GetName())
                nx = hist_nominal.FindBin(125.0)
                hist_nominal.SetBinContent(nx,0.0)
                hist_nominal.SetBinError(nx,0.0)
                hist_nominal.SetBinContent(nx+1,0.0)
                hist_nominal.SetBinError(nx+1,0.0)
                hist_nominal.SetBinContent(nx-1,0.0)
                hist_nominal.SetBinError(nx-1,0.0)
            hist_nominal.Write()
            hist_nominal_Blind =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassBlind"))
            checkbin(hist_nominal_Blind)
            nx = hist_nominal_Blind.FindBin(125.0)
            hist_nominal_Blind.SetBinContent(nx,0.0)
            hist_nominal_Blind.SetBinError(nx,0.0)
            hist_nominal_Blind.SetBinContent(nx+1,0.0)
            hist_nominal_Blind.SetBinError(nx+1,0.0)
            hist_nominal_Blind.SetBinContent(nx-1,0.0)
            hist_nominal_Blind.SetBinError(nx-1,0.0)
            hist_nominal_Blind.Write()
            
            hist_nominal_SR =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR"))
            hist_nominal_SR.SetName(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR")) 

            if((region == "FailSRv8p2") or (region == "FailSRv24")):
                hist_nominal_fit =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2Massfit"))
                checkbin(hist_nominal_fit)
                hist_nominal_fit.Write() 
                for  hist in hists_sys:
                    hist_fit =  hist.Clone(hist.GetName().replace("histJet2Mass", "histJet2Massfit"))
                    hist_fit.SetName(hist.GetName().replace("histJet2Mass", "histJet2Massfit"))
                    checkbin(hist_fit)      
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
                if(debug): print("debug 1:",idx,region)
                #data in SR is set to zero 
                for ibin in range(hist_nominal.GetNbinsX()):
                    hist_nominal_SR.SetBinContent(ibin+1,0)
                    hist_nominal_SR.SetBinError(ibin+1,0)
                 
            #elif (region != "FailSRv8p2") and (region != "FailSRv24"):
            else:
                #SR MC only include Higgs peak region
                for ibin in range(hist_nominal.GetNbinsX()):
                    if ((ibin<nx-2) or (ibin>nx)):
                        hist_nominal_SR.SetBinContent(ibin+1,0)
                        hist_nominal_SR.SetBinError(ibin+1,0)
                               
            checkbin(hist_nominal_SR)
                            
            hist_nominal_SR.Write() 
            
            ij=0
            for  hist in hists_sys:
                if(debug): print("for  hist in hists_sys ij ",ij)
                checkbin(hist)

                hist.Write()
                ij=ij+1

                hist_Blind =  hist.Clone(hist.GetName().replace("histJet2Mass", "histJet2MassBlind"))
                hist_Blind.SetName(hist.GetName().replace("histJet2Mass", "histJet2MassBlind"))
                checkbin(hist_Blind)
                nx = hist_Blind.FindBin(125.0)
                hist_Blind.SetBinContent(nx,0.0)
                hist_Blind.SetBinError(nx,0.0)
                hist_Blind.SetBinContent(nx+1,0.0)
                hist_Blind.SetBinError(nx+1,0.0)
                hist_Blind.SetBinContent(nx-1,0.0)
                hist_Blind.SetBinError(nx-1,0.0)
                hist_Blind.Write()
                
                hist_SR =  hist.Clone(hist.GetName().replace("histJet2Mass", "histJet2MassSR"))
                hist_SR.SetName(hist.GetName().replace("histJet2Mass", "histJet2MassSR"))

                if (idx==0) and (region != "FailSRv8p2") and (region != "FailSRv24") :
                    if(debug): print("debug syst 1:",idx,region)
                    for ibin in range(hist.GetNbinsX()):
                        hist_SR.SetBinContent(ibin+1,0)
                        hist_SR.SetBinError(ibin+1,0)
                #elif (region != "FailSRv8p2") and (region != "FailSRv24"):
                else:
                    if(debug): print("debug syst 2:",idx,region)
                    for ibin in range(hist.GetNbinsX()):
                        if ((ibin<nx-2) or (ibin>nx)):
                            hist_SR.SetBinContent(ibin+1,0)
                            hist_SR.SetBinError(ibin+1,0)
                checkbin(hist_SR)            
                hist_SR.Write() 

        inFile_this.Close()
    outFile.Close()

