import ROOT  as r
import sys
import math
import numpy as np

def checkbin(hist):
    for ibin in range(hist.GetNbinsX()):
        ibinv = hist.GetBinContent(ibin+1)
        if math.isnan(ibinv):
            print("bug isNaN")
        elif(ibinv==0 and ("Data" not in hist.GetName())):
            #print("empty bin, set to 0.000001")
            hist.SetBinContent(ibin+1, 0.000001)
            hist.SetBinError(ibin+1, 0.000001)
    return;            
                
if __name__ == "__main__":
    tag = sys.argv[1]
    tag_trig = sys.argv[2]
    vbdt = sys.argv[3]
    histdir = sys.argv[4]
    ttbar_bin1 = sys.argv[5]
    ttbar_bin1_up = "PNetp92"
    ttbar_bin1_down = "PNetp2"
    isblind = sys.argv[6]

    #obs = "__fatJet2MassSD"
    obs = "__fatJet2MassRegressed" #Regressed
    #print out information for debugging
    debug = False
    
    proc_file  =  ["data", "QCDggHVBF", "ttbar", "VH",  "ttH", "others",  "HHc1", "HHc0",  "HHc5", "HHc2p45", "HHVBFSM", "HHVBF0p511", "HHVBF1p511", "HHVBF110", "HHVBF112", "HHVBF121", "HHVBF101"]
    proc  = ["Data", "QCD", "TTJets", "VH",  "ttH", "others",  "ggHH_kl_1_kt_1_boost4b", "ggHH_kl_0_kt_1_boost4b",  
                  "ggHH_kl_5_kt_1_boost4b", "ggHH_kl_2p45_kt_1_boost4b", "qqHH_CV_1_C2V_1_kl_1_boost4b",
                  "qqHH_CV_0p5_C2V_1_kl_1_boost4b", "qqHH_CV_1p5_C2V_1_kl_1_boost4b", "qqHH_CV_1_C2V_1_kl_0_boost4b", 
                  "qqHH_CV_1_C2V_1_kl_2_boost4b", "qqHH_CV_1_C2V_2_kl_1_boost4b", "qqHH_CV_1_C2V_0_kl_1_boost4b"]

    #source of weight systematics name here should match that in the histogram name
    systs_weight = ["trigCorrHH2016", "trigCorrHH2017", "trigCorrHH2018", "pileupWeight","PNetShape","ttJetsCorr","BDT"+vbdt+"Shape", "triggerEffSF","PNetHbbScaleFactors","FSRPartonShower","ISRPartonShower","FSRPartonShower_Vjets","ISRPartonShower_Vjets","ggHHPDFacc","ggHHQCDacc","othersQCD","mHHTHunc"]
        
    #source of shape systematics 
    systs_shape = ["JER","JES_Abs","JES_Abs_2016","JES_Abs_2017","JES_Abs_2018","JES_BBEC1","JES_BBEC1_2016","JES_BBEC1_2017","JES_BBEC1_2018","JES_EC2","JES_EC2_2016","JES_EC2_2017","JES_EC2_2018","JES_FlavQCD","JES_HF","JES_HF_2016","JES_HF_2017","JES_HF_2018","JES_RelBal","JES_RelSample_2016","JES_RelSample_2017","JES_RelSample_2018","JMS","JMR","ttbarBin1Jet2PNetCut"] 
    
    outName = "HHTo4BPlots_Run2_BDT"+vbdt+tag+".root"
    
    outFile =  r.TFile(outName, "recreate")

    #ttbar yield: scale PNet9 yield to the actual yield with PNet>0.98
    inFile_this_ttbar_loose = r.TFile(histdir+tag+"_nominal/combine/ttbar.root",  "READ")
    ttbar_bin1_yield = inFile_this_ttbar_loose.Get("SRv8p2Bin1"+obs).Integral(2,18)
    print("ttbar yields nominal", ttbar_bin1_yield)
    print("ttbar yields loose nominal", inFile_this_ttbar_loose.Get("SRv8p2Bin1"+ttbar_bin1+obs).Integral(2,18))
    ratio_yield_ttbar = ttbar_bin1_yield/inFile_this_ttbar_loose.Get("SRv8p2Bin1"+ttbar_bin1+obs).Integral(2,18)
    print("ratio for ttbar yield", ratio_yield_ttbar)
                
    for idx in range(len(proc)):
        #if proc[idx]=="Data" or proc[idx]=="QCD":
        #        if(debug): print("proc: ",proc[idx], "continue")
        #        continue                    
        print("study process: ",proc_file[idx])
        
        #read histogram with nominal and weight syst
        inFile_this = r.TFile(histdir+tag+"_nominal/combine/"+proc_file[idx]+".root",  "READ")
        print("read file "+histdir+tag+"_nominal/combine/"+proc_file[idx]+".root")
        
        #read histogram for shape syst (jet related syst)
        inFile_systs_shape = []
        inFile2016_systs_shape = []
        inFile2017_systs_shape = []
        inFile2018_systs_shape = []
        for isysts_shape in systs_shape:
            if isysts_shape != "ttbarBin1Jet2PNetCut":
                if 'JES' in isysts_shape:
                    syst_name_up = isysts_shape[:3]+'Up'+isysts_shape[3:] 
                    syst_name_down =isysts_shape[:3]+'Down'+isysts_shape[3:]
                    #print(syst_name_up)
                    inFile_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_up+"/combine/"+proc_file[idx]+".root",  "READ"))
                    inFile_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_down+"/combine/"+proc_file[idx]+".root",  "READ"))
                    inFile2016_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_up+"/2016/"+proc_file[idx]+".root",  "READ"))
                    inFile2016_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_down+"/2016/"+proc_file[idx]+".root",  "READ"))
                    inFile2017_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_up+"/2017/"+proc_file[idx]+".root",  "READ"))
                    inFile2017_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_down+"/2017/"+proc_file[idx]+".root",  "READ"))
                    inFile2018_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_up+"/2018/"+proc_file[idx]+".root",  "READ"))
                    inFile2018_systs_shape.append(r.TFile(histdir+tag+"_"+syst_name_down+"/2018/"+proc_file[idx]+".root",  "READ"))
                else:
                    inFile_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Up/combine/"+proc_file[idx]+".root",  "READ"))
                    inFile_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Down/combine/"+proc_file[idx]+".root",  "READ"))
                    inFile2016_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Up/2016/"+proc_file[idx]+".root",  "READ"))
                    inFile2016_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Down/2016/"+proc_file[idx]+".root",  "READ"))
                    inFile2017_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Up/2017/"+proc_file[idx]+".root",  "READ"))
                    inFile2017_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Down/2017/"+proc_file[idx]+".root",  "READ"))
                    inFile2018_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Up/2018/"+proc_file[idx]+".root",  "READ"))
                    inFile2018_systs_shape.append(r.TFile(histdir+tag+"_"+isysts_shape+"_Down/2018/"+proc_file[idx]+".root",  "READ"))
            else:
                inFile_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/combine/"+proc_file[idx]+".root",  "READ"))
                inFile_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/combine/"+proc_file[idx]+".root",  "READ"))
                inFile2016_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/2016/"+proc_file[idx]+".root",  "READ"))
                inFile2016_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/2016/"+proc_file[idx]+".root",  "READ"))
                inFile2017_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/2017/"+proc_file[idx]+".root",  "READ"))
                inFile2017_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/2017/"+proc_file[idx]+".root",  "READ"))
                inFile2018_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/2018/"+proc_file[idx]+".root",  "READ"))
                inFile2018_systs_shape.append(r.TFile(histdir+tag+"_nominal_nosys/2018/"+proc_file[idx]+".root",  "READ"))
                
         
        #names of the analysis categories (3 BDT bins + 1 fail region for QCD fit)
        region_list = ["SRv8p2Bin1", "SRv8p2Bin2",  "SRv8p2Bin3", "FailSRv8p2"]
        
        #loop over analysis categories
        for iregion in region_list:
                
            inFile_this.cd()
            
            if (proc[idx]=='TTJets') and ('Bin1' in iregion):
                region = iregion+ttbar_bin1
            else:
                region = iregion
                
            print("hist name debug: ", region+obs)
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

                elif("PartonShower_Vjets" in sys):
                    sys_m = sys.replace("_Vjets","")
                    if ("others" in proc[idx]):
                        hist_Up = inFile_this.Get(region+sys_m+"Up"+obs)
                        hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                        hist_Down = inFile_this.Get(region+sys_m+"Down"+obs)
                        hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                    else:
                        #print("debug PartonShower:", sys,proc[idx])
                        hist_Up = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up")
                        hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up") 

                        hist_Down = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")               
                        hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down")
                        
                elif "PartonShower" in sys or "trigCorrHH" in sys:
                    if ("HH" in proc[idx]):
                        hist_Up = inFile_this.Get(region+sys+"Up"+obs)
                        hist_Down = inFile_this.Get(region+sys+"Down"+obs)
                        if "trigCorrHH" in sys:
                            downdiff = np.abs(hist_Down.Integral()-hist_nominal.Integral())
                            updiff = np.abs(hist_Up.Integral()-hist_nominal.Integral())
                            if (downdiff<updiff):
                                for ibin in range(hist_nominal.GetNbinsX()):
                                    hist_Down.SetBinContent(ibin+1, 2.*hist_nominal.GetBinContent(ibin+1)-hist_Up.GetBinContent(ibin+1))
                            else:
                                for ibin in range(hist_nominal.GetNbinsX()):
                                    hist_Up.SetBinContent(ibin+1, 2.*hist_nominal.GetBinContent(ibin+1)-hist_Down.GetBinContent(ibin+1))
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

                elif "othersQCD" in sys:
                    print("starting to cal othersQCD unc for ", proc[idx])
                    hist_Up = hist_nominal.Clone(region+sys+"Up"+obs)
                    hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Up") 

                    hist_Down = hist_nominal.Clone(region+sys+"Down"+obs)
                    hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"Down") 

                    if ("others" in proc[idx]): 

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
                
                elif "pileupWeight" in sys:
                    
                    hist_Up = inFile_this.Get(region+sys+"Up"+obs)
                    hist_Down = inFile_this.Get(region+sys+"Down"+obs)
                    
                    inFile2016_pileup = r.TFile(histdir+tag+"_nominal/2016/"+proc_file[idx]+".root",  "READ")
                    inFile2017_pileup = r.TFile(histdir+tag+"_nominal/2017/"+proc_file[idx]+".root",  "READ")
                    inFile2018_pileup = r.TFile(histdir+tag+"_nominal/2018/"+proc_file[idx]+".root",  "READ")
                  
                    hist_2016_pileup = inFile2016_pileup.Get(region+obs)
                    hist_2017_pileup = inFile2017_pileup.Get(region+obs)
                    hist_2018_pileup = inFile2018_pileup.Get(region+obs)
                    
                    hist_Up_2016_pileup = inFile2016_pileup.Get(region+sys+"Up"+obs)
                    hist_Up_2017_pileup = inFile2017_pileup.Get(region+sys+"Up"+obs)
                    hist_Up_2018_pileup = inFile2018_pileup.Get(region+sys+"Up"+obs)
                    hist_Down_2016_pileup = inFile2016_pileup.Get(region+sys+"Down"+obs)
                    hist_Down_2017_pileup = inFile2017_pileup.Get(region+sys+"Down"+obs)
                    hist_Down_2018_pileup = inFile2018_pileup.Get(region+sys+"Down"+obs)
                    
                    #print("full run 2: ", hist_nominal.Integral())
                    #print("2016: ", hist_2016_pileup.Integral(), " ", hist_Up_2016_pileup.Integral(), " ", hist_Down_2016_pileup.Integral())
                    #print("2017: ", hist_2017_pileup.Integral(), " ", hist_Up_2017_pileup.Integral(), " ", hist_Down_2017_pileup.Integral())
                    #print("2018: ", hist_2018_pileup.Integral(), " ", hist_Up_2018_pileup.Integral(), " ", hist_Down_2018_pileup.Integral())
                    
                    hist_Up_2016 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2016Up")
                    hist_Up_2016.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2016Up") 

                    hist_Down_2016 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2016Down")
                    hist_Down_2016.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2016Down") 
                    
                    hist_Up_2017 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2017Up")
                    hist_Up_2017.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2017Up") 

                    hist_Down_2017 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2017Down")
                    hist_Down_2017.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2017Down") 
                    
                    hist_Up_2018 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2018Up")
                    hist_Up_2018.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2018Up") 

                    hist_Down_2018 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2018Down")
                    hist_Down_2018.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+sys.replace(vbdt,"")+"2018Down") 

                    for ibin in range(hist_Up.GetNbinsX()): 
                        #print("hist_Up_2016_pileup.GetBinContent(ibin+1) bin",ibin,hist_Up_2016_pileup.GetBinContent(ibin+1))
                        hist_Up_2016.SetBinContent(ibin+1,hist_Up_2016_pileup.GetBinContent(ibin+1) + hist_2017_pileup.GetBinContent(ibin+1)+hist_2018_pileup.GetBinContent(ibin+1)) 
                        hist_Down_2016.SetBinContent(ibin+1,hist_Down_2016_pileup.GetBinContent(ibin+1) + hist_2017_pileup.GetBinContent(ibin+1)+hist_2018_pileup.GetBinContent(ibin+1)) 
                        
                        hist_Up_2017.SetBinContent(ibin+1,hist_Up_2017_pileup.GetBinContent(ibin+1) + hist_2016_pileup.GetBinContent(ibin+1)+hist_2018_pileup.GetBinContent(ibin+1)) 
                        hist_Down_2017.SetBinContent(ibin+1,hist_Down_2017_pileup.GetBinContent(ibin+1) + hist_2016_pileup.GetBinContent(ibin+1)+hist_2018_pileup.GetBinContent(ibin+1))
                        
                        hist_Up_2018.SetBinContent(ibin+1,hist_Up_2018_pileup.GetBinContent(ibin+1) + hist_2016_pileup.GetBinContent(ibin+1)+hist_2017_pileup.GetBinContent(ibin+1)) 
                        hist_Down_2018.SetBinContent(ibin+1,hist_Down_2018_pileup.GetBinContent(ibin+1) + hist_2016_pileup.GetBinContent(ibin+1)+hist_2017_pileup.GetBinContent(ibin+1)) 
                    
                    if proc[idx]=="TTJets" and ("SRv8p2Bin1" in region):
                        hist_Up_2016.Scale(ratio_yield_ttbar)
                        hist_Down_2016.Scale(ratio_yield_ttbar) 
                        hist_Up_2017.Scale(ratio_yield_ttbar)
                        hist_Down_2017.Scale(ratio_yield_ttbar)
                        hist_Up_2018.Scale(ratio_yield_ttbar)
                        hist_Down_2018.Scale(ratio_yield_ttbar)
                    
                    hists_sys.append(hist_Up_2016)
                    hists_sys.append(hist_Down_2016)
                    hists_sys.append(hist_Up_2017)
                    hists_sys.append(hist_Down_2017)
                    hists_sys.append(hist_Up_2018)
                    hists_sys.append(hist_Down_2018)
                        
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
            inFile2016 = r.TFile(histdir+tag+"_nominal_nosys/2016/"+proc_file[idx]+".root",  "READ")                                               
            inFile2017 = r.TFile(histdir+tag+"_nominal_nosys/2017/"+proc_file[idx]+".root",  "READ")                                               
            inFile2018 = r.TFile(histdir+tag+"_nominal_nosys/2018/"+proc_file[idx]+".root",  "READ")
                                  
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
                    
                    #split JMS, JMR into three years
                    if systs_shape[isysts_shape] == "JMR" or systs_shape[isysts_shape] == "JMS" or systs_shape[isysts_shape] == "JER":

                        print("splitted: ",systs_shape[isysts_shape])
                        
                        hist_JERUp_2016 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2016Up")
                        hist_JERUp_2016.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2016Up")
                
                        hist_JERDown_2016 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2016Down")
                        hist_JERDown_2016.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2016Down")
                    
                        hist_JERUp_2017 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2017Up")
                        hist_JERUp_2017.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2017Up") 

                        hist_JERDown_2017 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2017Down")
                        hist_JERDown_2017.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2017Down") 
                
                        hist_JERUp_2018 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2018Up")
                        hist_JERUp_2018.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2018Up") 

                        hist_JERDown_2018 = hist_nominal.Clone("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2018Down")
                        hist_JERDown_2018.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"2018Down") 
                                                
                        hist_y2016 = inFile2016.Get(region+obs)
                        hist_y2016.SetName("y2016")
                        hist_y2017 = inFile2017.Get(region+obs)
                        hist_y2017.SetName("y2017")
                        hist_y2018 = inFile2018.Get(region+obs)
                        hist_y2018.SetName("y2018")
                        hist_Up_y2016 = inFile2016_systs_shape[isysts_shape*2].Get(region+obs)
                        hist_Up_y2016.SetName("hist_Up_y2016")
                        hist_Up_y2017 = inFile2017_systs_shape[isysts_shape*2].Get(region+obs)
                        hist_Up_y2017.SetName("hist_Up_y2017")
                        hist_Up_y2018 = inFile2018_systs_shape[isysts_shape*2].Get(region+obs)
                        hist_Up_y2018.SetName("hist_Up_y2018")
                        hist_Down_y2016 =  inFile2016_systs_shape[isysts_shape*2+1].Get(region+obs) 
                        hist_Down_y2016.SetName("hist_Down_y2016")
                        hist_Down_y2017 =  inFile2017_systs_shape[isysts_shape*2+1].Get(region+obs) 
                        hist_Down_y2017.SetName("hist_Down_y2017")
                        hist_Down_y2018 =  inFile2018_systs_shape[isysts_shape*2+1].Get(region+obs) 
                        hist_Down_y2018.SetName("hist_Down_y2018")
     
                        print("full run 2: ", hist_nominal.Integral())
                        print("2016: ", hist_y2016.Integral(), " ", hist_Up_y2016.Integral(), " ", hist_Down_y2016.Integral())
                        print("2017: ", hist_y2017.Integral(), " ", hist_Up_y2017.Integral(), " ", hist_Down_y2017.Integral())
                        print("2018: ", hist_y2018.Integral(), " ", hist_Up_y2018.Integral(), " ", hist_Down_y2018.Integral())

                        for ibin in range(hist_nominal.GetNbinsX()): 
                            #print("hist_Up_2016.GetBinContent(ibin+1) bin",ibin,hist_Up_y2016.GetBinContent(ibin+1))
                            
                            hist_JERUp_2016.SetBinContent(ibin+1,hist_Up_y2016.GetBinContent(ibin+1) + hist_y2017.GetBinContent(ibin+1)+hist_y2018.GetBinContent(ibin+1)) 
                            hist_JERDown_2016.SetBinContent(ibin+1,hist_Down_y2016.GetBinContent(ibin+1) + hist_y2017.GetBinContent(ibin+1)+hist_y2018.GetBinContent(ibin+1)) 
                        
                            hist_JERUp_2017.SetBinContent(ibin+1,hist_Up_y2017.GetBinContent(ibin+1) + hist_y2016.GetBinContent(ibin+1)+hist_y2018.GetBinContent(ibin+1)) 
                            hist_JERDown_2017.SetBinContent(ibin+1,hist_Down_y2017.GetBinContent(ibin+1) + hist_y2016.GetBinContent(ibin+1)+hist_y2018.GetBinContent(ibin+1))
                        
                            hist_JERUp_2018.SetBinContent(ibin+1,hist_Up_y2018.GetBinContent(ibin+1) + hist_y2016.GetBinContent(ibin+1)+hist_y2017.GetBinContent(ibin+1)) 
                            hist_JERDown_2018.SetBinContent(ibin+1,hist_Down_y2018.GetBinContent(ibin+1) + hist_y2016.GetBinContent(ibin+1)+hist_y2017.GetBinContent(ibin+1))                             
                                
                        if proc[idx]=="TTJets" and ("SRv8p2Bin1" in region):
                    
                            hist_up_total = hist_JERUp_2016.Integral(2,18)#[50-220 GeV]
                            hist_JERUp_2016.Scale(ttbar_bin1_yield/hist_up_total)
                            hist_down_total = hist_JERDown_2016.Integral(2,18)
                            hist_JERDown_2016.Scale(ttbar_bin1_yield/hist_down_total) 
                            
                            hist_up_total = hist_JERUp_2017.Integral(2,18)
                            hist_JERUp_2017.Scale(ttbar_bin1_yield/hist_up_total)
                            hist_down_total = hist_JERDown_2017.Integral(2,18)
                            hist_JERDown_2017.Scale(ttbar_bin1_yield/hist_down_total)
                            
                            hist_up_total = hist_JERUp_2018.Integral(2,18)
                            hist_JERUp_2018.Scale(ttbar_bin1_yield/hist_up_total)
                            hist_down_total = hist_JERDown_2018.Integral(2,18)
                            hist_JERDown_2018.Scale(ttbar_bin1_yield/hist_down_total)
                    
                        print("splitted shape sys integral 2016 2017 2018")
                        print(hist_JERUp_2016.Integral(), hist_JERDown_2016.Integral())
                        print(hist_JERUp_2017.Integral(), hist_JERDown_2017.Integral())
                        print(hist_JERUp_2018.Integral(), hist_JERDown_2018.Integral())
                    
                        hists_sys.append(hist_JERUp_2016)
                        hists_sys.append(hist_JERDown_2016)
                        hists_sys.append(hist_JERUp_2017)
                        hists_sys.append(hist_JERDown_2017)
                        hists_sys.append(hist_JERUp_2018)
                        hists_sys.append(hist_JERDown_2018)
                        
                        #end of JMS JMR pileupWeight if condition
                    #end of else condition
                #end of shape syst loop
                
                hist_Up.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Up")
                hist_Down.SetName("histJet2Mass_"+outBinName+"_"+proc[idx]+"_"+systs_shape[isysts_shape]+"Down")
                
                if(debug): print("yields nominal", hist_nominal.Integral())
                
                if proc[idx]=="TTJets" and ("SRv8p2Bin1" in region):
                    
                    hist_up_total = hist_Up.Integral(2,18);#[50-220 GeV]
                    hist_Up.Scale(ttbar_bin1_yield/hist_up_total)
                    print("ttbar hist up integral ",hist_Up.Integral(2,18))
                    
                    hist_down_total = hist_Down.Integral(2,18)
                    hist_Down.Scale(ttbar_bin1_yield/hist_down_total)
                    print("ttbar hist down integral ",hist_Down.Integral(2,18))
                
                    hist_nominal_total = hist_nominal.Integral(2,18)
                    hist_nominal.Scale(ttbar_bin1_yield/hist_nominal_total)
                    
                hists_sys.append(hist_Up)
                hists_sys.append(hist_Down)
                                  

            outFile.cd()
            checkbin(hist_nominal)
            if (idx==0) and (region != "FailSRv8p2") and (region != "FailSRv24") and (isblind=="yes"):
                if(debug): print("debug blind data: ",region, " hist name ", hist_nominal.GetName())
                nx = hist_nominal.FindBin(125.0)
                hist_nominal.SetBinContent(nx,0.0)
                hist_nominal.SetBinError(nx,0.0)
                hist_nominal.SetBinContent(nx+1,0.0)
                hist_nominal.SetBinError(nx+1,0.0)
                hist_nominal.SetBinContent(nx-1,0.0)
                hist_nominal.SetBinError(nx-1,0.0)
            hist_nominal.Write()
            
            #_Blind histogram is the sideband region
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
            
            #_SR is just three bins around the 125 GeV. We don't use it for datacards though
            hist_nominal_SR =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR"))
            hist_nominal_SR.SetName(hist_nominal.GetName().replace("histJet2Mass", "histJet2MassSR")) 

            #_fit histJet2Massfit is a duplication of histJet2Mass for the fail region only. just a convention in the datacard, we can use histJet2Mass as well.
            if((region == "FailSRv8p2") or (region == "FailSRv24")):
                hist_nominal_fit =  hist_nominal.Clone(hist_nominal.GetName().replace("histJet2Mass", "histJet2Massfit"))
                checkbin(hist_nominal_fit)
                hist_nominal_fit.Write() 
                for  hist in hists_sys:
                    hist_fit =  hist.Clone(hist.GetName().replace("histJet2Mass", "histJet2Massfit"))
                    hist_fit.SetName(hist.GetName().replace("histJet2Mass", "histJet2Massfit"))
                    checkbin(hist_fit)      
                    hist_fit.Write()            
            
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
