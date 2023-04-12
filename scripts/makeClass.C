#include <stdlib.h>

void makeClass()
{
 TFile * file = new TFile("/eos/cms/store/group/phys_susy/razor/Run3Analysis/HH/HHTo4BNtupler/option5/nano/run3/combined/2022/ggHH_cHHH_1_TSG.root");

 TTree * tree = (TTree*)file->Get("tree");
 

 ofstream headerf;
 ofstream srcf;
 ofstream temp_h;
 ofstream temp_c1;
 ofstream temp_c2;
 headerf.open("hhtree.hh");
 srcf.open("hhtree.cc");
 temp_h.open("temp_h");
 temp_c1.open("temp_c1");
 temp_c2.open("temp_c2");
 //beginning of files
 headerf<<"#ifndef hhtree_HH"<<endl;
 headerf<<"#include <TBranch.h> \n#include <TTree.h> \n#include <TH1F.h> \n#include <TFile.h> \n#include <TBits.h> \n#include <vector> \n#include <unistd.h> \nusing namespace std; \n"<<endl;
 headerf<<"class hhtree{"<<endl;
 headerf<<" private:"<<endl;
 headerf<<" protected:"<<endl;
 headerf<<"   unsigned int index;"<<endl;

 srcf<<"#include \"hhtree.hh\""<<endl; 
 srcf<<"hhtree hh; \n"<<endl;
 srcf<<"void hhtree::Init(TTree *tree) \n{\n"<<endl;

 temp_c1<<"void hhtree::GetEntry(unsigned int idx)\n{\n"<<endl;
 temp_c1<<" index = idx;"<<endl;
 ////

 TSeqCollection *fullarray = tree->GetListOfAliases();
 bool have_aliases = true;
 if (!fullarray) {
   have_aliases = false;
   fullarray = tree->GetListOfBranches();
 }

 for (Int_t i = 0; i < fullarray->GetSize(); ++i) {
  TBranch *branch = 0;
  if (fullarray->At(i) == 0) continue;
  TString aliasname(fullarray->At(i)->GetName());
  if (have_aliases){
      branch = tree->GetBranch(tree->GetAlias(aliasname.Data()));
      aliasname = fullarray->At(i)->GetName();
  }
  else branch = (TBranch*)fullarray->At(i);
  if (branch == 0) continue;
  TString branchname(branch->GetName());
  TString branchtitle(branch->GetTitle());
  TString branchclass(branch->GetClassName());
  std::cout <<"branch: "<< branchtitle << std::endl;
 
  if (!branchname.BeginsWith("int") &&
        !branchname.BeginsWith("uint") &&
        !branchname.BeginsWith("ull") &&
        !branchname.BeginsWith("ullong") &&
        !branchname.BeginsWith("bool") &&
        !branchname.BeginsWith("float") &&
        !branchname.BeginsWith("double") &&
        !branchname.BeginsWith("uchars") &&
        !branchtitle.EndsWith("/F") &&
        !branchtitle.EndsWith("/I") &&
        !branchtitle.EndsWith("/i") &&
        !branchtitle.EndsWith("/l") &&
        !branchtitle.EndsWith("/L") &&
        !branchtitle.EndsWith("/O") &&
        !branchtitle.EndsWith("/D") &&
      !branchtitle.Contains("Trigger/flag bit (process: HLT)") &&
      !branchtitle.Contains("generator weight") &&
        !branchtitle.BeginsWith("TString") &&
        !branchtitle.BeginsWith("TBits") &&
        !branchclass.Contains("LorentzVector") &&
        !branchclass.Contains("int") &&
        !branchclass.Contains("uint") &&
        !branchclass.Contains("bool") &&
        !branchclass.Contains("float") &&
        !branchclass.Contains("double") &&
        !branchclass.Contains("string") &&
        !branchclass.Contains("TString") &&
        !branchclass.Contains("long long"))
        {
      continue;
        }
  TString type = "int";
  if(branchtitle.EndsWith("/F")) type = "float";
  if(branchtitle.EndsWith("/i")) type = "unsigned int";
  if(branchtitle.EndsWith("/l")) type = "unsigned long int";
  if(branchtitle.EndsWith("/L")) type = "long int";
  if(branchtitle.EndsWith("/O")) type = "bool";
  if(branchtitle.EndsWith("/Double")) type = "double";
  if(branchtitle.Contains("Trigger/flag bit (process: HLT)")) type = "bool";
  if(branchtitle.Contains("generator weight")) type = "float";
  srcf<<" "<<branchname<<"_branch = tree->GetBranch(\""<<branchname<<"\");"<<endl;

  if(branchtitle.Contains("["))
  {
    headerf<<"   "<<type<<" "<<branchname<<"_[10];"<<endl;
    headerf<<"   TBranch *"<<branchname<<"_branch;"<<endl;
    headerf<<"   bool "<<branchname<<"_isLoaded;"<<endl;
    temp_h<<"   const "<<type<<" *"<<branchname<<"();"<<endl;
    srcf<<" if("<<branchname<<"_branch) "<<branchname<<"_branch->SetAddress("<<branchname<<"_);"<<endl;
    temp_c2<<"const "<<type<<" * hhtree::"<<branchname<<"() \n{"<<endl;
  } 
  else
  {
    headerf<<"   "<<type<<" "<<branchname<<"_;"<<endl;
    headerf<<"   TBranch *"<<branchname<<"_branch;"<<endl;
    headerf<<"   bool "<<branchname<<"_isLoaded;"<<endl;
    temp_h<<"   const "<<type<<" &"<<branchname<<"();"<<endl;
    srcf<<" if("<<branchname<<"_branch) "<<branchname<<"_branch->SetAddress(&"<<branchname<<"_);"<<endl;
    temp_c2<<"const "<<type<<" &hhtree::"<<branchname<<"() \n{"<<endl;
  }
  temp_c1<<" "<<branchname<<"_isLoaded = false;"<<endl;

  temp_c2<<" if(not "<<branchname<<"_isLoaded)\n {"<<endl;
  temp_c2<<"   if("<<branchname<<"_branch != 0) "<<branchname<<"_branch->GetEntry(index);"<<endl;
  temp_c2<<"   else\n   {"<<endl;
  temp_c2<<"     printf(\"branch "<<branchname<<"_branch does not exist!\\n\");"<<endl;
  temp_c2<<"exit(1);\n   }"<<endl;
  temp_c2<<"   "<<branchname<<"_isLoaded = true;\n }"<<endl;
  temp_c2<<" return "<<branchname<<"_;\n}\n"<<endl;
  
 }
 temp_c1<<"}\n"<<endl;
 temp_h.close();
 temp_c1.close();
 temp_c2.close();

 headerf<<" \n "<<endl;
 headerf<<" public:"<<endl;
 headerf<<"   void Init(TTree * tree);"<<endl;
 headerf<<"   void GetEntry(unsigned int idx);"<<endl;
 headerf<<"   void LoadAllBranches();"<<endl;
 headerf.close();
 system("cat temp_h >> hhtree.hh");
 headerf.open("hhtree.hh", std::ofstream::app);
 headerf<<"};   \n"<<endl;
 headerf<<"#ifndef __CINT__   "<<endl;
 headerf<<"extern hhtree hh;   "<<endl;
 headerf<<"#endif   \n"<<endl;
 headerf<<"#endif   "<<endl;
 headerf.close();
 
 srcf<<"}\n"<<endl;
 srcf.close();
 system("cat temp_c1 >> hhtree.cc");
 system("cat temp_c2 >> hhtree.cc");

 system("rm temp*");

}
