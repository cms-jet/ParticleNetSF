#include "setTDRStyle.h"
#include <sstream>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include <vector>
#include "Rtypes.h"
#include "TColor.h"
#include "TVectorF.h"
#include <cstdlib>
#include <math.h>
#include "configuration.h"

TH2D *create2Dhisto(TString sample, TTree *tree,TString intLumi,TString cuts,
		    TString branchX,int binsX,float minX,float maxX,TString branchY,int binsY,float minY,float maxY,
		    bool useLog,TString name,bool data);
void makeTemplatesTop(TString path2file, TString era, TString cat, TString wpmin, TString wpmax, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cutmin="0==0", TString cutmax="0==0", TString suffix="none");
void makeMCHistosTop(TString name, TString path, std::vector<TString> processes, std::vector<TString> process_names, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
		     TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_);
void setTDRStyle();

double massScale(double mass, double scaleVal=1.05) { return scaleVal*mass; }

double massSmear(double mass, unsigned long lumi, unsigned long event, double sigma=0.1) {
  TRandom3 rnd((lumi << 10) + event);
  return rnd.Gaus(1, sigma)*mass; 
}


void make2DTemplates(TString sample, TString era, TString wpmin, TString wpmax) {
  
  conf::configuration(sample);
  
  if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") ) 
    { 
      std::cout << " In sample " << sample << " making 2D templates \n";
      makeTemplatesTop(sample,era,conf::category,wpmin,wpmax,conf::score_def,"200","1200");    
    }
}

void makeTemplatesTop(TString path2file, TString era, TString cat, TString wpmin, TString wpmax, TString score="fj_1_ParticleNetMD_XbbVsQCD", 
		      TString cutmin="0==0", TString cutmax="0==0", TString suffix="none") {

  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2(kTRUE);
  conf::configuration(path2file);

  // Processes, paths and lumi for each era
  vector<TString> processes     = conf::processes;
  vector<TString> process_names = conf::process_names;

  TString path;
  float intLumi;
  if (era == "2016") { path = conf::path_2016; intLumi= 19.52;  }
  //if (era == "2016") { path = conf::path_2016; intLumi= 16.81;  }
  if (era == "2017") { path = conf::path_2017; intLumi= 41.53; }
  if (era == "2018") { path = conf::path_2018; intLumi= 59.74; }
  ostringstream tmpLumi; tmpLumi << intLumi; TString lumi = tmpLumi.str();


  // Directory to store the templates   
  TString dirname1 = "templates2D";
  TString name0;
  if (score.Contains("ParticleNetMD"))      { name0 = "particlenetmd"; }
  else if (score.Contains("ParticleNet"))   { name0 = "particlenet"; }

  TString nameoutfile = conf::algo+"_tt1l_"+cat+"_"+wpmin+"to"+wpmax+"_"+era+"_"+cutmin+"to"+cutmax+"_templates";
  std::cout << " 2D templates name: " << nameoutfile << "\n";
  
  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }
 
  TFile *fout = new TFile("./"+dirname1+"/"+nameoutfile+".root","RECREATE");

  
  // Cuts and matching definition
  TString cut_ = "(passmetfilters && passMuTrig && fj_1_pt>="+cutmin+" && fj_1_pt<"+cutmax+")";

  TString c_base     = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && leptonicW_pt>150.) && ("+cut_+")";
  TString c_incl     = c_base+" && "+cut_;

  TString c_p3 = "( (fj_1_dr_T_Wq_max<0.8) && (fj_1_dr_T_b<0.8) )";
  TString c_p2 = "((fj_1_T_Wq_max_pdgId==0 && fj_1_dr_W_daus<0.8) || (fj_1_T_Wq_max_pdgId!=0 && fj_1_dr_T_b>=0.8 && fj_1_dr_T_Wq_max<0.8))";     
  TString c_p1 = "(!("+c_p3+" || "+c_p2+"))";


  std::vector<TString> cuts; cuts.clear();
  cuts.push_back(c_incl);
  cuts.push_back(c_incl+" && "+c_p3);
  cuts.push_back(c_incl+" && "+c_p2);
  cuts.push_back(c_incl+" && "+c_p1);

  // WP selection
  TString wp_val;

  TString c_p = "("+score+">"+wpmin+" && "+score+"<="+wpmax+")";
  TString c_f = "(!"+c_p+")";

  cuts.push_back(c_p);
  cuts.push_back(c_f);

  // Fit variables and axis ranges
  TString brX = conf::brX; int binsX = conf::binsX; float minX = conf::minX;  float maxX = conf::maxX;
  TString brY = conf::brY; int binsY = conf::binsY; float minY = conf::minY;  float maxY = conf::maxY;
  TString name = path2file+"_"+name0+"_"+cat+"_"+wpmin+"to"+wpmax+"_"+era;

  
  // Data histograms 
  TFile *f_data  = TFile::Open(path+"/data/singlemu_tree.root" , "READONLY");
  TTree *t_data  = (TTree*)f_data->Get("Events");
  TH2D *h_data_p = create2Dhisto(name,t_data,lumi,cuts[0]+" && "+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_p",true);
  TH2D *h_data_f = create2Dhisto(name,t_data,lumi,cuts[0]+" && "+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_f",true);
  fout->cd();
  h_data_p->Write("data_obs_pass");
  h_data_f->Write("data_obs_fail");
  
  // MC templates
  std::vector<TString> syst = conf::syst;
  TString name_, namesys_;
  for (int i0=0; i0<syst.size(); ++i0) {
    name_ = syst.at(i0); namesys_ = "";
    if (syst.at(i0)=="_")  { name_ = "nom"; namesys_ = "nom"; }
    if (syst.at(i0)=="pu") { namesys_ = "pu"; }
    
    if (syst.at(i0)=="_") { makeMCHistosTop(name,path,processes,process_names,name_,namesys_,lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout); }
    else {
      
      if(name_ == "lhescalemuf"){
          makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Up",lumi+"*(LHEScaleWeight[5]*LHEScaleWeightNorm[5])/(LHEScaleWeight[4]*LHEScaleWeightNorm[4])",cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
          makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Down",lumi+"*(LHEScaleWeight[3]*LHEScaleWeightNorm[3])/(LHEScaleWeight[4]*LHEScaleWeightNorm[4])",cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
      }
      else if(name_ == "lhescalemur"){
          makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Up",lumi+"*(LHEScaleWeight[7]*LHEScaleWeightNorm[7])/(LHEScaleWeight[4]*LHEScaleWeightNorm[4])",cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
          makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Down",lumi+"*LHEScaleWeight[1]*LHEScaleWeightNorm[1]/(LHEScaleWeight[4]*LHEScaleWeightNorm[4])",cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
      }                
      else{     
          makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Up",lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
          makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Down",lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
      }
    }
  }

  fout->Close();
  std::cout << "\n\n";
}


void makeMCHistosTop(TString name, TString path, std::vector<TString> processes, std::vector<TString> process_names, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
		     TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_) { 

  TString sys_type = "/";
  if ( sysType == "Up" ) { sys_type = "_up/"; } if ( sysType == "Down" ) { sys_type ="_down/"; }
		    
  TString sys_dir; 
  if      ( (sys == "nom") || (sys == "pu") || (sys == "jms") || (sys == "jmr") ) { sys_dir = "/mc_nom/"; }  
  else if ( (sys.Contains("lhe")) || (sys.Contains("ps")) ) { sys_dir = "/mc_sys/LHEWeight/"; }
  else                                                      { sys_dir = "/mc_sys/"+sys+sys_type; }

  if ( (sys == "nom") || (sys == "pu") ) { name = name+"_"+sys; }
  else                                   { name = name+"_"+sys+sysType; }

  TString name_b;
  if (sys == "nom")      { name_b = ""; }
  else if (sys == "pu")  { name_b = sysType; }
  else                   { name_b = sys+sysType; }
  
  std::vector<TH2D*>   h2ds;       h2ds.clear();
  std::vector<TString> h2ds_names; h2ds_names.clear();
  for (unsigned int i0=0; i0<processes.size(); ++i0) {
    
    TFile *f = TFile::Open(path+"/"+sys_dir+processes[i0]+"_tree.root","READONLY");
    TTree *t = (TTree*)f->Get("Events");
   
    if ( (process_names[i0] == "tt") || (process_names[i0] == "st") || (process_names[i0] == "ttv")) {
 
      TH2D *h_p3_p = create2Dhisto(name,t,wgts,cuts[1]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_p3_p",false); h2ds.push_back(h_p3_p); h2ds_names.push_back(process_names[i0]+"_p3_"+name_b+"_pass"); 
      TH2D *h_p2_p = create2Dhisto(name,t,wgts,cuts[2]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_p2_p",false); h2ds.push_back(h_p2_p); h2ds_names.push_back(process_names[i0]+"_p2_"+name_b+"_pass");
      TH2D *h_p1_p = create2Dhisto(name,t,wgts,cuts[3]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_p1_p",false); h2ds.push_back(h_p1_p); h2ds_names.push_back(process_names[i0]+"_p1_"+name_b+"_pass");
      
      TH2D *h_p3_f = create2Dhisto(name,t,wgts,cuts[1]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_p3_f",false); h2ds.push_back(h_p3_f); h2ds_names.push_back(process_names[i0]+"_p3_"+name_b+"_fail");
      TH2D *h_p2_f = create2Dhisto(name,t,wgts,cuts[2]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_p2_f",false); h2ds.push_back(h_p2_f); h2ds_names.push_back(process_names[i0]+"_p2_"+name_b+"_fail");
      TH2D *h_p1_f = create2Dhisto(name,t,wgts,cuts[3]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_p1_f",false); h2ds.push_back(h_p1_f); h2ds_names.push_back(process_names[i0]+"_p1_"+name_b+"_fail");
    }
    else {
      TH2D *h_p = create2Dhisto(name,t,wgts,cuts[0]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_p",false); h2ds.push_back(h_p); h2ds_names.push_back(process_names[i0]+"_"+name_b+"_pass");
      TH2D *h_f = create2Dhisto(name,t,wgts,cuts[0]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_f",false); h2ds.push_back(h_f); h2ds_names.push_back(process_names[i0]+"_"+name_b+"_fail");
    }
  } // end of looping over the trees
  
  
  // avoid zero bins in mc
  for (unsigned int i0=0; i0<h2ds.size(); ++i0) {
    for (unsigned int i1=0; i1<h2ds[i0]->GetNbinsX(); ++i1) {
      for (unsigned int i1y=0; i1y<h2ds[i0]->GetNbinsY(); ++i1y) {
        if (h2ds[i0]->GetBinContent(i1,i1y)<=0) { h2ds[i0]->SetBinContent(i1,i1y,0.001); h2ds[i0]->SetBinError(i1,i1y,0.001); }
      }
    }
  }

  // write histos to file
  f_->cd();
  for (unsigned int i0=0; i0<h2ds.size(); ++i0) { h2ds[i0]->Write(h2ds_names[i0]); }
  
} // end of makeMCHistosTop


TH2D *create2Dhisto(TString sample, TTree *tree,TString intLumi,TString cuts,
                    TString branchX,int binsX,float minX,float maxX,TString branchY,int binsY,float minY,float maxY,
                    bool useLog,TString name,bool data) {

  TH1::SetDefaultSumw2(kTRUE);
  
  TString puWgt;
  if (name.Contains("puUp"))        { puWgt = "puWeightUp"; }
  else if (name.Contains("puDown")) { puWgt = "puWeightDown"; }
  else                              { puWgt = "puWeight"; }


  TString genWgt = "xsecWeight*genWeight";

  // additinal weights
  TString wBosonWgt = "(w_qcdnlo_wgt_func(fj_1_pt)*w_ewknlo_wgt_func(fj_1_pt))"; //NLOQCD=1.35 for 2016
  TString zBosonWgt = "(z_qcdnlo_wgt_func(fj_1_pt)*z_ewknlo_wgt_func(fj_1_pt))"; //NLOQCD=1.45 for 2016
  std::cout << " Z W weights = " << zBosonWgt << " " << wBosonWgt << "\n";

  TString qcdWgt    = "1.";
  TString ttWgt     = "1."; if (sample.Contains("tt1L") || sample.Contains("ttbar1L") || sample.Contains("tt1l")) { ttWgt = "topptWeight"; };

  conf::configuration(sample);
  TString cut;
  if (data) { cut ="("+cuts+")"; } 
  else {
    if (name.Contains("qcd"))           { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+"*"+qcdWgt+")*("+cuts+")"; }
    else if (name.Contains("qcherwig")) { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+")*("+cuts+")"; }
    else if (name.Contains("wboson"))   { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+"*"+wBosonWgt+")*("+cuts+")"; }
    else if (name.Contains("zboson"))   { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+"*"+zBosonWgt+")*("+cuts+")"; }
    else if (name.Contains("tt"))       { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+"*"+ttWgt+")*("+cuts+")"; }
    else if (name.Contains("higgs"))    { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+")*("+cuts+")"; }
    else if (name.Contains("sm"))       { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+")*("+cuts+")"; }
    else                                { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+")*("+cuts+")"; }
  }
  
  std::cout << "\n";
  std::cout << "Process name: " << name << "\n";
  std::cout << "Cut = " << cut << "\n";
  std::cout << "\n";

  
  TH2D *hTemp = new TH2D(name,name,binsX,minX,maxX,binsY,minY,maxY); //hTemp->SetName(name);

  TString massScaleVal_ = "1.05"; if (name.Contains("Down")) { massScaleVal_ = "0.95"; }
  TString massSmearVal_ = "0.10"; if (name.Contains("Down")) { massSmearVal_ = "0."; }
  if (name.Contains("jms")) { std::cout << tree->Project(name,branchY+":(massScale("+branchX+","+massScaleVal_+"))",cut); }
  else if (name.Contains("jmr")) { std::cout << " In jmr \n"; tree->Project(name,branchY+":(massSmear("+branchX+",luminosityBlock,event,"+massSmearVal_+"))",cut); }
  else                           { std::cout << " In else \n"; tree->Project(name,branchY+":"+branchX,cut); }
  
  for (unsigned int i0x=0; i0x<hTemp->GetNbinsX(); ++i0x) {
    double error =0.; double integral = hTemp->IntegralAndError(i0x+1,i0x+1,hTemp->GetNbinsY(),hTemp->GetNbinsY()+1,error);
    hTemp->SetBinContent(i0x+1,hTemp->GetNbinsY(),integral);
    hTemp->SetBinError(i0x+1,hTemp->GetNbinsY(),error);
  } 
  
  return hTemp;
}
