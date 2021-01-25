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

void makeHistoHerwig(TString name_, TString path_, TString wgt_, std::vector<TString> cuts_, TString c_jet_, TString sys_, TFile *f_);
void makeHistoLHEPDF(TString name_, TString path_, TString wgt_, std::vector<TString> cuts_, TString c_jet_, TString sys_, TString var_, TFile *f_);
void makeHisto(TString path_, std::vector<TString> names_, std::vector<TString> wgt_, std::vector<float> scale_p, std::vector<float> scale_f, std::vector<TString> cuts_, 
	       std::vector<TString> c_p, std::vector<TString> c_f, TString sys_, TString var_, int bins_, float xmin_, float xmax_, TFile *fout_p, TFile *fout_f);
TH1D *create1Dhisto(TString sample, TTree *tree,TString intLumi,TString cuts,TString branch,int bins,float xmin,float xmax,
		    bool useLog,int color, int style,TString name,bool norm, bool data);
TH2D *create2Dhisto(TString sample, TTree *tree,TString intLumi,TString cuts,
		    TString branchX,int binsX,float minX,float maxX,TString branchY,int binsY,float minY,float maxY,
		    bool useLog,TString name,bool data);
void makeTemplates(TString path2file, TString era, TString cat, TString wp, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cutmin="0==0", TString cutmax="0==0", TString suffix="none");
void makeTemplatesTop(TString path2file, TString era, TString cat, TString wpmin, TString wpmax, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cutmin="0==0", TString cutmax="0==0", TString suffix="none");
TH2D *correct_2dext_norm(TH2D *h_nom, TH2D *h_ext);
void makeMCHistos(TString name, TString path, TString path2file, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
		  TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_);
void makeMCHistosTop(TString name, TString path, std::vector<TString> processes, std::vector<TString> process_names, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
		     TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_);
void makeMCHistosLHEPDFTop(TString name, TString path, std::vector<TString> processes, std::vector<TString> process_names, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
			   TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_);
std::vector<TH2D*> getLHPDFRmsAndSigma(std::vector<TH2D*> input);
void setTDRStyle();

// corrections to MC
/*
TH1D *h_w_ewknlo_wgt = nullptr;
float w_ewknlo_wgt_func(float x){ if (!h_w_ewknlo_wgt) return 1; return h_w_ewknlo_wgt->GetBinContent(h_w_ewknlo_wgt->FindFixBin(x)); }
TH1D *h_z_ewknlo_wgt = nullptr; 
float z_ewknlo_wgt_func(float x){ if (!h_z_ewknlo_wgt) return 1; return h_z_ewknlo_wgt->GetBinContent(h_z_ewknlo_wgt->FindFixBin(x)); }
*/

TH1D *h_z_ewknlo_wgt = nullptr; 
float z_ewknlo_wgt_func(float x){ if (!h_z_ewknlo_wgt) return 1; 
  //std::cout << "z_ewknlo = " << h_z_ewknlo_wgt->GetBinContent(h_z_ewknlo_wgt->FindFixBin(x)) << "\n";
  return h_z_ewknlo_wgt->GetBinContent(h_z_ewknlo_wgt->FindFixBin(x)); }

TH1D *h_z_qcdnlo_wgt = nullptr; 
float z_qcdnlo_wgt_func(float x){ if (!h_z_qcdnlo_wgt) return 1; 
  //std::cout << "z_qcdnlo = " << h_z_qcdnlo_wgt->GetBinContent(h_z_qcdnlo_wgt->FindFixBin(x)) << "\n";
  return h_z_qcdnlo_wgt->GetBinContent(h_z_qcdnlo_wgt->FindFixBin(x)); }

TH1D *h_w_ewknlo_wgt = nullptr; 
float w_ewknlo_wgt_func(float x){ if (!h_w_ewknlo_wgt) return 1; 
  //std::cout << "w_ewknlo = " << h_w_ewknlo_wgt->GetBinContent(h_w_ewknlo_wgt->FindFixBin(x)) << "\n";
  return h_w_ewknlo_wgt->GetBinContent(h_w_ewknlo_wgt->FindFixBin(x)); }

TH1D *h_w_qcdnlo_wgt = nullptr; 
float w_qcdnlo_wgt_func(float x){ if (!h_w_qcdnlo_wgt) return 1; 
  //std::cout << "w_qcdnlo = " << h_w_qcdnlo_wgt->GetBinContent(h_w_qcdnlo_wgt->FindFixBin(x)) << "\n";
  return h_w_qcdnlo_wgt->GetBinContent(h_w_qcdnlo_wgt->FindFixBin(x)); }

TH1D *h_qcd_incl_mass_wgt = nullptr;
float qcd_incl_mass_wgt(float x){ if (!h_qcd_incl_mass_wgt) return 1; return h_qcd_incl_mass_wgt->GetBinContent(h_qcd_incl_mass_wgt->FindFixBin(x)); }



// DDT
TH2D *hrwgt_dak8ddt_w_0p05 = nullptr;
double rewgtfuncDAK8DDT_w_0p05(double rho, double pt){
  if (!hrwgt_dak8ddt_w_0p05) return 1;
  return hrwgt_dak8ddt_w_0p05->GetBinContent(hrwgt_dak8ddt_w_0p05->GetXaxis()->FindBin(rho),hrwgt_dak8ddt_w_0p05->GetYaxis()->FindBin(pt));
}

TH2D *hrwgt_dak8ddt_w_0p10 = nullptr;
double rewgtfuncDAK8DDT_w_0p10(double rho, double pt){
  if (!hrwgt_dak8ddt_w_0p10) return 1;
  return hrwgt_dak8ddt_w_0p10->GetBinContent(hrwgt_dak8ddt_w_0p10->GetXaxis()->FindBin(rho),hrwgt_dak8ddt_w_0p10->GetYaxis()->FindBin(pt));
}

TH2D *hrwgt_dak8ddt_w_0p20 = nullptr;
double rewgtfuncDAK8DDT_w_0p20(double rho, double pt){
  if (!hrwgt_dak8ddt_w_0p20) return 1;
  return hrwgt_dak8ddt_w_0p20->GetBinContent(hrwgt_dak8ddt_w_0p20->GetXaxis()->FindBin(rho),hrwgt_dak8ddt_w_0p20->GetYaxis()->FindBin(pt));
}

TH2D *hrwgt_dak8ddt_w_0p30 = nullptr;
double rewgtfuncDAK8DDT_w_0p30(double rho, double pt){
  if (!hrwgt_dak8ddt_w_0p30) return 1;
  return hrwgt_dak8ddt_w_0p30->GetBinContent(hrwgt_dak8ddt_w_0p30->GetXaxis()->FindBin(rho),hrwgt_dak8ddt_w_0p30->GetYaxis()->FindBin(pt));
}

TH2D *hrwgt_dak8ddt_w_0p50 = nullptr;
double rewgtfuncDAK8DDT_w_0p50(double rho, double pt){
  if (!hrwgt_dak8ddt_w_0p50) return 1;
  return hrwgt_dak8ddt_w_0p50->GetBinContent(hrwgt_dak8ddt_w_0p50->GetXaxis()->FindBin(rho),hrwgt_dak8ddt_w_0p50->GetYaxis()->FindBin(pt));
}


double massScale(double mass, double scaleVal=1.05) { return scaleVal*mass; }

double massSmear(double mass, unsigned long lumi, unsigned long event, double sigma=0.1) {
  TRandom3 rnd((lumi << 10) + event);
  return rnd.Gaus(1, sigma)*mass; 
}


void make2DTemplates(TString sample, TString era, TString wpmin, TString wpmax) {
  
  conf::configuration(sample);
  
  if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") ) { 
    std::cout << " In sample " << sample << " making 2D templates \n";
    makeTemplatesTop(sample,era,conf::category,wpmin,wpmax,conf::score_def,"200","1200");    

  }

  if (sample == "zqq") { 
    std::cout << " In sample " << sample << " making 2D templates \n";
    makeTemplates(sample,era,conf::category,wpmin,conf::score_def,"200","1200"); 
  }
}

void makeTemplates(TString path2file, TString era, TString cat, TString wp,  TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cutmin="0==0", TString cutmax="0==0", TString suffix="none") {

  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2(kTRUE);
  conf::configuration(path2file);
  TString era_ = era; if (era == "2017" || era == "2018") { era_ = "2017"; }


  // Processes, paths and lumi for each era
  vector<TString> processes     = {"qcd","top","w-qq","w-qq","z-qq","z-qq","z-qq"};
  vector<TString> process_names = {conf::qcd.name,conf::tqq.name,conf::wcx.name,conf::wll.name,conf::zbb.name,conf::zcc.name,conf::zll.name};
  TString path;
  float intLumi;
  if (era == "2016") { 
    //    path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/ttbar1l/2016/"; 
    //    path = "/eos/uscms/store/user/pakontax/V2_Training_Official_nanoAODs_NEW_30Jan/2016/";
    path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200505/";
  //    path = "/eos/uscms/store/group/lpcjme/noreplica/NanoHRT/Trees/Apr08/muon/";
    intLumi= 36.8;
  }
  if (era == "2017") { 
    //path = "/eos/uscms/store/user/pakontax/V2_Training_Official_nanoAODs/2017_Muon_Channel/"; 
    path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2017/20201029_ak8/";
    intLumi= 44.98;
  }
  if (era == "2018") { 
    path = "/eos/uscms/store/user/pakontax/V2_Training_Official_nanoAODs/2018_Muon_Channel/"; 
    intLumi= 63.67;
  }
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();


  // Directory to store the templates   
  TString dirname1 = "templates2D";
  TString name0;
  if (score.Contains("tau21DDT"))      { name0 = "tau21ddt"; }
  if (score.Contains("DeepAK8_"))      { name0 = "dak8"; }
  if (score.Contains("DeepAK8MD"))     { name0 = "dak8md"; }
  if (score.Contains("DeepAK8DDT"))    { name0 = "dak8ddt"; }
  if (score.Contains("ParticleNetMD")) { name0 = "particlenetmd"; }
  TString nameoutfile = name0+"_zqq_"+cat+"_"+wp+"_"+era+"_"+cutmin+"to"+cutmax+"_templates";
  
  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }
 
  TFile *fout = new TFile("./"+dirname1+"/"+nameoutfile+".root","RECREATE");


  // corrections to mc
  TFile *f_v_ewknlo_wgt = TFile::Open("./corrections/kfactors.root" , "READONLY");
  TFile *f_w_qcdnlo_wgt = TFile::Open("./corrections/WJets_QCD_NLO.root" , "READONLY");
  TFile *f_z_qcdnlo_wgt = TFile::Open("./corrections/ZJets_QCD_NLO.root" , "READONLY");

  // Z bosons
  TH1D *h_z_qcdnlo = (TH1D*)f_v_ewknlo_wgt->Get("ZJets_012j_NLO/nominal");
  TH1D *h_z_ewknlo = (TH1D*)f_v_ewknlo_wgt->Get("EWKcorr/Z");
  TH1D *h_z_lo     = (TH1D*)f_v_ewknlo_wgt->Get("ZJets_LO/inv_pt");
  h_z_ewknlo_wgt   = (TH1D*)h_z_ewknlo->Clone("h_z_ewknlo_wgt"); h_z_ewknlo_wgt->Divide(h_z_qcdnlo);
  h_z_qcdnlo_wgt   = (TH1D*)f_z_qcdnlo_wgt->Get("Z_NLO_QCD_"+era_);

  // W bosons
  TH1D *h_w_qcdnlo = (TH1D*)f_v_ewknlo_wgt->Get("WJets_012j_NLO/nominal");
  TH1D *h_w_ewknlo = (TH1D*)f_v_ewknlo_wgt->Get("EWKcorr/W");
  TH1D *h_w_lo     = (TH1D*)f_v_ewknlo_wgt->Get("WJets_LO/inv_pt");
  h_w_ewknlo_wgt   = (TH1D*)h_w_ewknlo->Clone("h_w_ewknlo_wgt"); h_w_ewknlo_wgt->Divide(h_w_qcdnlo);
  h_w_qcdnlo_wgt   = (TH1D*)f_w_qcdnlo_wgt->Get("W_NLO_QCD_"+era_);


  // baseline selection
  TString wp_val = wp;
  TString c_p = "("+score+">="+wp+")";
  TString c_f = "(!"+c_p+")";
  TString c_p_ext = "("+score+">=0.54)";
  TString c_f_ext = "(!"+c_p_ext+")";

  TString cut_ = "(fj_1_pt>="+cutmin+" && fj_1_pt<"+cutmax+")";
  TString c_ht     = "ht>1000."; if (era != "2016") { c_ht     = "ht>1300."; }
  TString c_ht_qcd = "ht>1000."; if (era != "2016") { c_ht_qcd = "ht>1300."; }
  TString c_ht_ext = "ht>1000."; if (era != "2016") { c_ht_ext = "ht>700."; }
  TString c_base     = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=120000. && fj_1_sdmass>50. && fj_1_sdmass<250. && "+c_ht+" && met<150. && nlb_fj_pihalf==0 && nb_away_fj>=0) && ("+cut_+")";
  TString c_base_qcd = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=120000. && fj_1_sdmass>50. && fj_1_sdmass<250. && "+c_ht_qcd+" && met<150. && nlb_fj_pihalf==0 && nb_away_fj>=0) && ("+cut_+")";
  TString c_base_ext = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=120000. && fj_1_sdmass>50. && fj_1_sdmass<250. && ht>600. && nlb_fj_pihalf>=0) && ("+cut_+")";
  //TString c_base_ext = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=120000. && fj_1_sdmass>50. && fj_1_sdmass<250. && "+c_ht_ext+" && met<150. && nlb_fj_pihalf==0 && nb_away_fj>=0) && ("+cut_+")";
  TString c_incl     = c_base+" && "+cut_;
  TString c_incl_ext = c_base_ext+" && "+cut_;
  TString c_incl_qcd = c_base_qcd+" && "+cut_;
  std::vector<TString> cuts; cuts.clear();
  cuts.push_back(c_incl);
  cuts.push_back(c_incl_ext);
  cuts.push_back(c_p);     cuts.push_back(c_f);
  cuts.push_back(c_p_ext); cuts.push_back(c_f_ext);
  cuts.push_back(c_incl_qcd);

  // Fit variables and axis ranges
  TString brX = conf::brX; int binsX = conf::binsX; float minX = conf::minX;  float maxX = conf::maxX;
  TString brY = conf::brY; int binsY = conf::binsY; float minY = conf::minY;  float maxY = conf::maxY;
  TString name = path2file+"_"+name0+"_"+cat+"_"+wp+"_"+era;


  // Data histograms 
  TFile *f_data  = TFile::Open(path+"/jetht_tree.root" , "READONLY");
  TTree *t_data  = (TTree*)f_data->Get("Events");
  TH2D *h_data_p = create2Dhisto(name,t_data,lumi,cuts[0]+" && "+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_p",true);
  TH2D *h_data_f = create2Dhisto(name,t_data,lumi,cuts[0]+" && "+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_f",true);
  fout->cd();
  h_data_p->Write("data_obs_pass");
  h_data_f->Write("data_obs_fail");

  std:: cout << " Data = " << h_data_p->Integral() << "\n";

  // MC templates
  makeMCHistos(name,path,path2file,"nom","nom",lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
  //  makeMCHistos(name,path,path2file,"pu","puUp",lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
  //  makeMCHistos(name,path,path2file,"pu","puDown",lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);


  fout->Close();
  std::cout << "\n\n";
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

  // hack - it will not be needed in the new verison of ntuples
  if (score == "ak8_1_ParticleNetMD_XbbVsQCD") { score = "ak8_1_ParticleNetMD_Xbb/(ak8_1_ParticleNetMD_Xbb+ak8_1_ParticleNetMD_QCD)"; }

  // Processes, paths and lumi for each era
  //vector<TString> processes     = {"ttbar-powheg","singletop","ttv","w","diboson"};
  //vector<TString> process_names = {"tt","st","ttv","wqq","vv"};
  vector<TString> processes     = conf::processes;
  vector<TString> process_names = conf::process_names;

  TString path;
  float intLumi;
  if (era == "2016") { 
    path = conf::path_2016;
    intLumi= 36.8;
  }
  if (era == "2017") { 
    path = conf::path_2017;
    intLumi= 44.98;
  }
  if (era == "2018") { 
    path = conf::path_2018;
    intLumi= 63.67;
  }
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();


  // Directory to store the templates   
  TString dirname1 = "templates2D";
  TString name0;
  if      (score.Contains("tau21DDT"))      { name0 = "tau21ddt"; }
  else if (score.Contains("DeepAK8_"))      { name0 = "dak8"; }
  else if (score.Contains("DeepAK8MD"))     { name0 = "dak8md"; }
  else if (score.Contains("DeepAK8DDT"))    { name0 = "dak8ddt"; }
  else if (score.Contains("ParticleNetMD")) { name0 = "particlenetmd"; }
  else if (score.Contains("ParticleNet"))   { name0 = "particlenet"; }

  TString nameoutfile = conf::algo+"_tt1l_"+cat+"_"+wpmin+"to"+wpmax+"_"+era+"_"+cutmin+"to"+cutmax+"_templates";
  std::cout << " 2D templates name: " << nameoutfile << "\n";
  
  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }
 
  TFile *fout = new TFile("./"+dirname1+"/"+nameoutfile+".root","RECREATE");


  // Cuts and matching definition
  TString cut_ = "(passmetfilters && passMuTrig && fj_1_pt>="+cutmin+" && fj_1_pt<"+cutmax+")";
  TString c_base     = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && "+conf::brX+">"+conf::convertFloatToTString(conf::minX)+" && "+conf::brX+"<"+conf::convertFloatToTString(conf::maxX)+" && leptonicW_pt>150.) && ("+cut_+")";
  TString c_base_ext = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && "+conf::brX+">"+conf::convertFloatToTString(conf::minX)+" && "+conf::brX+"<"+conf::convertFloatToTString(conf::maxX)+" && leptonicW_pt>150.)";
  TString c_incl     = c_base+" && "+cut_;
  /*  
  TString c_p3 = "( (fj_1_dr_fj_top_wqmax<jetR) && (fj_1_dr_fj_top_b<jetR) )";
  TString c_p2 = "( (!"+c_p3+") && (fj_1_dr_fj_top_wqmax<jetR) && (fj_1_dr_fj_top_b>jetR) )";
  TString c_p1 = "(!("+c_p3+" || "+c_p2+"))";
  */
  TString c_p3 = "( (fj_1_dr_T<jetR) && (fj_1_dr_T_Wq_max<jetR) && (fj_1_dr_T_b<jetR) )";
  TString c_p2 = "( (!"+c_p3+") && (fj_1_T_Wq_max_pdgId<jetR) && (fj_1_dr_T_b>=jetR) )";
  TString c_p1 = "(!("+c_p3+" || "+c_p2+"))";

  std::vector<TString> cuts; cuts.clear();
  cuts.push_back(c_incl);
  cuts.push_back(c_incl+" && "+c_p3);
  cuts.push_back(c_incl+" && "+c_p2);
  cuts.push_back(c_incl+" && "+c_p1);


  // DDT maps ddt_dec2020
  //  TFile *f_dak8ddt_w_0p05 = TFile::Open("./corrections/WvsQCD_ddt_0p05rho.root" , "READONLY");
  TFile *f_dak8ddt_w_0p05 = TFile::Open("./corrections/ddt_dec2020/myDeepBoostedMap_0p05rho.root" , "READONLY");
  //TH2D *h_dak8ddt_w_0p05 = (TH2D*)f_dak8ddt_w_0p05->Get("DeepBoosted_WvsQCD_v_rho_v_pT_scaled_yx");
  TH2D *h_dak8ddt_w_0p05 = (TH2D*)f_dak8ddt_w_0p05->Get("DeepBoosted_WvsQCD_MD_v_rho_v_pT_scaled_yx_0p05");
  hrwgt_dak8ddt_w_0p05 = (TH2D*)h_dak8ddt_w_0p05->Clone("hrwgt_dak8ddt_w_0p05");

  TFile *f_dak8ddt_w_0p10 = TFile::Open("./corrections/ddt_dec2020/myDeepBoostedMap_0p10rho.root" , "READONLY");
  TH2D *h_dak8ddt_w_0p10 = (TH2D*)f_dak8ddt_w_0p10->Get("DeepBoosted_WvsQCD_MD_v_rho_v_pT_scaled_yx_0p10");
  hrwgt_dak8ddt_w_0p10 = (TH2D*)h_dak8ddt_w_0p10->Clone("hrwgt_dak8ddt_w_0p10");

  TFile *f_dak8ddt_w_0p20 = TFile::Open("./corrections/ddt_dec2020/myDeepBoostedMap_0p20rho.root" , "READONLY");
  TH2D *h_dak8ddt_w_0p20 = (TH2D*)f_dak8ddt_w_0p20->Get("DeepBoosted_WvsQCD_MD_v_rho_v_pT_scaled_yx_0p20");
  hrwgt_dak8ddt_w_0p20 = (TH2D*)h_dak8ddt_w_0p20->Clone("hrwgt_dak8ddt_w_0p20");

  // WP selection
  TString wp_val;

  TString c_p = "("+score+">"+wpmin+" && "+score+"<="+wpmax+")";
  TString c_f = "(!"+c_p+")";

  TString score_0p05, score_0p10, score_0p20;
  if (score.Contains("DeepAK8DDT")) {
    if (wpmin == "0p05") { score = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p05((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
    if (wpmin == "0p10") { score = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p10((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
    if (wpmin == "0p20") { score = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p20((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
    if (wpmin == "0p30") { score = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p30((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
    if (wpmin == "0p50") { score = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p50((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }

    c_p = "("+score+">=0)";
 
    if (wpmin == "0p10b") {
      score_0p05 = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p05((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))";
      score_0p10 = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p10((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))";
      c_p   = "("+score_0p10+">=0 && "+score_0p05+"<0)";     
    }

   if (wpmin == "0p05f") {
      score_0p05 = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p05((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))";
      score_0p20 = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p20((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))";
      c_p   = "("+score_0p20+">=0 && "+score_0p05+"<0)";     
   }

   if (wpmin == "0p10f") {
      score_0p10 = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p10((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))";
      score_0p20 = "(ak8_1_DeepAK8MD_WvsQCD-rewgtfuncDAK8DDT_w_0p20((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))";
      c_p   = "("+score_0p20+">=0 && "+score_0p10+"<0)";     
   }
    c_f = "(!"+c_p+")";

  }

  if (score.Contains("tau21DDT")) { 
    score = "((ak8_1_tau2/ak8_1_tau1) + 0.082*log(ak8_1_corr_sdmass/ak8_1_pt))";
    c_p = "("+score+"<0.43)";
    c_f = "(!"+c_p+")";
  }
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
      makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Up",lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
      makeMCHistosTop(name,path,processes,process_names,name_,namesys_+"Down",lumi,cuts,brX,binsX,minX,maxX,brY,binsY,minY,maxY,fout);
    }
  }

  /*
  // tt-herwig
  TFile *f_tt_h  = TFile::Open(path+"ttbar-herwig_tree.root" , "READONLY");
  TTree *t_tt_h  = (TTree*)f_tt_h->Get("Events");
  
  TH2D *h_tt_h_p1_p = create2Dhisto(name,t_tt_h,lumi,cuts[0]+"&&"+cuts[3]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p1_p",false);
  TH2D *h_tt_h_p2_p = create2Dhisto(name,t_tt_h,lumi,cuts[0]+"&&"+cuts[2]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p2_p",false);
  TH2D *h_tt_h_p3_p = create2Dhisto(name,t_tt_h,lumi,cuts[0]+"&&"+cuts[1]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p3_p",false);

  TH2D *h_tt_h_p1_f = create2Dhisto(name,t_tt_h,lumi,cuts[0]+"&&"+cuts[3]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p1_f",false);
  TH2D *h_tt_h_p2_f = create2Dhisto(name,t_tt_h,lumi,cuts[0]+"&&"+cuts[2]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p2_f",false);
  TH2D *h_tt_h_p3_f = create2Dhisto(name,t_tt_h,lumi,cuts[0]+"&&"+cuts[1]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p3_f",false);
  
  // avoid zero bins in mc
  for (unsigned int i0=0; i0<h_tt_h_p1_p->GetNbinsX(); ++i0) {
    if (h_tt_h_p1_p->GetBinContent(i0)<=0) { h_tt_h_p1_p->SetBinContent(i0,0.001);  h_tt_h_p1_p->SetBinError(i0,0.001); }
    if (h_tt_h_p2_p->GetBinContent(i0)<=0) { h_tt_h_p2_p->SetBinContent(i0,0.001);  h_tt_h_p2_p->SetBinError(i0,0.001); }
    if (h_tt_h_p3_p->GetBinContent(i0)<=0) { h_tt_h_p3_p->SetBinContent(i0,0.001);  h_tt_h_p3_p->SetBinError(i0,0.001); }

    if (h_tt_h_p1_f->GetBinContent(i0)<=0) { h_tt_h_p1_f->SetBinContent(i0,0.001);  h_tt_h_p1_f->SetBinError(i0,0.001); }
    if (h_tt_h_p2_f->GetBinContent(i0)<=0) { h_tt_h_p2_f->SetBinContent(i0,0.001);  h_tt_h_p2_f->SetBinError(i0,0.001); }
    if (h_tt_h_p3_f->GetBinContent(i0)<=0) { h_tt_h_p3_f->SetBinContent(i0,0.001);  h_tt_h_p3_f->SetBinError(i0,0.001); }  
  }
 
  // Write to root file
  h_tt_h_p1_p->Write("tt_herwig_p1_pass"); h_tt_h_p2_p->Write("tt_herwig_p2_pass"); h_tt_h_p3_p->Write("tt_herwig_p3_pass"); 
  h_tt_h_p1_f->Write("tt_herwig_p1_fail"); h_tt_h_p2_f->Write("tt_herwig_p2_fail"); h_tt_h_p3_f->Write("tt_herwig_p3_fail");
  */
  fout->Close();
  std::cout << "\n\n";
}


void makeMCHistos(TString name, TString path, TString path2file, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
		     TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_) { 

  conf::configuration(path2file);

  TString cc_match, bb_match, qq_match;
  bb_match = "(fj_1_nbhadrons>=1)";
  cc_match = "(fj_1_nbhadrons==0 && fj_1_nchadrons>=1)";
  qq_match = "(!("+bb_match+" || "+cc_match+"))";

  TString sys_type = "/";
  if ( sysType == "Up" ) { sys_type = "_up/"; } if ( sysType == "Down" ) { sys_type ="_down/"; }
		    
  TString sys_dir;  
  //if      ( (sys == "nom") || (sys == "pu") )               { sys_dir = "mc_nom/"; } 
  if      ( (sys == "nom") || (sys == "pu") )               { sys_dir = "/"; } 
  else if ( (sys.Contains("lhe")) || (sys.Contains("ps")) ) { sys_dir = "/LHEWeight/"; }
  else                                                      { sys_dir = "/"+sys+sys_type; }

  if ( (sys == "nom") || (sys == "pu") ) { name = name+"_"+sys; }
  else                                   { name = name+"_"+sys+sysType; }

  TString name_b;
  if (sys == "nom")      { name_b = ""; }
  else if (sys == "pu")  { name_b = sysType; }
  else                   { name_b = sys+sysType; }

  TFile *f_qcd   = TFile::Open(path+"/"+sys_dir+"qcd-mg_tree.root" , "READONLY");
  TFile *f_qcd_0 = TFile::Open(path+"/"+sys_dir+"qcd-mg-incl_tree.root" , "READONLY");
  TFile *f_qcd_1 = TFile::Open(path+"/"+sys_dir+"qcd-mg-benriched_tree.root" , "READONLY");
  TFile *f_qcd_2 = TFile::Open(path+"/"+sys_dir+"qcd-mg-bgenfilter_tree.root" , "READONLY");
  TFile *f_tt    = TFile::Open(path+"/"+sys_dir+"top_tree.root" , "READONLY");
  TFile *f_w     = TFile::Open(path+"/"+sys_dir+"w-qq_tree.root" , "READONLY");
  TFile *f_z     = TFile::Open(path+"/"+sys_dir+"z-qq_tree.root" , "READONLY");
  TFile *f_h     = TFile::Open(path+"/"+sys_dir+"ggH-bb_tree.root" , "READONLY");

  TTree *t_qcd   = (TTree*)f_qcd->Get("Events");
  TTree *t_qcd_0 = (TTree*)f_qcd_0->Get("Events");
  TTree *t_qcd_1 = (TTree*)f_qcd_1->Get("Events");
  TTree *t_qcd_2 = (TTree*)f_qcd_2->Get("Events");
  TTree *t_tt    = (TTree*)f_tt->Get("Events");
  TTree *t_w     = (TTree*)f_w->Get("Events");
  TTree *t_z     = (TTree*)f_z->Get("Events");

  TFile *f_qcd_her;
  TTree *t_qcd_her;
  if (sys == "nom") {
  f_qcd_her = TFile::Open(path+"/"+sys_dir+"qcd-herwig_tree.root" , "READONLY");
  t_qcd_her = (TTree*)f_qcd_her->Get("Events");
  }

  std::vector<TH2D*>   h2ds;       h2ds.clear();
  std::vector<TString> h2ds_names; h2ds_names.clear();

  // create pass templates
  TH2D *h_qcd_p      = create2Dhisto(name,t_qcd,wgts,cuts[6]+"&&"+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_p",false);
  TH2D *h_qcd_0_p    = create2Dhisto(name,t_qcd_0,wgts,cuts[6]+"&&"+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_0_p",false);
  TH2D *h_qcd_1_p    = create2Dhisto(name,t_qcd_1,wgts,cuts[6]+"&&"+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_1_p",false);
  TH2D *h_qcd_2_p    = create2Dhisto(name,t_qcd_2,wgts,cuts[6]+"&&"+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_2_p",false);
  TH2D *h_tt_p       = create2Dhisto(name,t_tt,wgts,cuts[0]+"&&"+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p",false);
  TH2D *h_wcx_p_i    = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[2]+" && fj_1_isW<0.8 && (fj_1_nchadrons>=1)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosoncc_p_i",false);
  TH2D *h_wll_p_i    = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[2]+" && fj_1_isW<0.8 && (fj_1_nchadrons==0)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosonll_p_i",false);
  TH2D *h_w_m_p_i    = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[2]+" && fj_1_isW<0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_p_i",false);
  TH2D *h_w_um_p_i   = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[2]+" && fj_1_isW>=0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_p_i",false);
  TH2D *h_zbb_m_p_i  = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ<0.8 && "+bb_match+" && "+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonbb_m_p_i",false);
  TH2D *h_zcc_m_p_i  = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ<0.8 && "+cc_match+" && "+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosoncc_m_p_i",false);
  TH2D *h_zqq_m_p_i  = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ<0.8 && "+qq_match+" && "+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonqq_m_p_i",false);
  TH2D *h_z_um_p_i   = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ>=0.8 && "+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_p_i",false);

  // take the boson templates using a relaxed selection
  TH2D *h_wcx_p    = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[4]+" && fj_1_isW<0.8 && (fj_1_nchadrons>=1)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosoncc_p",false);
  TH2D *h_wll_p    = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[4]+" && fj_1_isW<0.8 && (fj_1_nchadrons==0)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosonll_p",false);
  TH2D *h_w_m_p    = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[4]+" && fj_1_isW<0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_p",false);
  TH2D *h_w_um_p   = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[4]+" && fj_1_isW>=0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_p",false);
  TH2D *h_zbb_m_p  = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ<0.8 && "+bb_match+" && "+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonbb_m_p",false);
  TH2D *h_zcc_m_p  = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ<0.8 && "+cc_match+" && "+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosoncc_m_p",false);
  TH2D *h_zqq_m_p  = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ<0.8 && "+qq_match+" && "+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonqq_m_p",false);
  TH2D *h_z_um_p   = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ>=0.8 && "+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_p",false);

  // correct normalization
  h_wcx_p    = correct_2dext_norm(h_wcx_p_i   , h_wcx_p   ); 
  h_wll_p    = correct_2dext_norm(h_wll_p_i   , h_wll_p   );
  h_w_m_p    = correct_2dext_norm(h_w_m_p_i   , h_w_m_p   );
  h_w_um_p   = correct_2dext_norm(h_w_um_p_i  , h_w_um_p  ); 
  h_zbb_m_p  = correct_2dext_norm(h_zbb_m_p_i , h_zbb_m_p );
  h_zcc_m_p  = correct_2dext_norm(h_zcc_m_p_i , h_zcc_m_p );
  h_zqq_m_p  = correct_2dext_norm(h_zqq_m_p_i , h_zqq_m_p );
  h_z_um_p   = correct_2dext_norm(h_z_um_p_i  , h_z_um_p  );

  // store them in a vector
  std::cout << "name_b " << name_b << "\n"; 
  h2ds.push_back(h_qcd_p);   h2ds_names.push_back("h_"+conf::qcd.name+name_b+"_pass");
  h2ds.push_back(h_qcd_0_p); h2ds_names.push_back("h_qcd_0"+name_b+"_pass");
  h2ds.push_back(h_qcd_1_p); h2ds_names.push_back("h_qcd_1"+name_b+"_pass");
  h2ds.push_back(h_qcd_2_p); h2ds_names.push_back("h_qcd_2"+name_b+"_pass");
  h2ds.push_back(h_tt_p);    h2ds_names.push_back("h_"+conf::tqq.name+name_b+"_pass");
  h2ds.push_back(h_wcx_p);   h2ds_names.push_back("h_"+conf::wcx.name+name_b+"_pass");
  h2ds.push_back(h_wll_p);   h2ds_names.push_back("h_"+conf::wll.name+name_b+"_pass");
  h2ds.push_back(h_w_m_p);   h2ds_names.push_back("h_"+conf::wqq.name+name_b+"_pass");
  h2ds.push_back(h_w_um_p);  h2ds_names.push_back("h_"+conf::wll.name+name_b+"_pass_unmatched");
  h2ds.push_back(h_zbb_m_p); h2ds_names.push_back("h_"+conf::zbb.name+name_b+"_pass");
  h2ds.push_back(h_zcc_m_p); h2ds_names.push_back("h_"+conf::zcc.name+name_b+"_pass");
  h2ds.push_back(h_zqq_m_p); h2ds_names.push_back("h_"+conf::zll.name+name_b+"_pass");
  h2ds.push_back(h_z_um_p);  h2ds_names.push_back("h_"+conf::zll.name+name_b+"_pass_unmatched");


  std::cout << " qcd incl = " << h_qcd_p->Integral() 
	    << " qcd 0 = " << h_qcd_0_p->Integral()
	    << " qcd 1 = " << h_qcd_1_p->Integral()
	    << " qcd 2 = " << h_qcd_2_p->Integral() << "\n";
    
  // create fail template
  TH2D *h_qcd_f      = create2Dhisto(name,t_qcd,wgts,cuts[6]+"&&"+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_f",false);
  TH2D *h_qcd_0_f    = create2Dhisto(name,t_qcd_0,wgts,cuts[6]+"&&"+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_0_f",false);
  TH2D *h_qcd_1_f    = create2Dhisto(name,t_qcd_1,wgts,cuts[6]+"&&"+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_1_f",false);
  TH2D *h_qcd_2_f    = create2Dhisto(name,t_qcd_2,wgts,cuts[6]+"&&"+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_2_f",false);
  TH2D *h_tt_f       = create2Dhisto(name,t_tt,wgts,cuts[0]+"&&"+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_f",false);
  TH2D *h_wcx_f_i    = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[3]+" && fj_1_isW<0.8 && (fj_1_nchadrons>=1)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosoncc_f_i",false);
  TH2D *h_wll_f_i    = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[3]+" && fj_1_isW<0.8 && (fj_1_nchadrons==0)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosonll_f_i",false);
  TH2D *h_w_m_f_i    = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[3]+" && fj_1_isW<0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_f_i",false);
  TH2D *h_w_um_f_i   = create2Dhisto(name,t_w,wgts,cuts[0]+" &&"+cuts[3]+" && fj_1_isW>=0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_f_i",false);
  TH2D *h_zbb_m_f_i  = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ<0.8 && "+bb_match+" && "+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonbb_m_f_i",false);
  TH2D *h_zcc_m_f_i  = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ<0.8 && "+cc_match+" && "+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosoncc_m_f_i",false);
  TH2D *h_zqq_m_f_i  = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ<0.8 && "+qq_match+" && "+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonqq_m_f_i",false);
  TH2D *h_z_um_f_i   = create2Dhisto(name,t_z,wgts,cuts[0]+" && fj_1_isZ>=0.8 && "+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_f_i",false);

  // take the boson templates using a relaxed selection
  TH2D *h_wcx_f    = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[5]+" && fj_1_isW<0.8 && (fj_1_nchadrons>=1)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosoncc_f",false);
  TH2D *h_wll_f    = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[5]+" && fj_1_isW<0.8 && (fj_1_nchadrons==0)",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wbosonll_f",false);
  TH2D *h_w_m_f    = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[5]+" && fj_1_isW<0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_f",false);
  TH2D *h_w_um_f   = create2Dhisto(name,t_w,wgts,cuts[1]+" &&"+cuts[5]+" && fj_1_isW>=0.8",brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_f",false);
  TH2D *h_zbb_m_f  = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ<0.8 && "+bb_match+" && "+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonbb_m_f",false);
  TH2D *h_zcc_m_f  = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ<0.8 && "+cc_match+" && "+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosoncc_m_f",false);
  TH2D *h_zqq_m_f  = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ<0.8 && "+qq_match+" && "+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zbosonqq_m_f",false);
  TH2D *h_z_um_f   = create2Dhisto(name,t_z,wgts,cuts[1]+" && fj_1_isZ>=0.8 && "+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_f",false);

  // correct normalization
  h_wcx_f    = correct_2dext_norm(h_wcx_f_i   , h_wcx_f   ); 
  h_wll_f    = correct_2dext_norm(h_wll_f_i   , h_wll_f   );
  h_w_m_f    = correct_2dext_norm(h_w_m_f_i   , h_w_m_f   );
  h_w_um_f   = correct_2dext_norm(h_w_um_f_i  , h_w_um_f  ); 
  h_zbb_m_f  = correct_2dext_norm(h_zbb_m_f_i , h_zbb_m_f );
  h_zcc_m_f  = correct_2dext_norm(h_zcc_m_f_i , h_zcc_m_f );
  h_zqq_m_f  = correct_2dext_norm(h_zqq_m_f_i , h_zqq_m_f );
  h_z_um_f   = correct_2dext_norm(h_z_um_f_i  , h_z_um_f  );

  // store them in a vector
  h2ds.push_back(h_qcd_f);   h2ds_names.push_back("h_"+conf::qcd.name+name_b+"_fail");
  h2ds.push_back(h_qcd_0_f); h2ds_names.push_back("h_qcd_0"+name_b+"_fail");
  h2ds.push_back(h_qcd_1_f); h2ds_names.push_back("h_qcd_1"+name_b+"_fail");
  h2ds.push_back(h_qcd_2_f); h2ds_names.push_back("h_qcd_2"+name_b+"_fail");
  h2ds.push_back(h_tt_f);    h2ds_names.push_back("h_"+conf::tqq.name+name_b+"_fail");
  h2ds.push_back(h_wcx_f);   h2ds_names.push_back("h_"+conf::wcx.name+name_b+"_fail");
  h2ds.push_back(h_wll_f);   h2ds_names.push_back("h_"+conf::wll.name+name_b+"_fail");
  h2ds.push_back(h_w_m_f);   h2ds_names.push_back("h_"+conf::wqq.name+name_b+"_fail");
  h2ds.push_back(h_w_um_f);  h2ds_names.push_back("h_"+conf::wll.name+name_b+"_fail_unmatched");
  h2ds.push_back(h_zbb_m_f); h2ds_names.push_back("h_"+conf::zbb.name+name_b+"_fail");
  h2ds.push_back(h_zcc_m_f); h2ds_names.push_back("h_"+conf::zcc.name+name_b+"_fail");
  h2ds.push_back(h_zqq_m_f); h2ds_names.push_back("h_"+conf::zll.name+name_b+"_fail");
  h2ds.push_back(h_z_um_f);  h2ds_names.push_back("h_"+conf::zll.name+name_b+"_fail_unmatched");


  TH2D *h_qcd_her_p;
  TH2D *h_qcd_her_f;
  if (sys == "nom") {
    h_qcd_her_p = create2Dhisto(name,t_qcd_her,wgts,cuts[0]+"&&"+cuts[2],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_qcd_herwig_p",false); 
    h2ds.push_back(h_qcd_her_p); h2ds_names.push_back("h_qcd_herwig_pass");
    h_qcd_her_f = create2Dhisto(name,t_qcd_her,wgts,cuts[0]+"&&"+cuts[3],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_qcd_herwig_f",false); 
    h2ds.push_back(h_qcd_her_f); h2ds_names.push_back("h_qcd_herwig_fail");
  }


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

}




void makeMCHistosTop(TString name, TString path, std::vector<TString> processes, std::vector<TString> process_names, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
		     TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_) { 

  TString sys_type = "/";
  if ( sysType == "Up" ) { sys_type = "_up/"; } if ( sysType == "Down" ) { sys_type ="_down/"; }
		    
  TString sys_dir; 
  if      ( (sys == "nom") || (sys == "pu") || (sys == "jms") || (sys == "jmr") ) { sys_dir = "/mc/"; }  
  else if ( (sys.Contains("lhe")) || (sys.Contains("ps")) ) { sys_dir = "/LHEWeight/"; }
  else                                                      { sys_dir = "/"+sys+sys_type; }

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


void makeMCHistosLHEPDFTop(TString name, TString path, std::vector<TString> processes, std::vector<TString> process_names, TString sys, TString sysType, TString wgts, std::vector<TString> cuts,
			   TString brX, int binsX, float minX, float maxX, TString brY, int binsY, float minY, float maxY, TFile *f_) { 
  
  TString sys_dir =  "/LHEWeight/";
  name = name+"_"+sys;
  TString name_b = "lhepdf";
  
  std::vector<TH2D*>   h2ds;       h2ds.clear();
  std::vector<TString> h2ds_names; h2ds_names.clear();
 
  for (unsigned int i0=0; i0<processes.size(); ++i0) {
    
    TFile *f = TFile::Open(path+"/"+sys_dir+processes[i0]+"_tree.root","READONLY");
    TTree *t = (TTree*)f->Get("Events");
    
    std::vector<TH2D*> tmph2d_p;    tmph2d_p.clear();
    std::vector<TH2D*> tmph2d_f;    tmph2d_f.clear();
    std::vector<TH2D*> tmph2d_p3_p; tmph2d_p3_p.clear();
    std::vector<TH2D*> tmph2d_p2_p; tmph2d_p2_p.clear();
    std::vector<TH2D*> tmph2d_p1_p; tmph2d_p1_p.clear();
    std::vector<TH2D*> tmph2d_p3_f; tmph2d_p3_f.clear();
    std::vector<TH2D*> tmph2d_p2_f; tmph2d_p2_f.clear();
    std::vector<TH2D*> tmph2d_p1_f; tmph2d_p1_f.clear();
    for (unsigned i=0; i<30; ++i) {
 
      TString count = std::to_string(i);
      TString pdfwgt_ = "LHEPdfWeight["+count+"]*LHEPdfWeightNorm["+count+"]/(LHEPdfWeight[0]*LHEPdfWeightNorm[0])";
      TString wgts_ = wgts+"*"+pdfwgt_;
	
      if ( (process_names[i0] == "tt") || (process_names[i0] == "st") || (process_names[i0] == "ttv") ) {
 
	TH2D *h_p3_p = create2Dhisto(name,t,wgts_,cuts[1]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_p3_p",false); tmph2d_p3_p.push_back(h_p3_p);   
	TH2D *h_p2_p = create2Dhisto(name,t,wgts_,cuts[2]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_p2_p",false); tmph2d_p2_p.push_back(h_p2_p); 
	TH2D *h_p1_p = create2Dhisto(name,t,wgts_,cuts[3]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_p1_p",false); tmph2d_p1_p.push_back(h_p1_p);
      
	TH2D *h_p3_f = create2Dhisto(name,t,wgts_,cuts[1]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_p3_f",false); tmph2d_p3_f.push_back(h_p3_f);
        TH2D *h_p2_f = create2Dhisto(name,t,wgts_,cuts[2]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_p2_f",false); tmph2d_p2_f.push_back(h_p2_f);
        TH2D *h_p1_f = create2Dhisto(name,t,wgts_,cuts[3]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_p1_f",false); tmph2d_p1_f.push_back(h_p1_f);
      }
      else {
	TH2D *h_p = create2Dhisto(name,t,wgts_,cuts[0]+"&&"+cuts[4],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_p",false); tmph2d_p.push_back(h_p); 
	TH2D *h_f = create2Dhisto(name,t,wgts_,cuts[0]+"&&"+cuts[5],brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_"+process_names[i0]+"_"+count+"_f",false); tmph2d_f.push_back(h_f); 
      }
    } // end over the pdf variations

    // get mean and envelope

    if ( (process_names[i0] == "tt") || (process_names[i0] == "st") || (process_names[i0] == "ttv") ) {
      std::vector<TH2D*> h_p3_p_ = getLHPDFRmsAndSigma(tmph2d_p3_p); std::vector<TH2D*> h_p3_f_ = getLHPDFRmsAndSigma(tmph2d_p3_f);
      std::vector<TH2D*> h_p2_p_ = getLHPDFRmsAndSigma(tmph2d_p2_p); std::vector<TH2D*> h_p2_f_ = getLHPDFRmsAndSigma(tmph2d_p2_f);
      std::vector<TH2D*> h_p1_p_ = getLHPDFRmsAndSigma(tmph2d_p1_p); std::vector<TH2D*> h_p1_f_ = getLHPDFRmsAndSigma(tmph2d_p1_f);
      
      h2ds.push_back(h_p3_p_[0]); h2ds_names.push_back(process_names[i0]+"_p3_"+name_b+"Up_pass");
      h2ds.push_back(h_p3_p_[1]); h2ds_names.push_back(process_names[i0]+"_p3_"+name_b+"Down_pass");
      h2ds.push_back(h_p2_p_[0]); h2ds_names.push_back(process_names[i0]+"_p2_"+name_b+"Up_pass");
      h2ds.push_back(h_p2_p_[1]); h2ds_names.push_back(process_names[i0]+"_p2_"+name_b+"Down_pass");
      h2ds.push_back(h_p1_p_[0]); h2ds_names.push_back(process_names[i0]+"_p1_"+name_b+"Up_pass");
      h2ds.push_back(h_p1_p_[1]); h2ds_names.push_back(process_names[i0]+"_p1_"+name_b+"Down_pass");

      h2ds.push_back(h_p3_f_[0]); h2ds_names.push_back(process_names[i0]+"_p3_"+name_b+"Up_fail");
      h2ds.push_back(h_p3_f_[1]); h2ds_names.push_back(process_names[i0]+"_p3_"+name_b+"Down_fail");
      h2ds.push_back(h_p2_f_[0]); h2ds_names.push_back(process_names[i0]+"_p2_"+name_b+"Up_fail");
      h2ds.push_back(h_p2_f_[1]); h2ds_names.push_back(process_names[i0]+"_p2_"+name_b+"Down_fail");
      h2ds.push_back(h_p1_f_[0]); h2ds_names.push_back(process_names[i0]+"_p1_"+name_b+"Up_fail");
      h2ds.push_back(h_p1_f_[1]); h2ds_names.push_back(process_names[i0]+"_p1_"+name_b+"Down_fail");
    }

    else {
      std::vector<TH2D*> h_p_ = getLHPDFRmsAndSigma(tmph2d_p); std::vector<TH2D*> h_f_ = getLHPDFRmsAndSigma(tmph2d_f);
      
      h2ds.push_back(h_p_[0]); h2ds_names.push_back(process_names[i0]+"_"+name_b+"Up_pass");
      h2ds.push_back(h_p_[1]); h2ds_names.push_back(process_names[i0]+"_"+name_b+"Down_pass");

      h2ds.push_back(h_f_[0]); h2ds_names.push_back(process_names[i0]+"_"+name_b+"Up_fail");
      h2ds.push_back(h_f_[1]); h2ds_names.push_back(process_names[i0]+"_"+name_b+"Down_fail");
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
  
} // end of makeMCHistosLHEPDFTop


std::vector<TH2D*> getLHPDFRmsAndSigma(std::vector<TH2D*> input) {
  
  TH2D *htmp_up   = (TH2D*)input[0]->Clone("htmp_up"); 
  TH2D *htmp_down = (TH2D*)input[0]->Clone("htmp_down");
  
  for (unsigned int ibinX=0; ibinX<input[0]->GetNbinsX(); ++ibinX) {
    for (unsigned int ibinY=0; ibinY<input[0]->GetNbinsY(); ++ibinY) {
    float mean_ = 0.; float rms_ = 0.; float sigma_ = 0.;
    
    for (unsigned int ihisto=0; ihisto<input.size(); ++ihisto) {
      mean_ += input[ihisto]->GetBinContent(ibinX,ibinY);                                         
      rms_  += (input[ihisto]->GetBinContent(ibinX,ibinY))*(input[ihisto]->GetBinContent(ibinX,ibinY)); 
    }
    mean_ = mean_/input.size(); rms_ = sqrt(rms_/input.size()); sigma_ = sqrt((rms_*rms_) - (mean_*mean_)); 
    htmp_up->SetBinContent(ibinX,ibinY,(mean_ + sigma_)); 
    htmp_down->SetBinContent(ibinX,ibinY,(mean_ - sigma_));
    }
  }

  std::vector<TH2D*> output; output.clear();
  output.push_back(htmp_up);
  output.push_back(htmp_down);
  
  return output;
}


TH2D *correct_2dext_norm(TH2D *h_nom, TH2D *h_ext) {
  //std::cout << h_ext->GetNbinsX() << "  " << h_ext->GetNbinsY() << "\n";
  for (int iy=1; iy<=h_ext->GetNbinsY(); ++iy) {
    float int_nom = h_nom->Integral(1,h_nom->GetNbinsX(),iy,iy); //std::cout << " nom = " << int_nom << "\n";
    float int_ext = h_ext->Integral(1,h_ext->GetNbinsX(),iy,iy); //std::cout << " ext = " << int_ext << "\n";
    float scale = int_nom/int_ext;                               //std::cout << " scale = " << scale << "\n";

    for (int ix=1; ix<=h_ext->GetNbinsX(); ++ix) {
      h_ext->SetBinContent(ix,iy,scale*h_ext->GetBinContent(ix,iy));
      h_ext->SetBinError(ix,iy,scale*h_ext->GetBinError(ix,iy));
    }

  }
  
  return h_ext;
}


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
    //if (name.Contains("qcd"))           { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+"*"+qcdWgt+")*("+cuts+")"; }
    if (name.Contains("qcd_0"))         { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+"*"+qcdWgt+")*("+cuts+" && qcdSampleType==0)"; }
    if (name.Contains("qcd_1") || name.Contains("qcd_2") ) { cut = "("+intLumi+"*"+puWgt+"*"+genWgt+"*"+qcdWgt+")*("+cuts+" && qcdSampleType>0)"; }
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
  //if (name.Contains("jms"))      { std::cout << " In jms \n"; tree->Project(name,branchY+":(*"+branchX+")",cut); }
  else if (name.Contains("jmr")) { std::cout << " In jmr \n"; tree->Project(name,branchY+":(massSmear("+branchX+",luminosityBlock,event,"+massSmearVal_+"))",cut); }
  else                           { std::cout << " In else \n"; tree->Project(name,branchY+":"+branchX,cut); }

  // ad overflow bin
  for (unsigned int i0y=0; i0y<hTemp->GetNbinsY(); ++i0y) {
    double error =0.; double integral = hTemp->IntegralAndError(hTemp->GetNbinsX(),hTemp->GetNbinsX()+1,i0y+1,i0y+1,error);
    hTemp->SetBinContent(hTemp->GetNbinsX(),i0y+1,integral);
    hTemp->SetBinError(hTemp->GetNbinsX(),i0y+1,error);
  }

  for (unsigned int i0x=0; i0x<hTemp->GetNbinsX(); ++i0x) {
    double error =0.; double integral = hTemp->IntegralAndError(i0x+1,i0x+1,hTemp->GetNbinsY(),hTemp->GetNbinsY()+1,error);
    hTemp->SetBinContent(i0x+1,hTemp->GetNbinsY(),integral);
    hTemp->SetBinError(i0x+1,hTemp->GetNbinsY(),error);
  } 

  return hTemp;
}


void makeHisto(TString path_, std::vector<TString> names_, std::vector<TString> wgt_, std::vector<float> scalep_, std::vector<float> scalef_, std::vector<TString> cuts_, 
	       std::vector<TString> c_p, std::vector<TString> c_f, TString sys_, TString var_, int bins_, float xmin_, float xmax_, TFile *fout_p, TFile *fout_f) {
  

  for (unsigned int i=0; i<names_.size(); ++i) {

    TString fname;
    if (names_[i] == "wboson") { fname = "w-qq"; }
    if (names_[i] == "zboson") { fname = "z-qq"; }
    if (names_[i] == "higgs")  { fname = "ggH-bb"; }
    if (names_[i] == "tt")     { fname = "top"; }
    if (names_[i] == "qcd")    { fname = "qcd-mg"; }
    
    TFile *f_mc_; 
    if ((sys_=="pu")) { f_mc_ = TFile::Open(path_+"/"+fname+"_tree.root","READONLY"); }
    else if (sys_.Contains("lhe") || sys_.Contains("ps"))     { f_mc_ = TFile::Open(path_+"/mc_sys/LHEWeight/"+fname+"_tree.root","READONLY"); }
    else                               { f_mc_ = TFile::Open(path_+"/mc_sys/"+sys_+"_"+var_+"/"+fname+"_tree.root","READONLY"); }
    TTree *t_mc_ = (TTree*)f_mc_->Get("Events");

    TH1D *h_p_up = create1Dhisto(names_[i]+"_p_"+sys_+"Up",t_mc_,wgt_[i],cuts_[i]+"&&"+c_p[i],var_,bins_,xmin_,xmax_,false,1,1,"h_"+names_[i]+"_p_"+sys_+"Up",false,false); 
    h_p_up->Scale(scalep_[i]);
    TH1D *h_p_down = create1Dhisto(names_[i]+"_p_"+sys_+"Down",t_mc_,wgt_[i],cuts_[i]+"&&"+c_p[i],var_,bins_,xmin_,xmax_,false,1,1,"h_"+names_[i]+"_p_"+sys_+"Down",false,false); 
    h_p_down->Scale(scalep_[i]);
    TH1D *h_f_up = create1Dhisto(names_[i]+"_f_"+sys_+"Up",t_mc_,wgt_[i],cuts_[i]+"&&"+c_f[i],var_,bins_,xmin_,xmax_,false,1,1,"h_"+names_[i]+"_f_"+sys_+"Up",false,false);
    h_f_up->Scale(scalef_[i]);
    TH1D *h_f_down = create1Dhisto(names_[i]+"_f_"+sys_+"Down",t_mc_,wgt_[i],cuts_[i]+"&&"+c_f[i],var_,bins_,xmin_,xmax_,false,1,1,"h_"+names_[i]+"_f_"+sys_+"Down",false,false); 
    h_f_down->Scale(scalef_[i]);


    for (unsigned int ii=0; ii<h_p_up->GetNbinsX(); ++ii) {
      if (h_p_up->GetBinContent(ii)<=0) { h_p_up->SetBinContent(ii,0.001); h_p_up->SetBinError(ii,0.001); }
      if (h_p_down->GetBinContent(ii)<=0) { h_p_down->SetBinContent(ii,0.001); h_p_down->SetBinError(ii,0.001); }
      if (h_f_up->GetBinContent(ii)<=0) { h_f_up->SetBinContent(ii,0.001); h_f_up->SetBinError(ii,0.001); }
      if (h_f_down->GetBinContent(ii)<=0) { h_f_down->SetBinContent(ii,0.001); h_f_down->SetBinError(ii,0.001); }
    }

    fout_p->cd();
    h_p_up->Write(names_[i]+"_"+sys_+"Up");
    h_p_down->Write(names_[i]+"_"+sys_+"Down");
    
    fout_f->cd();
    h_f_up->Write(names_[i]+"_"+sys_+"Up");
    h_f_down->Write(names_[i]+"_"+sys_+"Down");
    
  }

} // end of makeHisto




TH1D *create1Dhisto(TString sample,TTree *tree,TString intLumi,TString cuts,TString branch,int bins,float xmin,float xmax,
		    bool useLog,int color, int style,TString name,bool norm, bool data) {
  TH1::SetDefaultSumw2(kTRUE);


  TString puWgt;
  if (name.Contains("puUp"))        { puWgt = "puWeightUp"; }
  else if (name.Contains("puDown")) { puWgt = "puWeightDown"; }
  else                              { puWgt = "puWeight"; }

  TString genWgt = "xsecWeight*genWeight";

  // additinal weights
  TString wBosonWgt = "(w_qcdnlo_wgt_func(fj_1_pt)*w_ewknlo_wgt_func(fj_1_pt))"; //NLOQCD=1.35 for 2016
  TString zBosonWgt = "(z_qcdnlo_wgt_func(fj_1_pt)*z_ewknlo_wgt_func(fj_1_pt))"; //NLOQCD=1.45 for 2016
  //TString qcdWgt    = "qcd_incl_mass_wgt(fj_1_sdmass)";
  //TString qcdWgt    = "qcd_incl_mass_wgt(fj_1_pt)";
  TString qcdWgt    = "1.";
  TString ttWgt     = "1."; if (sample.Contains("tt1L") || sample.Contains("ttbar1L")) { ttWgt = "topptWeight"; };
  // ttbar: add top-pt reweighting

  TString cut;
  if (data) { cut ="("+cuts+")"; } 
  //  if (data) { cut ="( (pass_mutrig || pass_eltrig) && "+cuts+")"; } 
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


  TH1D *hTemp = new TH1D(name,name,bins,xmin,xmax); //hTemp->SetName(name);
  tree->Project(name,branch,cut);

  hTemp->SetFillColor(0);
  hTemp->SetLineWidth(3);
  hTemp->SetMarkerSize(0);
  hTemp->SetLineColor(color);
  hTemp->SetFillColor(color);
  hTemp->SetLineStyle(style);

  // ad overflow bin             
  double error =0.; double integral = hTemp->IntegralAndError(bins,bins+1,error);
  hTemp->SetBinContent(bins,integral);
  hTemp->SetBinError(bins,error);

  if (norm) { hTemp->Scale(1./(hTemp->Integral())); }

  return hTemp;
} //~ end of create1Dhisto



void testarea() {

  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2(kTRUE);

  TString era = "2016";

  TFile *f_v_ewknlo_wgt = TFile::Open("./corrections/kfactors.root" , "READONLY");
  TFile *f_w_qcdnlo_wgt = TFile::Open("./corrections/WJets_QCD_NLO.root" , "READONLY");
  TFile *f_z_qcdnlo_wgt = TFile::Open("./corrections/ZJets_QCD_NLO.root" , "READONLY");


  // Z bosons
  TH1D *h_z_qcdnlo = (TH1D*)f_v_ewknlo_wgt->Get("ZJets_012j_NLO/nominal");
  TH1D *h_z_ewknlo = (TH1D*)f_v_ewknlo_wgt->Get("EWKcorr/Z");
  TH1D *h_z_lo     = (TH1D*)f_v_ewknlo_wgt->Get("ZJets_LO/inv_pt");
  h_z_ewknlo_wgt   = (TH1D*)h_z_ewknlo->Clone("h_z_ewknlo_wgt"); h_z_ewknlo_wgt->Divide(h_z_qcdnlo);
  h_z_qcdnlo_wgt   = (TH1D*)f_z_qcdnlo_wgt->Get("Z_NLO_QCD_"+era);

  // W bosons
  TH1D *h_w_qcdnlo = (TH1D*)f_v_ewknlo_wgt->Get("WJets_012j_NLO/nominal");
  TH1D *h_w_ewknlo = (TH1D*)f_v_ewknlo_wgt->Get("EWKcorr/W");
  TH1D *h_w_lo     = (TH1D*)f_v_ewknlo_wgt->Get("WJets_LO/inv_pt");
  h_w_ewknlo_wgt   = (TH1D*)h_w_ewknlo->Clone("h_w_ewknlo_wgt"); h_w_ewknlo_wgt->Divide(h_w_qcdnlo);
  h_w_qcdnlo_wgt   = (TH1D*)f_w_qcdnlo_wgt->Get("W_NLO_QCD_"+era);

  TString path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200505/";
  TFile *f = TFile::Open(path+"z-qq_tree.root" , "READONLY");
  TTree *t = (TTree*)f->Get("Events");
 
  float intLumi = 36.8;
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();

  TString c_match = "(fj_1_nbhadrons>=2)";
  TString c_p = "(fj_1_ParticleNetMD_XbbVsQCD>=0.93)";
  TString c_f = "(!"+c_p+")";
  TString c_p_ext = "(fj_1_ParticleNetMD_XbbVsQCD>=0.54)";

  TString cut_ = "(fj_1_pt>=200 && fj_1_pt<1200000)";
  TString c_base     = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=120000. && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>1000. && met<150000000. && nmb_fj_pihalf==0 && nb_away_fj>=0) && ("+cut_+")";
  TString c_base_ext = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=120000. && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>200. && nlb_fj_pihalf>=0) && ("+cut_+")";
  TString c_incl     = c_base+" && "+cut_;

  TString brX = "fj_1_sdmass"; int binsX=30; float minX = 50;  float maxX=200;
  //TString brX = "fj_1_sdmass"; int binsX=15; float minX = 50;  float maxX=200;
  TString brY = "fj_1_pt";     int binsY=40;  float minY = 200; float maxY=1200;
  
  TH1D *h_histo       = create1Dhisto("zboson",t,lumi,c_incl+"&&"+c_p+"&& fj_1_isZ<0.8 && ht>1000 && (fj_1_pt>=450 && fj_1_pt<500)",brX,binsX,minX,maxX,false,8,1,"h_zboson",false,false); h_histo->SetFillColor(0);
  h_histo->SetMarkerSize(0.);
  TH2D *h_histo2d_ext = create2Dhisto("zboson2d_ext",t,lumi,c_base_ext+"&& fj_1_isZ<0.8 &&"+c_p_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h2d_zboson_m_p_i",false);
  TH2D *h_histo2d   = create2Dhisto("zboson2d",t,lumi,c_incl+"&& fj_1_isZ<0.8 && ht>1000 && "+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h2d_zboson_m_p_b",false);
  /*
  TH1D *h_histo     = create1Dhisto("zboson",t,lumi,c_incl+"&&"+c_p+"&& ht>1000 && "+cut_,brX,binsX,minX,maxX,false,8,1,"h_zboson",false,false); h_histo->SetFillColor(0);
  TH2D *h_histo2d   = create2Dhisto("zboson2d",t,lumi,c_incl+"&& 0==0 &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h2d_zboson_m_p_i",false);
  TH2D *h_histo2d_b = create2Dhisto("zboson2d_b",t,lumi,c_incl+"&& ht>1000 && "+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h2d_zboson_m_p_b",false);  
  */

  //TH2D *h_histo2d_ext_b = correct_2dext_norm(h_histo2d, h_histo2d_ext);
  h_histo2d_ext = correct_2dext_norm(h_histo2d, h_histo2d_ext);

  TH1D *h1dfrom2d   = h_histo2d->ProjectionX("h1dfrom2d",h_histo2d->GetYaxis()->FindBin(450.),h_histo2d->GetYaxis()->FindBin(499.),"e");
  //TH1D *h1dfrom2d   = h_histo2d->ProfileX("h1dfrom2d",h_histo2d->GetYaxis()->FindBin(450.),h_histo2d->GetYaxis()->FindBin(500.),"e");
  h1dfrom2d->SetLineColor(2); h1dfrom2d->SetMarkerSize(0.);
  TH1D *h1dfrom2d_ext = h_histo2d_ext->ProjectionX("h1dfrom2d_ext",h_histo2d_ext->GetYaxis()->FindBin(450.),h_histo2d_ext->GetYaxis()->FindBin(499.),"e");
  h1dfrom2d_ext->SetLineColor(1); h1dfrom2d_ext->SetMarkerSize(0.);
  //h1dfrom2d_ext->Scale(h_histo2d->Integral()/h1dfrom2d_ext->Integral());
  //TH1D *h1dfrom2d_ext_b = h_histo2d_ext_b->ProjectionX("h1dfrom2d_ext_b",h_histo2d_ext_b->GetYaxis()->FindBin(450.),h_histo2d_ext_b->GetYaxis()->FindBin(499.),"e");
  //h1dfrom2d_ext_b->SetLineColor(4); h1dfrom2d_ext_b->SetMarkerSize(0.);


  std::cout << h_histo->Integral() << " " << h1dfrom2d->Integral() << " " << h1dfrom2d_ext->Integral() << " " << h_histo2d->Integral() <<  " " 
	    << h_histo2d->Integral(1,30,h_histo2d->GetYaxis()->FindBin(450.),h_histo2d->GetYaxis()->FindBin(499.)) << " " 
    
	    << "\n";

  TCanvas *c_ = new TCanvas("c_","c_",500,500);
  h_histo->Draw("HIST E0");
  h1dfrom2d_ext->Draw("HIST E0 sames");
  h1dfrom2d->Draw("HIST E0 sames");
  //h1dfrom2d_ext->Draw("HIST E0 sames");
}

