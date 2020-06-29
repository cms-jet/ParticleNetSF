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
void makeTemplatesTop(TString path2file, TString era, TString cat, TString wp, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cutmin="0==0", TString cutmax="0==0", TString suffix="none");
TH2D *correct_2dext_norm(TH2D *h_nom, TH2D *h_ext);
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

  TString cut_ = "(fj_1_pt>=200 && fj_1_pt<1200)";
  TString c_base     = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=1200. && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>1000. && met<150000000. && nmb_fj_pihalf==0 && nb_away_fj>=0) && ("+cut_+")";
  TString c_base_ext = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=1200. && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>200. && nlb_fj_pihalf>=0) && ("+cut_+")";
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


void mainfunction(TString sample) {
  if ( (sample == "tt1L") || (sample=="ttbar1L") ) { std::cout << sample << "\n";
    //    makeTemplatesTop(sample,"2016","bb","t","ak8_1_ParticleNetMD_XbbVsQCD","200","1200");
    makeTemplatesTop(sample,"2016","bb","t","ak8_1_ParticleNetMD_XbbVsQCD","200","1200"); 
  }
  if (sample == "zqq") { std::cout << sample<< "\n";
    //    makeTemplates(sample,"2016","bb","t","fj_1_ParticleNetMD_XbbVsQCD","200","1200");
    makeTemplates(sample,"2016","bb","t","fj_1_ParticleNetMD_XbbVsQCD","200","1200"); 
  }
}


void makeTemplates(TString path2file, TString era, TString cat, TString wp, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cutmin="0==0", TString cutmax="0==0", TString suffix="none") {

  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2(kTRUE);

  TString path;
  if (era == "2016") { path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200505/"; }
  if (era == "2017") { path = "/eos/uscms/store/group/lpcjme/LGOutputs/20190522_softb_2017/"; }
  if (era == "2018") { path = "/eos/uscms/store/group/lpcjme/LGOutputs/20190522_softb_2018/"; }

  TFile *f_data  = TFile::Open(path+"jetht_tree.root" , "READONLY");
  TFile *f_qcd   = TFile::Open(path+"qcd-mg_tree.root" , "READONLY");
  TFile *f_tt    = TFile::Open(path+"top_tree.root" , "READONLY");
  TFile *f_w     = TFile::Open(path+"w-qq_tree.root" , "READONLY");
  TFile *f_z     = TFile::Open(path+"z-qq_tree.root" , "READONLY");
  TFile *f_h     = TFile::Open(path+"ggH-bb_tree.root" , "READONLY");
  TFile *f_qcd_h = TFile::Open(path+"qcd-herwig_tree.root" , "READONLY");
  TFile *f_sm    = TFile::Open(path+"sm_tree.root" , "READONLY");

  TTree *t_data  = (TTree*)f_data->Get("Events");
  TTree *t_qcd   = (TTree*)f_qcd->Get("Events");
  TTree *t_tt    = (TTree*)f_tt->Get("Events");
  TTree *t_w     = (TTree*)f_w->Get("Events");
  TTree *t_z     = (TTree*)f_z->Get("Events");
  TTree *t_h     = (TTree*)f_h->Get("Events");
  TTree *t_qcd_h = (TTree*)f_qcd_h->Get("Events");
  TTree *t_sm    = (TTree*)f_sm->Get("Events");
  
  float intLumi;
  if (era == "2016") { intLumi= 36.8; }
  if (era == "2017") { intLumi= 44.98; }
  if (era == "2018") { intLumi= 63.67; }
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();

  // corrections to mc
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


  // baseline selection
  TString c_match;
  TString wp_val;
  if (cat == "bb") { 
    c_match = "(fj_1_nbhadrons>=2)";
    if (wp == "t") { wp_val = "0.93"; }
    if (wp == "m") { wp_val = "0.84"; }
    if (wp == "l") { wp_val = "0.54"; }
  }
  if (cat == "cc") { 
    c_match = "(fj_1_nbhadrons==0 && fj_1_nchadrons>=2)";
    if (wp == "t") { wp_val = "0.94"; }
    if (wp == "m") { wp_val = "0.90"; }
    if (wp == "l") { wp_val = "0.76"; }
  }  
  TString c_p = "("+score+">="+wp_val+")";
  TString c_f = "(!"+c_p+")";
  TString c_p_ext = "("+score+">=0.54)";
  TString c_f_ext = "(!"+c_p_ext+")";

  TString cut_ = "(fj_1_pt>="+cutmin+" && fj_1_pt<"+cutmax+")";
  TString c_base     = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=1200. && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>1000. && met<150. && nlb_fj_pihalf==0 && nb_away_fj>=0) && ("+cut_+")";
  TString c_base_ext = "(abs(fj_1_eta)<2.4 && fj_1_pt>=200. && fj_1_pt<=1200. && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>200. && nlb_fj_pihalf>=0) && ("+cut_+")";
  TString c_incl     = c_base+" && "+cut_;

  TString brX = "fj_1_sdmass"; int binsX=30; float minX = 50;  float maxX=200;
  //  TString brX = "fj_1_sdmass"; int binsX=15; float minX = 50;  float maxX=200;
  TString brY = "fj_1_pt";     int binsY=40;  float minY = 200; float maxY=1200;
  //  TString brX = "fj_1_sdmass"; int binsX=23; float minX = 40;  float maxX=201;
  //TString brY = "fj_1_pt";     int binsY=6;  float minY = 450; float maxY=1200;
  TString name = path2file+"_particlenet_"+cat+"_"+wp+"_"+era;


  // create pass templates
  TH2D *h_data_p   = create2Dhisto(name,t_data,lumi,c_incl+"&& ht>1000. && "+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_p",true);
  TH2D *h_qcd_p    = create2Dhisto(name,t_qcd,lumi,c_incl+"&&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_f",false);
  TH2D *h_tt_p     = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p",false);
  TH2D *h_tt_m_p   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_m_p",false);  
  TH2D *h_tt_um_p  = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_um_p",false);
  TH2D *h_w_p_i    = create2Dhisto(name,t_w,lumi,c_incl+" &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_p_i",false);
  TH2D *h_w_m_p_i  = create2Dhisto(name,t_w,lumi,c_incl+"&& fj_1_isW<0.8 &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_p_i",false);
  TH2D *h_w_um_p_i = create2Dhisto(name,t_w,lumi,c_incl+"&& fj_1_isW>=0.8 &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_p_i",false);
  TH2D *h_z_p_i    = create2Dhisto(name,t_z,lumi,c_incl+" &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_p_i",false);
  TH2D *h_z_m_p_i  = create2Dhisto(name,t_z,lumi,c_incl+"&& fj_1_isZ<0.8 &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_m_p_i",false);
  TH2D *h_z_um_p_i = create2Dhisto(name,t_z,lumi,c_incl+"&& fj_1_isZ>=0.8 &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_p_i",false);

  // take the boson templates using a relaxed selection
  TH2D *h_w_p    = create2Dhisto(name,t_w,lumi,c_base_ext+" && "+c_p_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_p",false);
  TH2D *h_w_m_p  = create2Dhisto(name,t_w,lumi,c_base_ext+" && fj_1_isW<0.8 &&"+c_p_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_p",false);
  TH2D *h_w_um_p = create2Dhisto(name,t_w,lumi,c_base_ext+" && fj_1_isW>=0.8 &&"+c_p_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_p",false);
  TH2D *h_z_p    = create2Dhisto(name,t_z,lumi,c_base_ext+" && "+c_p_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_p",false);
  TH2D *h_z_m_p  = create2Dhisto(name,t_z,lumi,c_base_ext+" && fj_1_isZ<0.8 &&"+c_p_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_m_p",false);
  TH2D *h_z_um_p = create2Dhisto(name,t_z,lumi,c_base_ext+" && fj_1_isZ>=0.8 &&"+c_p_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_p",false);

  // correct normalization
  h_z_p    = correct_2dext_norm(h_z_p_i    , h_z_p   ); 
  h_z_m_p  = correct_2dext_norm(h_z_m_p_i  , h_z_m_p );
  h_z_um_p = correct_2dext_norm(h_z_um_p_i , h_z_um_p);
  h_w_p    = correct_2dext_norm(h_w_p_i    , h_w_p   ); 
  h_w_m_p  = correct_2dext_norm(h_w_m_p_i  , h_w_m_p );
  h_w_um_p = correct_2dext_norm(h_w_um_p_i , h_w_um_p);
  /*
  float scale_w_p    = h_w_p_i->Integral()/h_w_p->Integral(); h_w_p->Scale(scale_w_p); 
  float scale_w_m_p  = h_w_m_p_i->Integral()/h_w_m_p->Integral(); h_w_m_p->Scale(scale_w_m_p);
  float scale_w_um_p = h_w_um_p_i->Integral()/h_w_um_p->Integral(); h_w_um_p->Scale(scale_w_um_p);

  float scale_z_p    = h_z_p_i->Integral()/h_z_p->Integral(); h_z_p->Scale(scale_z_p);
  //float scale_z_m_p  = h_z_m_p_i->Integral()/h_z_m_p->Integral(); //h_z_m_p->Scale(scale_z_m_p);
  float scale_z_m_p  = h_z_m_p_i->Integral()/h_z_m_p->Integral(); h_z_m_p->Scale(scale_z_m_p);
  float scale_z_um_p = h_z_um_p_i->Integral()/h_z_um_p->Integral(); h_z_um_p->Scale(scale_z_um_p);
  */

  // create fail template
  TH2D *h_data_f   = create2Dhisto(name,t_data,lumi,c_incl+"&& ht>1000. && "+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_f",true);
  TH2D *h_qcd_f    = create2Dhisto(name,t_qcd,lumi,c_incl+"&&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_f",false);
  TH2D *h_tt_f     = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_f",false);
  TH2D *h_tt_m_f   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_m_f",false);  
  TH2D *h_tt_um_f  = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_um_f",false);
  TH2D *h_w_f_i    = create2Dhisto(name,t_w,lumi,c_incl+" &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_f_i",false);
  TH2D *h_w_m_f_i  = create2Dhisto(name,t_w,lumi,c_incl+"&& fj_1_isW<0.8 &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_f_i",false);
  TH2D *h_w_um_f_i = create2Dhisto(name,t_w,lumi,c_incl+"&& fj_1_isW>=0.8 &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_f_i",false);
  TH2D *h_z_f_i    = create2Dhisto(name,t_z,lumi,c_incl+" &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_f_i",false);
  TH2D *h_z_m_f_i  = create2Dhisto(name,t_z,lumi,c_incl+"&& fj_1_isZ<0.8 &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_m_f_i",false);
  TH2D *h_z_um_f_i = create2Dhisto(name,t_z,lumi,c_incl+"&& fj_1_isZ>=0.8 &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_f_i",false);

  // take the boson templates using a relaxed selection
  TH2D *h_w_f    = create2Dhisto(name,t_w,lumi,c_base_ext+" &&"+c_f_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_f",false);
  TH2D *h_w_m_f  = create2Dhisto(name,t_w,lumi,c_base_ext+"&& fj_1_isW<0.8 &&"+c_f_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_m_f",false);
  TH2D *h_w_um_f = create2Dhisto(name,t_w,lumi,c_base_ext+"&& fj_1_isW>=0.8 &&"+c_f_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wboson_um_f",false);
  TH2D *h_z_f    = create2Dhisto(name,t_z,lumi,c_base_ext+" &&"+c_f_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_f",false);
  TH2D *h_z_m_f  = create2Dhisto(name,t_z,lumi,c_base_ext+"&& fj_1_isZ<0.8 &&"+c_f_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_m_f",false);
  TH2D *h_z_um_f = create2Dhisto(name,t_z,lumi,c_base_ext+"&& fj_1_isZ>=0.8 &&"+c_f_ext,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_zboson_um_f",false);

  // correct normalization
  h_z_f    = correct_2dext_norm(h_z_f_i    , h_z_f   ); 
  h_z_m_f  = correct_2dext_norm(h_z_m_f_i  , h_z_m_f );
  h_z_um_f = correct_2dext_norm(h_z_um_f_i , h_z_um_f);
  h_w_f    = correct_2dext_norm(h_w_f_i    , h_w_f   ); 
  h_w_m_f  = correct_2dext_norm(h_w_m_f_i  , h_w_m_f );
  h_w_um_f = correct_2dext_norm(h_w_um_f_i , h_w_um_f);

  
 /*
  float scale_w_f    = h_w_f_i->Integral()/h_w_f->Integral(); h_w_f->Scale(scale_w_f);
  float scale_w_m_f  = h_w_m_f_i->Integral()/h_w_m_f->Integral(); h_w_m_f->Scale(scale_w_m_f); 
  float scale_w_um_f = h_w_um_f_i->Integral()/h_w_um_f->Integral(); h_w_um_f->Scale(scale_w_um_f);
  
  float scale_z_f    = h_z_f_i->Integral()/h_z_f->Integral(); h_z_f->Scale(scale_z_f);
  float scale_z_m_f  = h_z_m_f_i->Integral()/h_z_m_f->Integral(); h_z_m_f->Scale(scale_z_m_f); 
  float scale_z_um_f = h_z_um_f_i->Integral()/h_z_um_f->Integral(); h_z_um_f->Scale(scale_z_um_f);
 */
  

  // Cross check yields
  std::cout << "pass:\n";
  std::cout << "N(data) = " << h_data_p->Integral() << " , N(sm) = " << h_qcd_p->Integral()+h_tt_p->Integral()+h_w_p->Integral()+h_z_p->Integral() << "\n";
  std::cout << "qcd = " << h_qcd_p->Integral() << " , ";
  std::cout << "tt  = " << h_tt_p->Integral() << " , ";
  std::cout << "zqq = " << h_z_m_p_i->Integral() << " " << h_z_m_p->Integral() << " " 
	    << h_z_m_p_i->Integral(1,30,h_z_m_p_i->GetYaxis()->FindBin(450.),h_z_m_p_i->GetYaxis()->FindBin(500.)) << " " 
	    << h_z_m_p->Integral(1,30,h_z_m_p->GetYaxis()->FindBin(450.),h_z_m_p->GetYaxis()->FindBin(500.))
	    << " , ";
  std::cout << "wqq = " << h_w_p->Integral() << " " << h_w_m_p->Integral() << "\n\n";
 
  std::cout << "fail:\n";
  std::cout << "N(data) = " << h_data_f->Integral() << " , N(sm) = " << h_qcd_f->Integral()+h_tt_f->Integral()+h_w_f->Integral()+h_z_f->Integral() << "\n";
  std::cout << "qcd = "<< h_qcd_f->Integral() << " , ";
  std::cout << "tt  = " << h_tt_f->Integral() << " , ";
  std::cout << "zqq = " << h_z_f->Integral() << " " << h_z_m_f->Integral() << " , ";
  std::cout << "wqq = " << h_w_f->Integral() << " " << h_w_m_f->Integral() << "\n\n";
  std::cout << "\n";
  

  // -------
  // systematics

  // herwig
  TH2D *h_qcd_h_f = create2Dhisto(name,t_qcd_h,lumi,c_incl+"&&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_h_f",false);
  TH2D *h_qcd_h_p = create2Dhisto(name,t_qcd_h,lumi,c_incl+"&&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_qcd_h_p",false);
 

  // avoid zero bins in mc
  for (unsigned int ii=0; ii<h_qcd_p->GetNbinsX(); ++ii) {
    if (h_qcd_p->GetBinContent(ii)<=0)   { h_qcd_p->SetBinContent(ii,0.001);   h_qcd_p->SetBinError(ii,0.001); }
    if (h_tt_p->GetBinContent(ii)<=0)    { h_tt_p->SetBinContent(ii,0.001);    h_tt_p->SetBinError(ii,0.001); }
    if (h_w_p->GetBinContent(ii)<=0)     { h_w_p->SetBinContent(ii,0.001);     h_w_p->SetBinError(ii,0.001); }
    if (h_w_m_p->GetBinContent(ii)<=0)   { h_w_m_p->SetBinContent(ii,0.001);   h_w_m_p->SetBinError(ii,0.001); }
    if (h_w_um_p->GetBinContent(ii)<=0)  { h_w_um_p->SetBinContent(ii,0.001);  h_w_um_p->SetBinError(ii,0.001); }
    if (h_z_p->GetBinContent(ii)<=0)     { h_z_p->SetBinContent(ii,0.001);     h_z_p->SetBinError(ii,0.001); }
    if (h_z_m_p->GetBinContent(ii)<=0)   { h_z_m_p->SetBinContent(ii,0.001);   h_z_m_p->SetBinError(ii,0.001); }
    if (h_z_um_p->GetBinContent(ii)<=0)  { h_z_um_p->SetBinContent(ii,0.001);  h_z_um_p->SetBinError(ii,0.001); }
    if (h_qcd_h_p->GetBinContent(ii)<=0) { h_qcd_h_p->SetBinContent(ii,0.001); h_qcd_h_p->SetBinError(ii,0.001); }

    if (h_qcd_f->GetBinContent(ii)<=0)   { h_qcd_f->SetBinContent(ii,0.001);   h_qcd_f->SetBinError(ii,0.001); }
    if (h_tt_f->GetBinContent(ii)<=0)    { h_tt_f->SetBinContent(ii,0.001);    h_tt_f->SetBinError(ii,0.001); }
    if (h_w_f->GetBinContent(ii)<=0)     { h_w_f->SetBinContent(ii,0.001);     h_w_f->SetBinError(ii,0.001); }
    if (h_w_m_f->GetBinContent(ii)<=0)   { h_w_m_f->SetBinContent(ii,0.001);   h_w_m_f->SetBinError(ii,0.001); }
    if (h_w_um_f->GetBinContent(ii)<=0)  { h_w_um_f->SetBinContent(ii,0.001);  h_w_um_f->SetBinError(ii,0.001); }
    if (h_z_f->GetBinContent(ii)<=0)     { h_z_f->SetBinContent(ii,0.001);     h_z_f->SetBinError(ii,0.001); }
    if (h_z_m_f->GetBinContent(ii)<=0)   { h_z_m_f->SetBinContent(ii,0.001);   h_z_m_f->SetBinError(ii,0.001); }
    if (h_z_um_f->GetBinContent(ii)<=0)  { h_z_um_f->SetBinContent(ii,0.001);  h_z_um_f->SetBinError(ii,0.001); }
    if (h_qcd_h_f->GetBinContent(ii)<=0) { h_qcd_h_f->SetBinContent(ii,0.001); h_qcd_h_f->SetBinError(ii,0.001); }
   
 }
 
  // make dir
  TString dirname1 = "particlenet_"+cat+"_"+wp+"_"+era;
  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }

  TString nameoutfile = "particlenet_"+cat+"_"+wp+"_"+era+"_"+cutmin+"to"+cutmax+"_templates"; 
  TFile *fout = new TFile("./"+dirname1+"/"+nameoutfile+".root","RECREATE");

  h_data_p->Write("data_obs_pass"); 
  h_qcd_p->Write("qcd_pass"); 
  h_tt_p->Write("tqq_pass"); h_tt_m_p->Write("tqq_pass_matched"); h_tt_um_p->Write("tqq_pass_unmatched");
  h_w_p->Write("wqq_pass"); h_w_m_p->Write("wqq_pass_matched"); h_w_um_p->Write("wqq_pass_unmatched");
  h_z_p->Write("zqq_pass"); h_z_m_p->Write("zqq_pass_matched"); h_z_um_p->Write("zqq_pass_unmatched");

  h_data_f->Write("data_obs_fail"); 
  h_qcd_f->Write("qcd_fail"); 
  h_tt_f->Write("tqq_fail"); h_tt_m_f->Write("tqq_fail_matched"); h_tt_um_f->Write("tqq_fail_unmatched");
  h_w_f->Write("wqq_fail"); h_w_m_f->Write("wqq_fail_matched"); h_w_um_f->Write("wqq_fail_unmatched");
  h_z_f->Write("zqq_fail"); h_z_m_f->Write("zqq_fail_matched"); h_z_um_f->Write("zqq_fail_unmatched");

  // systematics
  h_qcd_p->Write("qcd_herwig_pass"); h_qcd_f->Write("qcd_herwig_fail");
  
  fout->Close();
  std::cout << "\n\n";
}


void makeTemplatesTop(TString path2file, TString era, TString cat, TString wp, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cutmin="0==0", TString cutmax="0==0", TString suffix="none") {

  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2(kTRUE);

  TString path;
  if (era == "2016") { path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/ttbar1l/2016/"; }
  if (era == "2017") { path = "/eos/uscms/store/group/lpcjme/LGOutputs/20190522_softb_2017/"; }
  if (era == "2018") { path = "/eos/uscms/store/group/lpcjme/LGOutputs/20190522_softb_2018/"; }

  TFile *f_data  = TFile::Open(path+"singlemu_tree.root"       , "READONLY");
  TFile *f_tt    = TFile::Open(path+"ttbar-powheg_tree.root" , "READONLY");
  TFile *f_st    = TFile::Open(path+"singletop_tree.root" , "READONLY");
  TFile *f_w     = TFile::Open(path+"w_tree.root"         , "READONLY");
  TFile *f_vv    = TFile::Open(path+"diboson_tree.root"         , "READONLY");
  TFile *f_tt_h  = TFile::Open(path+"ttbar-herwig_tree.root" , "READONLY");
  TFile *f_sm    = TFile::Open(path+"sm_tree.root" , "READONLY");

  TTree *t_data  = (TTree*)f_data->Get("Events");
  TTree *t_tt    = (TTree*)f_tt->Get("Events");
  TTree *t_st    = (TTree*)f_st->Get("Events");
  TTree *t_w     = (TTree*)f_w->Get("Events");
  TTree *t_vv    = (TTree*)f_vv->Get("Events");
  TTree *t_tt_h  = (TTree*)f_tt_h->Get("Events");
  TTree *t_sm    = (TTree*)f_sm->Get("Events");
  
  float intLumi;
  if (era == "2016") { intLumi= 36.8; }
  if (era == "2017") { intLumi= 44.98; }
  if (era == "2018") { intLumi= 63.67; }
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();


  // selections
  TString wp_val;
  if (cat == "bb") {
    if (wp == "t") { wp_val = "0.93"; }
    if (wp == "m") { wp_val = "0.84"; }
    if (wp == "l") { wp_val = "0.54"; }
  }
  if (cat == "cc") {
    if (wp == "t") { wp_val = "0.94"; }
    if (wp == "m") { wp_val = "0.90"; }
    if (wp == "l") { wp_val = "0.76"; }
  }

  TString c_p = "("+score+">="+wp_val+")";
  TString c_f = "(!"+c_p+")";

  // matching definition
  TString c_p3 = "( (ak8_1_dr_fj_top_wqmax<0.8) && (ak8_1_dr_fj_top_b<0.8) )";
  TString c_p2 = "( (!"+c_p3+") && (ak8_1_dr_fj_top_wqmax<0.8) && (ak8_1_dr_fj_top_b>0.8) )";
  TString c_p1 = "(!("+c_p3+" || "+c_p2+"))";

  TString cut_ = "(ak8_1_pt>="+cutmin+" && ak8_1_pt<"+cutmax+")";
  TString c_base     = "(abs(ak8_1_eta)<2.4 && ak8_1_pt>=200. && ak8_1_mass>50. && ak8_1_mass<200.) && ("+cut_+")";
  TString c_base_ext = "(abs(ak8_1_eta)<2.4 && ak8_1_pt>=200. && ak8_1_mass>50. && ak8_1_mass<200.)";
  TString c_incl     = c_base+" && "+cut_;

  TString brX = "ak8_1_mass"; int binsX=30; float minX = 50;  float maxX=200;
  //TString brX = "ak8_1_mass"; int binsX=15; float minX = 50;  float maxX=200;
  TString brY = "ak8_1_pt";     int binsY=40;  float minY = 200; float maxY=1200;
  TString name = path2file+"_particlenet_"+cat+"_"+wp+"_"+era;


  // create pass templates
  TH2D *h_data_p    = create2Dhisto(name,t_data,lumi,c_incl+" && "+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_p",true);
  TH2D *h_tt_p1_p   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_p+"&&"+c_p1,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p1_p",false);
  TH2D *h_tt_p2_p   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_p+"&&"+c_p2,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p2_p",false);
  TH2D *h_tt_p3_p   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_p+"&&"+c_p3,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p3_p",false);
  TH2D *h_st_p1_p   = create2Dhisto(name,t_st,lumi,c_incl+"&&"+c_p+"&&"+c_p1,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_st_p1_p",false);
  TH2D *h_st_p2_p   = create2Dhisto(name,t_st,lumi,c_incl+"&&"+c_p+"&&"+c_p2,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_st_p2_p",false);
  TH2D *h_st_p3_p   = create2Dhisto(name,t_st,lumi,c_incl+"&&"+c_p+"&&"+c_p3,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_st_p3_p",false);
  TH2D *h_wqq_p     = create2Dhisto(name,t_w,lumi,c_incl+" &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wqq_p",false);
  TH2D *h_vv_p      = create2Dhisto(name,t_vv,lumi,c_incl+" &&"+c_p,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_vv_p",false);
  
  // create fail template
  TH2D *h_data_f    = create2Dhisto(name,t_data,lumi,c_incl+" && "+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_data_f",true);
  TH2D *h_tt_p1_f   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_f+"&&"+c_p1,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p1_f",false);
  TH2D *h_tt_p2_f   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_f+"&&"+c_p2,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p2_f",false);
  TH2D *h_tt_p3_f   = create2Dhisto(name,t_tt,lumi,c_incl+"&&"+c_f+"&&"+c_p3,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_p3_f",false);
  TH2D *h_st_p1_f   = create2Dhisto(name,t_st,lumi,c_incl+"&&"+c_f+"&&"+c_p1,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_st_p1_f",false);
  TH2D *h_st_p2_f   = create2Dhisto(name,t_st,lumi,c_incl+"&&"+c_f+"&&"+c_p2,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_st_p2_f",false);
  TH2D *h_st_p3_f   = create2Dhisto(name,t_st,lumi,c_incl+"&&"+c_f+"&&"+c_p3,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_st_p3_f",false);
  TH2D *h_wqq_f     = create2Dhisto(name,t_w,lumi,c_incl+" &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_wqq_f",false);
  TH2D *h_vv_f      = create2Dhisto(name,t_vv,lumi,c_incl+" &&"+c_f,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_vv_f",false);
 
  

  // Cross check yields
  std::cout << "pass: N(data) = " << h_data_p->Integral() << " , N(sm) = " << h_tt_p1_p->Integral()+h_tt_p2_p->Integral()+h_tt_p3_p->Integral()
    + h_st_p1_p->Integral() + h_st_p2_p->Integral() + h_st_p3_p->Integral() + h_wqq_p->Integral() + h_vv_p->Integral() << "\n";
  std::cout << "fail: N(data) = " << h_data_f->Integral() << " , N(sm) = " << h_tt_p1_f->Integral()+h_tt_p2_f->Integral()+h_tt_p3_f->Integral()
    + h_st_p1_f->Integral() + h_st_p2_f->Integral() + h_st_p3_f->Integral() + h_wqq_f->Integral() + h_vv_f->Integral() << "\n";
  std::cout << "\n";


  // -------
  // systematics

  // herwig
  TH2D *h_tt_h_p1_p = create2Dhisto(name,t_tt_h,lumi,c_incl+"&&"+c_p+"&&"+c_p1,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p1_p",false);
  TH2D *h_tt_h_p2_p = create2Dhisto(name,t_tt_h,lumi,c_incl+"&&"+c_p+"&&"+c_p2,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p2_p",false);
  TH2D *h_tt_h_p3_p = create2Dhisto(name,t_tt_h,lumi,c_incl+"&&"+c_p+"&&"+c_p3,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p3_p",false);

  TH2D *h_tt_h_p1_f = create2Dhisto(name,t_tt_h,lumi,c_incl+"&&"+c_f+"&&"+c_p1,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p1_f",false);
  TH2D *h_tt_h_p2_f = create2Dhisto(name,t_tt_h,lumi,c_incl+"&&"+c_f+"&&"+c_p2,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p2_f",false);
  TH2D *h_tt_h_p3_f = create2Dhisto(name,t_tt_h,lumi,c_incl+"&&"+c_f+"&&"+c_p3,brX,binsX,minX,maxX,brY,binsY,minY,maxY,false,"h_"+name+"_tt_h_p3_f",false);

  // avoid zero bins in mc
  for (unsigned int ii=0; ii<h_tt_p1_p->GetNbinsX(); ++ii) {
    if (h_tt_p1_p->GetBinContent(ii)<=0)   { h_tt_p1_p->SetBinContent(ii,0.001);    h_tt_p1_p->SetBinError(ii,0.001); }
    if (h_tt_p2_p->GetBinContent(ii)<=0)   { h_tt_p2_p->SetBinContent(ii,0.001);    h_tt_p2_p->SetBinError(ii,0.001); }
    if (h_tt_p3_p->GetBinContent(ii)<=0)   { h_tt_p3_p->SetBinContent(ii,0.001);    h_tt_p3_p->SetBinError(ii,0.001); }
    if (h_st_p1_p->GetBinContent(ii)<=0)   { h_st_p1_p->SetBinContent(ii,0.001);    h_st_p1_p->SetBinError(ii,0.001); }
    if (h_st_p2_p->GetBinContent(ii)<=0)   { h_st_p2_p->SetBinContent(ii,0.001);    h_st_p2_p->SetBinError(ii,0.001); }
    if (h_st_p3_p->GetBinContent(ii)<=0)   { h_st_p3_p->SetBinContent(ii,0.001);    h_st_p3_p->SetBinError(ii,0.001); }
    if (h_wqq_p->GetBinContent(ii)<=0)     { h_wqq_p->SetBinContent(ii,0.001);      h_wqq_p->SetBinError(ii,0.001); }
    if (h_vv_p->GetBinContent(ii)<=0)      { h_vv_p->SetBinContent(ii,0.001);       h_vv_p->SetBinError(ii,0.001); }
    if (h_tt_h_p1_p->GetBinContent(ii)<=0) { h_tt_h_p1_p->SetBinContent(ii,0.001);  h_tt_h_p1_p->SetBinError(ii,0.001); }
    if (h_tt_h_p2_p->GetBinContent(ii)<=0) { h_tt_h_p2_p->SetBinContent(ii,0.001);  h_tt_h_p2_p->SetBinError(ii,0.001); }
    if (h_tt_h_p3_p->GetBinContent(ii)<=0) { h_tt_h_p3_p->SetBinContent(ii,0.001);  h_tt_h_p3_p->SetBinError(ii,0.001); }

    if (h_tt_p1_f->GetBinContent(ii)<=0)   { h_tt_p1_f->SetBinContent(ii,0.001);    h_tt_p1_f->SetBinError(ii,0.001); }
    if (h_tt_p2_f->GetBinContent(ii)<=0)   { h_tt_p2_f->SetBinContent(ii,0.001);    h_tt_p2_f->SetBinError(ii,0.001); }
    if (h_tt_p3_f->GetBinContent(ii)<=0)   { h_tt_p3_f->SetBinContent(ii,0.001);    h_tt_p3_f->SetBinError(ii,0.001); }
    if (h_st_p1_f->GetBinContent(ii)<=0)   { h_st_p1_f->SetBinContent(ii,0.001);    h_st_p1_f->SetBinError(ii,0.001); }
    if (h_st_p2_f->GetBinContent(ii)<=0)   { h_st_p2_f->SetBinContent(ii,0.001);    h_st_p2_f->SetBinError(ii,0.001); }
    if (h_st_p3_f->GetBinContent(ii)<=0)   { h_st_p3_f->SetBinContent(ii,0.001);    h_st_p3_f->SetBinError(ii,0.001); }
    if (h_wqq_f->GetBinContent(ii)<=0)     { h_wqq_f->SetBinContent(ii,0.001);      h_wqq_f->SetBinError(ii,0.001); }
    if (h_vv_f->GetBinContent(ii)<=0)      { h_vv_f->SetBinContent(ii,0.001);       h_vv_f->SetBinError(ii,0.001); }
    if (h_tt_h_p1_f->GetBinContent(ii)<=0) { h_tt_h_p1_f->SetBinContent(ii,0.001);  h_tt_h_p1_f->SetBinError(ii,0.001); }
    if (h_tt_h_p2_f->GetBinContent(ii)<=0) { h_tt_h_p2_f->SetBinContent(ii,0.001);  h_tt_h_p2_f->SetBinError(ii,0.001); }
    if (h_tt_h_p3_f->GetBinContent(ii)<=0) { h_tt_h_p3_f->SetBinContent(ii,0.001);  h_tt_h_p3_f->SetBinError(ii,0.001); }  

}
 
  // make dir
  TString dirname1 = "particlenet_"+cat+"_"+wp+"_"+era;
  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }

  TString nameoutfile = "particlenet_tt1l_"+cat+"_"+wp+"_"+era+"_"+cutmin+"to"+cutmax+"_templates"; 
  TFile *fout = new TFile("./"+dirname1+"/"+nameoutfile+".root","RECREATE");

  h_data_p->Write("data_obs_pass");
  h_tt_p1_p->Write("tt_p1_pass"); h_tt_p2_p->Write("tt_p2_pass"); h_tt_p3_p->Write("tt_p3_pass");
  h_st_p1_p->Write("st_p1_pass"); h_st_p2_p->Write("st_p2_pass"); h_st_p3_p->Write("st_p3_pass"); 
  h_wqq_p->Write("wqq_pass");
  h_vv_p->Write("vv_pass");

  h_data_f->Write("data_obs_fail");
  h_tt_p1_f->Write("tt_p1_fail"); h_tt_p2_f->Write("tt_p2_fail"); h_tt_p3_f->Write("tt_p3_fail"); 
  h_st_p1_f->Write("st_p1_fail"); h_st_p2_f->Write("st_p2_fail"); h_st_p3_f->Write("st_p3_fail");
  h_wqq_f->Write("wqq_fail");
  h_vv_f->Write("vv_fail");
 

  // systematics  
  h_tt_h_p1_p->Write("tt_herwig_p1_pass"); h_tt_h_p2_p->Write("tt_herwig_p2_pass"); h_tt_h_p3_p->Write("tt_herwig_p3_pass"); 
  h_tt_h_p1_f->Write("tt_herwig_p1_fail"); h_tt_h_p2_f->Write("tt_herwig_p2_fail"); h_tt_h_p3_f->Write("tt_herwig_p3_fail");

  fout->Close();
  std::cout << "\n\n";
}


TH2D *correct_2dext_norm(TH2D *h_nom, TH2D *h_ext) {
  std::cout << h_ext->GetNbinsX() << "  " << h_ext->GetNbinsY() << "\n";
  for (int iy=1; iy<=h_ext->GetNbinsY(); ++iy) {
    float int_nom = h_nom->Integral(1,h_nom->GetNbinsX(),iy,iy);
    float int_ext = h_ext->Integral(1,h_ext->GetNbinsX(),iy,iy);
    float scale = int_nom/int_ext;

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


  TH2D *hTemp = new TH2D(name,name,binsX,minX,maxX,binsY,minY,maxY); //hTemp->SetName(name);
  tree->Project(name,branchY+":"+branchX,cut);

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
