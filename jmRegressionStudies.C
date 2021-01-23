#include "setTDRStyle.h"
#include <sstream>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include <vector>
#include "Rtypes.h"
#include "TColor.h"
#include "TVectorF.h"
#include <cstdlib>
#include "auxiliaryFunctions.h"
#include <string>
#include <unistd.h>


void setTDRStyle();
void jetMassVsTaggerScore(TString name, TTree *tree, TString match, std::vector<TString> wp, TString cut, TString branch, TString xname, TString yname, std::vector<TString> leg_);

void jmRegressionStudies() {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);

  TFile *f_ = new TFile("/eos/uscms/store/user/hqu/NanoTrees/20210109_VH_2018_0L/mc/zh-z2nn-h2cc_tree.root","READONLY");
  TTree *t_ = (TTree*)f_->Get("Events");


  // jet mass vs tagger score
  TString cut = "dr_ak15_hdaus<1.5 && ak15_pt>200. && ak15_pt<300";
  
  std::vector<TString> wp;  wp.clear();
  std::vector<TString> leg; leg.clear();
  wp.push_back("ak15_ParticleNetMD_HccVsQCD>0.99"); leg.push_back("ParticleNet-MD(cc)>0.99");
  wp.push_back("ak15_ParticleNetMD_HccVsQCD>0.96"); leg.push_back("ParticleNet-MD(cc)>0.96");
  wp.push_back("ak15_ParticleNetMD_HccVsQCD>0.90"); leg.push_back("ParticleNet-MD(cc)>0.90");

  TString gen_match = "0==0";
 
  jetMassVsTaggerScore("jms_vs_tagger_particlenet_cc",t_,gen_match,wp,cut,"ak15_regressed_mass","mass [GeV]","a. u.",leg);
 
}


void jetMassVsTaggerScore(TString name, TTree *tree, TString match, std::vector<TString> wp, TString cut, TString branch, TString xname, TString yname, std::vector<TString> leg_) {

  float xmin  = 50.;
  float xmax  = 200.;
  int   xbins = 30.;

  
  TLegend* leg = new TLegend(0.60,0.55,0.94,0.86);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetLineWidth(0);
  
  std::vector<TH1F*> h_; h_.clear();
  for (unsigned int i0=0; i0<wp.size(); ++i0) {
    
    ostringstream tmpIdx; tmpIdx << i0; TString idx = tmpIdx.str();
    TH1F *h_tmp_ = (TH1F*)aux::create1Dhisto("tmp",tree,"1.",cut+" && "+wp[i0]+" && "+match,branch,xbins,xmin,xmax,false,i0+1,i0+1,false,name+"_"+idx,true,false);
    h_.push_back(h_tmp_);
    leg->AddEntry(h_tmp_,leg_[i0],"L");
  }
  
  TCanvas *c_jms = (TCanvas*)aux::makeCanvasWithRatioMC(h_,name,xname,yname,leg);
  c_jms->Print("./plots/jms_"+name+".pdf");
  
  //return c_jms;
}
