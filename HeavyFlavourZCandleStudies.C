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
#include "configuration.h"
//#include <bits/stdc++.h>
#include <string>
#include <unistd.h>

void setTDRStyle();
std::vector<unsigned int> getColors();
std::vector<unsigned int> getStyles();
TH1F *create1Dhisto(TString sample,TTree *tree,TString intLumi,TString cuts,TString branch,int bins,float xmin,float xmax,
                    bool useLog,int color, int style,TString name,bool norm,bool data);
TGraphAsymmErrors *getSFGraph(TString dir,TString algo, TString cat, TString wp, TString era, std::vector<TString> ptrange, int ialgo, int nalgos, int color);
TH1F *createSovSqrtB(TString name, TH1F *h_s, TH1F *h_b);
TH1F *createRatioPassOvFail(TTree *tree, TString name, TString lumi, TString cut, TString addCut, TString var, int bins, float xmin, float xmax, TString xname,
                            int color, int style, bool log, bool norm, bool isdata);
TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, std::vector<TString> vars, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm);
TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm);
TCanvas *c_ratioOf1Dhistos(TString c_name, TString lumi, std::vector<TString> namest, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, 
			   std::vector<TString> namesc, std::vector<TString> cuts,
                           TString xname, TString yname, int nbins, float xmin, float xmax, bool uselog, bool norm);
TH1F *rescaleXaxis(TH1F *inputhisto, float xmin, float xmax);
void rescaleXaxis(TGraphAsymmErrors *inputhisto, double xmin, double scale);
TH1F *getDataMCratio(TGraphAsymmErrors *indata, TH1F *inMC);
void makePlotDataMC(TString name, TString dir,
                    std::vector<TTree*> trees, std::vector<TString> names, TString lumi,TString cut, TString addcut,
                    TString var,int nbins,float xmin,float xmax,TString xaxisname,
                    std::vector<TString> legends, std::vector<int> colors, bool logy, int norm);
void makePlotDataMC(TString name, TString dir,
                    std::vector<TTree*> trees, std::vector<TString> names, TString lumi,std::vector<TString> cuts, TString addcut,
                    TString var,int nbins,float xmin,float xmax,TString xaxisname,
                    std::vector<TString> legends, std::vector<int> colors, bool logy, int norm);
void makeDataMCPlotFromCombine(TString path2file, TString era, TString category, TString category_extra, TString wp,TString name, TString passOrFail, float xmin, float xmax, int nbins,TString xaxisname, bool log, TString sample);
void getSampleComposition(TString path2file, TString filename, TString score, TString year, TString ptrange, TString category);
TH1D *h1DHistoFrom2DTemplates(TString path2file,TString h2dname,TString name,double ymin, double ymax,int color,bool isdata);
TCanvas *makeCanvasWithRatio(std::vector<TH1D*> h1d, TString dir, TString name, TString xname, TString yname, TLegend *leg);
TCanvas *makeCanvasWithRatio(TH1D* hdata, std::vector<TH1D*> h1d, TString dir, TString name, TString xname, TString yname, TLegend *leg);
void makeDataMCFrom2DTemplates(TString path2file, TString nameoutfile, TString name, TString sample, double ymin, double ymax);
void makeDataMCFrom2DTemplatesTop(TString path2file, TString nameoutfile, TString name, TString sample, double ymin, double ymax);
TH1D *h1dShift(TString name, TH1D* inputh1d, int numOfBins, bool up);

TH1D *hrwgt = nullptr;

double rewgtfunc(double x){
  if (!hrwgt) return 1;
  return hrwgt->GetBinContent(hrwgt->FindFixBin(x));
}

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


// sample = tt1l, category = w, wp t or 0p05, .... 
void HeavyFlavourZCandleStudies(TString era, TString sample, TString category, TString wp, bool postfit=false, TString passOrFail="pass", TString category_extra="") {
  conf::configuration(sample);

  //TString era = (TString)era_;
  /*
  TString directory = "particlenet_bb_t_2016";
  getSampleComposition(directory,"particlenet_bb","t","2016","lowpt","p");
  getSampleComposition(directory,"particlenet_bb","t","2016","lowmedpt","p");
  getSampleComposition(directory,"particlenet_bb","t","2016","medpt","p");
  getSampleComposition(directory,"particlenet_bb","t","2016","medhipt","p");
  getSampleComposition(directory,"particlenet_bb","t","2016","hipt","p");
  */

  /*
  std::vector<std::string> name;
  // from ggH
  name.push_back("pt450to500"); name.push_back("pt500to550"); name.push_back("pt550to600"); name.push_back("pt600to675"); name.push_back("pt675to800"); name.push_back("pt800to1200");

  for (int i0=0; i0<name.size(); ++i0) { 
    makeDataMCPlotFromCombine("testarea/fitdir/",name[i0],"pass",40.,201.,23,"m_{SD} [GeV]",false);
    makeDataMCPlotFromCombine("testarea/fitdir/",name[i0],"fail",40.,201.,23,"m_{SD} [GeV]",false);
  }
  */
    
  std::vector<TString> name  = conf::name;
  std::vector<double>  ptmin = conf::ptmin;
  std::vector<double>  ptmax = conf::ptmax;    

  if ( (category=="bb") || (category=="cc") ) {
    /*
    name.push_back("pt200to300"); name.push_back("pt300to400"); name.push_back("pt400to800");
    ptmin.push_back(200.); ptmin.push_back(300.); ptmin.push_back(400.);
    ptmax.push_back(300.); ptmax.push_back(400.); ptmax.push_back(800.);*/
    //conf::name.push_back("pt200to450"); ptmin.push_back(200.); ptmax.push_back(450.);
    
    
    //    name.push_back("pt200to450"); name.push_back("pt450to600"); name.push_back("pt600to800"); name.push_back("pt800to1200");
    //    ptmin.push_back(200.); ptmin.push_back(450.); ptmin.push_back(600.); ptmin.push_back(800.);
    //    ptmax.push_back(450.); ptmax.push_back(600.); ptmax.push_back(800.); ptmax.push_back(1200.);
    //    name.push_back("pt450to600"); ptmin.push_back(450.); ptmax.push_back(600.);
  }

  if ( (category=="w") || (category=="top") ) {
    name.push_back("pt200to300"); name.push_back("pt300to400"); name.push_back("pt400to800");
    ptmin.push_back(200.); ptmin.push_back(300.); ptmin.push_back(400.);
    ptmax.push_back(300.); ptmax.push_back(400.); ptmax.push_back(800.);
    //name.push_back("pt200to450"); name.push_back("pt450to600"); name.push_back("pt600to800"); name.push_back("pt800to1200");
    //ptmin.push_back(200.); ptmin.push_back(450.); ptmin.push_back(600.); ptmin.push_back(800.);
    //ptmax.push_back(450.); ptmax.push_back(600.); ptmax.push_back(800.); ptmax.push_back(1200.);
    //name.push_back("pt450to600"); ptmin.push_back(450.); ptmax.push_back(600.);
  }


  if (postfit) {

    if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") || (sample == "comb")) { 
      std::cout << "postfit: " << sample << " " << category << " " << wp << " " << passOrFail << "\n";
      for (int i0=0; i0<name.size(); ++i0) { 
	//makeDataMCPlotFromCombine("dak8ddt",era,category,wp,name[i0],passOrFail,conf::minX,conf::maxX,conf::binsX,"m_{SD} [GeV]",false,sample);
	makeDataMCPlotFromCombine(conf::algo,era,category,category_extra,wp,name[i0],passOrFail,conf::minX,conf::maxX,conf::binsX,"m_{SD} [GeV]",false,sample); 

      }
    }
    


  }

  if (!postfit) {

  /*
  // from ggH
  name.push_back("pt200toInf"); 
  name.push_back("pt450to500"); name.push_back("pt500to550"); name.push_back("pt550to600"); name.push_back("pt600to675"); name.push_back("pt675to800"); name.push_back("pt800to1200");
  ptmin.push_back(200.);  ptmin.push_back(450.); ptmin.push_back(500.); ptmin.push_back(550.); ptmin.push_back(600.); ptmin.push_back(675.); ptmin.push_back(800.);
  ptmax.push_back(1200.); ptmax.push_back(500.); ptmax.push_back(550.); ptmax.push_back(600.); ptmax.push_back(675.); ptmax.push_back(800.); ptmax.push_back(1200.);
  */  
  //  name.push_back("pt200toInf"); 

  /*
  name.push_back("pt200to300"); name.push_back("pt300to400"); name.push_back("pt400to800");
  ptmin.push_back(200.); ptmin.push_back(300.); ptmin.push_back(400.);
  ptmax.push_back(300.); ptmax.push_back(400.); ptmax.push_back(800.); 
  */

    if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") ) { std::cout << sample << "\n";

      //      makeDataMCFrom2DTemplatesTop("templates2D/particlenet_tt1l_bb_l_2016_200to1200_templates.root","particlenet_tt1l_bb_l_2016","pt200to400","w",200,400);
      for (int i0=0; i0<name.size(); ++i0) {
      //makeDataMCFrom2DTemplatesTop("templates2D/particlenet_tt1l_bb_t_2016_200to1200_templates.root","particlenet_tt1l_bb_t_2016",name[i0],"top",ptmin[i0],ptmax[i0]);
      //      makeDataMCFrom2DTemplatesTop("templates2D/dak8_tt1l_bb_t_2018_200to1200_templates.root","dak8_tt1l_bb_t_2018",name[i0],"top",ptmin[i0],ptmax[i0]);
      //  makeDataMCFrom2DTemplatesTop("templates2D/dak8ddt_tt1l_bb_m_2018_200to1200_templates.root","dak8ddt_tt1l_bb_m_2018",name[i0],"top",ptmin[i0],ptmax[i0]);
	//makeDataMCFrom2DTemplatesTop(sample,"templates2D/particlenet_tt1l_bb_l_2016_200to1200_templates.root","particlenet_tt1l_bb_l_2016",name[i0],"top",ptmin[i0],ptmax[i0]);
	//makeDataMCFrom2DTemplatesTop("templates2D/dak8ddt_"+sample+"_"+category+"_"+wp+"_"+era+"_200to1200_templates.root",
	//			     "dak8ddt_"+sample+"_"+category+"_"+wp+"_"+era,name[i0],category,ptmin[i0],ptmax[i0]);  
	//	makeDataMCFrom2DTemplatesTop("templates2D/dak8md_tt1l_w_t_2016_200to1200_templates.root","dak8md_tt1l_w_t_2016",name[i0],"w",ptmin[i0],ptmax[i0]);
	
	//	makeDataMCFrom2DTemplatesTop("templates2D/dak8ddt_"+sample+"_"+category+"_"+wp+"_"+era+"_200to1200_templates.root",                      
	//			     "dak8ddt_"+sample+"_"+category+"_"+wp+"_"+era,name[i0],category,ptmin[i0],ptmax[i0]);


	makeDataMCFrom2DTemplatesTop("templates2D/"+conf::algo+"_"+sample+"_"+category+"_"+wp+"_"+era+"_200to1200_templates.root",                      
				     conf::algo+"_"+sample+"_"+category+"_"+wp+"_"+era,name[i0],category,ptmin[i0],ptmax[i0]);
      } 
    }
  
    if (sample == "zqq") { std::cout << sample<< "\n";
      for (int i0=0; i0<name.size(); ++i0) {
	//      makeDataMCFrom2DTemplates("particlenet_bb_t_2016/particlenet_bb_t_2016_450to1200_templates.root","particlenet_bb_t_2016",name[i0],ptmin[i0],ptmax[i0]);
	//      makeDataMCFrom2DTemplates("templates2D/particlenet_bb_t_2016_200to1200_templates.root","particlenet_bb_t_2016",name[i0],ptmin[i0],ptmax[i0]);
	makeDataMCFrom2DTemplates("templates2D/"+conf::algo+"_"+sample+"_"+category+"_"+wp+"_"+era+"_200to1200_templates.root",
				  conf::algo+"_"+sample+"_"+category+"_"+wp+"_"+era,name[i0],category,ptmin[i0],ptmax[i0]);
      } 
    }
  }

}

void makeDataMCFrom2DTemplates(TString path2file, TString nameoutfile, TString name, TString sample, double ymin, double ymax) {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
  conf::configuration(path2file);


  // prepare directories
  TString dirname1;
  if (path2file.Contains("tau21ddt"))      { dirname1 = "tau21ddt_sf"; }
  if (path2file.Contains("dak8md"))        { dirname1 = "dak8md_sf"; }
  if (path2file.Contains("dak8ddt"))       { dirname1 = "dak8ddt_sf"; }
  if (path2file.Contains("particlenet"))   { dirname1 = "particlenet_sf"; }
  if (path2file.Contains("particlenetmd")) { dirname1 = "particlenetmd_sf"; }
  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }
  nameoutfile = nameoutfile+"_"+name;

  // MC templates                                                                                                                                                                                                     
  std::vector<TString> syst = {"_","pu","jes","jer","met","lhescalemuf","lhescalemur","lhepdf"};
  //  std::vector<TString> syst = {"_"};                                                                                                                                                                              
  std::vector<TString> processes_in = {"tt_p3","st_p3","ttv_p3","tt_p2","st_p2","ttv_p2","tt_p1","st_p1","ttv_p1","wqq","vv"};


  TH1D *h_data_p = h1DHistoFrom2DTemplates(path2file,"data_obs_pass",name,ymin,ymax,1,true);
  //  TH1D *h_qcd_p  = h1DHistoFrom2DTemplates(path2file,"h_qcd_pass",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_0_p = h1DHistoFrom2DTemplates(path2file,"h_qcd_0_pass",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_1_p = h1DHistoFrom2DTemplates(path2file,"h_qcd_1_pass",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_2_p = h1DHistoFrom2DTemplates(path2file,"h_qcd_2_pass",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_p = (TH1D*)h_qcd_0_p->Clone("h_qcd_p"); h_qcd_p->Add(h_qcd_1_p); h_qcd_p->Add(h_qcd_2_p);
  //h_qcd_p->Smooth();
  TH1D *h_tqq_p  = h1DHistoFrom2DTemplates(path2file,"h_tt_pass",name,ymin,ymax,conf::tqq.color,false);
  TH1D *h_wcx_p  = h1DHistoFrom2DTemplates(path2file,"h_wbosoncx_pass",name,ymin,ymax,conf::wcx.color,false);
  TH1D *h_wll_p  = h1DHistoFrom2DTemplates(path2file,"h_wbosonll_pass",name,ymin,ymax,conf::wll.color,false);
  TH1D *h_zbb_p  = h1DHistoFrom2DTemplates(path2file,"h_zbosonbb_pass",name,ymin,ymax,conf::zbb.color,false);
  TH1D *h_zcc_p  = h1DHistoFrom2DTemplates(path2file,"h_zbosoncc_pass",name,ymin,ymax,conf::zcc.color,false);
  TH1D *h_zll_p  = h1DHistoFrom2DTemplates(path2file,"h_zbosonll_pass",name,ymin,ymax,conf::zll.color,false);

  std::cout << " Data = " << h_data_p->Integral() 
	    << " QCD 0 = " << h_qcd_0_p->Integral() 
	    << " QCD 1 = " << h_qcd_1_p->Integral()
	    << " QCD 2 = " << h_qcd_2_p->Integral()
	    << "\n";

  TH1D *h_data_f = h1DHistoFrom2DTemplates(path2file,"data_obs_fail",name,ymin,ymax,1,true);
  //TH1D *h_qcd_f  = h1DHistoFrom2DTemplates(path2file,"h_qcd_fail",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_0_f  = h1DHistoFrom2DTemplates(path2file,"h_qcd_0_fail",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_1_f  = h1DHistoFrom2DTemplates(path2file,"h_qcd_1_fail",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_2_f  = h1DHistoFrom2DTemplates(path2file,"h_qcd_2_fail",name,ymin,ymax,conf::qcd.color,false);
  TH1D *h_qcd_f = (TH1D*)h_qcd_0_f->Clone("h_qcd_f"); h_qcd_f->Add(h_qcd_1_f); h_qcd_f->Add(h_qcd_2_f);

  TH1D *h_tqq_f  = h1DHistoFrom2DTemplates(path2file,"h_tt_fail",name,ymin,ymax,conf::tqq.color,false);
  TH1D *h_wcx_f  = h1DHistoFrom2DTemplates(path2file,"h_wbosoncx_fail",name,ymin,ymax,conf::wcx.color,false);
  TH1D *h_wll_f  = h1DHistoFrom2DTemplates(path2file,"h_wbosonll_fail",name,ymin,ymax,conf::wll.color,false);
  TH1D *h_zbb_f  = h1DHistoFrom2DTemplates(path2file,"h_zbosonbb_fail",name,ymin,ymax,conf::zbb.color,false);
  TH1D *h_zcc_f  = h1DHistoFrom2DTemplates(path2file,"h_zbosoncc_fail",name,ymin,ymax,conf::zcc.color,false);
  TH1D *h_zll_f  = h1DHistoFrom2DTemplates(path2file,"h_zbosonll_fail",name,ymin,ymax,conf::zll.color,false);

  // add z and w contributions
  TH1D *h_wqq_p = (TH1D*)h_wcx_p->Clone("h_wbosonqq_pass"); h_wqq_p->Add(h_wll_p); h_wqq_p->SetLineColor(conf::wqq.color); h_wqq_p->SetFillColor(conf::wqq.color);
  TH1D *h_wqq_f = (TH1D*)h_wcx_f->Clone("h_wbosonqq_fail"); h_wqq_f->Add(h_wll_f); h_wqq_f->SetLineColor(conf::wqq.color); h_wqq_f->SetFillColor(conf::wqq.color);
  TH1D *h_zqq_p = (TH1D*)h_zbb_p->Clone("h_zbosonqq_pass"); h_zqq_p->Add(h_zcc_p); h_zqq_p->Add(h_zll_p); h_zqq_p->SetLineColor(conf::zqq.color); h_zqq_p->SetFillColor(conf::zqq.color);
  TH1D *h_zqq_f = (TH1D*)h_zbb_f->Clone("h_zbosonqq_fail"); h_zqq_f->Add(h_zcc_f); h_zqq_f->Add(h_zll_f); h_zqq_f->SetLineColor(conf::zqq.color); h_zqq_f->SetFillColor(conf::zqq.color);


  float qcdNormFactor = (h_data_p->Integral()+h_data_f->Integral())/(h_qcd_p->Integral()+h_qcd_f->Integral());
  std::cout << " QCD norm factor = " << qcdNormFactor << "\n";
  qcdNormFactor = 1.;
  h_qcd_p->Scale(qcdNormFactor);
  h_qcd_f->Scale(qcdNormFactor);

  float qcdScaleFactor = (h_data_p->Integral()/(h_data_p->Integral()+h_data_f->Integral()))/(h_qcd_p->Integral()/(h_qcd_p->Integral()+h_qcd_f->Integral()));
  std::cout << " QCD scale factor = " << qcdScaleFactor << "\n";
  qcdScaleFactor = 1.;
  h_qcd_p->Scale(qcdScaleFactor);
  
  //  float tqqNormFactor = 0.85;
  float tqqNormFactor = 1.;
  //  if (name.Contains("pt600to800"))  { tqqNormFactor = 0.85; }
  //  if (name.Contains("pt800to1200")) { tqqNormFactor = 0.80; }
  std::cout << " TQQ norm factor = " << tqqNormFactor << "\n";
  h_tqq_p->Scale(tqqNormFactor);
  h_tqq_f->Scale(tqqNormFactor);

  std::vector<TH1D*> histos_p; histos_p.clear();
  histos_p.push_back(h_data_p); histos_p.push_back(h_qcd_p); histos_p.push_back(h_tqq_p); histos_p.push_back(h_wqq_p); 
  // LG: this part may have to be revisited 
  // depends on the purity and if we need to add syst
  if (nameoutfile.Contains("bb")) { histos_p.push_back(h_zqq_p); }
  if (nameoutfile.Contains("cc")) { 
    histos_p.push_back(h_zqq_p);
    //    histos_p.push_back(h_zbb_p); 
    //    histos_p.push_back(h_zcc_p);
    //    histos_p.push_back(h_zoo_p);
  }
  

  std::vector<TH1D*> histos_f; histos_f.clear();
  histos_f.push_back(h_data_f); histos_f.push_back(h_qcd_f); histos_f.push_back(h_tqq_f); histos_f.push_back(h_wqq_f); 
  // LG: this part may have to be revisited 
  // depends on the purity and if we need to add syst
  if (nameoutfile.Contains("bb")) { histos_f.push_back(h_zqq_f); }
  if (nameoutfile.Contains("cc")) { 
    histos_f.push_back(h_zqq_f);
    //histos_f.push_back(h_zbb_f); 
    //histos_f.push_back(h_zcc_f);
    //histos_f.push_back(h_zoo_f);
  }

  TLegend* leg = new TLegend(0.75,0.55,0.92,0.88);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetLineWidth(0);
  leg->AddEntry(histos_p[0],"Data","PL");
  if (nameoutfile.Contains("bb")) { leg->AddEntry(histos_p[4],"Z","L"); }
  if (nameoutfile.Contains("cc")) {     
    leg->AddEntry(histos_p[4],"Z","L");
    //  leg->AddEntry(histos_p[4],"Z #rightarrow bb","L");
    //leg->AddEntry(histos_p[5],"Z #rightarrow cc","L");
    //leg->AddEntry(histos_p[6],"Z #rightarrow other","L"); 
  }
  leg->AddEntry(histos_p[3],"W","L");
  leg->AddEntry(histos_p[2],"Top","L");
  leg->AddEntry(histos_p[1],"QCD","L");


  TCanvas *c_pass = makeCanvasWithRatio(histos_p,"fit_"+sample+"_templates",name+"_pass","m_{SD} [GeV]","Enetries / bin",leg);  
  c_pass->Print(dirname1+"/"+nameoutfile+"_pass.pdf"); 

  TCanvas *c_fail = makeCanvasWithRatio(histos_f,"fit_"+sample+"_templates",name+"_fail","m_{SD} [GeV]","Enetries / bin",leg);
  c_fail->Print(dirname1+"/"+nameoutfile+"_fail.pdf");

  TCanvas *c_povf = new TCanvas("c_povf_"+name,"c_povf_"+name,500,500);
  TH1D *h_povf_qcd = (TH1D*)h_qcd_p->Clone("h_povf_qcd"); 
  h_povf_qcd->Divide(h_qcd_p,h_qcd_f);
  

  TF1  *f_pol3 = new TF1("f_pol3_"+name,"pol3",50.,250.);
  h_povf_qcd->Fit(f_pol3,"R");
  f_pol3->SetLineWidth(3);
  
  h_povf_qcd->GetYaxis()->SetTitle("Pass / Fail");
  h_povf_qcd->GetXaxis()->SetTitle("m_{SD} [GeV]");
  h_povf_qcd->SetFillColor(0);
  h_povf_qcd->Draw("HIST E0");
  f_pol3->Draw("L sames");
  c_povf->Print(dirname1+"/"+nameoutfile+"_pof_pol3.pdf");


  // correct QCD-pass using the Data-nonQCD in the fail region
  // and multiply by the pol3
  TH1D *h_data_nonQCD_f = (TH1D*)h_data_f->Clone("h_data_nonQCD_f");
  h_data_nonQCD_f->Add(h_tqq_f,-1.); h_data_nonQCD_f->Add(h_wqq_f,-1.); h_data_nonQCD_f->Add(h_zqq_f,-1.);
  h_data_nonQCD_f->SetMarkerSize(0);
  h_data_nonQCD_f->SetLineColor(92); h_data_nonQCD_f->SetFillColor(92);

  TH1D *h_qcd_pred_p = (TH1D*)h_data_nonQCD_f->Clone("h_qcd_pred_p"); 
  h_qcd_pred_p->Multiply(f_pol3);
  //h_qcd_pred_p->Multiply(h_povf_qcd);
  h_qcd_pred_p->SetMarkerSize(0);
  h_qcd_pred_p->SetLineColor(92); h_qcd_pred_p->SetFillColor(92);
  //h_qcd_pred_p->Scale(h_qcd_p->Integral()/h_qcd_pred_p->Integral());

  std::vector<TH1D*> histos_pred_p; histos_pred_p.clear();
  histos_pred_p.push_back(h_data_p); histos_pred_p.push_back(h_qcd_pred_p); histos_pred_p.push_back(h_tqq_p); histos_pred_p.push_back(h_wqq_p); histos_pred_p.push_back(h_zqq_p);
  TCanvas *c_pred_pass = makeCanvasWithRatio(histos_pred_p,"fit_templates",name+"_pred_pass","m_{SD} [GeV]","Enetries / bin",leg);  
  c_pred_pass->Print(dirname1+"/"+nameoutfile+"_pred_pass.pdf");


  // systematics

  // Herwig
  TH1D *h_qcd_h_p_up = h1DHistoFrom2DTemplates(path2file,"h_qcd_herwig_pass",name,ymin,ymax,92,false); h_qcd_h_p_up->SetName("h_qcd_h_p_up");
  TH1D *h_qcd_h_f_up = h1DHistoFrom2DTemplates(path2file,"h_qcd_herwig_fail",name,ymin,ymax,92,false); h_qcd_h_f_up->SetName("h_qcd_h_f_up");

  float qcdHerwigNormFactor = (h_data_p->Integral()+h_data_f->Integral())/(h_qcd_h_p_up->Integral()+h_qcd_h_f_up->Integral());
  std::cout << " QCD herwig norm factor = " << qcdHerwigNormFactor << "\n";
  //  qcdHerwigNormFactor = 1.;
  h_qcd_h_p_up->Scale(qcdHerwigNormFactor);
  h_qcd_h_f_up->Scale(qcdHerwigNormFactor);

  float qcdHerwigScaleFactor = (h_data_p->Integral()/(h_data_p->Integral()+h_data_f->Integral()))/(h_qcd_h_p_up->Integral()/(h_qcd_h_p_up->Integral()+h_qcd_h_f_up->Integral()));
  std::cout << " QCD herwig scale factor = " << qcdHerwigScaleFactor << "\n";
  //qcdHerwigScaleFactor = 1.;
  h_qcd_h_p_up->Scale(qcdHerwigScaleFactor);

  // symmetrize Herwig ; d = QCD- QCD_herwig(Up)
  TH1D *h_qcd_data_diff_p = (TH1D*)h_qcd_p->Clone("h_qcd_data_diff_p"); h_qcd_data_diff_p->Add(h_qcd_h_p_up,-1.);
  TH1D *h_qcd_data_diff_f = (TH1D*)h_qcd_f->Clone("h_qcd_data_diff_f"); h_qcd_data_diff_f->Add(h_qcd_h_f_up,-1.);

  TH1D *h_qcd_h_p_down = (TH1D*)h_qcd_h_p_up->Clone("h_qcd_h_p_down"); h_qcd_h_p_down->Add(h_qcd_data_diff_p,2.);
  TH1D *h_qcd_h_f_down = (TH1D*)h_qcd_h_f_up->Clone("h_qcd_h_f_down"); h_qcd_h_f_down->Add(h_qcd_data_diff_f,2.);

  // correct errors in down
  for (int i0=0; i0<h_qcd_h_p_down->GetNbinsX()+1; ++i0) { h_qcd_h_p_down->SetBinError(i0,h_qcd_h_p_up->GetBinError(i0)); }
  for (int i0=0; i0<h_qcd_h_f_down->GetNbinsX()+1; ++i0) { h_qcd_h_f_down->SetBinError(i0,h_qcd_h_f_up->GetBinError(i0)); }
    

  // jet mass scale
  TH1D *h_zqq_massscale_up_p   = h1dShift("h_zqq_massscale_up_p",h_zqq_p,1,true);  
  TH1D *h_zqq_massscale_down_p = h1dShift("h_zqq_massscale_down_p",h_zqq_p,1,false);  
  TH1D *h_zqq_massscale_up_f   = h1dShift("h_zqq_massscale_up_f",h_zqq_f,1,true);
  TH1D *h_zqq_massscale_down_f = h1dShift("h_zqq_massscale_down_f",h_zqq_f,1,false);

  TH1D *h_zbb_massscale_up_p   = h1dShift("h_zbb_massscale_up_p",h_zbb_p,1,true);  
  TH1D *h_zbb_massscale_down_p = h1dShift("h_zbb_massscale_down_p",h_zbb_p,1,false);  
  TH1D *h_zbb_massscale_up_f   = h1dShift("h_zbb_massscale_up_f",h_zbb_f,1,true);
  TH1D *h_zbb_massscale_down_f = h1dShift("h_zbb_massscale_down_f",h_zbb_f,1,false);

  TH1D *h_zcc_massscale_up_p   = h1dShift("h_zcc_massscale_up_p",h_zcc_p,1,true);  
  TH1D *h_zcc_massscale_down_p = h1dShift("h_zcc_massscale_down_p",h_zcc_p,1,false);  
  TH1D *h_zcc_massscale_up_f   = h1dShift("h_zcc_massscale_up_f",h_zcc_f,1,true);
  TH1D *h_zcc_massscale_down_f = h1dShift("h_zcc_massscale_down_f",h_zcc_f,1,false);

  TH1D *h_zll_massscale_up_p   = h1dShift("h_zll_massscale_up_p",h_zll_p,1,true);  
  TH1D *h_zll_massscale_down_p = h1dShift("h_zll_massscale_down_p",h_zll_p,1,false);  
  TH1D *h_zll_massscale_up_f   = h1dShift("h_zll_massscale_up_f",h_zll_f,1,true);
  TH1D *h_zll_massscale_down_f = h1dShift("h_zll_massscale_down_f",h_zll_f,1,false);
 
  TH1D *h_wqq_massscale_up_p   = h1dShift("h_wqq_massscale_up_p",h_wqq_p,1,true);  
  TH1D *h_wqq_massscale_down_p = h1dShift("h_wqq_massscale_down_p",h_wqq_p,1,false);  
  TH1D *h_wqq_massscale_up_f   = h1dShift("h_wqq_massscale_up_f",h_wqq_f,1,true);
  TH1D *h_wqq_massscale_down_f = h1dShift("h_wqq_massscale_down_f",h_wqq_f,1,false);

  TH1D *h_wcx_massscale_up_p   = h1dShift("h_wcx_massscale_up_p",h_wcx_p,1,true);  
  TH1D *h_wcx_massscale_down_p = h1dShift("h_wcx_massscale_down_p",h_wcx_p,1,false);  
  TH1D *h_wcx_massscale_up_f   = h1dShift("h_wcx_massscale_up_f",h_wcx_f,1,true);
  TH1D *h_wcx_massscale_down_f = h1dShift("h_wcx_massscale_down_f",h_wcx_f,1,false);

  TH1D *h_wll_massscale_up_p   = h1dShift("h_wll_massscale_up_p",h_wll_p,1,true);  
  TH1D *h_wll_massscale_down_p = h1dShift("h_wll_massscale_down_p",h_wll_p,1,false);  
  TH1D *h_wll_massscale_up_f   = h1dShift("h_wll_massscale_up_f",h_wll_f,1,true);
  TH1D *h_wll_massscale_down_f = h1dShift("h_wll_massscale_down_f",h_wll_f,1,false);

   
  // store pass and fail templates to files
  TFile *fout_p = new TFile("./"+dirname1+"/"+nameoutfile+"_templates_p.root","RECREATE");
  h_data_p->Write("data_obs");
  h_qcd_pred_p->Write("qcd");
  //h_qcd_p->Smooth();
  //h_qcd_p->Write("qcd");
  h_tqq_p->Write("tqq"); 
  h_wqq_p->Write("wqq"); 
  h_wcx_p->Write("wcx");
  h_wll_p->Write("wll"); 
  h_zqq_p->Write("zqq");
  h_zbb_p->Write("zbb");
  h_zcc_p->Write("zcc");
  h_zll_p->Write("zll");
  h_qcd_h_p_up->Write("qcd_herwigUp");
  h_qcd_h_p_down->Write("qcd_herwigDown");
  h_zqq_massscale_up_p->Write("zqq_jmsUp");
  h_zqq_massscale_down_p->Write("zqq_jmsDown");
  h_zbb_massscale_up_p->Write("zbb_jmsUp");
  h_zbb_massscale_down_p->Write("zbb_jmsDown");
  h_zcc_massscale_up_p->Write("zcc_jmsUp");
  h_zcc_massscale_down_p->Write("zcc_jmsDown");
  h_zll_massscale_up_p->Write("zll_jmsUp");
  h_zll_massscale_down_p->Write("zll_jmsDown");
  h_zqq_massscale_up_p->Write("zqq_jmspUp"); // was h_w here
  h_zqq_massscale_down_p->Write("zqq_jmspDown");
  h_wqq_massscale_up_p->Write("wqq_wjmsUp");
  h_wqq_massscale_down_p->Write("wqq_wjmsDown");
  h_wqq_massscale_up_p->Write("wqq_wjmspUp");
  h_wqq_massscale_down_p->Write("wqq_wjmspDown");
  h_wqq_massscale_up_p->Write("wqq_jmspUp");
  h_wqq_massscale_down_p->Write("wqq_jmspDown");
  h_wcx_massscale_up_p->Write("wcx_jmspUp");
  h_wcx_massscale_down_p->Write("wcx_jmspDown");
  h_wll_massscale_up_p->Write("wll_jmspUp");
  h_wll_massscale_down_p->Write("wll_jmspDown");
  fout_p->Close();
  
  TFile *fout_f = new TFile("./"+dirname1+"/"+nameoutfile+"_templates_f.root","RECREATE");
  h_data_f->Write("data_obs");
  //h_data_nonQCD_f->Scale(h_qcd_f->Integral()/h_data_nonQCD_f->Integral());
  //h_data_nonQCD_f->Write("qcd");
  h_qcd_f->Write("qcd");
  h_tqq_f->Write("tqq"); 
  h_wqq_f->Write("wqq"); 
  h_wcx_f->Write("wcx");
  h_wll_f->Write("wll"); 
  h_zqq_f->Write("zqq");
  h_zbb_f->Write("zbb");
  h_zcc_f->Write("zcc");
  h_zll_f->Write("zll");
  h_qcd_h_f_up->Write("qcd_herwigUp");
  h_qcd_h_f_down->Write("qcd_herwigDown");
  h_zqq_massscale_up_f->Write("zqq_jmsUp");
  h_zqq_massscale_down_f->Write("zqq_jmsDown");
  h_zbb_massscale_up_f->Write("zbb_jmsUp");
  h_zbb_massscale_down_f->Write("zbb_jmsDown");
  h_zcc_massscale_up_f->Write("zcc_jmsUp");
  h_zcc_massscale_down_f->Write("zcc_jmsDown");
  h_zll_massscale_up_f->Write("zll_jmsUp");
  h_zll_massscale_down_f->Write("zll_jmsDown");
  h_zqq_massscale_up_f->Write("zqq_jmspUp"); // was h_w here
  h_zqq_massscale_down_f->Write("zqq_jmspDown");
  h_wqq_massscale_up_f->Write("wqq_wjmsUp");
  h_wqq_massscale_down_f->Write("wqq_wjmsDown");
  h_wqq_massscale_up_f->Write("wqq_wjmspUp");
  h_wqq_massscale_down_f->Write("wqq_wjmspDown");
  h_wqq_massscale_up_f->Write("wqq_jmspUp");
  h_wqq_massscale_down_f->Write("wqq_jmspDown");
  h_wcx_massscale_up_f->Write("wcx_jmspUp");
  h_wcx_massscale_down_f->Write("wcx_jmspDown");
  h_wll_massscale_up_f->Write("wll_jmspUp");
  h_wll_massscale_down_f->Write("wll_jmspDown");
  fout_f->Close();

  histos_p.clear();
  histos_f.clear();
  histos_pred_p.clear();

  // jet mass resolution                   
  TString inFileNameP = "./"+dirname1+"/"+nameoutfile+"_templates_p.root";
  TString inFileNameF = "./"+dirname1+"/"+nameoutfile+"_templates_f.root";
  system("python jetMassSmearing.py "+inFileNameP+" zqq p");
  system("python jetMassSmearing.py "+inFileNameF+" zqq f");
  system("python jetMassSmearing.py "+inFileNameP+" zbb p");
  system("python jetMassSmearing.py "+inFileNameF+" zbb f");
  system("python jetMassSmearing.py "+inFileNameP+" zcc p");
  system("python jetMassSmearing.py "+inFileNameF+" zcc f");
  system("python jetMassSmearing.py "+inFileNameP+" wqq p");
  system("python jetMassSmearing.py "+inFileNameF+" wqq f");
  system("python jetMassSmearing.py "+inFileNameP+" wcx p");
  system("python jetMassSmearing.py "+inFileNameF+" wcx f");
  system("python jetMassSmearing.py "+inFileNameP+" wll p");
  system("python jetMassSmearing.py "+inFileNameF+" wll f");

}



void makeDataMCFrom2DTemplatesTop(TString path2file, TString nameoutfile, TString name, TString sample, double ymin, double ymax) {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
  conf::configuration(path2file);

  // prepare directories
  TString dirname1;
  if      (path2file.Contains("tau21ddt"))      { dirname1 = "tau21ddt_sf";      }
  else if (path2file.Contains("dak8md"))        { dirname1 = "dak8md_sf";        }
  else if (path2file.Contains("dak8ddt"))       { dirname1 = "dak8ddt_sf";       }
  else if (path2file.Contains("particlenetmd")) { dirname1 = "particlenetmd_sf"; }
  else if (path2file.Contains("particlenet"))   { dirname1 = "particlenet_sf";   }

  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }
  nameoutfile = nameoutfile+"_"+name;


  // Data templates 
  TH1D *h_data_p  = h1DHistoFrom2DTemplates(path2file,"data_obs_pass",name,ymin,ymax,1,true);
  TH1D *h_data_f  = h1DHistoFrom2DTemplates(path2file,"data_obs_fail",name,ymin,ymax,1,true);

  // MC templates
  //  std::vector<TString> syst = {"_","pu","jes","jer","met","lhescalemuf","lhescalemur","lhepdf"};
  //std::vector<TString> syst = {"_","pu","jes","jer","met"};
  std::vector<TString> syst = {"_","pu","jes","jer","met","lhescalemuf","lhescalemur"};
  //  std::vector<TString> syst = {"_"};
  std::vector<TString> processes_in = {"tt_p3","st_p3","ttv_p3","tt_p2","st_p2","ttv_p2","tt_p1","st_p1","ttv_p1","wqq","vv"};

  std::vector<TString> processes_out;
  if (sample == "top") { processes_out.push_back("tqq"); processes_out.push_back("tp2"); processes_out.push_back("tp1"); processes_out.push_back("other"); }
  if (sample == "w")   { processes_out.push_back("tp3"); processes_out.push_back("wqq"); processes_out.push_back("tp1"); processes_out.push_back("other"); }
  if (sample == "bb")  { processes_out.push_back("tqq"); processes_out.push_back("wqq"); processes_out.push_back("tp1"); processes_out.push_back("other"); }
  std::vector<int>     colors    = {conf::tp3.color,conf::tp2.color,conf::tp1.color,conf::other.color};
  std::vector<TString> legends   = {conf::tp3.legend_name,conf::tp2.legend_name,conf::tp1.legend_name,conf::other.legend_name};


  std::vector<TH1D*> hist_out_p; hist_out_p.clear();
  std::vector<TH1D*> hist_out_f; hist_out_f.clear();
  std::vector<TH1D*> hist_out_nom_p; hist_out_nom_p.clear();
  std::vector<TH1D*> hist_out_nom_f; hist_out_nom_f.clear();
  std::vector<TString> category; 
  category.push_back("pass"); 
  category.push_back("fail");
  
  for (unsigned int ic=0; ic<category.size(); ++ic) {
    TString catstr_; if (category[ic] == "pass") { catstr_ = "pass"; } else { catstr_ = "fail"; }
    
    for (unsigned int is=0; is<syst.size(); ++is) {

      int count = 0;
      //      for (unsigned int ip=0; ip<processes_in.size(); ip+=3) {
      for (unsigned int ip=0; ip<processes_in.size(); ip+=1) {
	if (syst[is]=="_") {
	  std::cout << processes_in[ip] << " " << ip << "\n";
	  TH1D *h_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
	  TH1D *h__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
	  h_->Add(h__);
	  if ( (processes_in[ip].Contains("tt_")) || (processes_in[ip].Contains("ttv_")) || (processes_in[ip].Contains("st_")) ) {
	    TH1D *h___ = h1DHistoFrom2DTemplates(path2file,processes_in[ip+2]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
	    h_->Add(h___);
	  }
	  h_->SetName(processes_out[count]); h_->SetLineColor(colors[count]); h_->SetFillColor(colors[count]);
	 
	  if (category[ic] == "pass") { hist_out_p.push_back(h_); hist_out_nom_p.push_back(h_); }
	  else                        { hist_out_f.push_back(h_); hist_out_nom_f.push_back(h_); } 
	}
      
	else {
	  
	  TH1D *h_up_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
	  TH1D *h_up__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
	  h_up_->Add(h_up__);
	  if ( (processes_in[ip].Contains("tt_")) || (processes_in[ip].Contains("ttv_")) || (processes_in[ip].Contains("st_")) ) {
	    TH1D *h_up___ = h1DHistoFrom2DTemplates(path2file,processes_in[ip+2]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
	    h_up_->Add(h_up___); 
	  }
	  h_up_->SetName(processes_out[count]+"_"+syst[is]+"Up"); h_up_->SetLineColor(colors[count]); h_up_->SetFillColor(colors[count]);
	  if (category[ic] == "pass") { hist_out_p.push_back(h_up_); } else { hist_out_f.push_back(h_up_); }

	  TH1D *h_down_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
	  TH1D *h_down__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
	  h_down_->Add(h_down__);
	  if ( (processes_in[ip].Contains("tt_")) || (processes_in[ip].Contains("ttv_")) || (processes_in[ip].Contains("st_")) ) {
	    TH1D *h_down___ = h1DHistoFrom2DTemplates(path2file,processes_in[ip+2]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
	    h_down_->Add(h_down___);
	  }
	  h_down_->SetName(processes_out[count]+"_"+syst[is]+"Down"); h_down_->SetLineColor(colors[count]); h_down_->SetFillColor(colors[count]);
	  if (category[ic] == "pass") { hist_out_p.push_back(h_down_); } else { hist_out_f.push_back(h_down_); }
	}
	if ( (processes_in[ip].Contains("tt_")) || (processes_in[ip].Contains("ttv_")) || (processes_in[ip].Contains("st_")) ) { ip = ip+2; }
	if ( (processes_in[ip].Contains("wqq")) || (processes_in[ip].Contains("vv")) )                                         { ip = ip+1; }
      
      ++count;
      }
    } // end of looping over the syst 
  } // end of category loop 


  
  // Additional systematics 
  /*
  // jet mass scale
  TH1D *h_tt_p3_massscale_up_p   = h1dShift("h_tt_p3_massscale_up_p",hist_out_p[0],1,true);    h_tt_p3_massscale_up_p->SetName(processes_out[0]+"_p3jmspUp");      
  TH1D *h_tt_p3_massscale_down_p = h1dShift("h_tt_p3_massscale_down_p",hist_out_p[0],1,false); h_tt_p3_massscale_down_p->SetName(processes_out[0]+"_p3jmspDown"); 
  TH1D *h_tt_p3_massscale_up_f   = h1dShift("h_tt_p3_massscale_up_f",hist_out_f[0],1,true);    h_tt_p3_massscale_up_f->SetName(processes_out[0]+"_p3jmsfUp");
  TH1D *h_tt_p3_massscale_down_f = h1dShift("h_tt_p3_massscale_down_f",hist_out_f[0],1,false); h_tt_p3_massscale_down_f->SetName(processes_out[0]+"_p3jmsfDown");

  TH1D *h_tt_p2_massscale_up_p   = h1dShift("h_tt_p2_massscale_up_p",hist_out_p[1],1,true);    h_tt_p2_massscale_up_p->SetName(processes_out[1]+"_vjmspUp"); 
  TH1D *h_tt_p2_massscale_down_p = h1dShift("h_tt_p2_massscale_down_p",hist_out_p[1],1,false); h_tt_p2_massscale_down_p->SetName(processes_out[1]+"_vjmspDown");
  TH1D *h_tt_p2_massscale_up_f   = h1dShift("h_tt_p2_massscale_up_f",hist_out_f[1],1,true);    h_tt_p2_massscale_up_f->SetName(processes_out[1]+"_vjmsfUp");
  TH1D *h_tt_p2_massscale_down_f = h1dShift("h_tt_p2_massscale_down_f",hist_out_f[1],1,false); h_tt_p2_massscale_down_f->SetName(processes_out[1]+"_vjmsfDown");

  TH1D *h_tt_p1_massscale_up_p   = h1dShift("h_tt_p1_massscale_up_p",hist_out_p[2],1,true);    h_tt_p1_massscale_up_p->SetName(processes_out[2]+"_p1jmspUp");  
  TH1D *h_tt_p1_massscale_down_p = h1dShift("h_tt_p1_massscale_down_p",hist_out_p[2],1,false); h_tt_p1_massscale_down_p->SetName(processes_out[2]+"_p1jmspDown");
  TH1D *h_tt_p1_massscale_up_f   = h1dShift("h_tt_p1_massscale_up_f",hist_out_f[2],1,true);    h_tt_p1_massscale_up_f->SetName(processes_out[2]+"_p1jmsfUp");
  TH1D *h_tt_p1_massscale_down_f = h1dShift("h_tt_p1_massscale_down_f",hist_out_f[2],1,false); h_tt_p1_massscale_down_f->SetName(processes_out[2]+"_p1jmsfDown");

  TH1D *h_tt_other_massscale_up_p   = h1dShift("h_tt_other_massscale_up_p",hist_out_p[3],1,true);    h_tt_other_massscale_up_p->SetName(processes_out[3]+"_ojmspUp"); 
  TH1D *h_tt_other_massscale_down_p = h1dShift("h_tt_other_massscale_down_p",hist_out_p[3],1,false); h_tt_other_massscale_down_p->SetName(processes_out[3]+"_ojmspDown");
  TH1D *h_tt_other_massscale_up_f   = h1dShift("h_tt_other_massscale_up_f",hist_out_f[3],1,true);    h_tt_other_massscale_up_f->SetName(processes_out[3]+"_ojmsfUp");
  TH1D *h_tt_other_massscale_down_f = h1dShift("h_tt_other_massscale_down_f",hist_out_f[3],1,false); h_tt_other_massscale_down_f->SetName(processes_out[3]+"_ojmsfDown");
  

  // store them
  hist_out_p.push_back(h_tt_p3_massscale_up_p);
  hist_out_p.push_back(h_tt_p3_massscale_down_p);
  hist_out_p.push_back(h_tt_p2_massscale_up_p);
  hist_out_p.push_back(h_tt_p2_massscale_down_p);
  hist_out_p.push_back(h_tt_p1_massscale_up_p);
  hist_out_p.push_back(h_tt_p1_massscale_down_p);
  hist_out_p.push_back(h_tt_other_massscale_up_p);
  hist_out_p.push_back(h_tt_other_massscale_down_p);
  
  hist_out_f.push_back(h_tt_p3_massscale_up_f);
  hist_out_f.push_back(h_tt_p3_massscale_down_f);
  hist_out_f.push_back(h_tt_p2_massscale_up_f);
  hist_out_f.push_back(h_tt_p2_massscale_down_f);
  hist_out_f.push_back(h_tt_p1_massscale_up_f);
  hist_out_f.push_back(h_tt_p1_massscale_down_f);
  hist_out_f.push_back(h_tt_other_massscale_up_f);
  hist_out_f.push_back(h_tt_other_massscale_down_f);
  */

  


  /*
  // Additional systematics: Herwig
  TH1D *h_tt_herwig_p1_p_up = h1DHistoFrom2DTemplates(path2file,"tt_herwig_p1_pass",name,ymin,ymax,4,false); h_tt_herwig_p1_p_up->SetName("h_tt_herwig_p1_p_up"); 
  TH1D *h_tt_herwig_p2_p_up = h1DHistoFrom2DTemplates(path2file,"tt_herwig_p2_pass",name,ymin,ymax,4,false); h_tt_herwig_p2_p_up->SetName("h_tt_herwig_p2_p_up");
  TH1D *h_tt_herwig_p3_p_up = h1DHistoFrom2DTemplates(path2file,"tt_herwig_p3_pass",name,ymin,ymax,4,false); h_tt_herwig_p3_p_up->SetName("h_tt_herwig_p3_p_up");
  h_tt_herwig_p1_p_up->Add(h_st_p1_p);
  h_tt_herwig_p2_p_up->Add(h_st_p2_p);
  h_tt_herwig_p3_p_up->Add(h_st_p3_p);

  TH1D *h_tt_herwig_p1_f_up = h1DHistoFrom2DTemplates(path2file,"tt_herwig_p1_fail",name,ymin,ymax,4,false); h_tt_herwig_p1_f_up->SetName("h_tt_herwig_p1_f_up"); 
  TH1D *h_tt_herwig_p2_f_up = h1DHistoFrom2DTemplates(path2file,"tt_herwig_p2_fail",name,ymin,ymax,4,false); h_tt_herwig_p2_f_up->SetName("h_tt_herwig_p2_f_up");
  TH1D *h_tt_herwig_p3_f_up = h1DHistoFrom2DTemplates(path2file,"tt_herwig_p3_fail",name,ymin,ymax,4,false); h_tt_herwig_p3_f_up->SetName("h_tt_herwig_p3_f_up");
  h_tt_herwig_p1_f_up->Add(h_st_p1_f);
  h_tt_herwig_p2_f_up->Add(h_st_p2_f);
  h_tt_herwig_p3_f_up->Add(h_st_p3_f);

  // symmetrize Herwig ; d = tt- tt_herwig(Up)
  TH1D *h_tt_p1_diff_p = (TH1D*)h_top_p1_p->Clone("h_tt_p1_diff_p"); h_tt_p1_diff_p->Add(h_tt_herwig_p1_p_up,-1.);
  TH1D *h_tt_p2_diff_p = (TH1D*)h_top_p2_p->Clone("h_tt_p2_diff_p"); h_tt_p2_diff_p->Add(h_tt_herwig_p2_p_up,-1.);
  TH1D *h_tt_p3_diff_p = (TH1D*)h_top_p3_p->Clone("h_tt_p3_diff_p"); h_tt_p3_diff_p->Add(h_tt_herwig_p3_p_up,-1.);

  TH1D *h_tt_p1_diff_f = (TH1D*)h_top_p1_f->Clone("h_tt_p1_diff_f"); h_tt_p1_diff_f->Add(h_tt_herwig_p1_f_up,-1.);
  TH1D *h_tt_p2_diff_f = (TH1D*)h_top_p2_f->Clone("h_tt_p2_diff_f"); h_tt_p2_diff_f->Add(h_tt_herwig_p2_f_up,-1.);
  TH1D *h_tt_p3_diff_f = (TH1D*)h_top_p3_f->Clone("h_tt_p3_diff_f"); h_tt_p3_diff_f->Add(h_tt_herwig_p3_f_up,-1.);

  TH1D *h_tt_herwig_p1_p_down = (TH1D*)h_tt_herwig_p1_p_up->Clone("h_tt_herwig_p1_p_down"); h_tt_herwig_p1_p_down->Add(h_tt_p1_diff_p,2.);
  TH1D *h_tt_herwig_p2_p_down = (TH1D*)h_tt_herwig_p2_p_up->Clone("h_tt_herwig_p2_p_down"); h_tt_herwig_p2_p_down->Add(h_tt_p2_diff_p,2.);
  TH1D *h_tt_herwig_p3_p_down = (TH1D*)h_tt_herwig_p3_p_up->Clone("h_tt_herwig_p3_p_down"); h_tt_herwig_p3_p_down->Add(h_tt_p3_diff_p,2.);  

  TH1D *h_tt_herwig_p1_f_down = (TH1D*)h_tt_herwig_p1_f_up->Clone("h_tt_herwig_p1_f_down"); h_tt_herwig_p1_f_down->Add(h_tt_p1_diff_f,2.);
  TH1D *h_tt_herwig_p2_f_down = (TH1D*)h_tt_herwig_p2_f_up->Clone("h_tt_herwig_p2_f_down"); h_tt_herwig_p2_f_down->Add(h_tt_p2_diff_f,2.);
  TH1D *h_tt_herwig_p3_f_down = (TH1D*)h_tt_herwig_p3_f_up->Clone("h_tt_herwig_p3_f_down"); h_tt_herwig_p3_f_down->Add(h_tt_p3_diff_f,2.);  

  // correct errors in down
  for (int i0=0; i0<h_tt_herwig_p1_p_down->GetNbinsX()+1; ++i0) { 
    h_tt_herwig_p1_p_down->SetBinError(i0,h_tt_herwig_p1_p_up->GetBinError(i0)); 
    h_tt_herwig_p2_p_down->SetBinError(i0,h_tt_herwig_p2_p_up->GetBinError(i0));
    h_tt_herwig_p3_p_down->SetBinError(i0,h_tt_herwig_p3_p_up->GetBinError(i0));
  }

  for (int i0=0; i0<h_tt_herwig_p1_f_down->GetNbinsX()+1; ++i0) { 
    h_tt_herwig_p1_f_down->SetBinError(i0,h_tt_herwig_p1_f_up->GetBinError(i0)); 
    h_tt_herwig_p2_f_down->SetBinError(i0,h_tt_herwig_p2_f_up->GetBinError(i0));
    h_tt_herwig_p3_f_down->SetBinError(i0,h_tt_herwig_p3_f_up->GetBinError(i0));
  }
  */  
  
  
  TLegend* leg = new TLegend(0.75,0.55,0.92,0.88);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetLineWidth(0);
  leg->AddEntry(h_data_p,"Data","PL");
  for (unsigned int i0=0; i0<hist_out_nom_p.size(); ++i0) { leg->AddEntry(hist_out_nom_p[i0],legends[i0],"L"); }

  TCanvas *c_pass = makeCanvasWithRatio(h_data_p,hist_out_nom_p,"fit_"+sample+"_templates",name+"_pass","m_{SD} [GeV]","Enetries / bin",leg);  
  c_pass->Print(dirname1+"/"+nameoutfile+"_pass.pdf"); 

  TCanvas *c_fail = makeCanvasWithRatio(h_data_f,hist_out_nom_f,"fit_"+sample+"_templates",name+"_fail","m_{SD} [GeV]","Enetries / bin",leg);
  c_fail->Print(dirname1+"/"+nameoutfile+"_fail.pdf");
  

  TFile *fout_p = new TFile("./"+dirname1+"/"+nameoutfile+"_templates_p.root","RECREATE");
  h_data_p->Write("data_obs");
  for (unsigned int i0=0; i0<hist_out_p.size(); ++i0) {
    hist_out_p[i0]->Write(hist_out_p[i0]->GetName()); std::cout << hist_out_p[i0]->GetName() << "\n";
  }
  // jms for softdrop corrections
  //  h_tt_p2_massscale_up_p->Write(processes_out[1]+"_vjmsUp");
  //  h_tt_p2_massscale_down_p->Write(processes_out[1]+"_vjmsDown");
  
  fout_p->Close();

  TFile *fout_f = new TFile("./"+dirname1+"/"+nameoutfile+"_templates_f.root","RECREATE");
  h_data_f->Write("data_obs");
  for (unsigned int i0=0; i0<hist_out_f.size(); ++i0) {
    hist_out_f[i0]->Write(hist_out_f[i0]->GetName());
  }
  // jms for softdrop corrections
  //h_tt_p2_massscale_up_f->Write(processes_out[1]+"_vjmsUp");
  //h_tt_p2_massscale_down_f->Write(processes_out[1]+"_vjmsDown");
  fout_f->Close();
  
  hist_out_p.clear(); hist_out_nom_p.clear();
  hist_out_f.clear(); hist_out_nom_f.clear();

  // jet mass resolution
  TString inFileNameP = "./"+dirname1+"/"+nameoutfile+"_templates_p.root";
  TString inFileNameF = "./"+dirname1+"/"+nameoutfile+"_templates_f.root";    
  //system("python jetMassSmearing.py "+inFileNameP+" "+processes_out[1]+" p");
  //system("python jetMassSmearing.py "+inFileNameF+" "+processes_out[1]+" f");

  // JMS and JMR
  TString smearVal = "0.15";
  TString scaleVal = "5."; // GeV
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[0]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 1 --suffix p3");
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[0]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF -1 --suffix p3");
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[1]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 1 --suffix v");
  system("python jetMassSandR.py --ifile "+inFileNameF+" --sample "+processes_out[1]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF -1 --suffix v");
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[2]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 1 --suffix p1");
  system("python jetMassSandR.py --ifile "+inFileNameF+" --sample "+processes_out[2]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF -1 --suffix p1");
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[3]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 1 --suffix o");
  system("python jetMassSandR.py --ifile "+inFileNameF+" --sample "+processes_out[3]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF -1 --suffix o");

  // m(SD) corrections
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[0]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 0 --suffix p3");
  system("python jetMassSandR.py --ifile "+inFileNameF+" --sample "+processes_out[0]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 0 --suffix p3");
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[1]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 0 --suffix v");
  system("python jetMassSandR.py --ifile "+inFileNameF+" --sample "+processes_out[1]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 0 --suffix v");
  system("python jetMassSandR.py --ifile "+inFileNameP+" --sample "+processes_out[2]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 0 --suffix p1");
  system("python jetMassSandR.py --ifile "+inFileNameF+" --sample "+processes_out[2]+" --smear "+smearVal+" --scale "+scaleVal+" --pOrF 0 --suffix p1");

  fout_p->Close();
  fout_f->Close();
  


}


TCanvas *makeCanvasWithRatio(std::vector<TH1D*> h1d, TString dir, TString name, TString xname, TString yname, TLegend *leg) {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);

  THStack *hs = new THStack("hs_"+name,"hs_"+name);
  for (int i0=1; i0<h1d.size(); ++i0) { hs->Add(h1d[i0]); }

  TH1D *h_tqq_lin = (TH1D*)h1d[2]->Clone("h_tqq_lin_"+name); h_tqq_lin->SetFillColor(0); h_tqq_lin->SetLineStyle(2); 
  TH1D *h_wqq_lin = (TH1D*)h1d[3]->Clone("h_wqq_lin_"+name); h_wqq_lin->SetFillColor(0); h_wqq_lin->SetLineStyle(2);
  TH1D *h_zqq_lin = (TH1D*)h1d[4]->Clone("h_zqq_lin_"+name); h_zqq_lin->SetFillColor(0); h_zqq_lin->SetLineStyle(2);  

  TH1D *h_sm = (TH1D*)h1d[1]->Clone("h_sm");
  h_sm->Reset("ICES");



  //  TString logstr = "lin"; if (logy) { logstr = "log"; }

  TCanvas *c_ = new TCanvas("c_"+name,"c_"+name,600,600); c_->SetName("c_"+name);
  TPad *pMain  = new TPad("pMain_"+name,"pMain+"+name,0.0,0.25,1.0,1.0);
  pMain->SetRightMargin(0.05);
  pMain->SetLeftMargin(0.17);
  pMain->SetBottomMargin(0.);
  pMain->SetTopMargin(0.1);
  TPad *pRatio = new TPad("pRatio_"+name,"pRatio_"+name,0.0,0.0,1.0,0.25);
  pRatio->SetRightMargin(0.05);
  pRatio->SetLeftMargin(0.17);
  pRatio->SetTopMargin(0.);
  pRatio->SetBottomMargin(0.37);
  pMain->Draw();
  pRatio->Draw();
  
  pMain->cd();
  for (int i0=0; i0<h1d.size(); ++i0) {
    if (i0==0) { 
      h1d[i0]->GetXaxis()->SetTitle(xname);
      h1d[i0]->GetYaxis()->SetTitle(yname);
      
      if (name.Contains("pass")) { h1d[i0]->GetYaxis()->SetRangeUser(0.,1.5*h1d[i0]->GetMaximum()); }
      else { h1d[i0]->GetYaxis()->SetRangeUser(0.,1.2*h1d[i0]->GetMaximum()); }
      h1d[i0]->Draw("P E0");
    }
    else { 
      h1d[i0]->Draw("HIST E0 sames");
      h_sm->Add(h1d[i0]);
    }
  }


  // data mc ratio
  TH1F *h_r = (TH1F*)h1d[0]->Clone("h_"+name+"_r"); h_r->Divide(h1d[0],h_sm);

  /*
  TList *histos = hs->GetHists();
  TIter next(histos);
  TH1F *hist;
  while ((hist =(TH1F*)next())) { hist->Scale(normVal); }
  */
  // continue drawing
  hs->Draw("HIST E0 sames");
  h_tqq_lin->Draw("HIST E0 sames");
  h_wqq_lin->Draw("HIST E0 sames");
  h_zqq_lin->Draw("HIST E0 sames");
  h1d[0]->Draw("P E0 sames");
  leg->Draw("sames");
  c_->RedrawAxis();
  pMain->RedrawAxis();

  pRatio->cd();
  gPad->SetBottomMargin(0.2);
  h_r->GetYaxis()->SetTitleOffset(0.9);
  h_r->GetYaxis()->SetTitleSize(0.1);
  h_r->GetYaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitleSize(0.1);
  h_r->GetXaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitle(xname);
  h_r->GetYaxis()->SetTitle("Data / MC");
  h_r->GetYaxis()->SetRangeUser(0.6,1.3);
  h_r->Draw("P E0");

  /*
  const int dir_err = system("mkdir -p ./"+dir);
  if (-1 == dir_err) {
    printf("Error creating directory!n");
    exit(1);
  }
  // c_->Print(dir+"/"+name+"_"+logstr+".png");
  // c_->Print(dir+"/"+name+"_"+logstr+".pdf");
  */
  return c_;
}


TCanvas *makeCanvasWithRatio(TH1D* hdata, std::vector<TH1D*> h1d, TString dir, TString name, TString xname, TString yname, TLegend *leg) {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);

  THStack *hs = new THStack("hs_"+name,"hs_"+name);
  std::vector<TH1D*> histos; histos.clear();
  for (int i0=0; i0<h1d.size(); ++i0) { 
    TString count = std::to_string(i0);
    hs->Add(h1d[i0]);
    TH1D* h_ = (TH1D*)h1d[i0]->Clone("h_"+count); h_->SetFillColor(0); h_->SetLineStyle(2);
    histos.push_back(h_);
  }

  TH1D *h_sm = (TH1D*)h1d[0]->Clone("h_sm");
  h_sm->Reset("ICES");


  //  TString logstr = "lin"; if (logy) { logstr = "log"; }
  TCanvas *c_ = new TCanvas("c_"+name,"c_"+name,600,600); c_->SetName("c_"+name);
  TPad *pMain  = new TPad("pMain_"+name,"pMain+"+name,0.0,0.25,1.0,1.0);
  pMain->SetRightMargin(0.05);
  pMain->SetLeftMargin(0.17);
  pMain->SetBottomMargin(0.);
  pMain->SetTopMargin(0.1);
  TPad *pRatio = new TPad("pRatio_"+name,"pRatio_"+name,0.0,0.0,1.0,0.25);
  pRatio->SetRightMargin(0.05);
  pRatio->SetLeftMargin(0.17);
  pRatio->SetTopMargin(0.);
  pRatio->SetBottomMargin(0.37);
  pMain->Draw();
  pRatio->Draw();
  
  pMain->cd();
  hdata->GetXaxis()->SetTitle(xname);
  hdata->GetYaxis()->SetTitle(yname);
  if (name.Contains("pass")) { hdata->GetYaxis()->SetRangeUser(0.,1.5*hdata->GetMaximum()); }
  else { hdata->GetYaxis()->SetRangeUser(0.,1.2*hdata->GetMaximum()); }
  hdata->Draw("P E0");

  for (int i0=0; i0<h1d.size(); ++i0) { h1d[i0]->Draw("HIST E0 sames"); h_sm->Add(h1d[i0]); }


  // data mc ratio
  TH1F *h_r = (TH1F*)hdata->Clone("h_"+name+"_r"); h_r->Divide(hdata,h_sm);

  /*
  TList *histos = hs->GetHists();
  TIter next(histos);
  TH1F *hist;
  while ((hist =(TH1F*)next())) { hist->Scale(normVal); }
  */
  // continue drawing
  hs->Draw("HIST E0 sames");
  for (unsigned int i0=0; i0<histos.size(); ++i0) { histos[i0]->Draw("HIST E0 sames"); }
  hdata->Draw("P E0 sames");
  leg->Draw("sames");
  c_->RedrawAxis();
  pMain->RedrawAxis();

  pRatio->cd();
  gPad->SetBottomMargin(0.2);
  h_r->GetYaxis()->SetTitleOffset(0.9);
  h_r->GetYaxis()->SetTitleSize(0.1);
  h_r->GetYaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitleSize(0.1);
  h_r->GetXaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitle(xname);
  h_r->GetYaxis()->SetTitle("Data / MC");
  h_r->GetYaxis()->SetRangeUser(0.6,1.3);
  h_r->Draw("P E0");

  /*
  const int dir_err = system("mkdir -p ./"+dir);
  if (-1 == dir_err) {
    printf("Error creating directory!n");
    exit(1);
  }
  // c_->Print(dir+"/"+name+"_"+logstr+".png");
  // c_->Print(dir+"/"+name+"_"+logstr+".pdf");
  */
  return c_;
}


TH1D *h1DHistoFrom2DTemplates(TString path2file,TString h2dname, TString name, double ymin, double ymax,int color,bool isdata) {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);


  TFile *f2d = TFile::Open(path2file,"READONLY");
  TH2D  *h2d = (TH2D*)f2d->Get(h2dname); h2d->SetDirectory(0);
  TH1D  *h1d = h2d->ProjectionX("h1d_"+h2dname,h2d->GetYaxis()->FindBin(ymin),h2d->GetYaxis()->FindBin(ymax),"e"); h1d->SetDirectory(0);
  //TH1D  *h1d = h2d->ProfileX("h1d_"+h2dname,h2d->GetYaxis()->FindBin(ymin),h2d->GetYaxis()->FindBin(ymax),"e");
  h1d->SetName(h2dname+"_"+name);
  //h1d->Rebin(2);
  
  h1d->SetLineColor(color);
  if (isdata) {
    h1d->SetFillColor(0); 
    h1d->SetMarkerStyle(20); h1d->SetMarkerSize(1.2); 
    h1d->SetLineWidth(3); h1d->SetLineColor(color);
  }
  else {
    h1d->SetFillColor(color);
    h1d->SetMarkerStyle(20); h1d->SetMarkerSize(0.);
    h1d->SetLineWidth(3); h1d->SetLineColor(color);
  }

  f2d->Close();

  return h1d;
}

/*
TH1D *h2DHistoFrom2DTemplates(TString path2file,TString h2dname, TString name, double ymin, double ymax,int color,bool isdata) {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);

  TFile *f2d = TFile::Open(path2file,"READONLY");
  TH2D  *h2d = (TH2D*)f2d->Get(h2dname);
  TH2D  *h2d = h2d->ProjectionX("h1d_"+h2dname,h2d->GetYaxis()->FindBin(ymin),h2d->GetYaxis()->FindBin(ymax),"e");
  h1d->SetName(h2dname+"_"+name);
  //h1d->Rebin(2);
  
  h1d->SetLineColor(color);
  if (isdata) {
    h1d->SetFillColor(0); 
    h1d->SetMarkerStyle(20); h1d->SetMarkerSize(1.2); 
    h1d->SetLineWidth(3); h1d->SetLineColor(color);
  }
  else {
    h1d->SetFillColor(color);
    h1d->SetMarkerStyle(20); h1d->SetMarkerSize(0.);
    h1d->SetLineWidth(3); h1d->SetLineColor(color);
  }

  return h1d;
}
*/

TH1D *h1dShift(TString name, TH1D* inputh1d, int numOfBins, bool up) {
  // need to be fixed for case with numOfBins>1 
 
  TH1D *hshifted = (TH1D*)inputh1d->Clone("hshifted_"+name);
  
  //hshifted->Reset("ICES");
  hshifted->SetName("hshifted_"+name);
  
  if (up) { 
    for (int i0=0; i0<=inputh1d->GetNbinsX()+1; ++i0) {
      if (i0<numOfBins) {
	hshifted->SetBinContent(i0,inputh1d->GetBinContent(i0));
	hshifted->SetBinError(i0,inputh1d->GetBinError(i0));
      }
      hshifted->SetBinContent(i0,inputh1d->GetBinContent(i0-numOfBins));
      hshifted->SetBinError(i0,inputh1d->GetBinError(i0-numOfBins));
    }
  }
  
  else {
    for (int i0=0; i0<=inputh1d->GetNbinsX()+1-numOfBins; ++i0) {
      hshifted->SetBinContent(i0,inputh1d->GetBinContent(i0+numOfBins));
      hshifted->SetBinError(i0,inputh1d->GetBinError(i0+numOfBins));
    }
  }

  // fix values of first bin
  if (up) { hshifted->SetBinContent(1,hshifted->GetBinContent(2)+0.0001); }

  return hshifted;
}


void compTemplateShapes() {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);

  TFile *f_w       = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200423/w-qq_tree.root" , "READONLY");
  TFile *f_w_ht200 = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200423/w-qq_ht200_tree.root" , "READONLY");
  TFile *f_z       = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200423/z-qq_tree.root" , "READONLY");
  TFile *f_z_ht200 = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200423/z-qq_ht200_tree.root" , "READONLY");
  TFile *f_h       = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200423/ggH-bb_tree.root" , "READONLY");
  TFile *f_top     = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200423/top_tree.root" , "READONLY");
  TFile *f_qcd     = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200423/qcd-mg_tree.root" , "READONLY");
  
  TFile *f_hb   = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200330/vhtobb_tree.root" , "READONLY");

  TTree *t_w       = (TTree*)f_w->Get("Events");
  TTree *t_w_ht200 = (TTree*)f_w_ht200->Get("Events");
  TTree *t_z       = (TTree*)f_z->Get("Events");
  TTree *t_z_ht200 = (TTree*)f_z_ht200->Get("Events");
  TTree *t_h       = (TTree*)f_h->Get("Events");
  TTree *t_top     = (TTree*)f_top->Get("Events");
  TTree *t_qcd     = (TTree*)f_qcd->Get("Events");
  
  TString c_inc = "(abs(fj_1_eta)<2.4 && fj_1_sdmass>30. && fj_1_sdmass<250. && ht>200. && fj_1_pt>400. && fj_1_ParticleNetMD_XbbVsQCD>0.)";

  std::vector<TTree*>  trees;  trees.clear();
  std::vector<TString> names;  names.clear();
  std::vector<bool>    isdata; isdata.clear();
  std::vector<int>     colors; colors.clear();
  std::vector<int>     styles; styles.clear();
  std::vector<TString> cuts;   cuts.clear();
  
  // bb working points
  std::vector<TString> c_bb; c_bb.clear();
  /*c_bb.push_back("fj_1_ParticleNetMD_XbbVsQCD>0.40");
    c_bb.push_back("fj_1_ParticleNetMD_XbbVsQCD>0.54"); c_bb.push_back("fj_1_ParticleNetMD_XbbVsQCD>0.84"); c_bb.push_back("fj_1_ParticleNetMD_XbbVsQCD>0.93");*/
  c_bb.push_back("fj_1_ParticleNetMD_XbbVsQCD<0.54");
  c_bb.push_back("fj_1_ParticleNetMD_XbbVsQCD<0.93 && fj_1_pt>600.");
  //c_bb.push_back("fj_1_pt>200 && fj_1_pt<300"); c_bb.push_back("fj_1_pt>600 && fj_1_pt<800"); c_bb.push_back("fj_1_pt>800 && fj_1_pt<1200");
  //c_bb.push_back("fj_1_pt>400 && fj_1_pt<500 && ht>200 && fj_1_ParticleNetMD_XbbVsQCD>0.93"); c_bb.push_back("fj_1_pt>400 && fj_1_pt<500 && ht>800 && fj_1_ParticleNetMD_XbbVsQCD>0.93");
  
  //trees.push_back(t_w_ht200); names.push_back("wboson"); isdata.push_back(false); colors.push_back(1); styles.push_back(1); cuts.push_back(c_inc+"&& fj_1_nbhadrons>=0 && (fj_1_isW<0.8)");
  //trees.push_back(t_z_ht200); names.push_back("zboson"); isdata.push_back(false); colors.push_back(1); styles.push_back(1); cuts.push_back(c_inc+"&& fj_1_nbhadrons>=0 && (fj_1_isZ<0.8)");


  trees.push_back(t_w_ht200); names.push_back("wboson"); isdata.push_back(false); colors.push_back(6); styles.push_back(1); cuts.push_back(c_inc+"&& fj_1_nbhadrons>=0 && (fj_1_isW<0.8)");
  trees.push_back(t_z_ht200); names.push_back("zboson"); isdata.push_back(false); colors.push_back(8); styles.push_back(2); cuts.push_back(c_inc+"&& fj_1_nbhadrons>=2 && (fj_1_isZ<0.8)");
  trees.push_back(t_h); names.push_back("higgs"); isdata.push_back(false); colors.push_back(2); styles.push_back(4); cuts.push_back(c_inc+"&& fj_1_nbhadrons>=2 && (fj_1_isH<0.8)");
  trees.push_back(t_top); names.push_back("top"); isdata.push_back(false); colors.push_back(4); styles.push_back(4); cuts.push_back(c_inc+"&& fj_1_nbhadrons>=0");
  trees.push_back(t_qcd); names.push_back("qcd"); isdata.push_back(false); colors.push_back(92); styles.push_back(1); cuts.push_back(c_inc+"&& fj_1_nbhadrons>=0");
  
  for (int it=0; it<trees.size(); ++it) {

    std::vector<TTree*> tr; tr.clear(); std::vector<TString> nm; nm.clear(); std::vector<bool> isd; isd.clear(); 
    std::vector<int> cl; cl.clear(); std::vector<int> ss; ss.clear(); 
    std::vector<TString> csp; csp.clear(); 
    std::vector<TString> csf; csf.clear(); 
    for (int ic=0; ic<c_bb.size(); ++ic) {
      ostringstream tmpIdx;
      tmpIdx << ic;
      TString idx = tmpIdx.str();
      tr.push_back(trees[it]); nm.push_back(names[it]+"_"+idx); isd.push_back(false); cl.push_back(ic+1); ss.push_back(ic+1); 
      csp.push_back(cuts[it]+" && "+c_bb[ic]);
      csf.push_back(cuts[it]+" && !("+c_bb[ic]+")");
    }

    //TCanvas *c_p = c_comp1Dhistos("c_massshape_bbloose_jetpt600to800_"+names[it]+"_p","36.1",nm,tr,isd,"fj_1_sdmass",csp,cl,ss,"m_{SD} [GeV]","a. u.",44,30.,250,false,true);
    //TCanvas *c_f = c_comp1Dhistos("c_massshape_bbloose_jetpt600to800_"+names[it]+"_f","36.1",nm,tr,isd,"fj_1_sdmass",csf,cl,ss,"m_{SD} [GeV]","a. u.",44,30.,250,false,true);
    TCanvas *c_f = c_comp1Dhistos("c_massshape_extrapolation_"+names[it]+"_f","36.1",nm,tr,isd,"fj_1_sdmass",csp,cl,ss,"m_{SD} [GeV]","a. u.",20,50.,250,false,true);
  }
 
}


void compDiscShapes(TString year, TString cat, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cut="0==0", TString suffix="none") {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);

  float intLumi;
  if (year == "2016") { intLumi = 36.8; }
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();

  
  TFile *f_qcd  = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200405/qcd-mg_tree.root" , "READONLY");
  TFile *f_tt   = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200405/top_tree.root" , "READONLY");
  TFile *f_w    = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200405/w-qq_tree.root" , "READONLY");
  TFile *f_z    = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200405/z-qq_tree.root" , "READONLY");
  TFile *f_h    = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200405/ggH-bb_tree.root" , "READONLY");

  TTree *t_qcd = (TTree*)f_qcd->Get("Events");
  TTree *t_tt  = (TTree*)f_tt->Get("Events");
  TTree *t_w   = (TTree*)f_w->Get("Events");
  TTree *t_z   = (TTree*)f_z->Get("Events");
  TTree *t_h   = (TTree*)f_h->Get("Events");


  TString c_match; 
  if (cat == "bb") { c_match = "(fj_1_nbhadrons>=2)"; } 
  if (cat == "cc") { c_match = "(fj_1_nbhadrons==0 && fj_1_nchadrons>=1)"; }
  
  TString c_inc = "(abs(fj_1_eta)<2.4 && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>1000.) && "+cut;
  

  std::vector<TTree*>  trees;  trees.size();
  std::vector<TString> names;  names.clear();
  std::vector<bool>    isdata; isdata.clear();
  std::vector<int>     colors; colors.size();
  std::vector<int>     styles; styles.size();
  std::vector<TString> cuts;   cuts.clear();
  
  //trees.push_back(t_h); names.push_back("higgs");  isdata.push_back(false); colors.push_back(2); styles.push_back(1); cuts.push_back(c_inc+"&&"+c_match+"&& (fj_1_isH<0.8)");
  trees.push_back(t_z); names.push_back("zboson"); isdata.push_back(false); colors.push_back(8); styles.push_back(1); cuts.push_back(c_inc+"&&"+c_match+"&& (fj_1_isZ<0.8)");
  trees.push_back(t_w); names.push_back("wboson"); isdata.push_back(false); colors.push_back(6); styles.push_back(1); cuts.push_back(c_inc+"&& fj_1_nchadrons>=1 && (fj_1_isW<0.8)");
                                                                                             
  TCanvas *c_fj_1_compDisc = c_comp1Dhistos("c_fj_1_compDisc_"+cat,lumi,names,trees,isdata,score,cuts,colors,styles,                                  
					    score,"a. u.",20,0.,1.,false,true);

}



void compDataMC(TString era, TString cat, TString wp, TString score="fj_1_ParticleNetMD_XbbVsQCD", TString cut="0==0", TString suffix="none") {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);
  conf::configuration("tt1l");

  // Processes, paths and lumi for each era
  vector<TString> processes      = {"ttbar-powheg","singletop","w","diboson"};
  vector<TString> process_names  = {"tt","st","wqq","vv"};
  vector<int>     process_color  = {conf::tqq.color, conf::tp2.color, 2, conf::other.color};
  vector<TString> process_legend = {conf::tqq.legend_name, "SingleTop", "W", conf::other.legend_name};

  TString path;
  float intLumi;
  if (era == "2016") {
    //path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/ttbar1l/2016/";
    //path = "/eos/uscms/store/user/pakontax/V2_Training_Official_nanoAODs_NEW_30Jan/2016/";
    //    path = "/uscmst1b_scratch/lpc1/3DayLifetime/loukas/particlenet_ak8_20200914_muon_2016/";
    path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/ttbar1l/20200914/particlenet_ak8_2016_20200914_muon/";
    //path = "/eos/uscms/store/user/pakontax/V2_Training_Official_nanoAODs/2016_Muon_Channel/";
    //    path = "/eos/uscms/store/group/lpcjme/noreplica/NanoHRT/Trees/Apr08/muon/";
    intLumi= 36.8;
  }
  if (era == "2017") {
    //path = "/eos/uscms/store/user/pakontax/V2_Training_Official_nanoAODs/2017_Muon_Channel/";
    path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/ttbar1l/20200914/particlenet_ak8_2017_20200914_muon/";
    intLumi= 44.98;
  }
  if (era == "2018") {
    //path = "/eos/uscms/store/user/pakontax/V2_Training_Official_nanoAODs/2018_Muon_Channel/";
    path = "/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/ttbar1l/20200914/particlenet_ak8_2018_20200914_muon/";
    intLumi= 63.67;
  }
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();


  // Directory to store the datamcplots
  TString name0;
  if (score.Contains("tau21DDT"))      { name0 = "tau21ddt"; }
  if (score.Contains("DeepAK8_"))      { name0 = "dak8"; }
  if (score.Contains("DeepAK8MD"))     { name0 = "dak8md"; }
  if (score.Contains("DeepAK8DDT"))    { name0 = "dak8ddt"; }
  if (score.Contains("ParticleNetMD")) { name0 = "particlenet"; }
  TString dirname1 = name0+"_sf";
  TString nameoutfile = name0+"_datamc_tt1l_"+cat+"_"+wp+"_"+era;

  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }


  // Cuts and matching definition
  TString cut_ = "(n_ak8>=1 && ak8_1_pt>=200 && passmetfilters && passMuTrig)";
  TString c_base     = "(abs(ak8_1_eta)<2.4 && ak8_1_pt>=200. && ak8_1_corr_sdmass>0.1 && ak8_1_corr_sdmass<2500. && ak8_1_mass>0.1) && ("+cut_+")";
  TString c_base_ext = "(abs(ak8_1_eta)<2.4 && ak8_1_pt>=200. && ak8_1_corr_sdmass>0.1 && ak8_1_corr_sdmass<2500. && ak8_1_mass>0.1)";
  TString c_inc     = c_base+" && "+cut_;

  // DDT maps
  TFile *f_dak8ddt_w_0p05 = TFile::Open("./corrections/WvsQCD_ddt_0p05rho.root" , "READONLY");
  TH2D *h_dak8ddt_w_0p05 = (TH2D*)f_dak8ddt_w_0p05->Get("DeepBoosted_WvsQCD_v_rho_v_pT_scaled_yx");
  hrwgt_dak8ddt_w_0p05 = (TH2D*)h_dak8ddt_w_0p05->Clone("hrwgt_dak8ddt_w_0p05");

  TFile *f_dak8ddt_w_0p10 = TFile::Open("./corrections/WvsQCD_ddt_0p10rho.root" , "READONLY");
  TH2D *h_dak8ddt_w_0p10 = (TH2D*)f_dak8ddt_w_0p10->Get("DeepBoosted_WvsQCD_v_rho_v_pT_scaled_yx");
  hrwgt_dak8ddt_w_0p10 = (TH2D*)h_dak8ddt_w_0p10->Clone("hrwgt_dak8ddt_w_0p10");

  TFile *f_dak8ddt_w_0p20 = TFile::Open("./corrections/WvsQCD_ddt_0p20rho.root" , "READONLY");
  TH2D *h_dak8ddt_w_0p20 = (TH2D*)f_dak8ddt_w_0p20->Get("DeepBoosted_WvsQCD_v_rho_v_pT_scaled_yx");
  hrwgt_dak8ddt_w_0p20 = (TH2D*)h_dak8ddt_w_0p20->Clone("hrwgt_dak8ddt_w_0p20");

  TFile *f_dak8ddt_w_0p30 = TFile::Open("./corrections/WvsQCD_ddt_0p30rho.root" , "READONLY");
  TH2D *h_dak8ddt_w_0p30 = (TH2D*)f_dak8ddt_w_0p30->Get("DeepBoosted_WvsQCD_v_rho_v_pT_scaled_yx");
  hrwgt_dak8ddt_w_0p30 = (TH2D*)h_dak8ddt_w_0p30->Clone("hrwgt_dak8ddt_w_0p30");


  TString var = score;
  std::cout << var << "\n";
  if (score.Contains("DeepAK8DDT")) {
    if (wp == "0p05") { var = "(ak8_1_DeepAK8_WvsQCD-rewgtfuncDAK8DDT_w_0p05((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
    if (wp == "0p10") { var = "(ak8_1_DeepAK8_WvsQCD-rewgtfuncDAK8DDT_w_0p10((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
    if (wp == "0p20") { var = "(ak8_1_DeepAK8_WvsQCD-rewgtfuncDAK8DDT_w_0p20((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
    if (wp == "0p30") { var = "(ak8_1_DeepAK8_WvsQCD-rewgtfuncDAK8DDT_w_0p30((2.*log(ak8_1_corr_sdmass/ak8_1_pt)),ak8_1_pt))"; }
  }

  //  var = "ak8_1_DeepAK8MD_WvsQCD";

  std::vector<TTree*> trees; trees.clear();
  std::vector<TString> names; names.clear();
  std::vector<bool>    isdata; isdata.clear();
  std::vector<int> colors; colors.clear();
  std::vector<int> styles; styles.clear();
  std::vector<TString> cuts; cuts.clear();
  std::vector<TString> legends; legends.clear();

  TFile *f_data  = TFile::Open(path+"/data/singlemu_tree.root" , "READONLY");
  TTree *t_data  = (TTree*)f_data->Get("Events");
  trees.push_back(t_data); names.push_back("data"); isdata.push_back(true); colors.push_back(1); styles.push_back(1); cuts.push_back(c_inc); legends.push_back("Data");

  for (int i0=0; i0<processes.size(); ++i0) {
    //TFile *f_ = TFile::Open(path+"/REPO_mc/mc_nom/"+processes[i0]+"_tree.root" , "READONLY");
    TFile *f_ = TFile::Open(path+"/mc_nom/"+processes[i0]+"_tree.root" , "READONLY");
    TTree *t_ = (TTree*)f_->Get("Events");
    trees.push_back(t_); names.push_back(process_names[i0]); isdata.push_back(false); colors.push_back(process_color[i0]); styles.push_back(1); cuts.push_back(c_inc); legends.push_back(process_legend[i0]);
  }

  //  TString cutB = "leptonicW_pt>175. && ak8_1_mass>50. && ak8_1_mass<250. && ("+var+">0.) && ak8_1_pt>400. && ak8_1_pt<800."; 
  TString cutB = "leptonicW_pt>175. && ak8_1_mass>50. && ak8_1_mass<250."; 
  makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,cutB,var,20,-1.,1.,"DeepAK8-DDT "+wp+"_b_"+era,legends,colors,false,0);
  //  makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,"0==0","leptonicW_pt",20,0.,500.,"leptonicW_pt_"+era,legends,colors,false,0);
  //  makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,cutB,"muon_pt",50,0.,500.,"muon_pt_b_"+era,legends,colors,false,0);
  //  makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,cutB,"ak8_1_pt",20,0.,1000.,"ak8_1_pt_b_"+era,legends,colors,false,0);
  //  makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,cutB,"ak8_1_mass",50,0.,250.,"ak8_1_mass_b_"+era,legends,colors,false,0);
  //  makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,cutB,"ak8_1_corr_sdmass",50,0.,250.,"ak8_1_corr_sdmass_b_"+era,legends,colors,false,0);
  //  makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,"0==0","ak8_1_eta",12,-3.,3.,"ak8_1_eta_"+era,legends,colors,false,0);
  // makePlotDataMC(nameoutfile,dirname1,trees,names,lumi,c_inc,"0==0",var,20,-1.,1.,"DeepAK8-DDT "+wp+"_"+era,legends,colors,false,0);
  
  // makePlotDataMC("particlenet_ht","testplotdir",trees,names,lumi,cuts,"0==0","ht",30,1000.,2500.,"H_{T} [GeV]",legends,colors,false,true);
  //makePlotDataMC(,"testplotdir",trees,names,lumi,cuts,"fj_1_ParticleNetMD_XccVsQCD>=0.76","fj_1_sdmass",30,50.,200.,"m_{SD} [GeV]",legends,colors,false,true);
  // makePlotDataMC("particlenet_xbbvsqcdt_pass","testplotdir",trees,names,lumi,cuts,"fj_1_ParticleNetMD_XccVsQCD>=0.94","fj_1_sdmass",30,50.,200.,"m_{SD} [GeV]",legends,colors,false,true);
  //makePlotDataMC("particlenet_xbbvsqcd_fail","testplotdir",trees,names,lumi,cuts,"fj_1_ParticleNetMD_XbbVsQCD<0.93","fj_1_sdmass",30,50.,200.,"m_{SD} [GeV]",legends,colors,false,true);
  /*
  TH1F *h_passovfail_qcd  = createRatioPassOvFail(trees[1],names[1],lumi,cuts[1],"fj_1_ParticleNetMD_XbbVsQCD>=0.93","fj_1_sdmass",30,50.,200.,"m_{SD} [GeV]",colors[1],styles[1],false,true,isdata[1]);
  TH1F *h_passovfail_data = createRatioPassOvFail(trees[0],names[0],lumi,cuts[0],"fj_1_ParticleNetMD_XbbVsQCD>=0.93","fj_1_sdmass",30,50.,200.,"m_{SD} [GeV]",colors[0],styles[0],false,true,isdata[0]);
  TCanvas *c_passovfail = new TCanvas("c_passovfail","c_passovfail",500,500);
  h_passovfail_qcd->GetYaxis()->SetTitle("N_{pass}/N_{fail}");
  h_passovfail_qcd->GetYaxis()->SetRangeUser(0.4,1.6);
  h_passovfail_qcd->GetXaxis()->SetTitle("m_{SD} [GeV]");
  h_passovfail_qcd->Draw("L E0");
  h_passovfail_data->Draw("P E0 sames");
  */
}

void getSF(TString dir, TString algo, TString cat, TString wp, TString era) {

  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);

  std::vector<TString> ptranges; ptranges.clear();
  ptranges.push_back("lowpt");
  ptranges.push_back("lowmedpt");
  ptranges.push_back("medpt");
  ptranges.push_back("medhipt");
  ptranges.push_back("hipt");
  //particlenet_bb_t_2016_nom/fitDiagnostics_particlenet_bb_t_2016_lowpt.root
  TGraphAsymmErrors *sfs = getSFGraph(dir,algo,cat,wp,era,ptranges,1,2,1);

}



TGraphAsymmErrors *getSFGraph(TString dir, TString algo, TString cat, TString wp, TString era, std::vector<TString> ptrange, int ialgo, int nalgos, int color) {
  
  double xval[20], xvalerr[20];
  double yval[20], yvalerrlo[20], yvalerrhi[20];
  
  for (unsigned int iptrange=0; iptrange<ptrange.size(); ++iptrange) {
    TFile *fdiag = TFile::Open("./"+dir+"/fitDiagnostics_"+algo+"_"+cat+"_"+wp+"_"+era+"_"+ptrange[iptrange]+".root", "READONLY" );
    TTree *t_   = (TTree*)fdiag->Get("tree_fit_sb");
    
    TString strCat, strCatLoErr, strCatHiErr;   
    strCat = "SF_z"; strCatLoErr = "SF_zLoErr"; strCatHiErr = "SF_zHiErr";
    TH1F *h_sf_cat      = new TH1F("h_sf_cat"      , "h_sf_cat"      , 10000,0.,10.);
    TH1F *h_sf_catLoErr = new TH1F("h_sf_catLoErr" , "h_sf_catLoErr" , 10000,0.,10.);
    TH1F *h_sf_catHiErr = new TH1F("h_sf_catHiErr" , "h_sf_catHiErr" , 10000,0.,10.);
    t_->Project("h_sf_cat",strCat);
    t_->Project("h_sf_catLoErr",strCatLoErr);
    t_->Project("h_sf_catHiErr",strCatHiErr);
    
    xval[iptrange]    = (iptrange+0.7) + 0.07*(ialgo - 1.*(nalgos-2));
    xvalerr[iptrange] = 0.03;
    yval[iptrange]    = h_sf_cat->GetMean();
    yvalerrlo[iptrange] = h_sf_catLoErr->GetMean();
    yvalerrhi[iptrange] = h_sf_catHiErr->GetMean();
    //   if (algos[ialgo]=="sdtau21" && ptrange[iptrange]=="low") { yvalerrlo[iptrange] = 0.01; yvalerrhi[iptrange] = 0.01; }
    
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    std::cout << h_sf_cat->GetMean() << "$^{+" << h_sf_catHiErr->GetMean() << "}_{-" << h_sf_catLoErr->GetMean() << "}$ & ";   
    
    h_sf_cat->Delete();
    h_sf_catLoErr->Delete();
    h_sf_catHiErr->Delete();
    
  } // end looping over the pt bins
  std::cout << "\n";
  
  TGraphAsymmErrors *gr_sf_tmp = new TGraphAsymmErrors(ptrange.size(),xval,yval,xvalerr,xvalerr,yvalerrlo,yvalerrhi); 
  gr_sf_tmp->SetLineColor(color);
  gr_sf_tmp->SetMarkerColor(color);
  return gr_sf_tmp;
}




void getSampleComposition(TString path2file, TString filename, TString score, TString year, TString ptrange, TString category) {
  TH1::SetDefaultSumw2(kTRUE);

  const int dir_err = system("mkdir -p ./"+(TString)path2file+"/plots_datamc");
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }
  //particlenet_bb_t_2016/particlenet_bb_t_2016_lowpt_templates_p.root
  TFile *fdiag = TFile::Open("./"+path2file+"/"+filename+"_"+score+"_"+year+"_"+ptrange+"_templates_"+category+".root", "READONLY" );

  // get prefit histograms
  TString prefitstr = "shapes_prefit";//h_particlenet_bb_t_2016_h_p
  TH1F *h_prefit_h_     = (TH1F*)fdiag->Get("h");     h_prefit_h_->SetName("h_h_");
  TH1F *h_prefit_zbb_   = (TH1F*)fdiag->Get("zbb");   h_prefit_zbb_->SetName("h_zbb_");
  TH1F *h_prefit_zcc_   = (TH1F*)fdiag->Get("zcc");   h_prefit_zcc_->SetName("h_zcc_");
  TH1F *h_prefit_zqq_   = (TH1F*)fdiag->Get("zqq");   h_prefit_zqq_->SetName("h_zqq_");
  TH1F *h_prefit_w_     = (TH1F*)fdiag->Get("w");     h_prefit_w_->SetName("h_w_");
  TH1F *h_prefit_tt_    = (TH1F*)fdiag->Get("tt");    h_prefit_tt_->SetName("h_tt_");
  TH1F *h_prefit_qcd_   = (TH1F*)fdiag->Get("qcd");   h_prefit_qcd_->SetName("h_qcd_");

  // get data
  TH1F *h_data_ = (TH1F*)fdiag->Get("data_obs");  h_data_->SetName("h_data_");

  std::cout << std::fixed;
  std::cout << std::setprecision(2);


  float totSM = 
    h_prefit_w_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))   + 
    h_prefit_zbb_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) + 
    h_prefit_zcc_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) + 
    h_prefit_zqq_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) +
    h_prefit_h_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))   +
    h_prefit_tt_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))  +
    h_prefit_qcd_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) ;
  
  // yields
  std::cout
    << h_prefit_w_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))   << " & "
    << h_prefit_zbb_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) << " & "
    << h_prefit_zcc_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) << " & "
    << h_prefit_zqq_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) << " & "
    << h_prefit_h_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))   << " & "
    << h_prefit_tt_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))  << " & "
    << h_prefit_qcd_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.)) << " & "
    << totSM                     << " & "
    << h_data_->Integral(h_data_->FindBin(85.),h_data_->FindBin(110.))        << " \\\\\n ";


  // purity
  std::cout 
    << h_prefit_w_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))/totSM << " & "
    << h_prefit_zbb_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))/totSM << " & "  
    << h_prefit_zcc_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))/totSM << " & "
    << h_prefit_zqq_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))/totSM << " & "
    << h_prefit_h_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))/totSM << " & "
    << h_prefit_tt_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))/totSM << " & "
    << h_prefit_qcd_->Integral(h_prefit_w_->FindBin(85.),h_prefit_w_->FindBin(110.))/totSM << " \\\\\n ";
  
}


void makeDataMCPlotFromCombine(TString path2file, TString era, TString category, TString category_extra, TString wp, TString name, TString passOrFail, 
			       float xmin, float xmax, int nbins,TString xaxisname, bool log, TString sample ="") {

  std::cout << " In makeDataMCPlotFromCombine\n";

  setTDRStyle();
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
  TH1::SetDefaultSumw2(kTRUE);
  conf::configuration(sample);

  if (xaxisname == "mass") { xaxisname = "m_{SD} [GeV]"; }

  const int dir_err = system("mkdir -p ./"+(TString)path2file+"_sf/plots_postfit");
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); } 

  TString fdiag_ = "./"+path2file+"_sf/fitdir/fitdiagnostics_"+path2file+"_"+sample+"_"+category+"_"+wp+"_"+era+"_"+name+".root";
  std::cout << " Opening fitdiagnostics: " << fdiag_ << "\n";
  TFile *fdiag = TFile::Open(fdiag_, "READONLY" );

  std::vector<TString> processes; processes.clear();
  std::vector<int> colors;        colors.clear();
  std::vector<TString> legends;   legends.clear();
  if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") ) {
    if (category == "w") { 
      std::cout << " sample = " << sample << " , category = " << category << "\n";
      processes.push_back("tp3"); processes.push_back("wqq"); processes.push_back("tp1"); processes.push_back("other"); processes.push_back("total");
      colors.push_back(conf::tp3.color); colors.push_back(conf::tp2.color); colors.push_back(conf::tp1.color); colors.push_back(conf::other.color); colors.push_back(602);
      legends.push_back(conf::tp3.legend_name); legends.push_back(conf::tp2.legend_name); legends.push_back(conf::tp1.legend_name); legends.push_back(conf::other.legend_name); legends.push_back("Total SM");
    }
  }

  if ( (category == "bb") && (category_extra == "") ){
    std::cout << " sample = " << sample << " , category = " << category << "\n";
    processes.push_back("zqq"); processes.push_back("wqq"); processes.push_back("tqq"); processes.push_back("qcd"); processes.push_back("total");
    colors.push_back(conf::zqq.color); colors.push_back(conf::wqq.color); colors.push_back(conf::tqq.color); colors.push_back(conf::qcd.color); colors.push_back(602);
    legends.push_back(conf::zqq.legend_name); legends.push_back(conf::wqq.legend_name); legends.push_back(conf::tqq.legend_name); legends.push_back(conf::qcd.legend_name); \
    legends.push_back("Total SM");
  }

  if ( (category == "bb") && (category_extra == "jmcor") ) {
    std::cout << " sample = " << sample << " , category = " << category << "\n";
    processes.push_back("tp1"); processes.push_back("wqq"); processes.push_back("tqq"); processes.push_back("other"); processes.push_back("total");
    colors.push_back(conf::tp1.color); colors.push_back(conf::wqq.color); colors.push_back(conf::tqq.color); colors.push_back(conf::other.color); colors.push_back(602);
    legends.push_back(conf::tp1.legend_name); legends.push_back(conf::wqq.legend_name); legends.push_back(conf::tqq.legend_name); legends.push_back(conf::other.legend_name); \
    legends.push_back("Total SM");
  }

  TString extra_ = "/"; if (sample == "comb") { extra_ = "/ch2_"; } 

  // get prefit histograms
  TString prefitstr = "shapes_prefit";

  TTree *tree = (TTree*)fdiag->Get("tree_fit_sb");
  TH1F *h_fit_status = new TH1F("h_fit_status_","h_fit_status_",20,-10.,10.); 
  tree->Project("h_fit_status_","fit_status");
  TString postfitstr = "shapes_fit_s"; if (h_fit_status->GetMean()<0.) { postfitstr = "shapes_prefit"; }

  std::vector<TH1F*> h_prefit; h_prefit.clear();
  std::vector<TH1F*> h_postfit; h_postfit.clear();
  for (unsigned int i0=0; i0<processes.size(); ++i0) {
    std::cout << prefitstr+extra_+passOrFail+"/"+processes[i0] << "\n";
    // get prefit histos
    TH1F *h_prefit_ = (TH1F*)fdiag->Get(prefitstr+extra_+passOrFail+"/"+processes[i0]); h_prefit_->SetName("h_pre_"+name+"_"+processes[i0]+"_"+passOrFail+"_");
    TH1F *h_prefit__ = rescaleXaxis(h_prefit_,conf::minX,conf::maxX); h_prefit_->SetName("h_pre_"+name+"_"+processes[i0]+"_"+passOrFail);
    h_prefit__->SetLineColor(colors[i0]); h_prefit__->SetLineStyle(2); h_prefit__->SetLineWidth(3); h_prefit__->SetMarkerSize(0);
    h_prefit.push_back(h_prefit__);

    // get postfit histograms
    TH1F *h_postfit_ = (TH1F*)fdiag->Get(postfitstr+extra_+passOrFail+"/"+processes[i0]); h_postfit_->SetName("h_post_"+name+"_"+processes[i0]+"_"+passOrFail+"_");
    TH1F *h_postfit__ = rescaleXaxis(h_postfit_,conf::minX,conf::maxX); h_postfit_->SetName("h_post_"+name+"_"+processes[i0]+"_"+passOrFail);
    h_postfit__->SetLineColor(colors[i0]); h_postfit__->SetLineStyle(1); h_postfit__->SetLineWidth(3); h_postfit__->SetMarkerSize(0);
    h_postfit.push_back(h_postfit__);
   
  }

  // get data
  TGraphAsymmErrors *h_data = (TGraphAsymmErrors*)fdiag->Get("shapes_prefit"+extra_+passOrFail+"/data"); h_data->SetName("h_data_"+name+"_"+passOrFail);
  rescaleXaxis(h_data, conf::minX, (conf::maxX-conf::minX)/(float)conf::binsX);
  h_data->SetLineColor(1); h_data->SetLineWidth(3); 
  h_data->SetMarkerColor(1); h_data->SetMarkerStyle(20); h_data->SetMarkerSize(1.2); 

  int numOfMC = h_prefit.size();
  TH1F *h_r_prefit = getDataMCratio(h_data,h_prefit[numOfMC-1]); h_r_prefit->SetName("h_r_prefit_"+name+"_"+passOrFail);
  h_r_prefit->SetMarkerSize(0); h_r_prefit->SetLineColor(42); h_r_prefit->SetLineStyle(2);
  
  TH1F *h_r_postfit = getDataMCratio(h_data,h_postfit[numOfMC-1]); h_r_postfit->SetName("h_r_postfit_"+name+"_"+passOrFail);
  h_r_postfit->SetMarkerColor(1); h_r_postfit->SetLineColor(1);

  TLegend* leg = new TLegend(0.60,0.55,0.94,0.86);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetLineWidth(0);
  //leg->AddEntry(h_postfit_h,"Higgs","L");
  for (unsigned int i0=0; i0<h_prefit.size(); ++i0) { leg->AddEntry(h_postfit[i0],legends[i0],"L"); }

  TLegend* legfit = new TLegend(0.4,0.7,0.65,0.86);
  legfit->SetFillStyle(0);
  legfit->SetFillColor(0);
  legfit->SetLineWidth(0);
  legfit->AddEntry(h_prefit[numOfMC-1],"Pre-fit","L");
  legfit->AddEntry(h_postfit[numOfMC-1],"Post-fit","L");

  TLegend* legfitr = new TLegend(0.35,0.80,0.55,0.93);
  legfitr->SetFillStyle(0);
  legfitr->SetFillColor(0);
  legfitr->SetLineWidth(0);
  legfitr->AddEntry(h_r_prefit,"Pre-fit","L");
  legfitr->AddEntry(h_r_postfit,"Post-fit","L");

  TPaveText *pt_cms = new TPaveText(0.11,0.77,0.4,0.9,"NDC");
  pt_cms->SetFillStyle(0);
  pt_cms->SetFillColor(0);
  pt_cms->SetLineWidth(0);
  pt_cms->AddText("CMS");
  pt_cms->SetTextSize(0.08);

  TPaveText *pt_preliminary = new TPaveText(0.2,0.63,0.4,0.9,"NDC");
  pt_preliminary->SetFillStyle(0);
  pt_preliminary->SetFillColor(0);
  pt_preliminary->SetLineWidth(0);
  pt_preliminary->AddText("Preliminary");
  pt_preliminary->SetTextFont(52);
  pt_preliminary->SetTextSize(0.06);

  TLatex pt_lumi;
  const char *longstring;
  if (path2file.Contains("2016")) { longstring = "35.9 fb^{-1} (13 TeV)"; }
  if (path2file.Contains("2017")) { longstring = "41.53 fb^{-1} (13 TeV)"; }
  if (path2file.Contains("2018")) { longstring = "59.74 fb^{-1} (13 TeV)"; }
  if (path2file.Contains("test")) { longstring = "35.9 fb^{-1} (13 TeV)"; }
  pt_lumi.SetTextSize(0.07);
  pt_lumi.SetTextFont(42);

  TString logstr = "lin"; if (log) { logstr = "log"; } 
  TCanvas *c = new TCanvas("c_"+name+"_"+category,"c_"+name+"_"+category,600,600); 
  c->SetName("c_"+name+"_"+category);

  TPad *pMain  = new TPad("pMain_"+name+"_"+category,"pMain+"+name+"_"+category,0.0,0.35,1.0,1.0);
  pMain->SetRightMargin(0.05);
  pMain->SetLeftMargin(0.17);
  pMain->SetBottomMargin(0.02);
  pMain->SetTopMargin(0.1);
  
  TPad *pRatio = new TPad("pRatio_"+name+"_"+category,"pRatio_"+name+"_"+category,0.0,0.0,1.0,0.35);
  pRatio->SetRightMargin(0.05);
  pRatio->SetLeftMargin(0.17);
  pRatio->SetTopMargin(0.02);
  pRatio->SetBottomMargin(0.4);
  pMain->Draw();
  pRatio->Draw();
 
  pMain->cd();
  float maxyld = h_postfit[numOfMC-1]->GetMaximum(); 
  if (h_prefit[numOfMC-1]->GetMaximum()>h_postfit[numOfMC-1]->GetMaximum()) { maxyld = h_prefit[numOfMC-1]->GetMaximum(); }
  if (log) { gPad->SetLogy(); h_prefit[numOfMC-1]->GetYaxis()->SetRangeUser(0.1,10.*maxyld); } else { h_prefit[numOfMC-1]->GetYaxis()->SetRangeUser(0.,1.8*maxyld); }
  h_prefit[numOfMC-1]->GetXaxis()->SetLabelSize(0.);
  h_prefit[numOfMC-1]->GetYaxis()->SetTitle("Events / bin");
  h_prefit[numOfMC-1]->GetXaxis()->SetTitle(xaxisname);
  h_prefit[numOfMC-1]->Draw("HIST E0");
  for (unsigned int i0=0; i0<h_prefit.size(); ++i0)  { h_prefit[i0]->Draw("HIST E0 sames"); }
  for (unsigned int i0=0; i0<h_postfit.size(); ++i0) { h_postfit[i0]->Draw("HIST E0 sames"); }
  h_data->Draw("P sames");
  leg->Draw("sames");
  pt_cms->Draw("sames");
  pt_preliminary->Draw("sames");
  std::cout << " i m here 12-d\n";
  pt_lumi.DrawLatexNDC(0.64,0.93,longstring);
  std::cout << " i m here 12-e\n";
  c->RedrawAxis();
  pRatio->cd();
  gPad->SetBottomMargin(0.2);
  h_r_postfit->GetYaxis()->SetTitleOffset(0.75);
  h_r_postfit->GetYaxis()->SetTitleSize(0.1);
  h_r_postfit->GetYaxis()->SetLabelSize(0.08);
  h_r_postfit->GetXaxis()->SetTitleSize(0.1);
  h_r_postfit->GetXaxis()->SetLabelSize(0.08);
  h_r_postfit->GetXaxis()->SetTitle(xaxisname);
  h_r_postfit->GetYaxis()->SetTitle("Data / Post-fit");
  h_r_postfit->GetYaxis()->SetRangeUser(0.8,1.2);
  h_r_postfit->Draw("P E0");
  h_r_prefit->Draw("HIST sames");
  c->RedrawAxis();
 
  if (log) {
    c->Print(path2file+"_sf/plots_datamc/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wp+"_"+passOrFail+"_log.pdf");
    c->Print(path2file+"_sf/plots_datamc/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wp+"_"+passOrFail+"_log.png");
  } 
  else {
    c->Print(path2file+"_sf/plots_postfit/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wp+"_"+passOrFail+"_lin.pdf");
    c->Print(path2file+"_sf/plots_postfit/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wp+"_"+passOrFail+"_lin.png");
  }

} // end of makeDataMCPlotsFromCombine



TH1F *createRatioPassOvFail(TTree *tree, TString name, TString lumi, TString cut, TString addCut, TString var, int bins, float xmin, float xmax, TString xname, 
			    int color, int style, bool log, bool norm, bool isdata) {
  
  TH1F *h_pass = create1Dhisto(name,tree,lumi,cut+"&&"+addCut,var,bins,xmin,xmax,log,color,style,"h_pass_"+name,norm,isdata);
  TH1F *h_fail = create1Dhisto(name,tree,lumi,cut+"&& !("+addCut+")",var,bins,xmin,xmax,log,color,style,"h_fail_"+name,norm,isdata);
  TH1F *h_pass_ov_fail = (TH1F*)h_pass->Clone("h_pass_ov_fail_"+name);
  h_pass_ov_fail->Divide(h_pass,h_fail);
  return h_pass_ov_fail;
}


TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, std::vector<TString> vars, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm) {

  std::vector<TH1F*> histos; histos.clear(); 
  for (int itree=0; itree<trees.size(); ++itree) {
    TH1F *h_tmp = create1Dhisto(names[itree],trees[itree],lumi,cuts[itree],vars[itree],nbins,xmin,xmax,uselog,colors[itree],styles[itree],"h_"+names[itree]+"_"+c_name,norm,isdata[itree]); 
    h_tmp->SetFillColor(0);
    histos.push_back(h_tmp);
  }

  TCanvas *c_comp1Dhistos = new TCanvas("c_comp1Dhistos_"+c_name,"c_comp1Dhistos_"+c_name,500,500);
  for (int ihisto=0; ihisto<histos.size(); ++ihisto) {
    if (ihisto==0) {
      histos[ihisto]->GetXaxis()->SetTitle(xname);
      histos[ihisto]->GetYaxis()->SetTitle(yname);
      histos[ihisto]->GetYaxis()->SetRangeUser(0.,1.);
      histos[ihisto]->Draw("HIST E0");
    }
    else {
      histos[ihisto]->Draw("HIST E0 sames");
    }
  }

  c_comp1Dhistos->Print("./plots/"+c_name+".pdf");
  return c_comp1Dhistos;
}


TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm) {

  std::vector<TH1F*> histos; histos.clear(); 
  for (int itree=0; itree<trees.size(); ++itree) {
    TH1F *h_tmp = create1Dhisto(names[itree],trees[itree],lumi,cuts[itree],var,nbins,xmin,xmax,uselog,colors[itree],styles[itree],"h_"+names[itree]+"_"+c_name,norm,isdata[itree]); 
    h_tmp->SetFillColor(0);
    histos.push_back(h_tmp);
  }
  
  // get histo max value
  float maxYval = 0.;
  for (int ihisto=0; ihisto<histos.size(); ++ihisto) { 
    float maxYval_tmp = histos[ihisto]->GetBinContent(histos[ihisto]->GetMaximumBin());
    if (maxYval_tmp>maxYval) { maxYval = maxYval_tmp; }
  }

  TCanvas *c_comp1Dhistos = new TCanvas("c_comp1Dhistos_"+c_name,"c_comp1Dhistos_"+c_name,500,500);
  for (int ihisto=0; ihisto<histos.size(); ++ihisto) {
    if (ihisto==0) {
      histos[ihisto]->GetXaxis()->SetTitle(xname);
      histos[ihisto]->GetYaxis()->SetTitle(yname);
      histos[ihisto]->GetYaxis()->SetRangeUser(0.,1.3*maxYval);
      histos[ihisto]->Draw("HIST E0");
    }
    else {
      histos[ihisto]->Draw("HIST E0 sames");
    }
  }
  //gSystem->cd("/uscms_data/d3/loukas/ParticleNetInData/CMSSW_10_2_18/src/PhysicsTools/NanoHRTTools/plotting/plots/particlenet/");
  c_comp1Dhistos->Print("./plots/"+c_name+".pdf");
  return c_comp1Dhistos;
}

TCanvas *c_ratioOf1Dhistos(TString c_name, TString lumi, std::vector<TString> namest, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, 
			   std::vector<TString> namesc, std::vector<TString> cuts,
			   TString xname, TString yname, int nbins, float xmin, float xmax, bool uselog, int norm) {

  std::vector<unsigned int> colors;
  colors.push_back(1);
  colors.push_back(2);
  colors.push_back(4);

  std::vector<unsigned int> lineStyles;
  lineStyles.push_back(1);
  lineStyles.push_back(2);
  lineStyles.push_back(3);

  std::vector<std::vector<TH1F*>> histos; histos.clear();
  for (int itree=0; itree<trees.size(); ++itree) {

    std::vector<TH1F*> histos_tmp; histos_tmp.clear();
    for (int icut=0; icut<cuts.size(); ++icut) {
      TH1F *h_tmp = create1Dhisto(namest[itree],trees[itree],lumi,cuts[icut],var,nbins,xmin,xmax,uselog,colors[itree],lineStyles[itree],"h_"+namest[itree]+"_"+namesc[icut]+"_"+c_name,norm,isdata[itree]);
      h_tmp->SetFillColor(0);
      histos_tmp.push_back(h_tmp);
    }
    histos.push_back(histos_tmp);
  }

  // calculate ratio wrt "0" histogram
  std::cout << " histos = " << histos.size() << "\n";
  std::vector<std::vector<TH1F*>> ratios; ratios.clear();
  for (int ihistos=0; ihistos<histos.size(); ++ihistos) {
    
    std::vector<TH1F*> ratios_tmp; ratios_tmp.clear();
    for (int ihisto=0; ihisto<histos[ihistos].size(); ++ihistos) {
      if (ihisto==0) { continue; }
      TString name = "h_r_"+(TString)histos[ihistos][ihisto]->GetName();
      TH1F* h_r = (TH1F*)histos[ihistos][ihisto]->Clone(name); 
      h_r->Divide(histos[ihistos][ihisto],histos[ihistos][0]);
      h_r->SetLineColor(colors[ihisto-1]); h_r->SetLineStyle(lineStyles[ihisto-1]);
      ratios_tmp.push_back(h_r);
    }
    ratios.push_back(ratios_tmp);
  }
  std::cout << ratios.size() << "\n";
  // plot them
  TCanvas *c_r = new TCanvas("c_ratioOf1Dhistos_"+c_name,"c_ratioOf1Dhistos_"+c_name,500,500);
  for (int ihistos=0; ihistos<ratios.size(); ++ihistos) {
    for (int ihisto=0; ihisto<ratios[ihistos].size(); ++ihisto) {
      
      if ( (ihistos==0) && (ihisto==0) ) {
	ratios[ihistos][ihisto]->GetXaxis()->SetTitle(xname);
	ratios[ihistos][ihisto]->GetYaxis()->SetTitle(yname);
	//ratios[ihistos][ihisto]->GetYaxis()->SetRangeUser(0.,2.);
	ratios[ihistos][ihisto]->Draw("HIST E0");
      }
      else { ratios[ihistos][ihisto]->Draw("HIST E0 sames"); }
    }
  }

  c_r->Print("./plots/"+c_name+".pdf");
  return c_r;
}


void makePlotDataMC(TString name, TString dir,
                    std::vector<TTree*> trees, std::vector<TString> names, TString lumi,TString cut, TString addcut,
                    TString var,int nbins,float xmin,float xmax,TString xaxisname,
                    std::vector<TString> legends, std::vector<int> colors, bool logy, int norm) {

  TH1::SetDefaultSumw2(kTRUE);
  std::cout << " working on " << name << " ...\n";

  TString sample;
  if (name.Contains("qcd"))     { sample = "qcd"; }
  if (name.Contains("photons")) { sample = "photons"; }


  TH1F *h_sm; TH1F *h_data; std::vector<TH1F*> mc_histos; mc_histos.clear();
  THStack *hs = new THStack("hs_"+name,"hs_"+name);

  TLegend* leg = new TLegend(0.65,0.55,0.92,0.88);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetLineWidth(0);

  TString logstr = "lin"; if (logy) { logstr = "log"; }

  TCanvas *c_ = new TCanvas(name,name,600,600); c_->SetName("c_"+name);
  TPad *pMain  = new TPad("pMain_"+name,"pMain+"+name,0.0,0.25,1.0,1.0);
  pMain->SetRightMargin(0.05);
  pMain->SetLeftMargin(0.17);
  pMain->SetBottomMargin(0.);
  pMain->SetTopMargin(0.1);

  TPad *pRatio = new TPad("pRatio_"+name,"pRatio_"+name,0.0,0.0,1.0,0.25);
  pRatio->SetRightMargin(0.05);
  pRatio->SetLeftMargin(0.17);
  pRatio->SetTopMargin(0.);
  pRatio->SetBottomMargin(0.37);
  pMain->Draw();
  pRatio->Draw();


  for (unsigned int i=0; i<trees.size(); ++i) {

    ostringstream tmpCount;
    tmpCount << i;
    TString count = tmpCount.str();

    bool data = false;  if (i==0) { data = true; }

    pMain->cd();
    TH1F *h = create1Dhisto(names[i],trees[i],lumi,cut+" && "+addcut,var,nbins,xmin,xmax,false,1,1,"h_"+name+"_"+count,false,data); h->SetName(names[i]);
    if (i==0) {
      h->SetMarkerSize(1.2); h->SetMarkerStyle(20); h->SetLineWidth(1); h->SetLineColor(colors[i]); h->SetFillColor(0);
      leg->AddEntry(h,legends[i],"P");
      h->GetYaxis()->SetTitleOffset(1.45);
      h->GetXaxis()->SetTitle(xaxisname);
      if (norm == 0) { h->GetYaxis()->SetTitle("Events / bin"); } else { h->GetYaxis()->SetTitle("a. u."); }
      if (logy) {
        gPad->SetLogy();
        h->GetYaxis()->SetRangeUser(1.,1000.*h->GetBinContent(h->GetMaximumBin()));
        if (var.Contains("tau") || var.Contains("eta") || var.Contains("ecf")) { h->GetYaxis()->SetRangeUser(1.,100000.*h->GetBinContent(h->GetMaximumBin())); }
      }
      else { h->GetYaxis()->SetRangeUser(0.,1.3*h->GetBinContent(h->GetMaximumBin())); }
      h->Draw("P E0");
      h_data = (TH1F*)h->Clone("h_"+name+"_data");
    }
    else {
      h->SetFillColor(colors[i]); h->SetMarkerSize(0); h->SetLineColor(colors[i]);
      leg->AddEntry(h,legends[i],"FL");
      h->Draw("HIST E0 sames");
    }

    if (i==1) { h_sm = (TH1F*)h->Clone("h_sm"); hs->Add(h); }
    if (i>1 ) { h_sm->Add(h); hs->Add(h); }

  } // end of looping over the trees

  // norm MC
  float normVal = h_data->Integral()/h_sm->Integral();
  std::cout << h_data->Integral() << " " << h_sm->Integral() << "\n";
  if (norm !=0) { h_sm->Scale(normVal); }

  // data mc ratio
  TH1F *h_r = (TH1F*)h_data->Clone("h_"+name+"_r"); h_r->Divide(h_data,h_sm);

  TList *histos = hs->GetHists();
  TIter next(histos);
  TH1F *hist;
  while ((hist =(TH1F*)next())) { hist->Scale(normVal); }

  // continue drawing
  pMain->cd();
  hs->Draw("HIST E0 sames");
  h_data->Draw("P E0 sames");
  leg->Draw("sames");
  c_->RedrawAxis();
  pMain->RedrawAxis();

  pRatio->cd();
  gPad->SetBottomMargin(0.2);
  h_r->GetYaxis()->SetTitleOffset(0.9);
  h_r->GetYaxis()->SetTitleSize(0.1);
  h_r->GetYaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitleSize(0.1);
  h_r->GetXaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitle(xaxisname);
  h_r->GetYaxis()->SetTitle("Data / MC");
  h_r->GetYaxis()->SetRangeUser(0.5,1.5);
  h_r->Draw("P E0");


  const int dir_err = system("mkdir -p ./"+dir);
  if (-1 == dir_err) {
    printf("Error creating directory!n");
    exit(1);
  }
  c_->Print(dir+"/"+name+"_"+logstr+".png");
  c_->Print(dir+"/"+name+"_"+logstr+".pdf");

}


void makePlotDataMC(TString name, TString dir,
                    std::vector<TTree*> trees, std::vector<TString> names, TString lumi,std::vector<TString> cuts, TString addcut,
                    TString var,int nbins,float xmin,float xmax,TString xaxisname,
                    std::vector<TString> legends, std::vector<int> colors, bool logy, int norm) {

  TH1::SetDefaultSumw2(kTRUE);
  std::cout << " working on b " << name << " ...\n";

  TString sample;
  if (name.Contains("qcd"))     { sample = "qcd"; }
  if (name.Contains("photons")) { sample = "photons"; }

  TH1F *h_sm; TH1F *h_data; std::vector<TH1F*> mc_histos; mc_histos.clear();
  THStack *hs = new THStack("hs_"+name,"hs_"+name);

  TLegend* leg = new TLegend(0.65,0.55,0.92,0.88);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetLineWidth(0);

  TString logstr = "lin"; if (logy) { logstr = "log"; }

  TCanvas *c_ = new TCanvas(name,name,600,600); c_->SetName("c_"+name);
  TPad *pMain  = new TPad("pMain_"+name,"pMain+"+name,0.0,0.25,1.0,1.0);
  pMain->SetRightMargin(0.05);
  pMain->SetLeftMargin(0.17);
  pMain->SetBottomMargin(0.);
  pMain->SetTopMargin(0.1);

  TPad *pRatio = new TPad("pRatio_"+name,"pRatio_"+name,0.0,0.0,1.0,0.25);
  pRatio->SetRightMargin(0.05);
  pRatio->SetLeftMargin(0.17);
  pRatio->SetTopMargin(0.);
  pRatio->SetBottomMargin(0.37);
  pMain->Draw();
  pRatio->Draw();


  std::vector<TH1F*> histoCol; histoCol.clear();
  for (unsigned int i=0; i<trees.size(); ++i) {

    ostringstream tmpCount;
    tmpCount << i;
    TString count = tmpCount.str();

    bool data = false;  if (i==0) { data = true; }

    pMain->cd();
    sample = names[i];
    TH1F *h = create1Dhisto(sample,trees[i],lumi,cuts[i]+" && "+addcut,var,nbins,xmin,xmax,false,1,1,"h_"+name+"_"+count,false,data); h->SetName(names[i]);
    std::cout << sample <<  " " << h->Integral(h->FindBin(75.),h->FindBin(105.)) << "\n";

    if (i==0) {
      h->SetMarkerSize(1.2); h->SetMarkerStyle(20); h->SetLineWidth(1); h->SetLineColor(colors[i]); h->SetFillColor(0);
      leg->AddEntry(h,legends[i],"P");
      h->GetYaxis()->SetTitleOffset(1.45);
      h->GetXaxis()->SetTitle(xaxisname);
      if (norm == 0) { h->GetYaxis()->SetTitle("Events / bin"); } else { h->GetYaxis()->SetTitle("a. u."); }
      if (logy) {
        gPad->SetLogy();
        h->GetYaxis()->SetRangeUser(1.,1000.*h->GetBinContent(h->GetMaximumBin()));
        if (var.Contains("tau") || var.Contains("eta") || var.Contains("ecf")) { h->GetYaxis()->SetRangeUser(1.,100000.*h->GetBinContent(h->GetMaximumBin())); }
      }
      else { h->GetYaxis()->SetRangeUser(0.,1.3*h->GetBinContent(h->GetMaximumBin())); }
      h->Draw("P E0");
      histoCol.push_back(h);
      h_data = (TH1F*)h->Clone("h_"+name+"_data");
    }
    else {
      h->SetFillColor(colors[i]); h->SetMarkerSize(0); h->SetLineColor(colors[i]);
      leg->AddEntry(h,legends[i],"FL");
      histoCol.push_back(h);
      h->Draw("HIST E0 sames");
    }

    if (i==1) { h_sm = (TH1F*)h->Clone("h_sm"); hs->Add(h); }
    if (i>1 ) { h_sm->Add(h); hs->Add(h); }

  } // end of looping over the trees

  // norm MC
  float normVal = h_data->Integral()/h_sm->Integral();
  if (norm !=0) { h_sm->Scale(normVal); }

  // data mc ratio
  TH1F *h_r = (TH1F*)h_data->Clone("h_"+name+"_r"); h_r->Divide(h_data,h_sm);

  TList *histos = hs->GetHists();
  TIter next(histos);
  TH1F *hist;
  while ((hist =(TH1F*)next())) { hist->Scale(normVal); }

  // continue drawing
  pMain->cd();
  hs->Draw("HIST E0 sames");
  h_data->Draw("P E0 sames");
  leg->Draw("sames");
  c_->RedrawAxis();
  pMain->RedrawAxis();

  pRatio->cd();
  gPad->SetBottomMargin(0.2);
  h_r->GetYaxis()->SetTitleOffset(0.9);
  h_r->GetYaxis()->SetTitleSize(0.1);
  h_r->GetYaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitleSize(0.1);
  h_r->GetXaxis()->SetLabelSize(0.08);
  h_r->GetXaxis()->SetTitle(xaxisname);
  h_r->GetYaxis()->SetTitle("Data / MC");
  h_r->GetYaxis()->SetRangeUser(0.8,1.2);
  h_r->Draw("P E0");


  TCanvas *c_signif = new TCanvas("c_signif_"+name,"c_signif_"+name,600,600); c_signif->SetName("c_signif_"+name);
  TPad *pUpper  = new TPad("pUpper_"+name,"pUpper+"+name,0.0,0.5,1.0,1.0);
  pUpper->SetRightMargin(0.05);
  pUpper->SetLeftMargin(0.17);
  pUpper->SetBottomMargin(0.);
  pUpper->SetTopMargin(0.1);

  TPad *pLower = new TPad("pLower_"+name,"pLower_"+name,0.0,0.0,1.0,0.5);
  pLower->SetRightMargin(0.05);
  pLower->SetLeftMargin(0.17);
  pLower->SetTopMargin(0.);
  pLower->SetBottomMargin(0.37);
  pUpper->Draw();
  pLower->Draw();

  // draw expected distributions but QCD
  pUpper->cd();
  histoCol[2]->GetYaxis()->SetTitle("Total SM-QCD");
  histoCol[2]->GetYaxis()->SetRangeUser(0.,1.2*histoCol[2]->GetBinContent(histoCol[2]->GetMaximumBin()));
  histoCol[2]->GetXaxis()->SetTitle(xaxisname);
  histoCol[2]->Draw("L E0"); //top
  histoCol[3]->Draw("L E0 sames"); // w
  histoCol[4]->Draw("L E0 sames"); // z->cat
  histoCol[6]->Draw("L E0 sames"); // h->cat
  pUpper->RedrawAxis();

  // draw S / sqrt(QCD)
  pLower->cd();
  TH1F *h_sovsqrtb_top = createSovSqrtB("h_sovsqrtb_top",histoCol[2],histoCol[1]); h_sovsqrtb_top->SetFillColor(0);
  h_sovsqrtb_top->GetYaxis()->SetTitle("S / #surd(QCD)");
  h_sovsqrtb_top->GetYaxis()->SetRangeUser(0.,1.2*h_sovsqrtb_top->GetBinContent(h_sovsqrtb_top->GetMaximumBin()));
  h_sovsqrtb_top->GetXaxis()->SetTitle(xaxisname);
  h_sovsqrtb_top->Draw("L"); //top
  TH1F *h_sovsqrtb_w = createSovSqrtB("h_sovsqrtb_w",histoCol[3],histoCol[1]); h_sovsqrtb_w->SetFillColor(0);
  h_sovsqrtb_w->Draw("L sames"); // w
  TH1F *h_sovsqrtb_zcat= createSovSqrtB("h_sovsqrtb_zcat",histoCol[4],histoCol[1]); h_sovsqrtb_zcat->SetFillColor(0);
  h_sovsqrtb_zcat->Draw("L sames"); // z->cat
  TH1F *h_sovsqrtb_h = createSovSqrtB("h_sovsqrtb_h",histoCol[6],histoCol[1]); h_sovsqrtb_h->SetFillColor(0);
  h_sovsqrtb_h->Draw("L sames"); // h->cat
  pLower->RedrawAxis();

  const int dir_err = system("mkdir -p ./"+dir);
  if (-1 == dir_err) {
    printf("Error creating directory!n");
    exit(1);
  }
  c_->Print(dir+"/"+name+"_"+logstr+".png");
  c_->Print(dir+"/"+name+"_"+logstr+".pdf");

}




TH1F *createSovSqrtB(TString name, TH1F *h_s, TH1F *h_b) {
  TH1F *h_sovsqrtb = (TH1F*)h_s->Clone(name);
  for (int i0=1; i0<=h_s->GetNbinsX(); ++i0) { 
    h_sovsqrtb->SetBinContent(i0, (h_s->GetBinContent(i0)/(sqrt(h_b->GetBinContent(i0)))) ); 
    h_sovsqrtb->SetBinError(i0,0.);
  }
  return h_sovsqrtb;
}




TH1F *create1Dhisto(TString sample_,TTree *tree,TString intLumi,TString cuts,TString branch,int bins,float xmin,float xmax,
                    bool useLog,int color, int style,TString name,bool norm,bool data) {
  TH1::SetDefaultSumw2(kTRUE);

  std::cout << "sample = " << sample_ << " " << name << "\n";
  TString sample;
  /* 
 if (sample_.Contains("signal"))     { sample = "qcd"; }
  if (sample_.Contains("higgs"))      { sample = "qcd"; }
  if (sample_.Contains("zboson"))     { sample = "qcd"; }
  if (sample_.Contains("wboson"))     { sample = "qcd"; }
  if (sample_.Contains("top"))        { sample = "qcd"; }
  if (sample_.Contains("tt"))        { sample = "qcd"; }
  if (sample_.Contains("st"))        { sample = "qcd"; }
  if (sample_.Contains("wqq"))        { sample = "qcd"; }
  if (sample_.Contains("vv"))        { sample = "qcd"; }
  if (sample_.Contains("sm"))         { sample = "qcd"; }
  if (sample_.Contains("znunu"))      { sample = "qcd"; }
  if (sample_.Contains("znunu-rwgt")) { sample = "qcd-rwgt"; }
  if (sample_.Contains("data"))       { sample = "qcd"; }
  if (sample_.Contains("qcd"))        { sample = "qcd"; }
  if (sample_.Contains("qcd-norwgt")) { sample = "qcd-norwgt"; }
*/

  TString cut;
  /*
  if (sample == "photons") {
    if (data) { cut ="(passPhoton165_HE10 && "+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*(passPhoton165_HE10 &&"+cuts+")"; }
  }
  */
  /*
  if (sample == "qcd") {
    if (data) { cut ="(passHTTrig && "+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*(passHTTrig &&"+cuts+")"; }
  }
  */
  /*
  if (sample == "qcd-norwgt") {
    if (data) { cut ="("+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }                                                                                                             
  }
  
  if (sample == "qcd") {
    if (data) { cut ="("+cuts+")"; }
    //    else      { cut ="(rewgtfunc(fj_1_pt)*xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*topptWeight*"+intLumi+")*("+cuts+")"; }
    //else      { cut ="(rewgtfunc(fj_1_sdmass)*xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }
    //else      { cut ="(rewgtfunc(sqrt((fj_1_sj1_eta-fj_1_sj2_eta)*(fj_1_sj1_eta-fj_1_sj2_eta) + (fj_1_sj1_eta-fj_1_sj2_phi)*(fj_1_sj1_eta-fj_1_sj2_phi)))*xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }                                                                                                             
  }
*/
   if (data) { cut ="("+cuts+")"; }
    //    else      { cut ="(rewgtfunc(fj_1_pt)*xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }
   else      { cut ="(genWeight*xsecWeight*puWeight*topptWeight*"+intLumi+")*("+cuts+")"; }
   //   else      { cut ="(xsecWeight*"+intLumi+")*("+cuts+")"; }
   /*
  if (sample == "signal") {
    if (data) { cut ="("+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }
  }
   */
   /*
   std::cout << " get xsec weight\n" << "\n";
   TH1F *hxsec = new TH1F("hxsec_"+name,"hxsec_"+name,10000,0.,100);
   tree->Project("hxsec_"+name,"xsecWeight",cut);
   std::cout << hxsec->GetMean() << "\n";
   */
   std::cout << " cut = " << cut << "\n";

  TH1F *hTemp = new TH1F(name,name,bins,xmin,xmax);
  tree->Project(name,branch,cut);

  hTemp->SetLineWidth(3);
  hTemp->SetMarkerSize(0);
  hTemp->SetLineColor(color);
  hTemp->SetFillColor(color);
  hTemp->SetLineStyle(style);

  double error =0.; double integral = hTemp->IntegralAndError(bins,bins+1,error);
  hTemp->SetBinContent(bins,integral);
  hTemp->SetBinError(bins,error);

  std::cout << "sample = " << sample_ << hTemp->Integral() << "\n";

  if (norm) { hTemp->Scale(1./(hTemp->Integral())); }

  return hTemp;
} //~ end of create1Dhisto


TH1F *rescaleXaxis(TH1F *inputhisto, float xmin, float xmax) {
  TH1::SetDefaultSumw2(kTRUE);

  int nbins = inputhisto->GetNbinsX();
  TString name = inputhisto->GetName();
  TH1F *outputhisto = new TH1F(name,name,nbins,xmin,xmax);

  for (unsigned int i0=0; i0<nbins; ++i0) {
    outputhisto->SetBinContent(i0+1,inputhisto->GetBinContent(i0+1));
    outputhisto->SetBinError(i0+1,inputhisto->GetBinError(i0+1));
  }

  return outputhisto;
}// end of rescaleXaxis



void rescaleXaxis(TGraphAsymmErrors *g, double xmin, double scaleX) {
  TH1::SetDefaultSumw2(kTRUE);

  int N=g->GetN();
  double *x=g->GetX();
  double *xerrh=g->GetEXhigh();
  double *xerrl=g->GetEXlow();
  int i=0; while(i<N) {
    x[i]=xmin+x[i]*scaleX;
    xerrh[i]=0.;
    xerrl[i]=0.;
    i=i+1;
  }
  g->GetHistogram()->Delete();
  g->SetHistogram(0);
}



std::vector<unsigned int> getColors() {
  std::vector<unsigned int> colors_; colors_.clear();
  colors_.push_back(92);
  colors_.push_back(4);
  colors_.push_back(6);
  colors_.push_back(182);
  colors_.push_back(8);
  colors_.push_back(2);
  return colors_;
}


std::vector<unsigned int> getStyles() {
  std::vector<unsigned int> styles_; styles_.clear();
  styles_.push_back(20);
  styles_.push_back(21);
  styles_.push_back(34);
  styles_.push_back(47);
  styles_.push_back(23);
  return styles_;
}

TH1F *getDataMCratio(TGraphAsymmErrors *indata, TH1F *inMC) {
  TH1::SetDefaultSumw2(kTRUE);

  TH1F *h_data = (TH1F*)inMC->Clone("h_data"); h_data->SetName("h_data");
  for (unsigned int i0=0; i0<inMC->GetNbinsX(); ++i0) {
    h_data->SetBinContent(i0+1, indata->GetY()[i0]);
    h_data->SetBinError(i0+1, (indata->GetEYlow()[i0]+indata->GetEYhigh()[i0])/2.);
  }

  TH1F *h_ratio = (TH1F*)h_data->Clone("h_ratio");
  h_ratio->Divide(h_data,inMC);
  h_ratio->SetMarkerStyle(20); h_ratio->SetMarkerSize(1);
  h_ratio->SetLineWidth(2); h_ratio->SetLineStyle(1);

  return h_ratio;
}

