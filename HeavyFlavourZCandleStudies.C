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
TCanvas *c_ratioOf1Dhistos(TString c_name, TString lumi, std::vector<TString> namest, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, 
			   std::vector<TString> namesc, std::vector<TString> cuts,
                           TString xname, TString yname, int nbins, float xmin, float xmax, bool uselog, bool norm);
TH1F *rescaleXaxis(TH1F *inputhisto, float xmin, float xmax);
void rescaleXaxis(TGraphAsymmErrors *inputhisto, double xmin, double scale);
TH1F *getDataMCratio(TGraphAsymmErrors *indata, TH1F *inMC);
void makeDataMCPlotFromCombine(TString path2file, TString era, TString category, TString category_extra, TString wpmin, TString wpmax, TString name, TString passOrFail, 
			       float xmin, float xmax, int nbins,TString xaxisname, bool log, TString sample);
TH1D *h1DHistoFrom2DTemplates(TString path2file,TString h2dname,TString name,double ymin, double ymax,int color,bool isdata);
TCanvas *makeCanvasWithRatio(TH1D* hdata, std::vector<TH1D*> h1d, TString dir, TString name, TString xname, TString yname, TLegend *leg);
void makeDataMCFrom2DTemplatesTop(TString path2file, TString nameoutfile, TString name, TString sample, double ymin, double ymax);
TH1D *h1dShift(TString name, TH1D* inputh1d, int numOfBins, bool up);

TH1D *hrwgt = nullptr;

double rewgtfunc(double x){
  if (!hrwgt) return 1;
  return hrwgt->GetBinContent(hrwgt->FindFixBin(x));
}

void HeavyFlavourZCandleStudies(TString era, TString sample, TString category, TString wpmin, TString wpmax, bool postfit=false, TString passOrFail="pass", TString category_extra="") 
{
  conf::configuration(sample);
  
  std::vector<TString> name  = conf::name;
  std::vector<double>  ptmin = conf::ptmin;
  std::vector<double>  ptmax = conf::ptmax;    
  
  if (postfit) 
    {
      if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") || (sample == "comb")) 
	{ 
	  std::cout << "postfit: " << sample << " " << category << " " << wpmin << "to" << wpmax << " " << passOrFail << "\n";
	  for (int i0=0; i0<name.size(); ++i0) 
	    { 
	      makeDataMCPlotFromCombine(conf::algo,era,category,category_extra,wpmin,wpmax,name[i0],passOrFail,conf::minX,conf::maxX,conf::binsX,"m_{SD} [GeV]",false,sample); 
	    }
	}
    }
  if (!postfit) 
    {
      if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") ) 
	{ 
	  std::cout << sample << "\n";
	  for (int i0=0; i0<name.size(); ++i0) 
	    { 
	      std::cout << sample<< "\n";
	      makeDataMCFrom2DTemplatesTop("templates2D/"+conf::algo+"_"+sample+"_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_200to1200_templates.root",                       
				     conf::algo+"_"+sample+"_"+category+"_"+wpmin+"to"+wpmax+"_"+era,name[i0],category,ptmin[i0],ptmax[i0]);
	    } 
	}
    }
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
  if (path2file.Contains("particlenetmd"))      { dirname1 = "particlenetmd_sf"; }
  else if (path2file.Contains("particlenet"))   { dirname1 = "particlenet_sf";   }

  const int dir_err = system("mkdir -p ./"+dirname1);
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }
  nameoutfile = nameoutfile+"_"+name;

  // Data templates 
  TH1D *h_data_p  = h1DHistoFrom2DTemplates(path2file,"data_obs_pass",name,ymin,ymax,1,true);
  TH1D *h_data_f  = h1DHistoFrom2DTemplates(path2file,"data_obs_fail",name,ymin,ymax,1,true);

  // MC templates
  std::vector<TString> syst = conf::syst;
  std::vector<TString> processes_in = conf::processes_in;

  std::vector<TString> processes_out;
  if (sample == "top") { processes_out.push_back("tp3"); processes_out.push_back("tp2"); processes_out.push_back("tp1"); processes_out.push_back("other"); }
  if (sample == "w")   { processes_out.push_back("tp3"); processes_out.push_back("tp2"); processes_out.push_back("tp1"); processes_out.push_back("other"); }
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
      for (unsigned int ip=0; ip<processes_in.size(); ip+=1) {
	if (syst[is]=="_") {
	  std::cout << processes_in[ip] << " " << ip << "\n";
          /////////////////////////////////////////////////////////////////////////
          
          TH1D *h_;
          if (ip>8){
               h_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
               TH1D *h__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
               h_->Add(h__);
          }
          else if(ip<=8){       
	       h_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
	       TH1D *h__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
	       TH1D *h___ = h1DHistoFrom2DTemplates(path2file,processes_in[ip+2]+"_"+syst[is]+catstr_,name,ymin,ymax,1,false);
	       h_->Add(h__); h_->Add(h___);
          }
          
          /////////////////////////////////////////////////////////////////////////
	  h_->SetName(processes_out[count]); h_->SetLineColor(colors[count]); h_->SetFillColor(colors[count]);
	 
	  if (category[ic] == "pass") { hist_out_p.push_back(h_); hist_out_nom_p.push_back(h_); }
	  else                        { hist_out_f.push_back(h_); hist_out_nom_f.push_back(h_); } 
	}
      
	else {

	  TString nameSyst = syst[is];
	  if ( (syst[is].Contains("jms")) || (syst[is].Contains("jmr")) ) {
	    nameSyst = processes_out[count]+syst[is];
	    if ( (processes_in[ip].Contains("wqq")) || (processes_in[ip].Contains("zqq")) || (processes_in[ip].Contains("tt_p2")) ) { nameSyst = "vqq"+syst[is]; }
	    if ( (processes_in[ip].Contains("wll")) || (processes_in[ip].Contains("qcd")) ) { nameSyst = "other"+syst[is]; } 
	  }
          if (ip>8){
          TH1D *h_up_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
          TH1D *h_up__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
          h_up_->Add(h_up__);
          h_up_->SetName(processes_out[count]+"_"+nameSyst+"Up"); h_up_->SetLineColor(colors[count]); h_up_->SetFillColor(colors[count]);
          if (category[ic] == "pass") { hist_out_p.push_back(h_up_); } else { hist_out_f.push_back(h_up_); }

          TH1D *h_down_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
          TH1D *h_down__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
          h_down_->Add(h_down__); 
          h_down_->SetName(processes_out[count]+"_"+nameSyst+"Down"); h_down_->SetLineColor(colors[count]); h_down_->SetFillColor(colors[count]);
          if (category[ic] == "pass") { hist_out_p.push_back(h_down_); } else { hist_out_f.push_back(h_down_); }

          }
	  else if (ip<=8){
	  TH1D *h_up_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
	  TH1D *h_up__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
	  h_up_->Add(h_up__);
	  TH1D *h_up___ = h1DHistoFrom2DTemplates(path2file,processes_in[ip+2]+"_"+syst[is]+"Up_"+catstr_,name,ymin,ymax,1,false);
	  h_up_->Add(h_up___); 
	  h_up_->SetName(processes_out[count]+"_"+nameSyst+"Up"); h_up_->SetLineColor(colors[count]); h_up_->SetFillColor(colors[count]);
	  if (category[ic] == "pass") { hist_out_p.push_back(h_up_); } else { hist_out_f.push_back(h_up_); }

	  TH1D *h_down_   = h1DHistoFrom2DTemplates(path2file,processes_in[ip]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
	  TH1D *h_down__  = h1DHistoFrom2DTemplates(path2file,processes_in[ip+1]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
	  h_down_->Add(h_down__); 
	  TH1D *h_down___ = h1DHistoFrom2DTemplates(path2file,processes_in[ip+2]+"_"+syst[is]+"Down_"+catstr_,name,ymin,ymax,1,false);
	  h_down_->Add(h_down___);
	  h_down_->SetName(processes_out[count]+"_"+nameSyst+"Down"); h_down_->SetLineColor(colors[count]); h_down_->SetFillColor(colors[count]);
	  if (category[ic] == "pass") { hist_out_p.push_back(h_down_); } else { hist_out_f.push_back(h_down_); }
	  }
        }
	ip = ip+2;
	++count;
      }
    } // end of looping over the syst 
  } // end of category loop 


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
  fout_p->Close();

  TFile *fout_f = new TFile("./"+dirname1+"/"+nameoutfile+"_templates_f.root","RECREATE");
  h_data_f->Write("data_obs");
  for (unsigned int i0=0; i0<hist_out_f.size(); ++i0) {
    hist_out_f[i0]->Write(hist_out_f[i0]->GetName());
  }
  fout_f->Close();
  
  hist_out_p.clear(); hist_out_nom_p.clear();
  hist_out_f.clear(); hist_out_nom_f.clear();

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
  h1d->SetName(h2dname+"_"+name);
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


TH1D *h1dShift(TString name, TH1D* inputh1d, int numOfBins, bool up) {
  // need to be fixed for case with numOfBins>1 
 
  TH1D *hshifted = (TH1D*)inputh1d->Clone("hshifted_"+name);
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

void makeDataMCPlotFromCombine(TString path2file, TString era, TString category, TString category_extra, TString wpmin, TString wpmax, TString name, TString passOrFail, 
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

  TString fdiag_ = "./"+path2file+"_sf/fitdir/fitdiagnostics_"+path2file+"_"+sample+"_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";
  std::cout << " Opening fitdiagnostics: " << fdiag_ << "\n";
  TFile *fdiag = TFile::Open(fdiag_, "READONLY" );

  std::vector<TString> processes; processes.clear();
  std::vector<int> colors;        colors.clear();
  std::vector<TString> legends;   legends.clear();
  if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") ) 
    {
      if (category == "w") 
	{ 
	  std::cout << " sample = " << sample << " , category = " << category << "\n";
	  processes.push_back("tp3"); processes.push_back("tp2"); processes.push_back("tp1"); processes.push_back("other"); processes.push_back("total");
	  colors.push_back(conf::tp3.color); colors.push_back(conf::tp2.color); colors.push_back(conf::tp1.color); colors.push_back(conf::other.color); colors.push_back(602);
	  legends.push_back(conf::tp3.legend_name); legends.push_back(conf::tp2.legend_name); legends.push_back(conf::tp1.legend_name); legends.push_back(conf::other.legend_name); legends.push_back("Total SM");
	}
      else if (category == "top")
	{ 
	  std::cout << " sample = " << sample << " , category = " << category << "\n";
	  processes.push_back("tp3"); processes.push_back("tp2"); processes.push_back("tp1"); processes.push_back("other"); processes.push_back("total");
	  colors.push_back(conf::tp3.color); colors.push_back(conf::tp2.color); colors.push_back(conf::tp1.color); colors.push_back(conf::other.color); colors.push_back(602);
	  legends.push_back(conf::tp3.legend_name); legends.push_back(conf::tp2.legend_name); legends.push_back(conf::tp1.legend_name); legends.push_back(conf::other.legend_name); legends.push_back("Total SM");
	}
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
  if (path2file.Contains("2016")) { longstring = "19.52 fb^{-1} (13 TeV)"; }
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
    c->Print(path2file+"_sf/plots_datamc/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wpmin+"to"+wpmax+"_"+passOrFail+"_log.pdf");
    c->Print(path2file+"_sf/plots_datamc/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wpmin+"to"+wpmax+"_"+passOrFail+"_log.png");
  } 
  else {
    c->Print(path2file+"_sf/plots_postfit/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wpmin+"to"+wpmax+"_"+passOrFail+"_lin.pdf");
    c->Print(path2file+"_sf/plots_postfit/"+path2file+"_"+sample+"_"+era+"_"+name+"_"+category+"_"+wpmin+"to"+wpmax+"_"+passOrFail+"_lin.png");
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
	ratios[ihistos][ihisto]->Draw("HIST E0");
      }
      else { ratios[ihistos][ihisto]->Draw("HIST E0 sames"); }
    }
  }

  c_r->Print("./plots/"+c_name+".pdf");
  return c_r;
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
  TString cut;
  if (data) { cut ="("+cuts+")"; }
  else      { cut ="(genWeight*xsecWeight*puWeight*topptWeight*"+intLumi+")*("+cuts+")"; }

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

