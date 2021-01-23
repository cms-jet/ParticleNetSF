#ifndef AUXILIARYFUNCTIONS_INCLUDE
#define AUXILIARYFUNCTIONS_INCLUDE

#include <vector>
#include <string>
#include <unistd.h>
#include <iostream>
#include <bits/stdc++.h>
#include <unistd.h>
#include "TString.h"

namespace aux {

  TH1F *create1Dhisto(TString sample_,TTree *tree,TString intLumi,TString cuts,TString branch,int bins,float xmin,float xmax,
		      bool useLog,int color, int style, bool fillHisto, TString name,bool norm,bool data) {
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
    if (fillHisto) { hTemp->SetFillColor(color); }
    hTemp->SetLineStyle(style);

    double error =0.; double integral = hTemp->IntegralAndError(bins,bins+1,error);
    hTemp->SetBinContent(bins,integral);
    hTemp->SetBinError(bins,error);

    std::cout << "sample = " << sample_ << hTemp->Integral() << "\n";

    if (norm) { hTemp->Scale(1./(hTemp->Integral())); }

    return hTemp;
  }


  TCanvas *makeCanvasWithRatioMC(std::vector<TH1F*> h1d, TString name, TString xname, TString yname, TLegend *leg) {

    TH1::SetDefaultSumw2(kTRUE);
    setTDRStyle();
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(1);

    std::vector<TH1F*> h_r; h_r.clear(); 
    TH1F *h_den = (TH1F*)h1d[0]->Clone(name+"_den");
    for (int i0=0; i0<h1d.size(); ++i0) { 
      
      ostringstream tmpIdx; tmpIdx << i0; TString idx = tmpIdx.str();      
      TH1F *h_r_ = (TH1F*)h1d[i0]->Clone(name+"_r_"+idx);
      h_r_->Divide(h_r_,h_den);
      h_r.push_back(h_r_);

    }

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
      std::cout << " Ploting main canvas " << i0 << "\n";
      if (i0==0) {
	h1d[i0]->GetXaxis()->SetTitle(xname);
	h1d[i0]->GetYaxis()->SetTitle(yname);
	h1d[i0]->GetYaxis()->SetRangeUser(0.,1.5*h1d[i0]->GetMaximum());
	h1d[i0]->Draw("HIST E0");
      }
      else { h1d[i0]->Draw("HIST E0 sames"); }
    }
    leg->Draw("sames");
    
    pRatio->cd();
    for (int i0=0; i0<h_r.size(); ++i0) {
      if (i0==0) {
        h_r[i0]->GetXaxis()->SetTitle(xname);
        h_r[i0]->GetYaxis()->SetTitle("Ratio");
        h_r[i0]->GetYaxis()->SetRangeUser(0.3,1.7);
        h_r[i0]->Draw("HIST");
      }
      else { h_r[i0]->Draw("HIST E0 sames"); }
    }

    return c_;
  }


}

#endif
