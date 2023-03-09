#include <string>
#include "configuration.h"

void makeOneDatacard(TString inputname, TString category, TString wpmin, TString wpmax, std::string name_, TString sample, TString era);
void makeOneDatacardTop(TString inputname, TString category, TString wpmin, TString wpmax,TString name, TString sample, TString era);

void makeDatacards(TString era, TString sample, TString category, TString wpmin, TString wpmax) {
  conf::configuration(sample);

  std::vector<TString> name  = conf::name;
  if (sample=="tt1L" || sample=="ttbar1L" || (sample=="ttbar1l") || (sample=="tt1l") ) 
    {
      for (int i0=0; i0<name.size(); ++i0) 
	{
	  makeOneDatacardTop(conf::algo+"_sf",category,wpmin,wpmax,name[i0],sample,era);
	}
    }
}


void makeOneDatacardTop(TString inputname, TString category, TString wpmin, TString wpmax, TString name_, TString sample, TString era) {
  conf::configuration(sample);

  std::cout << "makeData card : " << name_ << "\n"; 

  TString label0;
  TString name = (TString)name_;
  TString inputname_ = (TString)inputname;
  if (inputname_.Contains("tau21ddt"))                                                           { label0 = "tau21ddt"; }
  if (inputname_.Contains("dak8ddt"))                                                            { label0 = "dak8ddt"; }
  if (inputname_.Contains("dak8md"))                                                             { label0 = "dak8md"; }
  if (inputname_.Contains("dak8") && !(inputname_.Contains("md") || inputname_.Contains("ddt"))) { label0 = "dak8"; }
  if (inputname_.Contains("particlenetmd"))                                                      { label0 = "particlenetmd"; }
  if (inputname_.Contains("particlenet_"))                                                        { label0 = "particlenet"; }
  
  std::cout << label0 << "\n";
  
  const int dir_err = system("mkdir -p ./"+inputname_+"/fitdir/");
  if (-1 == dir_err) { printf("Error creating directory!n"); exit(1); }

  std::ofstream out("./"+inputname_+"/fitdir/datacard_"+label0+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt");
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::cout.rdbuf(out.rdbuf());
  TString p2name;
  TString p3name;
  TString othername;
  if (category == "top") {
    p2name = "tp2";
    p3name = "tp3";
    othername = "other";
  }
  if (category == "w") {
    p2name = "tp2";
    p3name = "tp3";
    othername ="other";
  }
  
  // pass 
  TFile *f_p = TFile::Open((TString)inputname+"/"+label0+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+"_templates_p.root","READONLY");
  TH1D  *h_obs_p    = (TH1D*)f_p->Get("data_obs"); 
  TH1D  *h_top_p1_p = (TH1D*)f_p->Get("tp1");
  TH1D  *h_top_p2_p = (TH1D*)f_p->Get(p2name);
  TH1D  *h_top_p3_p = (TH1D*)f_p->Get(p3name);
  TH1D  *h_other_p  = (TH1D*)f_p->Get("other");
  // fail
  TFile *f_f = TFile::Open((TString)inputname+"/"+label0+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+"_templates_f.root","READONLY");
  TH1D  *h_obs_f    = (TH1D*)f_f->Get("data_obs");   
  TH1D  *h_top_p1_f = (TH1D*)f_f->Get("tp1");
  TH1D  *h_top_p2_f = (TH1D*)f_f->Get(p2name);
  TH1D  *h_top_p3_f = (TH1D*)f_f->Get(p3name);
  TH1D  *h_other_f  = (TH1D*)f_f->Get("other");
  
  std::cout << "imax 2  number of channels\n";
  std::cout << "jmax 3  number of processes -1\n";
  std::cout << "kmax *  number of nuisance parameters (sources of systematical uncertainties)\n";
  std::cout << "------------\n";

  std::cout << "shapes  *  pass  " << (TString)inputname << "/"+label0+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_" << name << "_templates_p.root  $PROCESS $PROCESS_$SYSTEMATIC\n";
  std::cout << "shapes  *  fail  " << (TString)inputname << "/"+label0+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_" << name << "_templates_f.root  $PROCESS $PROCESS_$SYSTEMATIC\n";
  std::cout << "------------\n";


  std::cout << "bin             pass fail    \n";
  std::cout << "observation      " << h_obs_p->Integral(1,h_obs_p->GetNbinsX()) << " " << h_obs_f->Integral(1,h_obs_f->GetNbinsX()) << "\n";
  std::cout << "------------\n";
  std::cout << "# now we list the expected events for signal and all backgrounds in that bin\n";
  std::cout << "# the second 'process' line must have a positive number for backgrounds, and 0 for signal\n";
  std::cout << "# then we list the independent sources of uncertainties, and give their effect (syst. error)\n";
  std::cout << "# on each process and bin\n";
  std::cout << "bin             pass pass pass pass     fail fail fail fail  \n";
  std::cout << "process         " << p3name << " " << p2name << " tp1 other     " << p3name << " " << p2name << " tp1 other  \n";  
  std::cout << "process         4 3 6 -7      4 3 6 -7   \n";
  std::cout << "rate            -1 -1 -1 -1     -1 -1 -1 -1  \n"; 
  std::cout << "------------\n";


  std::cout << "lumi            lnN      1.025 1.025 1.025 1.025     1.025 1.025 1.025 1.025 \n";
  std::cout << "\n";
    
  std::cout << "tp3_xsec     lnN      1.05 - - -    1.05 - - - \n";
  //std::cout << "wqq_xsec     lnN      - 1.05 - -    - 1.05 - - \n";
  std::cout << "tp2_xsec     lnN      - 1.05 - -    - 1.05 - - \n";
  std::cout << "tp1_xsec     lnN      - - 1.05 -    - - 1.05 - \n";
  std::cout << "other_xsec   lnU      - - - 2.00    - - - 2.00 \n";
  
  // std::cout << "\n";
  // std::cout << "other_sf   lnN      - - - -    - - - 5.00 \n";
  

  std::cout << "\n";

  
  std::cout << "tp3jms       shapeU    1 - - -     1 - - - \n";  
  std::cout << "vqqjms       shapeU    - 1 - -     - 1 - - \n";
  std::cout << "tp1jms       shapeU    - - 1 -     - - 1 - \n";
  std::cout << "otherjms       shapeU    - - - 1     - - - 1 \n";

  std::cout << "tp3jmr       shapeU    1 - - -     1 - - - \n";  
  std::cout << "vqqjmr       shapeU    - 1 - -     - 1 - - \n";
  std::cout << "tp1jmr       shapeU    - - 1 -     - - 1 - \n";
  std::cout << "otherjmr       shapeU    - - - 1     - - - 1 \n";
  
  std::cout << "pu          shape    1 1 1 1     1 1 1 1 \n";
  std::cout << "jes         shape    1 1 1 1     1 1 1 1 \n";
  std::cout << "jer         shape    1 1 1 1     1 1 1 1 \n";
  std::cout << "met         shape    1 1 1 1     1 1 1 1 \n";
  std::cout << "lhescalemuf shape    1 1 1 1     1 1 1 1 \n";
  std::cout << "lhescalemur shape    1 1 1 1     1 1 1 1 \n";
  //std::cout << "lhepdf      shape    1 1 1 1     1 1 1 1 \n";
  
  //  std::cout << "herwig         shape    1 1 1 -     1 1 1 -\n";

  std::cout << "norm_top    rateParam    pass    " << p3name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    fail    " << p3name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    pass    " << p2name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    fail    " << p2name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    pass    tp1      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    fail    tp1      1   [0.,10.]\n";
  //std::cout << "norm_other    rateParam    pass    other      1   [0.,10.]\n";
  //std::cout << "norm_other    rateParam    fail    other      1   [0.,10.]\n";
  std::cout << "\n";

  std::cout << "*  autoMCStats  0\n";
  //std::cout << "dummy    lnN    1.001  1.001  1.001  1.001  1.001  1.001\n";

}

