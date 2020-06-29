#include <iostream>   // std::cout
#include <string>     // std::string, std::stof
void makeOneDatacard(std::string inputname, std::string name);
void makeOneDatacardTop(std::string inputname, std::string name, TString sample);


void makeDatacards(TString sample) {

  std::vector<std::string> name;
  // from ggH
  //  name.push_back("pt450to500"); name.push_back("pt500to550"); name.push_back("pt550to600"); name.push_back("pt600to675"); name.push_back("pt675to800"); name.push_back("pt800to1200");
  name.push_back("pt450to600"); name.push_back("pt600to800"); name.push_back("pt800to1200");

  if (sample=="zqq") { std::cout << sample << "\n";
    for (int i0=0; i0<name.size(); ++i0) { makeOneDatacard("../particlenet_sf",name[i0]); }  }

  if (sample=="tt1L" || sample=="ttbar1L") {
    makeOneDatacardTop("../particlenet_sf","pt200to400","w");
    for (int i0=0; i0<name.size(); ++i0) { makeOneDatacardTop("../particlenet_sf",name[i0],"top"); } }

}



void makeOneDatacard(std::string inputname, std::string name) {

  std::ofstream out("./fitdir/datacard_"+(TString)name+".txt");
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::cout.rdbuf(out.rdbuf());
  
  // ggH
  std::map<std::string,float> znormEW;
  znormEW["pt450to500"] = 1.05; znormEW["pt500to550"]= 1.05; znormEW["pt550to600"]= 1.07; znormEW["pt600to675"]= 1.07; znormEW["pt675to800"]= 1.07; znormEW["pt800to1200"]= 1.07;
  znormEW["pt450to600"] = 1.05; znormEW["pt500to600"] = 1.06; znormEW["pt600to800"]= 1.07; 
  std::map<std::string,float>::const_iterator znormEW_pos = znormEW.find(name);
  double znormEW_value = znormEW_pos->second; 

  std::map<std::string,float> wznormEW;
  wznormEW["pt450to500"] = 1.02; wznormEW["pt500to550"]= 1.02; wznormEW["pt550to600"]= 1.02; wznormEW["pt600to675"]= 1.06; wznormEW["pt675to800"]= 1.06; wznormEW["pt800to1200"]= 1.06;
  wznormEW["pt450to600"] = 1.02; wznormEW["pt500to600"] = 1.02; wznormEW["pt600to800"] = 1.06;
  std::map<std::string,float>::const_iterator wznormEW_pos = wznormEW.find(name);
  double wznormEW_value = wznormEW_pos->second;
  // ----

  // pass
  TFile *f_p = TFile::Open((TString)inputname+"/particlenet_bb_t_2016_"+(TString)name+"_templates_p.root","READONLY");
  TH1D  *h_obs_p = (TH1D*)f_p->Get("data_obs"); 
  TH1D  *h_qcd_p = (TH1D*)f_p->Get("qcd");
  TH1D  *h_tt_p  = (TH1D*)f_p->Get("tqq");
  TH1D  *h_w_p  = (TH1D*)f_p->Get("wqq");
  TH1D  *h_z_p  = (TH1D*)f_p->Get("zqq");
  
  TFile *f_f = TFile::Open((TString)inputname+"/particlenet_bb_t_2016_"+(TString)name+"_templates_f.root","READONLY");
  TH1D  *h_obs_f = (TH1D*)f_f->Get("data_obs");
  TH1D  *h_qcd_f = (TH1D*)f_f->Get("qcd");
  TH1D  *h_tt_f  = (TH1D*)f_f->Get("tqq");
  TH1D  *h_w_f  = (TH1D*)f_f->Get("wqq");
  TH1D  *h_z_f  = (TH1D*)f_f->Get("zqq");
  //TH1D  *h_h_f  = (TH1D*)f_f->Get("higgs");
  
  std::cout << "imax 2  number of channels\n";
  std::cout << "jmax 3  number of processes -1\n";
  std::cout << "kmax *  number of nuisance parameters (sources of systematical uncertainties)\n";
  std::cout << "------------\n";

  std::cout << "shapes  *  pass  " << (TString)inputname << "/particlenet_bb_t_2016_" << (TString)name << "_templates_p.root  $PROCESS $PROCESS_$SYSTEMATIC\n";
  std::cout << "shapes  *  fail  " << (TString)inputname << "/particlenet_bb_t_2016_" << (TString)name << "_templates_f.root  $PROCESS $PROCESS_$SYSTEMATIC\n";
  std::cout << "------------\n";


  std::cout << "bin             pass fail    \n";
  std::cout << "observation      " << h_obs_p->Integral(1,h_obs_p->GetNbinsX()) << " " << h_obs_f->Integral(1,h_obs_f->GetNbinsX()) << "\n";
  std::cout << "------------\n";
  std::cout << "# now we list the expected events for signal and all backgrounds in that bin\n";
  std::cout << "# the second 'process' line must have a positive number for backgrounds, and 0 for signal\n";
  std::cout << "# then we list the independent sources of uncertainties, and give their effect (syst. error)\n";
  std::cout << "# on each process and bin\n";
  std::cout << "bin             pass pass pass pass     fail fail fail fail  \n";
  std::cout << "process         zqq wqq tqq qcd     zqq wqq tqq qcd  \n";  
  std::cout << "process         -2 3 4 5      -2 3 4 5   \n";
  std::cout << "rate            -1 -1 -1 -1     -1 -1 -1 -1  \n"; 
  std::cout << "------------\n";


  std::cout << "lumi            lnN      1.025 1.025 1.025 1.025     1.025 1.025 1.025 1.025 \n";
  std::cout << "\n";
  
  // higher order corrections on w/z
  std::cout << "znormQ          lnN      1.10 1.10 - -     1.10 1.10 - - \n";
  std::cout << "znormEW         lnN      " << znormEW_value <<  " " << znormEW_value << " - -     " << znormEW_value << " " << znormEW_value << " - - \n";
  std::cout << "wznormEW        lnN      - " << wznormEW_value << " - -      - " << wznormEW_value << " - - \n";

  // xsections
  //  std::cout << "zqq_xsec        lnN      1.2 - - -     1.2 - - - \n";
  //std::cout << "wqq_xsec        lnN      - 1.2 - -     - 1.2 - - \n";
  //std::cout << "tqq_xsec        lnN      - - 1.10 -     - - 1.10 - \n";
  //std::cout << "qcd_xsec        lnN      - - - 1.1     - - - 1.1 \n";
  std::cout << "\n";

  //std::cout << "tqq_sf        lnN      -  - - -     - - 1.10 - \n";
  //std::cout << "wqq_sf        lnN      -  - - -     - 1.20 - - \n";

  
  //std::cout << "pu                                shape    1      1      1        1       1             1       1       1      1         1\n";
  //  std::cout << "zjmsp shape    1 - - -     1 - - -\n";
  std::cout << "wjms shape    - 1 - -     - 1 - -\n";
  //std::cout << "jmsp shape    1 1 - -     - - - -\n";
  //std::cout << "jmsf shape    - - - -     1 1 - -\n";
  
  //  std::cout << "herwig         shape    - - - 1     - - - 1\n";

  /*
  std::cout << "norm_wz    rateParam    pass    zqq      1   [0.,10.]\n";
  std::cout << "norm_wz    rateParam    fail    zqq      1   [0.,10.]\n";
  std::cout << "norm_wz    rateParam    pass    wqq      1   [0.,10.]\n";
  std::cout << "norm_wz    rateParam    fail    wqq      1   [0.,10.]\n";
  */
  /*
  std::cout << "norm_z    rateParam    pass    zqq      1   [0.,100.]\n";
  std::cout << "norm_z    rateParam    fail    zqq      1   [0.,100.]\n";
  std::cout << "norm_w    rateParam    pass    wqq      1   [0.,100.]\n";
  std::cout << "norm_w    rateParam    fail    wqq      1   [0.,100.]\n";
  */
  std::cout << "norm_top   rateParam    pass    tqq      1   [0.,10.]\n";
  std::cout << "norm_top   rateParam    fail    tqq      1   [0.,10.]\n";
  std::cout << "norm_qcd   rateParam    pass    qcd      1   [0.,10.]\n";
  std::cout << "norm_qcd   rateParam    fail    qcd      1   [0.,10.]\n";
  std::cout << "*  autoMCStats  0\n";
  //std::cout << "dummy    lnN    1.001  1.001  1.001  1.001  1.001  1.001\n";

}


void makeOneDatacardTop(std::string inputname, std::string name, TString sample) {

  std::ofstream out("./fitdir/datacard_tt1l_"+(TString)name+".txt");
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::cout.rdbuf(out.rdbuf());
  TString p2name;
  TString p3name;
  TString othername;
  if (sample == "top") {
    p2name = "tp2";
    p3name = "tqq";
    othername = "qcd";
  }
  if (sample == "w") {
    p2name = "wqq";
    p3name = "tp3";
    othername ="other";
  }

  // pass
  TFile *f_p = TFile::Open((TString)inputname+"/particlenet_tt1l_bb_t_2016_"+(TString)name+"_templates_p.root","READONLY");
  TH1D  *h_obs_p    = (TH1D*)f_p->Get("data_obs"); 
  TH1D  *h_top_p1_p = (TH1D*)f_p->Get("tp1");
  TH1D  *h_top_p2_p = (TH1D*)f_p->Get(p2name);
  TH1D  *h_top_p3_p = (TH1D*)f_p->Get(p3name);
  TH1D  *h_other_p  = (TH1D*)f_p->Get("other");
  
  TFile *f_f = TFile::Open((TString)inputname+"/particlenet_tt1l_bb_t_2016_"+(TString)name+"_templates_f.root","READONLY");
  TH1D  *h_obs_f    = (TH1D*)f_f->Get("data_obs");   
  TH1D  *h_top_p1_f = (TH1D*)f_f->Get("tp1");
  TH1D  *h_top_p2_f = (TH1D*)f_f->Get(p2name);
  TH1D  *h_top_p3_f = (TH1D*)f_f->Get(p3name);
  TH1D  *h_other_f  = (TH1D*)f_f->Get("other");
 
  std::cout << "imax 2  number of channels\n";
  std::cout << "jmax 3  number of processes -1\n";
  std::cout << "kmax *  number of nuisance parameters (sources of systematical uncertainties)\n";
  std::cout << "------------\n";

  std::cout << "shapes  *  pass  " << (TString)inputname << "/particlenet_tt1l_bb_t_2016_" << (TString)name << "_templates_p.root  $PROCESS $PROCESS_$SYSTEMATIC\n";
  std::cout << "shapes  *  fail  " << (TString)inputname << "/particlenet_tt1l_bb_t_2016_" << (TString)name << "_templates_f.root  $PROCESS $PROCESS_$SYSTEMATIC\n";
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
  std::cout << "process         2 3 -4 5      2 3 -4 5   \n";
  std::cout << "rate            -1 -1 -1 -1     -1 -1 -1 -1  \n"; 
  std::cout << "------------\n";


  std::cout << "lumi            lnN      1.025 1.025 1.025 1.025     1.025 1.025 1.025 1.025 \n";
  std::cout << "\n";
  
  std::cout << "tqq_xsec     lnN      1.10 - - -    1.10 - - - \n";
  std::cout << "wqq_xsec     lnN      - 1.10 - -    - 1.10 - - \n";
  std::cout << "tp1_xsec     lnN      - - 1.10 -    - - 1.10 - \n";
  std::cout << "other_xsec   lnN      - - - 1.2     - - - 1.2 \n";
  std::cout << "\n";

  std::cout << "other_sfp   lnN      - - - 1.3     - - - - \n";
  std::cout << "other_sff   lnN      - - - -     - - - 1.3 \n";

  std::cout << "wjms shape    - 1 - -     - 1 - -\n";
  //  std::cout << "jmsf shape    - - - -     - 1 - -\n";

  //std::cout << "pu                                shape    1      1      1        1       1             1       1       1      1         1\n";  
  std::cout << "herwig         shape    1 1 1 -     1 1 1 -\n";

  std::cout << "norm_top    rateParam    pass    " << p3name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    fail    " << p3name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    pass    " << p2name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    fail    " << p2name << "      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    pass    tp1      1   [0.,10.]\n";
  std::cout << "norm_top    rateParam    fail    tp1      1   [0.,10.]\n";
  std::cout << "\n";

  std::cout << "*  autoMCStats  0\n";
  //std::cout << "dummy    lnN    1.001  1.001  1.001  1.001  1.001  1.001\n";

}

