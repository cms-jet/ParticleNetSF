#ifndef CONFIGURATION_INCLUDE
#define CONFIGURATION_INCLUDE

#include <vector>
#include <string>
#include <unistd.h>
#include <iostream>
#include <bits/stdc++.h>
#include <unistd.h>
#include "TString.h"

namespace conf {

  struct process {
    
    TString name;
    TString legend_name;
    int     color;
    
  } zqq, zbb, zcc, zll, wqq, wcx, wll, tqq, qcd, tp2, tp3, tp1, other;
  
  TString brX; int binsX; float minX, maxX;
  TString brY; int binsY; float minY, maxY;
  TString algo;
  TString score_def; 
  TString category;
  
  std::vector<TString> name; 
  std::vector<double>  ptmin; 
  std::vector<double>  ptmax;
  std::vector<TString> processes;
  std::vector<TString> process_names;
  std::vector<TString> syst;
  std::vector<TString> processes_in;

  TString path_2016;
  TString path_2017;
  TString path_2018;

  TString jetCone;

  void configuration(TString sample) {

    TString jet_prefix;

    // =================== area to modify - tune ===================== //
   
    if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") || (sample == "tt1lw") ) {
      path_2016 = "/eos/uscms/store/user/pakontax/ParticleNET_UL_NanoV9_For_Copy/2016/from_Huilin/20211212_NanoAODv9_ak8_muon_2015/";
      //path_2016 = "/eos/uscms/store/user/pakontax/ParticleNET_UL_NanoV9_For_Copy/2016/from_Huilin/20211212_NanoAODv9_ak8_muon_2016/";
      path_2017 = "/eos/uscms/store/user/pakontax/ParticleNET_UL_NanoV9_For_Copy/2017/from_Huilin/20211212_NanoAODv9_ak8_muon_2017/";
      path_2018 = "/eos/uscms/store/user/pakontax/ParticleNET_UL_NanoV9_For_Copy/2018/from_Huilin/20211212_NanoAODv9_ak8_muon_2018/";
      
      jetCone    = "ak8";
      jet_prefix = "fj_1_";
      
      // brX = jet_prefix+"regressed_mass"; // brX = "ak8_1_mass" "ak8_corr_sdmass";
      brX = jet_prefix+"sdmass"; // brX = "ak8_1_mass" "ak8_corr_sdmass";
      brY = jet_prefix+"pt";
      category = "top";
      //category = "w";
 
     // for make2DTemplates
      processes.push_back("ttbar-powheg"); process_names.push_back("tt"); 
      processes.push_back("singletop");    process_names.push_back("st");
      processes.push_back("ttv");          process_names.push_back("ttv");
      processes.push_back("w");            process_names.push_back("wll");
      processes.push_back("diboson");      process_names.push_back("vv");
      //processes.push_back("qcd-mg");       process_names.push_back("qcd");

      // keep this format - pretify later - need to keep this order
      processes_in.push_back("tt_p3"); processes_in.push_back("st_p3"); processes_in.push_back("ttv_p3");
      processes_in.push_back("tt_p2"); processes_in.push_back("st_p2"); processes_in.push_back("ttv_p2"); 
      processes_in.push_back("tt_p1"); processes_in.push_back("st_p1"); processes_in.push_back("ttv_p1");
      processes_in.push_back("wll");   processes_in.push_back("vv");    //processes_in.push_back("qcd"); 

      // list of systematic uncertainties
      syst.push_back("_"); syst.push_back("pu"); syst.push_back("jes"); syst.push_back("jer"); syst.push_back("met"); 
      syst.push_back("jms"); syst.push_back("jmr");
      syst.push_back("lhescalemuf"); syst.push_back("lhescalemur"); //syst.push_back("lhepdf");
    }


    if (sample == "zqq") {
      jet_prefix = "fj_1_";
      brX = jet_prefix+"sdmass"; // brX = "ak8_1_corr_sdmass";
      brY = jet_prefix+"pt";
    }

    //algo      = "particlenetmd"; // deepak8ddt particlenetmd 
    algo      = "particlenet"; // deepak8ddt particlenetmd //leo
    score_def = jet_prefix+"ParticleNet_TvsQCD"; // DeepAK8DDT ParticleNetMD_XbbVsQCD
    //score_def = jet_prefix+"ParticleNet_WvsQCD"; // DeepAK8DDT ParticleNetMD_XbbVsQCD //leo
    //score_def = "(fj_1_ParticleNetMD_Xcc+fj_1_ParticleNetMD_Xqq)/(fj_1_ParticleNetMD_Xcc+fj_1_ParticleNetMD_Xqq+fj_1_ParticleNetMD_QCD)";
    binsX = 34; minX = 50;  maxX = 220.;
    binsY = 40; minY = 200; maxY = 1200.; //leo
       
     // Top
    //name.push_back("pt300to400"); ptmin.push_back(300.); ptmax.push_back(400.);
    //name.push_back("pt400to480"); ptmin.push_back(400.); ptmax.push_back(480.);  
    //name.push_back("pt480to600"); ptmin.push_back(480.); ptmax.push_back(600.);   
    name.push_back("pt600to1200"); ptmin.push_back(600.); ptmax.push_back(1200.);
    //name.push_back("200to1200"); ptmin.push_back(200.); ptmax.push_back(1200.);
     //W
    //name.push_back("pt200to300"); ptmin.push_back(200.); ptmax.push_back(300.);
    //name.push_back("pt300to400"); ptmin.push_back(300.); ptmax.push_back(400.);  
    //name.push_back("pt400to800"); ptmin.push_back(400.); ptmax.push_back(800.);   
    
    // =================== end of area to modify - tune ===================== //


    zqq.name        = "zboson";
    zqq.legend_name = "Z";
    zqq.color       = 8;

    zbb.name        = "zbosonbb";
    zbb.legend_name = "Z#rightarrow bb";
    zbb.color       = 8;

    zcc.name        = "zbosoncc";
    zcc.legend_name = "Z #rightarrow cc";
    zcc.color       = 843;
    
    zll.name        = "zbosonll";
    zll.legend_name = "Z #rightarrow ll";
    zll.color       = 432;

    wqq.name        = "wboson";
    wqq.legend_name = "W";
    wqq.color       = 6; 

    wcx.name        = "wbosoncx";
    wcx.legend_name = "W #rightarrow cX";
    wcx.color       = 6; 
    
    wll.name        = "wbosonll";
    wll.legend_name = "W #rightarrow ll";
    wll.color       = 606;
    
    tqq.name        = "tt";
    tqq.legend_name = "Top";
    tqq.color       = 4;
    
    qcd.name        = "qcd";
    qcd.legend_name = "QCD";
    qcd.color       = 92;
    
    tp2.name        = "tp2";
    tp2.legend_name = "W-merged";
    tp2.color       = 590;
    
    tp3.name        = "tp3";
    tp3.legend_name = "Top-merged";
    tp3.color       = 4;
    
    tp1.name        = "tp1";
    tp1.legend_name = "Non-merged";
    tp1.color       = 595;
    
    other.name        = "other";
    other.legend_name = "Other";
    other.color       = 6;
 
  }   

  TString convertFloatToTString(float input) {
    ostringstream tmpInputStr;
    tmpInputStr << input;
    TString inputStr = tmpInputStr.str();
    return inputStr;
  }

  std::string convertTStringToString(TString input) {
    std::string output;
    output = (string)input;
    return output;
  }

} 

#endif
