#ifndef CONFIGURATION_INCLUDE
#define CONFIGURATION_INCLUDE

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
  
  void configuration(TString sample) {

    TString jet_prefix;

    // =================== area to modify - tune ===================== //
   
    if ( (sample == "tt1L") || (sample=="ttbar1L") || (sample=="ttbar1l") || (sample == "tt1l") || (sample == "tt1lw") ) {
      jet_prefix = "ak8_1_";
      brX = jet_prefix+"corr_sdmass"; // brX = "ak8_1_mass";
      brY = jet_prefix+"pt";
    }

    if (sample == "zqq") {
      jet_prefix = "fj_1_";
      brX = jet_prefix+"sdmass"; // brX = "ak8_1_corr_sdmass";
      brY = jet_prefix+"pt";
    }

    algo      = "particlenetmd";
    score_def = jet_prefix+"ParticleNetMD_XbbVsQCD";
    binsX = 30; minX = 50;  maxX = 200.;
    binsY = 40; minY = 200; maxY = 1200.;
   
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

} 

#endif
