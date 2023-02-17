#include <bits/stdc++.h>
#include <string>
#include <unistd.h>
#include "configuration.h"

void makeOneFit(std::string era, std::string category, std::string wpmin, std::string wpmax, std::string name, std::string name1);
void makeOneFitTop(std::string era, std::string category, std::string wpmin, std::string wpmax, std::string name, std::string name1);
void makeCombFit(std::string era, std::string category, std::string wpmin, std::string wpmax, std::string name, std::string name1);

void makeFits(std::string era, std::string category, std::string wpmin, std::string wpmax, TString sample) {
  conf::configuration(sample);

  std::vector<TString> name  = conf::name;

  if (sample=="zqq") {
    for (int i0=0; i0<name.size(); ++i0) { makeOneFit(era,category,wpmin,wpmax,(std::string)name[i0],"particlenet"); } }
  
  if (sample=="tt1L" || sample=="ttbar1L" || (sample=="ttbar1l") || (sample=="tt1l")) {
    for (int i0=0; i0<name.size(); ++i0) { makeOneFitTop(era,category,wpmin,wpmax,(std::string)name[i0],(std::string)conf::algo); } 
  }

  if (sample=="comb") {
    for (int i0=0; i0<name.size(); ++i0) { makeCombFit(era,category,wpmin,wpmax,(std::string)name[i0],"particlenet"); } }

}


void makeOneFit(std::string era, std::string category, std::string wpmin, std::string wpmax, std::string name, std::string name1) {

  std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe "+name1+"_sf/fitdir/datacard_"+name1+"_zqq_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt --PO categories=zqq,wqq";//,qcd
  std::string multidimfit       = "combine -M MultiDimFit -m 125 "+name1+"_sf/fitdir/datacard_"+name1+"_zqq_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root  --algo=singles --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string fitdiagnostics    = "combine -M FitDiagnostics -m 125 "+name1+"_sf/fitdir/datacard_"+name1+"_zqq_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root --saveShapes --saveWithUncertainties --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string mvmultidimfitfile = "mv higgsCombineTest.MultiDimFit.mH125.root "+name1+"_sf/fitdir/multidimfit_"+name1+"_zqq_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";
  std::string mvfitdiagnostics  = "mv fitDiagnostics.root "+name1+"_sf/fitdir/fitdiagnostics_"+name1+"_zqq_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";

  
  const char *command_txt2workspace     = txt2workspace.c_str();
  const char *command_multidimfit       = multidimfit.c_str();
  const char *command_fitdiagnostics    = fitdiagnostics.c_str();
  const char *command_mvmultidimfitfile = mvmultidimfitfile.c_str(); 
  const char *command_mvfitdiagnostics  = mvfitdiagnostics.c_str();
  
  system(command_txt2workspace);
  system(command_multidimfit);
  system(command_fitdiagnostics);
  system(command_mvmultidimfitfile);
  system(command_mvfitdiagnostics);  
 
}

void makeOneFitTop(std::string era, std::string category, std::string wpmin, std::string wpmax, std::string name, std::string name1) {
  
  std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe "+name1+"_sf/fitdir/datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt --PO categories=tp3,tp2,tp1,other"; //For top-tag SFs
  //std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe "+name1+"_sf/fitdir/datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt --PO categories=tp2,tp3,tp1,other"; //For W-tag SFs
  
  std::string multidimfit       = "combine -M MultiDimFit -m 125 "+name1+"_sf/fitdir/datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root  --algo=singles --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  
  std::string fitdiagnostics    = "combine -M FitDiagnostics -m 125 "+name1+"_sf/fitdir/datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root --saveShapes --saveWithUncertainties --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  
  std::string mvmultidimfitfile = "mv higgsCombineTest.MultiDimFit.mH125.root "+name1+"_sf/fitdir/multidimfit_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";
  std::string mvfitdiagnostics  = "mv fitDiagnosticsTest.root "+name1+"_sf/fitdir/fitdiagnostics_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  const char *command_txt2workspace     = txt2workspace.c_str();
  const char *command_multidimfit       = multidimfit.c_str();
  const char *command_fitdiagnostics    = fitdiagnostics.c_str();
  const char *command_mvmultidimfitfile = mvmultidimfitfile.c_str(); 
  const char *command_mvfitdiagnostics  = mvfitdiagnostics.c_str();

  system(command_txt2workspace);
  system(command_multidimfit);
  system(command_fitdiagnostics);
  system(command_mvmultidimfitfile);
  system(command_mvfitdiagnostics);
  
  // calculate impacts
  std::string impacts_1 = "combineTool.py -M Impacts -d "+name1+"_sf/fitdir/datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root -m 125 --doInitialFit --robustFit 1 --exclude 'rgx{prop.*}'";
  std::string impacts_2 = "combineTool.py -M Impacts -d "+name1+"_sf/fitdir/datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root -m 125 --robustFit 1 --doFits --parallel 60 --exclude 'rgx{prop.*}'";
  std::string impacts_3 = "combineTool.py -M Impacts -d "+name1+"_sf/fitdir/datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root -m 125 -o impacts.json --exclude 'rgx{prop.*}'";
  std::string impacts_4 = "plotImpacts.py -i impacts.json -o impacts";
  std::string impacts_5 = "mv impacts.pdf "+name1+"_sf/fitdir/impacts_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".pdf";
  std::string impacts_6 = "mv impacts.json "+name1+"_sf/fitdir/impacts_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".json";
  std::string impacts_7 = "mv combine_logger.out "+name1+"_sf/fitdir/combine_logger_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".out";
  std::string impacts_8 = "python ../CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/data/tutorials/longexercise/diffNuisances.py "+name1+"_sf/fitdir/fitdiagnostics_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root -A -g plots.root";
  std::string impacts_9 = "mv plots.root "+name1+"_sf/fitdir/impact_plots_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";
 
  const char *command_impacts1 = impacts_1.c_str();
  const char *command_impacts2 = impacts_2.c_str();
  const char *command_impacts3 = impacts_3.c_str();
  const char *command_impacts4 = impacts_4.c_str();
  const char *command_impacts5 = impacts_5.c_str();
  const char *command_impacts6 = impacts_6.c_str();
  const char *command_impacts7 = impacts_7.c_str();
  const char *command_impacts8 = impacts_8.c_str();
  const char *command_impacts9 = impacts_9.c_str();

  system(command_impacts1);
  system(command_impacts2);
  system(command_impacts3);
  system(command_impacts4);
  system(command_impacts5);
  system(command_impacts6);
  system(command_impacts7);
  system(command_impacts8);
  system(command_impacts9);
  system("rm higgsCombine*.root");
  
  // post fit plots
  system("root -l -q \'HeavyFlavourZCandleStudies.C(\"'"+(TString)era+"'\",\"tt1l\",\"'"+(TString)category+"'\",\"'"+(TString)wpmin+"'\",\"'"+(TString)wpmax+"'\",true,\"pass\",\"jmcor\")\'");
  system("root -l -q \'HeavyFlavourZCandleStudies.C(\"'"+(TString)era+"'\",\"tt1l\",\"'"+(TString)category+"'\",\"'"+(TString)wpmin+"'\",\"'"+(TString)wpmax+"'\",true,\"fail\",\"jmcor\")\'");
  
}


void makeCombFit(std::string era, std::string category, std::string wpmin, std::string wpmax, std::string name, std::string name1) {

  std::string combinecards = "combineCards.py datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_pt200to450.txt datacard_"+name1+"_zqq_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt datacard_"+name1+"_tt1l_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt > datacard_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt";
  std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe "+name1+"_sf/fitdir/datacard_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt --PO categories=zqq,wqq,tqq,qcd";//,qcd
  std::string multidimfit       = "combine -M MultiDimFit -m 125 "+name1+"_sf/fitdir/datacard_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt  --algo=singles --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string fitdiagnostics    = "combine -M FitDiagnostics -m 125 "+name1+"_sf/fitdir/datacard_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root --saveShapes --saveWithUncertainties --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string mvmultidimfitfile = "mv higgsCombineTest.MultiDimFit.mH125.root "+name1+"_sf/fitdir/multidimfit_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";
  std::string mvfitdiagnostics  = "mv fitDiagnostics.root "+name1+"_sf/fitdir/fitdiagnostics_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";

  
  const char *command_combinecards      = combinecards.c_str();
  const char *command_txt2workspace     = txt2workspace.c_str();
  const char *command_multidimfit       = multidimfit.c_str();
  const char *command_fitdiagnostics    = fitdiagnostics.c_str();
  const char *command_mvmultidimfitfile = mvmultidimfitfile.c_str(); 
  const char *command_mvfitdiagnostics  = mvfitdiagnostics.c_str();
  
  std::cout << combinecards << "\n";
  std::cout << txt2workspace << "\n";
  std::string mv = name1+"_sf/fitdir";
  const char *command_mv = mv.c_str();
  chdir(command_mv);
  system(command_combinecards);
  chdir("../..");
  
  std::string mvdc = "mv datacard_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt "+name1+"_sf/fitdir/datacard_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".txt" ;
  const char *command_mvdc = mvdc.c_str();
  //  system(command_mvdc);
  system(command_txt2workspace);
  system(command_multidimfit);
  system(command_fitdiagnostics);
  system(command_mvmultidimfitfile);
  system(command_mvfitdiagnostics);  
  
  // calculate impacts
  std::string impacts_1 = "combineTool.py -M Impacts -d "+name1+"_sf/fitdir/multidimfit_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root -m 125 --doInitialFit --robustFit 1";
  std::string impacts_2 = "combineTool.py -M Impacts -d "+name1+"_sf/fitdir/multidimfit_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root -m 125 --robustFit 1 --doFits --parallel 60";
  std::string impacts_3 = "combineTool.py -M Impacts -d "+name1+"_sf/fitdir/multidimfit_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root -m 125 -o impacts.json -h --exclude {prop*}";
  std::string impacts_4 = "plotImpacts.py -i impacts.json -o impacts";
  std::string impacts_5 = "mv impacts.pdf particlenet_sf/fitdir/impacts_"+name1+"_comb_"+category+"_"+wpmin+"to"+wpmax+"_"+era+"_"+name+".root";

  const char *command_impacts1 = impacts_1.c_str();
  const char *command_impacts2 = impacts_2.c_str();
  const char *command_impacts3 = impacts_3.c_str();
  const char *command_impacts4 = impacts_4.c_str();
  const char *command_impacts5 = impacts_5.c_str();

  system(command_impacts1);
  system(command_impacts2);
  system(command_impacts3);
  system(command_impacts4);
  system(command_impacts5);
  
  // post fit plots
  system("root -l -q \'HeavyFlavourZCandleStudies.C(\"'"+(TString)era+"'\",\"comb\",\"'"+(TString)category+"'\",\"'"+(TString)wpmin+"'\",\"'"+(TString)wpmax+"'\",true,\"pass\")\'");
  system("root -l -q \'HeavyFlavourZCandleStudies.C(\"'"+(TString)era+"'\",\"comb\",\"'"+(TString)category+"'\",\"'"+(TString)wpmin+"'\",\"'"+(TString)wpmax+"'\",true,\"fail\")\'");
}
