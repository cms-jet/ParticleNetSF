#include <bits/stdc++.h>
#include <string>
#include <unistd.h>

void makeOneFit(std::string name);
void makeOneFitTop(std::string name);
void makeCombFit(std::string name);

void makeFits(TString sample) {
  
  std::vector<std::string> name;
  // from ggH
  //  name.push_back("pt450to500"); name.push_back("pt500to550"); name.push_back("pt550to600"); name.push_back("pt600to675"); name.push_back("pt675to800"); name.push_back("pt800to1200");
  name.push_back("pt450to600"); name.push_back("pt600to800"); name.push_back("pt800to1200");

  if (sample=="zqq") {
    for (int i0=0; i0<name.size(); ++i0) { makeOneFit(name[i0]); } }
  
  if (sample=="tt1L" || sample=="ttbar1L") {
    for (int i0=0; i0<name.size(); ++i0) { makeOneFitTop(name[i0]); } }

  if (sample=="comb") {
    for (int i0=0; i0<name.size(); ++i0) { makeCombFit(name[i0]); } }

}


void makeOneFit(std::string name) {
   
  std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe ./fitdir/datacard_"+name+".txt --PO categories=zqq,qcd,tqq";
  std::string multidimfit       = "combine -M MultiDimFit -m 125 ./fitdir/datacard_"+name+".root  --algo=singles --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string fitdiagnostics    = "combine -M FitDiagnostics -m 125 ./fitdir/datacard_"+name+".root --saveShapes --saveWithUncertainties --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string mvmultidimfitfile = "mv higgsCombineTest.MultiDimFit.mH125.root ./fitdir/multidimfit_"+name+".root";
  std::string mvfitdiagnostics  = "mv fitDiagnostics.root ./fitdir/fitdiagnostics_"+name+".root";
  
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
  std::string impacts_1 = "combineTool.py -M Impacts -d ./fitdir/datacard_"+name+".root -m 125 --doInitialFit --robustFit 1";
  std::string impacts_2 = "combineTool.py -M Impacts -d ./fitdir/datacard_"+name+".root -m 125 --robustFit 1 --doFits --parallel 60";
  std::string impacts_3 = "combineTool.py -M Impacts -d ./fitdir/datacard_"+name+".root -m 125 -o impacts.json";
  std::string impacts_4 = "plotImpacts.py -i impacts.json -o impacts";
  std::string impacts_5 = "mv impacts.pdf ./fitdir/impacts_"+name+".pdf";

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
  
}


void makeOneFitTop(std::string name) {
   
  std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe ./fitdir/datacard_tt1l_"+name+".txt --PO categories=tqq,wqq,tp1,other";
  std::string multidimfit       = "combine -M MultiDimFit -m 125 ./fitdir/datacard_tt1l_"+name+".root  --algo=singles --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string fitdiagnostics    = "combine -M FitDiagnostics -m 125 ./fitdir/datacard_tt1l_"+name+".root --saveShapes --saveWithUncertainties --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string mvmultidimfitfile = "mv higgsCombineTest.MultiDimFit.mH125.root ./fitdir/multidimfit_tt1l_"+name+".root";
  std::string mvfitdiagnostics  = "mv fitDiagnostics.root ./fitdir/fitdiagnostics_tt1l_"+name+".root";
  
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
  std::string impacts_1 = "combineTool.py -M Impacts -d ./fitdir/datacard_tt1l_"+name+".root -m 125 --doInitialFit --robustFit 1";
  std::string impacts_2 = "combineTool.py -M Impacts -d ./fitdir/datacard_tt1l_"+name+".root -m 125 --robustFit 1 --doFits --parallel 60";
  std::string impacts_3 = "combineTool.py -M Impacts -d ./fitdir/datacard_tt1l_"+name+".root -m 125 -o impacts.json";
  std::string impacts_4 = "plotImpacts.py -i impacts.json -o impacts";
  std::string impacts_5 = "mv impacts.pdf ./fitdir/impacts_tt1l_"+name+".pdf";

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
  
}


void makeCombFit(std::string name) {
  
  std::string combinecards      = "combineCards.py datacard_"+name+".txt datacard_tt1l_"+name+".txt datacard_tt1l_pt200to400.txt > datacard_comb_"+name+".txt";
  //std::string combinecards      = "combineCards.py datacard_"+name+".txt datacard_tt1l_"+name+".txt > datacard_comb_"+name+".txt";
  //std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe ./fitdir/datacard_comb_"+name+".txt --PO categories=zqq,wqq,qcd,tqq,tp1,other";
  std::string txt2workspace     = "text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe ./fitdir/datacard_comb_"+name+".txt --PO categories=zqq,wqq,qcd,tqq,tp1";
  std::string multidimfit       = "combine -M MultiDimFit -m 125 ./fitdir/datacard_comb_"+name+".root  --algo=singles --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string fitdiagnostics    = "combine -M FitDiagnostics -m 125 ./fitdir/datacard_comb_"+name+".root --saveShapes --saveWithUncertainties --robustFit=1 --cminDefaultMinimizerTolerance 5.";
  std::string mvmultidimfitfile = "mv higgsCombineTest.MultiDimFit.mH125.root ./fitdir/multidimfit_comb_"+name+".root";
  std::string mvfitdiagnostics  = "mv fitDiagnostics.root ./fitdir/fitdiagnostics_comb_"+name+".root";
  
  const char *command_combinecards      = combinecards.c_str();
  const char *command_txt2workspace     = txt2workspace.c_str();
  const char *command_multidimfit       = multidimfit.c_str();
  const char *command_fitdiagnostics    = fitdiagnostics.c_str();
  const char *command_mvmultidimfitfile = mvmultidimfitfile.c_str(); 
  const char *command_mvfitdiagnostics  = mvfitdiagnostics.c_str();
  
  chdir("fitdir");
  system(command_combinecards);
  chdir("..");
  system(command_txt2workspace);
  system(command_multidimfit);
  system(command_fitdiagnostics);
  system(command_mvmultidimfitfile);
  system(command_mvfitdiagnostics);  

  
  // calculate impacts
  std::string impacts_1 = "combineTool.py -M Impacts -d ./fitdir/datacard_comb_"+name+".root -m 125 --doInitialFit --robustFit 1";
  std::string impacts_2 = "combineTool.py -M Impacts -d ./fitdir/datacard_comb_"+name+".root -m 125 --robustFit 1 --doFits --parallel 60";
  std::string impacts_3 = "combineTool.py -M Impacts -d ./fitdir/datacard_comb_"+name+".root -m 125 -o impacts.json -h --exclude {prop*}";
  std::string impacts_4 = "plotImpacts.py -i impacts.json -o impacts";
  std::string impacts_5 = "mv impacts.pdf ./fitdir/impacts_comb_"+name+".pdf";

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
  
}
