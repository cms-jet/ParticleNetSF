#!/bin/bash 

object=$1
year=$2
version=$3

declare -a WPs_TopvsQCD
declare -a WPs_WvsQCD
declare -a WPs_Wvs_QCD_MD
declare -A WPs_FullVer_vs_QCD
declare -A WP_MDVer_vs_QCD
declare -a mist_rates

if [ ${object} == "T" ];
then
    category="top"
elif [ ${object} == "W" ];
then
    category="w"
fi

if [ ${object} == "T" ];
then
    mist_rates=("1p0" "0p5" "0p1") 
    #mist_rates=("1p0")
    if [ ${year} == 2018 ];
    then
      if [ ${version} == "Nominal" ]; 
      then 
        WPs_FullVer_vs_QCD=(["1p0"]="0.58" ["0p5"]="0.80" ["0p1"]="0.97") #0: 1p0, 1: 0p5, 2: 0p1
      fi
    elif [ ${year} == 2017 ];
    then
      if [ ${version} == "Nominal" ];
      then
        WPs_FullVer_vs_QCD=(["1p0"]="0.581" ["0p5"]="0.801" ["0p1"]="0.970") #0: 1p0, 1: 0p5, 2: 0p1
      fi
    elif [ ${year} == 2016 ];
    then
      if [ ${version} == "Nominal" ];
      then
	WPs_FullVer_vs_QCD=(["1p0"]="0.490" ["0p5"]="0.738" ["0p1"]="0.957") #0: 1p0, 1: 0p5, 2: 0p1 #2015
        #WPs_FullVer_vs_QCD=(["1p0"]="0.495" ["0p5"]="0.733" ["0p1"]="0.958") #0: 1p0, 1: 0p5, 2: 0p1 #2016
      fi    
    fi
elif [ ${object} == "W" ];
then   
    if [ ${year} == 2018 ];
    then
      if [ ${version} == "Nominal" ];
      then
        mist_rates=("5p0" "1p0" "0p5")
        WPs_FullVer_vs_QCD=(["5p0"]="0.70" ["1p0"]="0.94" ["0p5"]="0.98") #0: 5p0, 1: 1p0, 2: 0p5
      elif [ ${version} == "MD" ];     
      then
        mist_rates=("2p5" "1p0" "0p5")
        WPs_FullVer_vs_QCD=(["2p5"]="0.59" ["1p0"]="0.82" ["0p5"]="0.90") #0: 2p5, 1: 1p0, 2: 0p5
      fi 
    elif [ ${year} == 2017 ];
    then
      if [ ${version} == "Nominal" ];
      then
        mist_rates=("5p0" "1p0" "0p5")
        WPs_FullVer_vs_QCD=(["5p0"]="0.709" ["1p0"]="0.944" ["0p5"]="0.978") #0: 5p0, 1: 1p0, 2: 0p5
      elif [ ${version} == "MD" ];
      then
        mist_rates=("2p5" "1p0" "0p5")
        WPs_FullVer_vs_QCD=(["2p5"]="0.579" ["1p0"]="0.810" ["0p5"]="0.891") #0: 2p5, 1: 1p0, 2: 0p5
      fi
    elif [ ${year} == 2016 ];
    then
      if [ ${version} == "Nominal" ];
      then
        mist_rates=("5p0" "1p0" "0p5")
	WPs_FullVer_vs_QCD=(["5p0"]="0.677" ["1p0"]="0.935" ["0p5"]="0.974") #0: 5p0, 1: 1p0, 2: 0p5 #2015
        #WPs_FullVer_vs_QCD=(["5p0"]="0.668" ["1p0"]="0.934" ["0p5"]="0.974") #0: 5p0, 1: 1p0, 2: 0p5 #2016
      elif [ ${version} == "MD" ];
      then
        mist_rates=("2p5" "1p0" "0p5")
	WPs_FullVer_vs_QCD=(["2p5"]="0.637" ["1p0"]="0.845" ["0p5"]="0.910") #0: 2p5, 1: 1p0, 2: 0p5 #2015
        #WPs_FullVer_vs_QCD=(["2p5"]="0.642" ["1p0"]="0.842" ["0p5"]="0.907") #0: 2p5, 1: 1p0, 2: 0p5 #2016
      fi 
    fi
fi 

for era in ${year}
do
   for mistRate in "${mist_rates[@]}";
   do
      unset wpmin

      wpmin="${WPs_FullVer_vs_QCD[${mistRate}]}"
          
   
      cmd_templates2d=$(echo 'make2DTemplates.C("tt1l","'${era}'","'${wpmin}'","1.00")')
      cmd_templates1d=$(echo 'HeavyFlavourZCandleStudies.C("'${era}'","tt1l","'${category}'","'${wpmin}'","1.00",false,"pass")')
      cmd_datacards=$(echo 'makeDatacards.C("'${era}'","tt1l","'${category}'","'${wpmin}'","1.00")')
      cmd_makefits=$(echo 'makeFits.C("'${era}'","'${category}'","'${wpmin}'","1.00","tt1l")')
      #here
      #root -l -q ${cmd_templates2d}
      root -l -q ${cmd_templates1d}
      root -l -q ${cmd_datacards}
      root -l -q ${cmd_makefits} | tee ${object}_${year}_${version}_${wpmin}.txt
   done
done
