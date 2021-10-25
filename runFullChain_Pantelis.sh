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
    #mist_rates=("1p0" "0p5" "0p1") 
    mist_rates=("1p0")
    if [ ${year} == 2018 ];
    then
      if [ ${version} == "Nominal" ]; 
      then 
        WPs_FullVer_vs_QCD=(["1p0"]="0.59" ["0p5"]="0.81" ["0p1"]="0.97") #0: 1p0, 1: 0p5, 2: 0p1
      fi
    fi
elif [ ${object} == "W" ];
then   
    if [ ${year} == 2018 ];
    then
      if [ ${version} == "Nominal" ];
      then
        #mist_rates=("5p0" "1p0" "0p5")
        mist_rates=("5p0")
        WPs_FullVer_vs_QCD=(["5p0"]="0.72" ["1p0"]="0.95" ["0p5"]="0.98") #0: 5p0, 1: 1p0, 2: 0p5
      elif [ ${version} == "MD" ];     
      then
        mist_rates=("2p5" "1p0" "0p5")
        WPs_FullVer_vs_QCD=(["2p5"]="0.63" ["1p0"]="0.84" ["0p5"]="0.91") #0: 2p5, 1: 1p0, 2: 0p5
      fi 
    fi
fi 

#wpmin="0.90"
#wpmin="0.89"

#for era in 2016 2017 2018
for era in ${year}
do
   for mistRate in "${mist_rates[@]}";
   do
      unset wpmin

      wpmin="${WPs_FullVer_vs_QCD[${mistRate}]}"
          
   
      cmd_templates2d=$(echo 'make2DTemplates.C("tt1l","'${era}'","'${wpmin}'","1.00")')
      cmd_templates1d=$(echo 'HeavyFlavourZCandleStudies.C("'${era}'","tt1l","top","'${wpmin}'","1.00",false,"pass")')
      cmd_datacards=$(echo 'makeDatacards.C("'${era}'","tt1l","top","'${wpmin}'","1.00")')
      cmd_makefits=$(echo 'makeFits.C("'${era}'","top","'${wpmin}'","1.00","tt1l")')

      root -l -q ${cmd_templates2d}
      root -l -q ${cmd_templates1d}
      root -l -q ${cmd_datacards}
      root -l -q ${cmd_makefits}
   done
done
