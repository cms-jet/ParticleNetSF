Installing Combine (lpc)
------------------------
link: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/

voms-proxy-init --voms cms

setenv SCRAM_ARCH slc7_amd64_gcc700   (export SCRAM_ARCH=slc7_amd64_gcc700)

cmsrel CMSSW_10_2_13

'cd CMSSW_10_2_13/src'

'cmsenv'

'git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit'

'git clone https://github.com/cms-analysis/CombineHarvester.git'

'cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit'
'git fetch origin'
'git checkout v8.2.0'
'scramv1 b clean; scramv1 b' 

Clone specific branch from gitlab (top/W tagging framework):
-----------------------------------------------------------
'cd ../../../../'

'git clone -b ParticleNet_TopW_SFs_NanoV9 ssh://git@gitlab.cern.ch:7999/lpaizano/particlenet-top-w-calibration.git'

'cd particlenet-top-w-calibration/'

Run the full chain for SF extraction:
------------------------------------

./runFullChain.sh  [T|W]  [2016|2017|2018]  ["Nominal"|"MD"]

N.B.: Small tweaks are needed to be done inside the "configuration.h" and "makeFits.C" according to the object/version selected for the SF extraction.

gitlab link: https://gitlab.cern.ch/lpaizano/particlenet-top-w-calibration/-/tree/ParticleNet_TopW_SFs_NanoV9/

## Create 2D histograms pT(jet) vs m(jet) from the trees
 [this is to speed up the process]:

`root [0] .L make2DTemplates.C` 

`root [1] mainfunction("tt1l","2017","0.90","1.”)`

## Create the 1D templates:

`root [0] .L HeavyFlavourZCandleStudies.C` 

`root [1] HeavyFlavourZCandleStudies("2017","tt1l","bb","0.90","1.",false,"pass")`

## Create the datacards:
The Datacards will be produced under particlenet_sf/fitdir 
[if one keeps the default naming/settings]


`root [0] .L makeDatacards.C` 

`root [1] makeDatacards("2017","tt1l","bb","0.90","1.")`


## Do the fit:

`root [0] .L makeFits.C `

`root [1] makeFits("2017","bb","0.90","1.","tt1l")`

 or “zqq” or “tt1L” to fit only specific samples
