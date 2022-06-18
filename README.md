## Run the full chain for the SF extraction

./runFullChain.sh <object> [T|W] <year> [2016|2017|2018] <version> ["Nominal"|"MD"]

N.B.: Small tweaks are needed to be done inside the "configuration.h" and "makeFits.C" according to the object/version selected for the SF extraction.

## Setup CMSSW and combine



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
