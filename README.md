## Create 2D histograms pT(jet) vs m(jet) from the trees
 [this is to speed up the process]:

`root [0] .L make2DTemplates.C` 

`root [1] mainfunction("tt1L”) or mainfunction(“zqq”)`

## Create the 1D templates:

`root [0] .L HeavyFlavourZCandleStudies.C` 

`root [1] mainfunction("tt1L”) or mainfunction(“zqq”)`

## Create the datacards:
The Datacards will be produced under particlenet_sf/fitdir 
[if one keeps the default naming/settings]


`root [0] .L makeDatacards.C` 

`root [1] makeDatacards("tt1L”) or makeDatacards(“zqq”)`


## Do the fit:

`root [0] .L makeFits.C `

`root [1] makeFits("comb”)`

 or “zqq” or “tt1L” to fit only specific samples