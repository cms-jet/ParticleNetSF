#!/bin/bash 

wpmin="0p05"

for era in 2016 2017 2018
do
    cmd_templates2d=$(echo 'make2DTemplates.C("tt1l","'${era}'","'${wpmin}'","1.")')
    cmd_templates1d=$(echo 'HeavyFlavourZCandleStudies.C("'${era}'","tt1l","w","'${wpmin}'","1.",false,"pass")')
    cmd_datacards=$(echo 'makeDatacards.C("'${era}'","tt1l","w","'${wpmin}'","1.")')
    cmd_makefits=$(echo 'makeFits.C("'${era}'","w","'${wpmin}'","1.","tt1l")')

    #root -l -q ${cmd_templates2d}
    #root -l -q ${cmd_templates1d}
    root -l -q ${cmd_datacards}
    root -l -q ${cmd_makefits}
done
