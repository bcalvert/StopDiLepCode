#! /bin/bash
rm -rf PlotReadMe.txt
../makeRoot.sh GenerateReadMe.C
./GenerateReadMe -b -q >> PlotReadMe.txt

