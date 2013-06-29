#! /bin/bash
whichTTBar=$1
#for whichTTBar in 0 1 2 3 4 5 6 7 8
#for whichTTBar in 1 2 3 4 5 6 7 8
#do
 root -l -b -q "plotHaddersStopHack.C(1, ${whichTTBar}, 1)"
 root -l -b -q "plotHaddersStopHack.C(0, ${whichTTBar}, 1)"
 root -l -b -q "plotHaddersStopHack.C(1, ${whichTTBar}, 0)"
 root -l -b -q "plotHaddersStopHack.C(0, ${whichTTBar}, 0)"
#done
#hadd -f TTBar
#  LocalWords:  hadd
