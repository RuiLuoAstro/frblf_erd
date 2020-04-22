#!/bin/sh

for phis in 1e3 1e4
do 
    o1=./plots/simu/simdat_${phis}.eps
    o2=./plots/simu/simdat_${phis}_upper.eps
    in1=./nest_out/simu/simdat_${phis}
    in2=./nest_out/simu/simdat_${phis}_upper
    ./pltpost.py -f $in1 -o $o1 -title "Mock data" -up 0 -bo 1
    ./pltpost.py -f $in2 -o $o2 -title "Mock data" -up 1 -bo 1
done
exit  
