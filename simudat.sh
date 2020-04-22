#!/bin/sh

# FRB number
Nfrb=100

# LF parameters
alpha=-1.5
logls=45
logl0=40.0
dnu=1000

# Width parameters
mu=0.4
sigma=0.3

# Selection criteria
np=2
sn0=10

# Two surveys
g1=0.7
bw1=300
Ts1=30
fov1=0.55

g2=0.05
bw2=300
Ts2=100
fov2=30

# Host galaxy categories
galaxy_type=ALG_YMW16

for phis in 1e2 1e3 1e4
#for phis in 1e2
do 
    outputfile1=./simu/simdat_${phis}_${fov1}.txt
    outputfile2=./simu/simdat_${phis}_${fov2}.txt
    echo "python simufrb.py -ns $Nfrb -phis $phis -alpha $alpha -logls $logls -logl0 $logl0 -dnu $dnu -fgt $galaxy_type -mu ${mu} -sig ${sigma} -ga $g1 -np $np -bw $bw1 -ts $Ts1 -sn0 $sn0 -fov $fov1 -out $outputfile1"
    python simufrb.py -ns $Nfrb -phis $phis -alpha $alpha -logls $logls -logl0 $logl0 -dnu $dnu -fgt $galaxy_type -mu ${mu} -sig ${sigma} -ga $g1 -np $np -bw $bw1 -ts $Ts1 -sn0 $sn0 -fov $fov1 -out $outputfile1 &
    echo "python simufrb.py -ns $Nfrb -phis $phis -alpha $alpha -logls $logls -logl0 $logl0 -dnu $dnu -fgt $galaxy_type -mu ${mu} -sig ${sigma} -ga $g2 -np $np -bw $bw2 -ts $Ts2 -sn0 $sn0 -fov $fov2 -out $outputfile2"
    python simufrb.py -ns $Nfrb -phis $phis -alpha $alpha -logls $logls -logl0 $logl0 -dnu $dnu -fgt $galaxy_type -mu ${mu} -sig ${sigma} -ga $g2 -np $np -bw $bw2 -ts $Ts2 -sn0 $sn0 -fov $fov2 -out $outputfile2 &
done
exit
