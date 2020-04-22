export OMP_NUM_THREADS=2

fov1=0.55
fov2=30
galaxy_type=ALG_YMW16

for phis in 1e3 1e4
do 
    fout1=simdat_${phis}
    fout2=simdat_${phis}_upper
    fin1=./simu/simdat_${phis}_${fov1}.txt
    fin2=./simu/simdat_${phis}_${fov2}.txt
    echo "mpiexec.hydra -n 120 ./nest_simu.py -f1 $fin1 -f2 $fin2 -o $fout1 -g ${galaxy_type}"
    mpiexec.hydra -n 120 ./nest_simu.py -f1 $fin1 -f2 $fin2 -o $fout1 -g ${galaxy_type}
    echo "mpiexec.hydra -n 120 ./nest_simu.py -upper 1 -f1 $fin1 -f2 $fin2 -o $fout2 -g ${galaxy_type}"
    mpiexec.hydra -n 120 ./nest_simu.py -upper 1 -f1 $fin1 -f2 $fin2 -o $fout2 -g ${galaxy_type}
done
exit  
