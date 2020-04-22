export OMP_NUM_THREADS=2
#fgtype=ALG_YMW16
for fgtype in ALG_NE2001 ALG_YMW16
do
    echo "mpiexec.hydra -n 240 ./nest_samp.py -fc frb_cat.txt -fs tel_svy.txt -g ${fgtype} -upper 1 -o ${fgtype}_upper"
    mpiexec.hydra -n 240 ./nest_samp.py -fc frb_cat.txt -fs tel_svy.txt -g ${fgtype} -upper 1 -o ${fgtype}_upper 
    #echo "mpiexec.hydra -n 120 ./nest_samp.py -fc frb_cat.txt -fs tel_svy.txt -g ${fgtype} -upper 1 halo 1 -o ${fgtype}_upper_halo"
    #mpiexec.hydra -n 120 ./nest_samp.py -fc frb_cat.txt -fs tel_svy.txt -g ${fgtype} -upper 1 -halo 1 -o ${fgtype}_upper_halo
done
