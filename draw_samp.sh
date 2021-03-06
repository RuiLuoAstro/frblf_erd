#for fgtype in ETG_NE2001_upper ETG_YMW16_upper LTG_NE2001_upper LTG_YMW16_upper ALG_NE2001_upper ALG_YMW16_upper
#for fgtype in ETG_NE2001_upper_halo ETG_YMW16_upper_halo LTG_NE2001_upper_halo LTG_YMW16_upper_halo ALG_NE2001_upper_halo ALG_YMW16_upper_halo

for fgtype in ALG_NE2001_upper ALG_YMW16_upper
do
    echo "New plot of ${fgtype} has been finished."
    ./pltpost.py -f ./nest_out/previous/${fgtype} -o ./plots/previous/${fgtype}.eps -title "Real Sample" -up 1 -bo 1
done
