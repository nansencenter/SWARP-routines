# extracts only variables in $Vlist from ww3arctic files (originally they had everything)

D=/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC/
cd $D
Vlist="ice,hs,fp,dir"

for yr in 2015 2016
do
   for ddir in forecast analysis_m1 analysis_m2
   do
      cd $D/$yr/$ddir
      echo "************************************"
      pwd
      for f in SWARP*.nc
      do
         echo ncks -v $Vlist $f tmp.nc
         ncks -v $Vlist $f tmp.nc
         mv tmp.nc $f
      done
      echo "************************************"
   done
done
