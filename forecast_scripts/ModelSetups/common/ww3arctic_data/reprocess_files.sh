# extracts only variables in $Vlist from ww3arctic files (originally they had everything)
ME=`readlink -f $0`
D0=`dirname $ME`

D=/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC/
Vlist="ice,hs,fp,dir"               # conc,main wave variables
Vlist="$Vlist,uuss,vuss,utus,vtus"  # Stokes drift parameters

for yr in 2015 2016
do
   cd $D/$yr/originals
   rm ../analysis_m1/* ../analysis_m2/*
   for Cdate in $yr*
   do
      cdate=${Cdate:0:8} # Cdate=yyyymmdd00

      echo "************************************"
      echo "$D0/ww3_arctic_sort.sh $cdate 2 0"
      $D0/ww3_arctic_sort.sh $cdate 2 0
      echo "************************************"
      # exit
   done
done
