yr=2015
Edir=/work/shared/nersc/ECMWFR_T799/
echo " "

for vbl in D2M MSL T2M TCC U10M V10M
do
   # file to check
   ncfil=$Edir/ec_atmo_geo_la_${vbl}_$yr.nc
   echo Checking $ncfil ...

   # get number of recs
   ss=`ncdump -h $ncfil |grep UNLIMITED `
   ss=($ss)
   nrec=${ss[5]}
   nrec=${nrec#(}
   echo "Number of records          : $nrec"

   # desired number
   Nfc=$((8*4+3))
   jday=`date +%j`
   nrec0=$(($Nfc+4*$jday))
   echo "Expected number of records : $nrec0"
done
