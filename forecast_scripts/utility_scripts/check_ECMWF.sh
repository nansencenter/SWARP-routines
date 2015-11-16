Edir=/work/shared/nersc/ECMWFR_T799/
source $SWARP_ROUTINES/source_files/hex_vars.src
address=$FORECAST/fc_alert_email.txt
email=$(cat $address)

if [ $# -eq 0 ]
then
   alert=0
else
   alert=$1
fi

yr=`date +%Y`
correct=0
echo " " > txt

for vbl in D2M MSL T2M TCC U10M V10M
do
   # file to check
   ncfil=$Edir/ec_atmo_geo_la_${vbl}_$yr.nc
   echo Checking $ncfil ...                     >> txt

   # get number of recs
   ss=`$ncdump -h $ncfil |grep UNLIMITED `
   ss=($ss)
   nrec=${ss[5]}
   nrec=${nrec#(}
   echo "Number of records          : $nrec"    >> txt

   # desired number
   Nfc=$((8*4+3))
   jday=`date +%j`
   nrec0=$(($Nfc+4*$jday))
   echo "Expected number of records : $nrec0"   >> txt

   if [ $nrec -eq $nrec0 ]
   then
      correct=$((correct+1))
   fi
done

# EMAIL 
if [ $alert -eq 1 ]
then


   if [ $correct -lt 6 ]
   then
      mail -s "Missing records in ECMWF forcing" $email < txt
   fi

elif [ $alert -eq 2 ]
then

   if [ $correct -eq 6 ]
   then
      mail -s "ECMWF forcing OK now" $email < txt
   fi

else
   cat txt
   echo " "
fi
rm txt
