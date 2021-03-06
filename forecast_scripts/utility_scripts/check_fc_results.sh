for fctype in ice_only wavesice wavesice_ww3arctic
do
   for n in `seq 0 7`
   do

      fcdate=`date --date="today -${n}days" +%Y%m%d`
      dd=$TP4_REALTIME_RES/$fctype
      fil=$dd/$fcdate/final_output/*
      if [ -f $fil ]
      then
         echo " "
         echo "**************************************"
         echo "Latest TP4 $fctype forecast"
         echo $fcdate
         echo "**************************************"
         echo " "
         break
      fi
   done
done

for fctype in ice_only wavesice wavesice_ww3arctic
do
   for n in `seq 0 7`
   do

      fcdate=`date --date="today -${n}days" +%Y%m%d`
      dd=$BS1_REALTIME_RES/$fctype
      fil=$dd/$fcdate/final_output/*
      if [ -f $fil ]
      then
         echo " "
         echo "**************************************"
         echo "Latest BS1 $fctype forecast"
         echo $fcdate
         echo "**************************************"
         echo " "
         break
      fi
   done
done

for fctype in ice_only wavesice wavesice_ww3arctic
do
   for n in `seq 0 7`
   do

      fcdate=`date --date="today -${n}days" +%Y%m%d`
      dd=$FR1_REALTIME_RES/$fctype
      fil=$dd/$fcdate/final_output/*
      if [ -f $fil ]
      then
         echo " "
         echo "**************************************"
         echo "Latest FR1 $fctype forecast"
         echo $fcdate
         echo "**************************************"
         echo " "
         break
      fi
   done
done
