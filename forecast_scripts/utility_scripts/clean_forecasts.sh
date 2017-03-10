# clean older than 2 weeks
# - contents will have been done anyway probably
cday=`date --date="-14days" +%Y%m%d`

FCdirs=($TP4_REALTIME_RES/ice_only $TP4_REALTIME_RES/wavesice_ww3arctic $FR1_REALTIME_RES/ice_only $FR1_REALTIME_RES/wavesice_ww3arctic \
      $BS1_REALTIME_RES/ice_only $BS1_REALTIME_RES/wavesice_ww3arctic)

for FCdir in ${FCdirs[@]}
do

   cd $FCdir
   # pwd

   for dday in *
   do
      if [ ${#dday} -ne 8 ]
      then
         # leave unusual folders
         # - $dday should be a date
         continue
      fi
      if [ $dday -lt $cday ]
      then
         rm -r $dday
      fi
   done

done
