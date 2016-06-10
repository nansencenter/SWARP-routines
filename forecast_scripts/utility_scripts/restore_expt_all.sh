me=`readlink -f $0`  # full path of this script
here=`dirname $me`   # directory of this script
for reg in TP4 #FR1 BS1
do
   for FCtype in wavesice_ww3arctic #ice_only
   do
      $here/restore_expt.sh $reg $FCtype
   done
done
