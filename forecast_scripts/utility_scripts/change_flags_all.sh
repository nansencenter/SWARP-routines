me=`readlink -f $0`
here=`dirname $me`
fcm=$here/../ModelSetups
mod0=TP4a0.12

echo " "
echo "Print the following commands to compare flags for nested regions to TP4a0.12"
echo " "

for FCtype in ice_only wavesice_ww3arctic
do
   ss=$FCtype/inputs/flags
   f0=$fcm/$mod0/$ss
   # vim $f0
   for mod in FR1 BS1
   do
      f=$fcm/$mod/$ss
      echo "$here/compare_flags.sh $FCtype $mod"
   done
done
echo " "
