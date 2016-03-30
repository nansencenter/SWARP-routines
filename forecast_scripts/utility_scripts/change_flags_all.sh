me=`readlink -f $0`
here=`dirname $me`
fcm=$here/../ModelSetups
mod0=TP4a0.12

for FCtype in ice_only wavesice_ww3arctic
do
   ss=$FCtype/inputs/flags
   f0=$fcm/$mod0/$ss
   vim $f0
   for mod in FR1a0.03 BS1a0.045
   do
      f=$fcm/$mod/$ss
      vim -d $f0 $f
   done
done
