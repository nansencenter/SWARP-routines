me=`readlink -f $0`
here=`dirname $me`
fcm=$here/../ModelSetups
mod0=TP4a0.12

if [ $# -ne 1 ]
then
   echo "Usage: $me [forecast type]"
   echo "*forecast type: 'ice_only' or 'wavesice_ww3arctic'"
   echo "*region: 'FR1' or 'BS1'"
   echo "Compares flags for FR1 to the TP4 flags"
   echo "Uses same flags for BS1 as FR1"
   exit
else
   FCtype=$1
fi

mod=FR1a0.03
mod2=BS1a0.045

ss=$FCtype/inputs/flags
f0=$fcm/$mod0/$ss

f=$fcm/$mod/$ss
gvim -d $f $f0

f2=$fcm/$mod2/$ss
echo cp $f $f2
cp $f $f2
gvim -d $f2 $f0
