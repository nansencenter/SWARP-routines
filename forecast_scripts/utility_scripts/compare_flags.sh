me=`readlink -f $0`
here=`dirname $me`
fcm=$here/../ModelSetups
mod0=TP4a0.12

if [ $# -ne 2 ]
then
   echo "Usage: $me [forecast type] [region]"
   echo "*forecast type: 'ice_only' or 'wavesice_ww3arctic'"
   echo "*region: 'FR1' or 'BS1'"
   exit
else
   FCtype=$1
   mod1=$2
fi

if [ $mod1 == "FR1" ]
then
   mod=${mod1}a0.03
elif [ $mod1 == "BS1" ]
then
   mod=${mod1}a0.045
fi

ss=$FCtype/inputs/flags
f0=$fcm/$mod0/$ss

f=$fcm/$mod/$ss
gvim -d $f $f0
