# script to restore forecast expt directories
# if it has been cleaned

if [ $# -lt 2 ]
then
   echo Usage:
   echo restore_expt.sh Region Xno
   echo "where Region = TP4, FR1 or BS1"
   echo "where Xno = 1 (ice_only), 2 (waves-in-ice - WAM), 3 (waves-in-ice - WW3)"
   exit
else
   HYCOMreg=$1
   FCtype=$2
fi


echo " "
if [ $HYCOMreg == "TP4" ]
then
   R=TP4a0.12
elif [ $HYCOMreg == "FR1" ]
then
   R=FR1a0.03
elif [ $HYCOMreg == "BS1" ]
then
   R=BS1a0.045
else
   echo "Unknown region $HYCOMreg"
   exit
fi

if [ $FCtype == "ice_only" ]
then
   if [ $HYCOMreg == "TP4" ]
   then
      Xno=1
   else
      Xno=0
   fi
elif [ $FCtype == "wavesice_ww3arctic" ]
then
   if [ $HYCOMreg == "TP4" ]
   then
      Xno=3
   else
      Xno=1
   fi
elif [ $FCtype == "wavesice" ]
then
   if [ $HYCOMreg == "TP4" ]
   then
      Xno=2
   else
      echo "$FCtype is not set up for $R region"
      exit
   fi
else
   echo "Unknown forecast type $FCtype"
   exit
fi

echo " "
X=01.$Xno
E=01$Xno
cd $TP4_REALTIME/expt_$X
pwd

for f in bak/*
do
   if [ -f $f ]
   then
      echo cp $f .
      cp $f .
   fi
done


echo " "
echo Changing EXPT.src...
file=EXPT.src
cat bak/$file | sed \
-e "s/01.0/$X/g" \
-e "s/010/$E/g" \
> $file
echo " "


echo Changing blkdat.input...
file=blkdat.input
cat bak/$file | sed \
-e "s/010/$E/g" \
> $file
