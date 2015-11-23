# script to restore forecast expt directories
# if it has been cleaned

if [ $# -eq 0 ]
then
   echo Usage:
   echo restore_expt.sh Reg
   echo "where Reg = BS1 or FR1"
   exit
else
   Reg=$1
fi
Xno=0
X=01.$Xno
E=01$Xno

if [ $Reg == "BS1" ]
then
   rtdir=$BS1_REALTIME
   R=BS1a0.045
   nproc=112
elif [ $Reg == "FR1" ]
then
   rtdir=$FR1_REALTIME
   R=FR1a0.03
   nproc=51
fi
THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice_$Reg

echo " "
cd $rtdir/expt_$X
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
echo Changing REGION.src
file=REGION.src
cat ../bak/$file | sed \
-e "s/NATa1.00/$R/g" \
> ../$file


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
echo " "
