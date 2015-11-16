# script to restore forecast expt directories
# if it has been cleaned

if [ $# -eq 0 ]
then
   echo Usage:
   echo restore_expt.sh Xno
   echo "where Xno = 1 (ice_only), 2 (waves-in-ice - WAM), 3 (waves-in-ice - WW3)"
   exit
else
   X=01.$1
   E=01$1
fi

echo " "
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
echo " "
