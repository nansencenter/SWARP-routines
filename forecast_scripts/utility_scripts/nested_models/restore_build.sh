# script to restore forecast Build directories
# if it has been cleaned

if [ $# -eq 0 ]
then
   echo Usage:
   echo restore_build.sh reg
   echo "where reg = BS1 or FR1"
   exit
else
   reg=$1
fi
Xno=0
X=01.$Xno

if [ $reg == "BS1" ]
then
   rtdir=$BS1_REALTIME
   R=BS1a0.045
   nproc=112
elif [ $reg == "FR1" ]
then
   rtdir=$FR1_REALTIME
   R=FR1a0.03
   nproc=51
fi
THISFC=$SWARP_ROUTINES/forecast_scripts/nested_models/wavesice_$reg

echo " "
cd $rtdir/Build_V2.2.12_X$X
pwd

echo " "
echo Changing REGION.src
file=REGION.src
cat ../bak/$file | sed \
-e "s/NATa1.00/$R/g" \
> ../$file


if [ ! -f ../expt_$X/blkdat.input ]
then
   # need some files from expt directory
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/nested_models/restore_expt.sh $reg
fi

echo cp svn_Build/* .
cp svn_Build/* .

echo ./setup_patch.sh $nproc
./setuppatch.sh $nproc


echo "cp $THISFC/inputs/flags ."
echo " "
cp $THISFC/inputs/flags .
