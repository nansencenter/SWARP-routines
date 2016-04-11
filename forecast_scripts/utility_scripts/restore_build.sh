# script to restore forecast Build directories
# if it has been cleaned
source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -lt 2 ]
then
   echo Usage:
   echo restore_build.sh Region FCtype
   echo "where Region = TP4, FR1 or BS1"
   echo "where FCtype = 'ice_only', 'wavesice' (WAM), or 'wavesice_ww3arctic' (WW3 Arctic grid)"
   exit
else
   HYCOMreg=$1
   FCtype=$2
fi

echo " "
if [ $HYCOMreg == "TP4" ]
then
   R=TP4a0.12
   Nproc=133
elif [ $HYCOMreg == "FR1" ]
then
   R=FR1a0.03
   Nproc=51
elif [ $HYCOMreg == "BS1" ]
then
   R=BS1a0.045
   Nproc=112
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

X=01.$Xno
cd $RTmods/$R/Build_V2.2.12_X$X
pwd

# ls; exit

echo " "
echo Changing REGION.src
file=REGION.src
cat ../bak/$file | sed \
-e "s/NATa1.00/$R/g" \
-e "s/ECMWFR_T799/ECMWFR_T799_bak/g" \
> ../$file


if [ ! -f ../expt_$X/blkdat.input ]
then
   # need some files from expt directory
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_expt.sh $HYCOMreg $FCtype
fi

echo cp svn_Build/* .
cp svn_Build/* .

echo ./setup_patch.sh $Nproc
./setuppatch.sh $Nproc


THISFC=$SWARP_ROUTINES/forecast_scripts/ModelSetups/$R/$FCtype

echo "cp $THISFC/inputs/flags ."
echo " "
cp $THISFC/inputs/flags .

echo " "

if [ $Xno -eq 1 ] && [ $HYCOMreg == "TP4" ]
then
   # set up nesting
   cd ..
   ndir=nest_nersc/011/outer
   nsh=nest_nersc/bin/nest_outer.sh
   rm -f $ndir/SCRATCH/*
   $nsh $X $BS1_REALTIME 01.0 # BS1
   rm -f $ndir/SCRATCH/*
   $nsh $X $FR1_REALTIME 01.0 # FR1
fi
