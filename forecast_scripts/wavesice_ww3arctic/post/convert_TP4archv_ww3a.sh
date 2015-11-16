#!/bin/bash
# Runs hyc2proj on all files of archv type
# - to be run from data

source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/wavesice_ww3a
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC

# ==================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# =================================================================================

# In .bash_profile:
# 1. add MSCPROGS/bin to $PATH variable
#    (may need to compile MSCPROGS with "make" & "make install"
# 2. add this directory to $PATH variable
# 3. set environment variable $GIT_REPOS to the directory containing
#     this repository ("SWARP-routines")

# Info for hyc2proj
# - copied from MSCPROGS/Input & edited
h2p_in="$SWARP_ROUTINES/netcdf_production/Input"

# make a "working" directory
tday=$1
sdir=$THISFC2/$tday/netcdf
ddir=$THISFC2/$tday/bin
mkdir -p $sdir

cd $sdir
mkdir -p tmp
log=$sdir/tmp/convert_wav_log.txt
touch $log

mkdir -p tmp
cd tmp

# Info for hyc2proj
ln -s $h2p_in/proj.in .
ln -s $h2p_in/extract.archv_wav .
ln -s $h2p_in/extract.daily
ln -s $h2p_in/depthlevels.in .

# Info about grid from topo dir
tp4_input=/work/shared/nersc/msc/ModelInput/TOPAZ4/TP4a0.12/topo
ln -s $tp4_input/grid.info .
ln -s $tp4_input/regional.grid.a .
ln -s $tp4_input/regional.grid.b .
ln -s $tp4_input/regional.depth.a .
ln -s $tp4_input/regional.depth.b .

# CONVERT ARCH
if [ $(ls $ddir | tail -1) ]
then
   for f in $ddir/TP4archv_wav*.a
   do
      len=${#ddir}
      g=${f:$((len+1))}
      len=${#g}
      base=${g:0:$((len-2))}
      ln -s $f $ddir/$base.b .

      echo "********************************************************"
      echo "Running hyc2proj on $g..."
      echo "********************************************************"
      echo " "

      hyc2proj $g
      mv *.nc $sdir
   done
   echo "Binary files converted into .nc"       >> $log
else
   echo "NO archv* files in $ddir"              >> $log 
   mail -s "Convert_TP4archv_wav FAILED" $email < $log
fi

# CONVERT DAILY
if [ $(ls $ddir | tail -1) ]
then
   for f in $ddir/TP4DAILY*.a
   do
      len=${#ddir}
      g=${f:$((len+1))}
      len=${#g}
      base=${g:0:$((len-2))}
      ln -s $f $ddir/$base.b .

      echo "********************************************************"
      echo "Running hyc2proj on $g..."
      echo "********************************************************"
      echo " "

      hyc2proj $g
      mv *.nc $sdir
   done
   echo "Binary files converted into .nc"       >> $log
else
   echo "NO archv_wav* files in $ddir"              >> $log 
   mail -s "Convert_TP4archv_ww3a FAILED" $email < $log
fi

echo "********************************************************"
echo "Netcdf files in $sdir:"
ls -l $sdir
echo "********************************************************"

cd $sdir
rm -r tmp

