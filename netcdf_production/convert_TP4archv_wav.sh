#!/bin/bash
# Runs hyc2proj on all files of archv type
# - to be run from data

# EMAIL ADDRESS
address=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/fc_alert_email.txt
# ==================================================================================
email=$(cat $address)
# =================================================================================

# In .bash_profile:
# 1. add MSCPROGS/bin to $PATH variable
#    (may need to compile MSCPROGS with "make" & "make install"
# 2. add this directory to $PATH variable
# 3. set environment variable $GIT_REPOS to the directory containing
#     this repository ("SWARP-routines")

# Info for hyc2proj
# - copied from MSCPROGS/Input & edited
h2p_in="$GIT_REPOS/SWARP-routines/netcdf_production/Input"

# make a "working" directory
tday=$1
sdir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$tday/netcdf
ddir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work/$tday/bin
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
ln -s $h2p_in/depthlevels.in .

# Info about grid from topo dir
tp4_input=/work/shared/nersc/msc/ModelInput/TOPAZ4/TP4a0.12/topo
ln -s $tp4_input/grid.info .
ln -s $tp4_input/regional.grid.a .
ln -s $tp4_input/regional.grid.b .
ln -s $tp4_input/regional.depth.a .
ln -s $tp4_input/regional.depth.b .

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

echo "********************************************************"
echo "Netcdf files in $sdir:"
ls -l $sdir
echo "********************************************************"

cd $sdir
rm -r tmp

