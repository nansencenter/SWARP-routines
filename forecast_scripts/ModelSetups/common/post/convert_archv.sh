#!/bin/bash
# Runs hyc2proj on all files of archv type

# ===================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC
# ===================================================================================

# ==================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ==================================================================================

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
tday=$2
sdir=$THISFC2/$tday/netcdf
ddir=$THISFC2/$tday/bin
mkdir -p $sdir

cd $sdir
mkdir -p tmp
log=$sdir/tmp/convert_log.txt
touch $log
cd tmp

# Info for hyc2proj
ln -s $h2p_in/proj.in.polar_stereographic.$HYCOMreg proj.in
ln -s $h2p_in/extract.archv . 
ln -s $h2p_in/extract.daily . 
ln -s $h2p_in/depthlevels.in .

# Info about grid from topo dir
topo_input=$RTmods/$HYCOMreg/topo
ln -s $topo_input/grid.info .
ln -s $topo_input/regional.grid.a .
ln -s $topo_input/regional.grid.b .
ln -s $topo_input/regional.depth.a .
ln -s $topo_input/regional.depth.b .

# CONVERT ARCH
if [ $(ls $ddir | tail -1) ]
then
   for f in $ddir/${rungen}archv*.a
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

         echo hyc2proj $g
         echo cp *.nc $sdir
         echo " "
         hyc2proj $g
         cp *.nc $sdir
      done
else
   "No archv* files in $ddir" >> $log
   mail -s "convert_archv FAILED" $email < $log
fi

# CONVERT DAILY
if [ $(ls $ddir | tail -1) ]
then
   for f in $ddir/${rungen}DAILY*.a
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

         echo hyc2proj $g
         echo cp *.nc $sdir
         echo " "
         hyc2proj $g
         cp *.nc $sdir
      done
else
   "No DAILY* files in $ddir" >> $log
   mail -s "convert_archv FAILED" $email < $log
fi

cd $sdir
rm -r tmp

echo "********************************************************"
echo "Netcdf files in $sdir:"
ls -lh
echo "********************************************************"
