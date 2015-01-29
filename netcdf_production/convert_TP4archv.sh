#!/bin/bash
# Runs hyc2proj on all files of archv type
# - to be run from data

# In .bash_profile:
# 1. add MSCPROGS/bin to $PATH variable
#    (may need to compile MSCPROGS with "make" & "make install"
# 2. add this directory to $PATH variable
# 3. set environment variable $GIT_REPOS to the directory containing
#     this repository ("SWARP-routines")

# Info for hyc2proj
# - copied from MSCPROGS/Input & edited
h2p_in="$GIT_REPOS/SWARP-routines/netcdf_production/Input"
ln -s $h2p_in/proj.in .
ln -s $h2p_in/extract.archv .
ln -s $h2p_in/depthlevels.in .

# Info about grid from topo dir
topo=../../topo
ln -s $topo/grid.info .
ln -s $topo/regional.grid.a .
ln -s $topo/regional.grid.b .
ln -s $topo/regional.depths.a .
ln -s $topo/regional.depths.b .

odir=archv_netcdf
mkdir -p $odir
for f in TP4archv*.a
do
   hyc2proj $f
   #
   yr=${f:9:4}
   day=${f:14:3}
   dt=`jultodate $day $yr 1 1`
   #
   hr=${f:18:2}
   fout=${f:0:8}_${dt}_${hr}.nc
   mv $fout $odir
done

echo "Netcdf files in $odir:"
ll $odir
