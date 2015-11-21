#!/bin/bash
# This script will collect and archive the results of the local TP4 model

# =============================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# =============================================================================

# ===================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/ice_only         # scripts
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC
# ===================================================================================


# TODO use py_funs/forecast_tools.fc2png_all
ncfil=$THISFC2/final_output/SWARP*.nc
pdir=. # directory with pngs
gdir=. # directory to save gif to

for vbl in icec icetk
do
   $python nc2png_all.py -i $ncfil -v $vbl -odir $pdir -vec uice
   convert=/usr/bin/convert
   $convert -delay 15 -loop 0 *.png $vbl.gif
done
