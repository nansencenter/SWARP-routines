#!/bin/bash
# This script will collect and archive the results of the local TP4 model

# ===================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC=$SWARP_ROUTINES/forecast_scripts/ice_only         # scripts
THIS_SRC=$THISFC/inputs/THISFC.src
source $THIS_SRC
# ===================================================================================


# =============================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# =============================================================================

if [ $# -ne 1 ]
then
   echo "Usage: make_gifs.sh [date]"
   echo "-date in yyyymmdd format"
   exit
else
   cdate=$1
fi

ncfil=$THISFC2/$cdate/final_output/SWARPiceonly_forecast_start${cdate}T000000Z.nc
figdir=$THISFC2/$cdate/figures
mkdir -p $figdir

pdir=$figdir/pngs # directory with pngs
mkdir -p $pdir
gdir=$figdir/gifs # directory to save gifs to
mkdir -p $gdir

$python $THISFC/post/make_pngs.py --ncfile=$ncfil --outdir=$pdir

cd $pdir
for vdir in *
do
   cd $vdir
   $convert -delay 15 -loop 0 *.png $vdir.gif
   mv $vdir.gif $gdir
   cd ..
done
