#!/bin/bash
# Purpose of the script is to clean and redo the final products of wavesice *.nc
wdir=/work/timill/RealTime_Models/results/TP4a0.12/wavesice/work
wdir2=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/work
netdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/netcdf_production

# WAVESICE
cd $wdir
nfil=$(expr $(ls -1 | wc -l) - 1) #count the files -1 for the loop
echo $nfil
for f in $(seq 0 $nfil)
do
   dd=$(date --date="-$f day" '+%Y%m%d') #the count starts from today
   dd2=$(date --date="-$f day" '+%Y-%m-%d') #different format (needed by merge*.sh)
   cd $wdir/$dd/final_output
   rm *.nc
   /bin/bash $netdir/merge_TP4archv_wav.sh $dd $dd2 #create new final product
done

# ICE_ONLY
cd $wdir2
nfil2=$(expr $(ls -1 | wc -l) - 1) #count the files -1 for the loop
echo $nfil2
for f in $(seq 0 $nfil2)
do
   dd=$(date --date="-$f day" '+%Y%m%d')
   dd2=$(date --date="-$f day" '+%Y-%m-%d')
   cd $wdir2/$dd/final_output
   rm *.nc
   /bin/bash $netdir/merge_TP4archv.sh $dd $dd2
done

