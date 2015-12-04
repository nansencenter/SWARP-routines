source $SWARP_ROUTINES/source_files/hex_vars.src

rdir="$TP4_REALTIME/../results/TP4a0.12/wavesice/work"
cd $rdir

for cdate in *
do
   ofil=$cdate/final_output/*UTC*.nc
   rm -f $ofil
   $SWARP_ROUTINES/netcdf_production/merge_TP4archv_wav.sh $cdate
   ls $cdate/final_output
   echo " "
done
