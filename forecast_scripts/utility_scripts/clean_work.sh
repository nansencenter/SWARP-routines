# delete results that are older than 21 days
# (contents will have been already deleted anyway)
dt=`date --date="-21 days" +%Y%m%d`

SR="/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines"
source $SR/source_files/hex_vars.src

# old ice-only FCs
dd=$TP4_REALTIME/../results/TP4a0.12/ice_only/work
cd $dd
echo ' '
pwd
echo ' '

dirs=(*)
ndirs=${#dirs[@]}

for i in `seq 0 $((ndirs-1))`
do
   dti=${dirs[i]}
   dti2=${dti:0:8}
   if [ $dti2 -le $dt ]
   then
      echo rm -r $dti
      rm -r $dti
   fi
done

# old waves-ice FCs
dd=$TP4_REALTIME/../results/TP4a0.12/wavesice/work
cd $dd
echo ' '
pwd
echo ' '

dirs=(*)
ndirs=${#dirs[@]}

for i in `seq 0 $((ndirs-1))`
do
   dti=${dirs[i]}
   dti2=${dti:0:8}
   if [ $dti2 -le $dt ]
   then
      echo rm -r $dti
      rm -r $dti
   fi
done

# old WAMNSEA checking figures
dd=$TP4_REALTIME/../check_wamnsea
cd $dd
echo ' '
pwd
echo ' '

dirs=(*)
ndirs=${#dirs[@]}

for i in `seq 0 $((ndirs-1))`
do
   dti=${dirs[i]}
   dti2=${dti:0:8}
   if [ $dti2 -le $dt ]
   then
      echo rm -r $dti
      rm -r $dti
   fi
done
