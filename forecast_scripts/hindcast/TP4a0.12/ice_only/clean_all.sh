# script to set-up all experiments
# needs to be run from crontab since nesting is done to same folder
# -> can be crashes if files are opened by 2 expts at the same time
source $SWARP_ROUTINES/source_files/hex_vars.src
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/TP4
print_info=0
user=timill #TODO source hidden

P=$TP4_REALTIME
cd $P
refno=4
xref=$P/expt_01.$refno           # reference expt directory  (has restarts)
bref=$P/Build_V2.2.12_X01.$refno # reference Build directory (has hycom executable compiled)

E0=$((1$refno-1)) # increment this
X=01.$refno
Eno=0$E0
rno=0

md=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/
list=(`ls $md/2015/*gz`)
list2=(`ls $md/2016/*gz`)
Nlist=${#list[@]}
inext=Nlist
for loop_i in `seq 1 ${#list2[@]}`
do
   list[$inext]=${list2[$((loop_i-1))]}
   inext=$((inext+1))
   Nlist=$((Nlist+1))
done

for rno in `seq 1 $Nlist`
do
   tfil0=${list[$((rno-1))]}
   E0=$((E0+1))
   echo " "
   echo `basename $tfil0` $E0
   # continue

   Eno=0$E0
   X=0${E0:0:1}.${E0:1:1}
   xdir=$P/expt_$X            # reference expt directory  (has restarts)
   if [ -d $xdir ]
   then
      rm -f $xdir/log/mpijob.out
      cd $xdir/data
      rm -f TP4DAILY* TP4archv*
      rm -f $xdir/SCRATCH/Ecmwf.nc
      # deleting this and changing REGION.src makes preprocess.sh reset path to ECMWFR_T799_bak from ECMWFR_T799
   fi

done
