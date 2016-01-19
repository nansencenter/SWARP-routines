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

E0=1$refno # increment this
X=01.$refno
Eno=0$E0
rno=0

md=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/2015/
list=`ls $md/*gz`

for tfil0 in $list
do
   Eno=0$E0
   X=0${E0:0:1}.${E0:1:1}
   xdir=$P/expt_$X            # reference expt directory  (has restarts)
   if [ -d $xdir ]
   then
      rm -f $xdir/log/mpijob.out
      cd $xdir/data
      rm -f TP4DAILY* TP4archv*
      rm $xdir/SCRATCH/Ecmwf.nc
      # deleting this and changing REGION.src makes preprocess.sh reset path to ECMWFR_T799_bak from ECMWFR_T799
   fi

   # increment rno,E0
   rno=$((rno+1))
   E0=$((E0+1))
done
