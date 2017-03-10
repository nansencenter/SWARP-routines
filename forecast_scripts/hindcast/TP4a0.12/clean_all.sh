# script to set-up all experiments
# needs to be run from crontab since nesting is done to same folder
# -> can be crashes if files are opened by 2 expts at the same time
source $SWARP_ROUTINES/source_files/hex_vars.src
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/TP4
print_info=0
user=timill #TODO source hidden

P=$TP4_REALTIME
cd $P
refno=5
xref=$P/expt_01.$refno           # reference expt directory  (has restarts)
bref=$P/Build_V2.2.12_X01.$refno # reference Build directory (has hycom executable compiled)

E0=1$refno # don't clean ref dir
# E0=$((1$refno-1)) # do clean ref dir
X=01.$refno
Eno=0$E0
rno=0

cd $P
LIST=(expt_*)
Nexpts=${#LIST[@]}

for edir in ${LIST[@]}
do

   Bdir=Build_V2.2.12_X${edir:5:4}
   Enum1=10#${edir:5:2}
   Enum2=${edir:8:1}
   Enum=$((Enum1+0))$Enum2

   if [ $Enum -le $E0 ]
   then
      echo keeping $edir $Bdir
   else
      echo rm -r $edir $Bdir
      rm -r $edir $Bdir
   fi

done
