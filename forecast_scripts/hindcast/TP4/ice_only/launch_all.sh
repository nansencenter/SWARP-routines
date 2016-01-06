# script to set-up all experiments
# needs to be run from crontab since nesting is done to same folder
# -> can be crashes if files are opened by 2 expts at the same time
source $SWARP_ROUTINES/source_files/hex_vars.src
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/TP4
print_info=0
user=timill #TODO source hidden

lst1=(`$qstat | grep $user |grep TP4 |grep HC |grep " Q "`) # jobs in queue
lst2=(`$qstat | grep $user |grep TP4 |grep HC |grep " R "`) # running jobs
N1=${#lst1[@]}
N2=${#lst2[@]}

if [ $N1 -gt 0 ] || [ $N2 -gt 0 ]
then
   if [ $print_info -eq 1 ]
   then
      echo " "
      echo "Job in progress already - stopping"
      echo " "
   fi
   exit
elif [ $print_info -eq 1 ]
then
   echo $N1 $N2
   echo `$qstat | grep $user |grep TP4 |grep HC |grep " Q "`
   echo `$qstat | grep $user |grep TP4 |grep HC |grep " R "`
fi

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

logdir=$hd/ice_only/logs
mkdir -p $logdir
log=$logdir/hc.log
rm -f $log
echo `date` >> $log

for tfil0 in $list
do
   Eno=0$E0
   X=0${E0:0:1}.${E0:1:1}
   xdir=$P/expt_$X            # reference expt directory  (has restarts)

   cd $xdir
   mpj=log/mpijob.out

   tfil=`basename $tfil0`
   ryear=${tfil:10:4}
   rday=${tfil:15:3}
   rday2=10#$rday #base 10
   fday2=$((rday2+7))
   fday=`printf %3.3d $fday2`
   echo $tfil $ryear $rday $fday expt_$X
   # echo $X; exit # ????

   # ECMWF files in 2015 skip between
   # Tuesday 22-Aug 12:00 and Wednesday 23-Aug 06:00
   # ie 4 missing records
   ecmwf_missing=20150822
   jday_missing=10#`date --date=$ecmwf_missing +%j`
   jday_missing=$((jday_missing-1))
   if [ $jday_missing -gt $((rday2+0)) ] && [ $jday_missing -lt $fday2  ]
   then

      echo "can't run week ${ryear}_${rday}"
      echo "- no ECMWF forcing for $ecmwf_missing (${ryear}_${jday_missing})"
      launch=0

   elif [ $rday -eq 151 ]
   then

      echo "Possible problem with restart file $tfil"
      launch=0

   elif [ -f $mpj ]
   then

      # check if model crashed
      stat=`cat $mpj |grep "normal stop"`
      if [ ${#stat} -gt 0 ]
      then
         # didn't crash
         # - check if model output has the THERM_DIAG outputs
         stat=`cat $xdir/data/TP4archv.${ryear}_${rday}_000000Z.b |grep flx_ild`
         if [ ${#stat} -eq 0 ]
         then
            launch=1
         else
            echo "$xdir OK" >> $log
            launch=0
         fi

         if [ $print_info -eq 1 ] && [ $launch -eq 0 ]
         then
            echo "$xdir OK"
            echo " "
         fi
      else
         launch=1
      fi

   else
      launch=1
   fi

   if [ $launch -eq 1 ]
   then
      echo "$tfil $ryear $rday $fday"  >> $log
      echo "Launching job from `pwd`"  >> $log
      echo " "                         >> $log

      if [ $print_info -eq 1 ]
      then
         echo "Launching job from `pwd`"
         echo " "                       
      fi
      rm -f log/*
      $qsub pbsjob.sh
      exit
   fi

   # increment rno,E0
   rno=$((rno+1))
   E0=$((E0+1))

done
