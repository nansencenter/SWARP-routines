# script to set-up all experiments
# needs to be run from crontab since nesting is done to same folder
# -> can be crashes if files are opened by 2 expts at the same time
source $SWARP_ROUTINES/source_files/hex_vars.src
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/TP4a0.12
HCtype=wavesice

print_info=0
user=timill #TODO source hidden

ONEBYONE=0
if [ $ONEBYONE -eq 1 ]
   # check if another job is running
   # - should be OK in this case (no dump of nesting files)
   then
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
fi

P=$TP4_REALTIME
cd $P
refno=4
xref=$P/expt_01.$refno           # reference expt directory  (has restarts)
bref=$P/Build_V2.2.12_X01.$refno # reference Build directory (has hycom executable compiled)

E0=$((1$refno-1)) # increment this
X=01.$refno
Eno=0$E0
rno0=1

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

logdir=$hd/$HCtype/logs
mkdir -p $logdir
log=$logdir/hc.log
rm -f $log
echo `date` >> $log

#########################################################
# ECMWF files in 2015 skip between
# Saturday 22-Aug 12:00 and Sunday 23-Aug 06:00
# ie 4 missing records
bad_dates[0]=20150817 #Monday before this
bad_error[0]="Gap in ECMWF forcing 20150822-23"

# don't want to do last Monday's version 
# - winds are forecasts part of the week
bad_dates[1]="`date --date="last Monday" +%Y%m%d`"
bad_error[1]="Last Monday"

bad_dates[2]=20150601
bad_error[2]="Possible problem with restart"

# no WAM waves: 20150111 - no restart for then
# no WAM waves: 20150406
bad_dates[3]=20150406
bad_error[3]="No WAM waves on 20150406"

# no WAM waves: 20150612
bad_dates[4]=20150608
bad_error[4]="No WAM waves on 20150612"

# no WAM waves: 20151114
bad_dates[5]=20151109
bad_error[5]="No WAM waves on 20151114"

Nbad=${#bad_dates[@]}
# for loop_i in `seq 1 $Nbad`
# do
#    echo ${bad_dates[$((loop_i-1))]}
# done
#########################################################


if [ 1 -eq 0 ]
then
   # DO ALL
   E0_start=14
   E0_stop=2000
else
   # DO SOME
   # E0_start=36
   # E0_stop=56
   E0_start=18
   E0_stop=2000
fi


for rno in `seq 1 $Nlist`
do
   tfil0=${list[$((rno-1))]}
   E0=$((E0+1))
   Eno=0$E0
   X=0${E0:0:1}.${E0:1:1}
   xdir=$P/expt_$X            # reference expt directory  (has restarts)
   mpj=$xdir/log/mpijob.out

   tfil=`basename $tfil0`
   ryear=${tfil:10:4}
   rday=${tfil:15:3}
   rday2=10#$rday #base 10
   rday2=$((rday2+0))
   fday2=$((rday2+7))
   fday=`printf %3.3d $fday2`
   rdate=`date --date="${ryear}0101 +${rday2}days" "+%Y%m%d"`
   echo $tfil $rdate $rday $fday expt_$X

   # echo `basename $tfil0` $E0
   # continue

   if [ $E0 -lt $E0_start ] || [ $E0 -gt $E0_stop ]
   then
      echo "$E0: skipping"
      continue
   else
      echo "$E0: TODO"
   fi

   BAD_DAY=0
   for nbad in `seq 0 $((Nbad-1))`
   do
      if [ $rdate -eq ${bad_dates[$nbad]} ]
      then
         BAD_DAY=1
         BAD_ERROR=${bad_error[$nbad]}
         break
      fi
   done

   if [ $BAD_DAY -eq 1 ]
   then
      echo " "
      echo "Not launching week ${ryear}_${rday} ($rdate)"
      echo $BAD_ERROR
      echo " "
      continue
   fi

   if [ -f $mpj ]
   then

      # check if model crashed
      stat=`cat $mpj |grep "normal stop"`
      if [ ${#stat} -gt 0 ]
      then
         # didn't crash
         # - check if model output has the THERM_DIAG outputs
         stat=`cat $xdir/data/TP4archv_wav.${ryear}_${rday}_000000Z.b |grep taux_wav`
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
      cd $xdir
      echo "$tfil $ryear $rday $fday"  >> $log
      echo "Launching job from `pwd`"  >> $log
      echo " "                         >> $log

      if [ $print_info -eq 1 ]
      then
         echo "Launching job from `pwd`"
         echo " "                       
      fi

      rm -f log/*
      # cp $hd/$HCtype/inputs/archv.extract data
      cp $hd/$HCtype/inputs/extract.archv_wav data

      echo $qsub pbsjob.sh
      $qsub pbsjob.sh
      if [ $ONEBYONE -eq 1 ]
      then
         exit
      fi
   fi

done
