# redo some since winds are wrong
# - launch_all.sh will then redo them automatically when launched from cron

if [ 1 -eq 0 ]
then
   cd $TP4_REALTIME
   for E in `seq 39 56`
   do
      X=0${E:0:1}.${E:1:1}
      mpj=expt_$X/log/mpijob.out
      rm $mpj
   done
fi

if [ 1 -eq 1 ]
then
   # some of the 2015 060-067 results have got mixed up with other results somehow
   # - cleaning other data dir's
   for E in `seq 16 59`
   do
      X=0${E:0:1}.${E:1:1}
      cd $TP4_REALTIME/expt_$X/data
      pwd
      rm -f TP4DAILY_2015_060*
      rm -f TP4archv.2015_06*
   done
fi
