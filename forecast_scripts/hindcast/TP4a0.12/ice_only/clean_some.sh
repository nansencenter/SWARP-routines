# redo some since winds are wrong
# - launch_all.sh will then redo them automatically when launched from cron
cd $TP4_REALTIME
for E in `seq 39 56`
do
   X=0${E:0:1}.${E:1:1}
   mpj=expt_$X/log/mpijob.out
   rm $mpj
done
