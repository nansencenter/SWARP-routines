# $SWARP_ROUTINES defined in .bash_profile or crontab

# get other variables
source $SWARP_ROUTINES/source_files/hex_vars.src
WW3A="/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC"

if [ $# -eq 1 ]
then
   rerun=1
else
   rerun=0
fi

yr=`date "--date=yesterday" "+%Y"`
yday=`date "--date=yesterday" "+%Y%m%d"`
dd=$TP4_REALTIME/../check_ww3arctic/$yday
if [ -d $dd ]
then
   if [ $rerun -eq 1  ]
   then
      echo Script has already run today...
      echo "...re-running"
   else
      # echo Script has already run today...
      # echo "...quitting"
      exit
   fi
fi

ff="$WW3A/$yr/forecast/wam_nsea.fc.$yday.nc"
if [ ! -f $ff ]
then
   # no wave fc yet
   exit
fi

# launch check_ww3arctic.py
$python $FORECAST/alert_waves_ww3arctic/check_ww3arctic.py
