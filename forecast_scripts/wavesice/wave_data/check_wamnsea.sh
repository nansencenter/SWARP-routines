# $SWARP_ROUTINES defined in .bash_profile or crontab

# get other variables
source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -eq 1 ]
then
   rerun=1
else
   rerun=0
fi

tday=`date +%Y%m%d`
dd=$TP4_REALTIME/../check_wamnsea/$tday
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

tyear=`date +%Y`
ff="$wmsc/WAMNSEA/$tyear/forecasts/wam_nsea.fc.$tday.nc"
if [ ! -f $ff ]
then
   # no wave fc yet
   exit
fi

# load python and launch check_wamnsea.py
[ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
module load python/2.7.9-dso
$python $FORECAST/wavesice/wave_data/check_wamnsea.py
