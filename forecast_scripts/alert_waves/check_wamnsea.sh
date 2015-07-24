# $SWARP_ROUTINES defined in .bash_profile or crontab

# get other variables
source $SWARP_ROUTINES/source_files/hex_vars.src

tday=`date +%Y%m%d`
dd=$TP4_REALTIME/../check_wamnsea/$tday
if [ -d $dd ]
then
   # script has already run today
   exit
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
$python $FORECAST/alert_waves/check_wamnsea.py
