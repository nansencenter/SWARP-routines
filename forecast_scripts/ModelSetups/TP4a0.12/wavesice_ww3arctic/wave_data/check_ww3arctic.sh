# $SWARP_ROUTINES defined in .bash_profile or crontab

# get other variables
source $SWARP_ROUTINES/source_files/hex_vars.src
THISFC="$FORECAST/ModelSetups/TP4a0.12/wavesice_ww3arctic"
WW3A="$wmsc/WAVES_INPUT/WW3_ARCTIC"

if [ $# -eq 1 ]
then
   rerun=1
else
   rerun=0
fi

yr=`date "--date=yesterday" "+%Y"`
yday=`date "--date=yesterday" "+%Y%m%d"`
dd=$RTmods/check_ww3arctic/$yday
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

ff="$WW3A/$yr/forecast/SWARP_WW3_ARCTIC-12K_$yday.fc.nc"
if [ ! -f $ff ]
then
   # no wave fc yet
   exit
fi

subdir=$THISFC/wave_data
outdir=$RTmods/check_ww3arctic
mkdir -p $outdir

# load python and launch check_ww3arctic.py
[ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
module load python/2.7.9-dso
$python $subdir/check_ww3arctic.py --subdir=$subdir --outdir=$outdir
