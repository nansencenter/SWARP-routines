echo in check_wamnsea.sh
SWARP_ROUTINES=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/

# get other variables
source $SWARP_ROUTINES/source_files/hex_vars.src

# load python and launch check_wamnsea.py
# module load python/2.7.9-dso
$python FORECAST/alert_waves/check_wamnsea.py
