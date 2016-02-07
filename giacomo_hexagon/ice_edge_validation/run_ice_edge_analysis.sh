#!/bin/bash
# This script will download the model and osisaf files needed for the time MIZ series.
# It will then launch the time_series_miz.py and analyze the datasets
# Two inputs are needed for the python script:
# 1) MISSION = 'AOD','ICP','DFP','AARI'
# 2) DADATE = 'YYYYMMDD'
# 3) OUT_DIR = '/path/to/outputs'
# The format used will be: python time_series_miz.py MISSION DADATE OUT_DIR
# e.g. python time_series_miz.py AOD 20150515 /migrate/ice_edge_validation

mdldir=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/wavesice
osidir=/work/shared/nersc/msc/OSI-SAF
#wrkdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/validation/data
wrkdir=/home/charlie/Documents/SWARP-routines/giacomo_stuff/jiping
results=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/validation/outputs

startdate='19910101'
enddate='20141231'
sdj=$(date --date="$startdate" +%j)
edj=$(date --date="$enddate" +%j)
ndays=$(expr $edj - $sdj)

for n in $(seq 0 $ndays)
do
   hdate=$(date --date="$startdate + $n days" '+%Y%m%d')
   jdate=$(date --date="$hdate" '+%j')
   jdate=$(expr $jdate - 1)
   echo $jdate
   year=${hdate:0:4}
   python ice_edge_analysis.py "arctic" "${hdate}" "./outputs"
done
