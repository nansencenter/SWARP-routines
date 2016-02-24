#!/bin/bash
# 1) MISSION = 'OSI','MIC','MFD' as OSISAF, MODEL ICE CONC, MODEL FLOE size DISTRIBUTION
# 2) DADATE = 'YYYYMMDD'
# 3) OUT_DIR = '/path/to/outputs'
# The format used will be: python time_series_miz.py MISSION DADATE OUT_DIR
# e.g. python time_series_miz.py AOD 20150515 /work/banana_dir

mdldir=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/wavesice
osidir=/work/timill/GIACOMO/osisaf_repo
wrkdir=/work/timill/GIACOMO/MIZ_TIME_SERIES
results=/work/timill/GIACOMO/MIZ_TIME_SERIES

# NOTE THIS DATES WORK FOR THE OSISAF REPO!
startdate='20150302'
enddate='20150930'
sdj=$(date --date="$startdate" +%j)
edj=$(date --date="$enddate" +%j)
ndays=$(expr $edj - $sdj)

for n in $(seq 0 $ndays)
do
   python WIMchar.py ${hdate} "/work/timill/giacomo/WIMchar/output"
done
