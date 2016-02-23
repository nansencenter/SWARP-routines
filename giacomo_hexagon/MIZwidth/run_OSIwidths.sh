#!/bin/bash
# 1) MISSION = 'OSI','MIC','MFD' as OSISAF, MODEL ICE CONC, MODEL FLOE size DISTRIBUTION
# 2) DADATE = 'YYYYMMDD'
# 3) OUT_DIR = '/path/to/outputs'
# The format used will be: python time_series_miz.py MISSION DADATE OUT_DIR
# e.g. python time_series_miz.py AOD 20150515 /work/banana_dir

osidir=/work/timill/giacomo/osisaf_repo
wrkdir=/work/timill/giacomo/MIZ_TIME_SERIES

# NOTE THIS DATES WORK FOR THE OSISAF REPO!
startdate='19910101'
enddate='20141231'
ndays=$(( $(date -d "$startdate" +%s) - $(date -d "$enddate" +%s) ))

for n in $(seq 0 $ndays)
do
   # Downloading the data
   hdate=$(date --date="$startdate + $n days" '+%Y%m%d')
   jdate=$(date --date="$hdate" '+%j')
   echo $hdate
   o_obj=${osidir}/ice_conc_nh_polstere-100_multi_${hdate}1200.nc
   if [ -f ${o_obj} ]
   then
      python MIZwidth.py 'OSI' ${hdate} './outputs/OSI'
   else
      echo "Data not available"
   fi
done
