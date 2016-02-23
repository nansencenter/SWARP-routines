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
startdate='19910101'
enddate='20141231'
sdj=$(date --date="$startdate" +%j)
edj=$(date --date="$enddate" +%j)
ndays=$(expr $edj - $sdj)

for n in $(seq 0 $ndays)
do
   # Downloading the data
   hdate=$(date --date="$startdate + $n days" '+%Y%m%d')
   jdate=$(date --date="$hdate" '+%j')
   jdate=$(expr $jdate - 1)
   echo $jdate
   year=${hdate:0:4}
   o_obj=${osidir}/ice_conc_nh_polstere-100_multi_${hdate}1200.nc
   if [ -f ${o_obj} ]
   then
      echo "Found --> ${hdate}"
      cp ${o_obj} ${wrkdir}
   else   
      echo "NOT found --> ${hdate}"
   fi
   echo ""
   echo "Osisaf ${hdate} COMPLETE"
   echo ""
   if [ -f ${wrkdir}/ice_conc_nh_polstere-100_multi_${hdate}1200.nc ]
   then
      python MIZwidth.py 'OSI' '${hdate}' './OSI'
   else
      echo "Data not available"
   fi
done
