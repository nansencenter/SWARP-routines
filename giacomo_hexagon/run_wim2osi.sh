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
wrkdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/giacomo_hexagon/tmp
#wrkdir=/home/charlie/Documents/SWARP-routines/giacomo_stuff/jiping
results=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/validation/outputs

startdate='20150401'
enddate='20150930'
sdj=$(date --date="$startdate" +%j)
edj=$(date --date="$enddate" +%j)
ndays=$(expr $edj - $sdj)

for n in $(seq 0 $ndays)
do
   # Cleaning the data dir
   rm ${wrkdir}/OSI/*
   rm ${wrkdir}/MDL/*
   
   # Downloading the data
   hdate=$(date --date="$startdate + $n days" '+%Y%m%d')
   jdate=$(date --date="$hdate" '+%j')
   jdate=$(expr $jdate - 1)
   echo $jdate
   year=${hdate:0:4}
   o_obj=${osidir}/${year}_nh_polstere/ice_conc_nh_polstere-100_multi_${hdate}1200.nc
   m_obj=${mdldir}/${year}/SWARP_wavesice_forecast_${hdate}.tar.gz 
   if [ -f ${o_obj} ]
   then
      echo "OSISAF Found --> ${hdate}"
      cp ${o_obj} ${wrkdir}/OSI/
   else   
      echo "OSISAF NOT found --> ${hdate}"
   fi
   echo ""
   echo "Osisaf ${hdate} COMPLETE"
   echo ""
   if [ -f ${m_obj} ]
   then
      echo "MODEL Found --> ${hdate}"
      echo ""
      # NetCDF extraction
      steve=${hdate}/netcdf/TP4archv_wav_start${hdate}_000000Z_dump${hdate}_120000Z.nc
      echo "Unpacking..."
      tar -zxvf ${m_obj} ${steve} || echo "NetCDF not present"
      echo "Done."
      mv ${hdate} ${wrkdir}/MDL
      mv ${wrkdir}/MDL/${hdate}/netcdf/* ${wrkdir}/MDL
      rm -rf ${wrkdir}/MDL/${hdate}
   else   
      echo "MODEL NOT found --> ${hdate}"
   fi
   echo ""
   echo "Model ${hdate} COMPLETE"
   echo ""
   if [ -f ${wrkdir}/OSI/* ] && [ -f ${wrkdir}/MDL/*.nc ]
   then
      python wim2osi.py "${hdate}" "/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/giacomo_hexagon/outputs"
   else
      echo "Data not available"
   fi
done
