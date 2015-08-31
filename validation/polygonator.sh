#!/bin/bash
# This script will download the model and osisaf files needed for the time MIZ series.
# It will then launch the time_series_miz.py and analyze the datasets
# Two inputs are needed for the python script:
# 1) MISSION = 'AOD','ICP','DFP','AARI'
# 2) DADATE = 'YYYYMMDD'
# The format used will be: python time_series_miz.py MISSION DADATE
# e.g. python time_series_miz.py AOD 20150515

mdldir=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/wavesice
osidir=/work/shared/nersc/msc/OSI-SAF
wrkdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/validation/data
results=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/validation/outputs

startdate='20150601'
enddate='20150605'
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
      echo "Found --> ${hdate}"
      cp ${o_obj} ${wrkdir}/OSI/
   else   
      echo "NOT found --> ${hdate}"
   fi
   echo ""
   echo "Osisaf ${hdate} COMPLETE"
   echo ""
   if [ -f ${m_obj} ]
   then
      echo "Found --> ${hdate}"
      echo ""
#      # Binary extraction
#      steve=${hdate}/bin/TP4archv_wav.${year}_${jdate}_120000.a
#      alan=${hdate}/bin/TP4archv_wav.${year}_${jdate}_120000.b
#      echo "Unpacking..."
#      tar -zxvf ${m_obj} ${steve} || echo "Binary .a not present" 
#      tar -zxvf ${m_obj} ${alan} || echo "Binary .b not present" 
#      echo "Done."
#      mv ${hdate} ${wrkdir}/MDL
#      mv ${wrkdir}/MDL/${hdate}/bin/* ${wrkdir}/MDL
#      rm -rf ${wrkdir}/MDL/${hdate}
      # NetCDF extraction
      steve=${hdate}/netcdf/TP4archv_wav_start${hdate}_000000Z_dump${hdate}_120000Z.nc
      echo "Unpacking..."
      tar -zxvf ${m_obj} ${steve} || echo "NetCDF not present"
      echo "Done."
      mv ${hdate} ${wrkdir}/MDL
      mv ${wrkdir}/MDL/${hdate}/netcdf/* ${wrkdir}/MDL
      rm -rf ${wrkdir}/MDL/${hdate}
   else   
      echo "NOT found --> ${hdate}"
   fi
   echo ""
   echo "Model ${hdate} COMPLETE"
   echo ""
   if [ -f ${wrkdir}/OSI/* ] && ([ -f ${wrkdir}/MDL/*.a ] || [ -f ${wrkdir}/MDL/*.nc ])
   then
      # for ICP
      python time_series_miz.py "ICP" "${hdate}"
#      mv -f ${results}/ICP/* /work/users/timill/RESULTS/POLYGONS/ICP
      ## for DFP
      #python time_series_miz.py "DFP" "${hdate}"
      #mv -f ${results}/DFP/* /work/users/timill/RESULTS/POLYGONS/DFP
      ## for AOD
      #python time_series_miz.py "AOD" "${hdate}"
      #mv -f ${results}/AOD/* /work/users/timill/RESULTS/POLYGONS/AOD
   else
      echo "Data not available"
   fi
done
