if [ $# -lt 2 ]
then
   echo "Usage: wamnsea_download.sh [date] [type ie an or fc]"
   exit
fi

fdate=$1
fyear=${fdate:0:4}

ftype=$2

if [ $ftype ==  "an" ]
then
   cd /work/shared/nersc/msc/WAMNSEA/$fyear/analysis
   cat > ncftp.in<<EOF
open myocean
cd /Wave/METNO-WAM-NORDIC/
get wam_nsea.an.$fdate.nc
set confirm-close no
bye
EOF
else
   cd /work/shared/nersc/msc/WAMNSEA/$fyear/forecasts
   cat > ncftp.in<<EOF
open myocean
cd /Wave/METNO-WAM-NORDIC/
get wam_nsea.fc.$fdate.nc
set confirm-close no
bye
EOF
fi

ncftp < ncftp.in
rm ncftp.in
