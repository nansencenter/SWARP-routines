# do download from Ifremer ftp site
# either:
# - run after 1635 on $fdate
# or:
# - run after 0520 on $fdate +1day

if [ $# -lt 2 ]
then
   echo "Usage: ww3_arctic_download.sh [date] [cycle]"
   echo "*date in yyyymmdd format"
   echo "*cycle = 1 (1200 cycle) or 2 (0000 cycle)"
   exit
fi

fdate=$1
fyear=${fdate:0:4}

ftpsite="ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/wavewatch3/iowaga/HINDCAST/ARCTIC"
fdir="$ftpsite/$fyear/${fdate}00"

if [ $2 -eq 1 ]
then
   n2=6
else
   n2=7
fi


# cd /work/shared/nersc/msc/WAMNSEA/$fyear/analysis
cd /work/timill/tmp
for n in `seq -1 $n2`
do
   dd=`date --date="$fdate ${n}day" "+%Y%m%d"`
   wget $fdir/SWARP_WW3_ARCTIC-12K_${dd}.nc
   wget $fdir/SWARP_WW3_ARCTIC-12K_${dd}_ef.nc
done
