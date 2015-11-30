# do download from Ifremer ftp site
# either:
# - run after 1635 on $fdate
# or:
# - run after 0520 on $fdate +1day
source $SWARP_ROUTINES/source_files/hex_vars.src
# opts="--fl_fmt=netcdf4_classic"
opts="--fl_fmt=classic" #netcdf3 classic

[ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
module load nco

if [ $# -lt 2 ]
then
   echo "Usage: ww3_arctic_download.sh [date] [cycle]"
   echo "*date in yyyymmdd format"
   echo "*cycle = 1 (1200 cycle) or 2 (0000 cycle)"
   exit
fi

fdate=$1
fyear=${fdate:0:4}
fdir0=${fdate}00

cycle=$2
if [ $cycle -eq 3 ]
then
   # hindcasts
   ftpsite="ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/wavewatch3//HINDCAST/ARCTIC"
else
   # forecasts (cleaned each month)
   ftpsite="ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/wavewatch3/iowaga/HINDCAST/ARCTIC"
fi
fdir="$ftpsite/$fyear/$fdir0"

#main folders
WW3A="/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC"
ww3dir="$WW3A/${fyear}"
ww3dirEF="$WW3A/${fyear}_ef"
mkdir -p $ww3dir
mkdir -p $ww3dirEF

# log file
log=$WW3A/latest_download.txt
echo "** Update WW3_ARCTIC wave data from ftp.ifremer.fr **"   >  $log
echo " Today is $(date "+%Y%m%d") $(date "+%H:%M:%S")"         >> $log
echo ""                                                        >> $log
echo " Contents of folder ${fdate}00:"                         >> $log
echo ""                                                        >> $log

#subfolders
ww3dir0=$ww3dir/originals
mkdir -p $ww3dir0

#Ef subfolders
ww3dirEF0=$ww3dirEF/originals
mkdir -p $ww3dirEF0

if [ $cycle -eq 1 ]
then
   n2=6
elif [ $cycle -eq 2 ]
then
   n2=7
elif [ $cycle -eq 3 ]
then
   # just get hindcasts
   n2=0
fi


for n in `seq -1 $n2`   #all files
#for n in `seq -1 0`    #just analysis files
#for n in `seq 1 $n2`   #just forecast files
do
   dd=`date --date="$fdate ${n}day" "+%Y%m%d"`
   f1=SWARP_WW3_ARCTIC-12K_${dd}.nc
   f2=SWARP_WW3_ARCTIC-12K_${dd}_ef.nc

   # ==============================
   cd $ww3dir0
   echo "wget $fdir/$f1"
   wget $fdir/$f1
   if [ -f $f1 ]
   then
      echo $f1 >> $log
   else
      echo "$f1 MISSING" >> $log
   fi

   cd $ww3dirEF0
   echo "wget $fdir/$f2"
   wget $fdir/$f2
   if [ -f $f2 ]
   then
      echo $f2 >> $log
   else
      echo "$f2 MISSING" >> $log
   fi
   echo " " >> $log
   # ==============================

done
