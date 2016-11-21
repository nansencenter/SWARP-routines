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

if [ $# -lt 1 ]
then
   echo "Usage: ww3_arctic_download.sh [date]"
   echo "*date in yyyymmdd format"
   exit
fi

yday=$1                             # input yesterday's date
yyear=${yday:0:4}
fdate=`date --date="$yday +1day" "+%Y%m%d"`   # convert to today's date
fyear=${fdate:0:4}
fdir0=${yday}00
Fdir0=${fdate}T00

# forecasts
ftpsite="ftp://ftp.ifremer.fr/ifremer/ww3/FORECAST/ARCTIC/" # link updated 19/10/2016
fdir="$ftpsite/$fyear/$Fdir0"

#main folders
WW3A="/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC"
ww3dir="$WW3A/${yyear}"
ww3dirEF="$WW3A/${yyear}_ef"
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
mkdir -p $ww3dir0/$fdir0

#Ef subfolders
ww3dirEF0=$ww3dirEF/originals
mkdir -p $ww3dirEF0
mkdir -p $ww3dirEF0/$fdir0


for n in `seq -1 6`   #all files
#for n in `seq -1 0`    #just analysis files
#for n in `seq 1 6`   #just forecast files
do
   dd=`date --date="$fdate ${n}day" "+%Y%m%d"`
   f1=SWARP_WW3_ARCTIC-12K_${dd}.nc
   f2=SWARP_WW3_ARCTIC-12K_${dd}_ef.nc

   # ==============================
   cd $ww3dir0/$fdir0
   echo "wget $fdir/$f1"
   wget $fdir/$f1
   if [ -f $f1 ]
   then
      echo $f1 >> $log
   else
      echo "$f1 MISSING" >> $log
   fi

   cd $ww3dirEF0/$fdir0
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
