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
   echo "Usage: ww3_arctic_get_best_est.sh [date]"
   echo "*date in yyyymmdd format"
   exit
fi

# update best estimate from input date
dd=$1
dy=${dd:0:4}
fdir="$ftpsite/$dy/best_estimate"

# forecasts
ftpsite="ftp://ftp.ifremer.fr/ifremer/ww3/FORECAST/ARCTIC/" # link updated 19/10/2016
fdir="$ftpsite/$dy/best_estimate"

#main folders
WW3A="/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC"
ww3dir="$WW3A/${dy}"
ww3dirEF="$WW3A/${dy}_ef"
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
ww3dir0=$ww3dir/best_est_orig
mkdir -p $ww3dir0

#Ef subfolders
ww3dirEF0=$ww3dirEF/best_est_orig
mkdir -p $ww3dirEF0


# ==============================
# do downloads
cd $ww3dir0
pwd
f[0]=SWARP_WW3_ARCTIC-12K_${dd}.nc
echo "wget $fdir/${f[0]}"
wget $fdir/${f[0]}
F[0]=`readlink -f ${f[0]}`
D[0]=`readlink -f ../analysis_m1`

cd $ww3dirEF0
pwd
f[1]=SWARP_WW3_ARCTIC-12K_${dd}_ef.nc
echo "wget $fdir/${f[1]}"
wget $fdir/${f[1]}
F[1]=`readlink -f ${f[1]}`
D[1]=`readlink -f ../analysis_m1`
# ==============================

# convert
opts="--fl_fmt=classic" #netcdf3 classic
mkdir -p tmp
cd tmp
for n in `seq 0 1`
do
   f0=${f[$n]} 
   F0=${F[$n]} 
   g=U$f0

   #unpack, and change fill value
   echo ncpdq -U $F0 -o $g
   ncpdq -U $F0 -o $g
   ncatted -O -a _FillValue,,o,f,-32767 $g

   # repack, reformat
   echo ncpdq $opts -o $f0 $g
   ncpdq $opts -o $f0 $g

   # move to analysis_m1
   echo mv $f0 ${D[$n]}
   mv $f0 ${D[$n]}
done

# clean
cd ..
rm -r tmp
