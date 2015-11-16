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

ftpsite="ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/wavewatch3/iowaga/HINDCAST/ARCTIC"
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
ww3dir1=$ww3dir/analysis_m1
mkdir -p $ww3dir1
ww3dir2=$ww3dir/analysis_m2
mkdir -p $ww3dir2
ww3dir3=$ww3dir/forecast
mkdir -p $ww3dir3

#Ef subfolders
ww3dirEF0=$ww3dirEF/originals
mkdir -p $ww3dirEF0
ww3dirEF1=$ww3dirEF/analysis_m1
mkdir -p $ww3dirEF1
ww3dirEF2=$ww3dirEF/analysis_m2
mkdir -p $ww3dirEF2
ww3dirEF3=$ww3dirEF/forecast
mkdir -p $ww3dirEF3

if [ $2 -eq 1 ]
then
   n2=6
else
   n2=7
fi

# working directory
mkdir -p /work/timill/tmp
mkdir -p /work/timill/tmp/ww3a
cd /work/timill/tmp/ww3a

# lists of forecast files
list1=''
list2=''

for n in `seq -1 $n2`   #all files
#for n in `seq -1 0`    #just analysis files
#for n in `seq 1 $n2`   #just forecast files
do
   dd=`date --date="$fdate ${n}day" "+%Y%m%d"`
   f1=SWARP_WW3_ARCTIC-12K_${dd}.nc
   f2=SWARP_WW3_ARCTIC-12K_${dd}_ef.nc

   # ==============================
   wget $fdir/$f1
   if [ -f $f1 ]
   then
      echo $f1 >> $log
   else
      echo "$f1 MISSING" >> $log
   fi

   wget $fdir/$f2
   if [ -f $f2 ]
   then
      echo $f2 >> $log
   else
      echo "$f2 MISSING" >> $log
   fi
   echo " " >> $log
   # ==============================

   if [ $n -ge 1 ]
   then
      list1="$list1 $f1"
      list2="$list2 $f2"
   elif [ $n -eq 0 ]
   then
      #"yesterday's" file (-1 = "m1")
      yr=`date --date=$dd +%Y`

      # =========================================
      f=$f1
      g=U$f1
      dir0=$ww3dir0/$fdir0
      dir1=$ww3dir1

      # unpack & change fill value
      ncpdq -U $f -o $g
      ncatted -O -a _FillValue,,o,f,-32767 $g
      
      # move original
      mkdir -p $dir0
      echo mv $f $dir0
      mv $f $dir0

      # repack, reformat & move
      ncpdq $opts -o $f $g
      echo mv $f $dir1
      mv $f $dir1
      rm $g
      # =========================================

      # =========================================
      f=$f2
      g=U$f2
      dir0=$ww3dirEF0/$fdir0
      dir1=$ww3dirEF1

      # unpack & change fill value
      ncpdq -U $f -o $g
      ncatted -O -a _FillValue,,o,f,-32767 $g
      
      # move original
      mkdir -p $dir0
      echo mv $f $dir0
      mv $f $dir0

      # repack, reformat & move
      ncpdq $opts -o $f $g
      echo mv $f $dir1
      mv $f $dir1
      rm $g
      # =========================================

   elif [ $n -eq -1 ]
   then
      #"day-before-yesterday's" file (-2 = "m2")
      yr=`date --date=$dd +%Y`

      # =========================================
      f=$f1
      g=U$f1
      dir0=$ww3dir0/$fdir0
      dir1=$ww3dir2

      # unpack & change fill value
      ncpdq -U $f -o $g
      ncatted -O -a _FillValue,,o,f,-32767 $g
      
      # move original
      mkdir -p $dir0
      echo mv $f $dir0
      mv $f $dir0

      # repack, reformat & move
      ncpdq $opts -o $f $g
      echo mv $f $dir1
      mv $f $dir1
      rm $g
      # =========================================

      # =========================================
      f=$f2
      g=U$f2
      dir0=$ww3dirEF0/$fdir0
      dir1=$ww3dirEF2

      # unpack & change fill value
      ncpdq -U $f -o $g
      ncatted -O -a _FillValue,,o,f,-32767 $g
      
      # move original
      mkdir -p $dir0
      echo mv $f $dir0
      mv $f $dir0

      # repack, reformat & move
      ncpdq $opts -o $f $g
      echo mv $f $dir1
      mv $f $dir1
      rm $g
      # =========================================
   fi
done

if [ -z "$list1" ]
#if [ 1 -eq 1 ]
then
   echo No parameter files to process
else

   # ====================================================
   # Unpacking parameter files...
   List="$list1"
   dir0=$ww3dir0/$fdir0
   dir1=$ww3dir3
   tf=tmp.nc
   ofil1=SWARP_WW3_ARCTIC-12K_${fdate}.fc.nc

   list=""
   echo Unpacking parameter files...
   for f in $List
   do
      g=U$f
      ncpdq -U $f -o $g
      list="$list $g"
      #
      echo mv $f $dir0
      mv $f $dir0
   done
   
   # combine into single file
   echo ncrcat -o $tf $list
   ncrcat -o $tf $list
   
   # Packing/reformatting...
   echo Packing/reformatting...
   echo ncpdq $opts -o $ofil1 $tf
   ncatted -O -a _FillValue,,o,f,-32767 $tf
   ncpdq $opts -o $ofil1 $tf
   rm $list $tf
   
   # move forecast file to final location
   echo mv $ofil1 $dir1
   mv $ofil1 $dir1
   # ====================================================

fi

if [ -z "$list2" ]
#if [ 1 -eq 1 ]
then
   echo No Ef files to process
else

   # ====================================================
   # Unpacking Ef files...
   List="$list2"
   dir0=$ww3dirEF0/$fdir0
   dir1=$ww3dirEF3
   tf=tmp.nc
   ofil1=SWARP_WW3_ARCTIC-12K_${fdate}_ef.fc.nc

   list=""
   echo Unpacking Ef files...
   for f in $List
   do
      g=U$f
      ncpdq -U $f -o $g
      list="$list $g"
      #
      echo mv $f $dir0
      mv $f $dir0
   done
   
   # combine into single file
   echo ncrcat -o $tf $list
   ncrcat -o $tf $list
   
   # Packing/reformatting...
   echo Packing/reformatting...
   echo ncpdq $opts -o $ofil1 $tf
   ncatted -O -a _FillValue,,o,f,-32767 $tf
   ncpdq $opts -o $ofil1 $tf
   rm $list $tf
   
   # move forecast file to final location
   echo mv $ofil1 $dir1
   mv $ofil1 $dir1
   # ====================================================

fi
