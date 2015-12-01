# do download from Ifremer ftp site
# either:
# - run after 1635 on $fdate
# or:
# - run after 0520 on $fdate +1day
source $SWARP_ROUTINES/source_files/hex_vars.src
# opts="--fl_fmt=netcdf4_classic"
opts="--fl_fmt=classic" #netcdf3 classic
print_info=0

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

#main folders
WW3A="/work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC"
ww3dir="$WW3A/${fyear}"
ww3dirEF="$WW3A/${fyear}_ef"
mkdir -p $ww3dir
mkdir -p $ww3dirEF

# log file
log=$WW3A/latest_download.txt

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

if [ $cycle -eq 1 ]
then
   do_m1=1   #-1 result
   do_m2=1   #-2 result
   n2=6
elif [ $cycle -eq 2 ]
then
   do_m1=1   #-1 result
   do_m2=1   #-2 result
   n2=7
elif [ $cycle -eq 3 ]
then
   # just process hind cast
   do_m1=1   #-1 result
   do_m2=1   #-2 result
   n2=0
fi

dir0=$ww3dir0/${fdate}00
dir1=$ww3dirEF0/${fdate}00

# =================================================================================
# an-2:
if [ $do_m2 -eq 1 ]
then
   n=-1
   dd=`date --date="$fdate ${n}day" "+%Y%m%d"`
   ff=($dir0/SWARP_WW3_ARCTIC-12K_${dd}.nc $dir1/SWARP_WW3_ARCTIC-12K_${dd}_ef.nc)
   dir2s=($ww3dir2 $ww3dirEF2)

   for m in 0 1
   do
      dir2=${dir2s[$m]}
      f1=${ff[$m]}
      if [ $print_info -eq 1 ]
      then
         echo $f1
      fi

      # unpack & change fill value
      cd $dir2
      pwd

      f=`basename $f1`
      g=U$f
      if [ $print_info -eq 1 ]
      then
         echo "ncpdq -U $f1 -o $g"
      fi
      ncpdq -U $f1 -o $g
      ncatted -O -a _FillValue,,o,f,-32767 $g
            
      # repack, reformat
      if [ $print_info -eq 1 ]
      then
         echo "ncpdq $opts -o $f $g"
      fi
      ncpdq $opts -o $f $g
      rm $g
   done
fi
# # =================================================================================


# =================================================================================
# an-1:
if [ $do_m1 -eq 1 ]
then
   n=0
   dd=`date --date="$fdate ${n}day" "+%Y%m%d"`
   ff=($dir0/SWARP_WW3_ARCTIC-12K_${dd}.nc $dir1/SWARP_WW3_ARCTIC-12K_${dd}_ef.nc)
   dir2s=($ww3dir1 $ww3dirEF1)

   for m in 0 1
   do
      dir2=${dir2s[$m]}
      f1=${ff[$m]}
      if [ $print_info -eq 1 ]
      then
         echo $f1
      fi

      # unpack & change fill value
      cd $dir2
      pwd

      f=`basename $f1`
      g=U$f
      if [ $print_info -eq 1 ]
      then
         echo "ncpdq -U $f1 -o $g"
      fi
      ncpdq -U $f1 -o $g
      ncatted -O -a _FillValue,,o,f,-32767 $g
            
      # repack, reformat
      if [ $print_info -eq 1 ]
      then
         echo "ncpdq $opts -o $f $g"
      fi
      ncpdq $opts -o $f $g
      rm $g
   done
fi
# =================================================================================


# =================================================================================
# forecast
dir2s=($ww3dir3 $ww3dirEF3)
fcs=(SWARP_WW3_ARCTIC-12K_${fdate}.fc.nc SWARP_WW3_ARCTIC-12K_${fdate}_ef.fc.nc)
for n in `seq 1 $n2`
do
   dd=`date --date="$fdate ${n}day" "+%Y%m%d"`
   ff=($dir0/SWARP_WW3_ARCTIC-12K_${dd}.nc $dir1/SWARP_WW3_ARCTIC-12K_${dd}_ef.nc)

   for m in 0 1
   do
      dir2=${dir2s[$m]}
      f1=${ff[$m]}
      fout=${fcs[$m]}
      if [ $print_info -eq 1 ]
      then
         echo $f1
      fi

      # unpack & change fill value
      cd $dir2
      mkdir -p tmp
      cd tmp
      pwd

      f=`basename $f1`
      g=U$f
      if [ $print_info -eq 1 ]
      then
         echo "ncpdq -U $f1 -o $g"
      fi
      ncpdq -U $f1 -o $g
      ncatted -O -a _FillValue,,o,f,-32767 $g

      if [ $n -eq $n2 ]
      then
         # merge into 1 file
         gout=U$fout
         ncrcat *.nc $gout

         # repack, reformat
         if [ $print_info -eq 1 ]
         then
            echo "ncpdq $opts -o $fout $gout"
         fi
         ncpdq $opts -o $fout $gout

         # clean
         mv $fout ..
         cd ..
         rm -r tmp
      fi
   done
done
# =================================================================================
