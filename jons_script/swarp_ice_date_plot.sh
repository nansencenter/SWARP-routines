#!/bin/bash
# Read in date & run pyton scripts to plot fice hice validation fig
# source HYCOMreg & $FCtype

#THIS_SRC=`readlink -f $1` #get absolute path
#source $THIS_SRC
cdate=$3
#cdate=20151218
# figure base name
fign=Hindcast_

# Set directory to result file
hres=/work/timill/RealTime_Models/results_hindcasts

# set current day, year and topaz year day
cyear=$(date --date=$cdate +%Y)
cjday=10#$(date --date=$cday +%j)
cjday0=$(( $cjday - 1))
cjday=$(printf '%3.3d' $jday0)
# get last monday => topaz day!
sdate=$(date --date=$cdate +%Y%m%d)
wd=$(date --date=$sdate +%u)
while [ ! "$wd" == 1 ]
do
 wd=$(date --date=$sdate +%u)
 sdate=`date --date="$sdate -1 days" "+%Y%m%d"`
  #cday=`date --date="last Monday -7 days" "+%Y%m%d"`
 echo "date=$sdate Weekday=$wd"
done
echo "week day should be 1: wd=$wd"
# Set "start of week" day, year and topaz year day
syear=$(date --date=$sdate +%Y)
sjday=10#$(date --date=$sdate +%j)
sjday0=$(( $sjday - 1))
sjday=$(printf '%3.3d' $sjday0)
# local week folder
wdir=${syear}_${sjday}
echo "Plotting for date=${cdate} (TYD=${cjday})"
echo "Week folder=${wdir}"

# Daily binary .a files for cdate
adir=${hres}/{$HYCOMreg}/$FCtype/${syear}_GOOD/${wdir}/binaries/DAILY
afil=TP4DAILY_${syear}_${sjday}_${cyear}_${cjday}.a

# netcdf output file
ncdir=${hres}/{$HYCOMreg}/$FCtype/${syear}_GOOD/${wdir}/final_out
ncfil=SWARP_hindcast_ice_only_start${cdate}T000000Z.nc

# OSI-SAF data
osdir=/work-common/shared/nersc/msc/OSI-SAF/${cyear}_nh_polstere
osfil=ice_conc_nh_polstere-100_multi_${cdate}1200.nc

# SMOS data
smdir=/work/shared/nersc/msc/SMOS/netcdf/smos_sea_ice_thickness_v2.1_north
smfil=SMOS_Icethickness_v2.1_north_${cdate}.nc

echo "adir/afil ${adir}/${afil}"
if [ -n ${adir}/${afil} ]; then echo "${adir}/${afil} exist"; fi
echo "ncdir/ncfil ${ncdir}/${ncfil}"
if [ -n ${ncdir}/${ncfil} ]; then echo "${ncdir}/${ncfil} exist"; fi
echo "osdir/osfil ${osdir}/${osfil}"
if [ -n ${osdir}/${osfil} ]; then echo "${osdir}/${osfil} exist"; fi
echo "smdir/smfil ${smdir}/${smfil}"
if [ -n ${smdir}/${smfil} ]; then echo "${smdir}/${smfil} exist"; fi


