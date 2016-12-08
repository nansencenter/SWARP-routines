#!/bin/bash
# each day this script gets yesterday's OSISAF concentration from myocean

# get yesterday's date
ydate=`date -d "yesterday" '+%Y%m%d'`
ndays=5
if [ $# -ge 1 ]
then
   ydate=$1
fi
if [ $# -ge 2 ]
then
   ndays=$2
fi

# where the conc files are stored
hex_dir=/work/shared/nersc/msc/OSI-SAF
ME=`readlink -f $0`
Vdir=`dirname $ME`

CMEMS=1
if [ $CMEMS -eq 1 ]
then
   # where the files are stored in the CMEMS portal
   mo_dir0="Core/SEAICE_GLO_SEAICE_L4_NRT_OBSERVATIONS_011_001/METNO-GLO-SEAICE_CONC-NORTH-L4-NRT-OBS/"
   LNAME=cmems_ncftp # name of bookmark
else
   # where the files are stored in the myocean portal
   mo_dir0='SIW-TAC/SIW-OSISAF-GLO-SIT_SIE_SIC-OBS/conc/'
   LNAME=myocean # name of bookmark
fi

# check the last 5 days
for day in `seq 0 $ndays`
do
   cdate=`date -d "$ydate -${day}days" +%Y%m%d`
   cyear=`date -d "$cdate" '+%Y'`
   cmon=`date -d "$cdate" '+%m'`
   cday=`date -d "$cdate" '+%d'`

   # make yearly dir for file to go into
   dir0=$hex_dir/${cyear}_nh_polstere
   if [ ! -d dir0 ]
   then
      mkdir -p $dir0
      chmod -R g+r $dir0
      chmod -R g+w $dir0
   fi
   cd $dir0

   mo_dir=$mo_dir0/$cyear/$cmon
   fname=ice_conc_nh_polstere-100_multi_$cyear$cmon${cday}1200.nc
   if [ ! -f $fname ]
   then
      # if file isn't there, get from myocean portal
      LIST="$mo_dir/$fname"

#########################################################################
# get the file
# (can't handle too many at once it seems)
# make a text file (between "<<EOF" and "EOF")
  cat > ncftp.in<<EOF
open $LNAME
get $LIST
set confirm-close no
bye
EOF
         ncftp < ncftp.in # text file passed into ncftp
         rm ncftp.in
#########################################################################

   fi
done

#########################################################################
# Launch validation script
# - compares today's observation to relevant forecasts
$Vdir/ice_edge_OSISAF_1obs.sh $ydate
#########################################################################
