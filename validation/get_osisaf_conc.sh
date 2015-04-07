#!/bin/bash
# each day this script gets yesterday's OSISAF concentration from myocean

# get yesterday's date
yyear=`date -d "yesterday" '+%Y'`
ymon=`date -d "yesterday" '+%m'`
yday=`date -d "yesterday" '+%d'`
ydate=$yyear$ymon$yday

# where the conc files are stored
hex_dir=/work/shared/nersc/msc/OSI-SAF
cd $hex_dir
dir0=${yyear}_nh_polstere
if [ ! -d dir0 ]
then
   mkdir -p $dir0
   chmod -R g+r $dir0
   chmod -R g+w $dir0
fi
cd $dir0

# where the files are stored in the myocean portal 
mo_dir='SIW-TAC/SIW-OSISAF-GLO-SIT_SIE_SIC-OBS/conc/'
mo_dir=$mo_dir/$yyear/$ymon
for nday in `seq -w 01 $yday`
do
   fname=ice_conc_nh_polstere-100_multi_$yyear$ymon${nday}1200.nc
   if [ ! -f $fname ]
   then
      # if file isn't there, get from myocean portal
      LIST="$mo_dir/$fname"

#########################################################################
# get the file
# (can't handle too many at once it seems)
# make a text file (between "<<EOF" and "EOF")
  cat > ncftp.in<<EOF
open myocean 
get $LIST
set confirm-close no
bye
EOF
         ncftp < ncftp.in # text file passd into ncftp
         rm ncftp.in
#########################################################################

   fi

done
