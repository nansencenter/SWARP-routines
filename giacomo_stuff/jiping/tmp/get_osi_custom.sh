#!/bin/bash
# Get requested OSI-SAF file

echo "Insert date yyyymmdd"
read target
year=${target:0:4}
month=${target:4:2}
day=${target:6:2}

# where the conc files are stored
tmp_dir=/Documents/SWARP-routines/giacomo_stuff/jiping/tmp
cd $tmp_dir

# where the files are stored in the myocean portal 
mo_dir='SIW-TAC/SIW-OSISAF-GLO-SIT_SIE_SIC-OBS/conc/'
mo_dir=$mo_dir/$year/$month
for nday in `seq -w 01 $day`
do
   fname=ice_conc_nh_polstere-100_multi_$year$month${nday}1200.nc
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
