#!/bin/bash
# copy forecast file from hexagon to johansen

thour=`date +%H`
T0=2  # don't bother doing anything before 2am
T1=9  # email alert after this time
if [ $thour -lt $T0 ]
then
   exit
fi

if [ $# -eq 2 ]
then
   FCtype_hex=$1
   reg=$2
elif [ $# -eq 1 ]
then
   FCtype_hex=$1
   reg=TP4
else
   FCtype_hex=ice_only
   reg=TP4
fi

if [ $FCtype_hex == ice_only ]
then
   FCtype_joh="ice_only"               # type of FC on johansen
   FC_OUTPUT="SWARPiceonly_forecast"   # start of netcdf file

   # can rename here
   FC_OUTPUT2="SWARPiceonly_forecast"

elif [ $FCtype_hex == wavesice ]
then
   FCtype_joh="wavesice"               # type of FC on johansen
   FC_OUTPUT="SWARPwavesice_forecast"  # start of netcdf file

   # can rename here
   FC_OUTPUT2="SWARPwavesice_forecast"

elif [ $FCtype_hex == wavesice_ww3arctic ]
then
   FCtype_joh="wavesice"               # type of FC on johansen
   FC_OUTPUT="SWARPwavesice_WW3_forecast"   # start of netcdf file

   # can rename here
   FC_OUTPUT2="SWARPwavesice_forecast"

fi

if [ $reg == TP4 ]
then
   Reg=TP4a0.12
elif [ $reg == BS1 ]
then
   FCtype_joh="${FCtype_joh}_$reg"  # type of FC on johansen
   Reg=BS1a0.045
elif [ $reg == FR1 ]
then
   FCtype_joh="${FCtype_joh}_$reg"  # type of FC on johansen
   Reg=FR1a0.03
fi

print_info=1
if [ $print_info -eq 1 ]
then
   echo $reg $Reg
   echo $FCtype_hex $FCtype_joh
   echo $FC_OUTPUT $FC_OUTPUT2
fi

# ==================================================================================
# general variables
THREDDS="/mnt/10.11.12.232/Projects/SWARP/Thredds_data_dir" # final destination of forecast outputs
SRjoh="/Data/sim/tim/Projects/SWARP/SWARP-routines"         # location of SWARP routines on johansen
TP4rt="/work/timill/RealTime_Models/results/$Reg"           # location of forecast outputs on hexagon
hidden=$SRjoh/forecast_scripts/hidden                       # location of private info eg logs, ssh_info.src
FTP="/Data/FTPRoot/pub/Tim/SWARP_forecasts/$reg/"
FCemail=$hidden/FCemail.txt

# FC-specific variables
joh_dir=$THREDDS/$FCtype_joh  # directory on johansen where final file should go
tmp_dir=$joh_dir/tmp          # directory on johansen where final file + figures should go initially after copying
hex_dir=$TP4rt/$FCtype_hex    # directory on hexagon where final file ends up
# ==================================================================================

mkdir -p $tmp_dir

# get $keyname, $user
source $hidden/ssh_info.src
echo ho


# ==================================================================================
# EMAIL ADRESS
email=$(cat $FCemail)
# ==================================================================================

cplog=$hidden/logs/cp_${reg}_${FCtype_hex}.log
tday=$(date +%Y%m%d)

if [ -f $cplog ]
then
   if [ ! $tday == $(cat $cplog | sed '1!d') ]
   then
      mv $cplog $hidden/old_logs/$(cat $cplog | sed '1!d')
   fi
fi

if [ $(date +%A) == "Monday" ]
then
   rm -f $hidden/old_logs/*
fi

if [ ! -f "$cplog" ]
then
   touch $cplog
   echo $tday >> $cplog
fi

# finding the latest final product
# - check last $Nback days
Nback=0
wrn_count=0
for n in `seq 0 $Nback`
do
   echo "$(date +%H:%M) - looking for <<$FCtype_hex>> latest file" >> $cplog
   hdate=$(date --date="$n days ago" '+%Y%m%d')
   echo "Looking for product $hdate" >> $cplog
   hex_fil=${FC_OUTPUT}_start${hdate}T000000Z.nc   # final file on hexagon
   joh_fil=${FC_OUTPUT2}_start${hdate}T000000Z.nc  # final file on johansen

   # only do scp if file not present
   if [ ! -f $joh_dir/$joh_fil ]
   then
      echo "scp -i $HOME/.ssh/\$keyname \$user@hexagon.bccs.uib.no:$hex_dir/$hdate/final_output/$hex_fil $tmp_dir"
      scp -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hdate/final_output/$hex_fil $tmp_dir

      if [ -f $tmp_dir/$hex_fil ]
      then
         # if scp worked move it to THREDDS dir
         mv $tmp_dir/$hex_fil $joh_dir/$joh_fil
         chmod o+r $joh_dir/$joh_fil
         echo "Product found on $hdate!" >> $cplog
         echo "" >> $cplog
      else
         # if scp didn't work, give warning
         echo "No product on $hdate" >> $cplog
         wrn_count=$(expr $wrn_count + 1)
      fi

      #################################################################
      # gif's for website
      if [ $n -eq 0 ]
      then
         # copy gifs from latest forecast
         # TODO: mv to /Webdata and link to website
         echo scp -r -i $HOME/.ssh/\$keyname \$user@hexagon.bccs.uib.no:$hex_dir/$hdate/figures/gifs $tmp_dir
         scp -r -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hdate/figures/gifs $tmp_dir

         # copy to FTP
         chmod o+r $tmp_dir/gifs/*
         cp $tmp_dir/gifs/* $FTP
      fi
      #################################################################
   else
      echo "Product already on johansen" >> $cplog
      echo " " >> $cplog
   fi
done

if [ "$wrn_count" -gt 0 ] && [ $thour -gt $T1 ]
then
   echo ""                          >> $cplog
   mail -s "WARNING - Johansen Missing <<$FCtype_hex>> Product(s)" $email < $cplog
fi

# make key with:
# ssh-keygen -t dsa $HOME/.ssh/$keyname
# copy $HOME/.ssh/$keyname.pub to hexagon
# paste contents of $HOME/.ssh/$keyname to the bottom of ~/.ssh/authorised_keys (on target - hexagon)
