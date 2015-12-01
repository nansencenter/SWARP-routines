#!/bin/bash
# copy forecast file from hexagon to johansen
# NB use either this one OR copy_wavesice_ww3a.sh!!

# don't bother doing anything before 6am
thour=`date +%H`
T0=6
T1=11
if [ $thour -lt $T0 ]
then
   exit
fi

# ==================================================================================
# general variables
THREDDS="/mnt/10.11.12.232/Projects/SWARP/Thredds_data_dir" # final destination of forecast outputs
SRjoh="/Data/sim/tim/Projects/SWARP/SWARP-routines"         # location of SWARP routines on johansen
TP4rt="/work/timill/RealTime_Models/results/TP4a0.12"       # location of forecast outputs on hexagon
hidden=$SRjoh/forecast_scripts/hidden                       # location of private info eg logs, ssh_info.src

# FC-specific variables
FCtype_joh="wavesice"               # type of FC on johansen
FCtype_hex="wavesice"               # type of FC on johansen ("wavesice" or "wavesice_ww3arctic")
FC_OUTPUT="SWARPwavesice_forecast"  # start of netcdf file

joh_dir=$THREDDS/$FCtype_joh  #directory on johansen where final file should go
tmp_dir=$joh_dir/tmp          #directory on johansen where final file + figures should go initially after copying
hex_dir=$TP4rt/$FCtype_hex    #directory on hexagon where final file ends up
# ==================================================================================

mkdir -p $tmp_dir

# get $keyname, $user
source $hidden/ssh_info.src


# ==================================================================================
# EMAIL ADRESS
email=$(cat $FCemail)
# ==================================================================================

cplog=$hidden/cplog
tday=$(date +%Y%m%d)

if [ -f $cplog ] && [ $tday != $(cat $cplog | sed '1!d') ]
then
   mv $cplog $hidden/old_logs/$(cat $cplog | sed '1!d')
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
   hex_fil=${FC_output}_start${hdate}*  # final file

   # only do scp if file not present
   if [ ! -f $joh_dir/$hex_fil ]
   then
      echo scp -i $HOME/.ssh/\$keyname \$user@hexagon.bccs.uib.no:$hex_dir/$hdate/final_output/$hex_fil $tmp_dir
      scp -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hdate/final_output/$hex_fil $tmp_dir

      if [ -f $tmp_dir/$hex_fil ]
      then
         # if scp worked move it to THREDDS dir
         mv $tmp_dir/$hex_fil $joh_dir/
         chmod o+r $joh_dir/$hex_fil
         echo "Product found on $hdate!" >> $cplog
         echo "" >> $cplog

         # # =============================================================
         # # NB THIS ISN'T SET UP YET
         # if [ $n -eq 0 ]
         # then
         #    # copy gifs from latest forecast
         #    # TODO: mv to /Webdata and link to website
         #    echo scp -r -i $HOME/.ssh/\$keyname \$user@hexagon.bccs.uib.no:$hex_dir/$hdate/figures/gifs $tmp_dir
         #    scp -r -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hdate/figures/gifs $tmp_dir
         # fi
         # # =============================================================

      else
         # if scp didn't work, give warning
         echo "No product on $hdate" >> $cplog
         wrn_count=$(expr $wrn_count + 1)
      fi
   else
      echo "Product already on johansen" >> $cplog
      echo " " >> $cplog
   fi
done

if [ "$wrn_count" -gt 0 && $thour -gt $T1 ]
then
   echo ""                          >> $cplog
   mail -s "WARNING - Johansen Missing <<$FCtype_hex>> Product(s)" $email < $cplog
fi

# make key with:
# ssh-keygen -t dsa $HOME/.ssh/$keyname
# copy $HOME/.ssh/$keyname.pub to hexagon
# paste contents of $HOME/.ssh/$keyname to the bottom of ~/.ssh/authorised_keys (on target - hexagon)
