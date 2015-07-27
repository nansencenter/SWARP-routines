#!/bin/bash
# copy forecast file from hexagon to johansen

# don't bother doing anything before 2am
thour=`date +%H`
if [ $thour -lt 2 ]
then
   exit
fi

# final destination of forecast outputs
THREDDS="/mnt/10.11.12.232/Projects/SWARP/Thredds_data_dir"
# location of SWARP routines on johansen
SRjoh="/Data/sim/tim/Projects/SWARP/SWARP-routines"
# location of SWARP routines on hexagon
SRhex="/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines"
# location of forecast outputs on hexagon
TP4rt="/work/timill/RealTime_Models/results/TP4a0.12"

joh_dir=$THREDDS/ice_only #directory on johansen where final file should go
joh_dir_w=$THREDDS/wavesice #directory on johansen where final file should go
tmp_dir=$THREDDS/ice_only/tmp #directory on johansen where final file should go
tmp_dir_w=$THREDDS/wavesice/tmp #directory on johansen where final file should go
hex_dir=$TP4rt/ice_only/work # dir on hexagon where final file ends up
hex_dir_w=$TP4rt/wavesice/work # dir on hexagon where final file ends up
wrk_dir=$SRjoh/forecast_scripts/hidden # location of logs, ssh_info.src
# fc_dir=$SRhex/forecast_scripts

mkdir -p $tmp_dir
mkdir -p $tmp_dir_w

# get $keyname, $user
source $wrk_dir/ssh_info.src

# scp -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$fc_dir/fc_alert_email.txt $tmp_dir
cp $wrk_dir/../fc_alert_email.txt $tmp_dir

# EMAIL ADRESS
# ==================================================================================
email=$(cat $tmp_dir/fc_alert_email.txt)
# ==================================================================================

cplog=$wrk_dir/cplog
tday=$(date +%Y%m%d)

if [ -f $cplog ] && [ $tday != $(cat $cplog | sed '1!d') ]
then
   mv $cplog $wrk_dir/old_logs/$(cat $cplog | sed '1!d')
fi

if [ $(date +%A) == "Monday" ]
then
   rm -f $wrk_dir/old_logs/*
fi

if [ ! -f "$cplog" ]
then
   touch $cplog
   echo $tday >> $cplog
fi

# finding the latest final product - ICE_ONLY
# - check last $Nback days
Nback=4
wrn_count=0
for n in `seq 0 $Nback`
do
   echo "$(date +%H:%M) - looking for ICE_ONLY latest file" >> $cplog
   hdate=$(date --date="$n days ago" '+%Y%m%d')
   echo "Looking for product $hdate" >> $cplog
   hex_fil=SWARPiceonly_forecast_start${hdate}*  # final file

   # only do scp if file not present
   if [ ! -f $joh_dir/$hex_fil ]
   then
      echo scp -i $HOME/.ssh/\$keyname \$user@hexagon.bccs.uib.no:$hex_dir/$hdate/final_output/$hex_fil $tmp_dir
      scp -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hdate/final_output/$hex_fil $tmp_dir

      if [ -f $tmp_dir/$hex_fil ]
      then
         # if scp worked move it to THREDDS dir
         mv $tmp_dir/* $joh_dir/
         chmod o+r $joh_dir/$hex_fil
         echo "Product found on $hdate!" >> $cplog
         echo "" >> $cplog
      else
         # if scp didn't work, give warning
         echo "No product on $hdate" >> $cplog
         wrn_count=$(expr $wrn_count + 1)
      fi
   fi
done

if [ "$wrn_count" -gt 0 ]
then
   echo ""                          >> $cplog
   mail -s "WARNING - Johansen Missing ice_only Product(s)" $email < $cplog
fi

thour=`date +%H`
if [ $thour -lt 6 ]
then
   exit
fi
# finding the latest final product - WAVESICE
# - check last $Nback days
# TODO add check to see if file exists already, before downloading
wrn_count=0
for n in `seq 0 $Nback`
do
   echo "$(date +%H:%M) - looking for WAVESICE latest file" >> $cplog
   hdate=$(date --date="$n days ago" '+%Y%m%d')
   echo "Looking for product $hdate" >> $cplog
   hex_fil=SWARPwavesice_forecast_start${hdate}*  # final file

   # only do scp if final product is not present
   if [ ! -f $joh_dir_w/$hex_fil ]
   then
      echo scp -i $HOME/.ssh/\$keyname \$user@hexagon.bccs.uib.no:$hex_dir_w/$hdate/final_output/$hex_fil $tmp_dir_w
      scp -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir_w/$hdate/final_output/$hex_fil $tmp_dir_w

      if [ -f $tmp_dir_w/$hex_fil ]
      then
         # if scp worked move it to THREDDS dir
         mv $tmp_dir_w/* $joh_dir_w/
         chmod o+r $joh_dir_w/$hex_fil
         echo "Product found on $hdate!" >> $cplog
         echo "Latest product uploaded" >> $cplog
         echo "" >> $cplog
      else
         # if scp didn't work, give warning
         echo "No product on $hdate" >> $cplog
         wrn_count=$(expr $wrn_count + 1)
      fi
   fi
done

if [ "$wrn_count" -gt 0 ]
then
   echo ""                          >> $cplog
   mail -s "WARNING - Johansen Missing waves_ice Product(s)" $email < $cplog
fi

# make key with:
# ssh-keygen -t dsa $HOME/.ssh/$keyname
# copy $HOME/.ssh/$keyname.pub to hexagon
# paste contents of $HOME/.ssh/$keyname to the bottom of ~/.ssh/authorised_keys (on target - hexagon)
