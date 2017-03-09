#!/bin/bash
# copy forecast file from hexagon to johansen

print_info=0
thour=`date +%H`

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

   T0=2  # don't bother doing anything before 2am
   T1=9  # email alert after this time
   if [ $thour -lt $T0 ]
   then
      exit
   fi

elif [ $FCtype_hex == wavesice ]
then
   FCtype_joh="wavesice"               # type of FC on johansen
   FC_OUTPUT="SWARPwavesice_forecast"  # start of netcdf file

   # can rename here
   FC_OUTPUT2="SWARPwavesice_forecast"

   T0=6  # don't bother doing anything before 6am
   T1=11 # email alert after this time
   if [ $thour -lt $T0 ]
   then
      exit
   fi

elif [ $FCtype_hex == wavesice_ww3arctic ]
then
   FCtype_joh="wavesice"                     # type of FC on johansen
   FC_OUTPUT="SWARPwavesice_WW3_forecast"    # start of netcdf file

   # can rename here
   FC_OUTPUT2="SWARPwavesice_forecast"

   T0=6  # don't bother doing anything before 6am
   T1=11 # email alert after this time
   if [ $thour -lt $T0 ]
   then
      exit
   fi

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

# FTP directory
FTP="/Data/FTPRoot/pub/Tim/SWARP_forecasts/$reg/"

# SWARP web page
WEB="/WebData/swarp.nersc.no/SWARP_forecasts/$reg/"
FCemail=$hidden/FCemail.txt

# FC-specific variables
joh_dir=$THREDDS/$FCtype_joh  # directory on johansen where final file should go
tmp_dir=$joh_dir/tmp          # directory on johansen where final file + figures should go initially after copying
hex_dir=$TP4rt/$FCtype_hex    # directory on hexagon where final file ends up
# ==================================================================================

mkdir -p $tmp_dir

# get $keyname, $user
source $hidden/ssh_info.src

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
if [ $thour -eq 23 ]
then
   # only check previous days once a day
   # (plus cleaning of old FC's to save space)
   Nback=30
else
   Nback=0
fi
wrn_count=0
for n in `seq 0 $Nback`
do
   echo "$(date +%H:%M) - looking for <<$FCtype_hex>> latest file" >> $cplog
   hdate=$(date --date="$n days ago" '+%Y%m%d')
   echo "Looking for product $hdate" >> $cplog
   hex_fil=${FC_OUTPUT}_start${hdate}T000000Z.nc   # final file on hexagon
   joh_fil=${FC_OUTPUT2}_start${hdate}T000000Z.nc  # final file on johansen

   if [ $n -gt 7 ]
   then
      # clean older files to save space
      rm -f $joh_dir/$joh_fil
      continue
   fi

   # only try scp if file not present
   if [ ! -f $joh_dir/$joh_fil ]
   then
      # 1st check hexagon
      if [ ! "$(ssh -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no ls -A $hex_dir/$hdate/final_output/*.nc 2>/dev/null)" == "" ]
      then
         rm -f $tmp_dir/gifs/*
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
      elif [ $print_info -eq 1 ]
      then
         echo "Forecast file is not on hexagon"
      fi
   else
      echo "Netcdf product already on johansen" >> $cplog
      echo " " >> $cplog
   fi


   #################################################################
   # gif's for website
   # copy gif files and make comb figures if:
   # 1) the gif folder is cleared (when the .nc file is copied, see above) and
   # 2) the gif files exist on hexagon for the current day
   # if [ ! "$(ls -A $tmp_dir/gifs/*.gif )" ] && [ ! "$(ssh -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no ls -A $hex_dir/$hdate/figures/gifs/*.gif 2>/dev/null)" == "" ];

   # check for gifs' presence
   lst=($tmp_dir/gifs/*.gif)
   Ng=${#lst[@]}
   if [ $Ng -eq 1 ]
   then
      if [ ! -f ${lst[0]} ]
      then
         # change from 1 to 0
         Ng=0
      fi
   fi
   
   # try to copy if they are not here
   if [ $Ng -eq 0 ]
   then
      if [ $n -eq 0 ]
      then
         if [ ! "$(ssh -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no ls -A $hex_dir/$hdate/figures/gifs/*.gif 2>/dev/null)" == "" ]
         then
            # copy gifs from latest forecast
            # mv to /Webdata and link to website
            echo scp -r -i $HOME/.ssh/\$keyname \$user@hexagon.bccs.uib.no:$hex_dir/$hdate/figures/gifs $tmp_dir
            scp -r -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hdate/figures/gifs $tmp_dir

            # need to change permissions to make files accessible on FTP or SWARP area
            chmod o+r $tmp_dir/gifs/*
            chmod o+w $tmp_dir/gifs/*
            chmod g+r $tmp_dir/gifs/*
            chmod g+w $tmp_dir/gifs/*

            # copy to FTP
            cp $tmp_dir/gifs/* $FTP

            # make combined figures synced in time for WEB page
            if [ $FCtype_joh == "ice_only" ]; then
               # sea ice concentration and thickness
               convert \( $tmp_dir/gifs/icec.gif  -coalesce -append \) \
                       \( $tmp_dir/gifs/icetk.gif -coalesce -append \) \
                       +append -crop x600 +repage -set delay 15 -loop 0 $tmp_dir/gifs/IO1comb.gif
               # Sea ice and surface velocity speed
               convert \( $tmp_dir/gifs/uice.gif  -coalesce -append \) \
                       \( $tmp_dir/gifs/usurf.gif -coalesce -append \) \
                       +append -crop x600 +repage -set delay 15 -loop 0 $tmp_dir/gifs/IO2comb.gif
            fi
              
            if [ $FCtype_joh == "wavesice" ]; then
               # Max floe size and wave-in-ice significant wave height
               convert \( $tmp_dir/gifs/dmax.gif  -coalesce -append \) \
                       \( $tmp_dir/gifs/swh.gif -coalesce -append \) \
                       +append -crop x600 +repage -set delay 15 -loop 0 $tmp_dir/gifs/WIcomb.gif
            fi

            # need to change permissions to make files accessible on FTP or SWARP area
            chmod o+r $tmp_dir/gifs/*comb.gif
            chmod o+w $tmp_dir/gifs/*comb.gif
            chmod g+r $tmp_dir/gifs/*comb.gif
            chmod g+w $tmp_dir/gifs/*comb.gif

            # copy to SWARP web page area
            cp $tmp_dir/gifs/* $WEB
         elif [ $print_info -eq 1 ]
         then
            echo "No gifs on hexagon"
         fi
      fi
      #################################################################
#   else
#      echo "Product already on johansen" >> $cplog
#      echo " " >> $cplog
   fi
done

if [ "$wrn_count" -gt 0 ] && [ $thour -gt $T1 ]
then
   echo ""                                                                       >> $cplog
   mail -s "WARNING - Johansen Missing <<$FCtype_hex/$reg>> Product(s)" $email   <  $cplog
fi

# make key with:
# ssh-keygen -t dsa $HOME/.ssh/$keyname
# copy $HOME/.ssh/$keyname.pub to hexagon
# paste contents of $HOME/.ssh/$keyname to the bottom of ~/.ssh/authorised_keys (on target - hexagon)
