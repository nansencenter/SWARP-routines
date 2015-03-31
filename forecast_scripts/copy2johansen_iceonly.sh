# copy forecast file from hexagon to johansen
# NB template only - need to set $keyname and $user
joh_dir=/mnt/10.11.12.232/Projects/SWARP/Thredds_data_dir/ice_only #directory on johansen where final file should go
rdy_dir=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/website #directory on hexagon with the readytogo file
tmp_dir=/mnt/10.11.12.232/Projects/SWARP/Thredds_data_dir/ice_only/tmp #directory on johansen where final file should go
hex_dir=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/work # dir on hexagon where final file ends up
mkdir -p $tmp_dir

rm cplog.txt
touch cplog.txt

# cheching if everything on hexagon worked fine and if the file is ready
# let's start with a check every 20 mins
for tm in {0..10}
do
   cat > sftp.in<<EOF
cd $rdy_dir
get websitego $tmp_dir
rm websitego
bye
EOF
sftp -b sftp.in -oIdentityFile=$HOME/.ssh/$keyname $user@hexagon.bccs.uib.no
rm sftp.in
   if [ ! -f $tmp_dir/websitego ]
   then
      echo "Try number $tm (`date +%H:%M`) - files not ready" >> cplog.txt
      sleep 1200
   else
      rm $tmp_dir/websitego
      echo "Try number $tm (`date +%H:%M`)- files ready" >> cplog.txt
      # finding the latest final product
      for f in {0..100}
      do
       hdate=`date --date="$f days ago" '+%Y%m%d'`
       echo "Looking for product $hdate" >> cplog.txt
       hex_fil=SWARPiceonly_forecast_start*T000000Z.nc  # final file
       scp -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hdate/final_output/$hex_fil $tmp_dir
       if [ -f $tmp_dir/$hex_fil ]
       then
          mv $tmp_dir/* $joh_dir/
          echo "Product found on $hdate!" >> cplog.txt
          echo "" >> cplog.txt
          echo "Latest product uploaded" >> cplog.txt
          exit
       else
          echo "No product on $hdate" >> cplog.txt
       fi
    done
    echo "!WARNING! NO PRODUCT FOUND" >> cplog.txt
    exit
   fi
done

# make key with:
# ssh-keygen -t dsa $HOME/.ssh/$keyname
# copy $HOME/.ssh/$keyname.pub to hexagon
# paste contents of $HOME/.ssh/$keyname to the bottom of ~/.ssh/authorised_keys


