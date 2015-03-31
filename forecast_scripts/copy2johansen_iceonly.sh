# copy forecast file from hexagon to johansen
# NB only template - actual file (with $keyname and $user) in "hidden" (not part of repo)
joh_dir=/mnt/10.11.12.232/Projects/SWARP/Thredds_data_dir/ice_only #directory on johansen where final file should go
hex_dir= # dir on hexagon where final file ends up
hex_fil= # final file

# make key with:
# ssh-keygen -t dsa $HOME/.ssh/$keyname
# copy $HOME/.ssh/$keyname.pub to hexagon
# paste contents of $HOME/.ssh/$keyname to the bottom of ~/.ssh/authorised_keys
scp -i $HOME/.ssh/$keyname $user@hexagon.bccs.uib.no:$hex_dir/$hex_fil $joh_dir
