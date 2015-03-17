#!/bin/bash
#This is a script to both move and compress the TP4restart files.

#NEW VERSION: creating a TP4restart LIST that will register the name of every restart already archived
#I thought that reading through a text file would be faster than checking the whole work folder
clear

TP4rlog=$HOME/giacomo/restart_dir/TP4rlog #keeping a log, this can be manually cleared once in a while (monthly?) or by crontab
TP4rlist=$HOME/giacomo/restart_dir/TP4rlist #path to the LIST
echo "NEW OPERATION" >> $TP4rlog
date >> $TP4rlog #signing the time of the operation

rdir=/work/fanf/TOPAZ_RT #restart dir
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/2015 #backup dir

pckg=TP4restart* #looking for restarts
eye=`find ${rdir} -name $pckg` #gathering all the restarts
eye=( $eye ) #string into array
arch="" #empty archive variable
for el in "${eye[@]}" #analizing one by one the restart files
	do
	if [ -f $el ] #checking their existence
	then
		if grep -Fxq "$el" $TP4rlist #checking their presence in the LIST
		then
			echo "$el CHECKED" >> $TP4rlog #file present -> check confirmation (NO NEED TO PUT IN THE LOG, TEMP SAFETY)
		else 
			echo "$el" >> $TP4rlist #file not present -> update the LIST
			arch="$arch $el" #adding the file to the transfer variable
		fi
	else
		echo "MISSING RESTART FILES" >> $TP4rlog
		exit
	fi
done
if [ -z "$arch" ] #if arch is empty there are no new restart files
then
	echo "NO MODIFICATIONS - ARCHIVE UP TO DATE" >> $TP4rlog
else
	tar -zcvf ${bdir}/TP4restart_$(date +%Y)_$(date --date yesterday +%j).tar.gz $arch #creating the tar file (HYCOM julian-1 day)
	arch=( $arch) #defining as an array
	for el in "${arch[@]}"
		do
		echo "$el - ADDED" >> $TP4rlog
	done
	echo "FILES ADDED - ARCHIVE UP TO DATE" >> $TP4rlog
fi
cat $TP4rlog #printing the log
