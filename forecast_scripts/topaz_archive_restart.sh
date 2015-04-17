#!/bin/bash
#This is a script to both move and compress the TP4restart files.

#18-03-15 VERSION: creating a TP4restart LIST that will register the name of every restart already archived
#I thought that reading through a text file would be faster than checking the whole work folder


rdir=/work/fanf/TOPAZ_RT #restart dir
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts #backup dir
cyear=$(date -u "+%Y") # current day
ldir=$bdir/$cyear/log # where to put the log file and file list (one for each year)

mkdir -p $ldir # create the log dir if it doesn't exist already
TP4rlog=$ldir/TP4rlog #keeping a log (OK if this is for current year not restart year)

rm $TP4rlog
touch $TP4rlog

echo "NEW OPERATION" >> $TP4rlog
date >> $TP4rlog #signing the time of the operation

pckg=TP4restart*_mem001.a #looking for restarts (just check *.a files)
eye=`find ${rdir} -name $pckg` #gathering all the restarts
eye=( $eye ) #string into array
baselist="" #empty archive variable
for el in "${eye[@]}" #analizing one by one the restart files
	do
	if [ -f $el ] #checking their existence
	then
                afil=$(basename $el) # *.a filename without full path
                base=${afil%%_mem001.a} #cutting _mem001.a
		dcre=${base#TP4restart} #cutting TP4restart
                bfil=$rdir/${base}_mem001.b
                ufil=$rdir/${base}ICE.uf
                ryear=${dcre%%_*} #keeping the year, want to distinguish between ops year and file year
                TP4rlist=$bdir/$ryear/log/TP4rlist #path to the LIST

		if grep -Fxq "$base" $TP4rlist #checking their presence in the LIST (we treat the 3 files as a unique set)
		then
			echo "$base CHECKED" >> $TP4rlog #file present -> check confirmation
		else 
			echo "$base" >> $TP4rlist #file not present -> update the LIST                      
			baselist="$baselist $base" #adding the file group to the transfer variable
		fi
	else
		echo "MISSING RESTART FILES" >> $TP4rlog
		exit
	fi
done
if [ -z "$baselist" ] #if baselist is empty there are no new restart files
then
	echo "NO MODIFICATIONS - ARCHIVE UP TO DATE" >> $TP4rlog
else
        for base in $baselist
        do
           afil=${base}_mem001.a
           bfil=${base}_mem001.b
           ufil=${base}ICE.uf
	   dcre=${base#TP4restart} #cutting TP4restart
           ryear=${dcre%%_*} #keeping the year, we want to distinguish between ops year and file year
           tfil=${bdir}/$ryear/$base.tar.gz # tar file to create
           tar -zcvf $tfil -C ${rdir} $afil $bfil $ufil #remove /work/fanf/TOPAZ_RT from inside tar file
           echo "$base set -ADDED" >> $TP4rlog
        done
	echo "FILES ADDED - ARCHIVE UP TO DATE" >> $TP4rlog
fi
echo "Actual List:   "
echo ""
cat $TP4rlist #printing the list
echo "Today's Log:   "
echo ""
cat $TP4rlog #printing the log
