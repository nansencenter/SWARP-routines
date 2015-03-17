#!/bin/bash
#This is a script to both move and compress the TP4restart files.

#17-03-15 VERSION: creating a TP4restart LIST that will register the name of every restart already archived
#I thought that reading through a text file would be faster than checking the whole work folder

rdir=/work/fanf/TOPAZ_RT #restart dir
cyear=`date -ju "+%Y"` # current year
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/ #backup dir
mkdir -p $bdir/$cyear # create backup directory for current year if it doesn't exist already

ldir=$bdir/$cyear/log # where to put the log file and file list (one for each year)
mkdir -p $ldir # create the log dir if it doesn't exist already
TP4rlog=$ldir/TP4rlog #keeping a log (OK if this is for current year not restart year)
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
                base=${afil:0:21}
                bfil=$rdir/${base}_mem001.b
                ufil=$rdir/${base}ICE.uf
                #
                ryear=${base:10:4} # want to distinguish current year from year of restart file
                TP4rlist=$bdir/log/$ryear/TP4rlist #path to the LIST

		if grep -Fxq "$el" $TP4rlist #checking their presence in the LIST
		then
			echo "$el CHECKED" >> $TP4rlog #file present -> check confirmation (NO NEED TO PUT IN THE LOG, TEMP SAFETY)
		else 
			echo "$el" >> $TP4rlist #file not present -> update the LIST
                        echo "$bfil" >> $TP4rlist
                        echo "$ufil" >> $TP4rlist
                        #
			baselist="$baselist $base" #adding the file group to the transfer variable
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
        for base in $baselist
        do
           afil=$rdir/${base}_mem001.a
           bfil=$rdir/${base}_mem001.b
           ufil=$rdir/${base}ICE.uf
           ryear=${base:10:4} # want to distinguish current year from year of restart file
           tfil=${bdir}/$ryear/$base.tar.gz # tar file to create
           tar -zcvf $tfil $afil $bfil $ufil # TODO remove /work/fanf/TOPAZ_RT from inside tar file
           #
           echo "$afil -ADDED" $TP4rlog
           echo "$bfil -ADDED" $TP4rlog
           echo "$ufil -ADDED" $TP4rlog
        done
	# tar -zcvf ${bdir}/TP4restart_$(date +%Y)_$(date --date yesterday +%j).tar.gz $arch #creating the tar file (HYCOM julian-1 day)
	# arch=( $arch) #defining as an array
	# for el in "${arch[@]}"
	# 	do
	# 	echo "$el - ADDED" >> $TP4rlog
	# done
	echo "FILES ADDED - ARCHIVE UP TO DATE" >> $TP4rlog
fi
cat $TP4rlog #printing the log
