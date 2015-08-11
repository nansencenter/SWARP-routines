#!/bin/bash
#This is a script to both move and compress the TP4restart files.

#18-03-15 VERSION: creating a TP4restart LIST that will register the name of every restart already archived
#I thought that reading through a text file would be faster than checking the whole work folder

source $SWARP_ROUTINES/source_files/hex_vars.src

# EMAIL ADDRESS FOR THE WEEKLY UPDATE
email_list=$FORECAST/fc_alert_email.txt
# =========================================================================================================
email=$(cat $email_list)
# =========================================================================================================

# DIRECTORIES AND TIME DEFINITION
rdir=/work/fanf/TOPAZ_RT 
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts 
fcldir=$FORECAST/logs
cyear=$(date +%Y) 
pyear=$(expr $cyear - 1)
tday=$(date +%Y%m%d-%A)
tdd=$(date +%d)

# CREATING THE LOG DIRECTORY AND FILE
ldir=$bdir/$cyear/log 
mkdir -p $ldir
TP4rlog=$ldir/tp_archive_log.txt
TP4rlist=$ldir/tp_archive_list.txt
tmplist=$ldir/tmp_list.txt

if ! [ -f "$TP4rlist" ]
then
   cp $bdir/$pyear/$TP4rlist $ldir/
fi

touch $TP4rlog

# FETCHING OPERATION
echo $tday  >> $TP4rlog
echo ''     >> $TP4rlog


if [ "$(date +%A)" == "Monday" ]
then
   echo ""                    >> $TP4rlog
   echo "List used this past week:   " >> $TP4rlog
   echo ""                    >> $TP4rlog
   cat $TP4rlist              >> $TP4rlog                               # printing the list
   echo ""                    >> $TP4rlog
   echo ""                    >> $TP4rlog
fi

pckg=TP4restart*_mem001.a                                            # looking for restarts (just check *.a files)
eye=`find ${rdir} -name $pckg`                                       # gathering all the restarts
eye=( $eye )                                                         # string into array
baselist=""                                                          # empty archive variable
for el in "${eye[@]}"                                                # analizing one by one the restart files
   do
   if [ -f $el ]                                                # checking their existence
   then
      afil=$(basename $el)                                 # *.a filename without full path
      base=${afil%%_mem001.a}                              # cutting _mem001.a
      dcre=${base#TP4restart}                              # cutting TP4restart
      bfil=$rdir/${base}_mem001.b
      ufil=$rdir/${base}ICE.uf
      ryear=${dcre%%_*}                                    # keeping the year, want to distinguish between ops year and file year
      
      if grep -Fxq "$base" $TP4rlist                       # checking their presence in the LIST
      then
      	echo "$base CHECKED" >> $TP4rlog             # file present -> check confirmation
      else 
      	echo "$base" >> $TP4rlist                    # file not present -> update the LIST                      
                 sort $TP4rlist -o $TP4rlist                  # sort the list
      	baselist="$baselist $base"                   # adding the file group to the transfer variable
      fi
   else
      echo "MISSING RESTART FILES in $rdir" >> $TP4rlog
         echo "Check if the restarts are still uploaded on $rdir" >> $TP4rlog
         mail -s "WARNING - TOPAZ restart files" $email < $TP4rlog
      exit
   fi
done

# ARCHIVING OPERATION
if [ -z "$baselist" ]                                             # if baselist is empty there are no new restart files
then
   echo " NO MODIFICATIONS - ARCHIVE UP TO DATE " >> $TP4rlog
else
   for base in $baselist
   do
      afil=${base}_mem001.a
      bfil=${base}_mem001.b
      ufil=${base}ICE.uf
      dcre=${base#TP4restart}                                   # cutting TP4restart
      ryear=${dcre%%_*}                                         # keeping the year, we want to distinguish between ops year and file year
      tfil=${bdir}/$ryear/$base.tar.gz                          # tar file to create
      tar -zcvf $tfil -C ${rdir} $afil $bfil $ufil              # remove /work/fanf/TOPAZ_RT from inside tar file
      echo "$base -ADDED" >> $TP4rlog
   done
   echo "FILES ADDED - ARCHIVE UP TO DATE" >> $TP4rlog
fi

cp $TP4rlog $fcldir/

if [ "$(date +%A)" == "Monday" ]
then
   nol=$(cat $TP4rlist | wc -l)
   if [ $nol -gt 4 ]
   then
      tbr=$(cat $TP4rlist | sed '1!d')
      touch $tmplist
      echo "The following file will be removed from the list:  "  >> $TP4rlog
      echo $tbr                                                   >> $TP4rlog
      rm /work/timill/RealTime_Models/TP4a0.12/expt_01.1/data/${tbr}*
      sed '1d' $TP4rlist >> $tmplist
      mv $tmplist $TP4rlist
   fi
   weekn=$(expr $(date +%d) / 7)
   mail -s "Week $weekn - topaz_archive LOG" $email < $TP4rlog
   rm $TP4rlog 
fi

