#!/bin/bash
#This is a script to both move and compress the TP4restart files.

#18-03-15 VERSION: creating a TP4restart LIST that will register the name of every restart already archived
#I thought that reading through a text file would be faster than checking the whole work folder

source $SWARP_ROUTINES/source_files/hex_vars.src
mkdir -p /work/$USER/tmp
cd       /work/$USER/tmp

# =========================================================================================================
# EMAIL ADDRESS FOR THE WEEKLY UPDATE
email=$(cat $FCemail)
# =========================================================================================================

# DIRECTORIES AND TIME DEFINITION
rdir=/work/fanf/TOPAZ_RT 
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts 
fcldir=$FORECAST/logs
cyear=$(date +%Y) 
pyear=$(expr $cyear - 1)
tday=$(date +%Y%m%d-%A)
tdd=$(date +%d)
jday=$(date +%j)
if [ "$(date +%A)" == "Monday" ]
then
   #latest restart is today
   Rlatest=TP4restart${cyear}_$((jday-1))_00
else
   #latest restart is last Monday
   dt=`date --date="last Monday" "+%Y%m%d"`
   dty=`date --date="last Monday" "+%Y"`
   dtj=`date --date="last Monday" "+%j"`
   Rlatest=TP4restart${dty}_$((dtj-1))_00
fi

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

baselist=""                                                          # empty archive variable
have_latest=0
for el in $rdir/TP4restart*.a
   do
   afil=$(basename $el)                                 # *.a filename without full path
   if [ $afil==${Rlatest}_mem001.a ]
   then
      have_latest=1
   fi

   base=${afil%%_mem001.a}                              # cutting _mem001.a
   dcre=${base#TP4restart}                              # cutting TP4restart
   bfil=$rdir/${base}_mem001.b
   ufil=$rdir/${base}ICE.uf
   ryear=${dcre%%_*}                                    # keeping the year, want to distinguish between ops year and file year
   
   if grep -Fxq "$base" $TP4rlist                       # checking their presence in the list
   then
      echo "$base checked" >> $TP4rlog             # file present -> check confirmation
   else 
      echo "$base" >> $TP4rlist                    # file not present -> update the list                      
      sort $TP4rlist -o $TP4rlist                  # sort the list
      baselist="$baselist $base"                   # adding the file group to the transfer variable
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
      dcre=${base#TP4restart}                         # cutting TP4restart
      ryear=${dcre%%_*}                               # keeping the year, we want to distinguish between ops year and file year
      tfil=$base.tar.gz                               # tar file to create
      
      echo " "
      echo "creating $tfil"
      tar -zcvf $tfil -C ${rdir} $afil $bfil $ufil    # remove /work/fanf/TOPAZ_RT from inside tar file
      echo " "

      cp $tfil ${bdir}/$ryear
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

# =======================================================================
# WARNINGS
# 1. latest not in Francois's dir
if [ $have_latest -eq 0 ]
then
   tf=tmp.txt
   echo "Latest restart not present in $rdir" > $tf
   echo "$Rlatest" > $tf
   mail -s "WARNING: Latest restart not present in $rdir" $email < $tf
   rm $tf
# else
#    # test message
#    tf=tmp.txt
#    echo "Latest restart present in $rdir" > $tf
#    echo "$Rlatest" > $tf
#    mail -s "CONFIRM: Latest restart present in $rdir" $email < $tf
#    rm $tf
fi

# 2. latest not in migrate
if [ ! -f $bdir/$dty/$Rlatest.tar.gz ]
then
   tf=tmp.txt
   echo "Latest restart not present in $bdir/$dty" > $tf
   echo "$Rlatest" > $tf
   mail -s "WARNING: Latest restart not on /migrate" $email < $tf
   rm $tf
elif [ "$(date +%A)" == "Monday" ]
then
   # confirmation message
   tf=tmp.txt
   echo "Latest restart present in $bdir/$dty" > $tf
   echo "$Rlatest" > $tf
   mail -s "CONFIRM: Latest restart on /migrate" $email < $tf
   rm $tf
fi
# =======================================================================
