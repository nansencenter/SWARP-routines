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

print_info=1

# DIRECTORIES AND TIME DEFINITION
rdir=/work/fanf/TOPAZ_RT 
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts 
fcldir=$FORECAST/common/logs
mkdir -p $fcldir

cyear=$(date +%Y) 
pyear=$(expr $cyear - 1)
tday=$(date +%Y%m%d-%A)
cday=$(date +%Y%m%d)
tdd=$(date +%d)
jday=$(date +%j)
if [ "$(date +%A)" == "Monday" ]
then
   #latest restart is today
   Rlatest=TP4restart${cyear}_$((jday-1))_00
   dt=$cday
   dty=$cyear
   dtj=$jday
else
   #latest restart is last Monday
   dt=`date --date="last Monday" "+%Y%m%d"`
   dty=`date --date="last Monday" "+%Y"`
   dtj=`date --date="last Monday" "+%j"`
   Rlatest=TP4restart${dty}_$((dtj-1))_00
fi

if [ $print_info -eq 1 ]
then
   echo $Rlatest
fi

# CREATING THE LOG DIRECTORY AND FILE
ldir=$bdir/$cyear/log 
mkdir -p $ldir
TP4rlog=$ldir/tp_archive_log.txt
touch $TP4rlog

# FETCHING OPERATION
echo $tday  >> $TP4rlog
echo ''     >> $TP4rlog

baselist=""                               # empty archive variable
have_latest=0
for el in $rdir/TP4restart*.a
do
   afil=$(basename $el)                   # *.a filename without full path
   if [ $print_info -eq 1 ]
   then
      echo $afil
   fi

   if [ $afil==${Rlatest}_mem001.a ]
   then
      have_latest=1
   fi

   base=${afil%%_mem001.a}                # cutting _mem001.a
   dcre=${base#TP4restart}                # cutting TP4restart
   bfil=$rdir/${base}_mem001.b
   ufil=$rdir/${base}ICE.uf
   ryear=${dcre%%_*}                      # keeping the year, want to distinguish between ops year and file year
   
   echo $bdir/$ryear/$base.tar.gz
   if [ ! -f $bdir/$ryear/$base.tar.gz ]
   then
      baselist="$baselist $base"                   # adding the file group to the transfer variable
      if [ $print_info -eq 1 ]
      then
         echo "$base not present on /migrate"
         echo ""
      fi
   else
      if [ $print_info -eq 1 ]
      then
         echo "$base already present on /migrate"
         echo ""
      fi
   fi
done

if [ $print_info -eq 1 ]
then
   echo $baselist
   if [ $have_latest -eq 1 ]
   then
      echo "Latest restart present on $rdir"
      echo "Latest restart: $Rlatest"
   else
      echo "Latest restart not present on $rdir"
      echo "Latest restart: $Rlatest"
   fi
fi

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

      if [ $print_info -eq 1 ]
      then
         echo "cp $tfil ${bdir}/$ryear"
      fi
      cp $tfil ${bdir}/$ryear
      echo "$base ADDED" >> $TP4rlog
   done
fi

cp $TP4rlog $fcldir

if [ "$(date +%A)" == "Monday" ]
then
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
   echo "Latest restart not present in $rdir"                     > $tf
   echo "$Rlatest"                                                > $tf
   mail -s "WARNING: Latest restart not present in $rdir" $email  < $tf
   rm $tf
elif [ $print_info -eq 1 ]
then
   # Confirmation message
   tf=tmp.txt
   echo "Latest restart present in $rdir"                      > $tf
   echo "$Rlatest"                                             > $tf
   mail -s "CONFIRM: Latest restart present in $rdir" $email   < $tf
   rm $tf
fi

# 2. latest not in migrate
if [ ! -f $bdir/$dty/$Rlatest.tar.gz ]
then
   tf=tmp.txt
   echo "Latest restart not present in $bdir/$dty" > $tf
   echo "Latest restart: $Rlatest"                 > $tf
   echo " "                                        > $tf
   echo "Contents of $bdir/$dty:"                  > $tf
   ls -lh $bdir/$dty                               > $tf

   #email
   mail -s "WARNING: Latest restart not on /migrate" $email < $tf
   rm $tf
elif [ "$(date +%A)" == "Monday" ]
then
   # confirmation message
   tf=tmp.txt
   echo "Latest restart present in $bdir/$dty"     > $tf
   echo "Latest restart: $Rlatest"                 > $tf
   echo " "                                        > $tf
   echo "Contents of $bdir/$dty:"                  > $tf
   ls -lh $bdir/$dty                               > $tf

   #email
   mail -s "CONFIRMATION: Latest restart on /migrate" $email < $tf
   rm $tf
fi
# =======================================================================
