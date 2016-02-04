#!/bin/bash
#Get latest restart from the internal repo to the working dir

print_info=0 # print info for debugging

# ====================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC
# if [ $print_info -eq 1 ]
# then
#    echo ""
#    cat $THIS_SRC
#    echo ""
# fi
# ====================================================================================

# ====================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ====================================================================================

# DIRECTORIES AND DATELIST
datelist=$THISFC/logs/datelist.txt
rdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts      # directory with restarts
ddir=$TP4_REALTIME/expt_01.$Xno/data # location of TP4a0.12 directory (where forecast will be done)

# info dir
tday=$(cat $datelist | sed '1!d')
tyear=`date --date=$tday "+%Y"`
pyear=$(expr $tyear - 1)            # previous year
tday_j=`date --date=$tday "+%j"`
tday_j=$((tday_j-1))                # current day julian (1 Jan = 0)
idir=$THISFC2/$tday/info
mkdir -p $idir

# TEXTFILE AND LOG
logdir=$THISFC/logs
mkdir -p $logdir
out_restart=$logdir/last_restart.txt
log=$logdir/tp_get_log.txt
touch $log
echo $tday  >> $log
echo ""     >> $log

#####################################################################
# first look in data directory for last Monday's restart
# - this avoids migrate if possible
lastMon=`date --date="last Monday" "+%Y%m%d"`
if [ $lastMon -ge $tday ]
then
   lastMon=`date --date="$lastMon -7days" "+%Y%m%d"`
fi
lastMon_y=`date --date="$lastMon" "+%Y"`
lastMon_j=`date --date="$lastMon" "+%j"`
lastMon_j=$((lastMon_j-1))
rfil0=TP4restart${lastMon_y}_${lastMon_j}_00

if [[ -f $ddir/$rfil0.a && -f $ddir/$rfil0.b && -f $ddir/${rfil0}ICE.uf ]]
then
   echo Found the latest restart in data: $rfil0
   echo  $rfil0                > $out_restart
   cp    $out_restart            $idir
   #
   echo "Restart already in $ddir"  >> $log
   echo $afil                       >> $log
   echo $bfil                       >> $log
   echo $ufil                       >> $log
   echo ""                          >> $log
   exit
fi
#####################################################################

#####################################################################
# next check migrate

if [ -d $rdir ]
then
   # use this as a test for if migrate is down or not
   if [ $print_info -eq 1 ]
   then
      echo migrate OK
   fi

   f="unassigned"

   if [ -f $rdir/$lastMon_y/$rfil0.tar.gz ]
   then
      # quick check for last Monday's file
      if [ $print_info -eq 1 ]
      then
         echo " "
         echo "Last Monday's restart file is on /migrate"
      fi
      f=$rfil0
      F=$rdir/$lastMon_y/$f

   elif [ -d $rdir/$tyear ]
   then
      # find latest file from current year on migrate
      echo "Check current year's restart files"
      
      mdir=$rdir/$tyear
      cd $mdir

      rlist=(`ls *.gz`)
      nfil=`echo ${#rlist[@]}`
      if [ $nfil -gt 0 ]
      then
         f=${rlist[ $((nfil-1)) ]}
         F=$mdir/$f
      else
         echo "No restart files list found in $mdir"                       >> $log
         echo "Either check topaz_archive_restart or topaz_get_restart"    >> $log
      fi

   elif [ -d $rdir/$pyear ]
   then
      # find latest file from previous year on migrate
      echo "check previous year's restart files"

      mdir=$rdir/$pyear
      cd $mdir

      rlist=(`ls *.gz`)
      nfil=`echo ${#rlist[@]}`
      if [ $nfil -gt 0 ]
      then
         f=${rlist[ $((nfil-1)) ]}
         F=$mdir/$f
      else
         echo "No restart files list found in $mdir"                       >> $log
         echo "Either check topaz_archive_restart or topaz_get_restart"    >> $log
      fi
   fi

   if [ $f == 'unassigned' ]
   then
      # nothing found on migrate
      if [ $print_info -eq 1 ]
      then
         echo "No recent restarts"
         echo "Check topaz_archive_restart and topaz_get_restart ASAP"
      fi

      echo "No recent restarts"                                      >> $log
      echo "Check topaz_archive_restart and topaz_get_restart ASAP"  >> $log
      mail -s "WARNING - Restarts not found on /migrate" $email < $log
      exit
   else

      # check latest restart is not too old
      base=${f%%.tar.gz}
      baseday=${base#*TP4restart} #YYYY_JJJ_HH
      byr=${baseday:0:4}
      bdj=${baseday:5:3}

      if [ $byr -eq $tyear ]
      then
         J2=$(( 10#$tday_j ))                # convert to decimal
         J1=$(( 10#$bdj    ))                # convert to decimal
         age=$((J2-J1))
      else
         J2=$(( 10#$tday_j ))                # convert to decimal
         J0=`date --date="${byr}1231 +%j"`   # days in previous year
         J2=$((tday_j+J0))
         J1=$(( 10#$bdj    ))                # convert to decimal
         age=$((J2-J1))
      fi

      if [ $age -gt 13 ]
      then
         echo "Latest restart file on /migrate:"                  >> $log
         echo "$f"                                                >> $log
         echo "$age days old"                                     >> $log
         mail -s "WARNING - Restarts on /migrate too old" $email   < $log
      elif [ $print_info -eq 1 ]
      then
         echo "Age of restart $F"
         echo "$age - OK"
         echo " "
      fi
   fi

else
   # migrate is down
   # - check if any in data and that they are not too old
   if [ $print_info -eq 1 ]
   then
      echo "Migrate down" >> $log
   fi
   mail -s "WARNING - Migrate down" $email < $log

   cd $ddir
   rlist=(`ls TP4restart*.a`)
   nfil=`echo ${#rlist[@]}`
   if [ $nfil -gt 0 ]
   then
      f=$ddir/${rlist[ $((nfil-1)) ]}
   else
      echo "No restart files list found in $ddir"                       >> $log
      mail -s "Migrate down & no restarts in DATA" $email < $log
      exit
   fi

   # check latest restart is not too old
   base=${f%%.a}
   baseday=${base#TP4restart} #YYYY_JJJ
   byr=${baseday:0:4}
   bdj=${baseday:5:3}
   if [ $byr -eq $tyear ]
   then
      J2=$(( 10#$tday_j ))                # convert to decimal
      J1=$(( 10#$bdj    ))                # convert to decimal
      age=$((J2-J1))
   else
      J2=$(( 10#$tday_j ))                # convert to decimal
      J0=`date --date="${byr}1231 +%j"`   # days in previous year
      J2=$((tday_j+J0))
      J1=$(( 10#$bdj    ))                # convert to decimal
      age=$((J2-J1))
   fi

   if [ $age -gt 13 ]
   then
      echo "Latest restart file in $ddir:"                     >> $log
      echo "$f"                                                >> $log
      echo "$age days old"                                     >> $log
      mail -s "Migrate down & restarts in DATA too old" $email  < $log
      exit
   fi

   # if script comes here,
   # there are some recent-enough files in data
   # - make last_restart.txt
   echo  $base > $out_restart
   cp    $out_restart $idir
   #
   if [ $print_info -eq 1 ]
   then
      echo "Restart already in $ddir"
      echo $afil
      echo $bfil
      echo $ufil
   fi

   echo "Restart already in $ddir"  >> $log
   echo $afil                       >> $log
   echo $bfil                       >> $log
   echo $ufil                       >> $log
   echo ""                          >> $log

   exit
fi

# if script has reached this point,
# there is a restart file on migrate that needs to be retrieved

echo "Most recent restart: $f"
echo "Most recent restart: $f"                                    >> $log
echo "Unpacking..."                                               >> $log
echo " "                                                          >> $log
echo $f > $out_restart
cp    $out_restart $idir

# CREATING DAILY INFO DIR

# $f is now most recent restart
ryear=${f:10:4}	  # year of restart
rday_j=${f:15:3}  # julian day of restart
hr=${f:19:2}	  # hour of restart

# names of restart files we will finally use
f0=TP4restart${ryear}_${rday_j}_${hr}

# original names
afil0=${f0}_mem001.a
bfil0=${f0}_mem001.b
ufil0=${f0}ICE.uf

# final names
afil=${f0}.a
bfil=${f0}.b


# go to the data in $xdir
# if the most recent restart is NOT already there, then:
# - unpack restart there
# - rename files

cd $ddir
if [ 1 -eq 0 ]
then
   # don't delete old restarts
   # - they'll be deleted eventually
   # but will need them for previous week's hindcast
   rm -f TP4restart*
fi
cp $rdir/${ryear}/$f0.tar.gz .
tar -zxvf $f0.tar.gz
rm $f0.tar.gz
mv $afil0 $afil
mv $bfil0 $bfil

echo " "                      >> $log
echo "in $ddir:"              >> $log
echo mv $afil0 $afil          >> $log
echo mv $bfil0 $bfil          >> $log
echo "Don't rename $ufil0"    >> $log
echo " "                      >> $log                                              
