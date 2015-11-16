#!/bin/bash
#Get latest restart from the internal repo to the working dir

print_info=1 # print info for debugging

# ====================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=$1
source $THIS_SRC
if [ $print_info -eq 1 ]
then
   echo ""
   cat $THIS_SRC
   echo ""
fi
# gets:
# $THISFC   - scripts location
# $THISFC2  - outputs go here
# $Xno      - experiment number for running model
# ====================================================================================

# ====================================================================================
# EMAIL ADDRESS
address=$FORECAST/fc_alert_email.txt
email=$(cat $address)
# ====================================================================================

# DIRECTORIES AND DATELIST
datelist=$FORECAST/ice_only/logs/datelist.txt
rdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts      # directory with restarts
ddir=$TP4_REALTIME/expt_01.$Xno/data # location of TP4a0.12 directory (where forecast will be done)

# info dir
tday=$(cat $datelist | sed '1!d')
tyear=`date --date=$tday "+%Y"`
pyear=$(expr $tyear - 1)		              	# previous year
tday_j=$(cat $datelist | sed '6!d')                     # current day julian (1 Jan = 0)
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
   echo found the latest one:    $rfil0
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
         echo "Last Monday's restart file is there"
      fi
      f=$rfil0

   elif [ -d $rdir/$tyear ]
   then
      # find latest file from current year on migrate
      if [ $print_info -eq 1 ]
      then
         echo "check list of current year's restart files"
      fi
      lfil=$rdir/$tyear/log/tp_archive_list.txt
      if [ -f $lfil ]
      then
         f=`cat $lfil | grep "." | tail -1`
         echo "latest restart: $f"                    >> $log
         ryear=$tyear
      else
         echo "Restart file's list NOT FOUND in $rdir/$tyear"         >> $log
         echo "Either check topaz_archive_restart or topaz_get_restart" >> $log
      fi

   else
      # find latest file from previous year on migrate
      if [ $print_info -eq 1 ]
      then
         echo "check list of previous year's restart files"
      fi
      echo "No restarts in current year ($tyear)"     >> $log

      if [ -d $rdir/$pyear ]
      then
         lfil=$rdir/$pyear/log/tp_archive_list.txt                    #list of restart files
         if [ -f $lfil ]
         then
            sort $lfil -o $lfil
            f=$(cat $lfil | tail -1)
            echo "latest restart: $f"                 >> $log
            ryear=$pyear
         else
            echo "Restart file's list NOT FOUND in $rdir/$pyear"         >> $log
            echo "Either check topaz_archive_restart or topaz_get_restart" >> $log
         fi
      fi
   fi

   if [ $f == 'unassigned' ]
   then
      # nothing found on migrate
      if [ $print_info -eq 1 ]
      then
         echo "No recent restarts"
         echo "Couldn't find any restart's lists"
         echo "Check ASAP topaz_archive_restart and topaz_get_restart"
      fi

      echo "No recent restarts"                                      >> $log
      echo "Couldn't find any restart's lists"                       >> $log
      echo "Check ASAP topaz_archive_restart and topaz_get_restart"  >> $log
      mail -s "WARNING - Restart list NOT FOUND" $email < $log
      exit
   fi
else
   # migrate is down
   # - check if any in data and that they are not too old
   if [ $print_info -eq 1 ]
   then
      echo "Migrate down" >> $log
   fi
   mail -s "WARNING - Migrate down" $email < $log

   if [ -f $ddir/TP4restart*.a ]
   then
      for afil in $ddir/TP4restart*.a
      do
         n=$((n+1))
         base=${afil%%.a}
         baseday=${base#TPrestart} #YYYY_JJJ
         byr=${baseday:0:4}
         bdj=${baseday:5:3}
         echo $afil # $base $baseday $byr $bdj
         bdj=$((bdj+1))
         pyr=$((byr-1))
         basedate=`date --date='${pyr}1231' +%Y%m%d`
         bfil=${base}.b
         ufil=${base}ICE.uf

         if [[ -f $bfil && -f $afil && -f $ufil ]]
         then
            # check if .b file and .uf files are also there
            echo "Restarts already in $ddir" >> $log
            echo $afil                       >> $log
            echo $bfil                       >> $log
            echo $ufil                       >> $log
            echo " "                         >> $log

            if [ $n -eq 0 ]
            then
               latest_date=$basedate
               latest_file=$base
            elif [ $rdate -gt $latest_date ]
            then
               latest_date=$basedate
               latest_file=$base
            fi
         fi
      done
   else
      # no restart files in data
      echo ""                       >> $log
      echo "No restarts in $ddir"   >> $log
      mail -s "Migrate down & no restarts in DATA" $email < $log
      exit
   fi

   # if script comes here,
   # there are some files in data
   # - make last_restart.txt
   echo  $latest_file > $out_restart
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

   # give warning if restarts too old
   cutoff=`date -d '$tday - 13days' "+%Y%m%d"`
   if [ $latest_date -lt $cutoff ]
   then
      echo "Restarts too old"                               >> $log
      echo "Check topaz_archive_restart.sh"                 >> $log
      echo "This WARNING is in topaz_get_restart.sh"        >> $log
      mail -s "TP4 restarts are too old" $email < $log
   fi
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
ufil=${f0}ICE.uf


# go to the data in $xdir
# if the most recent restart is NOT already there, then:
# - unpack restart there
# - rename files

if [ ! -f $ddir/$afil ]
then
   cd $ddir
   rm -f TP4restart*
   cp $rdir/${ryear}/$f0.tar.gz .
   tar -zxvf $f0.tar.gz
   rm $f0.tar.gz
   mv $afil0 $ddir/$afil
   mv $bfil0 $ddir/$bfil
   mv $ufil0 $ddir/$ufil
   echo " "                                              >> $log
   echo mv $afil0 $ddir/$afil                            >> $log
   echo mv $bfil0 $ddir/$bfil                            >> $log
   echo mv $ufil0 $ddir/$ufil                            >> $log
   echo " "                                              >> $log                                              
else
   if [ $print_info -eq 1 ]
   then
      echo "Restart files already present in $ddir:"
   fi
   echo "Restart files already present in $ddir:"        >> $log
   echo $afil                                            >> $log
   echo $bfil                                            >> $log
   echo $ufil                                            >> $log
   echo " "                                              >> $log
fi

# check this test
# calculate days between current date
# and latest restart - if too old (>2 weeks?) send warning
dsr=$(expr $tday_j - $rday_j) #Days Since Restart

if [ $dsr -gt 13 ]
then
   echo "Restarts too old"                               >> $log
   echo "Check topaz_archive_restart.sh"                 >> $log
   echo "This WARNING is in topaz_get_restart.sh"        >> $log
   mail -s "TP4 restarts are too old" $email < $log
fi
