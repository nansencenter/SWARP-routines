#!/bin/bash
#Get latest restart from the internal repo to the working dir

# EMAIL ADDRESS
address=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/fc_alert_email.txt
# ====================================================================================
email=$(cat $address)
# ====================================================================================

# DIRECTORIES AND DATELIST
fcdir=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts
datelist=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/datelist.txt
rdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts      # directory with restarts
ddir=/work/timill/RealTime_Models/TP4a0.12/expt_01.1/data # location of TP4a0.12 directory (where forecast will be done)
logdir=$fcdir/logs
mkdir -p $logdir

# TEXTFILE AND LOG
out_restart=$fcdir/last_restart.txt
log=$logdir/tp_get_log.txt

if [ $(date +%A) == "Monday" ]
then
   mail -s "Weekly tp_get_restart log" $email < $log
   rm $log
fi

touch $log
echo $date  >> $log
echo ""     >> $log

cyear=$(cat $datelist | sed '3!d')			# current year
cmon=$(cat $datelist | sed '4!d')			# current month
cday=$(cat $datelist | sed '6!d')                       # current day julian (1 Jan = 1)
cday=$(expr $cday - 1)		                	# current day julian (1 Jan = 0)
pyear=$(expr $cyear - 1)		              	# previous year

f="unassigned"
if [ -d $rdir/$cyear ]
then
   lfil=$rdir/$cyear/log/tp_archive_list.txt                       #list of restart files
   if [ -f $lfil ]
   then
      f=`cat $lfil | grep "." | tail -1`
      echo "latest restart: $f"                    >> $log
      ryear=$cyear
   else
      echo "Restart file's list NOT FOUND in $rdir/$cyear"         >> $log
      echo "Either check topaz_archive_restart or topaz_get_restart" >> $log
      mail -s "WARNING - Restart list NOT FOUND" $email < $log
   fi
else
   echo "No restarts in current year ($cyear)"     >> $log
   #search previous year
   if [ -d $rdir/$pyear ]
   then
      lfil=$rdir/$pyear/log/tp_archive_list.txt                    #list of restart files
      if [ -f $lfil ]
      then
         f=`cat $lfil | grep "." | tail -1`
         echo "latest restart: $f"                 >> $log
         ryear=$pyear
      else
         echo "Restart file's list NOT FOUND in $rdir/$pyear"         >> $log
         echo "Either check topaz_archive_restart or topaz_get_restart" >> $log
         mail -s "WARNING - Restart list NOT FOUND" $email < $log
      fi
   fi
fi

if [ $f == 'unassigned' ]
then
   echo "No recent restarts"i                                     >> $log
   echo "Couldn't find any restart's lists"                       >> $log
   echo "Check ASAP topaz_archive_restart and topaz_get_restart"  >> $log
   mail -s "WARNING 2 - Restart list NOT FOUND" $email < $log
   exit
fi

echo "Most recent restart: $f"                                    >> $log
echo "Unpacking..."                                               >> $log
echo " "                                                          >> $log
echo $f > $out_restart

# CREATING DAILY INFO DIR
idir=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/work/$(cat $datelist | sed '1!d')/info
mkdir -p $idir
mv $out_restart $idir/

# $f is now most recent restart
ryear=${f:10:4}	# year of restart
jday=${f:15:3}	# julian day of restart
hr=${f:19:2}	# hour of restart

# names of restart files we will finally use
f0=TP4restart${ryear}_${jday}_${hr}

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
   cp $rdir/${ryear}/$f0.tar.gz .
   tar -zxvf $f0.tar.gz
   mv $afil0 $ddir/$afil
   mv $bfil0 $ddir/$bfil
   mv $ufil0 $ddir/$ufil
   echo " "                                              >> $log
   echo mv $afil0 $ddir/$afil                            >> $log
   echo mv $bfil0 $ddir/$bfil                            >> $log
   echo mv $ufil0 $ddir/$ufil                            >> $log
   echo " "                                              >> $log                                              
else
   echo "Restart files already present in $ddir:"        >> $log
   echo $afil                                            >> $log
   echo $bfil                                            >> $log
   echo $ufil                                            >> $log
   echo " "                                              >> $log
fi

# check this test
# calculate days between current date
# and latest restart - if too old (>2 weeks?) send warning
dsr=$(expr $cday - $jday) #Days Since Restart

if [ $dsr -gt 13 ]
then
   echo "Restarts too old"                               >> $log
   echo "Check topaz_archive_restart.sh"                 >> $log
   echo "This WARNING is in topaz_get_restart.sh"        >> $log
   mail -s "TP4 restarts are too old" $email < $log
fi

cp $log $fcdir/logs/

