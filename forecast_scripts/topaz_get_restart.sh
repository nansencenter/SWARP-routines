#!/bin/bash
#Get latest restart from the internal repo to the working dir

# textfile output:
out_restart=last_restart.txt

# fetch datelist
datelist=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/datelist
rdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts      # directory with restarts
ddir=/work/timill/RealTime_Models/TP4a0.12/expt_01.1/data # location of TP4a0.12 directory (where forecast will be done)

cyear=$(cat $datelist | sed '3!d')			# current year
cmon=$(cat $datelist | sed '4!d')			# current month
cday=$(cat $datelist | sed '6!d')                       # current day julian (1 Jan = 1)
cday=`expr $cday - 1`		                	# current day julian (1 Jan = 0)
pyear=`expr $cyear - 1`		                	# previous year

f="unassigned"
if [ -d $rdir/$cyear ]
then
   lfil=$rdir/$cyear/log/TP4rlist #list of restart files
   if [ -f $lfil ]
   then
      f=`cat $lfil | grep "." | tail -1`
      echo "latest restart: $f"
      ryear=$cyear
   else
      echo "No restarts in current year ($cyear)"
   fi
else
   echo "No restarts in current year ($cyear)"
fi

if [ f == 'unassigned' ]
then
   #search previous year
   if [ -d $rdir/$pyear ]
   then
      lfil=$rdir/$pyear/log/TP4rlist #list of restart files
      if [ -f $lfil ]
      then
         f=`cat $lfil | grep "." | tail -1`
         echo "latest restart: $f"
         ryear=$pyear
      else
         echo "No restarts in previous year ($pyear)"
      fi
   else
      echo "No restarts in previous year ($pyear)"
   fi
fi

if [ f == 'unassigned' ]
then
   echo "No recent restarts"
   exit
fi

echo "Most recent restart: $f"
echo "Unpacking..."
echo " "
echo $f > $out_restart
idir=/work/timill/RealTime_Models/results/TP4a0.12/ice_only/work/$(cat $datelist | sed '1!d')/info
mkdir -p $idir
mv $out_restart $idir

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
   echo " "
   echo mv $afil0 $ddir/$afil
   echo mv $bfil0 $ddir/$bfil
   echo mv $ufil0 $ddir/$ufil
   echo " "
else
   echo "Restart files already present in $ddir:"
   echo $afil
   echo $bfil
   echo $ufil
   echo " "
fi

# check this test
# calculate days between current date
# and latest restart - if too old (>2 weeks?) send warning
dsr=$(expr $cday - $jday) #Days Since Restart

if [ $dsr -gt 13 ]
then
   touch WARNING
   echo "Restarts too old" >> WARNING
   echo "Check topaz_archive_restart.sh" >> WARNING
   echo "This WARNING is in topaz_get_restart.sh" >> WARNING
   mail -s "RESTARTS WARNING" gcmdnt90@gmail.com < WARNING
   rm WARNING
fi

