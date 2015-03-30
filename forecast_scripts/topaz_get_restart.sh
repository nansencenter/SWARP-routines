#!/bin/bash
#Get latest restart from the internal repo to the working dir

# textfile output:
out_restart=last_restart.txt

rdir="/migrate/timill/restarts/TP4a0.12/SWARP_forecasts"      # directory with restarts

# no variable -> uses tp4dir ; #1 manually insert the dir ; >#2 error
if [ $# -eq 1 ]
then
   v1=$1
elif [ $# -eq 0 ]
then
   echo "Using default - proper location"
   v1=0
else
   echo "Wrong number of inputs to get_restart.sh"
   exit
fi

if [ $v1 -eq 0 ]
then
   #proper location
   tp4dir="/work/timill/RealTime_Models/TP4a0.12/" # location of TP4a0.12 directory (where forecast will be done)
   xdir="$tp4dir/expt_01.1"                     # location of expt directory
   ddir="$xdir/data"
else
   #test location
   xdir="$HOME/giacomo"
   mkdir -p $xdir/data
   ddir="$xdir/data"
fi

cyear=`date -u +%Y`			# current year
cmon=`date -u +%m`			# current month
cday=`date -d "yesterday" '+%j'`	# current day
pyear=`expr $cyear - 1`			# previous year
ceye=`find /${rdir}/${cyear} -name TP4restart${cyear}*`
ceye=( $ceye )
peye=`find /${rdir}/${pyear} -name TP4restart${pyear}*`
peye=( $peye )
cyfil=${rdir}/${cyear}/TP4restart${cyear}*
pyfil=${rdir}/${pyear}/TP4restart${pyear}*
shopt -s nullglob

#Looking for restarts in current year (${cyear}) or previous year..."
if [ ${#cyfil[@]} -gt 0 ]
then
   # loop over all files in current year
   # - last file is most recent date
   for f in ${ceye}
   do
	f=$(basename $f)
   done
elif [ ${#pyfil[@]} -gt 0 ]
then
   cd ${rdir}/${peye}
   # loop over all files in previous year
   # - last file is most recent date
   for f in ${pyfil}
   do
	f=$(basename $f)
   done
else
   echo "No recent restarts in ${rdir}"
   echo "(none from ${cyear} or ${pyear})."
   echo " "
   exit
fi

echo "Most recent restart: $f"
echo "Unpacking..."
echo " "
echo $f > $out_restart

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
   tar -zxvf $rdir/${ryear}/$f0.tar.gz
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

# calculate days between current date
# and latest restart - if too old (>2 weeks?) send warning
dsr=`expr $cday - $jday` #Days Since Restart

if [ $dsr -gt 13 ]
then
   echo "Restarts too old"
   exit
fi

# now we need to edit:
# - infile.in
# - pbsjob.sh
# then do:
# - qsub pbsjob.sh
