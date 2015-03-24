#!/bin/bash
#Get latest restart from the internal repo to the working dir

rdir="/migrate/timill/restarts/TP4a0.12/SWARP_forecasts"      # directory with restarts
echo "Type 0 for launch, 1 for test"
read vl
if [ $vl -eq 0 ]
then
   #proper location
   tp4dir="/work/timill/Model_Setups/TP4a0.12/" # location of TP4a0.12 directory (where forecast will be done)
   xdir="$tp4dir/expt_01.0"                     # location of expt directory
else
   #test location
   xdir="$HOME/giacomo"
   mkdir -p $xdir/data
fi

cyear=`date -u +%Y`			# current year
cmon=`date -u +%m`			# current month
cday=`date -d "yesterday" '+%j'`	# current day
pyear=`expr $cyear - 1`			# previous year

cd $rdir/${cyear}
#Looking for restarts in current year (${cyear}) or previous year..."
if [ -f TP4restart${cyear}*.tar.gz ]
then
   # loop over all files in current year
   # - last file is most recent date
   for f in TP4restart${cyear}_*.tar.gz
   do
      echo $f
   done
elif [ -f ${rdir}/${pyear}/TP4restart${pyear}_*.tar.gz ]
then
   cd ${rdir}/${pyear}
   # loop over all files in previous year
   # - last file is most recent date
   for f in TP4restart${pyear}_*.tar.gz
   do
      echo $f
   done
else
   echo "No recent restarts in ${rdir}"
   echo "(none from ${cyear} or ${pyear})."
   echo " "
   exit
fi

echo "Most recent restart: $f"
echo "Unpacking and copying to $xdir..."
echo " "

# $f is now most recent restart
ryea=${f:10:4}	# year of restart
jday=${f:15:3}	# julian day of restart
hr=${f:19:2}	# hour of restart

# names of restart files we will finally use
f0=TP4restart${ryea}_${jday}_${hr}
afil=${f0}_mem001.a
bfil=${f0}_mem001.b
ufil=${f0}ICE.uf


# go to the data in $xdir
# if the most recent restart is NOT already there, then:
# - unpack restart there
# - rename files

ddir="$xdir/data"

if [ ! -f $afil ]
then
   echo tar -xvf $rdir/$f
   tar -xvf $rdir/$f
   echo " "
   echo mv $afil $ddir
   echo mv $bfil $ddir
   echo mv $ufil $ddir
   echo " "
else
   echo "Restart files already present in $data_dir:"
   echo $afil
   echo $bfil
   echo $ufil
   echo " "
fi

# calculate days between current date
# and latest restart - if too old (>2 weeks?) send warning
dsr=`expr $cday - $jday` #Days Since Restart
if [ $dsr -ge 14 ]
then
   echo "Restarts too old"
   exit
fi

# now we need to edit:
# - infile.in
# - pbsjob.sh
# then do:
# - qsub pbsjob.sh
