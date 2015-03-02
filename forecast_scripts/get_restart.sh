####!bash
# get latest restart
# draft!! to be tested
P=`pwd`

restart_dir="$P/test"      # directory with restarts
tp4dir="$P"                # location of TP4a0.12 directory (where forecast will be done)
xdir="$tp4dir/expt_01.0"   # location of expt directory

date_now=`date -ju "+%Y%m%d"`
year_now=${date_now:0:4}      # current year
month_now=${date_now:4:2}     # current month
day_now=${date_now:6:2}       # current day
year_prev=$((year_now-1))     # previous year

cd $restart_dir
echo "Looking for restarts in current year (${year_now}) or previous year..."
echo " "
if [ -f TP4restart${year_now}_???_*.tar.gz ]
then
   # loop over all files in current year
   # - last file is most recent date
   for f in TP4restart${year_now}_???_*.tar.gz
   do
      echo $f
   done
elif [ -f TP4restart${year_prev}_???_*.tar.gz ]
then
   # loop over all files in previous year
   # - last file is most recent date
   for f in TP4restart${year_prev}_???_*.tar.gz
   do
      echo $f
   done
else
   echo "No recent restarts in $restart_dir"
   echo "(none from $year_now or $year_prev)."
   echo " "
   exit
fi

echo "Most recent restart: $f"
echo "Unpacking and copying to $data_dir..."
echo " "
cd $P

# $f is now most recent restart
restart_year=${f:10:4}   # year of restart
jday=${f:15:3} # julian day of restart
hr=${f:19:2}   # hour of restart
restart_date=`jul2date $jday $restart_year 1 1`
restart_mon=${restart_date:4:2}
restart_day=${restart_date:6:2}

# names of restart files we will finally use
f0=TP4restart${restart_year}_${jday}_${hr}
afil=${f0}.a
bfil=${f0}.b
ufil=${f0}ICE.uf

# go to $data_dir
# if the most recent restart is NOT already there, then:
# - unpack restart there
# - rename files
data_dir=$xdir/data
echo $data_dir
echo $xdir

cd $data_dir
if [ ! -f $afil ]
then
   echo tar -xvf $restart_dir/$f
   tar -xvf $restart_dir/$f
   echo " "
   echo mv ${f0}_mem001.a       $afil
   echo mv ${f0}_mem001.b       $bfil
   echo mv ${f0}ICE_mem001.uf   $ufil
   echo " "
   #
   mv ${f0}_mem001.a       $afil
   mv ${f0}_mem001.b       $bfil
   mv ${f0}ICE_mem001.uf   $ufil
else
   echo "Restart files already present in $data_dir:"
   echo $afil
   echo $bfil
   echo $ufil
   echo " "
fi

# TODO calculate days between current date
# and latest restart - if too old (>2 weeks?) send warning
days_since_restart=`date2jul $date_now $restart_year $restart_mon $restart_day`
if [ $days_since_restart -ge 15 ]
then
   echo "Warning old restart"
   # TODO send email
fi

# now we need to edit:
# - infile.in
# - pbsjob.sh
# then do:
# - qsub pbsjob.sh
