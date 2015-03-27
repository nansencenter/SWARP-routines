# get restart

# paths that should be in bash_profile - does crontab not recognise these?
SWARP_ROUTINES=$HOME/GITHUB_REPOSITORIES/SWARP_routines
TP4_REALTIME=/work/timill/RealTime_Models/TP4a0.12

$SWARP_ROUTINES/forecast_scripts/topaz_get_restart.sh # get latest restart file
# output year and date of forecast

#textfile outputs:
out_restart=last_restart.txt

year_today=`date -u +%Y`
rname=`cat $out_restart`
rgen=${rname:0:3}   # eg TP4
ryear=${rname:10:4} # year of restart file
rday=${rname:15:3}  # julian day of restart file (1 Jan = 0)

#################################################################
# make infile
jday0=`date -d "today" '+%j'` # julian day of today (1=1st Jan => need to change)
jday_today0=`expr $jday0 - 1`   # julian day of TOPAZ (0=1st Jan)
jday_today=`printf '%3.3d' $jday_today0`

# if last restart was in different year to this year:
if [ ! $year_today -eq $ryear ]
then
	jday_today=`expr $jday_today + $rday + 1` # integer
        jday_today=`printf '%3.3d' $jday_today`   # 3 digits (compare to rday)
fi

fc_days=5 # 5-day forecast check this works at/near year change
final_day=`expr $jday_today + $fc_days`

# print to screen - work out if last day of forecast is in a different year to current year
ndays=`date --date="${year_today}-12-31" +%j` # days in current year
if [ $final_day -gt $((ndays-1)) ]
then
	fc_final_day=`expr $final_day - $ndays`
	fc_year=`expr $year_today + 1`
else
	fc_year=$year_today
fi
echo "Restart files of ${ryear}_$rday"
echo "Forecast final day ${fc_year}_$fc_final_day"

make_infile4forecast.sh $rgen $ryear $rday $jday_today $final_day
xdir=$TP4_REALTIME/expt_01.1
infile=$xdir/infile.in
if [ $rday == $jday_today ]
then
   # delete "today" line
   echo "restart day is from today"
   echo "- editing infile.in"
   sed '17d' $infile >> infile.in.replace
   mv infile.in.replace $infile
fi
#################################################################

#################################################################
# Launch job
cd $xdir

# want to save archive files (TP4archv*.[ab]) every 3 hours
cp $SWARP_ROUTINES/forecast_scripts/inputs/blkdat.input.archv_3h blkdat.input

# clean data directory before run
rm data/TP4DAILY*
rm data/TP4archv*

# launch job
qsub pbsjob.sh
#################################################################

#################################################################
# NB not in this script (~2h later)
# TODO save results (initial restart file,TP4DAILY*,TP4archv) on migrate
# TODO convert files to netcdf with hyc2proj
#################################################################
