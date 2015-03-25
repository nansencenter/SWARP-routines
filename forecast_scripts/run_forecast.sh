# get restart
get_restart.sh # get restart file
# TODO output year and date of forecast

#################################################################
# make infile
jday=`date -d "today" '+%j'` # julian day of today (1=1st Jan => need to change)
jday_today=`expr($jday-1)`   # julian day of TOPAZ (0=1st Jan)
year_today= #TODO
if [ ! $year_today -eq $ryear ]
then
   ndays= # TODO (no of days in $ryear)
   $jday_today=$((jday_today+ndays))
fi

fc_days=5 # 5-day forecast TODO check this works at/near year change
final_day=$((jday_today+fc_days))

makeinfile4forecast.sh TP4 $ryear $rday $jday_today $final_day
#################################################################

#################################################################
# Launch job
xdir=$TP4_REALTIME/expt_01.1
cd $xdir
qsub pbsjob.sh
