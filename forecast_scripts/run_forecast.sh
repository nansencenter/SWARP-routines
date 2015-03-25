# get restart
topaz_get_restart.sh # get restart file
# output year and date of forecast

#################################################################
# make infile
jday0=`date -d "today" '+%j'` # julian day of today (1=1st Jan => need to change)
jday_yest=`expr $jday0 - 1`   # julian day of TOPAZ (0=1st Jan)
jday_today=`printf '%3.3d' $jday_yest`
year_today=`date -u +%Y`
if [ ! $year_today -eq $ryear ]
then
	jday_new=`expr $jday_today + $jday + 1`
	jday_today=`printf '%3.3d' $jday_new`
fi
fc_days=5 # 5-day forecast check this works at/near year change
final_day=`expr $jday_today + $fc_days`
if [ $final_day -gt 364 ]
then
	fc_final_day=`expr $final_day - 365`
	fc_year=`expr $ryear + 1`
else
	fc_year=$year_today
fi
echo "Restart files of $jday-$ryear"
echo "Forecast final day $fc_final_day-$fc_year"


makeinfile4forecast.sh TP4 $ryear $jday $jday_today $final_day
#################################################################

#################################################################
# Launch job
xdir=$TP4_REALTIME/expt_01.1
cd $xdir
qsub pbsjob.sh
