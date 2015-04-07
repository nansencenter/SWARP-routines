## Run loop over julian days to update forecast when it has been down for several days


MAINDIR=$HOME/Realtime

# Start time year, month, day in month
Y=2013
M=06
D=01

NRD=2
#Y2=2013
#M2=06
#D2=2

# julian start day
julstart=`$MAINDIR/subprogs/date2julday $Y $M $D`
# julian end date (last day to run forecast)
#julend=`$MAINDIR/subprogs/date2julday $Y2 $M2 $D2`

echo "Update forecast from"
echo "julstart: $julstart"
#echo "julend  : $julend"

for i in {1..2}
do
julday=`expr $julstart + $i - 1`

   echo "***************************************"
   echo "*** UPDATE FORECAST AND WEBPAGE FOR ***"
   echo "*** DATE    : $Y $M $D "
   echo "*** JULDAY  : $julday "
   echo ""

#$MAINDIR/barents_ice_forecast.sh barents $date $Y	
#$MAINDIR/barents_ice_forecast_WIM.sh barents $date $Y

   echo "*** UPDATE FINISHED "
   echo ""

done


