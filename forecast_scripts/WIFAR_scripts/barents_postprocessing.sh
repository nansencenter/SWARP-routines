#!/bin/bash
# Called by barents_forecast.sh
# Two inputs : julday, year
# ======================================================================
# 1. Check if run is completed
# 2. Call barents_ice_forecast.sh
# 3. Archive results
# ======================================================================
MAINDIR=/home/nersc/bergh/Realtime
WORKDIR=/work/bergh/BS1a0.045/expt_01.4/data/
#CRONDIR=/home/nersc/dany/.cron/
ARCHDIR=/bcmhsm/nersc/bergh/WIFAR/BS1a0.045/Operational/
region=barents
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/acml/4.3.0/gfortran64/lib:/opt/acml/4.3.0/ifort64/lib
# BELOW USED BEFORE :
##export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/acml/5.1.0/gfortran64/lib:/opt/acml/5.1.0/ifort64/lib
# INSTAD SAME AS IN barents_forecast ??
export PATH=/opt/torque/2.5.11/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/lib/intel64/:/opt/acml/5.1.0/gfortran64/lib:/opt/acml/5.1.0/ifort64/lib
## LOCATION='bergh@hexagon.bccs.uib.no:/bcmhsm/bergh/WIFAR/BS1a0.045/Operational/'

# To have access to the function hyc2proj
#export PATH=/home/nersc/intissar/bin:$PATH
#. ~/.bash_profile
#. ~/.bashrc

if [ $# -ne 1 ] ; then
  echo "Usage: $0 julday"
  exit
fi

echo "-------------------------------------------------------------------------"
echo " barents_postprocessing.sh"
echo "-------------------------------------------------------------------------"
echo ""

julday=$1

day_start=`$MAINDIR/subprogs/julday2dayinyear_out $julday`
day_start=`echo 00$day_start | tail -4c`
year_start=`$MAINDIR/subprogs/julday2year_out $julday`
#thedate=`echo $julday | $MAINDIR/subprogs/julday2date_short`

# Today
julday_2=`expr $julday + 1`
day_2=`$MAINDIR/subprogs/julday2dayinyear_out $julday_2`
day_2=`echo 00$day_2 | tail -4c`
year_2=`$MAINDIR/subprogs/julday2year_out $julday_2`

# Tomorrow
julday_3=`expr $julday + 2`
day_3=`$MAINDIR/subprogs/julday2dayinyear_out $julday_3`
day_3=`echo 00$day_3 | tail -4c`
year_3=`$MAINDIR/subprogs/julday2year_out $julday_3`

# End at day T+2 12:00 UTC
julday_end=`expr $julday + 3`
day_end=`$MAINDIR/subprogs/julday2dayinyear_out $julday_end`
day_end=`echo 00$day_end | tail -4c`
year_end=`$MAINDIR/subprogs/julday2year_out $julday_end`


# Check if the run is completed
echo " * Check if run is completed"
echo ""

if [ -s ${WORKDIR}/BS1restart${year_end}_${day_end}_00.a ] ; then
  echo " * File ${WORKDIR}/BS1restart${year_end}_${day_end}_00.a is present. Run completed."
  echo ""
else
  echo " WARNING : ${WORKDIR}/BS1restart${year_end}_${day_end}_00.a has not been created."
  echo " * Check if daily averages have been created"
  echo ""

   if [ -s ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_3}_${day_3}.a ] ; then
     echo " * File ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_3}_${day_3}.a is present."
   else
     echo " WARNING : ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_3}_${day_3}.a has not been created."
   fi

   if [ -s ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_2}_${day_2}.a ] ; then
     echo " * File ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_2}_${day_2}.a is present."
   else
     echo " WARNING : ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_2}_${day_2}.a has not been created."
   fi

   if [ -s ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_start}_${day_start}.a ] ; then
     echo " * File ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_start}_${day_start}.a is present."
   else
     echo " ERROR   : ${WORKDIR}/BS1DAILY_${year_start}_${day_start}_${year_start}_${day_start}.a has not been created."
     echo " * Run not completed. Quit."
     echo ""
     exit
   fi
fi


# Update the webpage for the next forecast
cd $MAINDIR
echo " * Launching barents_ice_forecast.sh $region $julday $year_start"
echo ""
$MAINDIR/barents_ice_forecast.sh $region $julday $year_start

# ---------------------------------------------------------------------
#                             ARCHIVING
# ---------------------------------------------------------------------

echo " * Archiving activated for daily best guess and once per month for restart best guess."
echo ""

cd $WORKDIR
# ----------------------------   Daily  -------------------------------

nice tar zcvf  Daily_BS1_${year_start}_${day_start}.tar.gz BS1DAILY_${year_start}_${day_start}_${year_start}_${day_start}*
scp -i $HOME/.ssh/scriptkey  Daily_BS1_${year_start}_${day_start}.tar.gz ${ARCHDIR}/Daily/
rm Daily_BS1_${year_start}_${day_start}.tar.gz

# ----------------------------  Restart -------------------------------
# Save only if we are the first day of the month
# This gives year month day
#thedate=`date +"%Y %m %d"`
theday=`date +"%d"`
if [ ${theday} -eq 14 ]; then
   nice tar zcvf Restart_BS1_${year_start}_${day_start}.tar.gz BS1restart${year_start}_${day_start}*
   scp -i $HOME/.ssh/scriptkey  Restart_BS1_${year_start}_${day_start}.tar.gz ${ARCHDIR}/Restart/
   rm  Restart_BS1_${year_start}_${day_start}.tar.gz
fi

cd $MAINDIR 