#!/bin/bash
# Called by framstrait_forecast.sh
# Two inputs : julday, year
# ======================================================================
# 1. Check if run is completed
# 2. Call framstrait_ice_forecast.sh
# 3. Archive results
# ======================================================================
rungen=TP4
MAINDIR=/home/nersc/bergh/Realtime
WORKDIR=/work/bergh/${rungen}a0.12/expt_01.1/data
#CRONDIR=/home/nersc/dany/.cron/
ARCHDIR=/migrate/bergh/TOPAZ/${rungen}a0.12/For_Operational_Runs/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/acml/4.3.0/gfortran64/lib:/opt/acml/4.3.0/ifort64/lib

# To have access to the function hyc2proj
#export PATH=/home/nersc/intissar/bin:$PATH
#. ~/.bash_profile
#. ~/.bashrc

if [ $# -ne 1 ] ; then
  echo "Usage: $0 julday"
  exit
fi

echo "-------------------------------------------------------------------------"
echo " topaz_postprocessing.sh"
echo "-------------------------------------------------------------------------"
echo ""

julday=$1

day_start=`$MAINDIR/subprogs/julday2dayinyear_out $julday`
day_start=`echo 00$day_start | tail -4c`
year_start=`$MAINDIR/subprogs/julday2year_out $julday`
#thedate=`echo $julday | $MAINDIR/subprogs/julday2date_short`

# Day not used 
julday_2=`expr $julday + 3`
day_2=`$julday_2 | $MAINDIR/subprogs/julday2dayinyear_out $julday_2`
day_2=`echo 00$day_2 | tail -4c`
year_2=`$MAINDIR/subprogs/julday2year_out $julday_2`



# End at day T+9 
julday_end=`expr $julday + 9`
day_end=`$MAINDIR/subprogs/julday2dayinyear_out $julday_end`
day_end=`echo 00$day_end | tail -4c`
year_end=`$MAINDIR/subprogs/julday2year_out $julday_end`

# Check if the run is completed
echo " * Check if run is completed"
echo ""

if [ -s ${WORKDIR}/${rungen}restart${year_end}_${day_end}_00.a ] ; then
   echo " * File ${WORKDIR}/${rungen}restart${year_end}_${day_end}_00.a is present. Run completed."
   echo ""
else
   echo " * File ${WORKDIR}/${rungen}restart${year_end}_${day_end}_00.a is NOTpresent."
   echo " * Run not completed. Quit."
   echo ""
   exit
fi



# ---------------------------------------------------------------------
#                             ARCHIVING
# ---------------------------------------------------------------------
####### Need to save nesting files +restart occasionnaly!!!!
echo " * Archiving activated for daily best guess and once per month for restart best guess."
echo ""

cd $WORKDIR
# ----------------------------   Daily  -------------------------------

#nice tar zcvf  Daily_${rungen}_${year_start}_${day_start}.tar.gz ${rungen}DAILY_${year_start}_${day_start}_${year_start}_${day_start}*
#scp -i $HOME/.ssh/scriptkey  Daily_${rungen}_${year_start}_${day_start}.tar.gz ${ARCHDIR}/Daily/
#rm Daily_${rungen}_${year_start}_${day_start}.tar.gz
#rm ${rungen}DAILY_${year_start}_*

# ----------------------------  Restart -------------------------------
nice tar zcvf Restart_${rungen}_${year_start}_${day_start}.tar.gz ${rungen}restart${year_start}_${day_start}*
scp -i $HOME/.ssh/scriptkey  Restart_${rungen}_${year_start}_${day_start}.tar.gz ${ARCHDIR}/Restart/
rm  Restart_${rungen}_${year_start}_${day_start}.tar.gz
#rm  ${rungen}restart${year_2}_${day_2}*

cd $MAINDIR 
