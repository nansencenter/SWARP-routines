#! /bin/bash
# This is the first script of the forecast script set for Fram Strait.
# It does not require any input and is called by cron.
# Cron settings can be changed by editing 'mycronfile' and updating it
# using 'crontab mycronfile'.
# ======================================================================
# 1. Check if restart file present
# 2. Generate infile.in
# 3. Submit the job and wait until it's finished
# 4. Launch framstrait_forecast_postprocessing.sh
# ======================================================================
MAINDIR=/home/nersc/bergh/Realtime
WORKDIR=/work/bergh/FR1a0.03/expt_01.1/
EMAIL="jon.bergh@nersc.no"
export PATH=/opt/torque/torque-2.3.4-snap.200809221601/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/074/lib/intel64/:/opt/acml/4.3.0/gfortran64/lib:/opt/acml/4.3.0/ifort64/lib


if [ $# -ne 0 ]
then
  echo "Usage: $0 no arguments"
  exit
fi

echo "-------------------------------------------------------------------------"
echo " framstrait_forecast.sh"
echo "-------------------------------------------------------------------------"

thedate=`date +"%Y %m %d"`
julday=`$MAINDIR/subprogs/date2julday $thedate`
#julday=`expr $julday`
julday=`expr $julday + 1`

# Force start date
#julday=22686

echo 'Runnning framstrait_forecast.sh for julday+1='

echo "julday = $julday"
echo ""
day=`$MAINDIR/subprogs/julday2dayinyear_out $julday`
day=`echo 00$day | tail -4c`
year=`$MAINDIR/subprogs/julday2year_out $julday`

day1_hycom=`expr $day + 1`  # HYCOM days > 365 in case of new year shift
day2_hycom=`expr $day + 2`
day3_hycom=`expr $day + 3`

day1_hycom=`echo 00$day1_hycom | tail -4c`
day2_hycom=`echo 00$day2_hycom | tail -4c`
day3_hycom=`echo 00$day3_hycom | tail -4c`

echo " * Simulation will run from $year day $day to day $day3_hycom"
echo ""

if [  -s ${WORKDIR}/data/FR1restart${year}_${day}_00.a ]
then
  echo " * FR1restart${year}_${day}_00.a exists"
  echo ""
else
  errormessage=" ERROR : framstrait_forecast.sh : FR1restart${year}_${day}_00.a does not exist"
  echo $errormessage
  exit
fi

echo " * Generating infile.in"
echo ""
${MAINDIR}/makeinfile4forecast.sh FR1 $year $day $day1_hycom $day2_hycom $day3_hycom

if [ -s ${MAINDIR}/infiles/infile.in ]
then
  mv ${MAINDIR}/infiles/infile.in $WORKDIR/
else
  errormessage=" ERROR : framstrait_forecast.sh : did not generate infile.in"
  echo $errormessage
  exit
fi

echo " * Submitting job"
echo ""


#echo " * Preparing forcing files"
#echo ""
#$HOME/bin/forfun_nersc ecnc old

echo " * Integrating 72 hours forecast. Will process results when run completed"
echo ""
cd ${WORKDIR}
/opt/torque/torque-2.3.4-snap.200809221601/bin/qsub pbsjob_FR1_RT.sh


while  ! [ -f ./data/FR1restart${year}_${day3_hycom}_00.a ]
do
  echo " * Waiting 20 minutes ..."
  echo ""
  sleep 1200
done

# Proceed with the postprocessing
echo " * Launching framstrait_forecast_postprocessing.sh $julday"
echo ""
#cd ${MAINDIR}
#cat pbsjob_template_postprocess.sh | sed -e "s/julday/${julday}" >pbs_a.sh
#cat pbs_a.sh | sed -e "s/region/framstrait">pbsjob_FR_postprocess.sh
#rm pbs_a.sh
#/opt/torque/torque-2.3.4-snap.200809221601/bin/qsub pbsjob_FR_postprocess.sh
#Call the following routine
${MAINDIR}/framstrait_postprocessing.sh $julday
