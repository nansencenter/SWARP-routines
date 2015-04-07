#! /bin/bash
# This is the first script of the forecast script set for Barents Sea.
# It does not require any input and is called by cron.
# Cron settings can be changed by editing 'mycronfile' and updating it
# using 'crontab mycronfile'.
# ======================================================================
# 1. Check if restart file present
# 2. Generate infile.in
# 3. Submit the job and wait until it's finished
# 4. Launch barents_forecast_postprocessing.sh
# ======================================================================
MAINDIR=/home/nersc/bergh/Realtime
WORKDIR=/work/bergh/BS1a0.045/expt_01.4/
#EMAIL="jon.bergh@nersc.no"
#export PATH=/opt/torque/torque-2.3.4-snap.200809221601/bin:$PATH
export PATH=/opt/torque/2.5.11/bin:$PATH
#export LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/074/lib/intel64/:/opt/acml/4.3.0/gfortran64/lib:/opt/acml/4.3.0/ifort64/lib
export LD_LIBRARY_PATH=/opt/intel/lib/intel64/:/opt/acml/5.1.0/gfortran64/lib:/opt/acml/5.1.0/ifort64/lib

if [ $# -ne 0 ]
then
  echo "Usage: $0 no arguments"
  exit
fi

echo "-------------------------------------------------------------------------"
echo " barents_forecast.sh"
echo "-------------------------------------------------------------------------"

thedate=`date +"%Y %m %d"`
julday=`$MAINDIR/subprogs/date2julday $thedate`
#julday=`expr $julday`
julday=`expr $julday + 1`
###### RESTART CHANGE TEMPORARELY ###
## julday=`expr $julday - 1`
#################################
julday_end=`expr $julday + 3`
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

day_end=`$MAINDIR/subprogs/julday2dayinyear_out $julday_end`
year_end=`$MAINDIR/subprogs/julday2year_out $julday_end`
day_end=`echo 00$day_end | tail -4c`
#year_end=`echo 00$year_end | tail -4c`


echo " * Simulation will run from $year day $day to day $day3_hycom"
echo ""
echo " * Last restart file will be ${WORKDIR}/data/BS1restart${year_end}_${day_end}_00.a"
echo ""

if [  -s ${WORKDIR}/data/BS1restart${year}_${day}_00.a ]
then
  echo " * BS1restart${year}_${day}_00.a exists"
  echo ""
else
  errormessage=" ERROR : barents_forecast.sh : BS1restart${year}_${day}_00.a does not exist"
  echo $errormessage
  exit
fi

echo " * Generating infile.in"
echo ""
${MAINDIR}/makeinfile4forecast.sh BS1 $year $day $day1_hycom $day2_hycom $day3_hycom

if [ -s ${MAINDIR}/infiles/infile.in ]
then
  mv ${MAINDIR}/infiles/infile.in $WORKDIR/
else
  errormessage=" ERROR : barents_forecast.sh : did not generate infile.in"
  echo $errormessage
  exit
fi

echo " * Submitting job"
echo ""


## echo " * Preparing forcing files"
## echo ""
## echo " * Update ECMWFR files in work/shared/nersc direct from TOPAZ on Ve"
## rsync -vzr -e "ssh -i /home/nersc/bergh/.ssh/nansenkey" fanf@ve.hpc.ntnu.no:/home/ntnu/fanf/ECMWFR/ /work/shared/nersc/ECMWFR_T799

#$HOME/bin/forfun_nersc ecnc old

echo " * Integrating 72 hours forecast. Will process results when run completed"
echo ""
cd ${WORKDIR}
#/opt/torque/torque-2.3.4-snap.200809221601/bin/qsub pbsjob_BS1_RT.sh
/opt/torque/2.5.11/bin/qsub pbsjob_BS1_RT.sh

count=0
while  ! [ -f ./data/BS1restart${year_end}_${day_end}_00.a ]
do
  if [  $count -ge 36 ]
  then
    echo " ***********************************************************************************"
    echo " ******* PBS JOB NOT COMPLETED AFTER 6 HOURS - EXIT barents_forecast.sh ************"
    echo " ***********************************************************************************"
    exit 
  fi
  count=`expr $count + 1`
  echo " * Count = $count PBS JOB running - Waiting 10 minutes ..."
  echo ""
  sleep 600
done

# Proceed with the postprocessing
${MAINDIR}/barents_postprocessing.sh $julday
echo " * Launching barents_forecast_postprocessing.sh $julday in jobname: RT-BS1 "
echo ""

#cd ${MAINDIR}
#cat pbs_postprocess.sh | sed \
#   -e "s/region/barents/g" \
#   -e "s/julday/$julday/g" \
#   > pbs_post_process_barents.sh
#
#/opt/torque/torque-2.3.4-snap.200809221601/bin/qsub pbs_postprocess_barents.sh 


