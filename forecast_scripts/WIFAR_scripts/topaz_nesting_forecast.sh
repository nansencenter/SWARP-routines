#! /bin/bash
# This is the first script of the forecast script set for Fram Strait.
# It does not require any input and is called by cron.
# Cron settings can be changed by editing 'mycronfile' and updating it
# using 'crontab mycronfile'.
# ======================================================================
# 1. Download restart file from myocean.met.no  
# 2. Make sue the downloading worked
# 3. Generate infile.in
# 3. Submit the job and wait until it's finished
# 4. Launch topaz_forecast_postprocessing.sh
# ======================================================================
MAINDIR=/home/nersc/bergh/Realtime
WORKDIR=/work/bergh/TP4a0.12/expt_01.1/
EMAIL="jon.bergh@nersc.no"
export PATH=/opt/torque/torque-2.3.4-snap.200809221601/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/074/lib/intel64/:/opt/acml/4.3.0/gfortran64/lib:/opt/acml/4.3.0/ifort64/lib


if [ $# -ne 0 ]
then
  echo "Usage: $0 no arguments"
  exit
fi

echo "-------------------------------------------------------------------------"
echo " topaz_nesting_forecast.sh"
echo "-------------------------------------------------------------------------"

thedate=`date +"%Y %m %d"`

julday=`$MAINDIR/subprogs/date2julday $thedate`
#julday=`expr $julday`
julday=`expr $julday + 1`
#### force from defined julday relative to todays date ###
## julday=`expr $julday - 6`
##################################
todate=`$MAINDIR/subprogs/julday2date_human $julday`
echo "Run at julday = $julday"
echo "Run at date   = $todate"
echo ""

day=`$MAINDIR/subprogs/julday2dayinyear_out $julday`
day=`echo 00$day | tail -4c`
year=`$MAINDIR/subprogs/julday2year_out $julday`

################################################################
# Get restart files from myocean/met.no/ running at Vilje
################################################################
# HUSK! Remember to also set days in topaz_postprocessing.sh
# check for -days get into last year
echo "day: $day"
echo "year: $year"
echo "julday: $julday"

# Set start, end, and restart output dates
day1_julian=`expr $julday - 2`
#day2_julian=`expr $julday + 5`
#day3_julian=`expr $julday + 9`

day1_hycom=`$MAINDIR/subprogs/julday2dayinyear_out $day1_julian`
year=`$MAINDIR/subprogs/julday2year_out $day1_julian`
day2_hycom=`expr $day1_hycom + 7`
day3_hycom=`expr $day1_hycom + 11`
#
day1_hycom=`echo 00$day1_hycom | tail -4c`
day2_hycom=`echo 00$day2_hycom | tail -4c`
day3_hycom=`echo 00$day3_hycom | tail -4c`
#
echo " * Simulation will run from $year day $day1_hycom to day $day3_hycom"
echo ""

#if [  -s ${WORKDIR}/data/Metno/TP4restart${year}_${day1_hycom}_00.a ]
#then
#  echo " *Metno TP4restart${year}_${day1_hycom}_00.a exists"
#  echo ""
#  mv ${WORKDIR}/data/Metno/TP4restart${year}_${day1_hycom}_00* ${WORKDIR}/data/
#else
#  errormessage=" ERROR : topaz_nesting_forecast.sh : Metno TP4restart${year}_${day1_hycom}_00.a does not exist"
#  echo $errormessage
#  exit
#fi

##########################################################################
# Download restart file from myocean.met.no
#########################################################################

echo " * Download restart files from myocean.met.no if exists, else use old restart files"

# check if file exist in ./data, if so mv it to tmprestart (create tmprestart if needed)
if [  -s ${WORKDIR}/data/TP4restart${year}_${day1_hycom}_00.a ]
 then
  echo " Moving old restart files to ./data/tmprestart"
#    if ! [ -s ${WORKDIR}/data/tmprestart ]
#      echo " ./data/tmprestart doesn't exist, make dir tmprestart"
#      mkdir ${WORKDIR}/data/tmprestart
#    fi
    if [  -s ${WORKDIR}/data/tmprestart ]
     then
      echo " The directory - ${WORKDIR}/data/tmprestart - exist"
     else
      echo " The directory - ${WORKDIR}/data/tmprestart - DOES NOT EXIST"
      echo " Creating new directory: ${WORKDIR}/data/tmprestart"
      mkdir ${WORKDIR}/data/tmprestart
    fi
      mv ${WORKDIR}/data/TP4restart${year}_${day1_hycom}_00* ${WORKDIR}/data/tmprestart
fi

# go to data dir and download restart file from ftp server using ncftp
cd ${WORKDIR}/data
  cat > ncftp.in<<EOF
open myocean
cd ARC-MFC/ARC-METNO-ARC-TOPAZ4_2_PHYS-FOR/NERSC-MSC/
get TP4restart${year}_${day1_hycom}_00_mem001.a TP4restart${year}_${day1_hycom}_00_mem001.b TP4restart${year}_${day1_hycom}_00ICE.uf
set confirm-close no
bye
EOF
   ncftp < ncftp.in
   rm ncftp.in

mv TP4restart${year}_${day1_hycom}_00_mem001.a TP4restart${year}_${day1_hycom}_00.a
mv TP4restart${year}_${day1_hycom}_00_mem001.b TP4restart${year}_${day1_hycom}_00.b

cd ${WORKDIR}


###########################################################
# Download restart file direct from TOPAZ run on Ve #######
###########################################################

#VEDIR=/prod/forecast/sea/TOPAZ4
##VEDIR=/home/ntnu/fanf/ECMWFR
#echo " * Download restart files direct from TOPAZ on Ve"
#cd ${WORKDIR}/data
#scp -p -i ~/.ssh/nansenkey fanf@ve.hpc.ntnu.no:$VEDIR/TP4restart${year}_${day1_hycom}_00_mem001.a ${WORKDIR}/data 
#scp -p -i ~/.ssh/nansenkey fanf@ve.hpc.ntnu.no:$VEDIR/TP4restart${year}_${day1_hycom}_00_mem001.b ${WORKDIR}/data
#scp -p -i ~/.ssh/nansenkey fanf@ve.hpc.ntnu.no:$VEDIR/TP4restart${year}_${day1_hycom}_00ICE.uf ${WORKDIR}/data
#
#mv TP4restart${year}_${day1_hycom}_00_mem001.a TP4restart${year}_${day1_hycom}_00.a
#mv TP4restart${year}_${day1_hycom}_00_mem001.b TP4restart${year}_${day1_hycom}_00.b
#
#cd ${WORKDIR}


# check if new restart file exist in ./data, else cp old file from tmprestart and proceed
if [  -s ${WORKDIR}/data/TP4restart${year}_${day1_hycom}_00.a ]
 then
   echo " * Restart file TP4restart${year}_${day1_hycom}_00.a copied to ./data, proceed forecast ..."
 else
   echo " ***** Restart file TP4restart${year}_${day1_hycom}_00.a doesn't exist ***** "
   echo " ******************************************************************** "
#    if [  -s ${WORKDIR}/data/tmprestart/TP4restart${year}_${day1_hycom}_00.a ]
#     then
#       echo " * Use old restart files from end of last forecast, proceed forecast ..."
#       mv ${WORKDIR}/data/tmprestart/TP4restart${year}_${day1_hycom}_00* ${WORKDIR}/data
#    fi
   exit
fi

##############################################################################
#But so far we run the Intissar's Free run version of TOPAZ4
##############################################################################


#day1_hycom=`expr $day + 3`  # HYCOM days > 365 in case of new year shift
#day2_hycom=`expr $day + 7`
#day3_hycom=`expr $day + 8`

#day1_hycom=`echo 00$day1_hycom | tail -4c`
#day2_hycom=`echo 00$day2_hycom | tail -4c`
#day3_hycom=`echo 00$day3_hycom | tail -4c`

#echo " * Simulation will run from $year day $day1_hycom to day $day3_hycom"
#echo ""

#if [  -s ${WORKDIR}/data/TP4restart${year}_${day}_00.a ]
#then
#  echo " *Free run TP4restart${year}_${day}_00.a exists"
#  echo ""
#else
#  errormessage=" ERROR : topaz_nesting_forecast.sh : Free run TP4restart${year}_${day}_00.a does not exist"
#  echo $errormessage
#  exit
#fi

echo " * Generating infile.in"
echo ""
#${MAINDIR}/makeinfile4forecast.sh TP4 $year $day $day1_hycom $day2_hycom $day3_hycom
${MAINDIR}/makeinfile4topazforecast.sh TP4 $year $day1_hycom $day2_hycom $day3_hycom

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

cd ${WORKDIR}


echo " * Integrating forecast"
echo ""
##/opt/torque/torque-2.3.4-snap.200809221601/bin/qsub pbsjob_TP4_RT.sh
/opt/torque/2.5.11/bin/qsub pbsjob_TP4_RT.sh

while  ! [ -f ${WORKDIR}/data/TP4restart${year}_${day3_hycom}_00.a ]
do
  echo " * Waiting 20 minutes ..."
  echo ""
  sleep 1200
done

# Proceed with the postprocessing
echo " * Save Topaz restart files from $hycom_day1 - 7"
echo ""
${MAINDIR}/topaz_postprocessing.sh $julday
