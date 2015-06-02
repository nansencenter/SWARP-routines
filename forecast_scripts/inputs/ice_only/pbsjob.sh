#! /bin/bash
#
#  Make sure I use the correct shell.
#
#PBS -S /bin/bash
#
# give the job a name
#PBS -N "TP4x012fc"
#
#  Specify the project the job belongs to
#
#PBS -A nn2993k
#
#  We want 24 hours on 51 cpu's:
#
#PBS -l walltime=01:10:00,mppwidth=133
#
#  The job needs 1 GB memory per cpu:
##PBS -l mppmem=1000mb
#
#  Send me an email on  a=abort, b=begin, e=end
#
#PBS -m a
#
#  Use this email address (check that it is correct):
#PBS -M timothy.williams@nersc.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o log/mpijob.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e log/mpijob.err
#

# ------------------- Fetch Environment ------------------------------
# -------- these are needed in preprocess scripts---------------------
NMPI=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
NOMP=`qstat -f $PBS_JOBID | awk '/mppdepth/ {print $3}'`
[ -z "$NOMP" ] && NOMP=0
export NOMP NMPI
echo "nmpi=$NMPI"
echo "nomp=$NOMP" # Not really used yet
echo "PBS_O_WORKDIR= $PBS_O_WORKDIR "
# -----------------End Fetch Environment -----------------------------

# Normal model executable:
EXEC=hycom
# Model executable if want to use PAT (performance analysis tool):
# EXEC=hycom+pat
# EXEC=hycom+apa

# SWARP post-processing option
# 0 - do nothing after postprocess.sh
# 1 - run ice-only PP:
#     /home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/process_FCresults.sh
# 2 - run waves-in-ice PP:
SWARP_PP=1

# Enter directory from where the job was submitted
cd $PBS_O_WORKDIR       ||  { echo "Could not go to dir $PBS_O_WORKDIR  "; exit 1; }

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }

# Transfer data files to scratch - must be in "expt_XXX" dir for this script
./preprocess.sh         ||  { echo "Preprocess had fatal errors "; exit 1; }

# Enter Scratch/run dir and Run model
cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
echo 'TW: launching job from:'
echo `pwd`
aprun -n $NMPI -m 1000M ./${EXEC}  # Run hycom

# Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
./postprocess.sh 

# extra processing for SWARP forecasts
if [ $SWARP_PP -eq 1 ]
then
   # SWARP post-processing - ice only version
   /home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/ice_only/process_FCresults.sh
elif [ $SWARP_PP -eq 2 ]
then
   # SWARP post-processing - waves version
   /home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/wavesice/process_FCresults_wav.sh
fi

datelist=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/datelist.txt
cday=$(cat $datelist | sed '1!d')
cyear=$(cat $datelist | sed '3!d')

# LAUNCH THE WAVE JOB ONLY IF WE HAVE ALREADY DOWNLOADED THE WAMNSEA PRODUCT
if [ -f /work/shared/nersc/msc/WAMNSEA/$cyear/forecasts/wam_nsea.fc.$cday.nc ]
then
   /home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts/wavesice/run_forecast_wav.sh
fi

exit $?
