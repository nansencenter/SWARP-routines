#! /bin/bash
#
#  Make sure I use the correct shell.
#
#PBS -S /bin/bash
#
# give the job a name
#PBS -N "MIZ_BS1_FR1"
#
#  Specify the project the job belongs to
#
#PBS -A nn2993k
#
#  We want .5 hours on 1 cpu:
#
#PBS -l walltime=00:30:00,mppwidth=1
#
#  The job needs 1 GB memory per cpu:
##PBS -l mppmem=1000mb
#
#  Send me an email on  a=abort, b=begin, e=end
#
#PBS -m a
#
#  Use this email address (check that it is correct):
#PBS -M gcmdnt90@gmail.com
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

# Enter directory from where the job was submitted
cd $PBS_O_WORKDIR       ||  { echo "Could not go to dir $PBS_O_WORKDIR  "; exit 1; }

mkdir -p log
EXEC=run_MIZ_BS1_FR1.sh
echo 'Launching job from:'
echo `pwd`
aprun -n $NMPI -m 1000M ./${EXEC}  # Run hycom

exit $?
