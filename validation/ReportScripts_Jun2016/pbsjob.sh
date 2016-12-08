#! /bin/bash
#
#  Make sure I use the correct shell.
#
#PBS -S /bin/bash
#
# give the job a name
###PBS -N "AODs_all"
#PBS -N "MIZmap_all"
#
#  Specify the project the job belongs to
#
#PBS -A nn2993k
#
#  We want .5 hours on 1 cpu:
#
#PBS -l walltime=00:20:00,mppwidth=1
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

module load python/2.7.9-dso
module load spatialindex/1.8.5
# python="/work/apps/python/2.7.9-dso-gnu/bin/python"

OPT=2
echo $PYTHONPATH

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

echo 'Launching job from:'
echo `pwd`

if [ $OPT -eq 0 ]
then
   pyscript="test_import.py"
elif [ $OPT -eq 1 ]
then
   pyscript="AODs_all.py"
elif [ $OPT -eq 2 ]
then
   pyscript="MIZmap_all.py"
fi

echo "Running $pyscript"
# aprun -n $NMPI -m 1000M $python $pyscript
aprun -B python $pyscript

echo "Exit qsub"
exit $?
