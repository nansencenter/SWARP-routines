#! /bin/bash
# This is the first script of the forecast script set for Fram Strait.
# It does not require any input and is called by cron.
# Cron settings can be changed by editing 'mycronfile' and updating it
# using 'crontab mycronfile'.
# ======================================================================
# 1. Update ECMWFR atmospheric forecast files from met.no at Ve 
# ======================================================================
MAINDIR=/home/nersc/bergh/Realtime
SHOST=fanf@ve.hpc.ntnu.no
SDIR=/home/ntnu/fanf/ECMWFR/
LDIR=/work/shared/nersc/ECMWFR_T799
LKEY=/home/nersc/bergh/.ssh/nansenkey

echo 'Update ECMWFR atmospheric forecast from met.no via ve.hpc.ntnu.no'

time rsync -vzrh -e "ssh -i $LKEY" $SHOST:$SDIR $LDIR


