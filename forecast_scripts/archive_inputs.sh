#!/bin/bash
# One script to execute both topaz_archive_restart.sh and wamnsea_update.sh
FCDIR=/home/nersc/timill/GITHUB-REPOSITORIES/SWARP-routines/forecast_scripts

$FCDIR/topaz_archive_restart.sh
$FCDIR/wamnsea_update.sh

