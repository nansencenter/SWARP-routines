#!/bin/bash
# This script will collect and archive the results of the local TP4 model

echo "Collecting data produced in date `date +%d/%m/%Y`"

# defining all the dir that will be used
RTM=/work/timill/RealTime_Models
DFDIR=$RTM/TP4a0.12/expt_01.1/data
NCDIR=$RTM/post_processing/archv_netcdf
FODIR=$RTM/results/TP4a0.12/ice_only/final_output
WKDIR=$RTM/results/TP4a0.12/ice_only/work


