#!/bin/bash
# Check if backup folders exist before moving 
#  SWARP forecast from /scratch to /project/ area
# INPUT
# 1) HYCOMreg, tex TP4a0.12, ti get correct path and FCtype
# 2) forecast version, same as folder name, tex "ice_only"

proj=NS2993K
HYCOMreg=$1
FCtype=$2
backuptype=$3
path1=/projects/$proj/SWARP_FC/$backuptype
path2=/projects/$proj/SWARP_FC/$backuptype/$HYCOMreg
path3=/projects/$proj/SWARP_FC/$backuptype/$HYCOMreg/$FCtype
[ ! -d $path1 ] && mkdir -p $path1 
[ ! -d $path2 ] && mkdir -p $path2
[ ! -d $path3 ] && mkdir -p $path3
