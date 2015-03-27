#!/bin/bash
# This script will collect and archive the results of the local TP4 model

WDIR="/work/timill/RealTime_Models/TP4a0.12/expt_01.1/data"
SDIR="/migrate/timill/MODELS/TP4a0.12/expt_01.0/data"

cd $WDIR

eye=( TP4DAILY*.a )
for el in ${eye[@]}
do
   iy=${el:9:4}
   id=${el:14:3}
   fd=${el:23:3}
   fy=${el:18:4}
   adfil=${el}
   bdfil=TP4DAILY_${iy}_${id}_${fy}_${fd}.b
   pfil=TP4archv.${iy}_${fd}*
   tfil=TP4DAILYarchv.$iy_$id.tar.gz
   touch $tfil
   tar -zcvf $tfil $pfil $adfil $bdfil
   mkdir -p $SDIR/${yy}
   mv $tfil $SDIR/${yy}
   echo "Daily and archv files from $fd/$fy"
   echo "stored in $SDIR"
done
