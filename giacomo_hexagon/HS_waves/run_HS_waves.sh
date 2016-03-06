#!/bin/bash

startdate='20150302'
enddate='20150930'
sdj=$(date --date="$startdate" +%j)
edj=$(date --date="$enddate" +%j)
ndays=$(expr $edj - $sdj)

for n in $(seq 0 $ndays)
do
   hdate=$(date --date="$startdate + $n days" '+%Y%m%d')
   python HS_waves.py ${hdate}
done
