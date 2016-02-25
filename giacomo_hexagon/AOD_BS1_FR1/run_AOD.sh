#!/bin/bash

startdate='20150308'
enddate='20150930'
sdj=$(date --date="$startdate" +%j)
edj=$(date --date="$enddate" +%j)
ndays=$(expr $edj - $sdj)

for n in $(seq 0 $ndays)
do
   hdate=$(date --date="$startdate + $n days" '+%Y%m%d')
   python AOD_BS1_FR1.py ${hdate} 'multi'
done
