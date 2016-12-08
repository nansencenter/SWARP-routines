#!/usr/bin/env python
## This script downloads atmospheric forcing data
##
## Follow instruction here and register to ECMWF and get API key etc
##  https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets
## API key from https://api.ecmwf.int/v1/key/
## paste into: $HOME/.ecmwfapirc
##
## To get other fields, go to ECMWF web page: 
## http://apps.ecmwf.int/datasets/data/interim-full-daily
## select time period, fields and select "View the MARS request" to get python code 
## add 'format'    : "netcdf" to get NetCDF format, and
## rename "target": "CHANGEME",

from ecmwfapi import ECMWFDataServer
# server = ECMWFDataServer()

ATM   = 1

#set resolution
if 1:
   grid  = "0.5/0.5"
else:
   grid  = "1/1"

# Set date and time 
if 0:
   year  = "2016"
   date  = str(year+"-01-01/to/"+year+"-03-31")
else:
   year  = "2015"
   date  = str(year+"-01-01/to/"+year+"-12-31")

#area="Custom" "N90,W-180,S10,W180" # ????
# TODO THIS AREA LIMIT SEEMS NOT TO WORK
# area="Northern Hemisphere"
# select forecast or analised fields 
# forecast step = 0 and time 00:00/06:00 etc and type= "fc"
# analysed step = 3,6,9 or 12 ??? and type = "an"
type="an"
step="0"
time="00:00:00/06:00:00/12:00:00/18:00:00"
#type="fc"
#step="3/6/9/12"
#time="00:00:00/12:00:00"

## OPERATIONAL VARIABLES
if ATM==1:
   # atmospheric fields only
   varnos   = ['168.128','167.128','165.128','166.128','151.128','164.128']
   varnames = ["D2M"    ,"T2M"    ,"U10M"   ,"V10M"   ,"MSL"    ,"TCC"    ]
elif ATM==2:
   # atmospheric fields+sea ice fields
   varnos   = ['168.128','167.128','165.128','166.128','151.128','164.128','31.128','34.128','141.128']
   varnames = ["D2M"    ,"T2M"    ,"U10M"   ,"V10M"   ,"MSL"    ,"TCC"    ,"CICE"  ,"SST"   ,"HSNW"   ]
elif ATM==0:
   # sea ice fields
   varnos   = ['31.128','34.128','141.128']
   varnames = ["CICE","SST","HSNW"]

for i,varname in enumerate(varnames):
   target   = varname+"_"+year+".nc"
   varno    = str(varnos[i])

   server = ECMWFDataServer()

   print ""
   print i
   print varname
   print varno
   print target
   print ""


   server.retrieve({
      "class": "ei",
      "dataset": "interim",
      "date": date,
      "expver": "1",
      "grid": grid,
      "levtype": "sfc",
      "param": varno,
      "step": step,
      "stream": "oper",
      "time": time,
      "type": type,
      "format": "netcdf",
      "target": target,
   })
