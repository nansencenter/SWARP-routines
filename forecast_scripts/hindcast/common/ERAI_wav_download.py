#!/usr/bin/env python
## This script download atmosphere, sea ice, ocean and wave data
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


#set resolution
if 0:
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

type="an"
step="0"
time="00:00:00/06:00:00/12:00:00/18:00:00"


## WAVE VARIABLES
varnos   = ['230.140','232.140','229.140']
varnames = ["MWD","MWP","SWH"]

for i,varname in enumerate(varnames):
   target   = varname+"_"+year+".nc"
   varno    = str(varnos[i])
   
   server   = ECMWFDataServer()
   
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
      "stream": "wave",
      "time": "00:00:00/06:00:00/12:00:00/18:00:00",
      "type": type,
      "format": "netcdf",
      "target": target,
   })
