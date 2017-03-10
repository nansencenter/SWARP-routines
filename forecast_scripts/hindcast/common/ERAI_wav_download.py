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
   grid  = "0.5/0.5" # 0.5 degree grid
else:
   grid  = "1/1" # 1 degree grid

# dictionary of variable id's
varnos   = {'MWD':'230.140',
            'MWP':'232.140',
            'SWH':'229.140'}

varnames = ["MWD","MWP","SWH"] # variables to get

server   = ECMWFDataServer()

for nyear in range(2014,2015):
   # Set date and time 
   year  = str(nyear)
   if year == 2016:
      date  = str(year+"-01-01/to/"+year+"-03-31")
   else:
      date  = str(year+"-01-01/to/"+year+"-12-31")

   type="an"
   step="0"
   time="00:00:00/06:00:00/12:00:00/18:00:00"


   ## WAVE VARIABLES
   for varname in varnames:
      target   = varname+"_"+year+".nc"
      varno    = str(varnos[varname])

      print ""
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
