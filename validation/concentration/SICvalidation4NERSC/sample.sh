#!/bin/bash

fimexbin="/usr/bin"
Rbin="/usr/bin"
IOpath="."
WORKpath="."
SRCpath="."

day1="20150413"
date1="2015-04-13"
date2="2015-04-19"


# Ice chart results:

ic_file=${WORKpath}/"iceChartResults.nc"
ic_clock="14:00:00"
icDefFile=${IOpath}/iceChart.ncml

test -e $ic_file $ && /bin/rm $ic_file

varList_OSISAF="time xc yc lon lat crs ice_concentration"
varString_OSISAF=""
for varName in ${varList_OSISAF}
do
  varString_OSISAF=${varString_OSISAF}" --extract.selectVariables="${varName}
done

ncString_OSISAF=" --input.type=ncml --output.type=nc4"
fimexArgs="--input.file="${icDefFile}${ncString_OSISAF}${varString_OSISAF}" --output.file="${ic_file}
${fimexbin}/fimex ${fimexArgs} --extract.reduceTime.start="${date1}\ ${ic_clock}" --extract.reduceTime.end="${date2}\ ${ic_clock}"

echo "...got ice chart results"


# TOPAZ results:

tp_file=${WORKpath}/"TOPAZresults.nc"
tp_clock="00:00:00"

test -e $tp_file $ && /bin/rm $tp_file

varList_TOPAZ="time x y longitude latitude polar_stereographic fice"
varString_TOPAZ=""
for varName in ${varList_TOPAZ}
do
  varString_TOPAZ=${varString_TOPAZ}" --extract.selectVariables="${varName}
done

tp_href="http://thredds.met.no/thredds/dodsC/topaz/dataset-topaz4-arc-myoceanv2-"${day1}
ncString_TOPAZ=" --input.type=netcdf --output.type=netcdf"
fimexArgs="--input.file="${tp_href}${ncString_TOPAZ}${varString_TOPAZ}" --output.file="${tp_file}
${fimexbin}/fimex ${fimexArgs} --extract.reduceTime.start="${date1}\ ${tp_clock}" --extract.reduceTime.end="${date2}\ ${tp_clock}"

echo "...got TOPAZ results"


# Perform validation:

bdate=${date1}
export IOpath SRCpath WORKpath bdate
${Rbin}/R CMD BATCH ${SRCpath}/SIvalidate.R

echo "...done"


exit
