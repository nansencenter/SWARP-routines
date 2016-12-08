if [ $# -ne 2 ]
then
   me=`basename $0`
   echo "Usage"
   echo "$me infile outfile"
   echo where infile and outfile are netcdf files 
   exit
fi
ncfil=$1
out=$2
vlist=model_depth,stereographic,icec,icetk,uice,vice,usurf,vsurf
ncks -d x,350,500 -d y,300,450 -v $vlist $ncfil.nc $out
