# clean up files from panoply reading grb2 files

dirs[0]="eg_grib2"
dirs[1]="out"

ndirs=${#dirs[@]}

for n in `seq 0 $(($ndirs-1))`
do
   dd=${dirs[$n]}
   rm $dd/*.ncx
   rm $dd/*.gbx9
done
