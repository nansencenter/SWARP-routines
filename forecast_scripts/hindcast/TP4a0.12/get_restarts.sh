# script to copy all restart files to $TP4_REALTIME/expt_01.4/data
# - other expt's will link to them
source $SWARP_ROUTINES/source_files
cd $TP4_REALTIME/expt_01.4/data

md=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/2015/
list=(`ls $md/*gz`)
list2=(`ls $md/../2016/*gz`)
Nlist=${#list[@]}
inext=Nlist
for loop_i in `seq 1 ${#list2[@]}`
do
   list[$inext]=${list2[$((loop_i-1))]}
   inext=$((inext+1))
   Nlist=$((Nlist+1))
done

for rno in `seq 1 $Nlist`
do
   tfil0=${list[$((rno-1))]}
   tfil=`basename $tfil0`
   year=${tfil:10:4}
   day=${tfil:15:3}
   echo $tfil $year $day
   cp $tfil0 .
   tar -xf $tfil
   rm $tfil
done

for afil in *_mem001.a
do
   bfil=${afil%.a}.b
   base=${afil%_mem001.a}
   mv $afil $base.a
   mv $bfil $base.b
done
