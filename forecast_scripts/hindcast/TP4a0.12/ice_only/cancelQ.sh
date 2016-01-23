lst1=(`qstat | grep timill |grep TP4 |grep HC |grep " Q "`) # jobs in queue
N1=${#lst1[@]}

for n in `seq 0 6  $((N1-1))`
do
   jno=${lst1[$n]}
   echo "qdel $jno"
   qdel $jno
done
