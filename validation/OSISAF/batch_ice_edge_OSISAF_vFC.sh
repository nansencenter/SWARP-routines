# cdate=20160109 #start day
# fdate=20170306 #final day

cdate=20160801 #start day
fdate=20160831  #final day

while [ $cdate -le $fdate ]
do
   $SWARP_ROUTINES/validation/OSISAF/ice_edge_OSISAF_1obs.sh $cdate 1 0
   cdate=`date --date="$cdate +1days" +%Y%m%d`
done
