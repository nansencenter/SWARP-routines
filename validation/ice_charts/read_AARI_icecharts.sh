indir="/work/shared/nersc/msc/AARI_icecharts"


for reg in Bar Gre
do
   for year in 2015 2016
   do
      for crit in FA_only FB_only
      do
         if [ $crit=="FA_only" ]
         then
            outdir=$indir/$reg/MIZ_FA
            mkdir -p $outdir
         else
            outdir=$indir/$reg/MIZ_FB
            mkdir -p $outdir
         fi
         python read_icechart.py --indir=$indir/$reg/sigrid/$year --outdir=$outdir --MIZ_criteria=$crit --chart_source=AARI
      done
   done
done
