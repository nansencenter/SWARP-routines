# analyse all hindcasts:
source $SWARP_ROUTINES/source_files/hex_vars.src
myfile=`readlink -f $0`
pydir=`dirname $myfile`

R=TP4a0.12
rungen=${R:0:3}
HCtype=wavesice
HCtype2=$HCtype
#
print_info=1
#
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/$R
user=timill #TODO source hidden
hd2=$hd/$HCtype
outdir0=$TP4_REALTIME/../results_hindcasts/$R/$HCtype

###########################################
# Make one directory for all binary files:
# - then link to weekly directories
outdirA=$outdir0/analysis
mkdir -p $outdirA

outdir1=$outdirA/binaries_all/
outdirD=$outdir1/DAILY

if [ 1 -eq 0 ]
then
   # redo linking
   rm -rf $outdir1
   mkdir -p $outdir1
   mkdir $outdirD


   for YEAR in 2015 2016
   do
      outdir=$outdir0/${YEAR}_GOOD
      cd $outdir

      for subdir in ${YEAR}_*
      do
         cd $outdir
         echo $subdir
         SUBDIR=`readlink -f $subdir` # full path
         echo $SUBDIR


         cd $outdirD
         echo ''
         echo '*********************************************'
         pwd
         echo ''
         echo 'ln -s $SUBDIR/binaries/DAILY/$rungen*.a .'
         ln -s $SUBDIR/binaries/DAILY/$rungen*.a .
         echo ' '
         echo 'ln -s $SUBDIR/binaries/DAILY/$rungen*.b .'
         ln -s $SUBDIR/binaries/DAILY/$rungen*.b .
         echo '*********************************************'
         echo ' '

      done
   done
fi
###########################################

# analyse ice edge distance in DAILY files
# could also look at MIZ width (fice,hice)
cd $outdirD
plotting=False

for afil in $rungen*.a
do
   yr=${afil:18:4}
   jday=$((0+10#${afil:23:3}))
   pyr=$((yr-1))
   cdate=`date --date="${pyr}-12-31 +${jday}days" +%Y%m%d`

   # plotting=True
   cdate0=20150701
   if [ $cdate -lt $cdate0 ]
   then
      echo jumping to $cdate0
      continue
   fi

   datedir=$outdirA/DAILY
   mkdir -p $datedir
   datedir=$datedir/${cdate}
   if [ ! -d $datedir ]
   then
      # only run analysis if not already done
      mkdir $datedir
      echo "$python $pydir/analyse_DAILY.py --infile=$afil --outdir=$datedir --plotting=$plotting"
      $python $pydir/analyse_DAILY.py --infile=$afil --outdir=$datedir --plotting=$plotting
   else
      echo "Skipping $datedir - already done"
   fi
done
