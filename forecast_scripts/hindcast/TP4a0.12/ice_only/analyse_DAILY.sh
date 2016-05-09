# analyse all hindcasts:
source $SWARP_ROUTINES/source_files/hex_vars.src
myfile=`readlink -f $0`
pydir=`dirname $myfile`

R=TP4a0.12
rungen=${R:0:3}
HCtype=ice_only
HCtype2=$HCtype
#
print_info=1
#
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/$R
user=timill #TODO source hidden
hd2=$hd/$HCtype
outdir0=$TP4_REALTIME/../results_hindcasts/$R/$HCtype

# bug or bugfix version
outdir0=$outdir0/bugfix2
# outdir0=$outdir0/bugfix
# outdir0=$outdir0/bug
# ls $outdir0

# time span
# YEARS="2015 2016"
YEARS=2015

###########################################
# Make one directory for all binary files:
# - then link to weekly directories
outdirA=$outdir0/analysis
mkdir -p $outdirA

outdir1=$outdirA/binaries_all/
outdirD=$outdir1/DAILY

redo_links=0
if [ ! -d $outdirD ] || [ $redo_links -eq 1 ]
then
   # redo linking
   rm -rf $outdir1
   mkdir -p $outdir1
   mkdir $outdirD


   for YEAR in $YEARS
   do
      outdir=$outdir0/${YEAR}
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

# plotting=False
plotting=True # make plots
cdate0=20150101 # don't process files from before this date

for afil in $rungen*.a
do
   yr=${afil:18:4}
   jday=$((0+10#${afil:23:3}))
   cdate=`date --date="${yr}-01-01 +${jday}days" +%Y%m%d`

   if [ 1 -eq 1 ]
   then
      echo $afil
      echo $yr
      echo $cdate
      # exit
   fi

   if [ $cdate -lt $cdate0 ]
   then
      echo jumping to $cdate0
      continue
   elif [[ ! $YEARS == *$yr* ]]
   then
      echo "Not processing $yr"
      continue
   fi

   datedir=$outdirA/DAILY
   mkdir -p $datedir
   datedir=$datedir/${cdate}
   # if [ 1 -eq 1 ]
   if [ ! -d $datedir ]
   then
      # only run analysis if not already done
      mkdir -p $datedir
      echo "$python $pydir/analyse_DAILY.py --infile=$afil --outdir=$datedir --plotting=$plotting"
      $python $pydir/analyse_DAILY.py --infile=$afil --outdir=$datedir --plotting=$plotting
   else
      echo "Skipping $datedir - already done"
   fi
done
