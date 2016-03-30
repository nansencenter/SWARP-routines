# compile all FC models
me=`readlink -f $0`
here=`dirname $me`

source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -eq 0 ]
then
   echo Usage: $0 restore
   echo "*restore=1: do want to restore Build directories"
   echo "(also rerun nesting if using TP4 domain)"
   echo "*restore=0: don't want to restore Build directories"
   exit
else
   restore=$1
fi


# restore Build directory if needed
if [ $restore -eq 1 ]
then
   # restore deleted files in Build
   # - this also gets flags from the proper "inputs"
   RBS=$SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_build.sh
   for reg in TP4 BS1 FR1
   do
      for FCtype in ice_only wavesice_ww3arctic
      do
         echo $RBS $reg $FCtype
         $RBS $reg $FCtype   || echo "Build set-up error: $reg $FCtype"
      done
   done
fi


# recompile
logdir=$here/logs
mkdir -p $logdir
logfile=$logdir/compile.log
rm -r $logfile
touch $logfile

for Xno in 1 3
do
   for rdir in $TP4_REALTIME
   do
      cd $rdir/Build_V2.2.12_X01.$Xno
      pwd                                       >> $logfile
      echo " "                                  >> $logfile
      python $here/print_flags.py               >> $logfile
      echo " "                                  >> $logfile
      make clean
      make || echo "!!**COMPILATION ERROR**!!"  >> $logfile
      echo " "                                  >> $logfile
   done
done

for Xno in 0 1
do
   for rdir in $BS1_REALTIME $FR1_REALTIME
   do
      cd $rdir/Build_V2.2.12_X01.$Xno
      pwd                                       >> $logfile
      echo " "                                  >> $logfile
      python $here/print_flags.py               >> $logfile
      echo " "                                  >> $logfile
      make clean
      make || echo "!!**COMPILATION ERROR**!!"  >> $logfile
      echo " "                                  >> $logfile
   done
done
