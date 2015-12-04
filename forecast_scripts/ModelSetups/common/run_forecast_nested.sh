#!/bin/bash
# This script will: 
# 1) Create today's datelist
# 2) Get the restarts - topaz_get_restart.sh
# 3) Prepare the infile.in - make_infile4forecast.sh
# 4) Run the model - pbsjob.sh

# get variables
# ($SWARP_ROUTINES defined in crontab or .bash_profile)
source $SWARP_ROUTINES/source_files/hex_vars.src

print_info=0 # print info to screen (or email in crontab)
test_pre=0
SWARP_PP=1

if [ $# -lt 1 ]
then
   echo "Usage: $0 <<Source file>>"
   exit
else
   THIS_SRC=`readlink -f $1` #get absolute path
   source $THIS_SRC
fi

manual=0
if [ $# -eq 2 ]
then
   echo ""
   echo "Running manually (over-riding checks)"
   echo ""
   manual=1
fi

if [ $print_info -eq 1 ]
then
   echo "Source file    : $THIS_SRC"
   echo "HYCOM region   : $HYCOMreg"
   echo "FC type        : $FCtype"
fi

xdir=$RTmods/$HYCOMreg/expt_01.$Xno
if [ $FCtype == "ice_only" ]
then
   # always use reanalysis restart
   restart_OPT=1
else
   # else can have option of using ice-only FC
   restart_OPT=$TP4restart_OPT
fi

if [ $print_info -eq 1 ]
then
   echo ""
   echo "Source file    : $THIS_SRC"
   echo "HYCOM region   : $HYCOMreg"
   echo "FC type        : $FCtype"
   echo "Expt no        : $Xno"
   echo "PBS job name   : $JOBNAME"
   echo "expt dir       : $xdir"
   echo "Results to     : $THISFC2"
   echo "Email          : $FCemail"
   echo ""
fi


## ===========================================================
# Initial check before run
# (do this before datelist.txt is changed)
# 1. check if forecast is already running
msg=`$qstat | grep ${rungen}x01${Xno}fc`
if [ ! -z "${msg}" ] && [ $manual -eq 0 ]
then
   if [ $print_info -eq 1 ]
   then
      echo "$FCtype_long is already running - stopping"
      echo "pbs job message:"
      echo $msg
   fi
   exit
fi
## ===========================================================


# CREATING THE LOG
logdir=$THISFC/logs
mkdir -p $logdir
log=$logdir/run_forecast_log.txt
if [ -f "$log" ]
then
   rm $log
fi


#=======================================================================
# DATE INFO - TODAY
if [ 1 -eq 1 ]
then
   cday=$(date +%Y%m%d)
else
   # test last Tuesday (restart and dump day are the same
   # for $TP4restart_OPT=2)
   cday=$(date --date="last Tuesday" +%Y%m%d)
   echo " "
   echo "******************************************************"
   echo "WARNING: testing artificial start date:" $cday
   echo "******************************************************"
   echo " "
fi
cyear=$(date --date=$cday +%Y)
jday0=$(date --date=$cday +%j)
jday_today0=$(( $jday0 - 1))              # julian day of TOPAZ (0=1st Jan)
jday_today=$(printf '%3.3d' $jday_today0) # 3 digits

rundir=$THISFC2/$cday                     # where the results will end up
mkdir -p $rundir
mkdir -p $rundir/info

# CREATING THE DATELIST
datelist=$rundir/info/datelist.txt
echo $cday                             >  $datelist
echo $(date --date=$cday +%Y-%m-%d)    >> $datelist
echo $cyear                            >> $datelist
echo $(date --date=$cday +%m)          >> $datelist
echo $(date --date=$cday +%d)          >> $datelist
echo $jday_today                       >> $datelist
cp $datelist $logdir

if [ $print_info -eq 1 ]
then
   echo Datelist:
   echo $datelist
   echo Contents:
   cat $datelist
   echo " "
fi

# sts need yesterday's info as well
yday=`date --date="$cday -1day" +%Y%m%d`
yyear=`date --date="$yday" +%Y`
ydj=`date --date="$yday" +%j`
ydj=$((ydj-1))
#=======================================================================

if [ $restart_OPT -eq 2 ]
then
   # ======================================================================
   # DATE INFO - YESTERDAY 
   # - need to get restart for yesterday
   rdate=$(date --date="$cday -1day" +%Y%m%d)
   ryear=$(date --date=$rdate +%Y)
   rwday=$(date --date=$rdate +%A)
   rday_j=$(date --date=$rdate +%j)
   rday=$((rday_j-1))                        # 1 Jan = day 0
   rdir=$RTmods/$HYCOMreg/expt_01.1/data
   rname=${rungen}restart${ryear}_${rday}_00
   if [ $rwday == 'Monday' ]
   then
      # yesterday's restart file is topaz reanalysis file
      restart_OPT=1
   fi
   # ======================================================================
fi
# =========================================================================##


# =========================================================================
if [ $manual -eq 0 ]
then
   # Checks before run:
   # 2. check if TP4 ice-only forecast has run
   #    - needed for nesting
   oTP4=$RTres/TP4a0.12/ice_only/$cday/final_output/SWARPiceonly_forecast_start${cday}T000000Z.nc
   ls $oTP4
   if [ ! -f $oTP4 ]
   then
      echo hi
      if [ $print_info -eq 1 ]
      then
         echo "TP4 ice-only FC has not yet run - stopping"
      fi
      exit
   fi

   # 3. check if forecast has already run
   if [ -f $rundir/final_output/SWARP*.nc ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "$FCtype_long has already run - stopping"
      fi
      exit
   fi

   # 4. check if waves FC file is there
   do_check=0
   if [ $FCtype == "wavesice" ]
   then
      fc_fil=$wamnsea/$cyear/forecasts/wam_nsea.fc.$cday.nc
      do_check=1
   elif [ $FCtype == "wavesice_ww3arctic" ]
   then
      fc_fil=$ww3_arctic/$yyear/forecast/SWARP_WW3_ARCTIC-12K_$yday.fc.nc
      do_check=1
   fi

   if [ ! -f $fc_fil ] && [ $do_check -eq 1 ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "$FCtype_long file is not yet present - stopping"
         echo $fc_fil
         date
      fi
      exit
   fi

   if [ $restart_OPT -eq 2 ] && [ ! $FCtype == "ice_only" ]
   then
      # 5. check if local ice-only forecast has run
      if [ ! -f $RTres/$HYCOMreg/ice_only/final_output/SWARP*.nc ]
      then
         if [ $print_info -eq 1 ]
         then
            echo "Ice-only FC has not yet run - stopping"
         fi
         exit
      fi
   fi
fi
# =========================================================================


# =========================================================================
if [ $restart_OPT -eq 1 ]
then
   # get reanalysis restart
   echo "Launching topaz_curviint_restart_reanalysis.sh @ $(date)"   > $log
   $FCcommon/pre/topaz_curviint_restart_reanalysis.sh $THIS_SRC            # get latest restart file

   # GETTING INFO FROM LAST_RESTART
   cd $rundir  # just in case we've changed dir in script
   out_restart=info/last_restart.txt
   rname=$(cat $out_restart)
   ryear=${rname:10:4}  # year of restart file
   rday=${rname:15:3}   # julian day of restart file (1 Jan = 0)
   rdate=`date --date="$ryear-01-01 +${rday}days" +%Y%m%d`
else
   afil=$rdir/$rname.a
   bfil=$rdir/$rname.b
   ufil=$rdir/${rname}ICE.uf
   if [ -f $afil ] && [ -f $bfil ] && [ -f $ufil ]
   then
      echo "cp $afil $xdir/data" > $log
      echo "cp $bfil $xdir/data" > $log
      echo "cp $ufil $xdir/data" > $log
      cp $afil $xdir/data
      cp $bfil $xdir/data
      cp $ufil $xdir/data
      #
      out_restart=$rundir/info/last_restart.txt
      echo $rname > $out_restart
      cp $out_restart $logdir
   else
      echo "Yesterday's restart files from ice-only run not present:"
      echo "$afil"
      echo "$bfil"
      echo "$ufil"
      exit
   fi
fi
# =========================================================================


# =========================================================================
# DATE INFO - LAST DAY OF FORECAST
# $final_day=julian day relative to $ryear (NB 3 digits)
fin_day=`date --date="$cday +${FCdays}days" "+%Y%m%d"`
fin_year=$(date --date=$fin_day +%Y)
fin_day_j=$(date --date=$fin_day +%j)
fin_day_j0=$((fin_day_j-1))
if [ $ryear -eq $fin_year ]
then
   final_day=$(printf '%3.3d' $fin_day_j0) # 3 digits
else
   final_day=$((r_ndays+fin_day_j0))
fi

if [ $ryear -eq $yyear ]
then
   day2=$(printf '%3.3d' $ydj) # 3 digits
else
   day2=$((r_ndays+ydj))
fi
# =========================================================================


# =========================================================================
# MAKE INFILE 

# =========================================================================
# want to save next Monday's restart
Nday=`date --date="next Monday" +%Y%m%d`
Nyear=$(date --date=$Nday +%Y)
Nday_j=$(date --date=$Nday +%j)
Nday_j=$((Nday_j-1))
if [ $ryear -eq $Nyear ]
then
   Rday=$(printf '%3.3d' $Nday_j) # 3 digits
else
   Rday=$((r_ndays+Nday_j))
fi
# =========================================================================

# print to screen
echo "Restart date   : $rdate"
echo "Current date   : $cday"
echo "Final date     : $fin_day"
echo "Final hour     : $FCfinal_hour"
echo "Restart files of $rname"
echo "Forecast final day ${fin_year}_$(printf '%3.3d' $fin_day_j0)_$FCfinal_hour (${ryear}_${final_day}_$FCfinal_hour)"

if [ $FCtype == "ice_only" ] && [ $TP4restart_OPT -eq 2 ] && [ $day2 -ne $rday ]
then
   if [ $Rday -ge $final_day ]
   then
      days="$rday $day2 $final_day"
      Ropts="F T F"
   else
      days="$rday $day2 $Rday $final_day"
      Ropts="F F T F"
   fi
else
   if [ $Rday -ge $final_day ]
   then
      days="$rday $final_day"
      Ropts="F F"
   else
      days="$rday $Rday $final_day"
      Ropts="F T F"
   fi
fi

ftmp='tmp.txt'
echo "rungen         $rungen"             >> $ftmp
echo "expt_dir       $xdir"               >> $ftmp
echo "refyear        $ryear"              >> $ftmp
echo "days           $days"               >> $ftmp
echo "restart_opts   $Ropts"              >> $ftmp
echo "final_hour     $FCfinal_hour"       >> $ftmp
echo "nest_outer     $nesting_outer"      >> $ftmp
echo "nest_inner     $nesting_inner"      >> $ftmp


echo "Launching make_infile.py @ $(date)" >> $log
echo ""                                   >> $log
cat $ftmp                                 >> $log
echo ""                                   >> $log

# load python and launch make_infile.py
[ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
module load python/2.7.9-dso
$python $FCcommon/pre/make_infile.py --infile=$ftmp
rm $ftmp

infile=$xdir/infile.in
if [ $print_info -eq 1 ]
then
   echo ""
   echo "Infile: $infile"
   echo ""
   cat $infile
   echo ""
fi
cp $infile $rundir/info
# =========================================================================


# =========================================================================
# Get other inputs
cp $THISFC/inputs/blkdat.input      $xdir
cp $FCcommon/inputs/preprocess.sh   $xdir
cp $FCcommon/inputs/postprocess.sh  $xdir

JF=$xdir/pbsjob.sh
if [ $rungen == "BS1" ]
then
   Nproc=112 # no of CPUs
elif [ $rungen == "FR1" ]
then
   Nproc=51 # no of CPUs
fi

if [ $SWARP_PP -eq 1 ]
then
   # set up post-processing
   cat $FCcommon/inputs/pbsjob.sh | sed \
    -e "s/JOBNAME/$JOBNAME/g" \
    -e "s/MPPWIDTH/$Nproc/g" \
    -e "s/exit \$?//g" \
    > $JF

   PP=$FCcommon/post/process_FCresults.sh
   echo "# SWARP post-processing:"  >> $JF
   echo "PP=$PP"                    >> $JF
   echo "THIS_SRC=$THIS_SRC"        >> $JF
   echo '$PP $THIS_SRC'             >> $JF
   echo " "                         >> $JF
   echo " "                         >> $JF
   echo 'exit $?'                   >> $JF
else
   cat $FCcommon/inputs/pbsjob.sh | sed \
    -e "s/JOBNAME/$JOBNAME/g" \
    -e "s/MPPWIDTH/$Nproc/g" \
    > $JF
fi


if [ $print_info -eq 1 ]
then
   echo "cp $THISFC/inputs/blkdat.input      $xdir"
   echo "cp $FCcommon/inputs/preprocess.sh   $xdir"
   echo "cp $FCcommon/inputs/postprocess.sh  $xdir"
fi

#choose variables to extract (-DARCHIVE_SELECT)
exfil=$THISFC/inputs/archv.extract
if [ -f $exfil ]
then
   cp $exfil $xdir/data
   if [ $print_info -eq 1 ]
   then
      echo cp $exfil $xdir/data
   fi
fi
# =========================================================================


# =========================================================================
# clean data directory before run
rm -f $xdir/data/${rungen}DAILY*
rm -f $xdir/data/${rungen}archv*

# clean log file - else mpijob.out gets too big
rm -f $xdir/log/*

if [ ! $test_pre -eq 1 ]
then
   # launch job
   echo "Launching pbsjob @ $(date)"   >> $log
   cd $xdir
   $qsub pbsjob.sh
fi
# =========================================================================
