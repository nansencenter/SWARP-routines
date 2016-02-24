#!/bin/bash
# This script will: 
# 1) Create today's datelist
# 2) Get the restarts - topaz_get_restart.sh
# 3) Prepare the infile.in - make_infile4forecast.sh
# 4) Run the model - pbsjob.sh

# get variables
# ($SWARP_ROUTINES defined in crontab or .bash_profile)
source $SWARP_ROUTINES/source_files/hex_vars.src

print_info=1 # print info to screen (or email in crontab)
test_pre=1
SWARP_PP=1

if [ $# -lt 1 ]
then
   echo "Usage: $0 <<Source file>>"
   exit
else
   # 1st input:
   # different source file for wavewatch3,ice-only,or WAM
   # ! Use thishc.src for hindcast HC runs
   THIS_SRC=`readlink -f $1` #get absolute path
   source $THIS_SRC
fi

xdir=$RTmods/$HYCOMreg/expt_01.$Xno
#if [ $FCtype == "ice_only" ]
#then
#   # always use reanalysis restart
#   restart_OPT=1
#else
#   # else can have option of using ice-only FC
#   restart_OPT=$TP4restart_OPT
#fi

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
# 1. check if hindcast HC is already running
msg=`$qstat | grep ${rungen}x01${Xno}hc`
if [ 1 -eq 1 ]
then
   if [ ${#msg} -ne 0 ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "$FCtype_long is already running - stopping"
         echo "pbs job message:"
         echo $msg
      fi
      exit
   fi
else
   echo "WARNING: already-running check disabled in script"
fi
## ===========================================================


manual=0
if [ $# -eq 2 ]
then
   echo ""
   echo "Running manually (over-riding checks)"
   echo ""
   manual=1
fi


# CREATING THE LOG
logdir=$THISFC/logs
mkdir -p $logdir
log=$logdir/run_hindcast_log.txt
if [ -f "$log" ]
then
   rm $log
fi

#if[ $# -eq 3 ]
#then
   echo ""
   echo "Running in hindcast mode, set cdate to previous Monday"
   echo "Run from Monday to Sunday midnight, +7 days"
   echo "TODO assure that TP4 restart from MET.NO forecast is used"
   echo "" 

#fi
#=======================================================================
# DATE INFO - RUNNING HINDCAST FROM LAST WEEKS MONDAY TO SUNDAY (+7 DAYS)
# SETTING cday TO PREVIOUS WEEK'S MONDAY
#
cday=`date --date="last Monday -7 days" "+%Y%m%d"`
#
#THIS DOES NOT WORK: cday=`date --date="$lastMon -7 days" "+%Y%m%d"`
 

#=======================================================================


##=======================================================================
## DATE INFO - TODAY
##if [ 1 -eq 1 ]
#then
#   cday=$(date +%Y%m%d)
#else
#   # test last Tuesday (restart and dump day are the same
#   # for $TP4restart_OPT=2)
#   cday=$(date --date="last Tuesday" +%Y%m%d)
#   echo " "
#   echo "******************************************************"
#   echo "WARNING: testing artificial start date:" $cday
#   echo "******************************************************"
#   echo " "
#fi
cyear=$(date --date=$cday +%Y)
jday0=10#$(date --date=$cday +%j)
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

# # sts need yesterday's info as well
# yday=`date --date="$cday -1day" +%Y%m%d`
# yyear=`date --date="$yday" +%Y`
# ydj=10#`date --date="$yday" +%j`
# ydj=$((ydj-1))
#=======================================================================

#if [ $restart_OPT -eq 2 ]
#then
#   # ======================================================================
#   # DATE INFO - CURRENT HINDCAST START DAY 
#   # cdate == rdate
   rdate=$(date --date=$cday +%Y%m%d)
   ryear=$(date --date=$rdate +%Y)
   rwday=$(date --date=$rdate +%A)
   rday_j=10#$(date --date=$rdate +%j)
   rday=$((rday_j-1))                        # 1 Jan = day 0
   rday=$(printf '%3.3d' $rday) # 3 digits padding zero
   # USE met.no backup files insterad! rdir=$RTmods/$HYCOMreg/expt_01.1/data
   rname=${rungen}restart${ryear}_${rday}_00
#   if [ $rwday == 'Monday' ]
#   then
#      # yesterday's restart file is topaz reanalysis file
#      restart_OPT=1
#   fi
#   # ======================================================================
#fi
## =========================================================================##

# =========================================================================
echo "Manual = $manual"
if [ $manual -eq 0 ]
then
   # Checks before run:
   # 2. check if forecast has already run
   if [ -f $rundir/final_output/SWARP*.nc ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "$FCtype_long has already run - stopping"
      fi
      exit
   fi
#fi

# TODO check if wave data is avaialable, use only analysis:
#   /work/shared/nersc/msc/WAVES_INPUT/WW3_ARCTIC/2016/analysis_m1/..

#   3. check if waves FC file is there
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
#
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
#
   if [ ! "$HYCOMreg" == "TP4a0.12" ]
   then
      # TODO check output is collected into filename "hindcast" ?!!!
      # 4.a check if outer model TP4 hindcast has run?
      # 4.b check if nesting files are present?
      # outer model?
      if [ ! -f $RTres/TP4a0.12/ice_only/final_output/SWARPiceonly_hindcast_start${rday}T000000Z.nc ]
      then
         if [ $print_info -eq 1 ]
         then
            echo "Ice-only FC has not yet run - stopping"
         fi
         exit
      fi
      ## MOVED BELOW final_day TODO check if last days Nesting files exit? (assume previous exists & saved in nest_out_${HYCOMreg}_hindcast)      
      #if [ ! -f /work/timill/RealTime_Models/${HYCOMreg}/expt_01.$X0/nest_out_${HYCOMreg}_hindcast/nest_${fin_year}_${fin_day_j0}_00.hdr ]
      #then
      #   if [ $print_info -eq 1 ]
      #   then
      #      echo "No nesting hindcast files for end date nest_${fin_year}_${fin_day_j0}_00.hdr - stopping"
      #   fi
      #   exit
      #fi
   fi
fi
# =========================================================================


# =========================================================================
#if [ $restart_OPT -eq 1 ]
#then
   # get reanalysis restart
   # DONT NEED THIS NEED TO BE A RESTART THE SAME DAY AS cday (here current hindcast day)
   #echo "Launching topaz_get_restart_reanalysis.sh @ $(date)"   > $log
   #$FCcommon/pre/topaz_get_restart_reanalysis.sh $THIS_SRC            # get latest restart file

   # GETTING INFO FROM LAST_RESTART
   # cd $rundir  # just in case we've changed dir in script
   #out_restart=info/last_restart.txt
   #rname=$(cat $out_restart)
   #ryear=${rname:10:4}  # year of restart file
   #rday=${rname:15:3}   # julian day of restart file (1 Jan = 0)
   #rdate=`date --date="$ryear-01-01 +${rday}days" +%Y%m%d`
   #r_ndays=10#`date --date="$ryear-12-31" +%j`
#else
# TODO Change from /migrate to /work/shared ? or to /Norstore/ ...?
   #cd $rundir  # just in case we've changed dir in script
   #
if [ "$HYCOMreg" == "TP4a0.12" ] 
then
   echo " USE for now /migrate met.no backup restart files"
   mdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/${ryear} 
   tarfil=$mdir/$rname.tar.gz
   if [ -f $tarfil ]
   then
      echo "Untar file and rename from *mem001.a to *.a"
      echo "$mdir/$tarfil"
      rm $xdir/data/$rname*
      cp $tarfil $xdir/data
      cd $xdir/data
      tar xzf $rname.tar.gz
      mv $rname*.a $rname.a
      mv $rname*.b $rname.b
      rm $rname.tar.gz
      cd $rundir
   fi
   afil=$xdir/data/$rname.a
   bfil=$xdir/data/$rname.b
   ufil=$xdir/data/${rname}ICE.uf
   if [ -f $afil ] && [ -f $bfil ] && [ -f $ufil ]
   then
      echo "Untar in $xdir/data"
      echo "untar & cp $tarfil" > $log
      echo "to  $afile" > $log
      echo "to  $bfile" > $log
      echo "to  $ufile" > $log
      #cp $afil $xdir/data
      #cp $bfil $xdir/data
      #cp $ufil $xdir/data
      #
      #out_restart=$rundir/info/last_restart.txt
      #echo $rname > $out_restart
      #cp $out_restart $logdir
   else
      echo "TP4 Restart met.no files for $rday from ice-only run not present:"
      echo "$tarfil"
      #echo "$bfil"
      #echo "$ufil"
      exit
   fi
else
   # Use restart file from regional (BS1 or FR1) ice_only run
   afil=$xdir/data/$rname.a
   bfil=$xdir/data/$rname.b
   ufil=$xdir/data/${rname}ICE.uf
   if [ -f $afil ] && [ -f $bfil ] && [ -f $ufil ]
   then
      echo "Regional restart files are present in:"
      echo "$xdir/data/"
   else
      echo echo "Regional restart files are NOT present in: EXIT"
      echo "$xdir/data/"
      exit
   fi
fi
#fi
# =========================================================================


# =========================================================================
# DATE INFO - LAST DAY OF FORECAST
# $final_day=julian day relative to $ryear (NB 3 digits)
fin_day=`date --date="$rdate +${FCdays}days" "+%Y%m%d"`
fin_year=$(date --date=$fin_day +%Y)
fin_day_j=10#$(date --date=$fin_day +%j) # 3 digits already - need to convert to base 10
fin_day_j0=$((fin_day_j-1))
fin_day_j0=$(printf '%3.3d' $fin_day_j0) # 3 digits
r_ndays=10#`date --date="$ryear-12-31" +%j`
if [ $cyear -eq $fin_year ]
then
   final_day=$fin_day_j0
else
   final_day=$((r_ndays+fin_day_j0))
   final_day=$(printf '%3.3d' $final_day) # 3 digits
fi
#
#if [ $ryear -eq $yyear ]
#then
#   day2=$(printf '%3.3d' $ydj) # 3 digits
#else
#   day2=$((r_ndays+ydj))
#fi
# =========================================================================

# Check if last days Nesting files exit? (assume previous exists & saved in nest_out_${HYCOMreg}_hindcast)      
if [ ! "$HYCOMreg" == "TP4a0.12" ]
then
   if [ ! -f /work/timill/RealTime_Models/${HYCOMreg}/expt_01.$Xno/nest_out_${HYCOMreg}_hindcast/nest_${fin_year}_${fin_day_j0}_00.hdr ]
   then
      if [ $print_info -eq 1 ]
      then
         echo "No nesting hindcast files for end date nest_${fin_year}_${fin_day_j0}_00.hdr - stopping"
      fi 
      exit
   fi
fi
# =========================================================================
# MAKE INFILE 

# print to screen
echo ""
echo "Restart date   : $rdate"
echo "Hindcast date  : $cday"
echo "Final date     : $fin_day"
echo "Final hour     : $FCfinal_hour"
echo "Restart files of $rname"
echo "Forecast final day ${fin_year}_${fin_day_j0}_$FCfinal_hour (${ryear}_${final_day}_$FCfinal_hour)"

#if [ $FCtype == "ice_only" ] && [ $TP4restart_OPT -eq 2 ] && [ $day2 -ne $rday ]
#then
#   # echo "Launching make_infile_3days.sh @ $(date)"                  >> $log
#   # inputs="$THIS_SRC $rungen $ryear $rday $day2 $final_day $FCfinal_hour $nesting_outer $nesting_inner"
#   # echo "$FCcommon/make_infile_3days.sh $inputs"
#   # $FCcommon/pre/make_infile_3days.sh $inputs
#   days="$rday $day2 $final_day"
#   Ropts="F T F"
#else
#   # echo "Launching make_infile_2days.sh @ $(date)"                  >> $log
#   # inputs="$THIS_SRC $rungen $ryear $rday $final_day $FCfinal_hour $nesting_outer $nesting_inner"
#   # echo "$FCcommon/make_infile_2days.sh $inputs"
#   # $FCcommon/pre/make_infile_2days.sh $inputs
#   days="$rday $final_day"
#   Ropts="F T"
#fi

# If hindcast save last days restart
days="$rday $final_day"
Ropts="F T"
# Hindcast run from cday to final_day

ftmp='tmp.txt'
echo "rungen         $rungen"             >  $ftmp
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
# =========================================================================

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
# Get other inputs (!!for hindcast use blkdat.hindcast.input)
cp $THISFC/inputs/blkdat.hindcast.input      $xdir/blkdat.input
cp $FCcommon/inputs/preprocess.sh   $xdir
cp $FCcommon/inputs/postprocess.sh  $xdir

JF=$xdir/pbsjob.sh
Nproc=133 # no of CPUs

if [ $SWARP_PP -eq 1 ]
then
   # set up post-processing
   cat $FCcommon/inputs/pbsjob.sh | sed \
    -e "s/JOBNAME/$JOBNAME/g" \
    -e "s/MPPWIDTH/$Nproc/g" \
    -e "s/exit \$?//g" \
    > $JF

   PP=$FCcommon/post/process_HCresults.sh
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
   echo "cp $THISFC/inputs/blkdat.hindcast.input      $xdir/blkdat.input"
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
# TODO Do we need this?? clean data directory before run 
# rm -f $xdir/data/${rungen}DAILY*
# rm -f $xdir/data/${rungen}archv*

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
