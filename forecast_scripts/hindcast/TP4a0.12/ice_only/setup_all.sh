# script to copy all restart files to $TP4_REALTIME/expt_01.4/data
# - other expt's will link to them
source $SWARP_ROUTINES/source_files/hex_vars.src
NE=$SWARP_ROUTINES/forecast_scripts/utility_scripts/NewExp.sh
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/TP4a0.12
Done=0
FULL_RESET=0

P=$TP4_REALTIME
cd $P
refno=4
xref=$P/expt_01.$refno           # reference expt directory  (has restarts)
bref=$P/Build_V2.2.12_X01.$refno # reference Build directory (has hycom executable compiled)

if [ 1 -eq 1 ]
then
   rno0=1      # start at the beginning 1 (ref = 0: do that manually)
   do_num=1000 # stop after this number done
else
   rno0=30     # start at no $((rno0+1))
   do_num=1 # stop after this number done
fi
E0=1$refno # increment this
X=01.$refno
Eno=0$E0

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

logdir=$hd/ice_only/logs
mkdir -p $logdir
log=$logdir/hc.log
rm -f $log
touch $log

#########################################################
# ECMWF files in 2015 skip between
# Saturday 22-Aug 12:00 and Sunday 23-Aug 06:00
# ie 4 missing records
bad_dates[0]=20150817 #Monday before this
bad_error[0]="Gap in ECMWF forcing 20150822-23"

# don't want to do last Monday's version 
# - winds are forecasts part of the week
bad_dates[1]="`date --date="last Monday" +%Y%m%d`"
bad_error[1]="Last Monday"

Nbad=${#bad_dates[@]}
# for loop_i in `seq 1 $Nbad`
# do
#    echo ${bad_dates[$((loop_i-1))]}
# done
#########################################################

for rno in `seq $((rno0+1)) $Nlist`
do
   tfil0=${list[$((rno-1))]}
   E0=$((E0+1))
   Eno=0$E0
   X=0${E0:0:1}.${E0:1:1}

   tfil=`basename $tfil0`
   ryear=${tfil:10:4}
   rday=${tfil:15:3}
   rday2=10#$rday #base 10
   rday2=$((rday2+0))
   fday2=$((rday2+7))
   fday=`printf %3.3d $fday2`

   dref=$ryear-01-01
   rdate=`date --date="$dref +${rday2}days" +%Y%m%d`
   echo " "
   echo $tfil $ryear $rday $fday $rdate

   BAD_DAY=0
   for nbad in `seq 0 $((Nbad-1))`
   do
      if [ $rdate -eq ${bad_dates[$nbad]} ]
      then
         BAD_DAY=1
         BAD_ERROR=${bad_error[$nbad]}
         break
      fi
   done

   if [ $BAD_DAY -eq 1 ]
   then
      echo " "
      echo "Not setting up week ${ryear}_${rday} ($rdate)"
      echo $BAD_ERROR
      echo " "
      continue
   fi

   xdir=$P/expt_$X            # new expt directory  (has restarts)
   bdir=$P/Build_V2.2.12_X$X  # new Build directory (has hycom executable compiled)
   if [ ! -d $xdir ] || [ ! -d $bdir ]
   then
      echo "Setting up `basename $xdir` `basename $bdir`"
      $NE 01$refno $Eno          # run NewExp.sh script
   elif [ $FULL_RESET -eq 0 ]
   then
      echo " "
      echo "Skipping `basename $xdir` `basename $bdir`"
      echo " "
      continue
   fi

   cd $bdir
   echo ln -s $bref/hycom .
   ln -s $bref/hycom .

   cd $xdir/data
   mkdir -p log
   echo ln -s $xref/data/TP4restart${ryear}_${rday}* .  # link to restart 
   ln -s $xref/data/TP4restart${ryear}_${rday}* .  # link to restart 
   cd $xdir

   # archv.extract
   id=$hd/ice_only/inputs
   cp $id/archv.extract data

   # blkdat.input
   cat $id/blkdat.input  | sed \
    -e "s/IEX/$Eno/g" \
    > blkdat.input

   #pbsjob.sh
   JOBNAME="TP4X${Eno}HC"
   cat $id/pbsjob.sh | sed \
    -e "s/JOBNAME/$JOBNAME/g" \
     > pbsjob.sh

   ###################################################################
   # make infile
   rungen=TP4
   days="$rday $fday"
   Ropts="F F"
   FCfinal_hour="00"
   nesting_outer="T"
   nesting_inner="F"

   ftmp='tmp.txt'
   echo "rungen         $rungen"             >> $ftmp
   echo "expt_dir       $xdir"               >> $ftmp
   echo "refyear        $ryear"              >> $ftmp
   echo "days           $days"               >> $ftmp
   echo "restart_opts   $Ropts"              >> $ftmp
   echo "final_hour     $FCfinal_hour"       >> $ftmp
   echo "nest_outer     $nesting_outer"      >> $ftmp
   echo "nest_inner     $nesting_inner"      >> $ftmp

   # load python and launch make_infile.py
   [ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
   module load python/2.7.9-dso
   $python $FCcommon/pre/make_infile.py --infile=$ftmp
   rm $ftmp

   # change forcing
   ftmp=infile.in.tmp
   mv infile.in $ftmp
   cat $ftmp | sed \
    -e "s/ecncF/ecnc2/g" \
     > infile.in
   rm $ftmp
   ###################################################################


   ###################################################################
   #setup nesting
   rm -f nesting.in
   cd $P
   # use non-interactive version of nest_nersc/bin/nest_outer.sh
   rm -f nest_nersc/$Eno/outer/SCRATCH/*

   $hd/ice_only/nest_outer.sh $X $FR1_REALTIME 01.0
   rm -f nest_nersc/$Eno/outer/SCRATCH/*

   $hd/ice_only/nest_outer.sh $X $BS1_REALTIME 01.0
   rm -f nest_nersc/$Eno/outer/SCRATCH/*
   ###################################################################

   Done=$((Done+1))
   if [ $Done -eq $do_num ]
   then
      exit
   fi

done
