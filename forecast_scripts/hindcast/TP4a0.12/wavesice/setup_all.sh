# script to set up weeks after 2015_060 based on expt 1.4
# - links to executable from that expt
# - links to restarts from that expt
# - makes pbsjob.sh
# - makes blkdat.input
# - makes infile.in
# - DOESN'T GET archv.extract - USE archv FILES FROM ice_only
source $SWARP_ROUTINES/source_files/hex_vars.src
NE=$SWARP_ROUTINES/forecast_scripts/utility_scripts/NewExp.sh
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/TP4a0.12
HCtype=wavesice
Done=0
FULL_RESET=1

P=$TP4_REALTIME
cd $P
refno=4
xref=$P/expt_01.$refno           # reference expt directory  (has restarts)
bref=$P/Build_V2.2.12_X01.$refno # reference Build directory (has hycom executable compiled)

E0=$((1$refno-1)) # increment this
X=01.$refno
Eno=0$E0

md=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/
list=(`ls $md/2015/*gz`)
list2=(`ls $md/2016/*gz`)
Nlist=${#list[@]}
inext=Nlist
for loop_i in `seq 1 ${#list2[@]}`
do
   list[$inext]=${list2[$((loop_i-1))]}
   inext=$((inext+1))
   Nlist=$((Nlist+1))
done

logdir=$hd/$HCtype/logs
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


if [ 1 -eq 1 ]
then
   # DO ALL
   E0_start=15 # do 1st (14) manually
   E0_stop=2000
else
   # DO SOME
   # E0_start=36
   # E0_stop=56
   E0_start=57
   E0_stop=2000
fi


for rno in `seq 1 $Nlist`
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

   # echo `basename $tfil0` $E0
   # continue

   if [ $E0 -lt $E0_start ] || [ $E0 -gt $E0_stop ]
   then
      echo "$E0: skipping"
      continue
   else
      echo "$E0: TODO"
   fi

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

   # use same executable as expt 1.4
   cd $bdir
   echo "rm hycom; ln -s $bref/hycom ."
   rm hycom
   ln -s $bref/hycom .

   cd $xdir/data
   mkdir -p log
   echo "rm TP4restart*; ln -s $xref/data/TP4restart${ryear}_${rday}* ."  # link to restart 
   rm TP4restart*
   ln -s $xref/data/TP4restart${ryear}_${rday}* .  # link to restart 
   cd $xdir

   id=$hd/$HCtype/inputs
   # # archv.extract
   # cp $id/archv.extract data

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
   nesting_outer="F"
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
   echo "$python $FCcommon/pre/make_infile.py --infile=$ftmp"
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


   ## ###################################################################
   ## NO NESTING
   ## #setup nesting
   ## rm -f nesting.in
   ## cd $P
   ## # use non-interactive version of nest_nersc/bin/nest_outer.sh
   ## rm -f nest_nersc/$Eno/outer/SCRATCH/*

   ## $hd/$HCtype/nest_outer.sh $X $FR1_REALTIME 01.0
   ## rm -f nest_nersc/$Eno/outer/SCRATCH/*

   ## $hd/ice_only/nest_outer.sh $X $BS1_REALTIME 01.0
   ## rm -f nest_nersc/$Eno/outer/SCRATCH/*
   ## ###################################################################

   Done=$((Done+1))
done
