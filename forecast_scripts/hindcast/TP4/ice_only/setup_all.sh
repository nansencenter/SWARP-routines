# script to copy all restart files to $TP4_REALTIME/expt_01.4/data
# - other expt's will link to them
source $SWARP_ROUTINES/source_files/hex_vars.src
NE=$SWARP_ROUTINES/forecast_scripts/utility_scripts/NewExp.sh
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/TP4
Done=0

P=$TP4_REALTIME
cd $P
refno=4
xref=$P/expt_01.$refno           # reference expt directory  (has restarts)
bref=$P/Build_V2.2.12_X01.$refno # reference Build directory (has hycom executable compiled)

if [ 1 -eq 1 ]
then
   rno0=0      # start at the beginning 1=$((rno0+1))
   do_num=1000 # stop after this number done
else
   rno0=30     # start at no $((rno0+1))
   do_num=1 # stop after this number done
fi
E0=1$refno # increment this
X=01.$refno
Eno=0$E0
rno=0

md=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts/2015/
list=`ls $md/*gz`

logdir=$hd/ice_only/logs
mkdir -p $logdir
log=$logdir/hc.log
rm -f $log
touch $log

for tfil0 in $list
do
   Eno=0$E0
   X=0${E0:0:1}.${E0:1:1}

   if [ $rno -gt $rno0 ]
   then
      tfil=`basename $tfil0`
      ryear=${tfil:10:4}
      rday=${tfil:15:3}
      rday2=10#$rday #base 10
      fday2=$((rday2+7))
      fday=`printf %3.3d $fday2`
      echo $tfil $ryear $rday $fday
      # echo $X; exit # ????

      # ECMWF files in 2015 skip between
      # Tuesday 22-Aug 12:00 and Wednesday 23-Aug 06:00
      # ie 4 missing records
      ecmwf_missing=20150822
      jday_missing=10#`date --date=$ecmwf_missing +%j`
      jday_missing=$((jday_missing-1))
      if [ $jday_missing -gt $((rday2+0)) ] && [ $jday_missing -lt $fday2  ]
      then
         echo "can't run week ${ryear}_${rday}"
         echo "- no ECMWF forcing for $ecmwf_missing (${ryear}_${jday_missing})"
      else

         $NE 01$refno $Eno          # run NewExp.sh script
         xdir=$P/expt_$X            # reference expt directory  (has restarts)
         bdir=$P/Build_V2.2.12_X$X  # reference Build directory (has hycom executable compiled)

         cd $bdir
         ln -s $bref/hycom .

         cd $xdir
         mkdir -p log
         ln -s $xref/data/TP4restart${ryear}_${rday}* data  # link to restart 

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
      fi # check forcing is there

   fi # check not in $xref

   # increment rno,E0
   rno=$((rno+1))
   E0=$((E0+1))

done
