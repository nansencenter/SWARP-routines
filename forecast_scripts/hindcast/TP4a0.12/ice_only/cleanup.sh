source $SWARP_ROUTINES/source_files/hex_vars.src
refno=5

E0_start=$((1+1$refno))
E0_stop=99
P=$TP4_REALTIME

for E0 in `seq $E0_start $E0_stop`
do
   E=0$E0
   X=0${E0:0:1}.${E0:1:1}
   xdir=$P/expt_$X            # new expt directory  (has restarts)
   bdir=$P/Build_V2.2.12_X$X  # new Build directory (has hycom executable compiled)
   # echo $E $X; continue

   if [ ! -d $xdir ] && [ ! -d $bdir ]
   then
      echo "($xdir,$bdir) not present"
      continue
   fi
   # echo $E $X; continue

   # delete expt and build dir's
   echo rm -rf $xdir $bdir
   rm -rf $xdir $bdir
   # exit

   # clean force etc
   cd $P/force/rivers
   rm $E
   pwd
   ls -lh
   echo " "

   cd $P/force/nersc_era40
   rm $E
   pwd
   ls -lh
   echo " "

   cd $P/relax
   rm $E
   pwd
   ls -lh
   echo " "

   cd $P/nest_nersc
   rm $E
   pwd
   ls -lh
   echo " "

   cd $P/tides_nersc
   rm $E
   pwd
   ls -lh
   echo " "

done
