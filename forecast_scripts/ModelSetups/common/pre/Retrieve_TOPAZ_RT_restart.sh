me=`readlink -f $0`
here=`dirname $me`
if [ $# -ne 1 ]
then
   echo `basename $me` rbase
   echo rbase = eg TP4restart2016_101_00
   exit
else
   rbase=$1
fi

# =========================================================================================================
# EMAIL ADDRESS FOR THE WEEKLY UPDATE
email=$(cat $FCemail)
# =========================================================================================================

# vilje info
vdir=/prod/forecast/work/sea/TOPAZ4/Forecast07/
vkey=$HOME/.ssh/key_TWhex2vil #TODO put in hidden
vid=laurentb@vilje.hpc.ntnu.no #TODO put in hidden

# hex info
hdir=/work/timill/TOPAZ_RT
cd $hdir
gfil=$rbase.tar.gz

# cp to migrate
bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts 
yr=${rbase:10:4}
Bdir=$bdir/$yr
mkdir -p $Bdir


if [ -f $Bdir/$gfil ]
then
   echo $rbase already on /migrate
   mail -s "CONFIRMATION: Latest restart on /migrate" $email < $gfil
   exit
else
   # get from vilje
   lst="${rbase}_mem001.a ${rbase}_mem001.b ${rbase}ICE.uf"
   count=0
   for f in $lst
   do
      # echo $f
      scp -i $vkey $vid:$vdir/$f .
      if [ -f $f ]
      then
         count=$((count+1))
      fi
   done

   if [ $count -ne 3 ]
   then
      mail -s "WARNING: Latest restart not on vilje" $email < $gfil
      exit
   else
      echo $gfil > toto
      echo copying to /migrate >> toto

      mail -s "CONFIRMATION: Latest restart on vilje" $email < toto
      rm toto
   fi

   # tar up
   echo tar -zcvf $gfil $lst
   tar -zcvf $gfil $lst

   # copy to /migrate
   echo cp $gfil $Bdir
   cp $gfil $Bdir

   # cp to data x01.1
   ddir=$TP4_REALTIME/expt_01.1/data
   echo cp ${rbase}_mem001.a $ddir/${rbase}.a
   echo cp ${rbase}_mem001.b $ddir/${rbase}.b
   echo cp ${rbase}ICE.uf    $ddir
   cp ${rbase}_mem001.a $ddir/${rbase}.a
   cp ${rbase}_mem001.b $ddir/${rbase}.b
   cp ${rbase}ICE.uf    $ddir
fi
