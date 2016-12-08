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

# vilje info
vdir=/prod/forecast/work/sea/TOPAZ4/Forecast07/
vkey=$HOME/.ssh/key_TWhex2vil #TODO put in hidden
vid=laurentb@vilje.hpc.ntnu.no #TODO put in hidden

# hex info
hdir=/work/timill/TOPAZ_RT
cd $hdir
gfil=$rbase.tar.gz

# scp from vilje and tar up files
if [ 1 -eq 1 ]
then
   lst="${rbase}_mem001.a ${rbase}_mem001.b ${rbase}ICE.uf"
   for f in $lst
   do
      # echo $f
      scp -i $vkey $vid:$vdir/$f .
   done
   tar -zcvf $gfil $lst
fi

bdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts 
yr=${rbase:10:4}

# cp to data x01.1
ddir=$TP4_REALTIME/expt_01.1/data
echo cp ${rbase}_mem001.a $ddir/${rbase}.a
echo cp ${rbase}_mem001.b $ddir/${rbase}.b
echo cp ${rbase}ICE.uf    $ddir
cp ${rbase}_mem001.a $ddir/${rbase}.a
cp ${rbase}_mem001.b $ddir/${rbase}.b
cp ${rbase}ICE.uf    $ddir

# cp to migrate
if [ 1 -eq 1 ]
then
   echo cp $gfil $bdir/$yr
   cp $gfil $bdir/$yr
fi
