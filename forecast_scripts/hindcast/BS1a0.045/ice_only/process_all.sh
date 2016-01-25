source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -eq 0 ]
then
   echo "`basename $0` [batch number]"
   exit
else
   batch=$1
fi
R=BS1a0.045
T=04
HCtype=ice_only
rungen=${R:0:3}
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/$R/$HCtype

ddir=$BS1_REALTIME/expt_01.2/data
ndir=$BS1_REALTIME/expt_01.2/nest_out_${R}_${T}
rdir0="$BS1_REALTIME/../results_hindcasts/$R/"
rdir="$rdir0/$HCtype"
mkdir -p $rdir0
mkdir -p $rdir


if [ $batch -eq 1 ]
then
   year1=2015
   year2=2015
   day1=67
   day2=151
elif [ $batch -eq 2 ]
then
   year1=2015
   year2=2016
   day1=242
   day2=17
fi


Day1=`printf %3.3d $day1`
Day2=`printf %3.3d $day2`
Rdir0=$rdir/batch_${year1}_${Day1}_${year2}_${Day2}
Rdir=$Rdir0/binaries
mkdir -p $Rdir0
mkdir -p $Rdir

echo " "
echo "Base dir: $Rdir0"
echo " "

bdir1=$Rdir/DAILY
bdir2=$Rdir/restarts
bdir3=$Rdir/archv
bdir4=$Rdir/nesting
mkdir -p $bdir1
mkdir -p $bdir2
mkdir -p $bdir3
mkdir -p $bdir4

figdir=$Rdir0/figures
mkdir -p $figdir
figdir=$Rdir0/figures/DAILY
mkdir -p $figdir
figdir=$Rdir0/figures/DAILY/pngs
mkdir -p $figdir
gifdir=$Rdir0/figures/DAILY/gifs
mkdir -p $gifdir

idir=$Rdir0/info
mkdir -p $idir
cp $hd/inputs/flags                 $idir
cp $hd/inputs/blkdat.input          $idir
cp $hd/inputs/infile.in.batch$batch $idir/infile.in
# cp $ddir/../log/mpijob.out          $idir
# exit


DAILY=1 # copy DAILY files
RESTS=1 # copy restarts
ARCHV=1 # copy archive files
NESTS=1 # copy nesting
PLOTS=1 # plots for rough qual control

if [ $DAILY -eq 1 ]
then
   # cp DAILY avg's
   for afil in $ddir/${rungen}DAILY*.a
   do
      Bdir=$bdir1

      f=`basename $afil`
      bfil=${afil%.a}.b
      #
      YEAR1=${f:9:4}
      DAY1=10#${f:14:3}
      YEAR=${f:18:4}
      DAY=10#${f:23:3}
      day=$((DAY+0))
      # echo $YEAR $DAY $day

      use_file=0
      if [ $year2 -eq $year1 ]
      then
         if [ $day -lt $day2 ] && [ $day -ge $day1 ]
         then
            use_file=1
         fi
      else
         if [ $YEAR -eq $year1 ] && [ $day -ge $day1 ]
         then
            use_file=1
         elif [ $YEAR -eq $year2 ] && [ $day -lt $day2 ]
         then
            use_file=1
         elif [ $YEAR -gt $year2 ] && [ $YEAR -lt $year2 ]
         then
            use_file=1
         fi
      fi
      if [ $use_file -eq 1 ]
      then
         echo " "
         echo "cp $afil $Bdir"
         echo "cp $bfil $Bdir"
         echo " "
         cp $afil $Bdir
         cp $bfil $Bdir
      fi

   done
fi

if [ $RESTS -eq 1 ]
then
   # cp restarts
   for afil in $ddir/${rungen}restart*.a
   do
      Bdir=$bdir2

      f=`basename $afil`
      bfil=${afil%.a}.b
      ufil=${afil%.a}ICE.uf
      #
      YEAR=${f:10:4}
      DAY=10#${f:15:3}
      day=$((DAY+0))
      # echo $YEAR $DAY $day

      use_file=0
      if [ $year2 -eq $year1 ]
      then
         if [ $day -ge $day1 ] && [ $day -le $day2 ]
         then
            use_file=1
         fi
      else
         if [ $YEAR -eq $year1 ] && [ $day -ge $day1 ]
         then
            use_file=1
         elif [ $YEAR -eq $year2 ] && [ $day -le $day2 ]
         then
            use_file=1
         elif [ $YEAR -gt $year2 ] && [ $YEAR -le $year2 ]
         then
            use_file=1
         fi
      fi

      if [ $use_file -eq 1 ]
      then
         echo " "
         echo "cp $afil $Bdir"
         echo "cp $bfil $Bdir"
         echo "cp $ufil $Bdir"
         echo " "
         cp $afil $Bdir
         cp $bfil $Bdir
         cp $ufil $Bdir
      fi

   done
fi

if [ $ARCHV -eq 1 ]
then
   # cp archv files
   for afil in $ddir/${rungen}archv.*.a
   do
      Bdir=$bdir3

      f=`basename $afil`
      bfil=${afil%.a}.b
      #
      YEAR=${f:9:4}
      DAY=10#${f:14:3}
      day=$((DAY+0))
      # echo $YEAR $DAY $day

      use_file=0
      if [ $year2 -eq $year1 ]
      then
         if [ $day -ge $day1 ] && [ $day -le $day2 ]
         then
            use_file=1
         fi
      else
         if [ $YEAR -eq $year1 ] && [ $day -ge $day1 ]
         then
            use_file=1
         elif [ $YEAR -eq $year2 ] && [ $day -le $day2 ]
         then
            use_file=1
         elif [ $YEAR -gt $year2 ] && [ $YEAR -le $year2 ]
         then
            use_file=1
         fi
      fi

      if [ $use_file -eq 1 ]
      then
         echo " "
         echo "cp $afil $Bdir"
         echo "cp $bfil $Bdir"
         echo " "
         cp $afil $Bdir
         cp $bfil $Bdir
         # exit
      fi

   done
fi

if [ $NESTS -eq 1 ]
then
   # cp nesting files
   for afil in $ndir/nest*.hdr
   do
      Bdir=$bdir4

      f=`basename $afil`
      bfil=${afil%.hdr}
      # echo $afil $bfil $f
      
      YEAR=${f:5:4}
      DAY=10#${f:10:3}
      HOUR=10#${f:14:2}
      day=$((DAY+0))
      hour=$((HOUR+0))
      # echo $YEAR $DAY $day

      use_file=0
      if [ $year2 -eq $year1 ]
      then
         if [ $day -ge $day1 ] && [ $day -le $day2 ]
         then
            use_file=1
         fi
      else
         if [ $YEAR -eq $year1 ] && [ $day -ge $day1 ]
         then
            use_file=1
         elif [ $YEAR -eq $year2 ] && [ $day -le $day2 ]
         then
            use_file=1
         elif [ $YEAR -gt $year2 ] && [ $YEAR -le $year2 ]
         then
            use_file=1
         fi
      fi

      if [ $YEAR -eq $year2 ] && [ $day -eq $day2 ] && [ $hour -gt 0 ]
      then
         use_file=0
      fi

      if [ $use_file -eq 1 ]
      then
         echo " "
         echo "cp ${bfil}* $Bdir"
         echo " "
         cp $bfil* $Bdir
         # exit
      fi

   done
fi

if [ $PLOTS -eq 1 ]
then
   # test plots from DAILY avg
   $python $hd/make_pngs.py --indir=$bdir1 --outdir=$figdir
   cd $figdir
   for direc in *
   do
      echo convert -delay 15 -loop 0 $direc/+.png $direc.gif
      convert -delay 15 -loop 0 $direc/*.png $direc.gif
      mv $direc.gif $gifdir
   done
fi
