source $SWARP_ROUTINES/source_files/hex_vars.src
if [ $# -eq 1 ]
then
   cdate=$1
else
   echo "Usage: `basename $0` [date]"
   echo "- where date is in yyyymmdd format"
   exit
fi

regions="TP4a0.12 BS1a0.045 FR1a0.03"
# regions="TP4a0.12"
# regions="BS1a0.045 FR1a0.03"
FCtypes="ice_only"
basedir="$RTmods/results/"
valdir="$RTmods/validation/"
mkdir -p $valdir

ME=`readlink -f $0`
Vdir=`dirname $ME`

# load python
[ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
module load python/2.7.9-dso
module load geos
export PYTHONPATH=$PYTHONPATH:$SWARP_ROUTINES/py_funs

for reg in $regions
do
   for FCtype in $FCtypes
   do
      mkdir -p $valdir/$reg
      mkdir -p $valdir/$reg/$FCtype

      outdir=$valdir/$reg/$FCtype/$cdate #where results will go (1 dir for each forecast)
      mkdir -p $outdir
      outdir=$outdir/OSISAF #where results will go (1 dir for each forecast)
      mkdir -p $outdir
      figdir=$outdir/figs
      mkdir -p $figdir

      FC_rootdir=$basedir/$reg/$FCtype # where data from forecast is
      echo " "
      echo $reg - $FCtype
      echo "Data from   : $FC_rootdir"
      echo "Results to  : $outdir"
      echo " "
      if [ 1 -eq 1 ]
      then
         echo $python $Vdir/ice_edge_OSISAF_1obs.py --date=$cdate --outdir=$outdir --FC_rootdir=$FC_rootdir
         $python $Vdir/ice_edge_OSISAF_1obs.py --date=$cdate --outdir=$outdir --FC_rootdir=$FC_rootdir
      fi


      # ==============================================
      # make gifs of ice edge plots
      cd $outdir

      # loop over figs from each FC day
      # - order them in chronological order
      ddirs=(FC*days)
      Ndirs=${#ddirs[@]}
      N=$Ndirs

      for dd in FC*days
      do
         f=$dd/*IceEdge_OSISAF*.png
         # echo $f
         if [ -f $f ]
         then
            # echo ln -s $f fig$N.png
            ln -s $f fig$N.png
         fi
         N=$((N-1))
      done

      # ls -lh
      gfil=IceEdge_OSISAF${cdate}.gif
      convert -delay 75 -loop 0 fig*.png $gfil
      mv $gfil $figdir
      rm fig*.png
      # ==============================================


      # ==============================================
      # make gifs of conc anomaly plots
      cd $outdir

      # loop over figs from each FC day
      # - order them in chronological order
      ddirs=(FC*days)
      Ndirs=${#ddirs[@]}
      N=$Ndirs

      for dd in FC*days
      do
         f=$dd/conc_anomaly_OSISAF*.png
         if [ -f $f ]
         then
            ln -s $f fig$N.png
         fi
         N=$((N-1))
      done

      gfil=conc_anomaly_OSISAF${cdate}.gif
      convert -delay 75 -loop 0 fig*.png $gfil
      mv $gfil $figdir
      rm fig*.png
      # ==============================================

   done
done
