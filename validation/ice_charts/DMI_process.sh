source $SWARP_ROUTINES/source_files/hex_vars.src
if [ $# -eq 1 ]
then
   cyear=$1
   Indir=/work/shared/nersc/msc/DMI_icecharts/sigrid/$cyear
else
   Indir=`pwd`
fi

outdir=/work/shared/nersc/msc/DMI_icecharts/MIZ_FA

ME=`readlink -f $0`
Vdir=`dirname $ME`

# load python
[ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
module load python/2.7.9-dso
module load geos
export PYTHONPATH=$PYTHONPATH:$SWARP_ROUTINES/py_funs

# extract MIZ from ice chart
pyscript=$Vdir/read_icechart.py
for indir in $Indir/*
do
   echo $python $pyscript --indir=$indir --outdir=$outdir --chart_source=DMI --MIZ_criteria="FA_only" --overwrite=False
   $python $pyscript --indir=$indir --outdir=$outdir --chart_source=DMI --MIZ_criteria="FA_only" --overwrite=False
done

# collect time series for the year:
$python $Vdir/collect_MIZ_widths.py $outdir/MIZ_polys_classified/$cyear
