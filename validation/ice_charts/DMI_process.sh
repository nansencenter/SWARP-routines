source $SWARP_ROUTINES/source_files/hex_vars.src
if [ $# -eq 1 ]
then
   cyear=$1
else
   echo "Usage: `basename $0` [year]"
   exit
fi

indir=`pwd`
outdir=/work/shared/nersc/msc/DMI_icecharts/sigrid/$cyear

ME=`readlink -f $0`
Vdir=`dirname $ME`

# load python
[ -f /etc/bash.bashrc ] && . /etc/bash.bashrc
module load python/2.7.9-dso
export PYTHONPATH=$PYTHONPATH:$SWARP_ROUTINES/py_funs

# extract MIZ from ice chart
pyscript=$Vdir/read_icechart.py
$python --indir=$indir --outdir=$outdir --chart_source=DMI --MIZ_criteria="FA_only" --overwrite=False

# TODO MIZ width/area?
