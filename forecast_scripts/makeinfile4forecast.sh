#!/usr/bin/ksh
# ======================================================================
# Six arguments:
# 1. run_version
# 2. reference_year
# 3. from_day (-1 at 00:00UTC)
# 4. to_day (+0)
# 5. to_day (+2 at 00:00 and 12:00UTC)
# Other infile parameters must be changed in infile.mal.outer (TP4) or infile.mal
# ======================================================================

#MAINDIR=/home/nersc/bergh/Realtime
#MAINDIR=$HOME/giacomo/SWARP-routines/forecast_scripts

# set SWARP_ROUTINES with this line in ~/.bash_profile
# (NB set to correct path where routines are)
# export SWARP_ROUTINES=$HOME/GITHUB-REPOSITORIES/SWARP-routines
MAINDIR=$SWARP_ROUTINES/forecast_scripts

# set place of realtime model in .bash_profile also
tp4dir=$TP4_REALTIME
xdir=$tp4dir/expt_01.1
#xdir=$MAINDIR/infiles

if [ $# -ne 5 ]
then
  echo "Usage: $0 <5 inputs> (see script)"
  exit
else

  echo "-------------------------------------------------------------------------"
  echo " makeinfile4forecast.sh $1"
  echo "-------------------------------------------------------------------------"

  rungen=$1
  thr=`printf '%3.3d' $3`	#three digits modification
  fou=`printf '%3.3d' $4`	#three digits modification
  fiv=`printf '%3.3d' $5`	#three digits modification
  cd ${MAINDIR}/infiles
  if [ -z "${rungen}" -o "${rungen}" == "TP4" ]
  then
     file=infile.mal.outer
  else
     file=infile.mal
  fi
  
  if [ -f $file ]
  then
    echo "Run version:    $1"
    echo "Reference year: $2"
    echo "From day:       $thr"
    echo "To day 1:       $fou"
    echo "To day 2:       $fiv"
    if [ -f $xdir/infile.in ]
    then
      rm -f $xdir/infile.in
    fi

    cat $file | sed \
    -e "s/nd1/$thr/g" \
    -e "s/nd2/$fou/g" \
    -e "s/nd3/$fiv/g" \
    -e "s/ref_year/$2/g" \
    -e "s/run_gen/$1/g" \
    > $xdir/infile.in

    echo " $xdir/infile.in created"
  else
    echo " ERROR : infile.mal does not exist" 
    exit
  fi
fi

