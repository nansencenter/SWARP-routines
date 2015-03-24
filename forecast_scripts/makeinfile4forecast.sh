#!/usr/bin/ksh
# ======================================================================
# Six arguments:
# 1. run_version
# 2. reference_year
# 3. from_day (-1 at 00:00UTC)
# 4. to_day (+0)
# 5. to_day (+1)
# 6. to_day (+2 at 00:00 and 12:00UTC)
# Other infile parameters must be changed in infile.mal
# ======================================================================
#MAINDIR=/home/nersc/bergh/Realtime
MAINDIR=$HOME/SWARP-routines/forecast_scripts

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
    if [ -f infile.in ]
    then
      rm -f infile.in
    fi

    cat $file | sed \
    -e "s/nd1/$thr/g" \
    -e "s/nd2/$fou/g" \
    -e "s/nd3/$fiv/g" \
    -e "s/ref_year/$2/g" \
    -e "s/run_gen/$1/g" \
    > infile.in

    echo " * infile.in created"
  else
    echo " ERROR : infile.mal does not exist" 
    exit
  fi
fi

