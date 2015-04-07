#!/bin/ksh
# ======================================================================
# Five arguments:
# 1. run_version
# 2. reference_year
# 3. from_day (-7 at 00:00UTC)
# 4. to_day (+7)
# 5. to_day (+9 at 00:00 and 12:00UTC)
# Other infile parameters must be changed in infile.mal
# ======================================================================
MAINDIR=/home/nersc/bergh/Realtime

if [ $# -ne 5 ]
then
  echo "Usage: $0 <5 inputs> (see script)"
  exit
else

  echo "-------------------------------------------------------------------------"
  echo " makeinfile4topazforecast.sh $1"
  echo "-------------------------------------------------------------------------"

  rungen=$1
  
  cd ${MAINDIR}/infiles
 
    file=infile.mal.topaz
  
  if [ -f $file ]
  then
    echo "Run version:    $1"
    echo "Reference year: $2"
    echo "From day:       $3"
    echo "To day 1:       $4"
    echo "To day 2:       $5"
    if [ -f infile.in ]
    then
      rm -f infile.in
    fi

    cat $file | sed \
    -e "s/nd1/$3/g" \
    -e "s/nd2/$4/g" \
    -e "s/nd3/$5/g" \
    -e "s/reference_year/$2/g" \
    -e "s/run_version/$1/g" \
    > infile.in

    echo " * infile.in created"
  else
    echo " ERROR : infile.mal does not exist" 
    exit
  fi
fi

