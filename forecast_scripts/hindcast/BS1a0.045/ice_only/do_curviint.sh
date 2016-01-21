# #! /bin/bash

# run in base directory of local model
rm -rf curviint
BASE=/work/timill/RealTime_Models
if [ $# -eq 1 ]
then
   batch=$1
else
   echo "Usage: do_curviint.sh [batch_number]"
   echo " batch_number=1 starts 2015_060"
   echo " batch_number=2 starts 2015_060"
   exit
fi

############################################################
#
tp4_path="$BASE/TP4a0.12"
X_out="01.4"                                          #expt no for TP4
if [ $batch -eq 1 ]
then
   rfil0="TP4restart2015_067_00"                         #name of TP4 restart file
else
   rfil0="TP4restart2015_242_00"                         #name of TP4 restart file
fi
rfil="${tp4_path}/expt_${X_out}/data/${rfil0}.a"      #full path to TP4 restart file

if [ 1 -eq 1 ]
then
   #BS1a0.045 model
   rungen="BS1"
   X_inn="01.2"                                       #expt no for BS1
   E_inn="012"                                       #expt no for FR1
   path_inn="$BASE/BS1a0.045"
else
   #FR1a0.03 model
   rungen="FR1"
   X_inn="01.2"                                       #expt no for FR1
   E_inn="012"                                       #expt no for FR1
   path_inn="$BASE/FR1a0.03"
fi
# don't need to change anything below
############################################################

echo " "
echo "infile: $rfil0"
echo " "

cd ${path_inn}
pwd
rm -f curviint/$E_inn/SCRATCH/*

# make the restart file by interpolation from TP4:
echo bin/curviint.sh ${X_inn} ${X_out} ${tp4_path} ${rfil}
other_nersc/bin/curviint.sh ${X_inn} ${X_out} ${tp4_path} ${rfil}


# get them, rename them and copy them to data directory:
E=${X_inn:0:2}${X_inn:3:4}       #eg 01.0 -> 010
d0=$path_inn/curviint/$E         #current location
d1=$path_inn/expt_$X_inn/data    #desired destination

echo ""
echo "**Renaming & copying outputs of curviint.sh..."
echo "...final location of restart files:"
echo ""

cd $d0
for f in $rfil0*
do
   lf=${#f}
   f2=$rungen${f:3:$lf}
   cp $f $d1/$f2
   echo $d1/$f2
done
echo ""
