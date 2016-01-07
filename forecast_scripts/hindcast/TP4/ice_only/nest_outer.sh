#!/bin/bash
#
# Routine sets up pivot points needed by the outer model when interpolating 
# model values to an inner nested model.
# Programs used are "nestpivots" - this routine was a part of nersc hycom before-
# now moved outside for easier debugging / error checking
#
# Final files are placed in $BASEDIR/nest_nersc/$E/
#
# REQUIRES: a working topography for experiment $E (regional.grid/regional.depth)


#set -x
pget=cp
pput=cp
# Experiment also needs topography (through experiment number)
if [ $# -ne 3 ] ; then
   echo "This routine will set up nesting from the outer model "
   echo "You need to specify the experiment number of the outer "
   echo "mode (this model), and the path to the base directory "
   echo "of the inner model, along with the experiment number of"
   echo "the inner model"
   echo
   echo "Example:"
   echo "  $(basename $0) 01.0 /work/knutali/h22_tst/FRMa0.05/ 01.1"
   echo "where 01.0 is the experiment number of the outer model and"
   echo "01.1 is the experiment number of the inner model. The path"
   echo "to the top-level directory of nested model is in this example "
   echo "/work/knutali/h22_tst/FRMa0.05/"
   exit
fi
export X=$1
PATHINNER=$2
IX=$3

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly

#########################################################################################################
# TW modified
# NB need to launch from BASEDIR eg TP4a0.12
# cd $(dirname $0)/../../ 
export BASEDIR=$(pwd)/  
#########################################################################################################

cd -
echo "Sourcing Outer model"
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }
echo X=$X E=$E T=$T R=$R

# Check that pointer to MSCPROGS is set
if [ -z ${MSCPROGS} ] ; then
   echo "MSCPROGS Environment not set "
   exit
else
   if [ ! -d ${MSCPROGS} ] ; then
      echo "MSCPROGS not properly set up"
      echo "MSCPROGS not a directory at ${MSCPROGS}"
      exit
   fi
fi

# Rename these variables, cause we need to get the same from the outer model
OX=$X
OE=$E
OT=$T
OR=$R


# Now source inner model setup
echo
X=$IX
echo "Sourcing Inner model"
source ${PATHINNER}/REGION.src || { echo "Could not source ${PATHINNER}/REGION.src" ; exit 1 ; }
source ${PATHINNER}/expt_$X/EXPT.src || { echo "Could not source ${PATHINNER}/expt_$X/EXPT.src" ; exit 1 ; }
echo X=$X E=$E T=$T R=$R

# Rename these variables, cause we need to get the same from the outer model
IX=$X
IE=$E
IT=$T
IR=$R

#Unset these variables to avoid unwanted effects
unset  -v X E T R

#
#
# --- Create nesting fields for outer model (should be compatible with standard hycom)
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
D=$BASEDIR/nest_nersc/$OE/outer/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

# Get regional files and grid.info from out model
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b grid.info
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b grid.info
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b    || \
   { echo "Could not get outer regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a    || \
   { echo "Could not get outer regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${OR}_${OT}.a regional.depth.a || \
   { echo "Could not get outer regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${OR}_${OT}.b regional.depth.b || \
   { echo "Could not get outer regional.depth.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/grid.info grid.info || \
   { echo "Could not get outer grid.info file " ; exit 1 ; }

# Make nest directory in SCRATCH - Then get regional files from inner model
mkdir -p nest 
cd nest
${pget} ${PATHINNER}/topo/regional.grid.b regional.grid.b    || \
   { echo "Could not get inner regional.grid.b file " ; exit 1 ; }
${pget} ${PATHINNER}/topo/regional.grid.a regional.grid.a    || \
   { echo "Could not get inner regional.grid.a file " ; exit 1 ; }
${pget} ${PATHINNER}/topo/depth_${IR}_${IT}.a regional.depth.a || \
   { echo "Could not get inner regional.depth.a file " ; exit 1 ; }
${pget} ${PATHINNER}/topo/depth_${IR}_${IT}.b regional.depth.b || \
   { echo "Could not get inner regional.depth.b file " ; exit 1 ; }
   # add the field necessary in preprocess
${pget} ${PATHINNER}/topo/regional.grid.b ${D}/regional.grid_${IR}_${IT}.b   || \
   { echo "Could not get inner regional.grid.b file " ; exit 1 ; }
${pget} ${PATHINNER}/topo/regional.grid.a ${D}/regional.grid_${IR}_${IT}.a   || \
   { echo "Could not get inner regional.grid.a file " ; exit 1 ; }
${pget} ${PATHINNER}/topo/depth_${IR}_${IT}.a ${D}/depth_${IR}_${IT}.a || \
   { echo "Could not get inner depth_${IR}_${IT}.a file " ; exit 1 ; }
${pget} ${PATHINNER}/topo/depth_${IR}_${IT}.b ${D}/depth_${IR}_${IT}.b || \
   { echo "Could not get inner depth_${IR}_${IT}.b file " ; exit 1 ; }


####################################################################################################
#Now we dump the path of nesting folder into nesting.in
# TW modified
if ! [ -f ${BASEDIR}/expt_${OX}/nesting.in ] 
then
   touch ${BASEDIR}/expt_${OX}/nesting.in
   echo "${PATHINNER}/expt_${IX}/nest_out_${IR}_${IT}/"> ${BASEDIR}/expt_${OX}/nesting.in
   echo "the nesting will be dump in ${PATHINNER}/expt_${IX}/nest_out_${IR}_${IT}/"
else
   echo "Already have a file ${BASEDIR}/expt_${OX}/nesting.in"
   echo "The model is already dumping nesting for the following model:"
   cat ${BASEDIR}/expt_${OX}/nesting.in
   echo "${PATHINNER}/expt_${IX}/nest_out_${IR}_${IT}/">> ${BASEDIR}/expt_${OX}/nesting.in
fi

echo "you are now dumping nesting for the following nested models:"
cat ${BASEDIR}/expt_${OX}/nesting.in
echo "Check the permissions !!!!"
echo "if you are not happy with the destination edit ${BASEDIR}/expt_${OX}/nesting.in"
cd ../
####################################################################################################

#Create nesting relaxation mask - nestpivots must be set up properly with no gotchas
#Inner model grid files are now in local nest directory
${MSCPROGS}/src/Nesting-2.2/nestpivots ./nest/

# Move to experiment directory - prepend with id and topo version of inner model
cd nest
for i in pivots* ; do
   mv $i $D/${IR}_${IT}_$i
done

echo
echo "All looks ok - Nesting file and diagnostic file in $D"
echo "Nesting file names are prepended with region=$IR and topo version=$IT (from inner model)"
