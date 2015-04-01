#!/bin/bash
# Called by framstrait_forecast_postprocessing.sh
# Three inputs : region, julday, year_start
# ======================================================================
# The script seeks for file named iceforecast.$region and uses values therein
# to determine some aspects of the output data.
# 1. Create images with ice_forecast_maps.m
# 2. Copy and archive images to the webpage server
# 3. Call gen_wp.sh to update the webpage
# ======================================================================

# Set up directory variables. WORKDIR is temporary work directory. 
# CRONDIR is the location of this script and its input files (see below)
MAINDIR=$HOME/Realtime
WORKDIR=$MAINDIR/workdir
jobdir=workdir
# WEBDIR=$MAINDIR/webpage
EXPDIR=/work/bergh/FR1a0.03/expt_01.1/data/
#AMSRDIR=/Network/Servers/sverdrup-e3.nersc.no/Data/ssmi/NewSystem
# do not work since 24/04/2011 
# AMSRDIR=/Data/ssmi/NewSystem
#AMSRDIR1=/Data/ArcticROOS/norsex/NORSEX/AMSR/
#AMSRDIR2=/Data/ArcticROOS/web_arctic/norsex/DATA/

# This sets global variables used in function call "alert_message"
# the maillist variable is the person(s) to receive alerts per mail (when
# CRON=1)
CRON=1 # Set to 1 when run in cron, anything else if run interactively
MAILLIST=jon.bergh@nersc.no

# Variables used when calling "alert_message" - indicates success or failure
#PROCID=$(basename $0) # PROCID is the name of this routine withouth path
PROCID=$0 # PROCID is the name of this routine withouth path
SUCCESS="SUCCESS : $PROCID"
ERROR="ERROR : $PROCID"

region=$1
julday=$2
year_start=$3

# Quit if no arguments are supplied to this routine
#[ $# -ne 3 ] && alert_message "$ERROR"  "Routine needs 3 arguments" 
if [ $# -ne 3 ] ; then
  echo "Usage: $0 <3 inputs> (see script)"
  exit
fi


# This function is used to send messages per mail when run CRON=1, otherwise
# it just gives a message to terminal
alert_message() {
   SUBJECT=$1
   MESG=$2
   if [ $CRON -eq 1 ] ; then # Submitted through cron 
      echo "$MESG"
      echo "$MESG" | mail  -s "$SUBJECT" $MAILLIST
   else
      echo "$MESG"
   fi
   exit 0
}

echo "-------------------------------------------------------------------------"
echo " framstrait_ice_forecast.sh"
echo "-------------------------------------------------------------------------"
echo ""

hycom_day=`$MAINDIR/subprogs/julday2dayinyear_out $julday`
day=`echo 00$hycom_day | tail -4c`
year=`$MAINDIR/subprogs/julday2year_out $julday`

hycom_day1=`expr $day + 1`  # HYCOM days > 365 in case of new year shift
hycom_day2=`expr $day + 2`
hycom_day3=`expr $day + 3`

year1=$year_start
year2=$year_start

if [ $hycom_day1 -gt 364 ] ; then
  day1=`expr $hycom_day1 - 365`
  year1=`expr $year_start + 1`
else
  day1=$hycom_day1
fi
day1=`echo 00$day1 | tail -4c`

if [ $hycom_day2 -gt 364 ] ; then
  day2=`expr $hycom_day2 - 365`
  year2=`expr $year_start + 1`
else
  day2=$hycom_day2
fi
day2=`echo 00$day2 | tail -4c`

# We will look for input files called iceforecast.$region.
# Exit if not found.
#[ ! -f $MAINDIR/iceforecast.${region} ] && alert_message "$ERROR" "iceforecast.$region does not exist"
if [ ! -f $MAINDIR/iceforecast.${region} ] ; then
   echo " ERROR : iceforecast.$region does not exist"
   exit
else
   echo " * The file iceforecast.$region exists"
   echo ""
fi

# Get variables from input file iceforecast.$region
# proj_eval - The mapping to use for the matlab m_map function.
# TARGETDIR - This denotes the base directory of the web page to create.
# ARCHIVE   - An archive directory for old ice forecast images.
# NTIMES    - How many dates to produce forecasts for. We process the NTIMES
#             latest dates in the OpenDAP dataset.
proj_eval=$(head -n 2 $MAINDIR/iceforecast.${region} | tail -n1 )
TARGETDIR=$(head -n 3 $MAINDIR/iceforecast.${region} | tail -n1 )
ARCHIVE=$(head -n 4 $MAINDIR/iceforecast.${region} | tail -n1 )
datfile=$TARGETDIR/forecast.dat

echo " * Reading forecast.$region : "
echo ""
echo " PROJECTION = $proj_eval"
echo " TARGET     = $TARGETDIR"
echo " ARCHIVE    = $ARCHIVE"
echo ""

# Make temporary work directory if not already present.
# Brutally clean up that directory before proceeding.
[ ! -d $WORKDIR ] && mkdir -p $WORKDIR

echo " * Cleaning the workdir"
echo ""
cd $WORKDIR && rm -f $WORKDIR/*

# Put all necessary files in WORKDIR
ln -sf $EXPDIR/FR1DAILY_${year_start}_${day}_${year_start}_${day}.*  $WORKDIR
ln -sf $EXPDIR/FR1DAILY_${year_start}_${day}_${year1}_${day1}.*  $WORKDIR
ln -sf $EXPDIR/FR1DAILY_${year_start}_${day}_${year2}_${day2}.*  $WORKDIR

#ln -sf $EXPDIR/latlon.dat        $WORKDIR
#ln -sf $EXPDIR/ndepths400x320.uf $WORKDIR
#ln -sf $EXPDIR/newpos.uf         $WORKDIR
ln -sf $EXPDIR/regional.grid.*   $WORKDIR
ln -sf $EXPDIR/regional.depth.*  $WORKDIR
ln -sf $EXPDIR/grid.info         $WORKDIR

echo " * Launching Matlab and ice_forecast_maps.m"
echo ""

# This seems to be the safest way to run Matlab when no display is set.
# It starts a "fake" X window session in display number 9
#XVFB=`which Xvfb`
#$XVFB :9 -screen 0 1024x768x24 &

# The bulk of this routine is done by matlab script ice_forecast_maps.
# It is run inside the virtual frame buffer.
#env DISPLAY=:9.0 matlab -nodesktop -nosplash <<EOF
#warning off MATLAB:FINITE:obsoleteFunction;
#$proj_eval;
#ice_forecast_maps;
#exit;
#EOF
cat $MAINDIR/pbs_postprocess.sh | sed \
   -e "s/region/framstrait/g" \
   -e "s/type/forecast/g" \
   -e "s/jobdir/$jobdir/g" \
   > matlab_framstrait.sh

/opt/torque/torque-2.3.4-snap.200809221601/bin/qsub  matlab_framstrait.sh 




while  ! [ -f $WORKDIR/forecast.dat ]
do
  echo " * Waiting 10 minutes ..."
  echo ""
  sleep 300
done




# Replace FR1DAILY* files with those needed for validation
rm FR1DAILY_${year_start}_${day}*

valday1=`expr $hycom_day - 6`
valday2=`expr $hycom_day - 5`
valday3=`expr $hycom_day - 4`

valyear1=$year_start
valyear2=$year_start
valyear3=$year_start

if [ $valday1 -lt 0 ] ; then
  valday1=`expr 365 + $valday1`
  valyear1=`expr $year_start - 1`
fi
valday1=`echo 00$valday1 | tail -4c`
ln -sf $EXPDIR/FR1DAILY_${valyear1}_${valday1}_${valyear1}_${valday1}.*  $WORKDIR

if [ $valday2 -lt 0 ] ; then
  valday2=`expr 365 + $valday2`
  valyear2=`expr $year_start - 1`
fi
valday2=`echo 00$valday2 | tail -4c`
ln -sf $EXPDIR/FR1DAILY_${valyear2}_${valday2}_${valyear2}_${valday2}.* $WORKDIR

if [ $valday3 -lt 0 ] ; then
  valday3=`expr 365 + $valday3`
  valyear3=`expr $year_start - 1`
fi
valday3=`echo 00$valday3 | tail -4c`
ln -sf $EXPDIR/FR1DAILY_${valyear3}_${valday3}_${valyear3}_${valday3}.* $WORKDIR

valjulday=`expr $julday - 6`
valjulday1=`expr $julday - 5`
valjulday2=`expr $julday - 4`

valdate=`echo $valjulday | $MAINDIR/subprogs/julday2date_short`
valdate1=`echo $valjulday1 | $MAINDIR/subprogs/julday2date_short`
valdate2=`echo $valjulday2 | $MAINDIR/subprogs/julday2date_short`

#echo ""
#echo " * Retrieving AMSR-E data"
#echo ""
#AMSDIR not updated since 24-04-2011
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR/NORSEX/AMSR/AMSRNH${valdate}.uf .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR/NORSEX/AMSR/AMSRNH${valdate1}.uf .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR/NORSEX/AMSR/AMSRNH${valdate2}.uf .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR/DATA/psn12lons_v2.dat .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR/DATA/psn12lats_v2.dat .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR1/AMSRNH${valdate}.uf .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR1/AMSRNH${valdate1}.uf .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR1/AMSRNH${valdate2}.uf .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR2/psn12lons_v2.dat .
#scp -i ~/.ssh/nansenkey topaz@nansen.nersc.no:$AMSRDIR2/psn12lats_v2.dat .

echo ""
echo " * Retrieving OSI SAF data from myocean.met.no"
echo ""
## the file format: ice_conc_nh_201202051200.nc
## ssmifile=`ice_conc_nh_${valdate}1200.nc`
# read in to ncftp.in file and execute
# Using the ncftp bookmark, need to run the full version first as
# ncftp -u username -p password ftp.myocean.met.no
# When asked, save bookmark as myocean.
# Bookmarks will be saved in $HOME/.ncftp/bookmarks
# 
#get ice_conc_nh_${valdate}1200.nc
#get ice_conc_nh_${valdate1}1200.nc
#get ice_conc_nh_${valdate2}1200.nc

LIST=""
for i in $valdate $valdate1 $valdate2
do
nr=$i"1200"
if [ ! -e ice_conc_nh_$nr.nc ] ; then
  LIST=$LIST" ice_conc_nh_$nr.nc"
fi
done

  cat > ncftp.in<<EOF
open myocean 
cd SIW-TAC/SIW-OSISAF-GLO-SIT_SIE_SIC-OBS/conc/
get $LIST
set confirm-close no
bye
EOF
   ncftp < ncftp.in
   rm ncftp.in


echo " * Launching Matlab and ice_validation_maps.m"
echo ""
cat $MAINDIR/pbs_postprocess.sh | sed \
   -e "s/region/framstrait/g" \
   -e "s/type/validation/g" \
   -e "s/jobdir/$jobdir/g" \
   > matlab_framstrait.sh
/opt/torque/torque-2.3.4-snap.200809221601/bin/qsub matlab_framstrait.sh 

occ=0
while  ! [ -f $WORKDIR/validation.dat ]
do
  echo " * Waiting 4 minutes ..."
  echo ""
  let occ=occ+1
  sleep 240
done

# The matlab script above will create several jpeg files. These files have names
# of this type: "ice-forecast-05-Dec-2008_variable.jpg"
# where variable is either drift, hice or speed.
# The file forecast.dat is also created, which consists of one line for each
# time processed, with each line containing a date and the jpeg files created for
# this date. The script proceeds by using these files to set up the ice forecast
# webpage.

# Move the existing forecast jepg images into the ARCHIVE directory.
echo ""
echo " * Archiving images"
echo ""
echo " * Moving previous images to $ARCHIVE using SSH"
echo ""
command1=`ssh -i ~/.ssh/nansenkey topaz@nansen.nersc.no "mv $TARGETDIR/ice-*.jpg $ARCHIVE"`
#command1=`ssh topaz@nansen.nersc.no "mv $TARGETDIR/ice-*.jpg $ARCHIVE"`
echo "$command1"

# Now copy the newly created ice forecast images to the TARGETDIR directory.
echo " * Copying newly created images to $TARGETDIR using SSH"
echo ""
scp -p -i ~/.ssh/nansenkey $WORKDIR/ice-forecast*.jpg topaz@nansen.nersc.no:$TARGETDIR
#scp -p $WORKDIR/ice-forecast*.jpg topaz@nansen.nersc.no:$TARGETDIR 
scp -p -i ~/.ssh/nansenkey $WORKDIR/forecast.dat topaz@nansen.nersc.no:$TARGETDIR
#scp -p $WORKDIR/forecast.dat topaz@nansen.nersc.no:$TARGETDIR
scp -p -i ~/.ssh/nansenkey $WORKDIR/ice-validation*.jpg topaz@nansen.nersc.no:$TARGETDIR
#scp -p -i ~/.ssh/id_dsa_pub $WORKDIR/ice-validation*.jpg topaz@nansen.nersc.no:$TARGETDIR
scp -p -i ~/.ssh/nansenkey $WORKDIR/validation.dat topaz@nansen.nersc.no:$TARGETDIR
#scp -p -i ~/.ssh/id_dsa_pub $WORKDIR/validation.dat topaz@nansen.nersc.no:$TARGETDIR
command2=`ssh -i ~/.ssh/nansenkey topaz@nansen.nersc.no "chmod 644 $TARGETDIR/*.jpg"`
#command2=`ssh topaz@nansen.nersc.no "chmod 644 $TARGETDIR/*.jpg"`
echo "$command2"

# This runs the gen_wp.sh script which sets up the forecast web page for the
# region. Note the location of the gen_wp.sh. Input to the gen_wp.sh script is 
# the name of the region and the forecast.dat file. The Template directory 
# contains some template files which are needed to set up the final web page.
# See that script for more info. 
# Note that at this point we are in the target directory, whereas the template
# files are in the same location as gen_wp.sh
echo " * Launching gen_wp.sh for $region"
echo ""
command3=`ssh -i ~/.ssh/nansenkey topaz@nansen.nersc.no "$TARGETDIR/gen_wp.sh $region $datfile"`
#command3=`ssh topaz@nansen.nersc.no "$TARGETDIR/gen_wp.sh $region $datfile"`
echo "$command3"

# Everything should be ok here - notify people of our success
#MAILSENT=1
#alert_message "$SUCCESS" "Web page updated for Fram Strait : $(ls *.jpg)"

echo " * Web page updated for Fram Strait"
echo ""
