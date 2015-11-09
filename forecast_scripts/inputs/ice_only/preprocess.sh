#!/bin/bash
# Error counter for preprocessing

# TW: Executable to get from Build (set in pbsjob.sh):
# Usually "hycom", 
#  but can also be eg "hycom+pat" or "hycom+apa"
#  if using PAT (Performance Analysis Tool);
# Made it an argument of preprocess.sh so only need to
#  change pbsjob.sh if change executable to run
if [ $# -eq 0 ] ; then
   # Default (no argument to script) is hycom
   EXEC="hycom"
else
   EXEC=$1
fi

#Check environment
cd $(dirname 0)/../
export BASEDIR=$(pwd)
cd -
echo BASEDIR=$BASEDIR
if [ -z ${BASEDIR} ] ; then
   tellerror "preprocess.sh: BASEDIR Environment not set "
   exit
fi

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src "; exit 1; }


logfile=$P/log/pre
[ ! -d $P/log ] && mkdir $P/log
[ -f $logfile ] && rm $logfile

#Redirect stderr to log
exec 2>$logfile.err

numerr=0

tellerror () {
  echo "[FATAL  ] $1" 
  let numerr=$numerr+1; 
  echo "[FATAL  ] $1" >> $logfile
}
tellwarn () {
  echo "[WARNING] $1" 
  echo "[WARNING] $1" >> $logfile
}


# Create data and scratch dirs
if [ ! -d $D ] ; then
   mkdir -p $D   || { echo "Could not create data dir : $S " ; exit 1 ; }
fi
if [ ! -d $S ] ; then
   mkdir -p $S   || { echo "Could not create data dir : $S " ; exit 1 ; }
fi


# Set integration limits
# ./subprogs/setlimits.py
# pyt=/work/apps/python/2.7.2-dso-gnu/bin/python
$python subprogs/setlimits.py

# Remove init file after limits are set
[ -f INITIALIZE ]  && rm INITIALIZE

# Source info from current experiment - listed in "THIS_EXPT.src"
cd       $S || { tellerror "no scratch dir $S" ;  exit 1 ;}


#
#
#
# --- Set up time limits for this segment
#
#
mv $P/limits limits             || tellerror "no limits file created by setlimits.py"  
cp $P/blkdat.input blkdat.input || tellerror "No blkdat.input file" 
cp $P/infile.in infile.in       || tellerror "no infile.in" 
cp $P/infile2.in infile2.in     || tellerror "no infile2.in" 


touch ports.input tracer.input infile.evp infile.icestate
rm    ports.input tracer.input infile.evp infile.icestate
[ -s  $P/ports.input ] && cp $P/ports.input ports.input
[ -s  $P/tracer.input ] && cp $P/tracer.input tracer.input
[ -s  $P/infile.evp ] && cp $P/infile.evp infile.evp
[ -s  $P/infile.icestate ] && cp $P/infile.icestate infile.icestate

#
# --- check that iexpt from blkdat.input agrees with E from this script.
# --- fetch some things from blkdat.input 
#
export EB=`grep "'iexpt ' =" blk* | awk '{printf("%03d", $1)}'`
export PRIVER=`grep "'priver' =" blk* | awk '{printf("%1d", $1)}'`
export TRCFLG=`grep "'trcflg' =" blk* | awk '{printf("%1d", $1)}'`
export YRFLAG=`grep "'yrflag' =" blk* | awk '{printf("%1d", $1)}'`
export JERLV=`grep "'jerlv0' =" blk* | awk '{printf("%1d", $1)}'`
export SSSRLX=`grep "'sssflg' =" blk* | awk '{printf("%1d", $1)}'`
export SSTRLX=`grep "'sstflg' =" blk* | awk '{printf("%1d", $1)}'`
export RLX=`grep "'relax ' =" blk* | awk '{printf("%1d", $1)}'`
export TRCRLX=`grep "'trcrlx' =" blk* | awk '{printf("%1d", $1)}'`
export K=`grep "'kdm   ' =" blkdat.input | awk '{printf("%1d", $1)}'`
[ "$EB" == "$E" ] || tellerror " blkdat.input iexpt ${EB} different from this experiment ${E} set in EXPT.src"
echo "Fetched from blkdat.input:"
echo "--------------------------"
echo "EB     is $EB    "
echo "PRIVER is $PRIVER"
echo "TRCFLG is $TRCFLG"
echo "YRFLAG is $YRFLAG"
echo "JERLV  is $JERLV "
echo "SSS    is $SSSRLX"
echo "SST    is $SSTRLX"
echo "RLX    is $RLX"
echo "TRCRLX is $TRCRLX"
echo


#
# --- Fetch some things from infile.in
#
# TODO iday may be ordinal day > days in year
rungen=$(head -n 2 infile.in  | tail -n1 | cut -c1-3)
refyear=$(head -n 3 infile.in  | tail -n1 | cut -c1-4)
refyear=$(printf "%4.4d" $refyear)
iday=$(head -n 4 infile.in  | tail -n1 | cut -c1-4)
iday=`echo 00$iday | tail -4c`
ihour=$(head -n 4 infile.in  | tail -n1 | cut -c7-8)
ihour=`echo 0$ihour | tail -3c`
frc=$(head -n 6 infile.in  | tail -n1 | cut -c1-5 | tr -d " ")
clm=$(head -n 6 infile.in  | tail -n1 | cut -c7-11 | tr -d " ") 
nestoflag=$(head -n 11 infile.in | tail -n1 | awk '{printf("%1s", $1)}')
nestiflag=$(head -n 12 infile.in | tail -n1 | awk '{printf("%1s", $1)}')
tideflag=$(head -n 13 infile.in | tail -n1 | cut -c1 )
tidechoice=$(head -n 13 infile.in | tail -n1 | cut -c3-5 )
gpflag=$(head -n 14 infile.in | tail -n1 | awk '{printf("%1s", $1)}')
echo "Fetched from infile.in:"
echo "-----------------------"
echo "rungen  is              : $rungen"
echo "refyear is              : $refyear"
echo "iday    is              : $iday"
echo "ihour   is              : $ihour"
echo "Forcing option is       : ${frc}"
echo "Climatology option is   : ${clm}"
echo "Outer nest flag is      : $nestoflag"
echo "Inner nest flag is      : $nestiflag"
echo "Tide flag and choice is : $tideflag $tidechoice"
echo "Gridpoint flag is       : $gpflag "
echo


#
# --- Fetch some things from limits
#
tstart=$(cat limits | tr " " -s | sed "s/^[ ]*//" | cut -d " " -f1)
tstop=$(cat limits | tr " " -s | sed "s/^[ ]*//" | cut -d " " -f2)

# Get init flag from start time
init=0
[ `echo $tstart | awk '{print ($1 < 0.0)}'` == 1  ] && init=1
echo "Fetched from limits:"
echo "--------------------"
echo "init is $init"
echo


#################################################################################
# TW: - Have enabled the possibility of dumping snapshots of the ice variables
#        eg conc, thickness, stresses and strains, ice pressure, floe size (if using waves)
#     - Use the flags -DICE_DYN_DIAG and -DICE_DYN_DIAG_DUMP to
#        dump a netcdf file into the directory $IDDdir (in SCRATCH)
IDDdir="$S/ICE_DYN_DIAG"
if [ ! -d $IDDdir ]; then
   mkdir $IDDdir
fi
#################################################################################

#################################################################################
# TW: Waves stuff:

# Check if waves-only run,
# and if so, whether ice conc/thickness come from daily average files or observations;
wav_infil="infile.waves"
latest_WIV="1.2"           #latest version of infile.waves
waves_only=0               #variable used by preprocess_waves.sh and below (don't need restart file, or to use inner/outer nesting if it = 1)
waves_ext_data=0           #variable used by preprocess_waves.sh

if [ -f "${P}/${wav_infil}" ]; then 

   cp $P/$wav_infil .
   WIV=$(head -n 1 $wav_infil  | tail -n1 | awk '{printf("%3s", $1)}')
   wav_flag1=$(head -n 2 $wav_infil | tail -n1 | awk '{printf("%1s", $1)}')
   wav_flag2=$(head -n 3 $wav_infil | tail -n1 | awk '{printf("%1s", $1)}')
   # wav_flag3=$(head -n 6 $wav_infil | tail -n1 | awk '{printf("%1s", $1)}')
   # wav_flag4=$(head -n 7 $wav_infil | tail -n1 | awk '{printf("%1s", $1)}')

   if [ ! "${WIV}" = "${latest_WIV}" ] ; then
      tellerror "preprocess.sh: update version of ${wav_infil} to ${latest_WIV}"
      exit
   fi

   echo "Fetched from $wav_infil:"
   echo "-------------------------------"
   echo "Version of infile.waves : $WIV"
   echo "Waves-only run          : $wav_flag1"
   echo "EXTERNAL c/h DATA?      : $wav_flag2"
   # echo "c DATA OPTION           : $wav_flag3"
   # echo "h DATA OPTION           : $wav_flag4"
   echo

   if [ "$wav_flag1" = "T" ] ; then
      waves_only=1
      if [ "$wav_flag2" = "T" ] ; then
         #Using external data
         #(Alternative is daily average model data)
         waves_ext_data=1

         # #Determine source of concentration data;
         # if [ "$wav_flag3" = "1" ] ; then
         #    #AMSR-E concentration data (6.5km)
         #    c_option=1
         # elif [ "$wav_flag3" = "2" ] ; then
         #    #OSISAF concentration data (25km)
         #    c_option=2
         # elif [ "$wav_flag3" = "3" ] ; then
         #    #AMSR2 concentration data (3.5km)
         #    c_option=3
         # fi

         # #Determine source of thickness data;
         # if [ "$wav_flag4" = "1" ] ; then
         #    #SMOS concentration data (~25km)
         #    h_option=1
         # fi
      fi
   fi
fi

prep_waves=${P}/preprocess_waves.sh
if [ -f "$prep_waves" ]; then 
   # Set paths to wave data and link to them
   echo "Running waves script: $prep_waves"
   echo $prep_waves $waves_only $waves_ext_data
   $prep_waves $waves_only $waves_ext_data
   echo ""
else
   echo "No waves script $prep_waves to run"
fi

#################################################################################
#C
#C --- turn on detailed debugging.
#C
#touch PIPE_DEBUG

#
# --- pget, pput "copy" files between scratch and permanent storage.
# --- Can both be cp if the permanent filesystem is mounted locally.
# --- KAL - we use pget=pput=cp
#Fanf: I/O is too slow on login node we link 
#export pget=/bin/cp
#export pput=/bin/cp
export pget='/bin/ln -sf'
export pput='/bin/ln -sf'

echo
echo "Initialization complete - now copying necessary files to scratch area"

#
#
# --- Set up MPI partition for this segment
#
#
[  "$NMPI" == "" ] && NMPI=0  # Init if unset
if [ $NMPI != 0 ] ; then
  export NPATCH=`echo $NMPI | awk '{printf("%04d", $1)}'`
  /bin/rm -f patch.input
  /bin/cp $BASEDIR/topo/partit/depth_${R}_${T}.${NPATCH}  patch.input || \
     tellerror "Could not get patch.input  ($BASEDIR/topo/partit/depth_${R}_${T}.${NPATCH})"
# /bin/cp $BASEDIR/topo/partit/depth_${R}_${T}.${NPATCH}u patch.input
fi

#
# --- input files from file server.
#
echo "**Retrieving grid and topography files"
${pget} $BASEDIR/topo/regional.grid.a regional.grid.a || telerror "no grid file regional.grid.a" 
${pget} $BASEDIR/topo/regional.grid.b regional.grid.b || telerror "no grid file regional.grid.a" 
${pget} $BASEDIR/topo/depth_${R}_${T}.a regional.depth.a || telerror "no topo file depth_${R}_${T}.a" 
${pget} $BASEDIR/topo/depth_${R}_${T}.b regional.depth.b || telerror "no topo file depth_${R}_${T}.b" 
${pget} $BASEDIR/topo/grid.info  grid.info || telerror "no grid.info file " 


#
# --- Check forcing and climatology option in infile.in
#
if [ "$clm" == "era40" ] ; then
   CLMDIR=$BASEDIR/force/nersc_era40/$E/
elif [ "$clm" == "old" ] ; then
   CLMDIR=$BASEDIR/force/nersc_old/$E/
elif [ "$clm" == "ncepr" ] ; then
   CLMDIR=$BASEDIR/force/nersc_ncepr/$E/
elif [ "$clm" != "prep" ] ; then
   tellerror "Unknown climate option $clm"
   CLMDIR=""
fi

   

#
#
# --- Climatology atmospheric forcing - always copied
#
#
if [ "$CLMDIR" != "" -a "$clm" != "prep" ] ; then
   echo "**Retrieving climatology forcing files"
   for i in tauewd taunwd wndspd radflx shwflx vapmix \
      airtmp precip uwind vwind clouds relhum slp ; do
      # TODO - doesnt handle offsets
      ${pget} $CLMDIR/$i.a      forcing.$i.a || tellerror "Could not fetch $CLMDIR/$i.a"
      ${pget} $CLMDIR/$i.b      forcing.$i.b || tellerror "Could not fetch $CLMDIR/$i.b"
   done

   # TODO - check flags for these
   for i in seatmp surtmp ; do
      ${pget} $CLMDIR/$i.a forcing.$i.a || tellwarn "Could not fetch $CLMDIR/$i.a"
      ${pget} $CLMDIR/$i.b forcing.$i.b || tellwarn "Could not fetch $CLMDIR/$i.b"
   done
fi

#
# --- Synoptic forcing option - for now just sets up links
#
if [ "$frc" != "month" -a  "$clm" != "prep" ] ; then
   echo "**Setting up INLINE_FORCING synoptic forcing"

   # pathvar is name of variable
   pathvar=""
   if [ "${frc}" == "era40" ] ; then
      pathvar="ERA40_PATH"
   elif [ "${frc}" == "era-i" ] ; then
      pathvar="ERAI_PATH"
   elif [ "${frc}" == "ncepr" ] ; then
      pathvar="NCEP_PATH"
   elif [ "${frc}" == "ecmwf" ] ; then
      pathvar="ECMWF_PATH"
      pathvar2="./Ecmwfr"
   elif [ "${frc}" == "metno" ] ; then
      pathvar="METNO_PATH"
      pathvar2="./Met.no"
   elif [ "${frc}" == "ecnc" ] ; then
      pathvar="ECNC_PATH"
      pathvar2="./Ecmwf.nc"
   elif [ "${frc}" == "ecnc2" ] ; then
      pathvar="ECNC_PATH"
      pathvar2="./Ecmwf.nc"
   elif [ "${frc}" == "ecncF" ] ; then
      pathvar="ECNC_PATH"
      pathvar2="./Ecmwf.nc"
   elif [ "${frc}" == "ncepr" ] ; then
      pathvar="NCEP_PATH"
   else
      tellerror "No method for forcing $frc (did you set it up?)"
   fi

   if [ "$pathvar" != "" ] ; then
      pathval="${!pathvar}" # pathval is value of variable named pathvar
      echo "Forcing is $frc, environment variable $pathvar is $pathval"
      if [ "${pathval}" == "" -o ! -d ${pathval:=tomtomtom} ] ; then  #tomtomtom is substituted if pathval is empty
         tellerror "$pathvar not set or not a directory"
      else # Could be present
         if [ "$pathvar2" != ""  ] ; then
            ln -s ${pathval} $pathvar2 || tellwarn "Could not link $pathvar ${pathval} as $pathvar2 (it might be present already)"
         else
            ln -s ${pathval} . || tellwarn "Could not link $pathvar ${pathval} (it might be present)"
         fi
      fi
   else
      tellerror "No pathvar $pathvar set up for forcing $frc. Set it in REGION.src "
   fi
#else
#   tellerror "Forcing variable frc undefined"
fi


#
# --- Pre-prepared forcing option - for now just sets up links
#
if [ "$clm" == "prep" ] ; then
   echo "**Setting up pre-prepared synoptic forcing from force/synoptic/$E"

   DIR=$BASEDIR/force/synoptic/$E/

   for i in tauewd taunwd wndspd radflx shwflx vapmix \
      airtmp precip uwind vwind clouds relhum slp ; do

      echo "|--> $i"

      cp $DIR/${i}.a forcing.${i}.a ||  tellerror "Could not fetch $DIR/$i.a"
      cp $DIR/${i}.b forcing.${i}.b ||  tellerror "Could not fetch $DIR/$i.b"
   done
fi

#
# --- time-invarent heat flux offset
#
#TODO handle offsets
#setenv OF ""
#setenv OF "_413"
if [ "$CLMDIR" != "" ] ; then
   if [  "$OF" != "" ];  then
     touch  forcing.offlux.a  forcing.offlux.b
     if [ ! -s  forcing.offlux.a ] ; then
        ${pget} $CLMDIR/offlux${OF}.a forcing.offlux.a &
     fi
     if [ ! -s  forcing.offlux.b ] ; then
        ${pget} $CLMDIR/offlux${OF}.b forcing.offlux.b &
     fi 
   fi 
fi 


#
# --- river forcing
# --- KAL: rivers are experiment-dependent
#
if [ $PRIVER -eq 1 ] ; then
   echo "**Setting up river forcing"
   ${pget} $BASEDIR/force/rivers/$E/rivers.a forcing.rivers.a || tellerror "Could not get river .a file"
   ${pget} $BASEDIR/force/rivers/$E/rivers.b forcing.rivers.b || tellerror "Could not get river .b file"
fi


#
# --- nutrient river forcing
# --- KAL: rivers are experiment-dependent
#
if [ $TRCFLG -eq 9 ] ; then
   echo "**Setting up river nutrient forcing"
   ${pget} $BASEDIR/force/rivers/$E/rivnitr.a forcing.rivnitr.a || tellerror "Could not get nitrate river .a file"
   ${pget} $BASEDIR/force/rivers/$E/rivnitr.b forcing.rivnitr.b || tellerror "Could not get nitrate river .b file"
   ${pget} $BASEDIR/force/rivers/$E/rivphos.a forcing.rivphos.a || tellerror "Could not get phosphate river .a file"
   ${pget} $BASEDIR/force/rivers/$E/rivphos.b forcing.rivphos.b || tellerror "Could not get phosphate river .b file"
   ${pget} $BASEDIR/force/rivers/$E/rivsili.a forcing.rivsili.a || tellerror "Could not get silicate river .a file"
   ${pget} $BASEDIR/force/rivers/$E/rivsili.b forcing.rivsili.b || tellerror "Could not get silicate river .b file"
   ${pget} $BASEDIR/force/rivers/$E/rivdono.a forcing.rivdono.a || tellerror "Could not get DON river .a file"
   ${pget} $BASEDIR/force/rivers/$E/rivdono.b forcing.rivdono.b || tellerror "Could not get DON river .b file"
fi


#
# --- kpar forcing
#
if [ $JERLV -eq 0 ] ; then
   echo "**Setting up kpar forcing"
   ${pget} $BASEDIR/force/seawifs/kpar.a forcing.kpar.a || tellerror "Could not get kpar.a file"
   ${pget} $BASEDIR/force/seawifs/kpar.b forcing.kpar.b || tellerror "Could not get kpar.b file"
fi


#
# --- relaxation
#
if [ $SSSRLX -eq 1 -o $SSTRLX -eq 1 -o $RLX -eq 1 -o $init -eq 1 ] ; then
   echo "**Setting up relaxation"
   for i in saln temp intf ; do
      j=$(echo $i | head -c3)
      ${pget} $BASEDIR/relax/${E}/relax_$j.a relax.$i.a  || tellerror "Could not get relax.$i.a"
      ${pget} $BASEDIR/relax/${E}/relax_$j.b relax.$i.b  || tellerror "Could not get relax.$i.b"
   done
   if [ $TRCFLG -eq 9 ] ; then
     for i in nitr phos sili oxyg ; do
        j=$(echo $i | head -c3)
        ${pget} $BASEDIR/relax/${E}/relax_$j.a relax.$i.a  || tellerror "Could not get relax.$i.a"
        ${pget} $BASEDIR/relax/${E}/relax_$j.b relax.$i.b  || tellerror "Could not get relax.$i.b"
     done
   fi
fi
if [ $RLX -eq 1 ] ; then
   echo "**Setting up relaxation masks"
   ${pget} $BASEDIR/relax/${E}/relax_rmu.a relax.rmu.a  || tellerror "Could not get relax.rmu.a"
   ${pget} $BASEDIR/relax/${E}/relax_rmu.b relax.rmu.b  || tellerror "Could not get relax.rmu.b"
fi
if [ $TRCRLX -eq 1 ] ; then
   echo "**Setting up relaxation masks"
   ${pget} $BASEDIR/relax/${E}/relax_rmu.a relax.rmutr.a  || tellerror "Could not get relax.rmu.a"
   ${pget} $BASEDIR/relax/${E}/relax_rmu.b relax.rmutr.b  || tellerror "Could not get relax.rmu.b"
fi

#
# - thermobaric reference state?
# 
echo "**Setting up various forcing files"
touch tbaric.a tbaric.b
${pget} $BASEDIR/topo/tbaric.a tbaric.a  || tellwarn "Could not get tbaric.a"
${pget} $BASEDIR/topo/tbaric.b tbaric.b  || tellwarn "Could not get tbaric.b"

#
# - Mystery field ....
#
#${pget} $BASEDIR/relax/${E}/iso_sigma.a iso.sigma.a  || tellwarn "Could not get iso.sigma.a"
#${pget} $BASEDIR/relax/${E}/iso_sigma.b iso.sigma.b  || tellwarn "Could not get iso.sigma.b"

#
# - 
#
${pget} ${BASEDIR}/force/other/iwh_tabulated.dat . || tellerror "Could not get iwh_tabulated.dat"

#C
#touch iso.top.a
#touch iso.top.b
#if (-z iso.top.a) then
#   ${pget} ${D}/../../relax/${E}/iso_top.a iso.top.a  &
#endif
#if (-z iso.top.b) then
#   ${pget} ${D}/../../relax/${E}/iso_top.b iso.top.b  &
#endif
#C
#touch thkdf4.a
#touch thkdf4.b
#if ( -z thkdf4.a ) then
#   ${pget} ${D}/../../relax/${E}/thkdf4.a thkdf4.a || tellwarn "Could not get thkdf4.a"
#endif
#if (-z thkdf4.b) then
#   ${pget} ${D}/../../relax/${E}/thkdf4.b thkdf4.b  &
#   ${pget} ${D}/../../relax/${E}/clim_tran.txt  clim_tran.txt || tellwarn "Could not get clim_tran.txt"
#   ${pget} ${D}/../../relax/${E}/cm_tran* . || tellwarn "Could not get clim_tran.txt"
#endif

${pget} ${D}/../../relax/${E}/thkdf4.a thkdf4.a || tellwarn "Could not get thkdf4.a"
${pget} ${D}/../../relax/${E}/thkdf4.b thkdf4.b || tellwarn "Could not get thkdf4.b"
${pget} ${D}/../../relax/${E}/clim_tran.txt  clim_tran.txt || tellwarn "Could not get clim_tran.txt"
${pget} ${D}/../../relax/${E}/cm_tran* . || tellwarn "Could not get clim_tran.txt"

# Retrieve infile_gp if gp is active
if [ "$gpflag" == "T" ] ; then
   ${pget} ${P}/infile_gp.in . || tellerror "Could not get infile_gp.in"
fi

# Retrieve nesting.in if outer nesting is active
# if [ "$nestoflag" == "T" ] ; then
if [ "$nestoflag" == "T" -a "${waves_only}" -eq 0 ] ; then
   #only create nesting if not using waves-only mode [TW: 20131115]

   echo "**Setting up outer nesting"

   # The nesting files will be dumped in a subdirectory of the outer model scratch dir. This 
   # is a new approach. After the run is finished, data should be copied to the data dir. The naming
   # scheme of the output directory is ./nest_out_${IR}_${IT}/
   ${pget} ${P}/nesting.in . || tellerror "Could not get nesting.in"

   # Due to the naming scheme we can retrieve the inner region and topo version like this
   cat nesting.in | \
   while read NESTDIR ; do
      #get the last name part of the path i.e. nest_out_${IR}_${IT} and etract IR and IT
      #multiple path is possible
      #example: /work/dany/FR1b0.12/expt_01.1/nest_out_FR1b0.12_01/
      II=`echo $NESTDIR | awk --field-separator=/  '{print $(NF-1)}'`
      IT=`echo $II | awk --field-separator=_  '{print $(NF)}'`
      IR=`echo $II | awk --field-separator=_  '{print $(NF-1)}'`

      # Set up nest dir
      mkdir -p $NESTDIR || tellerror "Could not create directory $NESTDIR"
      ${pget} $BASEDIR/nest_nersc/${E}/outer/depth_${IR}_${IT}.b $NESTDIR/regional.depth.b  || \
         tellerror "Could not get $BASEDIR/nest_nersc/${E}/outer/depth_${IR}_${IT}.b"
      ${pget} $BASEDIR/nest_nersc/${E}/outer/depth_${IR}_${IT}.a $NESTDIR/regional.depth.a  || \
         tellerror "Could not get $BASEDIR/nest_nersc/${E}/outer/depth_${IR}_${IT}.a"
      ${pget} $BASEDIR/nest_nersc/${E}/outer/regional.grid_${IR}_${IT}.a $NESTDIR/regional.grid.a  || \
         tellerror "Could not get $BASEDIR/nest_nersc/${E}/outer/regional.grid_${IR}_${IT}.a"
      ${pget} $BASEDIR/nest_nersc/${E}/outer/regional.grid_${IR}_${IT}.b $NESTDIR/regional.grid.b  || \
         tellerror "Could not get $BASEDIR/nest_nersc/${E}/outer/regional.grid_${IR}_${IT}.b"
      ${pget} $BASEDIR/nest_nersc/${E}/outer/${IR}_${IT}_pivots.uf $NESTDIR/pivots.uf || \
         tellerror "Could not get $BASEDIR/nest_nersc/${E}/outer/${IR}_${IT}_pivots.uf"
   done
      
fi


# Retrieve inner nesting files if active
# if [ "$nestiflag" == "T" ] ; then
if [ "$nestiflag" == "T" -a "${waves_only}" -eq 0 ] ; then
   #don't need nesting if not using waves-only mode
   echo "**Setting up inner nesting"

   # We need one output nesting directory - to be linked in as "nest" in SCRATCH dir
   [ ! -d $BASEDIR/expt_${X}/nest_out_${R}_${T}  -a ! -L $BASEDIR/expt_${X}/nest_out_${R}_${T}  ]  && \
      tellerror "No nest input directory $BASEDIR/expt_${X}/nest_out_${R}_${T} exists for this model"

   [ -L nest ] && rm nest  
   ln -s $BASEDIR/expt_${X}/nest_out_${R}_${T}/  nest

#   [ ! -d nest ] && mkdir nest  
#   [ ! -d nest ] &&  tellerror "Could not create nesting subdirectory"
   ${pget} ${BASEDIR}/nest_nersc/$E/inner/ports.input.nest ports.input || tellerror "Could not get nesting ports.input"
   ${pget} ${BASEDIR}/nest_nersc/$E/inner/rmu.a nest/ || tellerror "Could not get nesting rmu.a"
   ${pget} ${BASEDIR}/nest_nersc/$E/inner/rmu.b nest/ || tellerror "Could not get nesting rmu.b"

#   # Also copy regional.grid.[ab] and regional.depth.[ab] to nesting directory
#   ${pget} $BASEDIR/topo/regional.grid.a nest/regional.grid.a || telerror "no grid file regional.grid.a" 
#   ${pget} $BASEDIR/topo/regional.grid.b nest/regional.grid.b || telerror "no grid file regional.grid.a" 
#   ${pget} $BASEDIR/topo/depth_${R}_${T}.a nest/regional.depth.a || telerror "no topo file depth_${R}_${T}.a" 
#   ${pget} $BASEDIR/topo/depth_${R}_${T}.b nest/regional.depth.b || telerror "no topo file depth_${R}_${T}.b" 


   # Number of files in nest dir
   numnest=$(ls -1 nest/ | grep \.hdr$ | wc -l)

   if [ $numnest -eq 0 ]  ; then
      tellerror "nest directory does not contain .hdr files"
      
   else

      # Sample file
      filenest=$(ls -1 nest/ | grep \.hdr$ | head -n1)

      # Very crude sanity check on input nest dir
      nestk=$(head -n 1 nest/${filenest}  | sed "s/.*://" | awk '{printf("%d", $2)}')

      [ $nestk -ne $K ] && tellerror "Nest files and blkdat.input have different vertical dim $K vs $nestk "
   fi

fi


# Retrieve nersc tidal data set if active
if [ "$tideflag" == "T" ] ; then
   ${pget} ${BASEDIR}/tides_nersc/$E/${tidechoice}obc_elev.dat . || tellerror "Could not get tidal data ${tidechoice}obc_elev.dat "
fi
   



echo
#
# --- try to get restart input from various areas
#
if [ $YRFLAG -ne 3 -a $init -eq 1 ] ; then
   echo "No restart needed"
elif [ ${waves_only} -eq 1 ] ; then
   echo "**Waves-only mode: no restart needed"
else
   filename=${rungen}restart${refyear}_${iday}_${ihour}

   # Special case if CURVIINT flag set
   if [ -f $P/CURVIINT ] ; then
      if [ -f $BASEDIR/curviint/$E/${filename}.a -a $BASEDIR/curviint/$E/${filename}.b ] ; then
         echo "using CURVIINT restart files"
         for i in $BASEDIR/curviint/$E/${filename}* ; do
            echo "Using restart file $i"
            nname=$(basename $i)
            cp $i $nname
         done
         rm $P/CURVIINT
      else
         tellerror "CURVIINT set but could not find CURVIINT restart file ${filename}.[ab]"
      fi

   # Try to fetch restart from data dir $D
   elif [ -f $D/${filename}.a -a $D/${filename}.b ] ; then
      echo "using restart files ${filename}.[ab] from data dir $D"
      if [ -L $D/${filename}.a -a $S/${filename}.b ]; then
         echo "The file is a link, we copy" 
         cp $D/${filename}* .
      else
         ${pget} $D/${filename}* .
      fi
   # - scratch is the last place we look
   elif [ -f $S/${filename}.a -a $S/${filename}.b ] ; then
      echo "Using restart file ${filename}.[ab]  already in scratch dir $S"
        

   else
      tellerror "Could not find restart file ${filename}.[ab] (for init make sur yrflag < 3)"
   fi

fi




#
# --- model executable
#
echo "Getting executable ${EXEC} from Build"
rm -f hycom
/bin/cp ${BASEDIR}/Build_V${V}_X${X}/hycom . || tellerror "Could not get hycom executable"



#
# --- summary printout
#
[ ! -d old ] && mkdir -p old
touch   summary_out
/bin/mv summary_out old/summary_old


#
# --- heat transport output
#
touch   flxdp_out.a flxdp_out.b
/bin/mv flxdp_out.a old/flxdp_old.a
/bin/mv flxdp_out.b old/flxdp_old.b
#
touch   ovrtn_out
/bin/mv ovrtn_out old/ovrtn_old
#
# --- clean up old archive files, typically from batch system rerun.
#


if [ ! -d KEEP ]
then
    mkdir KEEP
fi
touch archv.dummy.b
for f  in arch* ; do
  /bin/mv $f KEEP/$f
done

# copy archv.extract from data to SCRATCH
# - if using ARCHIVE_SELECT flag,
#   this determines fields to dump
if [ -f $D/archv.extract ]
then
   echo cp $D/archv.extract $S
   cp $D/archv.extract $S
fi


#
# --- let all file copies complete.
#
wait
#
# --- zero file length means no rivers.
#
if [ ! -s forcing.rivers.a ] ; then
   /bin/rm forcing.rivers.[ab]
fi


#C
#C --- Nesting input archive files for next segment.
#C
#if (-e ./nest) then
#  cd ./nest
#  touch archv_${NB}.tar
#  if (-z archv_${NB}.tar) then
#    ${pget} ${D}/nest/archv_${NB}.tar archv_${NB}.tar &
#  endif
#  cd ..
#endif
#C
#chmod ug+x hycom
#/bin/ls -laFq
#C
#if (-e ./nest) then
#  ls -laFq nest
#endif


echo
if [ $numerr -ne 0 ] ; then
   echo "Some errors were fatal - rectify. You may find more info in these files"
   echo "$logfile"  
   echo "$logfile.err"
else
   echo "Some nonfatal errors errors occured. You may find more info here:"
   echo "$logfile"  
   echo "$logfile.err"
   echo
   echo "Model now set up to run in $S"
fi

# Tell where stuff ended up

exit $numerr # Fails if any fatal errors occured



