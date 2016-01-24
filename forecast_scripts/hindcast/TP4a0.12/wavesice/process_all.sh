# script to set-up all experiments
# needs to be run from crontab since nesting is done to same folder
# -> can be crashes if files are opened by 2 expts at the same time
source $SWARP_ROUTINES/source_files/hex_vars.src

R=TP4a0.12
rungen=${R:0:3}
HCtype=wavesice
HCtype2=$HCtype
#
print_info=1
#
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/$R
user=timill #TODO source hidden
hd2=$hd/$HCtype
outdir0=$TP4_REALTIME/../results_hindcasts/$R/$HCtype
mkdir -p $outdir0


P=$TP4_REALTIME
cd $P
refno=4
xref=$P/expt_01.$refno           # reference expt directory  (has restarts)
bref=$P/Build_V2.2.12_X01.$refno # reference Build directory (has hycom executable compiled)

E0=$((1$refno-1)) # increment this
X=01.$refno
Eno=0$E0
rno=0

md=/migrate/timill/restarts/$R/SWARP_forecasts
list=(`ls $md/2015/*gz`)
list2=(`ls $md/2016/*gz`)
Nlist=${#list[@]}
inext=Nlist
for loop_i in `seq 1 ${#list2[@]}`
do
   list[$inext]=${list2[$((loop_i-1))]}
   inext=$((inext+1))
   Nlist=$((Nlist+1))
done

logdir=$hd2/logs
mkdir -p $logdir
log=$logdir/hc_proc.log
rm -f $log
echo `date` >> $log

#########################################################
# ECMWF files in 2015 skip between
# Saturday 22-Aug 12:00 and Sunday 23-Aug 06:00
# ie 4 missing records
bad_dates[0]=20150817 #Monday before this
bad_error[0]="Gap in ECMWF forcing 20150822-23"

# don't want to do last Monday's version 
# - winds are forecasts part of the week
bad_dates[1]="`date --date="last Monday" +%Y%m%d`"
bad_error[1]="Last Monday"

bad_dates[2]=20150601
bad_error[2]="Possible problem with restart"

# no WAM waves: 20150111 - no restart for then
# no WAM waves: 20150406
bad_dates[3]=20150406
bad_error[3]="No WAM waves on 20150406"

# no WAM waves: 20150612
bad_dates[4]=20150608
bad_error[4]="No WAM waves on 20150612"

# no WAM waves: 20151114
bad_dates[5]=20151109
bad_error[5]="No WAM waves on 20151114"

Nbad=${#bad_dates[@]}
# for loop_i in `seq 1 $Nbad`
# do
#    echo ${bad_dates[$((loop_i-1))]}
# done
#########################################################

if [ 1 -eq 1 ]
then
   # DO ALL
   E0_start=14
   E0_stop=2000
else
   # DO SOME
   # E0_start=36
   # E0_stop=56
   E0_start=57
   E0_stop=2000
fi
HYC2PROJ=1
NCMERGE=1

DONE=0
for rno in `seq 1 $Nlist`
do
   tfil0=${list[$((rno-1))]}
   E0=$((E0+1))
   Eno=0$E0
   X=0${E0:0:1}.${E0:1:1}
   xdir=$P/expt_$X            # reference expt directory  (has restarts)
   ddir=$xdir/data

   # echo `basename $tfil0` $E0
   # continue

   if [ $E0 -lt $E0_start ] || [ $E0 -gt $E0_stop ]
   then
      echo "$E0: skipping"
      continue
   else
      echo "$E0: TODO"
   fi


   ##############################################################################
   # hyc2proj inputs
   cd $ddir
   cp $xdir/../topo/grid.info .
   cp $hd2/inputs/extract.archv_wav .
   cp $hd2/inputs/extract.daily .
   HDIR=$SWARP_ROUTINES/netcdf_production/Input
   cp $HDIR/depthlevels.in .
   cp $HDIR/proj.in.polar_stereographic.$R proj.in
   ##############################################################################


   ##############################################################################
   # DATE
   tfil=`basename $tfil0`
   ryear=${tfil:10:4}
   rday=${tfil:15:3}
   rday2=10#$rday #base 10
   fday2=$((rday2+7))
   fday=`printf %3.3d $fday2`
   rdate=`date --date="$ryear-01-01 +${rday}days" "+%Y%m%d"`
   rdate_long=`date --date="$rdate" "+%Y-%m-%d"`
   echo $tfil $ryear $rday $fday $E0

   outdir=${outdir0}/${ryear}_GOOD
   mkdir -p $outdir

   BAD_DAY=0
   for nbad in `seq 0 $((Nbad-1))`
   do
      if [ $rdate -eq ${bad_dates[$nbad]} ]
      then
         BAD_DAY=1
         BAD_ERROR=${bad_error[$nbad]}
         break
      fi
   done

   if [ $BAD_DAY -eq 1 ]
   then
      echo " "
      echo "Not launching week ${ryear}_${rday} ($rdate)"
      echo $BAD_ERROR
      echo " "
      continue
   fi
   ##############################################################################


   #############################################################
   # make directories for results
   OUTDIR=$outdir/${ryear}_${rday}
   mkdir -p $OUTDIR
   # echo $OUTDIR
   # continue

   BDIR=$OUTDIR/binaries
   mkdir -p $BDIR
   mkdir -p $BDIR/DAILY
   mkdir -p $BDIR/archv
   mkdir -p $BDIR/archv_wav

   NDIR=$OUTDIR/netcdf
   mkdir -p $NDIR
   mkdir -p $NDIR/DAILY
   mkdir -p $NDIR/archv
   mkdir -p $NDIR/archv_wav
   #############################################################

   # echo CONTINUE;continue
   dinf=$OUTDIR/info
   mkdir -p $dinf
   echo cp $hd2/inputs/flags    $dinf
   cp $hd2/inputs/flags    $dinf

   OP="cp"
   echo $OP $xdir/log/mpijob.out $dinf
   $OP $xdir/log/mpijob.out $dinf

   if [ $HYC2PROJ -eq 1 ]
   then
      #############################################################
      # run hyc2proj & move files to output location
      for afil in ${rungen}DAILY*.a
      do
         hyc2proj $afil
         $OP $afil $BDIR/DAILY

         bfil=${afil%.a}.b
         $OP $bfil $BDIR/DAILY

         mv ${rungen}DAILY*.nc $NDIR/DAILY
      done


      for afil in ${rungen}archv.*.a
      do
         hyc2proj $afil
         $OP $afil $BDIR/archv

         bfil=${afil%.a}.b
         $OP $bfil $BDIR/archv

         mv ${rungen}archv*.nc $NDIR/archv
      done


      for afil in ${rungen}archv_wav.*.a
      do
         hyc2proj $afil
         $OP $afil $BDIR/archv_wav

         bfil=${afil%.a}.b
         $OP $bfil $BDIR/archv_wav

         mv ${rungen}archv_wav*.nc $NDIR/archv_wav
      done
   fi # HYC2PROJ
   #############################################################
   

   #############################################################
   # MAKE FINAL NETCDF FILE
   if [ $NCMERGE -eq 1 ]
   then
      if [ $HCtype == "ice_only" ]
      then
         # merge archv nc files
         ADIR=$NDIR/archv
      else
         # merge archv_wav nc files
         ADIR=$NDIR/archv_wav
      fi

      # Unpack nc files so merge is done correctly
      cd $ADIR
      rm -r tmp.nc tmp
      mkdir -p tmp
      for f in *.nc
      do
         echo "Unpacking $f"
         ncpdq -U $f tmp/$f
         cd tmp

         # set missing value to be -32767
         ncatted -O -a _FillValue,,o,f,-32767 $f
         ncatted -O -a missing_value,,o,f,-32767 $f

         cd $ADIR
      done

      # do merge
      ncrcat tmp/*.nc tmp.nc

      ofil=SWARP_hindcast_${HCtype2}_start${rdate}T000000Z.nc

      # make final file by repacking tmp.nc
      echo " "                                              >> $log
      echo "Making $ofil (ncpdq)..."                        >> $log
      ncpdq tmp.nc $ofil
      
      # clean and move to final_output
      rm -r tmp tmp.nc
      odir=$OUTDIR/final_output # weekly file of merged netcdfs
      mkdir -p $odir
      mv $ofil $odir
      cd $odir
      #############################################################


      #########################################################################
      # FINALLY ADD INFO

      # change variable names to standard abbreviations in GRIB2 tables
      ncrename -v fice,icec   $ofil #ice cover
      ncrename -v hice,icetk  $ofil #ice thickness

      # change time:units
      ncatted -O -h -a units,time,m,c,"seconds since 1970-01-01T00:00:00Z"  $ofil

      #####################################################################################
      # most variable attributes set in hyc2proj:
      # - just add shape of earth to "int stereographic"
      #   (radius from hyc2proj - mod_toproj.F90)
      ncatted -O -a "semi_major_axis",stereographic,c,f,6378273.     $ofil
      ncatted -O -a "semi_minor_axis",stereographic,c,f,6378273.     $ofil
      #####################################################################################

      #####################################################################################
      # Global attributes for  files
      # (o or c = overwrite/create)
      ncatted -O -h -a software_version,global,c,c,"NERSC-HYCOM (TOPAZ)"            $ofil
      ncatted -O -h -a references,global,c,c,"www.nersc.no"                         $ofil
      ncatted -O -h -a comment,global,c,c," "                                       $ofil

      if [ $rungen == "TP4" ]
      then
         ncatted -O -h -a area_name,global,c,c,"TP4a0.12"                                 $ofil
         ncatted -O -h -a area_resolution,global,c,c,"12.5km"                             $ofil
         ncatted -O -h -a area_description,global,c,c,"Arctic and North Atlantic Oceans"  $ofil
      elif [ $rungen == "BS1" ]
      then
         ncatted -O -h -a area_name,global,c,c,"BS1a0.045"                    $ofil
         ncatted -O -h -a area_resolution,global,c,c,"4.5km"                  $ofil
         ncatted -O -h -a area_description,global,c,c,"Barents and Kara Sea"  $ofil
      elif [ $rungen == "FR1" ]
      then
         ncatted -O -h -a area_name,global,c,c,"FR1a0.03"            $ofil
         ncatted -O -h -a area_resolution,global,c,c,"3.5km"         $ofil
         ncatted -O -h -a area_description,global,c,c,"Fram Strait"  $ofil
      fi

      if [ $HCtype == "wavesice" ]
      then

         ncatted -O -h -a title,global,o,c,"SWARP waves-in-ice hindcast"            $ofil # o=overwrite/create, c=format (also f=float)
         ncatted -O -h -a field_type,global,c,c,"6-hourly"                          $ofil
         ncatted -O -h -a hindcast_range,global,c,c,"7 day hindcast"                $ofil
         ncatted -O -h -a wave_forcing,global,c,c,"WAM North Sea Arctic (met.no)"   $ofil
         ncatted -O -h -a wave_forcing_contact,global,c,c,"bruce.hackett@met.no"    $ofil

      elif [ $HCtype == "wavesice_ww3arctic" ]
      then

         ncatted -O -h -a title,global,o,c,"SWARP waves-in-ice hindcast"            $ofil # o=overwrite/create, c=format (also f=float)
         ncatted -O -h -a field_type,global,c,c,"3-hourly"                          $ofil
         ncatted -O -h -a hindcast_range,global,c,c,"7 day hindcast"                $ofil
         ncatted -O -h -a wave_forcing,global,c,c,"WW3 Arctic"                      $ofil
         #
         wref="ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/wavewatch3/iowaga/HINDCAST/ARCTIC"
         ncatted -O -h -a wave_forcing_reference,global,c,c,$wref                   $ofil

      elif [ $HCtype == "ice_only" ]
      then

         ncatted -O -h -a title,global,o,c,"SWARP ice-only hindcast"                $ofil # o=overwrite/create, c=format (also f=float)
         ncatted -O -h -a field_type,global,c,c,"3-hourly"                          $ofil
         ncatted -O -h -a hindcast_range,global,c,c,"7 day hindcast"                $ofil
         ncatted -O -h -a wave_forcing,global,c,c,"None"                            $ofil

      fi

      ncatted -O -h -a hindcast_start_date,global,c,c,"${rdate_long}T00:00:00Z"      $ofil
      ncatted -O -h -a institution,global,c,c,"NERSC"                               $ofil
      ncatted -O -h -a institution_references,global,c,c,"http://www.nersc.no/"     $ofil
      ncatted -O -h -a data_centre,global,c,c,"NERSC"                               $ofil
      ncatted -O -h -a data_centre_references,global,c,c,"www.nersc.no"             $ofil
      ncatted -O -h -a contact,global,c,c,"timothy.williams@nersc.no"               $ofil
      ncatted -O -h -a project,global,c,c,"SWARP"                                   $ofil
      ncatted -O -h -a project_references,global,c,c,"swarp.nersc.no"               $ofil
      ncatted -O -h -a distribution_statement,global,c,c,"No restrictions"          $ofil
      ncatted -O -h -a operational_status,global,c,c,"test"                         $ofil
      #
      ncatted -O -h -a title,global,o,c,"SWARP waves-in-ice hindcast"               $ofil # o=overwrite/create, c=format (also f=float)
      # ncatted -O -h -a history,global,o,c,"NERSC-HYCOM output->hyc2proj->ncrcat"    $ofil

      # Restart file date
      reformat=0
      if [ $reformat -eq 1 ]
      then
         # reformat bulletin date
         # eg 2015-05-09T00:00:00Z -> 2015-05-09 UTC 00:00:00
         bdate_info=$(ncinfo $ofil | grep "bulletin_date")
         split=(${bdate_info//\ / })   # make array with delimiter space "\ "
         bdate_time=${split[1]}        # eg 2015-05-09T00:00:00Z
         split=(${bdate_time//T/ })    # make array with delimiter "T"
         bdate=${split[0]}             # eg 2015-05-09
         btime=${split[1]}             # eg 00:00:00Z
         btime=${btime:0:8}            # remove "Z"
         ncatted -O -h -a bulletin_date,global,o,c,"$bdate UTC $btime" $ofil
      fi

      # Rename bulletin_date to restart_file_date (clearer)
      ncrename -a global@bulletin_date,restart_file_date    $ofil

      # delete old attribute(s)
      ncatted -a field_date,global,d,,                      $ofil
      ncatted -a history,global,d,,                         $ofil
      #####################################################################################
   fi # NCMERGE

   DONE=$((DONE+1))
done
