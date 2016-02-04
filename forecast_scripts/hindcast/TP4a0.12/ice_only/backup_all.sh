# backup all hindcasts:
source $SWARP_ROUTINES/source_files/hex_vars.src

R=TP4a0.12
rungen=${R:0:3}
HCtype=ice_only
HCtype2=$HCtype
#
print_info=1
#
hd=$SWARP_ROUTINES/forecast_scripts/hindcast/$R
user=timill #TODO source hidden
hd2=$hd/$HCtype
outdir0=$TP4_REALTIME/../results_hindcasts/$R/$HCtype

migdir0=/migrate/timill/RESULTS/$R/SWARP_hindcasts
mkdir -p $migdir0
migdir0=$migdir0/$HCtype
mkdir -p $migdir0


for YEAR in 2015 2016
do
   outdir=$outdir0/${YEAR}_GOOD
   cd $outdir

   fp=final_output
   mkdir -p $fp

   migdir=$migdir0/$YEAR
   mkdir -p $migdir

   for subdir in ${YEAR}_*
   do
      echo $subdir
      chmod -R g+r $subdir
      cp $subdir/$fp/*.nc $fp
      tfil=$subdir.tar.gz 
      tar -zcvf $tfil $subdir
      cp $tfil $migdir
      echo " "
   done

   chmod -R g+r $fp
   chmod -R o+r $fp
   Tfil=${fp}_$YEAR.tar.gz
   tar -zcvf $Tfil $fp
done
