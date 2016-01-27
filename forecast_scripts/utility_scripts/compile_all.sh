# compile all FC models

source $SWARP_ROUTINES/source_files/hex_vars.src

if [ $# -gt 0 ]
then
   # restore deleted files in Build
   # - this also gets flags from the proper "inputs"
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_build.sh TP4 ice_only            || echo "Build set-up error: TP4 ice_only          "
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_build.sh TP4 wavesice_ww3arctic  || echo "Build set-up error: TP4 wavesice_ww3arctic"
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_build.sh BS1 ice_only            || echo "Build set-up error: BS1 ice_only          "
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_build.sh BS1 wavesice_ww3arctic  || echo "Build set-up error: BS1 wavesice_ww3arctic"
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_build.sh FR1 ice_only            || echo "Build set-up error: FR1 ice_only          "
   $SWARP_ROUTINES/forecast_scripts/utility_scripts/restore_build.sh FR1 wavesice_ww3arctic  || echo "Build set-up error: FR1 wavesice_ww3arctic"
fi


# for Xno in 1 3
for Xno in 3
do
   cd $TP4_REALTIME/Build_V2.2.12_X01.$Xno
   pwd
   make clean
   make || echo "compilation error"
done

for Xno in 0 1
do
   cd $BS1_REALTIME/Build_V2.2.12_X01.$Xno
   pwd
   make clean
   make || echo "compilation error"

   cd $FR1_REALTIME/Build_V2.2.12_X01.$Xno
   pwd
   make clean
   make || echo "compilation error"
done
