DMIdir="/work/shared/nersc/msc/DMI_icecharts/sigrid"
outdirA="/work/shared/nersc/msc/DMI_icecharts/MIZ_FA"
outdirB="/work/shared/nersc/msc/DMI_icecharts/MIZ_FB"
mkdir -p $outdirA $outdirB

for year in 2015 2016
do
   python read_icechart.py --indir=$DMIdir/$year --outdir=$outdirA --MIZ_criteria=FA_only --chart_source=DMI
   python read_icechart.py --indir=$DMIdir/$year --outdir=$outdirB --MIZ_criteria=FB_only --chart_source=DMI
done
