webdir="http://wdc.aari.ru/datasets/d0004/"
ToDir=/work/shared/nersc/msc/AARI_icecharts

for region in Gre Bar
do
   cd $ToDir
   mkdir -p $region
   cd $region

   for year in 2015 2016
   do
      wget -r -nH -e robots=off --no-parent --cut-dirs=5 -R "index.html*" -R "*.zip" $webdir/$region/sigrid/$year/
   done
done
