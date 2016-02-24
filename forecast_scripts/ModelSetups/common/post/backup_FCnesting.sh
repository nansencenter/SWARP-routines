# Run from hindcast to ackup nested models (BS1a0.045 and FR1a0.03) nesting files for from prev.
# Monday to Sunday (7days)!
# INPUT
# $1=$THIS_SRC (link to "This forecast source file")
# $2=rdate YYYYMMDD "start date"
# $2=fdate YYYYMMDD "final date"

huser=timill
# Source model and forecast settings/links 
if [ ! $# -eq 3 ]
then
   echo "Backup of nesting files need 3 inputs 1) THIS_SRC 2) rdate 3) fdate "
   exit
fi
source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC
# Hexgon dir
keyname=/home/nersc/timill/.ssh/key_hex2store
nestdir=$RTmods/$HYCOMreg/expt_01.$Xno/SCRATCH/nest
# Norstore dir
commonstore=/norstore_osl/home/timill/SWARP-routines/forecast_scripts/ModelSetups/common/NorStore
destdir=/projects/NS2993K/SWARP_FC/nesting/$HYCOMreg/$FCtype
datelist=$THISFC/logs/datelist.txt
# read input start and end date
rdate=$2
fdate=$3
# Set start and end date, year and topaz year day
rdate=`date -d "$rdate" +%Y%m%d`
fdate=`date -d "$fdate" +%Y%m%d`
# OUTPUT TAR FILE
nestout=FC_${HYCOMreg}_${rdate}_${fdate}.tar.gz
# Gether nesting files from one week into "nestlist", loop over all days
cd $nestdir
nestlist=()
cdate=`date -d "$rdate" +%Y%m%d`
while [[ $cdate -le $fdate ]]; do
   echo $cdate
   year=`date -d "$cdate" +%Y`
   YD=`date -d "$cdate" +%j`
   # Topaz year days start at 000, so -1 from YD
   TYD=`expr $YD - 1`
   TYD=`printf %03d $TYD`
   # echo "year TYD = $year $TYD"
   addlist=()
   addlist=$(ls nest_${year}_${TYD}*)
   # echo "addlist: ${addlist[*]})"
   if [ -n "$addlist" ]; then
      nestlist+=("${addlist[*]}")
   else
      errormessage= "!! No nesting files for TOPAZ day $year $TYD, though continue !!"
      echo $errormessage 1>&2
   fi
   cdate=`date -d "$cdate + 1 days" +%Y%m%d`
done
if [ -n "$nestlist" ]; then
   echo "tar zcvfh $nestout nestlist"
# echo "   tar zcvfh $nestout ${nestlist[*]}"
else
   errormessage="!! No nesting files at all for $dm9l to $dm2l !!"
   echo $errormessage 1>&2 
   # exit
fi
# check if remote directory exist, if not make one
echo "ssh -i $keyname $huser@login3.norstore.uio.no "sh $commonstore/swarp_fc_mkdir.sh $HYCOMreg $FCtype nesting"
ssh -i $keyname $huser@login3.norstore.uio.no "sh $commonstore/swarp_fc_mkdir.sh $HYCOMreg $FCtype nesting""
# send to Norstore /projects..../nesting  directory
scp -i $keyname $nestout $huser@login3.norstore.uio.no:$destdir
#echo "scp -i $keyname $nestout $huser@login3.norstore.uio.no:$destdir"
# remove .tar file from /work
rm $nestdir/$nestout

