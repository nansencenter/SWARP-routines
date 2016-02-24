# Run on to backup nested models (BS1a0.045 and FR1a0.03) restart files for cdate in datelist!
# Run after hindcast an assumes restart files still are in /work/..../data/ 
# INPUT
# $1=$THIS_SRC (link to "This forecast source file")
# $2=cdate YYYYMMDD (the forecast day, if not defined will be set to current day)
huser=timill
keyname=/home/nersc/timill/.ssh/key_hex2store
# TOPAZ dir
resdir=$RTmods/$HYCOMreg/expt_01.$Xno/data
# Norstore dir
commonstore=/norstore_osl/home/timill/SWARP-routines/forecast_scripts/ModelSetups/common/NorStore
destdir=/projects/NS2993K/SWARP_FC/restarts/$HYCOMreg/$FCtype
#
# Source model and forecast settings/links 
source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC
eatelist=$THISFC/logs/datelist.txt
destdir=/projects/NS2993K/SWARP_FC/restarts/$HYCOMreg/$FCtype
# Use cdate in datelist (hindcast cdate==rdate)
cdate=$(cat $datelist | sed '1!d')
Y=`date -d "$cdate" +%Y`
YD=`date -d "$cdate" +%j`
# Topaz year day (start at day 000)
TYD=`expr $YD - 1`
# Filename
resfilebase=${rungen}restart${Y}_$TYD
resfile=${resfilebase}.tar.gz
cd $resdir
tar zcvf $resfile ${resfilebase}_00.a ${resfilebase}_00.b ${resfilebase}_00ICE.uf
# check if directory exist, if not make one
ssh -i $keyname $huser@login3.norstore.uio.no "sh $commonstore/swarp_fc_mkdir.sh $HYCOMreg $FCtype restarts"
# send to Norstore /projects..../restart  directory
scp -i $keyname $resfile $huser@login3.norstore.uio.no:$destdir
# remove .tar file from /work
rm $resdir/$resfile

