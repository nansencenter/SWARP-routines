#!/bin/bash
#Get latest restart from the internal repo to the working dir
#Interpolate with curviint

print_info=1 # print info for debugging

# ====================================================================================
source $SWARP_ROUTINES/source_files/hex_vars.src
THIS_SRC=`readlink -f $1`
source $THIS_SRC
# if [ $print_info -eq 1 ]
# then
#    echo ""
#    cat $THIS_SRC
#    echo ""
# fi
# ====================================================================================

# ====================================================================================
# EMAIL ADDRESS
email=$(cat $FCemail)
# ====================================================================================

# DIRECTORIES AND DATELIST
datelist=$THISFC/logs/datelist.txt
xdir=$RTmods/$HYCOMreg/expt_01.$Xno
ddir=$xdir/data   # where forecast will be done

# TP4 locations
rdir=/migrate/timill/restarts/TP4a0.12/SWARP_forecasts   # directory with TP4 restarts
Ddir=$TP4_REALTIME/expt_01.1/data # location of TP4a0.12 directory (where forecast will be done)

# info dir
tday=$(cat $datelist | sed '1!d')
tyear=`date --date=$tday "+%Y"`
pyear=$(expr $tyear - 1)            # previous year
tday_j=`date --date=$tday "+%j"`
tday_j=$((tday_j-1))                # current day julian (1 Jan = 0)
idir=$THISFC2/$tday/info
mkdir -p $idir

# TEXTFILE AND LOG
logdir=$THISFC/logs
mkdir -p $logdir
out_restart=$logdir/last_restart.txt
log=$logdir/tp_get_log.txt
touch $log
echo $tday  >> $log
echo ""     >> $log

# ====================================================================================
# first look in data directory for last Monday's restart
# - this avoids curviint if possible
lastMon=`date --date="last Monday" "+%Y%m%d"`
if [ $lastMon -ge $tday ]
then
   lastMon=`date --date="$lastMon -7days" "+%Y%m%d"`
fi
lastMon_y=`date --date="$lastMon" "+%Y"`
lastMon_j=`date --date="$lastMon" "+%j"`
lastMon_j=$((lastMon_j-1))
rfil0=${rungen}restart${lastMon_y}_${lastMon_j}_00

if [[ -f $ddir/$rfil0.a && -f $ddir/$rfil0.b && -f $ddir/${rfil0}ICE.uf ]]
then
   echo Found the latest restart in data: $rfil0
   echo  $rfil0                > $out_restart
   cp    $out_restart            $idir
   #
   echo "Restart already in $ddir"  >> $log
   echo $rfil.a                     >> $log
   echo $rfil.b                     >> $log
   echo ${rfil}ICE.uf               >> $log
   echo ""                          >> $log
   exit
fi
# ====================================================================================


# ====================================================================================
# next look in TP4 data directory for last Monday's restart
# - then run curviint
lastMon=`date --date="last Monday" "+%Y%m%d"`
if [ $lastMon -ge $tday ]
then
   lastMon=`date --date="$lastMon -7days" "+%Y%m%d"`
fi
lastMon_y=`date --date="$lastMon" "+%Y"`
lastMon_j=`date --date="$lastMon" "+%j"`
lastMon_j=$((lastMon_j-1))
Rfil0=TP4restart${lastMon_y}_${lastMon_j}_00

if [[ -f $Ddir/$Rfil0.a && -f $Ddir/$Rfil0.b && -f $Ddir/${Rfil0}ICE.uf ]]
then
   echo Found the latest restart in data: $Rfil0
   echo  $Rfil0         > $out_restart
   cp    $out_restart   $idir
   #
   echo "TP4 restart located in $Ddir" >> $log
   echo "Running curviint..."          >> $log

   # Run curviint
   X=01.$Xno
   E=01$Xno
   cd $xdir/..
   pwd
   inputs="$X 01.1 $TP4_REALTIME $Ddir/$Rfil0.a"
   echo "other_nersc/bin/curviint.sh $inputs"
   other_nersc/bin/curviint.sh $inputs

   # get results, rename and move:
   cdir=`readlink -f curviint/$E`
   cp $cdir/$Rfil0.a $ddir/$rfil0.a
   cp $cdir/$Rfil0.b $ddir/$rfil0.b
   cp $cdir/${Rfil0}ICE.uf $ddir/${rfil0}ICE.uf
   echo ""
else
   echo "Latest TP4 restart not present"  >> $log
   echo "$Rfil0"                          >> $log
   echo "Latest TP4 restart not present"
   echo "$Rfil0"
fi
# ====================================================================================
