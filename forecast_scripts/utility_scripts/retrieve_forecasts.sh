me=$0
if [ $# -ne 2 ]
then
   echo "`basename $me` [region] [FCtype]"
   exit
fi
reg=$1
FCtype=$2

if [ $reg == "TP4" ]
then
   wdir=$TP4_REALTIME_RES/$FCtype
   mdir=/migrate/timill/RESULTS/TP4a0.12/SWARP_forecasts/$FCtype
fi

cd $wdir
for Ydir in $mdir/*
do
   # yearly directory
   year=`basename $Ydir`
   if [ $year -lt 2016 ]
   then
      echo Skipping $year
      continue
   fi

   for Gfil in $Ydir/*
   do
      # loop over daily forecasts
      gfil=`basename $Gfil`
      gg=${gfil%*.tar.gz}
      N=${#gg}
      cdate=${gg:$((N-8)):$N}
      echo $cdate

      # if [ ${cdate:4:2} -lt 12 ]
      # then
      #    continue
      # fi

      if [ -d $cdate ]
      then
         # don't overwrite directories
         echo Skipping $cdate
         continue
      else
         echo Getting $cdate
         # continue
         mkdir $cdate
         cd $cdate
      fi

      echo cp $Gfil .
      cp $Gfil .

      echo tar -zxf $gfil
      tar -zxf $gfil
      rm $gfil
      cd $wdir
   done
done
