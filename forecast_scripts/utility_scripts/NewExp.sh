# run from regional directory eg TP4a0.12
exp=1
lnk=1
Bld=1

if [ $# -ne 2 ] ; then
   echo "This script will create a new experiment subdirectory based on"
   echo "an already existing experiment directory. Some basic variables"
   echo "will be set in files EXPT.src and blkdat.input - the rest must be "
   echo "set by user. Script will stop if the new experiment directory "
   echo "exists."
   echo "Usage:"
   echo "   $(basename $0) old_experiment_number new_experiment_number"
   echo
   echo "Example:"
   echo "   $(basename $0) 010 011"
   echo "   makes expt_01.1 from expt_01.0"
   echo
   exit 1
fi


P=`pwd`
in1=$1
Xold=${in1:0:2}.${in1:2:1}
eold=$P/expt_$Xold
E=$2
X=${E:0:2}.${E:2:1}
enew=$P/expt_$X

if [ ! -d $eold ]
then
   echo " "
   echo "Directory $eold does not exist"
   echo "- launch from root model directory eg 'TP4a0.12'"
   echo " "
   exit
else
   source REGION.src       # need R variable
   rungen=${R:0:3}
fi

if [ $rungen == 'TP4' ]
then
   nmpi=133
elif [ $rungen == 'BS1' ]
then
   nmpi=112
elif [ $rungen == 'FR1' ]
then
   nmpi=51
fi
JOBNAME=${rungen}HC$E

mkdir -p $enew
mkdir -p $enew/data

if [ $exp -eq 1 ]
then
   # copy linked dirs:
   cp -a $eold/subprogs $enew
   cp -a $eold/bak      $enew

   # copy files:
   # NB still need to change blkdat.input and EXPT.src
   for f in $eold/*
   do
      if [ -f $f ]
      then
         cp $f $enew
      fi
   done

   cd $enew
   mkdir -p log

   # Change EXPT.src
   mv EXPT.src EXPT.src.tmp
   cat EXPT.src.tmp | \
      sed "s/^X=.*/X=\"$X\"                                                          # X based on dir name (expt_$X)/" | \
      sed "s/^E=.*/E=\"$E\"                                                           # E is X without \".\"/  " \
      > EXPT.src
   rm EXPT.src.tmp

   # Set up new blkdat.input
   mv blkdat.input blkdat.input.tmp
   ss=`cat blkdat.input.tmp | grep iexpt`
   cat blkdat.input.tmp | \
      sed "s/$ss/ $E\t  'iexpt ' = experiment number x10/" \
      > blkdat.input
   rm blkdat.input.tmp

   # change pbsjob.sh
   cat $FCcommon/inputs/pbsjob.sh | sed \
        -e "s/JOBNAME/$JOBNAME/g" \
        -e "s/MPPWIDTH/$nmpi/g" \
        -e "s/exit \$?//g" \
        > pbsjob.sh

   # nesting
   if [ -d $eold/nest_out* ]
   then
      ln -s $eold/nest_out* .
   fi
   
fi

if [ $lnk -eq 1 ]
then
   # make links
   cd $P/force/rivers
   ln -sf $1 $E
   pwd
   ls -lh
   echo " "

   cd $P/force/nersc_era40
   ln -sf $1 $E
   pwd
   ls -lh
   echo " "

   cd $P/relax
   ln -sf $1 $E
   pwd
   ls -lh
   echo " "

   cd $P/nest_nersc
   ln -sf $1 $E
   pwd
   ls -lh
   echo " "

   cd $P/tides_nersc
   ln -sf $1 $E
   pwd
   ls -lh
   echo " "
fi

if [ $Bld -eq 1 ]
then
   # set up Build
   cd $P
   source REGION.src       # need R variable
   source $enew/EXPT.src   # need V variable

   B1=`readlink -f Build_V${V}_X$Xold`
   B2=`readlink -f Build_V${V}_X$X`
   mkdir -p $B2
   cd $B2

   # 
   cp    $B1/flags         $B2
   cp    $B1/Makefile      $B2
   cp    $B1/dependencies  $B2
   cp    $B1/setuppatch.sh $B2
   cp -a $B1/svn_Build     $B2

   if [ -f $B1/dimensions_nersc.h ]
   then
      cp    $B1/dimensions_nersc.h $B2
   else
      # run setuppatch.sh (needs blkdat.input correct)
     ./setuppatch $nmpi
   fi

   # cp svn_Build/* .
   pwd
   ls -lh
   echo " "

fi

# # compiler linked files
# CLF=/home/nersc/timill/hycom/hycom-TWaddons/compiler_linked_files/$R
# if [ ! -d "$CLF/X01.$2" ]
# then
#    echo cp -r $CLF/X01.$1 $CLF/X01.$2
#    cp -r $CLF/X01.$1 $CLF/X01.$2
# fi

# cd $P/expt_01.$2
# pwd
# echo ln -s $CLF/X01.$2 bak
# ln -s $CLF/X01.$2 bak
# source EXPT.src



echo ""
echo "Also (TW comments):"
echo "1. Now need to get restart files and maybe nesting_out directory"
echo "2. Check links to flags, dependencies and Makefile are OK"
