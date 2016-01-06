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
   echo "   $(basename $0) 0 1"
   echo "   makes expt_01.1 from expt_01.0"
   echo
   exit 1
fi


P=`pwd`
eold=$P/expt_01.$1
X=01.$2
E=01$2
enew=$P/expt_$X

if [ ! -d $eold ]
then
   echo "Launch from root model directory eg 'TP4a0.12'"
   exit
fi

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
   ln -sf 01$1 01$2
   pwd
   ls -lh
   echo " "

   cd $P/force/nersc_era40
   ln -sf 01$1 01$2
   pwd
   ls -lh
   echo " "

   cd $P/relax
   ln -sf 01$1 01$2
   pwd
   ls -lh
   echo " "

   cd $P/nest_nersc
   ln -sf 01$1 01$2
   pwd
   ls -lh
   echo " "

   cd $P/tides_nersc
   ln -sf 01$1 01$2
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

   B1=`readlink -f Build_V${V}_X01.$1`
   B2=`readlink -f Build_V${V}_X01.$2`
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

      if [ ${R:0:3} == "TP4" ]
      then
        ./setuppatch 133
      elif [ ${R:0:3} == "BS1" ]
      then
        ./setuppatch 112
      elif [ ${R:0:3} == "FR1" ]
      then
        ./setuppatch 51
      fi
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
echo "3. Run setuppatch.sh in Build to make dimensions_nersc.h"
echo "   eg ./setuppatch.sh 133 for TP4"
