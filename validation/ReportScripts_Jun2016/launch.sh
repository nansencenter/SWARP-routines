# run from /work
# - can't import modules from /home
# so copy them and prepend location to $PYTHONPATH

PF=$SWARP_ROUTINES/py_funs
mkdir -p py_funs
cp $PF/*.py py_funs

# export PYTHONPATH=py_funs:$PYTHONPATH
# echo $PYTHONPATH

mkdir -p log
rm -f log/*

# N=1
N=60
for n in `seq 1 $N`
do
   qsub pbsjob.sh
done
